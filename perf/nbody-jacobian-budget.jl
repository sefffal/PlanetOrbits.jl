# Go/no-go budget for an analytic-Jacobian AHL21 gradient path.
#
# Run:
#   julia --project=perf perf/nbody-jacobian-budget.jl
#
# The question this answers, *before* anyone writes the 500 lines of derivative
# algebra upstream's `compute_jacobian_gamma!` represents:
#
#   An analytic path costs (local Jacobian assembly) + (chain-rule matmuls).
#   The chain-rule matmuls are irreducible — they are the same for any
#   assembly, analytic or not. So if the matmuls *alone* already cost more
#   than ForwardDiff's entire marginal (Dual{P} step minus Float64 step), the
#   analytic route cannot win at that (N, P) and no algebra will save it.
#
# The matmul cost is measured on a synthetic propagation with the real
# structure and the real shapes — same layout, same gather/scatter, same
# sequence of blocks as `ahl21_step` — but with the local Jacobian blocks
# supplied as constants rather than computed. That is deliberately a *lower
# bound*: a free assembly.
#
# Layout of the propagated sensitivity S = ∂(x,v)/∂θ, (6N) × P:
#   rows      1 .. 3N   x components, body-major
#   rows 3N+1 .. 6N     v components, body-major
# plus dmdθ = ∂m/∂θ (N × P) for the explicit mass dependence of each block.

using StaticArrays
using LinearAlgebra: I
using BenchmarkTools
using ForwardDiff
using ForwardDiff: Dual, Partials
using Printf
using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2

const NS = (2, 3, 4, 5)          # bodies
const PS = (7, 13, 20, 29, 40)   # parameters

# ---------------------------------------------------------------- workloads
# A compact, hierarchical system in AU / day / M⊙ — the regime where AHL21 is
# actually reached for (γ ≈ 0.3–1) per step. Values only set the γ band and
# the Newton iteration counts; the cost model does not depend on them further.
function demo_state(N::Int)
    m = [1.0; fill(3e-5, N - 1)]
    x = zeros(3, N); v = zeros(3, N)
    for i in 2:N
        a = 0.1 * 1.6^(i - 2)                  # 0.1 AU and out, near-resonant
        vc = sqrt(G_au3_day2 * m[1] / a)
        ph = 0.7 * i
        x[1, i] = a * cos(ph); x[2, i] = a * sin(ph); x[3, i] = 0.01a
        v[1, i] = -vc * sin(ph); v[2, i] = vc * cos(ph); v[3, i] = 0.0
    end
    # zero the barycentre so the star carries the reflex
    for k in 1:3
        x[k, 1] = -sum(m[i] * x[k, i] for i in 2:N) / m[1]
        v[k, 1] = -sum(m[i] * v[k, i] for i in 2:N) / m[1]
    end
    P_in = 2π * sqrt(0.1^3 / (G_au3_day2 * m[1]))
    return (st = NBodyState(SMatrix{3,N}(x), SMatrix{3,N}(v)),
            Gm = SVector{N}(G_au3_day2 .* m),
            h = P_in / 20)
end

# ------------------------------------------------- synthetic chain-rule step
# One pair block: 12 state rows in, 12 out, plus the two masses of the pair.
@inline function chain_pair!(S::MMatrix{M,P,T}, idx::SVector{12,Int},
                             J::SMatrix{12,12,T}, Jm::SMatrix{12,2,T},
                             dmdθ::MMatrix{N,P,T}, i::Int, j::Int) where {M,P,T,N}
    @inbounds for c in 1:P
        t = SVector{12,T}(ntuple(k -> S[idx[k], c], Val(12)))
        u = J * t + Jm * SVector{2,T}(dmdθ[i, c], dmdθ[j, c])
        for k in 1:12
            S[idx[k], c] = u[k]
        end
    end
    return nothing
end

# The corrector: Δv depends on x and the masses only (never on v), so its
# block is (3N)×(3N) onto the velocity rows plus (3N)×N of mass columns.
@inline function chain_phi!(S::MMatrix{M,P,T}, Jvx::SMatrix{K,K,T},
                            Jvm::SMatrix{K,N,T}, dmdθ::MMatrix{N,P,T}) where {M,P,T,K,N}
    @inbounds for c in 1:P
        xs = SVector{K,T}(ntuple(k -> S[k, c], Val(K)))
        ms = SVector{N,T}(ntuple(k -> dmdθ[k, c], Val(N)))
        dv = Jvx * xs + Jvm * ms
        for k in 1:K
            S[K+k, c] += dv[k]
        end
    end
    return nothing
end

@inline function chain_drift!(S::MMatrix{M,P,T}, h::T, ::Val{K}) where {M,P,T,K}
    @inbounds for c in 1:P, k in 1:K
        S[k, c] += h * S[K+k, c]
    end
    return nothing
end

# The full step, structurally identical to `ahl21_step`:
#   drift(h/2) ∘ pairs(forward) ∘ ϕ†(h) ∘ pairs(reverse) ∘ drift(h/2)
function chain_step!(S::MMatrix{M,P,T}, dmdθ::MMatrix{N,P,T}, b, h::T) where {M,P,T,N}
    K = 3N
    chain_drift!(S, 0.5h, Val(K))
    for p in eachindex(b.idx)
        chain_pair!(S, b.idx[p], b.J[p], b.Jm[p], dmdθ, b.i[p], b.j[p])
    end
    chain_phi!(S, b.Jvx, b.Jvm, dmdθ)
    for p in reverse(eachindex(b.idx))
        chain_pair!(S, b.idx[p], b.J[p], b.Jm[p], dmdθ, b.i[p], b.j[p])
    end
    chain_drift!(S, 0.5h, Val(K))
    return S
end

# Variant C: parameters in the *outer* loop, so one sensitivity column lives in
# registers for the whole step and the only memory traffic is streaming the
# Jacobian blocks (which are shared across columns). The gather/scatter of the
# matrix-at-a-time version above is the obvious thing to blame for its cost, so
# this is the fair upper bound on what any chain-rule layout can do.
@inline function chain_column(s::MVector{M,T}, m::SVector{N,T}, b, h::T, ::Val{K}) where {M,T,N,K}
    @inbounds begin
        for k in 1:K
            s[k] += 0.5h * s[K+k]
        end
        for p in eachindex(b.idx)
            idx = b.idx[p]
            t = SVector{12,T}(ntuple(k -> s[idx[k]], Val(12)))
            u = b.J[p] * t + b.Jm[p] * SVector{2,T}(m[b.i[p]], m[b.j[p]])
            for k in 1:12
                s[idx[k]] = u[k]
            end
        end
        xs = SVector{K,T}(ntuple(k -> s[k], Val(K)))
        dv = b.Jvx * xs + b.Jvm * m
        for k in 1:K
            s[K+k] += dv[k]
        end
        for p in reverse(eachindex(b.idx))
            idx = b.idx[p]
            t = SVector{12,T}(ntuple(k -> s[idx[k]], Val(12)))
            u = b.J[p] * t + b.Jm[p] * SVector{2,T}(m[b.i[p]], m[b.j[p]])
            for k in 1:12
                s[idx[k]] = u[k]
            end
        end
        for k in 1:K
            s[k] += 0.5h * s[K+k]
        end
    end
    return s
end

function chain_step_colwise!(S::MMatrix{M,P,T}, dmdθ::MMatrix{N,P,T}, b, h::T) where {M,P,T,N}
    K = M ÷ 2
    @inbounds for c in 1:P
        s = MVector{M,T}(ntuple(k -> S[k, c], Val(M)))
        m = SVector{N,T}(ntuple(k -> dmdθ[k, c], Val(N)))
        chain_column(s, m, b, h, Val(K))
        for k in 1:M
            S[k, c] = s[k]
        end
    end
    return S
end

# Variant R: exploit the pair block's rank. Δ(x_i,v_i,x_j,v_j) depends on the
# pair only through the *relative* state (x_i−x_j, v_i−v_j), and is written
# back with ±mass fractions — so the 12×12 block upstream materializes and
# BLAS-multiplies actually has rank 6. Applying it as
#   relative (6 subs) → 6×6 → distribute (12 madds)
# is ~78 madds per column against 168 for the dense 12×12. Upstream's own
# `jac_ij` is the dense form, so modelling only that would have been unfair to
# the analytic route.
@inline function chain_column_rank6(s::MVector{M,T}, m::SVector{N,T}, b, h::T,
                                    ::Val{K}) where {M,T,N,K}
    @inbounds begin
        for k in 1:K
            s[k] += 0.5h * s[K+k]
        end
        for pass in 1:2
            rng = pass == 1 ? eachindex(b.idx) : reverse(eachindex(b.idx))
            for p in rng
                idx = b.idx[p]
                rel = SVector{6,T}(ntuple(k -> s[idx[k]] - s[idx[k+6]], Val(6)))
                d = b.J6[p] * rel + b.Jm6[p] * SVector{2,T}(m[b.i[p]], m[b.j[p]])
                fi = b.fi[p]; fj = b.fj[p]
                for k in 1:6
                    s[idx[k]]   += fj * d[k]
                    s[idx[k+6]] -= fi * d[k]
                end
            end
        end
        xs = SVector{K,T}(ntuple(k -> s[k], Val(K)))
        dv = b.Jvx * xs + b.Jvm * m
        for k in 1:K
            s[K+k] += dv[k]
        end
        for k in 1:K
            s[k] += 0.5h * s[K+k]
        end
    end
    return s
end

function chain_step_rank6!(S::MMatrix{M,P,T}, dmdθ::MMatrix{N,P,T}, b, h::T) where {M,P,T,N}
    K = M ÷ 2
    @inbounds for c in 1:P
        s = MVector{M,T}(ntuple(k -> S[k, c], Val(M)))
        m = SVector{N,T}(ntuple(k -> dmdθ[k, c], Val(N)))
        chain_column_rank6(s, m, b, h, Val(K))
        for k in 1:M
            S[k, c] = s[k]
        end
    end
    return S
end

pairs_of(N) = Tuple((i, j) for i in 1:N-1 for j in i+1:N)

function make_blocks(N::Int, ::Type{T}) where {T}
    K = 3N
    prs = pairs_of(N)
    idx = map(((i, j),) -> SVector{12,Int}(vcat(3(i-1)+1:3i, 3(j-1)+1:3j,
                                                K+3(i-1)+1:K+3i, K+3(j-1)+1:K+3j)), prs)
    return (idx = idx,
            i   = map(first, prs),
            j   = map(last, prs),
            J   = map(_ -> SMatrix{12,12,T}(I + 1e-3 .* randn(12, 12)), prs),
            Jm  = map(_ -> SMatrix{12,2,T}(1e-3 .* randn(12, 2)), prs),
            J6  = map(_ -> SMatrix{6,6,T}(1e-3 .* randn(6, 6)), prs),
            Jm6 = map(_ -> SMatrix{6,2,T}(1e-3 .* randn(6, 2)), prs),
            fi  = map(_ -> T(0.3), prs),
            fj  = map(_ -> T(0.7), prs),
            Jvx = SMatrix{K,K,T}(1e-3 .* randn(K, K)),
            Jvm = SMatrix{K,N,T}(1e-3 .* randn(K, N)))
end

# ---------------------------------------------------------------- ForwardDiff
# Cost of one `ahl21_step` with P partials attached to the state and masses.
function dual_step_cost(N::Int, P::Int, w)
    Tg = ForwardDiff.Tag{typeof(identity),Float64}
    D = Dual{Tg,Float64,P}
    seed(x) = D(x, Partials(ntuple(_ -> randn(), P)))
    stD = NBodyState(SMatrix{3,N,D}(seed.(w.st.x)), SMatrix{3,N,D}(seed.(w.st.v)))
    GmD = SVector{N,D}(seed.(w.Gm))
    b = @benchmark ahl21_step($stD, $GmD, $(w.h)) samples=400 evals=3
    return minimum(b).time * 1e-3   # µs
end

function main()
    println("AHL21 analytic-Jacobian budget — is the chain rule alone affordable?")
    println("All times µs per `ahl21_step`-equivalent. 1-core sandbox, minimum of 400 samples.\n")
    @printf("%3s %4s | %8s %8s %8s | %8s %8s %8s %7s | %s\n",
            "N", "P", "primal", "Dual{P}", "FD marg", "chain-M", "chain-C", "chain-R",
            "best/FD", "budget for assembly")
    println("-"^112)
    for N in NS
        w = demo_state(N)
        bp = @benchmark ahl21_step($(w.st), $(w.Gm), $(w.h)) samples=400 evals=3
        tprimal = minimum(bp).time * 1e-3
        blocks = make_blocks(N, Float64)
        for P in PS
            tdual = dual_step_cost(N, P, w)
            marg = tdual - tprimal
            S = MMatrix{6N,P,Float64}(randn(6N, P))
            dm = MMatrix{N,P,Float64}(randn(N, P))
            bm = @benchmark chain_step!($S, $dm, $blocks, $(w.h)) samples=400 evals=3
            bc = @benchmark chain_step_colwise!($S, $dm, $blocks, $(w.h)) samples=400 evals=3
            br = @benchmark chain_step_rank6!($S, $dm, $blocks, $(w.h)) samples=400 evals=3
            tm = minimum(bm).time * 1e-3
            tc = minimum(bc).time * 1e-3
            tr = minimum(br).time * 1e-3
            best = min(tm, tc, tr)
            budget = marg - best      # what the analytic assembly may cost, µs/step
            @printf("%3d %4d | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f %7.2f | %+7.3f µs = %+6.1f%% of primal\n",
                    N, P, tprimal, tdual, marg, tm, tc, tr, best / marg,
                    budget, 100 * budget / tprimal)
        end
        println()
    end
    println("`budget` is what one step's worth of analytic local Jacobians (all pairs +")
    println("corrector) may cost for the analytic path to break even with ForwardDiff.")
    println("Negative means it loses before any assembly is written at all.")
end

main()
