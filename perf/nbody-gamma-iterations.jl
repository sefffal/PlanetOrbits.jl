# How many Newton iterations does the universal-variable γ solve actually take?
#
#   julia --project=perf perf/nbody-gamma-iterations.jl
#
# This is the number that decides whether an `EnzymeRules` rule for
# `_solve_gamma` is worth writing (§10.6.8). Enzyme has no such rule — the
# implicit rule is ForwardDiff `Dual` method dispatch, invisible to it — so
# Enzyme differentiates the loop while ForwardDiff short-circuits it to one
# correction. With k iterations and W partials, the available multiplier on the
# γ portion is:
#
#   batched forward:  k·W  →  k + W    (tangents through every iteration)
#   reverse:          ~3k  →  k + 3    (per-iteration cost is P-independent)
#
# so the arithmetic win is larger for forward, and both scale with k. At k = 3
# the rule is barely worth writing; at k = 8–10 it plausibly reorders the modes.
#
# Method: the counting kernels below are line-for-line copies of
# `_kepler_driftij!` and `_gamma_iterate`, but they call the *real*
# `_delxv_gamma` for the state update, so the trajectory is the real one and
# only the iteration count is added. `main` asserts the stepped state is `===`
# `ahl21_step`'s, which is what makes the counts trustworthy: if the copies ever
# drift from the originals, the gate fails rather than reporting a wrong k.

using StaticArrays
using Printf
using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2,
                    _delxv_gamma, _gamma_F_dF, _cubic1, _dot3, _comp_sum,
                    _drift!, _phisalpha!

const M3   = (1.071, 1.36e-5, 2.34e-5)     # Kepler-36-like
const SMA3 = (0.0, 0.1153, 0.1283)
const PH3  = (0.0, 0.4, 2.1)

# Exact copy of `_gamma_iterate`'s initial guess + Newton loop, instrumented.
function gamma_iterations(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
    if !iszero(zeta)
        zinv = 6 / zeta
        gamma = _cubic1(0.5 * eta * sqb * zinv, r0 * signb * beta0 * zinv,
                        -h * signb * beta0 * sqb * zinv)
    elseif !iszero(eta)
        reta = r0 / eta
        disc = reta^2 + 2h / eta
        gamma = disc > 0 ? sqb * (-reta + sqrt(disc)) : h * r0inv * sqb
    else
        gamma = h * r0inv * sqb
    end
    gamma1 = 2gamma
    gamma2 = 3gamma
    iter = 0
    step = zero(gamma)
    while iter < 20
        gamma2 = gamma1
        gamma1 = gamma
        F, dF = _gamma_F_dF(gamma, gm, r0, beta0, signb, sqb, zeta, eta, h)
        step = F / dF
        gamma -= step
        iter += 1
        (gamma == gamma2 || gamma == gamma1) && break
    end
    # Relative size of the last Newton correction. At a converged root this is
    # ~1 ulp; if the cap is hit because the iterate is in a >2-cycle at the
    # ulp level it stays ~1e-16, and if the solve genuinely failed it does not.
    return iter, gamma, abs(step) / max(abs(gamma), eps(gamma))
end

# Exact copy of `_delxv_gamma`'s prologue, up to the γ solve.
function count_pair!(counts::Vector{Pair{Int,Float64}}, x0::SVector{3}, v0::SVector{3}, gm, h, drift_first::Bool)
    rtmp = x0 - h * v0
    r0 = drift_first ? sqrt(_dot3(rtmp)) : sqrt(_dot3(x0))
    r0inv = inv(r0)
    beta0 = 2 * gm * r0inv - _dot3(v0)
    signb = sign(beta0)
    sqb = sqrt(signb * beta0)
    zeta = gm - r0 * beta0
    eta = drift_first ? _dot3(rtmp, v0) : _dot3(x0, v0)
    k, _, res = gamma_iterations(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
    push!(counts, k => res)
    return nothing
end

# Exact copy of `_kepler_driftij!`, plus the count.
@inline function count_driftij!(counts, x, v, xerr, verr, Gm, i::Int, j::Int, h, drift_first::Bool)
    @inbounds x0 = SVector(x[1,i]-x[1,j], x[2,i]-x[2,j], x[3,i]-x[3,j])
    @inbounds v0 = SVector(v[1,i]-v[1,j], v[2,i]-v[2,j], v[3,i]-v[3,j])
    gm = Gm[i] + Gm[j]
    iszero(gm) && return nothing
    count_pair!(counts, x0, v0, gm, h, drift_first)
    delx, delv = _delxv_gamma(x0, v0, gm, h, drift_first)
    gminv = inv(gm)
    mi = Gm[i] * gminv
    mj = Gm[j] * gminv
    @inbounds for k in 1:3
        x[k,i], xerr[k,i] = _comp_sum(x[k,i], xerr[k,i],  mj*delx[k])
        x[k,j], xerr[k,j] = _comp_sum(x[k,j], xerr[k,j], -mi*delx[k])
        v[k,i], verr[k,i] = _comp_sum(v[k,i], verr[k,i],  mj*delv[k])
        v[k,j], verr[k,j] = _comp_sum(v[k,j], verr[k,j], -mi*delv[k])
    end
    return nothing
end

# Exact copy of `ahl21_step`'s composition.
function count_step(counts, st::NBodyState{N,T,L}, Gm::SVector{N}, h) where {N,T,L}
    x = MMatrix(st.x); xerr = MMatrix(st.xerr)
    v = MMatrix(st.v); verr = MMatrix(st.verr)
    h2 = 0.5 * h
    _drift!(x, xerr, v, h2)
    for i in 1:N-1, j in i+1:N
        count_driftij!(counts, x, v, xerr, verr, Gm, i, j, h2, true)
    end
    _phisalpha!(x, v, verr, Gm, h)
    for i in N-1:-1:1, j in N:-1:i+1
        count_driftij!(counts, x, v, xerr, verr, Gm, i, j, h2, false)
    end
    _drift!(x, xerr, v, h2)
    return NBodyState{N,T,L}(SMatrix(x), SMatrix(v), SMatrix(xerr), SMatrix(verr))
end

function demo_state(::Type{T}, m::NTuple{N,<:Real}, sma, ph) where {T,N}
    x = zeros(T, 3, N); v = zeros(T, 3, N)
    for i in 2:N
        a = T(sma[i]); vc = sqrt(T(G_au3_day2) * T(m[1]) / a); p = T(ph[i])
        x[1,i] = a*cos(p); x[2,i] = a*sin(p); x[3,i] = a*T(0.002)
        v[1,i] = -vc*sin(p); v[2,i] = vc*cos(p)
    end
    for k in 1:3
        x[k,1] = -sum(T(m[i])*x[k,i] for i in 2:N)/T(m[1])
        v[k,1] = -sum(T(m[i])*v[k,i] for i in 2:N)/T(m[1])
    end
    return NBodyState(SMatrix{3,N,T}(x), SMatrix{3,N,T}(v))
end

bitexact(a::NBodyState, b::NBodyState) =
    a.x === b.x && a.v === b.v && a.xerr === b.xerr && a.verr === b.verr

function run(label, m, sma, ph, h, nstep)
    Gm = SVector{length(m),Float64}(G_au3_day2 .* m)
    st = demo_state(Float64, m, sma, ph)
    ck = st                       # counted march
    ref = st                      # real march
    counts = Pair{Int,Float64}[]
    for _ in 1:nstep
        ck = count_step(counts, ck, Gm, h)
        ref = ahl21_step(ref, Gm, h)
    end
    gate = bitexact(ck, ref) ? "bit-identical" : "*** DIVERGED — counts untrustworthy ***"

    hist = zeros(Int, 21)
    for (k, _) in counts; hist[k+1] += 1; end
    kbar = sum(first, counts) / length(counts)
    capped = [r for (k, r) in counts if k == 20]
    worst = maximum(last, counts)
    @printf("\n%s\n  gate vs ahl21_step: %s   (%d γ solves over %d steps)\n",
            label, gate, length(counts), nstep)
    print("  k histogram: ")
    for k in 0:20
        hist[k+1] > 0 && @printf("k=%d:%.1f%%  ", k, 100*hist[k+1]/length(counts))
    end
    @printf("\n  mean k = %.2f, max = %d\n", kbar, maximum(first, counts))
    @printf("  final |Δγ|/γ: worst over all solves %.2e", worst)
    isempty(capped) ? println() :
        @printf("; over the %d capped solves %.2e … %.2e\n",
                length(capped), minimum(capped), maximum(capped))
    for W in (7, 21)
        @printf("    W=%2d  forward multiplier kW/(k+W) = %.2f×   reverse 3k/(k+3) = %.2f×\n",
                W, kbar*W/(kbar+W), 3kbar/(kbar+3))
    end
    return kbar
end

function main()
    P_in = 2π * sqrt(SMA3[2]^3 / (G_au3_day2 * M3[1]))
    @printf("inner period %.4f d\n", P_in)
    for div in (20, 60, 200)
        run(@sprintf("3-body Kepler-36-like, h = P/%d = %.5f d", div, P_in/div),
            M3, SMA3, PH3, P_in/div, 200)
    end
    # An eccentric, wider pair: γ per step is larger, so this is the
    # unfavourable end for the initial guess.
    run("3-body, wider + eccentric outer (h = P_in/20)",
        (1.0, 1e-3, 5e-4), (0.0, 0.1, 0.9), (0.0, 0.0, 2.8), P_in/20, 200)
end

main()
