# Does the analytic-Jacobian analysis change for Enzyme instead of ForwardDiff?
#
# Run one mode per process (a hard crash must not take the others with it):
#   for M in fwd rev fwd_rta rev_rta; do
#     julia --project=<env-with-Enzyme> perf/nbody-enzyme.jl $M
#   done
# or just `perf/run-enzyme-probe.sh`.
#
# Three questions with different answers:
#
#  1. Does Enzyme run the AHL21 kernels, and where does it stop? Bisected:
#     `_delxv_gamma` is a pure SVector→SVector function, while `ahl21_step`
#     builds **MMatrix scratch inside the differentiated region**. That is
#     precisely the shape §10.4.2 bisected to as the cause of the last
#     `EnzymeRuntimeActivityError` (an in-region allocation defeating static
#     activity analysis), so it is the prime suspect here too. Running both
#     tells us whether the AHL21 blocker is the kernels' *mathematics* or just
#     their mutable scratch — and the latter is fixable, since `ahl21_step` is
#     already documented as a pure state → state map.
#
#  2. Do our analytic rules apply? No. `_solve_gamma`'s implicit-function rule
#     and the G₃/H₁/H₂ rules are ForwardDiff `Dual` *method dispatch*
#     (§10.4.4); Enzyme never sees them and instead differentiates the Newton
#     iteration and the G/H series, whose termination is exact float equality
#     on the primal (`g3 == g31`) with nothing checking the derivative has
#     converged. The accuracy column is where that shows up.
#
#  3. Does the architecture conclusion change? For **reverse** mode, yes —
#     its cost is independent of the parameter count, which is the analytic
#     Jacobian's only advantage, so it subsumes rather than competes. For
#     **batched forward** the analytic comparison is unchanged in kind, since
#     batched forward still scales with the parameter count exactly as
#     ForwardDiff does.
#
# `set_runtime_activity` is a release valve, not a configuration flag — it
# disables static activity analysis and costs ~2× (§10.4.2). Reported
# separately and never mixed into the headline numbers.

using StaticArrays
using ForwardDiff
using BenchmarkTools
using Printf
using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2, _delxv_gamma
using Enzyme

setprecision(BigFloat, 256)

include(joinpath(@__DIR__, "nbody-ref-step.jl"))

const MODE = isempty(ARGS) ? "fwd" : ARGS[1]
const NSTEP = 5
const MASSES = (1.071, 1.36e-5, 2.34e-5)
const SMA    = (0.0, 0.1153, 0.1283)
const PHASE  = (0.0, 0.4, 2.1)
const H = 2π * sqrt(SMA[2]^3 / (G_au3_day2 * MASSES[1])) / 20

function demo_state(::Type{T}, m::SVector{3,T}) where {T}
    N = 3
    x = zeros(T, 3, N); v = zeros(T, 3, N)
    for i in 2:N
        a = T(SMA[i]); vc = sqrt(T(G_au3_day2)*m[1]/a); ph = T(PHASE[i])
        x[1,i] = a*cos(ph); x[2,i] = a*sin(ph); x[3,i] = a*T(0.002)
        v[1,i] = -vc*sin(ph); v[2,i] = vc*cos(ph)
    end
    for k in 1:3
        x[k,1] = -sum(m[i]*x[k,i] for i in 2:N)/m[1]
        v[k,1] = -sum(m[i]*v[k,i] for i in 2:N)/m[1]
    end
    return NBodyState(SMatrix{3,N,T}(x), SMatrix{3,N,T}(v))
end

# (A) Pure SVector kernel — no mutable scratch anywhere.
function pair_scalar(u::SVector{7,T}) where {T}
    x0 = SVector{3,T}(u[1], u[2], u[3])
    v0 = SVector{3,T}(u[4], u[5], u[6])
    dx, dv = _delxv_gamma(x0, v0, u[7], T(H), true)
    return dx[1] + 2dx[2] + 3dx[3] + 4dv[1] + 5dv[2] + 6dv[3]
end

# (B) Full step — allocates MMatrix scratch inside the region.
function step_scalar(m::SVector{3,T}) where {T}
    st = demo_state(T, m)
    Gm = SVector{3,T}(T(G_au3_day2) .* m)
    for _ in 1:NSTEP
        st = ahl21_step(st, Gm, T(H))
    end
    return st.x[1,2] + 2*st.x[2,3] + 3*st.v[3,1]
end

# Same map over plain Matrix, so a BigFloat reference is reachable (MMatrix
# refuses non-isbits eltypes). Asserted equal to `step_scalar` in Float64.
function step_scalar_ref(m::AbstractVector{T}) where {T}
    st = demo_state(T, SVector{3,T}(m))
    Gm = collect(T(G_au3_day2) .* m)
    x = Matrix{T}(st.x); v = Matrix{T}(st.v)
    xerr = zeros(T, 3, 3); verr = zeros(T, 3, 3)
    for _ in 1:NSTEP
        ref_step(x, v, xerr, verr, Gm, T(H))
    end
    return (x[1,2]+xerr[1,2]) + 2*(x[2,3]+xerr[2,3]) + 3*(v[3,1]+verr[3,1])
end

# (C) The same integration over the propagator's *full* parameter set —
# 18 state components plus 3 masses. This is the case that actually
# discriminates the modes: (B) has 3 parameters, where forward mode is
# structurally favoured and reverse mode's P-independence cannot show.
function state_scalar(u::AbstractVector{T}) where {T}
    x = SMatrix{3,3,T}(ntuple(l -> u[l], Val(9)))
    v = SMatrix{3,3,T}(ntuple(l -> u[9+l], Val(9)))
    Gm = SVector{3,T}(T(G_au3_day2) * u[19], T(G_au3_day2) * u[20],
                      T(G_au3_day2) * u[21])
    st = NBodyState(x, v)
    for _ in 1:NSTEP
        st = ahl21_step(st, Gm, T(H))
    end
    return st.x[1,2] + 2*st.x[2,3] + 3*st.v[3,1]
end

function state_scalar_ref(u::AbstractVector{T}) where {T}
    x = Matrix{T}(reshape(u[1:9], 3, 3)); v = Matrix{T}(reshape(u[10:18], 3, 3))
    Gm = T.(T(G_au3_day2) .* u[19:21])
    xerr = zeros(T, 3, 3); verr = zeros(T, 3, 3)
    for _ in 1:NSTEP
        ref_step(x, v, xerr, verr, Gm, T(H))
    end
    return (x[1,2]+xerr[1,2]) + 2*(x[2,3]+xerr[2,3]) + 3*(v[3,1]+verr[3,1])
end

const U7 = SVector{7,Float64}(0.09, 0.02, 0.001, -0.5, 2.3, 0.05,
                              G_au3_day2 * (MASSES[1] + MASSES[2]))
const M3 = SVector{3,Float64}(MASSES)
const ST0 = demo_state(Float64, M3)
const U21 = SVector{21,Float64}(ST0.x..., ST0.v..., MASSES...)

mode_obj() = MODE in ("rev", "rev_rta") ? Reverse : Forward
with_rta() = endswith(MODE, "_rta")
function mode_arg()
    m = mode_obj()
    return with_rta() ? Enzyme.set_runtime_activity(m) : m
end

function attempt(label, f, ref, scale)
    print(rpad(label, 26), " "); flush(stdout)
    local g
    try
        g = f()
    catch e
        println("FAILED  ", replace(sprint(showerror, e), '\n' => ' ')[1:min(end, 150)])
        return
    end
    gv = collect(g isa Tuple ? g[1] : g)
    err = Float64(maximum(abs, BigFloat.(gv) .- ref) / scale)
    @printf("ok   relerr %9.2e", err)
    b = @benchmark $f() samples=100 evals=1
    @printf("   %9.2f µs\n", minimum(b).time * 1e-3)
end

function main()
    @printf("mode = %s   (runtime activity: %s)\n\n", MODE, with_rta() ? "ON — ~2× valve" : "off")

    rp = ForwardDiff.gradient(pair_scalar, SVector{7,BigFloat}(big.(U7)))
    rs = ForwardDiff.gradient(step_scalar_ref, SVector{3,BigFloat}(big.(MASSES)))
    sp = maximum(abs, rp); ss = maximum(abs, rs)
    @printf("reference map reproduces ahl21_step in Float64: %s\n\n",
            step_scalar_ref(M3) === step_scalar(M3) ? "bit-identical" :
            "differs by $(abs(step_scalar_ref(M3) - step_scalar(M3)))")

    println("(A) `_delxv_gamma` — pure SVector, no mutable scratch")
    attempt("  ForwardDiff", () -> ForwardDiff.gradient(pair_scalar, U7), rp, sp)
    attempt("  Enzyme $(MODE)", () -> Enzyme.gradient(mode_arg(), pair_scalar, U7), rp, sp)

    println("\n(B) `ahl21_step` ×$NSTEP, P=3 (masses) — MMatrix scratch in-region")
    attempt("  ForwardDiff", () -> ForwardDiff.gradient(step_scalar, M3), rs, ss)
    attempt("  Enzyme $(MODE)", () -> Enzyme.gradient(mode_arg(), step_scalar, M3), rs, ss)

    ru = ForwardDiff.gradient(state_scalar_ref, SVector{21,BigFloat}(big.(U21)))
    su = maximum(abs, ru)
    println("\n(C) `ahl21_step` ×$NSTEP, P=21 (full state + masses) — the case that")
    println("    discriminates the modes; reverse mode should be flat in P here")
    attempt("  ForwardDiff", () -> ForwardDiff.gradient(state_scalar, U21), ru, su)
    attempt("  Enzyme $(MODE)", () -> Enzyme.gradient(mode_arg(), state_scalar, U21), ru, su)
end

main()
