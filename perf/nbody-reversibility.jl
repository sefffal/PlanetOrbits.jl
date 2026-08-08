# Is a tape-free reverse mode admissible for AHL21?
#
#   julia --project=perf perf/nbody-reversibility.jl
#
# Reverse-mode AD over an n-step integration needs the intermediate states
# q₀…q_{n-1} during the reverse sweep, because step k's local Jacobian is
# evaluated at q_k. Three ways to get them:
#
#   (1) tape      — store them (and every scalar inside each step). O(n) memory,
#                   and the cache traffic is what makes Enzyme reverse 1.5×
#                   slower than ForwardDiff at P=21 (§10.6.6).
#   (2) checkpoint — store every √n, recompute *forward* inside a window.
#                   O(√n) memory, ~2× primal time, and **exact**: forward
#                   recomputation reproduces the primal bit-for-bit.
#   (3) reverse the map — store only q_n and recover q_{k-1} from q_k by
#                   stepping the velocity-reversed state. O(1) memory, no
#                   recomputation. This is the "tape-free" option, and it is
#                   admissible only if the recovered q_k is close enough to the
#                   real one that the adjoint is still the adjoint of the
#                   problem we actually solved.
#
# AHL21 is a palindromic composition, so (3) is exact in exact arithmetic and
# `_march!` already relies on it to reach epochs before t0. In floating point
# it is not, for a specific reason: the Kepler drift solves γ by Newton with
# exact-float-equality termination, so the reversed solve is a *different*
# root-find, not a sign flip of the forward one.
#
# The measurement below is a round trip — n steps forward, reverse v, n steps
# forward, reverse v — against two controls that separate the two ways it can
# go wrong:
#
#   * 2-body: the map is exact for one pair (φ† vanishes), so the system is
#     integrable and non-chaotic. Whatever error appears here is *reversal
#     roundoff* alone.
#   * 1-ulp probe: perturb the initial state in the last bit and march forward
#     n steps. This is the intrinsic amplification of a roundoff-sized
#     perturbation. If the round trip tracks it, reversal is as good as the
#     problem allows and the limit is chaos, which bounds every method equally
#     — including the taped one, whose *gradient* is then just as ill-
#     conditioned. If the round trip is much worse, reversal itself is lossy.

using StaticArrays
using Printf
using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2, _reverse_v

# Kepler-36-like: two close, near-resonant super-Earths — deliberately the
# chaotic end of what this propagator is used for.
const M3   = (1.071, 1.36e-5, 2.34e-5)
const SMA3 = (0.0, 0.1153, 0.1283)
const PH3  = (0.0, 0.4, 2.1)

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

march(st, Gm, h, n) = (for _ in 1:n; st = ahl21_step(st, Gm, h); end; st)

# n steps forward, then n steps of the velocity-reversed state, then undo the
# velocity flip. Exact arithmetic would return `st` unchanged.
roundtrip(st, Gm, h, n) = _reverse_v(march(_reverse_v(march(st, Gm, h, n)), Gm, h, n))

# Relative separation of two states, measured against the trajectory's own
# scale (never element-by-element against a near-zero component).
function sep(a::NBodyState, b::NBodyState)
    sx = maximum(abs, a.x); sv = maximum(abs, a.v)
    return max(maximum(abs, a.x .- b.x) / sx, maximum(abs, a.v .- b.v) / sv)
end

bitexact(a::NBodyState, b::NBodyState) =
    a.x === b.x && a.v === b.v && a.xerr === b.xerr && a.verr === b.verr

# Perturb one position component by 1 ulp — a roundoff-sized kick.
function ulp_kick(st::NBodyState{N,T,L}) where {N,T,L}
    x = MVector(st.x[:]); x[4] = nextfloat(x[4])
    return NBodyState{N,T,L}(SMatrix{3,N,T,L}(x), st.v, st.xerr, st.verr)
end

function report(label, m, sma, ph, h, ns)
    Gm = SVector{length(m),Float64}(G_au3_day2 .* m)
    st0 = demo_state(Float64, m, sma, ph)
    println("\n", label)
    @printf("  %6s  %12s  %12s  %10s\n", "steps", "round trip", "1-ulp kick", "bit-exact")
    for n in ns
        rt = roundtrip(st0, Gm, h, n)
        kick = sep(march(st0, Gm, h, n), march(ulp_kick(st0), Gm, h, n))
        @printf("  %6d  %12.2e  %12.2e  %10s\n", n, sep(st0, rt), kick,
                bitexact(st0, rt) ? "yes" : "no")
    end
end

function main()
    P_in = 2π * sqrt(SMA3[2]^3 / (G_au3_day2 * M3[1]))
    h = P_in / 20
    ns = (1, 10, 100, 1_000, 10_000)

    @printf("inner period %.4f d, h = P/20 = %.5f d\n", P_in, h)

    report("3-body (Kepler-36-like, near-resonant — chaotic)", M3, SMA3, PH3, h, ns)
    report("2-body (map is exact for one pair — integrable control)",
           (M3[1], M3[2]), (SMA3[1], SMA3[2]), (PH3[1], PH3[2]), h, ns)
end

main()
