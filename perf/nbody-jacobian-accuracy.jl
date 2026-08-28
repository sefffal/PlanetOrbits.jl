# Does the ForwardDiff pair Jacobian suffer the cancellations Agol's extra
# H-functions were introduced to avoid?
#
# Run:
#   julia --project=perf perf/nbody-jacobian-accuracy.jl
#
# Context. §10.4.3 of the design doc examined Zach Langford's cancellation
# warning at the level of the *individual* special functions (∂G₃/∂γ and
# friends) and found it did not bite. But upstream's analytic Jacobian
# introduces five special functions that appear **nowhere in the primal**:
#
#     H₃ = G₁G₂ − 3G₃            H₆ = 2G₂² − 3G₁G₃
#     H₅ = G₁G₂ − (2+G₀)G₃       H₇ = G₁G₂(1−2G₀) + 3G₀²G₃
#     H₈ = G₁G₂ − 3G₀G₃
#
# Every one is a difference of same-order terms — at small γ the leading
# powers cancel and the result is two orders in γ smaller than its operands.
# They exist *only* in `compute_jacobian_gamma!`. So the cancellation question
# was never really about ∂G₃/∂γ; it is about whether the composite Jacobian
# entries lose digits, and that is what this measures.
#
# Two referees, doing different jobs:
#   (a) ALGEBRA — Float64 AD vs BigFloat(256) central differences of the
#       primal. Catches a wrong rule. Good to ~1e-20 with the step used here.
#   (b) ROUNDOFF — Float64 AD vs BigFloat(256) AD of the *same expressions*.
#       Identical algebra, 256-bit arithmetic, so any disagreement is purely
#       Float64 conditioning: exactly the cancellation question.
#
# Errors are reported against the Jacobian's own norm, never element-by-element
# against a near-zero entry (§10.4.2's lesson).

using StaticArrays
using ForwardDiff
using Printf
using PlanetOrbits: _delxv_gamma, _solve_gamma, G_au3_day2

setprecision(BigFloat, 256)

# The pair kernel as a flat R⁸ → R⁶ map: (x0, v0, gm, h) ↦ (Δx, Δv).
function pair_kernel(u::AbstractVector{T}, drift_first::Bool) where {T}
    x0 = SVector{3,T}(u[1], u[2], u[3])
    v0 = SVector{3,T}(u[4], u[5], u[6])
    delx, delv = _delxv_gamma(x0, v0, u[7], u[8], drift_first)
    return SVector{6,T}(delx[1], delx[2], delx[3], delv[1], delv[2], delv[3])
end

# Central differences in BigFloat. The step is 1e-18 of the *group* scale, not
# of the component — a component that is exactly zero (and the sweep has
# several) would otherwise give δ = 0 and a NaN column.
function fd_jacobian(u::Vector{BigFloat}, drift_first::Bool)
    sx = maximum(abs, @view u[1:3]); sv = maximum(abs, @view u[4:6])
    scale = (sx, sx, sx, sv, sv, sv, abs(u[7]), abs(u[8]))
    J = zeros(BigFloat, 6, 8)
    for c in 1:8
        δ = scale[c] * big"1e-18"
        up = copy(u); up[c] += δ
        um = copy(u); um[c] -= δ
        J[:, c] .= (pair_kernel(up, drift_first) .- pair_kernel(um, drift_first)) ./ (2δ)
    end
    return J
end

# A pair state on the conic selected by `ecc` (< 1 elliptic, > 1 hyperbolic),
# with h solved so the universal anomaly over the step actually *is* γ_target.
# Bisecting rather than using the small-γ estimate matters: γ ≈ sqb·h/r is only
# good below ~0.3, and the interesting band is where the series/closed-form
# switch and the H-function cancellations live.
function pair_geometry(ecc::Float64)
    gm = G_au3_day2 * 1.0
    a = 0.1
    r = a * abs(1 - ecc)                  # periapsis distance, either conic
    x0 = SVector(r, 0.0, 0.0)
    vmag = sqrt(gm * (1 + ecc) / r)       # vis-viva at periapsis, both conics
    v0 = SVector(0.0, vmag * 0.9, vmag * 0.1)
    return (x0 = x0, v0 = v0, gm = gm, r = r)
end

function pair_case(γ_target::Float64, ecc::Float64, drift_first::Bool)
    g = pair_geometry(ecc)
    γ(h) = achieved_gamma((x0 = g.x0, v0 = g.v0, gm = g.gm, h = h), drift_first)
    beta0 = 2g.gm / g.r - sum(abs2, g.v0)
    hi = γ_target * g.r / sqrt(abs(beta0))   # small-γ estimate as a starting bracket
    lo = 0.0
    for _ in 1:200                           # grow the bracket, then bisect
        γ(hi) >= γ_target && break
        lo = hi; hi *= 1.5
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        γ(mid) < γ_target ? (lo = mid) : (hi = mid)
    end
    return (x0 = g.x0, v0 = g.v0, gm = g.gm, h = 0.5 * (lo + hi))
end

function achieved_gamma(c, drift_first)
    rtmp = c.x0 - c.h * c.v0
    r0 = drift_first ? sqrt(sum(abs2, rtmp)) : sqrt(sum(abs2, c.x0))
    beta0 = 2c.gm / r0 - sum(abs2, c.v0)
    signb = sign(beta0); sqb = sqrt(signb * beta0)
    zeta = c.gm - r0 * beta0
    eta = drift_first ? sum(rtmp .* c.v0) : sum(c.x0 .* c.v0)
    return _solve_gamma(c.gm, r0, inv(r0), beta0, signb, sqb, zeta, eta, c.h)
end

function run(drift_first::Bool)
    println(drift_first ? "\n### drift_first = true  (−drift then Kepler)" :
                          "\n### drift_first = false (Kepler then −drift)")
    @printf("%9s %6s %8s | %10s %10s | %10s %10s\n",
            "γ target", "e", "γ actual", "algebra", "worst col", "roundoff", "worst col")
    println("-"^78)
    for ecc in (0.0, 0.6, 0.95, 1.4)
        for γt in (1e-3, 1e-2, 0.1, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0)
            c = pair_case(γt, ecc, drift_first)
            u64 = [c.x0..., c.v0..., c.gm, c.h]
            ubig = BigFloat.(u64)
            J64 = ForwardDiff.jacobian(u -> pair_kernel(u, drift_first), u64)
            Jbig_ad = ForwardDiff.jacobian(u -> pair_kernel(u, drift_first), ubig)
            Jbig_fd = fd_jacobian(ubig, drift_first)
            scale = maximum(abs, Jbig_ad)
            ealg = Float64.(abs.(J64 .- Jbig_fd) ./ scale)
            ernd = Float64.(abs.(J64 .- Jbig_ad) ./ scale)
            γa = achieved_gamma(c, drift_first)
            @printf("%9.0e %6.2f %8.4f | %10.2e %10d | %10.2e %10d\n",
                    γt, ecc, γa, maximum(ealg), argmax(ealg)[2],
                    maximum(ernd), argmax(ernd)[2])
        end
        println()
    end
end

println("Pair-kernel Jacobian ∂(Δx,Δv)/∂(x0,v0,gm,h), errors relative to max|J|.")
println("Columns: 1-3 x0, 4-6 v0, 7 gm, 8 h.")
run(true)
run(false)
