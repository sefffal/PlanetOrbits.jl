# The one cancellation upstream engineered around by hand — does it bite us?
#
# Run:
#   julia --project=perf perf/nbody-mass-derivative.jl
#
# In `kepler_driftij_gamma!` NbodyGradient computes the mass column of the pair
# Jacobian two different ways, and says why:
#
#     jac_ij[k, 7] = jac_mass[k]*m[j]                                  # ∂/∂m_i
#     jac_ij[k,14] = mi*delxv[k]*mijinv + GNEWT*mj*jac_kepler[k,7]     # ∂/∂m_j
#     # "Compute the mass jacobian separately since otherwise cancellations
#     #  happen in kepler_driftij_gamma"
#
# The pair kernel applies  x_i += mj·Δx(gm)  with  mj = m_j/(m_i+m_j),
# gm = G(m_i+m_j). So
#
#     ∂(mj·Δx)/∂m_i = −m_j/(m_i+m_j)² · Δx  +  mj·G·∂Δx/∂gm
#
# — two terms of opposite sign. To leading order Δx is linear in gm, which
# makes the second term ≈ +m_j/(m_i+m_j)²·Δx: the leading order cancels
# identically. That is why upstream derived a closed form for it. The ∂/∂m_j
# direction has no such cancellation (both terms are positive), which is
# exactly why upstream leaves *that* one naive.
#
# Our ForwardDiff path forms the cancelling combination in Float64 for the
# ∂/∂m_i direction. This measures how many digits that costs, against a
# BigFloat(256) evaluation of the identical expressions.
#
# Two ratios, doing different jobs:
#   cancel  = |operand| / |result|  — digits *lost* to cancellation
#   relerr  = |Float64 − BigFloat| / |result|  — digits *remaining*
# A large `cancel` with a small `relerr` means the cancellation is benign.

using StaticArrays
using ForwardDiff
using Printf
using PlanetOrbits: _delxv_gamma, _solve_gamma, G_au3_day2

setprecision(BigFloat, 256)

# The increment the pair kernel adds to body i, as a function of m_i alone.
function incr_i(mi_, mj_, x0::SVector{3}, v0::SVector{3}, h, drift_first)
    T = promote_type(typeof(mi_), typeof(mj_), eltype(x0), typeof(h))
    gm = T(G_au3_day2) * (mi_ + mj_)
    delx, delv = _delxv_gamma(SVector{3,T}(x0), SVector{3,T}(v0), gm, T(h), drift_first)
    w = mj_ / (mi_ + mj_)
    return SVector{6,T}(w * delx[1], w * delx[2], w * delx[3],
                        w * delv[1], w * delv[2], w * delv[3])
end

# The two cancelling operands, at the working precision — for the `cancel` ratio.
function operand_scale(mi_::T, mj_::T, x0, v0, h, drift_first) where {T<:Real}
    gm = T(G_au3_day2) * (mi_ + mj_)
    delx, delv = _delxv_gamma(SVector{3,T}(x0), SVector{3,T}(v0), gm, T(h), drift_first)
    s = mj_ / (mi_ + mj_)^2
    return SVector{6,T}(s * delx[1], s * delx[2], s * delx[3],
                        s * delv[1], s * delv[2], s * delv[3])
end

function pair_case(γ_target::Float64, ecc::Float64, mratio::Float64)
    mi = 1.0                       # host
    mj = mratio                    # companion
    gm = G_au3_day2 * (mi + mj)
    a = 0.1
    r = a * (1 - ecc)
    x0 = SVector(r, 0.0, 0.0)
    vmag = sqrt(gm * (1 + ecc) / r)
    v0 = SVector(0.0, vmag * 0.9, vmag * 0.1)
    beta0 = 2gm / r - sum(abs2, v0)
    h = γ_target * r / sqrt(abs(beta0))
    return (mi = mi, mj = mj, x0 = x0, v0 = v0, h = h)
end

function run(drift_first::Bool)
    println(drift_first ? "\n### drift_first = true" : "\n### drift_first = false")
    @printf("%7s %6s %8s | %10s %10s | %10s %10s\n",
            "γ", "e", "m_j/m_i", "∂/∂m_i", "cancel", "relerr", "digits left")
    println("-"^74)
    for mratio in (1e-3, 1e-5, 1e-8)
        for ecc in (0.0, 0.6), γt in (1e-2, 0.1, 0.5, 1.0, 2.0)
            c = pair_case(γt, ecc, mratio)
            f64(m) = incr_i(m, c.mj, c.x0, c.v0, c.h, drift_first)
            fbig(m) = incr_i(m, big(c.mj), BigFloat.(c.x0), BigFloat.(c.v0),
                             big(c.h), drift_first)
            d64 = ForwardDiff.derivative(f64, c.mi)
            dbig = ForwardDiff.derivative(fbig, big(c.mi))
            res = maximum(abs, dbig)
            opd = maximum(abs, operand_scale(big(c.mi), big(c.mj), c.x0, c.v0,
                                             c.h, drift_first))
            cancel = Float64(opd / res)
            relerr = Float64(maximum(abs, BigFloat.(d64) .- dbig) / res)
            @printf("%7.2f %6.2f %8.0e | %10.3e %10.1e | %10.2e %10.1f\n",
                    γt, ecc, mratio, Float64(res), cancel, relerr,
                    relerr > 0 ? -log10(relerr) : 17.0)
        end
        println()
    end
end

println("Mass column of the pair Jacobian: ∂(m_j/(m_i+m_j) · Δx)/∂m_i.")
println("`cancel` = |leading operand| / |result| (digits lost to cancellation).")
println("`relerr` = Float64 ForwardDiff vs BigFloat(256), relative to the result.")
run(true)
run(false)
