# ---------------------------------------------------
# Analytic Jacobian of the AHL21 pair kernel — NOT part of the package.
#
# This file is deliberately **not** `include`d by `PlanetOrbits.jl`. It is the
# analytic derivative machinery of NbodyGradient.jl v0.2.3 (`jac_delxv_gamma!`
# + `compute_jacobian_gamma!`, © 2021 Eric Agol, David Hernandez & Zach
# Langford, MIT), ported to our SMatrix state so its accuracy and its cost can
# be measured against the shipped ForwardDiff path. It is kept for future
# comparison — see §10.6 of design/planetorbits-v2-nbody-migration.md for the
# measurement that decided against shipping it.
#
# Load it with `include`, e.g. from perf/ or test/:
#   include(joinpath(pkgdir(PlanetOrbits), "src", "nbody-jacobian.jl"))
#
# What it computes, for one pair over one (sub)step:
#
#   J = ∂(Δx, Δv) / ∂(x0₁,x0₂,x0₃, v0₁,v0₂,v0₃, k, h)      (6×8)
#   jac_mass = the cancellation-free form of the ∂/∂m_i column (6)
#
# where k = G(m_i + m_j). The 3×3 sub-blocks are a diagonal plus two rank-1
# outer products — which is the structural reason a matrix-free JVP (what
# ForwardDiff does) is cheaper than forming and applying J.
# ---------------------------------------------------

using StaticArrays
using PlanetOrbits: _G3, _H1, _H2, _solve_gamma, _dot3, GH_SERIES_CUTOFF

# ---------------------------------------------------
# The four special functions that appear ONLY in the Jacobian. Each is a
# difference of same-order terms — at small γ the leading powers cancel and
# the result is two orders in γ below its operands, which is exactly why
# upstream gave each one its own series rather than evaluating the difference.
#
#   H₃ = G₁G₂ − 3G₃          H₆ = 2G₂² − 3G₁G₃
#   H₅ = G₁G₂ − (2+G₀)G₃     H₈ = G₁G₂ − 3G₀G₃
#
# Series coefficients transcribed from NbodyGradient `src/utils.jl`.
# ---------------------------------------------------

_JH3(γ, β, sqb; gc=GH_SERIES_CUTOFF) =
    γ < gc  ? _JH3_series(γ, β, sqb) :
    β >= 0  ? (4sin(γ) - sin(γ)*cos(γ) - 3γ) / (β*sqb) :
              (4sinh(γ) - sinh(γ)*cosh(γ) - 3γ) / (β*sqb)

function _JH3_series(γ, β, sqb)
    x2 = -sign(β) * γ^2
    term = one(γ)/30; h3 = one(γ)/10
    h31 = 2h3; h32 = 2h3; n = 0; iter = 0; four2n = 4*one(γ)
    while true
        h32 = h31; h31 = h3; n += 1
        term *= x2; term /= (2n+4)*(2n+5); four2n *= 4
        h3 += term * (four2n - 1)
        iter += 1
        (iter >= 100 || h3 == h32 || h3 == h31) && break
    end
    return h3 * (-x2^2 * γ / (β * sqb))
end

_JH5(γ, β, sqb; gc=GH_SERIES_CUTOFF) =
    γ < gc  ? _JH5_series(γ, β, sqb) :
    β >= 0  ? (3sin(γ) - 2γ - γ*cos(γ)) / (β*sqb) :
              (3sinh(γ) - 2γ - γ*cosh(γ)) / (β*sqb)

function _JH5_series(γ, β, sqb)
    x2 = -sign(β) * γ^2
    term = one(γ)/60; h5 = term
    h51 = 2h5; h52 = 2h5; n = 0; iter = 0
    while true
        h52 = h51; h51 = h5; n += 1
        term *= x2*(n+1); term /= (2n+5)*(2n+4)*n
        h5 += term
        iter += 1
        (iter >= 100 || h5 == h52 || h5 == h51) && break
    end
    return h5 * (-x2^2 * γ / (β * sqb))
end

_JH6(γ, β; gc=GH_SERIES_CUTOFF) =
    γ < gc  ? _JH6_series(γ, β) :
    β >= 0  ? (9 - 8cos(γ) - cos(2γ) - 6γ*sin(γ)) / (2β^2) :
              (9 - 8cosh(γ) - cosh(2γ) + 6γ*sinh(γ)) / (2β^2)

function _JH6_series(γ, β)
    x2 = -sign(β) * γ^2
    term = one(γ)/360; h6 = one(γ)/40
    h61 = 2h6; h62 = 2h6; n = 0; iter = 0; four2n = 16*one(γ)
    while true
        h62 = h61; h61 = h6; n += 1
        term *= x2; term /= (2n+5)*(2n+6); four2n *= 4
        h6 += term * (four2n - 3n - 7)
        iter += 1
        (iter >= 100 || h6 == h62 || h6 == h61) && break
    end
    return h6 * (-x2^3 / β^2)
end

_JH8(γ, β, sqb; gc=GH_SERIES_CUTOFF) =
    γ < gc  ? _JH8_series(γ, β, sqb) :
    β >= 0  ? (-3γ*cos(γ) + sin(γ) + sin(2γ)) / (β*sqb) :
              (-3γ*cosh(γ) + sinh(γ) + sinh(2γ)) / (β*sqb)

function _JH8_series(γ, β, sqb)
    x2 = -sign(β) * γ^2
    term = one(γ)/120; h8 = 3*one(γ)/20
    h81 = 2h8; h82 = 2h8; n = 0; iter = 0; four2n = 32*one(γ)
    while true
        h82 = h81; h81 = h8; n += 1
        term *= x2; term /= (2n+4)*(2n+5); four2n *= 4
        h8 += term * (four2n - 14 - 6n)
        iter += 1
        (iter >= 100 || h8 == h82 || h8 == h81) && break
    end
    return h8 * (x2^2 * γ / (β * sqb))
end

# ---------------------------------------------------
# Primal pass, returning every intermediate the Jacobian needs. Mirrors
# `PlanetOrbits._delxv_gamma` exactly — the returned Δx/Δv must be
# bit-identical to it, which the gate asserts with `===`.
# ---------------------------------------------------
function _delxv_gamma_full(x0::SVector{3}, v0::SVector{3}, k, h, drift_first::Bool)
    rtmp = x0 - h * v0
    r0 = drift_first ? sqrt(_dot3(rtmp)) : sqrt(_dot3(x0))
    r0inv = inv(r0)
    beta = 2 * k * r0inv - _dot3(v0)
    betainv = inv(beta)
    signb = sign(beta)
    sqb = sqrt(signb * beta)
    zeta = k - r0 * beta
    eta = drift_first ? _dot3(rtmp, v0) : _dot3(x0, v0)
    gamma = _solve_gamma(k, r0, r0inv, beta, signb, sqb, zeta, eta, h)
    # Spelled exactly as `_delxv_gamma` spells it — `2*signb*sx^2*betainv` is
    # not bit-identical to `_g12`'s `2*sign(beta)*sx^2/beta`, and the gate
    # asserts the primal here is `===` the shipped one.
    xx = 0.5 * gamma
    if beta > 0
        sx, cx = sincos(xx)
    else
        sx = sinh(xx); cx = exp(-xx) + sx
    end
    g1 = 2 * sx * cx / sqb
    g2 = 2 * signb * sx^2 * betainv
    g0 = 1 - beta * g2
    g3 = _G3(gamma, beta, sqb)
    r = r0 * g0 + eta * g1 + k * g2
    rinv = inv(r)
    dfdt = -k * g1 * rinv * r0inv
    if drift_first
        h1 = zero(gamma); h2 = zero(gamma)
        fm1 = -k * r0inv * g2
        gmh = k * r0inv * (h * g2 - r0 * g3)
        dgdtm1 = k * r0inv * rinv * (h * g1 - r0 * g2)
    else
        h1 = _H1(gamma, beta); h2 = _H2(gamma, beta, sqb)
        fm1 = k * rinv * (g2 - k * r0inv * h1)
        gmh = k * rinv * (r0 * h2 + eta * h1)
        dgdtm1 = -k * rinv * g2
    end
    delx = fm1 * x0 + gmh * v0
    delv = dfdt * x0 + dgdtm1 * v0
    return (; delx, delv, gamma, g0, g1, g2, g3, h1, h2, dfdt, fm1, gmh, dgdtm1,
            r0, r, r0inv, rinv, k, h, beta, betainv, eta, sqb, zeta)
end

# ---------------------------------------------------
# Assembly. Every 3×3 block of J is
#
#     coeff·I + x0 ⊗ (Axx·x0 + Axv·v0) + v0 ⊗ (Bxx·x0 + Bxv·v0)
#
# (upstream's doubly-nested i/j loop, written as what it is). `A` carries the
# derivatives of the x0-coefficient of the output, `B` those of the
# v0-coefficient; `dia` is the coefficient itself, appearing on the diagonal.
# ---------------------------------------------------
struct _CoeffDerivs{T}
    xx::T; xv::T; vx::T; vv::T; k::T; h::T
end

@inline function _jac_rows(x0::SVector{3,T}, v0::SVector{3,T}, dia_x, dia_v,
                           A::_CoeffDerivs, B::_CoeffDerivs) where {T}
    # Returns a 3×8 block: columns 1-3 ∂/∂x0, 4-6 ∂/∂v0, 7 ∂/∂k, 8 ∂/∂h.
    return SMatrix{3,8,T}(ntuple(Val(24)) do lin
        j = (lin - 1) % 3 + 1
        c = (lin - 1) ÷ 3 + 1
        if c <= 3
            i = c
            (i == j ? dia_x : zero(T)) +
                (A.xx*x0[i] + A.xv*v0[i]) * x0[j] + (B.xx*x0[i] + B.xv*v0[i]) * v0[j]
        elseif c <= 6
            i = c - 3
            (i == j ? dia_v : zero(T)) +
                (A.vx*x0[i] + A.vv*v0[i]) * x0[j] + (B.vx*x0[i] + B.vv*v0[i]) * v0[j]
        elseif c == 7
            A.k * x0[j] + B.k * v0[j]
        else
            A.h * x0[j] + B.h * v0[j]
        end
    end)
end

"""
    _delxv_gamma_jac(x0, v0, k, h, drift_first) -> (delx, delv, J, jac_mass)

Analytic Jacobian of the universal-variable Kepler drift for one pair.
`J` is 6×8: `∂(Δx, Δv) / ∂(x0, v0, k, h)`. `jac_mass` is the cancellation-free
spelling of the ∂/∂m_i column described in `nbody-mass-derivative.jl`, scaled
so the caller applies it as `jac_ij[:, m_i] = jac_mass * m_j` (upstream's
convention); it carries a factor G² and so takes `G` as a keyword.
"""
function _delxv_gamma_jac(x0::SVector{3,T}, v0::SVector{3,T}, k, h,
                          drift_first::Bool; G=one(T)) where {T}
    s = _delxv_gamma_full(x0, v0, k, h, drift_first)
    (; gamma, g0, g1, g2, g3, h1, h2, dfdt, fm1, gmh, dgdtm1,
       r0, r, r0inv, rinv, beta, betainv, eta, sqb, zeta) = s
    r0inv2 = r0inv^2; r0inv3 = r0inv2 * r0inv
    rinv2 = rinv^2;   rinv3 = rinv2 * rinv
    hsq = h^2; ksq = k^2

    h6 = _JH6(gamma, beta)
    h3 = _JH3(gamma, beta, sqb)
    h5 = _JH5(gamma, beta, sqb)
    h8 = -2h3 + 3h5

    d  = (h + eta*g2 + 2*k*g3) * betainv
    c1 = d - r0*g3
    c2 = eta*g0 + g1*zeta
    c3 = d*k + g1*r0^2
    c17 = r0 - r - g2*k
    c20 = k*(g2*k + r) - g0*r0*zeta
    c12 = g0*h - g1*r0

    if drift_first
        c13 = g1*h - g2*r0
        c9  = 2*g2*h - 3*g3*r0
        c10 = k*r0inv2^2 * (-g2*r0*h + k*c9*betainv - c3*c13*rinv)
        c24 = r0inv3 * (r0*(2*k*r0inv - beta)*betainv - g1*c3*rinv/g2)

        A = _CoeffDerivs(
            fm1*c24,                                            # dfm1dxx
            -fm1*(g1*rinv + h*c24),                             # dfm1dxv
            -fm1*(g1*rinv + h*c24),                             # dfm1dvx
            fm1*rinv*(-r0*g2 + k*h6*betainv/g2 + h*(2*g1 + h*r*c24)),
            fm1*(1/k + g1*c1*rinv*r0inv/g2 - 2*betainv*r0inv),  # dfm1dk
            fm1*(g1*rinv*(1/g2 + 2*k*r0inv - beta) - eta*c24))  # dfm1dh
        B = _CoeffDerivs(
            c10,                                                # dgmhdxx
            -g2*k*c13*rinv*r0inv - h*c10,                       # dgmhdxv
            -g2*k*c13*rinv*r0inv - h*c10,                       # dgmhdvx
            2*g2*h*k*c13*rinv*r0inv + hsq*c10 +
                k*betainv*rinv*r0inv*(r0^2*h8 - beta*h*r0*g2^2 + (h*k + eta*r0)*h6),
            r0inv*(k*c1*c13*rinv*r0inv + g2*h - g3*r0 - k*c9*betainv*r0inv),
            g2*k*r0inv + k*c13*rinv*r0inv +
                g2*k*(2*k*r0inv - beta)*c13*rinv*r0inv - eta*c10)

        c21 = (g2*k - r0)*(beta*c3 - k*g1*r)*betainv*rinv2*r0inv3/g1 +
              eta*g1*rinv*r0inv2 - 2r0inv2
        c22 = rinv*(-g1 - g0*g2/g1 + g2*c2*rinv)
        c25 = k*rinv*r0inv2*(-g2 + k*(c13 - g2*r0)*betainv*r0inv2 - c13*r0inv -
                             c12*c3*rinv*r0inv2 + c13*c2*c3*rinv2*r0inv2 -
                             c13*c20*betainv*rinv*r0inv2)
        c26 = k*rinv2*r0inv*(-g2*c12 - g1*c13 + g2*c13*c2*rinv)
        c34 = (-beta*eta^2*g2^2 - eta*k*h8 - h6*ksq - 2beta*eta*r0*g1*g2 +
               (g2^2 - 3*g1*g3)*beta*k*r0 - beta*g1^2*r0^2)*betainv*rinv2 +
              (eta*g2^2)*rinv/g1 + (k*h8)*betainv*rinv/g1
        # `h2` is zero on this branch in the primal; the Jacobian needs it.
        h2j = _H2(gamma, beta, sqb)
        c33 = d*k*rinv3*r0inv*k*(h*g2 - r0*g3) +
              k*(-eta*k*g1*g2^2 - g1*g2*g3*ksq - r0*eta*beta*g1*g2^2 -
                 r0*k*g1*h2j - beta*g2^2*g0*r0^2)*betainv*rinv2*r0inv

        C = _CoeffDerivs(
            dfdt*c21,                                           # ddfdtdxx
            dfdt*(c22 - h*c21),                                 # ddfdtdxv
            dfdt*(c22 - h*c21),                                 # ddfdtdvx
            dfdt*(c34 - 2*h*c22 + hsq*c21),                     # ddfdtdvv
            dfdt*(1/k - betainv*r0inv - c17*betainv*rinv*r0inv -
                  c1*(g1*c2 - g0*r)*rinv2*r0inv/g1),            # ddfdtdk
            dfdt*(g0*rinv/g1 - c2*rinv2 - (2*k*r0inv - beta)*c22 - eta*c21))
        D = _CoeffDerivs(
            c25,
            c26 - h*c25,
            c26 - h*c25,
            c33 - 2*h*c26 + hsq*c25,
            rinv*r0inv*(-k*(c13 - g2*r0)*betainv*r0inv + c13 -
                        k*c13*c17*betainv*rinv*r0inv + k*c1*c12*rinv*r0inv -
                        k*c1*c2*c13*rinv2*r0inv),
            g1*k*rinv*r0inv + k*c12*rinv2*r0inv - k*c2*c13*rinv3*r0inv -
                (2*k*r0inv - beta)*c26 - eta*c25)

        # Cancellation-free mass column.
        h4 = -_H1(gamma, beta)*beta
        dfm1dk2 = r0*h4 + k*h6
        dgmhdk2 = h6*g3*ksq + eta*r0*(h6 + g2*h4) + r0^2*g0*h5 + k*eta*g2*h6 +
                  (g1*h6 + g3*h4)*k*r0
        ddfdtdk2 = -(g2*k - r0)*(beta*r0*(g3 - g1*g2) - beta*eta*g2^2 + k*h3) *
                   betainv*rinv2*r0inv
        ddgdtdk2 = k*betainv*rinv2*r0inv*(
            -beta*eta^2*g2^4 + eta*g2*(g1*g2^2 + g1^2*g3 - 5*g2*g3)*k +
            g2*g3*h3*ksq + 2eta*r0*beta*g2^2*(g3 - g1*g2) +
            (4g3 - g0*g3 - g1*g2)*(g3 - g1*g2)*r0*k +
            beta*(2g1*g3*g2 - g1^2*g2^2 - g3^2)*r0^2)
        jm_x = (G*r0inv)^2 * betainv * rinv .* (dfm1dk2 .* x0 .- dgmhdk2 .* v0)
        jm_v = G^2 * r0inv * rinv .* (ddfdtdk2 .* x0 .+ ddgdtdk2 .* v0)
    else
        c9  = g2*r - h1*k
        c14 = r0*g2 - k*h1
        c15 = eta*h1 + h2*r0
        c16 = eta*h2 + g1*gamma*r0/sqb
        c19 = 4*eta*h1 + 3*h2*r0
        c23 = h2*k - r0*g1

        A = _CoeffDerivs(
            k*rinv3*betainv*r0inv^4*(k*h1*r^2*r0*(beta - 2*k*r0inv) +
                beta*c3*(r*c23 + c14*c2) + c14*r*(k*(r - g2*k) + g0*r0*zeta)),
            k*rinv2*r0inv*(k*(g2*h2 + g1*h1) - 2g1*g2*r0 + g2*c14*c2*rinv),
            k*rinv2*r0inv*(k*(g2*h2 + g1*h1) - 2g1*g2*r0 + g2*c14*c2*rinv),
            k*r0inv*rinv2*betainv*(2eta*k*(g2*g3 - g1*h1) + (3g3*h2 - 4h1*g2)*ksq +
                beta*g2*r0*(3h1*k - g2*r0) +
                c14*rinv*(-beta*g2^2*eta^2 + eta*k*(2g0*g3 - h2) - h6*ksq +
                          (-2eta*g1*g2 + k*(h1 - 2g1*g3))*beta*r0 - beta*g1^2*r0^2)),
            rinv*r0inv*(4*h1*ksq*betainv*r0inv - k*h1 - 2*g2*k*betainv + c14 -
                k*c14*c17*betainv*rinv*r0inv + k*(g1*r0 - k*h2)*c1*rinv*r0inv -
                k*c14*c1*c2*rinv2*r0inv),
            (g1*k - h2*ksq*r0inv - k*c14*c2*rinv*r0inv)*rinv2)
        B = _CoeffDerivs(
            k*rinv*r0inv*(h2 + k*c19*betainv*r0inv2 - c16*c3*rinv*r0inv2 +
                c2*c3*c15*(rinv*r0inv)^2 - c15*c20*betainv*rinv*r0inv2),
            k*rinv2*(h1*r - g2*c16 - g1*c15 + g2*c2*c15*rinv),
            k*rinv2*(h1*r - g2*c16 - g1*c15 + g2*c2*c15*rinv),
            k*betainv*rinv2*(2*eta^2*(g1*h1 - g2*g3) + eta*k*(4g2*h1 - 3h2*g3) +
                r0*eta*(4g0*h1 - 2g1*g3) + 3r0*k*((g1 + beta*g3)*h1 - g3*g2) +
                (g0*h8 - beta*g1*(g2^2 + g1*g3))*r0^2 -
                c15*rinv*(beta*g2^2*eta^2 + eta*k*h8 + h6*ksq +
                          (2eta*g1*g2 - k*(g2^2 - 3g1*g3))*beta*r0 + beta*g1^2*r0^2)),
            rinv*(k*c1*c16*rinv*r0inv + c15 - k*c15*c17*betainv*rinv*r0inv -
                  k*c19*betainv*r0inv - k*c1*c2*c15*rinv2*r0inv),
            k*rinv3*(r*c16 - c2*c15))

        C = _CoeffDerivs(
            dfdt*(eta*g1*rinv - 2 - g0*c3*rinv*r0inv/g1 + c2*c3*r0inv*rinv2 -
                  k*(k*g2 - r0)*betainv*rinv*r0inv)*r0inv2,
            -dfdt*(g0*g2/g1 + (r0*g1 + eta*g2)*rinv)*rinv,
            -dfdt*(g0*g2/g1 + (r0*g1 + eta*g2)*rinv)*rinv,
            -k*rinv3*r0inv*betainv*((beta*eta*g2^2 + k*h8)*(r0*g0 + k*g2) +
                g1*(-h6*ksq + (-2eta*g1*g2 + (h1 - 2g1*g3)*k)*beta*r0 -
                    beta*g1^2*r0^2)),
            dfdt*(1/k + c1*(r0 - g2*k)*r0inv*rinv2/g1 - betainv*r0inv*(1 + c17*rinv)),
            dfdt*(r0 - g2*k)*rinv2/g1)
        D = _CoeffDerivs(
            rinv2*r0inv3*((eta*g2 + g1*r0)*k*c3*rinv +
                          g2*k*(k*(g2*k - r) - g0*r0*zeta)*betainv),
            k*g2*rinv3*(r*g1 + r0*g1 + eta*g2),
            k*g2*rinv3*(r*g1 + r0*g1 + eta*g2),
            k*betainv*rinv3*(eta^2*beta*g2^3 - eta*k*g2*h3 + 3r0*eta*beta*g1*g2^2 +
                r0*k*(-g0*h6 + 3beta*g1*g2*g3) + beta*g2*(g0*g2 + g1^2)*r0^2),
            rinv*r0inv*(-r0*g2 + g2*k*(r + r0 - g2*k)*betainv*rinv - k*g1*c1*rinv +
                        k*g2*c1*c2*rinv2),
            k*rinv3*(g2*c2 - r*g1))

        h7 = beta*g1*g2^2 - g0*h8
        dfm1dk2 = betainv*r0inv*rinv2*(
            r*(2eta*k*(g1*h1 - g3*g2) + (4g2*h1 - 3g3*h2)*ksq - eta*r0*beta*g1*h1 +
               (g3*h2 - 4g2*h1)*beta*k*r0 + g2*h1*beta^2*r0^2) -
            c14*(-eta^2*beta*g2^2 - k*eta*h8 - ksq*h6 -
                 eta*r0*beta*(g1*g2 + g0*g3) + 2*(h1 - g1*g3)*beta*k*r0 -
                 (g2 - beta*g1*g3)*beta*r0^2))
        dgmhdk2 = betainv*rinv2*(
            r*(2eta^2*(g3*g2 - g1*h1) + eta*k*(3g3*h2 - 4g2*h1) +
               r0*eta*(beta*g3*(g1*g2 + g0*g3) - 2g0*h6) +
               (-h6*(g1 + beta*g3) + g2*(2g3 - h2))*r0*k +
               (h7 - beta^2*g1*g3^2)*r0^2) -
            c15*(-beta*eta^2*g2^2 + eta*k*(-h2 + 2g0*g3) - h6*ksq -
                 r0*eta*beta*(h2 + 2g0*g3) + 2beta*(2*h1 - g2^2)*r0*k +
                 beta*(beta*g1*g3 - g2)*r0^2))
        ddfdtdk2 = (r0 - g2*k)*betainv*r0inv*rinv2 *
                   (-eta*beta*g2^2 + h3*k + (g3 - g1*g2)*beta*r0)
        ddgdtdk2 = betainv*rinv2*(-beta*eta^2*g2^3 + eta*k*g2*h3 +
                   eta*r0*beta*g2*(g3 - 2g1*g2) + (h6 - beta*g2^3)*r0*k +
                   beta*g1*(g3 - g1*g2)*r0^2)
        jm_x = G^2 * rinv * r0inv .* (dfm1dk2 .* x0 .+ dgmhdk2 .* v0)
        jm_v = G^2 * rinv * r0inv .* (ddfdtdk2 .* x0 .+ ddgdtdk2 .* v0)
    end

    top = _jac_rows(x0, v0, fm1, gmh, A, B)
    bot = _jac_rows(x0, v0, dfdt, dgdtm1, C, D)
    J = vcat(top, bot)
    jac_mass = SVector{6,T}(jm_x[1], jm_x[2], jm_x[3], jm_v[1], jm_v[2], jm_v[3])
    return s.delx, s.delv, J, jac_mass
end
