# ---------------------------------------------------
# Hyperbolic Kepler equation:  M = e·sinh(H) − H,  e > 1
#
# `H` is the hyperbolic anomaly (the code keeps the `EA` naming used by the
# elliptical solvers). Roots.jl can solve this too — and did in v1 — but it
# allocates, which would put every hyperbolic system outside the
# allocation-free contract the rest of the package holds to. This is a fixed
# iteration count of Halley steps from a closed-form guess, so it is
# branch-light, allocation-free, and differentiable by the same implicit-rule
# trick used for the elliptical solvers.
# ---------------------------------------------------

"""
    PlanetOrbits.HyperbolicHalley()

Allocation-free Halley solver for the hyperbolic Kepler equation
(`M = e·sinh(H) − H`, `e > 1`). Selected automatically by `Auto()` when
`e > 1`; the elliptical solvers do not apply there.
"""
struct HyperbolicHalley <: AbstractSolver end

# Halley converges cubically; 12 steps is far past machine precision across
# the tested (M, e) grid and keeps the loop trip count compile-time constant
# so it unrolls rather than branching per iteration.
const _HYP_ITERS = 12

@inline function kepler_solver(_MA::Real, e::Real, ::HyperbolicHalley)
    T = float(promote_type(typeof(_MA), typeof(e)))
    MA = T(_MA)
    ecc = T(e)
    # asinh(M/e) is the large-|M| asymptote of the equation and stays finite
    # and well-scaled through M = 0, so it works as a global starting point.
    H = asinh(MA / ecc)
    for _ in 1:_HYP_ITERS
        sh = sinh(H)
        ch = cosh(H)
        f = ecc * sh - H - MA
        f1 = ecc * ch - one(T)          # f′  (> 0 for e > 1)
        f2 = ecc * sh                   # f″
        # Halley: x -= 2ff′ / (2f′² − ff″)
        denom = 2 * f1 * f1 - f * f2
        H -= (2 * f * f1) / denom
    end
    return H
end

# Implicit differentiation, hyperbolic branch. From M = e·sinh(H) − H,
#   dM = (e·cosh H − 1)·dH + sinh(H)·de
# so dH = (dM − sinh(H)·de) / (e·cosh H − 1). Note the sign on the `de` term
# differs from the elliptical rule — it is not the same formula with cos→cosh.
@inline function kepler_solver(M::Dual{Tg}, e::Dual{Tg}, method::HyperbolicHalley) where {Tg}
    H = kepler_solver(value(M), value(e), method)
    sh, ch = sinh(H), cosh(H)
    invu = inv(value(e) * ch - 1)
    return Dual{Tg}(H, invu * partials(M) - (sh * invu) * partials(e))
end
@inline function kepler_solver(M::Dual{Tg}, e::Real, method::HyperbolicHalley) where {Tg}
    H = kepler_solver(value(M), e, method)
    invu = inv(e * cosh(H) - 1)
    return Dual{Tg}(H, invu * partials(M))
end
@inline function kepler_solver(M::Real, e::Dual{Tg}, method::HyperbolicHalley) where {Tg}
    H = kepler_solver(M, value(e), method)
    sh, ch = sinh(H), cosh(H)
    invu = inv(value(e) * ch - 1)
    return Dual{Tg}(H, -(sh * invu) * partials(e))
end
