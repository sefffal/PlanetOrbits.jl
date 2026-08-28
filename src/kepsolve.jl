# ---------------------------------------------------
# Kepler Equation Solvers
# ---------------------------------------------------
abstract type AbstractSolver end

"""
    PlanetOrbits.Auto()

Automatic choice of Kepler solver algorithm.
Currently defaults to PlanetOrbits.Markley()
"""
struct Auto <: AbstractSolver end

include("kepsolve-markley.jl")
include("kepsolve-hyperbolic.jl")

"""
    PlanetOrbits.RootsMethod(method::Roots.AbstractUnivariateZeroMethod, kwargs...)

Wraps a root finding method from Roots.jl. You can also pass keyword arguments
that will be forwarded to Roots to control the tolerance.

Examples:
```julia
method = PlanetOrbits.RootsMethod(Roots.Newton())
method = PlanetOrbits.RootsMethod(Roots.Thukral5B())
method = PlanetOrbits.RootsMethod(Roots.Bisection())
method = PlanetOrbits.RootsMethod(Roots.A42())
method = PlanetOrbits.RootsMethod(Roots.Newton(), rtol=1e-3, atol=1e-3)
```
"""
struct RootsMethod{M,K} <: AbstractSolver
    method::M
    kwargs::K
end
RootsMethod(method; kwargs...) = RootsMethod(method, kwargs)
include("kepsolve-roots.jl")

# Fallback kepler solver function.
# If algorithm is unspecified, select the best one here.
kepler_solver(MA, e) = kepler_solver(MA, e, Auto())
function kepler_solver(MA, e, ::Auto)
    # On the primal: the concrete solvers' `Dual` methods below re-check the
    # branch with `value(e)`, so classifying the `Dual` here would send
    # e = Dual(1.0, +∂) to `Markley` (lexicographically < 1 is false, so
    # `>= 1` picks hyperbolic — and the reverse for a negative partial) while
    # the method it lands in disagrees. One classification, on the value.
    if _primal(e) < 1
        kepler_solver(MA, e, Markley())
    else
        kepler_solver(MA, e, HyperbolicHalley())
    end
end

# ---------------------------------------------------
# Implicit differentiation of the Kepler solve.
#
# The eccentric anomaly is defined implicitly by Kepler's equation
# E - e sin E = M, so once the primal E is known its derivatives have the
# closed form dE = (dM + sin(E) de) / (1 - e cos E). Solving in the value
# domain and applying this rule is both faster and more accurate than
# grinding the iterative solver through Dual/tracked number types.
# ---------------------------------------------------

# ChainRules consumers (Zygote, etc.)
using ChainRulesCore
@scalar_rule PlanetOrbits.kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u, sin(Ω) / u)

# ForwardDiff Duals: dispatch on the argument types so Dual-valued solves
# run the solver once on primal values only. Defined per concrete solver
# type (not Auto, whose body branches on e and re-dispatches here; and not a
# Union, which would be ambiguous against each solver's own Real method).
# Only valid for the elliptical branch (E - e sin E = M); hyperbolic Duals
# fall through to the generic iterative path.
#
# The branch test is `_primal(e)`, not `value(e)`: one `value` strips a single
# Dual level, so under a nested Dual (a Hessian) it is still a `Dual` being
# compared to 1 — lexicographically, which is exactly the ordering this test
# exists to avoid. `_primal` recurses to the Float64 underneath, so every
# level of nesting classifies the same conic.
for S in (:Markley, :RootsMethod)
    @eval begin
        @inline function kepler_solver(M::Dual{Tg}, e::Dual{Tg}, method::$S) where {Tg}
            if _primal(e) >= 1
                return invoke(kepler_solver, Tuple{Real,Real,typeof(method)}, M, e, method)
            end
            E = kepler_solver(value(M), value(e), method)
            s, c = sincos(E)
            invu = inv(1 - value(e) * c)
            return Dual{Tg}(E, invu * partials(M) + (s * invu) * partials(e))
        end
        @inline function kepler_solver(M::Dual{Tg}, e::Real, method::$S) where {Tg}
            if _primal(e) >= 1
                return invoke(kepler_solver, Tuple{Real,Real,typeof(method)}, M, e, method)
            end
            E = kepler_solver(value(M), e, method)
            c = cos(E)
            invu = inv(1 - e * c)
            return Dual{Tg}(E, invu * partials(M))
        end
        @inline function kepler_solver(M::Real, e::Dual{Tg}, method::$S) where {Tg}
            if _primal(e) >= 1
                return invoke(kepler_solver, Tuple{Real,Real,typeof(method)}, M, e, method)
            end
            E = kepler_solver(M, value(e), method)
            s, c = sincos(E)
            invu = inv(1 - value(e) * c)
            return Dual{Tg}(E, (s * invu) * partials(e))
        end
    end
end

# The same rule in the form the propagator actually consumes — sin E and cos E
# rather than E itself. Given the *primal* root's sincos pair (`sE`, `cE`) and
# the Duals it came from, the Dual pair costs no transcendental at all:
#
#     ∂E = (∂M + sin E ∂e) / (1 − e cos E)
#     sin E = Dual(sin E, +cos E ∂E)      cos E = Dual(cos E, −sin E ∂E)
#
# Splitting it out this way removes a second full `sincos` that the propagator
# was paying: `kepler_solver`'s Dual methods above already evaluate sincos(E) to
# build ∂E, and `_anomaly_sincos` then called `sincos` on the Dual root they
# returned — recomputing a transcendental whose value was already in hand.
#
# It is also the seam that lets the primal root come from anywhere, including
# the SIMD batch kernel (kepsolve-simd.jl), since nothing here touches the
# solver. `e` carries the eccentricity's partials — zero for a Float64 row.
@inline function _dual_sincosE(MA::Dual{Tg,V,N}, e::Dual{Tg,V,N},
                               sE::V, cE::V) where {Tg,V,N}
    invu = inv(one(V) - value(e) * cE)
    ∂E = invu * partials(MA) + (sE * invu) * partials(e)
    return Dual{Tg}(sE, cE * ∂E), Dual{Tg}(cE, -sE * ∂E)
end
