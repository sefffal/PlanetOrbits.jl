module PlanetOrbitsForwardDiffExt

using PlanetOrbits
using ForwardDiff: ForwardDiff, Dual, value, partials
using Roots
# Current primal value takes 60ns to calculate.
# Derivative tracing through function takes 94ns.
# This seems too slow! Once the primal is known, we should only
# need a single sin or cos and some division/multiplication.


# Hacky: we have an annoying method abiguity between ::Dual being specific and ::SomeSolver being specific.
# the compiler doesn't know which one to choose! We want it to choose these ones.
# The best way I have found so far is to iterate at compile time through all subtypes and define gradient 
# rules for each. This won't work for new solvers added by a user.

function define_partials_for_solver(T::Symbol)
    @eval function PlanetOrbits.kepler_solver(M::Dual{T}, e::Real, method::PlanetOrbits.$T) where T
        if value(e) >= 1
            @error "diff rules need to be updated for unbound orbits. Review implicit derivative of hyperbolic keplers eqn."
        end
        EA = PlanetOrbits.kepler_solver(value(M),e,method)
        temp = 1 - e*cos(EA)
        return Dual{T}(
            EA,
            partials(M)/temp
        )
    end
    @eval function PlanetOrbits.kepler_solver(M::Real, e::Dual{T}, method::PlanetOrbits.$T) where T
        if value(e) >= 1
            @error "diff rules need to be updated for unbound orbits. Review implicit derivative of hyperbolic keplers eqn."
        end
        EA = PlanetOrbits.kepler_solver(M,value(e),method)
        sea, cea = sincos(EA)
        temp = 1 - value(e)*cea
        return Dual{T}(
            EA,
            partials(e)*sea/temp
        )
    end
    @eval function PlanetOrbits.kepler_solver(M::Dual{T}, e::Dual{T}, method::PlanetOrbits.$T) where T
        if value(e) >= 1
            @error "diff rules need to be updated for unbound orbits. Review implicit derivative of hyperbolic keplers eqn."
        end
        EA = PlanetOrbits.kepler_solver(value(M),value(e),method)
        sea, cea = sincos(EA)
        temp = 1 - value(e)*cea
        return Dual{T}(
            EA,
            partials(M)/temp + partials(e)*sea/temp,
        )
    end
end


define_partials_for_solver(:Goat)
define_partials_for_solver(:RootsMethod)

# define_partials_for_solver(:Markley)
# Shocker! Currently it's faster to diff through the Markley algorithm than it is to run it and then compute a 
# a single `sincos` call. SIMD is amazing! 
# We leave this implementation here for future in case these peformance tradeoffs change.


end