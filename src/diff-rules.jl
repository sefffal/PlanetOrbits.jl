
# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)


# TODO: need to test to 100% verify these are right
# using DiffRules
# DiffRules.@define_diffrule PlanetOrbits.kepler_solver(M, e)  = :(u = 1 - $e*cos(kepler_solver($M, $e)); 1 / u), :(Ω=kepler_solver($M, $e); u = 1 - $e*cos(Ω); sin(Ω) / u)


# We have analytic gradients for these already calculated. But with the above defintion
# they don't have much of an effect.
# @scalar_rule raoff(o::OrbitSolution) pmra(o)
# @scalar_rule decoff(o::OrbitSolution) pmdec(o)
# @scalar_rule propmotionanom(o::OrbitSolution) acceleration(o)
# @scalar_rule orbitsolve_ν(elem::VisualOrbit{T}, t; tref=58849)
