
# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# We have analytic gradients for these already calculated. But with the above defintion
# they don't have much of an effect.
# @scalar_rule raoff(o::OrbitSolution) pmra(o)
# @scalar_rule decoff(o::OrbitSolution) pmdec(o)
# @scalar_rule propmotionanom(o::OrbitSolution) acceleration(o)
# @scalar_rule orbitsolve_ν(elem::KeplerianElements{T}, t; tref=58849)

# @scalar_rule OrbitSolution(xang, yang, ẋang, ẏang, żcart, ẍang, ÿang) OrbitSolution(ẋang, ẏang, ẍang, ÿang, nothing, nothing, nothing)