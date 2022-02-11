
# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)
