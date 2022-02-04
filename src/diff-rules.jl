
# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# We try to support symbolic manipulation using Symbolics.jl, but it's
# not reasonable to use `remp2pi` on a symbolic variable.
# We therefore have a special fallback method for that case. We 
# define it when both packages get loaded by the user using Requires.
# @inline rem2pi_safe(x) = rem2pi(x, RoundNearest)
# Define a scale rule to allow autodiff to diff through rem2pi
# @scalar_rule rem2pi_safe(x) x