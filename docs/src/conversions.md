# Converting Between Orbit Types

You can convert between several of the supported orbit types.

Examples:
```@example 1
using PlanetOrbits, Plots


# Specify a Visual{KepOrbit}
orb_vis = Visual{KepOrbit}(M=1, e=0.4, a=1, i=2, Ω=3, ω=1, tp=0, plx=10.0);

# Convert to Thiele-Innes
orb_ti = ThieleInnesOrbit(orb_vis)

# Convert back to Visual{KepOrbit}
orb_vis2 = Visual{KepOrbit}(orb_ti)

# Convert to a CartesianOrbit (specified by position and velocity)
# We have to solve the orbit at a particular time (or true anomally, mean anomally, etc)
# Then we can use that solution to construct a CartesianOrbit
orb_vis_sol = orbitsolve(orb_vis,0)
orb_cart = CartesianOrbit(orb_vis_sol; tol=1e-4) # default is 1e-8

# Solve each orbit at the same date
time = mjd("2023-01-01")

sol_vis =   orbitsolve(orb_vis, time)
sol_ti =    orbitsolve(orb_ti, time)
sol_vis2 =  orbitsolve(orb_vis2, time)
sol_cart =  orbitsolve(orb_cart, time)

plot( aspectratio=1, legend=:none,)
xlims!(-1.5,1.5)
ylims!(-1.5,1.5)
zlims!(-1.5,1.5)

plot!(sol_vis,  color=1, kind=(:x,:y,:z))
plot!(sol_ti,   color=2, kind=(:x,:y,:z))
plot!(sol_vis2, color=3, kind=(:x,:y,:z))
plot!(sol_cart, color=4, kind=(:x,:y,:z))

scatter!([0], [0], [0], marker=:circle,  color=:white, ms=6)
```

When converting to a [`CartesianOrbit`](@ref), the `tol` parameter controls how near-equitorial and near-circular orbits are handled. 
If eccentricity is below `tol`, then it is zeroed and the orbit is treated as circular (changing how `ω` is set).
If the absolute value of inclination is below `tol`, then it is zeroed and the orbit is treated as equatorial (changing how `Ω` is set).