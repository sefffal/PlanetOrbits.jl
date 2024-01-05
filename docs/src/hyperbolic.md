# Hyperbolic Orbits

PlanetOrbits.jl has preliminary support for hyperbolic orbits.
They are currently supported with [`KepOrbit`](@ref) and [`CartesianOrbit`](@ref) but not [`ThieleInnesOrbit`](@ref).


```@example 1
using PlanetOrbits, Plots

# Specify orbit with a Campbell parameters (KepOrbit)
orb = orbit(M=1, e=1.1, a=1, i=2, Ω=3, ω=1, tp=mjd("2024-01-01"));
sol = orbitsolve(orb, mjd("2024-3-01"))
plot(sol, tspan=mjd("2024-01-01") .+ (-300,100))

```


```@example 1
using PlanetOrbits, Plots

# Specify orbit with a state vector (CartesianOrbit)
orb = orbit(
    x = 1.0,
    y = 0.3,
    z = 0.001,
    vx = 0,
    vy = 9,
    vz = 0.0,
    M = 1,
    tref = mjd("2024-01-01")
)
sol = orbitsolve(orb, mjd("2024-3-01"))
plot(sol, tspan=mjd("2024-01-01") .+ (-300,100))
```