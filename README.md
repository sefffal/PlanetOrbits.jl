# DirectOrbits.jl

Tools for displaying, solving Keplerian orbits in the context of direct imaging.

See also [DirectImages.jl](//github.com/sefffal/DirectImages.jl)


# Usage
```julia
orbit = Orbit(
    a = 1.0,
    i = deg2rad(10),
    e = 0.0,
    τ = 0.0,
    μ = 1.0,
    ω = deg2rad(19),
    Ω = deg2rad(19),
    plx = 35.
)

# Display one full period of the orbit (requires `using Plots` beforehand)
plot(orbit)

# Get an Cartesian coordinate at a given epoch as an SVector{Float64}()
pos = xyz(orbit, 1.0) # at time in MJD (modified Julian days)

```
