# DirectOrbits.jl

Tools for displaying, solving Keplerian orbits in the context of direct imaging.

See also [DirectImages.jl](//github.com/sefffal/DirectImages.jl)


# Usage
```julia
orbit = Orbit(
    a = 1.0,
    i = deg2rad(45),
    e = 0.25,
    τ = 0.0,
    μ = 1.0,
    ω = deg2rad(0),
    Ω = deg2rad(120),
    plx = 35.
)

# Display one full period of the orbit (requires `using Plots` beforehand)
plot(orbit label="My Planet")
```
![Orbit Plot](docs/orbit-sample.png)

```julia
# Get an Cartesian coordinate at a given epoch as an SVector{Float64}()
julia> pos = xyz(orbit, 1.0) # at time in MJD (modified Julian days)
3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
  0.02003012254093835
  0.01072871196981525
 -0.019306398386215368
```
