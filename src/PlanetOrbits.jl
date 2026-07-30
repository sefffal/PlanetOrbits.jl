"""
# PlanetOrbits

A package for calculating the orbits of hierarchical systems of bodies in the
context of direct imaging, astrometry, and radial velocity.

Construct a `System` of `Body`s nested into `Orbit`s, solve it at a set of
epochs with `orbitsolve` (producing a `Trajectory`), and query observables
between any pair of references (bodies, barycentres, photocentres):

```julia
A = Body(mass=1.1, name=:A)
b = Body(mass=5mjup, name=:b)
sys = System(Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0); plx=24.5)
traj = orbitsolve(sys, epochs)
sol = traj[1]
raoff(sol, b, A)                    # relative astrometry [mas]
radvel(sol, A, barycentre(sys))     # stellar reflex RV [m/s]
```
"""
module PlanetOrbits

using LinearAlgebra
using StaticArrays
import ForwardDiff
using ForwardDiff: Dual, Partials, value, partials

include("constants.jl")
include("kepsolve.jl")
include("body.jl")
include("frame.jl")
include("system.jl")
include("trajectory.jl")
include("keplerian-approx.jl")
include("kepsolve-simd.jl")
include("nbody-kernels.jl")
include("ahl21.jl")
include("observables.jl")
include("sugar.jl")
include("time.jl")

end # module
