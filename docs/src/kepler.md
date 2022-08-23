# Kepler Solvers

The heart of this package is being able to take a set of Keplerian elements and output relative positions, velocities, etc.

This normaly requires solving Kepler's equation numerically. This package supports a multitude of solver algorithms that can be passed to [`orbitsolve`](@ref):

* [`PlanetOrbits.Auto`](@ref)
* [`PlanetOrbits.Markley`](@ref)
* [`PlanetOrbits.Goat`](@ref)
* [`PlanetOrbits.RootsMethod`](@ref)

The last of these `RootsMethod`, allows one to substitute any algorithm from the Roots.jl package. These include many different classical and modern root finding algorithms.chosen precision, including artibrary precision BigFloats. Using big floats with, for example, `Roots.PlanetOrbits.Thukral5B` and a tight tolerenace, allows you to solve orbits up to arbitrary precision.

The default choice is `Auto`, which currently selects `Markley` for all cases. The Markley algorithm is very fast, reasonably accurate, and always converges, making it a good default choice.


The Markley algorithm is slightly tweaked version of the algorithm in [AstroLib.jl](http://juliaastro.github.io/AstroLib.jl/stable/ref/#AstroLib.kepler_solver).
On my laptop, this solves for a single eccentric anomaly in just 71 ns.
Since it is implemented in pure Julia, there is no overhead from calling into a C or Cython compiled function and no need for vectorization.