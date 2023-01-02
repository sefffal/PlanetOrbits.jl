<img height=150 src="https://github.com/sefffal/PlanetOrbits.jl/blob/master/docs/src/assets/logo.png"/>

# PlanetOrbits.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sefffal.github.io/PlanetOrbits.jl/dev)
[![codecov](https://codecov.io/gh/sefffal/PlanetOrbits.jl/branch/master/graph/badge.svg?token=QLTCBWVV98)](https://codecov.io/gh/sefffal/PlanetOrbits.jl)

Tools for solving simple Keplerian orbits. 
The primary use case is mapping orbital elements into e.g. Cartesian coordinates at different times.
A Plots.jl recipe is included for easily plotting orbits.
One can for instance calculate an orbit around a star in 3D, a projected position in the sky, a radial velocity curve, or stellar astrometric accleration over time.

It's a great tool for visualizing different orbits (see examples) and generating nice animations (e.g. with Plots or Luxor.jl).
This package has been designed for good performance and composability with a wide range of packages in the Julia ecosystem, including ForwardDiff. 
It forms the backbone of [DirectDetections.jl](https://github.com/sefffal/DirectDetections.jl), a modelling framework for all kinds of exoplanet data.

See documentation at https://sefffal.github.io/PlanetOrbits.jl/dev
