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
It forms the backbone of [Octofitter.jl](https://github.com/sefffal/Octofitter.jl), a modelling framework for all kinds of exoplanet data.

See documentation at https://sefffal.github.io/PlanetOrbits.jl/dev

## N-body integration and attribution

PlanetOrbits v2 includes the AHL21 symplectic N-body propagator
(`method=AHL21(h=…)`), whose numerical kernels are derived from
[NbodyGradient.jl](https://github.com/ericagol/NbodyGradient) (© 2021 Eric
Agol, David Hernandez & Zach Langford, MIT license), merged into this package
with the authors' agreement. If you publish results computed with the AHL21
propagator, please cite:

> Agol, Hernandez & Langford (2021), MNRAS 507, 1582
> ([arXiv:2106.02188](https://arxiv.org/abs/2106.02188))

The hierarchical-system formalism (hierarchy matrix / A-matrix) follows
Hamers & Portegies Zwart (2016).
