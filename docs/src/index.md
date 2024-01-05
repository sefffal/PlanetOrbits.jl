
# PlanetOrbits.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/sefffal/PlanetOrbits.jl)


Tools for solving Keplerian orbits, especially in the context of exoplanet detection and modelling.
The functions of this package allow one to propagate two-body orbits using a variety of orbital basis sets (classic Campbell, Thiele-Innes, and Cartesian state-vectors.).
A fully featured A Plots.jl recipe is included for easily plotting different orbit properties in space or time, or against other variables.

Bound circular and elliptical orbits are fully supported. Support for hyperbolic orbits is experimental.

Among other uses, it can be used to calculates the projected positions of planets, radial velocity, and proper motion anomaly.
It is also a great tool for visualizing different orbits (see examples) and generating nice animations (e.g. with Plots or Luxor.jl).

This package has been designed for good performance and composability with a wide range of packages in the Julia ecosystem.
Automatic differentiation with ForwardDiff, Enzyme, and Zygote are supported. 

A variety of Kepler solvers are provided. Arbitrary precision can be achieved by specifying orbits using Julia's built in `BigFloat` datatype and using a solver with user-specified tolerance.

To fit orbits to observations, see [Octofitter.jl](https://github.com/sefffal/Octofitter.jl).

See also [AstroImages.jl](https://github.com/JuliaAstro/AstroImages.jl).

```@raw html
<video src="assets/51-eri-orbit.mp4" autoplay loop width=300 height=300>
```


### Tutorials
```@contents
Pages = ["introdcution.md", "plots.md", "image-warping.md"]
Depth = 5
```

### Documentation
```@contents
Pages = ["api.md", "conventions.md", "kepler.md"]
Depth = 5
```


## Attribution

If you find this package useful in your research, please cite the following [paper](https://dx.doi.org/10.3847/1538-3881/acf5cc) (open-access link).


This software package contains calculations that are adapted from various open source packages, including:
* NASA/JPL SPICE (public domain)
* keplerorbit.py by Spencer Wallace (MIT license)
* PoliaAstro (MIT license)
* Orbitize by Blunt et al. (BSD 3-Clause License)
* RadVel by Fulton et al. (MIT license)

These codes were useful references in the development of this package but are not redistributed.