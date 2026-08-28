
# PlanetOrbits.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/sefffal/PlanetOrbits.jl)

Tools for calculating orbits of hierarchical systems of bodies, in the context
of direct imaging, astrometry, and radial velocity.

The package is built around one idea:

> A propagator produces per-body barycentric Cartesian states at requested
> epochs, and every observable is a function of a pair of *references* —
> weighted combinations of body states — plus system-level frame information.

Everything else follows from it. Observables are pairwise queries between any
bodies, barycentres or photocentres; propagators are swappable behind a single
output format; and the whole construct → solve → query path is
allocation-free, so it can sit inside an MCMC hot loop. Automatic
differentiation with ForwardDiff flows through construction, both propagators,
and the observables.

To fit orbits to observations, see
[Octofitter.jl](https://github.com/sefffal/Octofitter.jl).

## Installation

These docs describe **PlanetOrbits v2**. If you are coming from v1, see the
[migration guide](@ref migration). Install into a fresh project environment:

```julia
using Pkg
Pkg.activate("my-orbits")     # a fresh, dedicated environment
Pkg.add("PlanetOrbits")
```

To fit orbits to observations, install
[Octofitter.jl](https://github.com/sefffal/Octofitter.jl) v9, which is built on
PlanetOrbits v2 and installs it automatically.

## Quick start

```julia
using PlanetOrbits
import PlanetOrbits as PO       # System, Body and Orbit are not exported

A = PO.Body(mass=1.1, name=:A)             # [M⊙]
b = PO.Body(mass=5mjup, name=:b)

sys = PO.System((A, b), (
    PO.Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
); plx=24.5)

traj = orbitsolve(sys, [58900.0, 59000.0, 59100.0])   # sorted epochs [MJD]

raoff(traj[1], b, A)                    # separation of b from A [mas]
radvel(traj[1], A, barycentre(sys))     # stellar reflex velocity [m/s]
```

Every orbit names its reference explicitly with `about=`. That is what lets
one API describe a lone planet, a moon, a circumbinary planet or a 2+2
quadruple — see [Hierarchical systems](@ref).

!!! note "`System`, `Body` and `Orbit` are not exported"
    Octofitter.jl owns those three names unqualified, and its versions build
    these. Either qualify them (`import PlanetOrbits as PO`) or import them
    explicitly (`using PlanetOrbits: Body, Orbit, System`). Everything else —
    `orbitsolve`, the observables, `bodies`, `barycentre`, `photocentre` — is
    exported normally.

## What is supported

- Bound circular and elliptical orbits, and **unbound hyperbolic orbits**
  (`e > 1`), including velocities.
- Hierarchical systems of any shape: Jacobi chains, astrocentric sets, moons,
  circumbinary planets, 2+2 quadruples, and mixtures of conventions.
- Two propagators — `KeplerianApprox` (each hierarchy row on an independent
  Keplerian orbit) and `AHL21` (full symplectic N-body) — behind one interface,
  so switching costs one keyword.
- Arbitrary precision, by building the system from `BigFloat` values and
  choosing a solver with a user-specified tolerance.
- Plotting through a [Makie](https://docs.makie.org) extension: orbit tracks,
  calendar-date epoch axes, and shared axis-label metadata.

## Where to go next

If you are new to orbital elements and the two-body problem, the free online
text [Orbital Mechanics & Astrodynamics](https://orbital-mechanics.space/) by
Bryan Weber is an excellent companion; the tutorials below link into the
relevant sections as they come up.

### Tutorials
```@contents
Pages = ["introduction.md", "hierarchies.md", "parametrizations.md", "nbody.md", "plots.md"]
Depth = 5
```

### Documentation
```@contents
Pages = ["api.md", "conventions.md", "precision.md", "kepler.md", "migration.md"]
Depth = 5
```

## Attribution

If you find this package useful in your research, please cite the following
[paper](https://dx.doi.org/10.3847/1538-3881/acf5cc) (open-access link).

The N-body propagator is the symplectic map of Agol, Hernandez & Langford
(2021), MNRAS 507, 1582 ([arXiv:2106.02188](https://arxiv.org/abs/2106.02188)),
whose kernels are merged into this package with the authors' agreement from
[NbodyGradient.jl](https://github.com/ericagol/NbodyGradient.jl) (MIT,
© Eric Agol). Please cite that paper when publishing results computed with
`AHL21`. The hierarchy-matrix formalism follows Hamers & Portegies Zwart
(2016).

This software package contains calculations that are adapted from various open
source packages, including:
* NASA/JPL SPICE (public domain)
* keplerorbit.py by Spencer Wallace (MIT license)
* PoliaAstro (MIT license)
* Orbitize by Blunt et al. (BSD 3-Clause License)
* RadVel by Fulton et al. (MIT license)

These codes were useful references in the development of this package but are
not redistributed.

```@raw html
<video src="assets/51-eri-orbit.mp4" autoplay loop width=300 height=300>
```
