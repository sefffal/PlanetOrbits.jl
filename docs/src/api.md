# API Documentation

## Constructing systems

```@docs
PlanetOrbits.Body
PlanetOrbits.Orbit
PlanetOrbits.System
```

### Convention constructors

```@docs
Jacobi
Astrocentric
```

### Alternative parametrizations

```@docs
ThieleInnes
thieleinnes
```

### Two-body conveniences

Both are deliberately unexported; opt in with
`using PlanetOrbits: orbit, rvorbit`.

```@docs
PlanetOrbits.orbit
PlanetOrbits.rvorbit
```

## References

Anywhere an observable expects a reference you may pass a `BodyRef` from
`bodies(sys)`, a named `Body` value, or a `Symbol`. Barycentres and
photocentres are `WeightedPoint`s, and so is any weight vector a likelihood
builds for itself — see [Blended sources & photocentres](@ref).
`framedirection` is the one reference that is a *direction* rather than a point in
space, which is what makes an observer-aware read absolute rather than
relative.

```@docs
bodies
barycentre
photocentre
fluxes
framedirection
PlanetOrbits.BodyRef
PlanetOrbits.WeightedPoint
```

## Solving

```@docs
orbitsolve
orbitsolve!
```

### Propagators

```@docs
KeplerianApprox
AHL21
```

## Observables

All observables take a solution and two references, read as
`f(sol, of, relative_to)`.

| function | units | requires |
|---|---|---|
| `posx`, `posy`, `posz` | AU | — |
| `velx`, `vely`, `velz` | AU / julian yr | — |
| `radvel` | m/s | — |
| `raoff`, `decoff` | mas | `plx` |
| `projectedseparation` | mas | `plx` |
| `posangle` | rad | — |
| `pmra`, `pmdec` | mas/yr | `plx` |

`raoff`, `decoff`, `projectedseparation` and `posangle` also take a fourth
argument — the observer's barycentric position [ICRS, AU] — which turns on the
annual–orbital (Kopeikin) coupling and exact parallax factors. Those forms
require a full absolute frame and a trajectory solved with
`observing_geometry=true`; see [Precision opt-outs](@ref).

The frame quantities `frame_ra`, `frame_dec`, `frame_pmra`, `frame_pmdec` and
`frame_rv` give the propagated barycentre frame at the solution epoch, and
require a full absolute frame.

## Kepler solvers

```@docs
PlanetOrbits.Auto
PlanetOrbits.Markley
PlanetOrbits.Goat
PlanetOrbits.HyperbolicHalley
PlanetOrbits.RootsMethod
```

## System properties

```@docs
period
totalmass
semimajoraxis
eccentricity
inclination
meanmotion
periastron
distance
```

## Plotting utilities

Backend-independent metadata and epoch grids — no plotting package required:

```@docs
plotinfo
plotlabel
paraminfo
orbit_track_epochs
plot_epochs
orbit_phase
```

Added by the Makie extension, which activates when any Makie backend is
loaded. See [Plotting](@ref).

```@docs
orbitlines!
orbit_theme
MJDConversion
add_mjd_axis!
```

## Index

```@index
```
