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

```@docs
bodies
barycentre
photocentre
fluxes
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
| `accra`, `accdec` | mas/yr² | `plx` |

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

## Index

```@index
```
