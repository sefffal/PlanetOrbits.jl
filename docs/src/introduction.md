# Introduction

This page walks through the whole path: build a system, solve it at a set of
epochs, and query observables.

## Bodies

A [`PlanetOrbits.Body`](@ref) is a point mass. Masses are in solar masses; the
constants `msun`, `mjup` and `mearth` are exported so you can write the unit
you mean.

```@example intro
using PlanetOrbits
import PlanetOrbits as PO

A = PO.Body(mass=1.1, name=:A)
b = PO.Body(mass=5mjup, name=:b)
nothing # hide
```

The `name` is a label, and it is what observables use to find a body later.
Give every body one.

## Orbits

An [`PlanetOrbits.Orbit`](@ref) is one Keplerian relationship: the orbit of
something **about** something else. The reference is always explicit.

```@example intro
o = PO.Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
nothing # hide
```

Angles are radians, `a` is in AU, and `tp` — the epoch of periastron passage —
is an MJD. There is no `M` keyword: the gravitating mass of an orbit is the
total mass of the bodies it binds, taken from the `Body` values themselves.

You can give `P` [days] instead of `a`, and there are alternative
parametrizations for the eccentricity and the orbital phase — see
[Parametrizations](@ref).

## Systems

A [`PlanetOrbits.System`](@ref) is a list of bodies plus the orbits relating
them. There must be exactly one orbit fewer than there are bodies: each orbit
supplies one relative coordinate, and the remaining degree of freedom is the
system barycentre.

```@example intro
sys = PO.System((A, b), (o,); plx=24.5)
```

The keyword arguments define the system's frame:

| given | effect |
|---|---|
| nothing | physical-unit observables only (`posx`, `radvel`, …) |
| `plx` [mas] | angular observables in mas (`raoff`, `pmra`, …) |
| `ra`, `dec`, `plx`, `pmra`, `pmdec`, `rv`, `ref_epoch` | full absolute frame |

The absolute frame adds 3D-motion compensation and light-travel-time
correction, which matter for nearby, high-proper-motion stars.

## Solving

`orbitsolve` takes the system and a **sorted** vector of epochs, and returns a
`Trajectory` — per-body barycentric states at every epoch.

```@example intro
epochs = collect(range(58800.0, 59600.0, length=5))
traj = orbitsolve(sys, epochs)
sol = traj[1]
nothing # hide
```

Indexing a trajectory gives a per-epoch solution, which is a cheap view rather
than a copy. Solving all epochs at once is the primary entry point: it is what
N-body integration requires, and it lets the Kepler solves batch across epochs
under SIMD. There is a single-epoch convenience, `orbitsolve(sys, t)`, for
interactive use.

## Observables

Every observable is a query between two references. A reference can be a body
(by value, by `Symbol`, or by a handle from `bodies`), a barycentre, or a
photocentre.

```@example intro
raoff(sol, b, A)          # RA offset of b relative to A [mas]
```

```@example intro
radvel(sol, A, barycentre(sys))   # reflex radial velocity of the star [m/s]
```

The argument order reads as "of, relative to". These are all equivalent:

```@example intro
refs = bodies(sys)
(raoff(sol, b, A), raoff(sol, :b, :A), raoff(sol, refs.b, refs.A))
```

`bodies(sys)` returns persistent, `isbits` handles. Named `Body` values and
`Symbol`s resolve to them by name, and because the name is a type parameter
that resolution constant-folds away — so all three spellings cost the same in
a hot loop.

The available observables:

| function | units | needs |
|---|---|---|
| `posx`, `posy`, `posz` | AU | — |
| `velx`, `vely`, `velz` | AU / julian yr | — |
| `radvel` | m/s | — |
| `raoff`, `decoff` | mas | `plx` |
| `projectedseparation` | mas | `plx` |
| `posangle` | rad | — |
| `pmra`, `pmdec` | mas/yr | `plx` |
| `accra`, `accdec` | mas/yr² | `plx` |

Asking for an angular observable on a system with no parallax is an error
rather than a silently wrong number.

## Barycentres and photocentres

`barycentre(sys)` is the whole-system barycentre; `barycentre(sys, members...)`
is the barycentre of a subsystem. `photocentre(sys; band)` is the flux-weighted
point — what a blended astrometric source actually measures.

```@example intro
Al = PO.Body(mass=1.0, flux=(G=1.0,), name=:A)
bl = PO.Body(mass=0.2, flux=(G=1.0,), name=:b)
lum = PO.System((Al, bl),
    (PO.Orbit(bl, about=Al; a=4.0, e=0.1, i=0.5, tp=58849.0),); plx=30.0)
lsol = orbitsolve(lum, 59000.0)

# equal brightness ⇒ the photocentre sits midway between the two bodies
mid = (raoff(lsol, :A, barycentre(lum)) + raoff(lsol, :b, barycentre(lum))) / 2
(raoff(lsol, photocentre(lum), barycentre(lum)), mid)
```

Bodies carry per-band fluxes, so `photocentre(sys; band=:G)` selects a weight
set. Setting the host's flux to 1.0 makes every other body's flux a contrast
ratio.

## Orbit properties

```@example intro
(period(sys), semimajoraxis(sys), eccentricity(sys), totalmass(sys), distance(sys))
```

For a system with more than one orbit, pass a row index: `period(sys, 2)`.

## Allocation-free use

For hot loops, preallocate the trajectory once and reuse it. `orbitsolve!`
performs no allocation of its own, so with caller-owned storage the whole
construct → solve → query path is allocation-free — including under
ForwardDiff `Dual`s.

```@example intro
traj_buf = Trajectory(sys, epochs)

function total_sep(θ, buf)
    A = PO.Body(mass=θ[1], name=:A)
    b = PO.Body(mass=θ[2], name=:b)
    s = PO.System((A, b),
        (PO.Orbit(b, about=A; a=θ[3], e=θ[4], i=θ[5], tp=58849.0),); plx=24.5)
    orbitsolve!(buf, s)
    sum(raoff(x, :b, :A) for x in buf)
end

θ = [1.1, 5mjup, 8.0, 0.1, 0.5]
total_sep(θ, traj_buf)          # warm up
@allocated total_sep(θ, traj_buf)
```

Note the system is rebuilt from scratch on every call: `System` values are
cheap, immutable and `isbits`, and that is the intended usage under MCMC.

Gradients flow through the same path:

```@example intro
import ForwardDiff
ForwardDiff.gradient(t -> total_sep(t, Trajectory{eltype(t)}(
    PO.System((PO.Body(mass=t[1], name=:A), PO.Body(mass=t[2], name=:b)),
        (PO.Orbit(PO.Body(mass=t[2], name=:b), about=PO.Body(mass=t[1], name=:A);
                  a=t[3], e=t[4], i=t[5], tp=58849.0),); plx=24.5), epochs)), θ)
```
