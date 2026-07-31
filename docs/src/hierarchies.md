# Hierarchical systems

A system is a flat list of bodies plus a flat list of orbits. Each orbit says
what orbits what, and nothing is inferred. That one decision is what lets the
same API describe a lone planet, a moon, a circumbinary planet, or a 2+2
quadruple.

## The reference grammar

Both endpoints of an orbit — the thing orbiting and the `about=` it orbits —
take the same grammar:

| spec | meaning |
|---|---|
| `A` | body `A` |
| `(A, b)` | the barycentre of `A` and `b` |

That is all of it. A body, or a set of bodies meaning their barycentre.

```@example hier
using PlanetOrbits
import PlanetOrbits as PO

A = PO.Body(mass=1.1, name=:A)
b = PO.Body(mass=8mjup, name=:b)
c = PO.Body(mass=2mjup, name=:c)
nothing # hide
```

## Jacobi vs. astrocentric

For a star with two planets there are two natural conventions, and they are
*different models*, not different spellings:

```@example hier
els1 = (; a=2.5, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
els2 = (; a=8.0, e=0.3, i=0.6, ω=0.4, Ω=2.0, tp=57000.0)

jacobi = PO.System((A, b, c), (
    PO.Orbit(b, about=A;      els1...),
    PO.Orbit(c, about=(A, b); els2...),   # c orbits the A+b barycentre
); plx=50.0)

astro = PO.System((A, b, c), (
    PO.Orbit(b, about=A; els1...),
    PO.Orbit(c, about=A; els2...),        # c orbits the star itself
); plx=50.0)
nothing # hide
```

Under the default `KeplerianApprox` propagator, each orbit is solved
independently — so **the orbits *are* the model**. The two systems above give
observably different predictions from identical element values:

```@example hier
ts = [58000.0, 59000.0, 60000.0]
tj = orbitsolve(jacobi, ts); ta = orbitsolve(astro, ts)
maximum(abs(raoff(tj[k], :c, :A) - raoff(ta[k], :c, :A)) for k in 1:3)  # mas
```

They also imply different gravitating masses per orbit — Jacobi's outer orbit
is bound by all three bodies, astrocentric's by only two:

```@example hier
(jacobi.rows[2].M, astro.rows[2].M)
```

There is deliberately no default and no inference from semi-major axis. Pick
one and say so.

!!! tip "Under `AHL21` the convention is pure bookkeeping"
    The N-body propagator uses the orbits only to build initial conditions, so
    any invertible set describing the same configuration integrates
    identically. Convention only changes the *model* under `KeplerianApprox`.

## Convention constructors

`Jacobi` and `Astrocentric` expand to the explicit spellings above, which is
easier to read for long chains:

```@example hier
jacobi2 = PO.System((A, b, c), Jacobi(A, b => els1, c => els2); plx=50.0)
jacobi2.Ainv ≈ jacobi.Ainv
```

`show` reports the convention it resolved:

```@example hier
astro
```

## Moons

A moon orbits its host, not the system barycentre. In v1 — and in any
tree-shaped representation — this is unreachable, because the host would have
to appear in two places at once. Here it is just another `about=`:

```@example hier
star = PO.Body(mass=1.0, name=:A)
planet = PO.Body(mass=10mjup, name=:b)
moon = PO.Body(mass=1mearth, name=:m)

sysm = PO.System((star, planet, moon), (
    PO.Orbit(planet, about=star;   a=5.2,   e=0.05, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
    PO.Orbit(moon,   about=planet; a=0.007, e=0.0,  i=0.4, ω=0.0, Ω=1.0, tp=58849.0),
); plx=50.0)

solm = orbitsolve(sysm, 58900.0)
(projectedseparation(solm, :m, :b), projectedseparation(solm, :b, :A))
```

## Circumbinary planets

A planet orbiting a tight binary is the barycentre spelling:

```@example hier
Aa = PO.Body(mass=0.6, name=:Aa)
Ab = PO.Body(mass=0.4, name=:Ab)
p  = PO.Body(mass=1mjup, name=:b)

cb = PO.System((Aa, Ab, p), (
    PO.Orbit(Ab, about=Aa;        a=0.2, e=0.05, i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
    PO.Orbit(p,  about=(Aa, Ab);  a=3.0, e=0.1,  i=0.4, ω=1.0, Ω=2.0, tp=58000.0),
); plx=80.0)
nothing # hide
```

## Set exteriors: 2+2 quadruples

Both endpoints take the grammar, so the *orbiting* side can be a barycentre
too. A 2+2 quadruple is two tight pairs orbiting each other:

```@example hier
Ba = PO.Body(mass=0.8, name=:Ba)
Bb = PO.Body(mass=0.7, name=:Bb)

quad = PO.System((Aa, Ab, Ba, Bb), (
    PO.Orbit(Ab, about=Aa; a=0.5, e=0.1, i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
    PO.Orbit(Bb, about=Ba; a=0.6, e=0.2, i=0.4, ω=0.3, Ω=0.2, tp=58849.0),
    PO.Orbit((Ba, Bb), about=(Aa, Ab); a=50.0, e=0.3, i=0.5, ω=0.4, Ω=0.3, tp=58000.0),
); plx=20.0)
```

## Mixed conventions

Mixing is legal, and `show` says so rather than pretending the system has one
convention. A moon spelled about its host, inside a chain that is otherwise
Jacobi, is a mixture — and that is a real modelling statement, not a
formatting detail.

## Validation

The orbit set has to determine every body's position. When it does not, the
error names the rows involved:

```@example hier
try
    PO.System((A, b, c), (PO.Orbit(b, about=A; a=1.0),))
catch err
    println(sprint(showerror, err))
end
```

```@example hier
try
    PO.System((A, b, c), (PO.Orbit(b, about=A; a=1.0), PO.Orbit(b, about=A; a=2.0)))
catch err
    println(sprint(showerror, err))
end
```

## Zero-mass bodies

Test particles are supported throughout. A massless set has no mass-weighted
barycentre, so its limit — the members' geometric centre, which for a single
particle is just its own position — is used instead. The star sits exactly at
the barycentre of a system whose companions are massless:

```@example hier
z = PO.System(
    (PO.Body(mass=1.0, name=:A), PO.Body(mass=0.0, name=:b)),
    (PO.Orbit(PO.Body(mass=0.0, name=:b), about=PO.Body(mass=1.0, name=:A);
              a=5.0, e=0.2, i=0.4, tp=58849.0),); plx=40.0)
abs(raoff(orbitsolve(z, 58900.0), :A, barycentre(z)))
```

## How it works

Internally each orbit contributes one row to a hierarchy matrix `H`, which maps
absolute barycentric positions to the relative coordinates the orbits describe.
Inverting `H` and dropping the barycentre column gives `sys.Ainv`, which turns
per-orbit relative states back into per-body absolute states.

This is one code path for every convention. It is a little slower than the
closed form a pure Jacobi chain admits — around 20 ns against a 25 µs
likelihood evaluation — and in exchange there is no second path in which an
astrocentric system can be silently evaluated with the Jacobi formula.
