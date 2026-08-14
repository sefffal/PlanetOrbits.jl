# Parametrizations

An orbit needs six numbers. Which six is up to you: `Orbit` accepts several
standard parametrizations, organised into groups. **Supply exactly one
alternative from each group.**

| group | alternatives | default |
|---|---|---|
| size | `a` [AU] \| `P` [days] | *required* |
| shape | (`e`, `ω`) \| (`secosω`, `sesinω`) \| (`ecosω`, `esinω`) | `e=0, ω=0` |
| phase | `tp` \| `M0` + `epoch` \| `θ` + `epoch` | `tp=0` |
| orientation | `i`, `Ω` | `0`, `0` |

Or replace all of them at once with [Cartesian initial conditions](@ref), or
replace size and orientation together with [Thiele-Innes constants](@ref).

Mistakes are mechanical:

```@example param
using PlanetOrbits
import PlanetOrbits as PO
A = PO.Body(mass=1.1, name=:A); b = PO.Body(mass=5mjup, name=:b)

try
    PO.Orbit(b, about=A; a=5.0, P=3000.0)
catch err
    println(sprint(showerror, err))
end
```

## Size: `a` or `P`

```@example param
sysa = PO.System((A, b), (PO.Orbit(b, about=A; a=8.0, e=0.1),))
sysp = PO.System((A, b), (PO.Orbit(b, about=A; P=period(sysa), e=0.1),))
(semimajoraxis(sysp), period(sysp))
```

`P` is converted with the orbit's own gravitating mass, so under different
conventions the same `P` gives a different `a` — as it should.

For a bound orbit both must be positive and finite: `a` sizes the conic and the
derived constants divide by it, so `a = 0` sends the mean motion `√(GM/a³)` to
`Inf` and `a = Inf` sends `J = 2πa/P` to `NaN`, in either case making every
observable `NaN` rather than failing. Both ends are reachable from an ordinary
prior — `a ~ Uniform(0, 100)` inverse-transforms to *exactly* `0.0` once the
sampler proposes a sufficiently negative unconstrained coordinate — so they are
rejected at construction with a
[`PlanetOrbits.OrbitDomainError`](@ref) naming the
value, which a likelihood catches as a quiet `-Inf`. A non-positive or
non-finite `P` is rejected the same way. (Hyperbolic orbits are the exception
that proves the rule: there `a < 0` is the convention — see below.)

!!! warning "`P` is in days"
    `period(sys)` returns days, so the two round-trip. Imaging users who think
    in years get a plausible-looking 365× error rather than a crash, so `show`
    prints the period in both units.

## Shape

`secosω` = √e·cosω and `ecosω` = e·cosω sample the eccentricity *disc* rather
than the half-plane, which removes the ω degeneracy as e → 0. That is usually
what you want when sampling.

```@example param
e, ω = 0.3, 1.1
s1 = PO.System((A, b), (PO.Orbit(b, about=A; a=5.0, e=e, ω=ω),))
s2 = PO.System((A, b), (PO.Orbit(b, about=A; a=5.0, secosω=√e*cos(ω), sesinω=√e*sin(ω)),))
(eccentricity(s1), eccentricity(s2))
```

## Phase

`tp` is the canonical element — the epoch of periastron passage, in MJD. The
alternatives are measured at an epoch you supply:

- `M0` — mean anomaly at `epoch` [rad].
- `θ` — sky-plane position angle at `epoch` [rad].

```@example param
base = PO.System((A, b), (PO.Orbit(b, about=A; a=5.0, e=0.4, ω=1.1, i=0.5, Ω=2.2, tp=59000.0),); plx=25.0)
ep = 59211.0
θ_obs = posangle(orbitsolve(base, ep), :b, :A)

viaθ = PO.System((A, b),
    (PO.Orbit(b, about=A; a=5.0, e=0.4, ω=1.1, i=0.5, Ω=2.2, θ=θ_obs, epoch=ep),); plx=25.0)
(periastron(base), periastron(viaθ))
```

`θ` is often the best-constrained quantity for a directly imaged planet with a
short observational arc, which is why it is worth parametrizing on directly.
Recovering `tp` from it needs no mass and no `a`: the deprojection cancels the
radius factor, and only the mean motion needs the period the constructor
already has.

!!! note "`τ` is not accepted"
    The dimensionless phase `τ` needs hidden period and reference-epoch state
    and has no clean meaning under N-body integration. Use `tp` or `M0`.

## Cartesian initial conditions

Give the full relative state instead of elements. This determines everything,
so it replaces every group above:

```@example param
sol = orbitsolve(base, 59123.0)
cart = PO.System((A, b), (PO.Orbit(b, about=A;
    x  = posx(sol, :b, :A), y  = posy(sol, :b, :A), z  = posz(sol, :b, :A),
    vx = velx(sol, :b, :A), vy = vely(sol, :b, :A), vz = velz(sol, :b, :A),
    epoch = 59123.0),); plx=25.0)
(semimajoraxis(cart), eccentricity(cart), periastron(cart))
```

Positions are AU and velocities AU / julian year — the same frame and units as
`posx` and `velx` — so a state read out of a solution reconstructs the same
orbit exactly, as above.

Because `a` comes from vis-viva and `e` from the eccentricity vector, unbound
states need no special handling: they simply come out with `a < 0` and `e > 1`.

!!! tip "Elements ↔ state vector"
    The conversion in both directions — six elements to a position-and-velocity
    pair and back, via the eccentricity vector and the vis-viva equation — is
    worked through step by step in
    [Classical Orbital Elements and the State Vector](https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html).

## Hyperbolic orbits

`e > 1` is supported, with velocities. The semi-major axis is negative by
convention; a positive value is taken as |a|.

```@example param
hyp = PO.System((A, PO.Body(mass=0.0, name=:b)),
    (PO.Orbit(PO.Body(mass=0.0, name=:b), about=A;
              a=-5.0, e=1.4, i=0.5, ω=1.1, Ω=2.2, tp=59000.0),); plx=25.0)
(period(hyp), semimajoraxis(hyp))
```

The period is `Inf`, and `P=` is rejected for unbound orbits. `e == 1` exactly
is also rejected — the elements are degenerate for parabolae, so use Cartesian
initial conditions there. (For why `e = 1` is a genuinely singular case rather
than an implementation gap, see
[Parabolic Trajectories](https://orbital-mechanics.space/the-orbit-equation/parabolic-trajectories.html);
the bound and unbound cases are
[Elliptical](https://orbital-mechanics.space/the-orbit-equation/elliptical-orbits.html)
and
[Hyperbolic](https://orbital-mechanics.space/the-orbit-equation/hyperbolic-trajectories.html).)

## Thiele-Innes constants

`ThieleInnes` returns the size and orientation elements to splat in, so it
composes with the shape and phase groups:

```@example param
ti = thieleinnes(base)
```

```@example param
viaTI = PO.System((A, b),
    (PO.Orbit(b, about=A; ThieleInnes(; ti...)..., e=0.4, tp=59000.0),); plx=25.0)
(semimajoraxis(viaTI), inclination(viaTI))
```

Pass `plx` if your constants are in mas rather than AU. This parametrization
has no coordinate singularity at `e → 0` or `i → 0`, which is why Gaia-only
astrometric fits use it.

!!! note "The node is ambiguous by ±180°"
    `(ω, Ω)` and `(ω+π, Ω+π)` give *identical* Thiele-Innes constants, so the
    inverse is genuinely two-valued. The two branches share a sky-plane track
    and have opposite line-of-sight motion — the standard astrometric node
    ambiguity, which only radial velocities break. `ThieleInnes` returns the
    branch with `Ω ∈ [0, π)`.

## Radial-velocity-only fits

Radial velocities constrain `m·sin i`, never `m` and `i` separately, and give
no handle on `Ω` or the parallax. The convention is `i = π/2`, `Ω = 0`, and no
parallax — so angular observables are unavailable rather than silently wrong:

```@example param
using PlanetOrbits: rvorbit
rv = rvorbit(M=1.1, msini=2mjup, P=400.0, e=0.2, ω=1.0, tp=59000.0)
radvel(orbitsolve(rv, 59100.0), :A, barycentre(rv))
```

Under this convention the secondary's mass *is* its m·sin i; reading it as a
true mass gives a lower bound.

## The `M=` compatibility override

An orbit's gravitating mass normally comes from the bodies it binds. `M=`
overrides it:

```@example param
ov = PO.System((A, b), (PO.Orbit(b, about=A; a=8.0, M=2.5),))
(ov.rows[1].M, period(ov))
```

This is **physically inconsistent** — it decouples Kepler's third law from the
masses that drive the reflex amplitudes — and exists only to reproduce
published fits bit-for-bit and to match orbitize!/RadVel conventions. `show`
labels it as a compatibility override. It is not a modelling choice.
