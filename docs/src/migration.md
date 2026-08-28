# [Migrating from v1](@id migration)

This page is for people with working v1 code. If you are learning the package
now, skip it — nothing here is needed to use v2.

## The one idea that changed

In v1 an **orbit was the model**. `Visual{KepOrbit}(M=1.1, a=8.0, plx=24.5, …)`
carried the total mass, the elements and the distance in a single object, and
which of the seven orbit types you had determined which observables you could
ask for.

In v2 a **system is the model**. Mass lives on bodies, distance and sky
position live on the system's frame, and orbits say only what orbits what:

```julia
A = PO.Body(mass=1.1, name=:A)
b = PO.Body(mass=10mjup, name=:b)
sys = PO.System((A, b), (PO.Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2,
                                  tp=58849.0),); plx=24.5)
```

Everything below follows from that. There is one type instead of seven, the
capability tables are gone (a system either has a parallax or it does not), and
adding a third body is adding a third body rather than reaching for a different
representation.

## Constructors

| v1 | v2 |
|---|---|
| `KepOrbit(M=…, a=…, e=…, i=…, ω=…, Ω=…, tp=…)` | `System((A, b), (Orbit(b, about=A; a=…, e=…, i=…, ω=…, Ω=…, tp=…),))` |
| `Visual{KepOrbit}(…, plx=…)` | the same, plus `; plx=…` on `System` |
| `AbsoluteVisual{KepOrbit}(…, ra=…, dec=…, pmra=…, pmdec=…, rv=…, ref_epoch=…)` | the same, plus the full frame keyword set on `System` |
| `RadialVelocityOrbit(M=…, a=…, e=…, ω=…, tp=…)` | `PlanetOrbits.rvorbit(M=…, msini=…, a=…, e=…, ω=…, tp=…)` |
| `ThieleInnesOrbit(A=…, B=…, F=…, G=…, e=…, tp=…, plx=…)` | `Orbit(b, about=A; ThieleInnes(; A, B, F, G, plx)..., e=…, tp=…)` |
| `CartesianOrbit(x=…, …, vz=…, M=…, tref=…)` | `Orbit(b, about=A; x=…, …, vz=…, epoch=…)` |
| `orbit(a=…, M=…, plx=…)` | `PlanetOrbits.orbit(a=…, M=…, plx=…)` — same behaviour, no longer exported |

`Body`, `Orbit` and `System` are **not exported**, because Octofitter owns
those three names unqualified. Either qualify them (`import PlanetOrbits as PO`)
or import them explicitly (`using PlanetOrbits: Body, Orbit, System`).

The v1 orbit-type conversions (`ThieleInnesOrbit(orb_vis)`,
`CartesianOrbit(sol)`) have no v2 equivalent, because there is nothing to
convert *between*. Their two useful directions survive as functions:
[`thieleinnes(sys)`](@ref) returns the constants of a system, and a Cartesian
state read out of a solution reconstructs the same orbit exactly — see
[Cartesian initial conditions](@ref).

## Masses

`M` is no longer an orbital element. Bodies carry mass, and an orbit's
gravitating mass is the total mass of the bodies it binds. So the v1 pattern of
passing the total mass to the orbit and then the planet mass separately to each
observable is gone.

Masses are uniformly in solar masses, with `msun`, `mjup` and `mearth` exported
as plain multipliers: `mass=10mjup`.

`M=` survives on `Orbit` as a labelled compatibility override for reproducing
published fits bit-for-bit. It is physically inconsistent — it decouples
Kepler's third law from the masses that drive the reflex amplitudes — and
`show` labels it as such. See [The `M=` compatibility override](@ref).

## Solving

v1 solved one epoch at a time. v2 solves a **sorted vector** of epochs and
returns a `Trajectory`; indexing it gives the per-epoch solution v1 code
expects.

```julia
traj = orbitsolve(sys, epochs)   # sorted Vector, the primary entry point
sol  = traj[1]
sol  = orbitsolve(sys, 59000.0)  # single-epoch convenience, for interactive use
```

Solving all epochs together is what N-body integration requires and is what
lets the Kepler solves batch under SIMD, so it is worth restructuring loops to
use it rather than calling the single-epoch form repeatedly.

`orbitsolve_ν` and `orbitsolve_meananom` have no v2 equivalent.

## Observables

Every observable is now a query between two references, read as
`f(sol, of, relative_to)`:

| v1 | v2 |
|---|---|
| `raoff(sol)` | `raoff(sol, b, A)` |
| `radvel(sol)` | `radvel(sol, b, A)` |
| `radvel(sol, M_planet)` — the star's reflex | `radvel(sol, A, barycentre(sys))` |
| `pmra(sol, M_planet)` | `pmra(sol, A, barycentre(sys))` |
| — | `radvel(sol, b, A)` now includes the Einstein term; `velz` is the kinematic quantity |
| — | `raoff(sol, b, A, obs_pos)` — from an observer that is not at the SSB |

That third row is the substantive one. In v1 the star's reflex motion was
obtained by scaling the relative orbit by a mass ratio you supplied at the call
site; in v2 the bodies already carry mass, so the reflex is just a query
against the barycentre and there is no ratio to get wrong.

The one-argument forms (`raoff(sol)`) still work, but only for two-body
systems, where there is exactly one thing they could mean.

**Not carried over:** `accra`, `accdec`, `meananom`, `trueanom`, `eccanom`,
`semiamplitude`. For orbital phase use [`orbit_phase`](@ref), which returns the
mean-anomaly phase in `[0, 2π)`.

## Angles and epochs

Unchanged: angles are radians, epochs are MJD, `a` is AU, `P` and `period` are
days. The dimensionless phase `τ` is not a constructor keyword. If you have one
— from an older PlanetOrbits release, or from orbitize! — expressed as a
fraction of a period past a reference epoch `tref`, spell it as the equivalent
mean anomaly at that epoch:

```julia
Orbit(b, about=A; a=…, M0=2π * τ, epoch=tref)
```

See [Parametrizations](@ref) for the other phase and shape alternatives, which
are constructor *groups* rather than aliases.

## Plotting moved to Makie

v1's Plots.jl recipes — `plot(orb, kind=(:x, :y))`, `kind=:radvel`,
`body=(:primary, :secondary)` — are gone. v2 has a Makie extension instead:
`lines(sys, :b, :A)` for a quick look, `orbitlines!` for a phase-coloured
track, and plain Makie over the observables for anything else. The `kind`
keyword has no equivalent because you now name the observables you want
directly. See [Plotting](@ref).

## Numerical differences

Results will not be bit-identical to v1. In rough order of how likely you are
to notice:

- **The barycentric light-travel sign fix.** v1 had the barycentric
  subtraction the wrong way round, which inverted the sign of the whole
  correction. This is a bug fix rather than a precision tier, so it applies in
  every mode and there is no flag that restores the old behaviour. It affects
  systems built with a full absolute frame; parallax-only and frameless
  systems are untouched by it.
- **Observing geometry is modelled, and on by default** — viewing-direction
  rotation, differential (per-body) light-travel time, line-of-sight
  projection, and a per-body rather than per-system AU→mas scale.
  `observing_geometry=false` selects exactly the v1 geometry. See
  [Precision opt-outs](@ref) for how large each correction is and when it is
  worth declining.

  Concretely, in the regression suite that replays v1 fixtures: with the
  defaults, results agree with v1 to a relative 3e-3 (frameless and
  parallax-only) and 3e-2 (full absolute frame). With
  `observing_geometry=false` the frameless and parallax-only cases return to
  agreement at 1e-13 — the observing pass really is the only thing that
  changed for them — while the absolute-frame cases keep a residual, which is
  the sign fix above.
- **`radvel` is now the spectroscopic velocity.** v1 returned the kinematic
  line-of-sight velocity; v2 adds the **Einstein term** — the second-order
  Doppler and gravitational-redshift difference between the two references.
  There is no flag, because no reduction pipeline can have removed the
  orbit-varying part of it (it depends on the sampled orbit) and the constant
  part is absorbed by the fitted offset either way. `velz` (in AU/julian yr)
  is unchanged and remains the kinematic quantity, so
  `velz(sol, b, A) * au2m / year2sec_julian` reproduces v1's number exactly.

  How much your results move depends entirely on which pair you difference,
  and it is three orders of magnitude:

  | use | Einstein term | *varying* part, which is what shifts a fit |
  |---|---|---|
  | stellar reflex, planet companion — `radvel(sol, A, barycentre(sys))` | 0.03–5.7 cm/s | ≤ 0.3 cm/s |
  | stellar reflex, stellar/BD companion | 0.8–1.7 m/s | m/s if eccentric |
  | relative RV, Jupiter at 1 AU, e = 0 | 4.4 m/s | 0 (absorbed by the offset) |
  | relative RV, hot Jupiter at 0.05 AU | 89 m/s | 0 (absorbed by the offset) |
  | relative RV, Jupiter at 1 AU, **e = 0.4** | 2.7–8.4 m/s | **5.6 m/s** |
  | relative RV, imaged brown dwarf, 30 AU, e = 0.3 | 0.10–0.22 m/s | 0.12 m/s |

  So an ordinary reflex-RV planet fit will not notice. A **relative**-RV fit of
  an eccentric companion may, and if it does, the v1 answer was the biased one.
  Two further consequences: masses now enter RV predictions (and their
  gradients) through the potential, and the fitted γ absorbs a slightly
  different set of constants — `v_sys²/2c` and the star's own *surface*
  redshift, which are still not modelled. Full tables in
  [Precision opt-outs](@ref).
- **Thiele-Innes inversion.** v1's `ThieleInnesOrbit` inverse carried documented
  π errors for `Ω ≥ π` and for `ω + Ω > 2π`. v2 inverts through half-angle sums,
  which have no quadrant fixups and no near-zero divisions, so those cases are
  now correct.
- **Hierarchy conventions are now explicit.** A v1 orbit was always a two-body
  relationship, so whatever a multi-planet model meant was decided by the code
  that assembled the orbits. v2 makes you write `about=A` or `about=(A, b)`,
  and the two are different models that give observably different predictions
  from the same element values. If you are porting a multi-planet fit, work out
  which one your v1 setup assumed — see [Jacobi vs. astrocentric](@ref).

## What v2 can do that v1 could not

Worth knowing about before you port a workaround:

- **Moons, circumbinary planets and 2+2 quadruples.** A v1 orbit described
  exactly two bodies, so a body that is simultaneously a host and a companion
  had no representation at all. See [Hierarchical systems](@ref).
- **Photocentres.** Blended catalogue sources are a first-class reference, so a
  partially resolved pair needs no hand-rolled flux weighting. See
  [Blended sources & photocentres](@ref).
- **N-body integration.** `AHL21` swaps in for the Keplerian approximation with
  one keyword and every observable works unchanged. See
  [N-body integration](@ref).
- **Hyperbolic orbits with velocities**, allocation-free. See
  [Hyperbolic orbits](@ref).
- **Observers that are not at the solar-system barycentre.** `raoff` and
  `decoff` take an optional per-epoch observer position, which turns on the
  annual–orbital (Kopeikin 1995) coupling and exact parallax factors by exact
  geometry rather than a series. Naming [`framedirection`](@ref) as the reference
  gives absolute astrometry the full parallax ellipse with no separately
  computed factors. See [Precision opt-outs](@ref).
- **Allocation-free construction and solving under ForwardDiff `Dual`s**, which
  is what makes the whole path usable inside an MCMC hot loop.
