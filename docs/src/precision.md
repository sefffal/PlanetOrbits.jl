# Precision opt-outs

`orbitsolve` models source-side geometry well below the microarcsecond. Most
data does not need all of it, and the corrections that go unused are not free.
Two keyword arguments let you decline them:

```julia
orbitsolve(sys, epochs; observing_geometry=true, barycentric_lighttime=true)
```

Both default to `true`: PlanetOrbits cannot see how precise your measurements
are, so the accurate answer is the default and opting out is an explicit
assertion by the caller. They are keywords on `orbitsolve`/`orbitsolve!` rather
than fields of `System` because they describe the precision of the
*observations*, not the physical system.

**They gate different corrections, and one does not imply the other.** The
common mistake is to treat "high precision data" as a single axis. It is two.

| | `observing_geometry` | `barycentric_lighttime` |
|---|---|---|
| gates | how the barycentric states are *viewed* | the emission epoch `t_em` |
| scales with | angular extent **ρ** of the system on sky | proximity × proper motion |
| units of the error | µas of position | seconds of time |
| turn off when | ρ is small, or precision is coarse | the system is distant or slow |

A wide, nearby binary observed with GRAVITY needs the first and not the
second. A transiting system measured by transit timing needs the second and not
the first — its planets are unresolved, so ρ ≈ 0 and every angular correction
vanishes, while the timing correction does not.

## `observing_geometry`

Four corrections separate a barycentric state from what an observer at the
solar-system barycentre sees: rotation into the viewing direction at each
epoch, per-body (differential) light-travel time, a per-body rather than
per-system AU→mas scale, and the line-of-sight projection. All four are gated
by this one flag; the first three scale with the angular excursion actually
observed, ρ, and the fourth does not.

| correction | size | observable |
|---|---|---|
| viewing-direction rotation | 4.85e-3 · ρ[mas] · μ[″/yr] · T[yr] µas | position |
| differential (per-body) light-travel time | 0.099 · ρ[mas] · √(M/a[AU]) µas | position |
| depth scaling | 4.85e-6 · ρ[mas]² µas | position |
| line-of-sight projection | 0.023 · μ[″/yr] · a[AU] m/s | **radial velocity** |

For **absolute** astrometry ρ is the photocentre reflex, not the relative
orbit. A Jupiter analog at 10 pc gives ρ = 0.475 mas and therefore 0.005 /
0.021 / 1e-6 µas — four orders of magnitude below a 30–100 µas per-epoch
precision. For a resolved system at GRAVITY+ or CRIRES+ precision, ρ is the
full separation and the same formulas give tens of µas.

For the first three rows, then, the test is a single comparison of
`max(ρ) ×` the coefficients above against your declared measurement precision.
It is not a distance cut: distance enters only through ρ.

!!! warning "The fourth row does not scale with ρ, and is not in µas"
    Two bodies at different sky positions are seen along *different* unit
    vectors, so a velocity common to both — above all the barycentre's
    transverse space velocity — projects differently onto each. That
    contributes ≈ `0.023 · μ[″/yr] · a[AU]` m/s to a **relative** radial
    velocity: 24 cm/s for a Barnard-like host with a 1 AU companion, and
    **≈ 239 m/s** for the same host with a companion at 1000 AU. It is
    independent of distance, and phase-locked at the orbital period rather
    than secular.

    A cm/s-precision RV user who read only the µas table would turn this flag
    off and lose a term well above their noise floor. If you are fitting
    radial velocities of a high-proper-motion system — especially a wide pair
    — size this row, not the first three. It is also why no reduction
    pipeline can remove the term: it depends on the fitted barycentre frame
    and on each body's fitted direction.

`observing_geometry=false` selects the cheap geometry — one shared AU→mas scale
per epoch, no rotation, no retardation, no line-of-sight projection.

## `barycentric_lighttime`

The system as a whole moves in three dimensions, so its distance changes, so
light observed at `t_obs` was emitted at a different `t_em`. `orbitsolve` finds
`t_em` by solving

```
t_em + (d(t_em) − d_ref) / c = t_obs
```

and evaluates the orbit there. `barycentric_lighttime=false` takes
`t_em = t_obs` instead.

The raw size of this shift is large — a tenth of a day to well over a day — but
most of it is unobservable. The constant part is degenerate with `tp` and the
linear part with the period, so what actually changes a fit is only the
*curvature*, driven by the perspective acceleration v_tan²/d:

| system | raw shift | after removing the degenerate linear part |
|---|---|---|
| 125 pc, slow | 0.061 d | 0.01 s |
| 41 pc, 0.1″/yr | 0.30 d | 0.04 s |
| 8 pc, 1.7″/yr | 0.73 d | 2.0 s |
| 1.8 pc, 10.4″/yr (Barnard-like) | 1.34 d | 15.9 s |

Translate that into your own observable before deciding. For astrometry, a
timing error δt displaces a companion by roughly `ρ · 2π · δt/P`: even the
Barnard-like extreme is ~0.002 µas at ρ = 1 mas and P = 1 yr. For transit
timing, δt *is* the observable, and seconds are within reach of current
facilities — though note the shift is common-mode across every body in the
system, so it is partly absorbed by the fitted periods and relative timings
between planets are insensitive to it at first order.

!!! note
    `barycentric_lighttime=false` does **not** stop the frame being propagated.
    `frame_pmra`, `frame_pmdec` and `frame_rv` are still evaluated at each
    epoch, including their perspective acceleration — that curvature *is* the
    absolute-astrometry signal, and deleting it would silently change your
    model rather than merely cheapen it. The frame-level equivalent of "do not
    propagate at all" is not a solve-time flag; it is building the `System`
    with a `Parallax` frame (pass only `plx`) instead of an `AbsoluteFrame`
    (pass `ra`, `dec`, `pmra`, `pmdec`, `rv`, `ref_epoch`).

### The two settings are two conventions, and your data picks one

Hipparcos and Gaia catalog proper motions are **apparent** quantities — the
rate of the light-time-affected direction against the time of light *arrival*
(Butkevich & Lindegren 2014; the Gaia documentation says so explicitly) — and
the catalogs were reduced with the light-time-free standard model.

- `barycentric_lighttime=false` **is** that standard model: the catalog
  values propagated linearly as they stand. Use it whenever the quantities
  you compare against are catalog-convention absolute astrometry — catalog
  proper motions, Hipparcos–Gaia position differences, published abscissae.
  This is not an approximation being tolerated; it is the reduction
  convention of the data (B&L 2014, Sects. 5.5 and 6.1).
- `barycentric_lighttime=true` is the rigorous apparent path: the catalog
  proper motions are internally de-Dopplered to the true space velocity
  (`μ_true = μ_app · (1 + v_r/c)`, the same observed→inertial step as ERFA's
  `starpv`), the emission epoch is solved along that worldline, and the
  *angular* frame readouts — position, `frame_pmra`, `frame_pmdec` — are
  apparent quantities, `d/dt_obs`.

`frame_rv` deliberately stays in the spectroscopic convention on both paths:
it is the coordinate radial rate at the emission event, which is what a
Doppler measurement of that event corresponds to, not the rate at which the
light-time-affected distance changes against arrival time. That mirrors the
catalogs themselves, which publish apparent proper motions beside a
spectroscopic radial velocity — and it means `frame_rv` at `ref_epoch` is the
`rv` you passed in, whichever way the flag is set.

Both conventions reproduce the catalog values identically at `ref_epoch`, and
they differ away from it only by the genuine second-order light-time terms —
measured over ±25 yr at three of the most extreme catalog kinematics there
are: 0.10 mas / 0.008 mas/yr for Barnard's star, 0.16 mas / 0.013 mas/yr for
Kapteyn's star (faster and receding at 245 km/s, so the larger of the two),
and below a microarcsecond for a garden-variety host like ups And. (B&L's own
verdict for Barnard: light-time effects are negligible at the 1 mas level for
baselines ≪ 114 yr.) The flag is therefore a *convention* choice for catalog
work, and a physics switch only for absolute timing.
The `test/lighttime-bl2014.jl` oracle gates both paths against an independent
rigorous implementation.

## Rules of thumb

Turn `observing_geometry` **on** for resolved systems at µas precision —
GRAVITY-class relative astrometry, or anything where the separation reaches
tens of mas and the quoted uncertainty is below ~10 µas.

Turn `barycentric_lighttime` **on** for nearby, high-proper-motion systems, and
for any transit-timing or other absolute-timing data.

Turn both **off** for a large survey of distant unresolved systems at
Gaia-class precision, which is where the savings actually matter — but not if
that survey includes cm/s radial velocities of high-proper-motion hosts; see
the warning above.

When in doubt, leave them on and measure whether it costs you anything — the
defaults are the accurate ones. Octofitter automates exactly that measurement:
its `System` takes `observing_geometry=:auto`, which resolves the flag once at
build by comparing both settings against your own uncertainties over draws
from the priors, and reports what it decided.

## The Einstein term in `radvel` (modelled, always)

[`radvel`](@ref) returns the **spectroscopic** radial velocity — what a
spectrograph reports. `velz` is the kinematic quantity, for dynamics. That
distinction is *kinematic vs. spectroscopic*, not "coordinate vs. apparent":
`velz` already carries the line-of-sight projection.

The difference between them is the **Einstein term**: the second-order Doppler
(time-dilation) shift plus the gravitational redshift, differenced between the
two references. Per body,

```
Ein_i = ( ½|v_tot,i|² + Σ_{j≠i} G·mⱼ / r_ij ) / c
radvel(sol, t, r) = velz-difference + (Ein(t) − Ein(r))
```

with `v_tot` the body's **total** barycentric velocity — orbital plus the
frame's space velocity when the system has an `AbsoluteFrame`, orbital alone
for `Parallax`/`NoFrame`. Nothing emits from a barycentre, so `Ein` there is
zero and the stellar-reflex case keeps the star's own term in full; a
photocentre blends its members'.

**There is no flag for this, by design.** A precision opt-out is an assertion
about the *data* — "my reduction already removed this" or "my errors are
larger than this". Neither applies: the orbit-varying part of the Einstein
term depends on the sampled orbit (e, the masses, r(t)), so no reduction
pipeline can ever have removed it, and the constant part is absorbed by the
instrument offset whether or not it is modelled. It is therefore computed on
both settings of `observing_geometry` — that flag chooses the precision of the
geometry and is not permitted to change what `radvel` means.

Three consequences worth stating plainly:

1. **Masses now enter radial-velocity predictions**, and their gradients, via
   Φ. In a fit that couples the mass to the RV amplitude anyway this changes
   little; in a relative-RV fit of a directly-imaged companion it is a new,
   genuine constraint.
2. **`γ` means something slightly different.** The fitted systemic offset now
   absorbs `v_sys²/2c` and the star's own surface potential (see below)
   instead of those plus the terms now modelled.
3. **v0.11 results shift.** See the migration guide; the tables below are the
   size of the shift.

Their combined size for body 1 of a pair is

```
Δv = (1/c) · [ v₁²/2 + G·M₂/r ]
```

with a natural scale of `G·M₂/(a·c)` = **2.96 m/s · (M₂/M⊙) / (a/AU)**.

Which body you ask about changes the answer by three orders of magnitude, so
the two common uses of `radvel` have to be priced separately. These are the
tables to read when deciding whether the change from v0.11 moves your results,
and whether the term is worth anything to you.

**Stellar reflex —** `radvel(sol, A, barycentre(sys))`. Here `M₂` is the
companion mass and `v₁` the star's reflex, so both terms are tiny, and most of
what remains is constant:

| term | size | observable? |
|---|---|---|
| star's own *surface* gravitational redshift | 636 m/s (solar) | **not modelled** — needs a radius; constant, so absorbed by γ |
| `v_sys²/2c` | 1.5 m/s at 30 km/s | constant — absorbed by γ |
| `v_sys · v_orb / c` | 3 mm/s (30 m/s reflex), **3 m/s** (30 km/s reflex) | **varies at the orbital period** |
| `G·M₂/r + v_orb²/2c` | 0.03–5.7 cm/s (planets), **0.8–1.7 m/s** (stellar/BD) | constant if circular; varies with `e` |

So for a planet-mass companion these sit below current precision — a 1 AU,
e = 0.4 Jupiter varies by 0.3 cm/s over its orbit, and a circular hot Jupiter's
5.7 cm/s is constant and therefore absorbed into γ. For a **stellar or
brown-dwarf companion on an eccentric orbit they reach the m/s level** and are
phase-locked to the orbit, which is where they can bias elements rather than
merely shift γ. The eccentric variation scales as `2e/(1−e²) · G·M₂/(a·c)`.

**Relative RV —** `radvel(sol, b, A)`, the observable for directly-imaged
companions and transmission spectroscopy. Now the roles swap: the *companion*
sits deep in the *star's* potential and moves fast, so its Einstein term
carries `G·M_A/r` and its own orbital speed, and it dominates the difference:

| system | relative Einstein term | variation over the orbit |
|---|---|---|
| Jupiter-mass, 1 AU, e = 0 | 4.4 m/s | 0 (constant) |
| hot Jupiter, 0.05 AU, e = 0 | 89 m/s | 0 (constant) |
| Jupiter analog, 5.2 AU, e = 0.05 | 0.80–0.91 m/s | 0.11 m/s |
| Jupiter, 1 AU, **e = 0.4** | 2.7–8.4 m/s | **5.6 m/s** |
| imaged brown dwarf, 30 AU, e = 0.3 | 0.10–0.22 m/s | 0.12 m/s |

The constant column is absorbed by whatever velocity offset the reduction
already fits, so it changes nothing you fit. The **variation** is not, and for
a close-in eccentric companion it is several m/s — comparable to the precision
high-resolution cross-correlation is now reaching on directly-detected
companions. If you fit relative RVs of an eccentric companion, this is the
term on this page most likely to matter to you, and the one whose absence in
v0.11 most likely biased your elements.

## Observer-aware observables (opt-in, per read)

Everything above puts the observer at the solar-system barycentre. The
observer-aware forms of [`raoff`](@ref) and `decoff` take an explicit observer
position instead:

```julia
raoff(sol, b, A, obs_pos)      # relative astrometry from `obs_pos`
raoff(sol, A, framedirection, obs_pos)   # absolute: parallax ellipse included
```

`obs_pos` is the observer's barycentric position in ICRS Cartesian
coordinates, in AU, at that epoch — the Earth's, Gaia's at L2, or `(0, 0, 0)`
for the SSB, in which case they reduce exactly to the zero-argument forms.

The geometry is exact rather than a series in ϖ: each body's apparent
direction is computed from where the observer actually is, so the
annual–orbital (Kopeikin 1995) coupling and the exact per-body parallax
factors fall out of the same expression. Which part you get is decided by
which reference you name, not by a flag:

- against another **body** or a **barycentre** — points at the same distance —
  the first-order parallax is common to both and cancels, leaving only the
  differential (Kopeikin) term, `Δθ ≈ 4.85 · z[AU] / d[pc]²` µas per AU of
  observer displacement, with `z` the **line-of-sight** separation. Unlike
  every entry in the `observing_geometry` table this does *not* scale with the
  sky separation ρ: a face-on orbit at 1″ has almost none of it, an edge-on
  orbit of the same size has all of it. Worked values: 0.48 µas for a
  companion 10 AU deep at 10 pc, 3.9 µas at 20 AU and 5 pc, 1.4 µas for a
  Barnard-like host at 1.8 pc with a 1 AU companion, and ~48 µas for a wide
  1000 AU system at 10 pc.
- against [`framedirection`](@ref) — a *direction*, which no observer displacement
  moves — the target keeps its parallax factor in full, so you get the
  complete parallax ellipse plus the orbit. This is absolute astrometry, and
  it needs no separately-computed parallax factors.

For pulsar timing the Kopeikin coupling is the standard route to `i` and `Ω`;
for astrometry it matters only for nearby, wide, inclined systems, and then it
is not small.

### Which reference for which data

Both are designed uses, so neither is refused — but they answer different
questions, and with no observer argument they give nearly identical numbers,
so the choice is easy to make by accident. Pick from what your instrument
measured:

| your data | reference | what the observer argument adds |
|---|---|---|
| relative astrometry: a companion's offset from its host (imaging, interferometry) | the host body, e.g. `raoff(sol, b, A, obs)` | the differential (Kopeikin) term only |
| a resolved pair's separation from its own barycentre or photocentre | `barycentre(sys)` / `photocentre(sys)` | the differential term only |
| absolute astrometry: a source's position against the sky (Gaia, Hipparcos, HST FGS) | [`framedirection`](@ref) | the full parallax ellipse, plus the orbit |
| a blended catalog source's absolute position | `photocentre(sys)` as the **target**, `framedirection` as the reference | the full parallax ellipse of the blended point |

The rule underneath: name a reference at a **finite distance** when the
measurement is differential, because the observer's displacement moves both
endpoints and cancels; name the **direction** when the measurement is against
the sky, because nothing cancels. If you get it wrong the zero-argument
answers are unchanged and only the observer-aware ones move — by the parallax,
which for a nearby target is orders of magnitude larger than anything else on
this page, so it will not be subtle once you look.

Two requirements, both enforced with an error rather than a silent
degradation. The system needs an `AbsoluteFrame` — `obs_pos` is ICRS, so
placing it relative to the target requires the target's ICRS direction — and
the trajectory must have been solved with `observing_geometry=true`, since the
cheap path stores no per-epoch viewing triad and a µas-level observer coupling
computed on top of it would not be coherent anyway.

`orbitsolve` gains nothing from any of this. There is no ambient observer
state: an observer is visible at every call site that has one, so nothing can
change silently, and one trajectory serves a model containing several
observers (Hipparcos, Gaia, a ground spectrograph) at once. Ephemerides stay
the caller's business — PlanetOrbits takes positions, not spacecraft.

!!! note "Timescales"
    Epochs passed to `orbitsolve` are barycentric (BJD\_TDB-like MJD), and no
    timescale machinery enters PlanetOrbits. A constant-rate rescaling between
    TCB and TDB (1.5e-8) is absorbed exactly by the fitted period, so it costs
    a paragraph rather than code.

## What is not modelled at all

Everything in this section is outside every setting above: `orbitsolve` does
not compute it, and no keyword or argument turns it on. Several of these terms
are *larger* than corrections the flags do gate, so read this section before
assuming a sub-microarcsecond or sub-decimetre-per-second budget is met.

### Observer-side geometry

- **Differential stellar aberration** — 9.9 µas at 100 mas separation, 99 µas
  at 1″. Depends on the observer's velocity, so it is removed by the data
  reduction or belongs in the instrument layer.

- **The annual radial-velocity term.** The same swing in viewing direction
  that produces the Kopeikin coupling reprojects the target's transverse space
  velocity onto the line of sight, giving an annual signal of amplitude
  `v_t · ϖ` — 1.5 cm/s for 30 km/s at 10 pc, but 24 cm/s for a Barnard-like
  90 km/s at 1.8 pc, above the stability floor of current spectrographs. It is
  deliberately absent rather than merely unimplemented: a barycentric-velocity
  correction of Wright & Eastman (2014) grade already removes it given the
  catalog parallax and proper motion, which every RV pipeline applies, so
  modelling it here would double-count. There are no observer-aware `radvel`
  forms for that reason.

### Relativistic terms

- **The body's own surface gravitational redshift** — 636 m/s for a solar
  photosphere. It needs a radius, which `Body` does not carry, and it is
  constant, so it is absorbed by the fitted velocity offset γ. (The
  *companion's* potential at the body, which is not constant, **is**
  modelled — see above.)

- **Shapiro delay** — tens of µs for a solar-mass companion. Irrelevant to RV,
  relevant to edge-on timing.

- **Post-Newtonian corrections to the orbit itself**, such as relativistic
  periastron precession.

!!! note
    Of the usual timing triad, the Rømer terms *are* modelled — the common one
    by `barycentric_lighttime` and the per-body one by `observing_geometry` —
    and so now is the Einstein term. Only Shapiro is missing.


## Cost

150 epochs, one companion, absolute frame, single core. Every combination is
allocation-free.

| `observing_geometry` | `barycentric_lighttime` | `orbitsolve!` |
|---|---|---|
| `true` | `true` | 15.8 µs |
| `true` | `false` | 10.3 µs |
| `false` | `true` | 12.0 µs |
| `false` | `false` | 6.6 µs |

For reference, the Kepler solve itself is 3.8 µs of that, so with both
corrections declined the solver is most of what remains.

## Checking the difference

The honest way to choose is to solve both ways and compare against your own
uncertainties. Doing it one flag at a time also shows that the two are
independent:

```@example prec
using PlanetOrbits
using PlanetOrbits: Body, Orbit, System

A = Body(mass=0.16, name=:A)
b = Body(mass=1.5e-6, name=:b)
orb = Orbit(b, about=A; a=1.0, e=0.2, i=0.9, ω=1.1, Ω=2.2, tp=57000.0)
# Barnard-like: 1.83 pc, ~10.4"/yr — deliberately the hardest case
sys = System((A, b), (orb,); plx=546.5, ra=45.0, dec=20.0,
             pmra=7323.0, pmdec=7323.0, rv=-110e3, ref_epoch=57388.0)

epochs = collect(range(57388.0, 57388.0 + 3652.5, length=25))
full = orbitsolve(sys, epochs)

Δ(t) = maximum(abs(raoff(full[k], b, A) - raoff(t[k], b, A))
               for k in eachindex(epochs)) * 1000   # µas

(geometry_off = Δ(orbitsolve(sys, epochs; observing_geometry=false)),
 lighttime_off = Δ(orbitsolve(sys, epochs; barycentric_lighttime=false)),
 both_off = Δ(orbitsolve(sys, epochs; observing_geometry=false,
                         barycentric_lighttime=false)))
```

Two things to read carefully here. The separation is ~0.5″, so this is a case
where the observing geometry genuinely matters — at a tenth of that separation
it would fall by a factor of ten or more.

And the light-travel figure is measured **at fixed elements**, which is the
pessimistic reading: in an actual fit the constant and linear parts of the
shift are absorbed by `tp` and the period, leaving only the curvature from the
table above. Use the fixed-element number to decide whether to bother
investigating; use the curvature to decide whether it changes your posterior.
