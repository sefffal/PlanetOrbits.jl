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
epoch, per-body (differential) light-travel time, line-of-sight projection, and
a per-body rather than per-system AU→mas scale. All of them scale with the
angular excursion actually observed, ρ:

| correction | size |
|---|---|
| viewing-direction rotation | 4.85e-3 · ρ[mas] · μ[″/yr] · T[yr] µas |
| differential (per-body) light-travel time | 0.099 · ρ[mas] · √(M/a[AU]) µas |
| depth scaling | 4.85e-6 · ρ[mas]² µas |

For **absolute** astrometry ρ is the photocentre reflex, not the relative
orbit. A Jupiter analog at 10 pc gives ρ = 0.475 mas and therefore 0.005 /
0.021 / 1e-6 µas — four orders of magnitude below a 30–100 µas per-epoch
precision. For a resolved system at GRAVITY+ or CRIRES+ precision, ρ is the
full separation and the same formulas give tens of µas.

So the test is a single comparison of `max(ρ) ×` the coefficients above against
your declared measurement precision. It is not a distance cut: distance enters
only through ρ.

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

## Rules of thumb

Turn `observing_geometry` **on** for resolved systems at µas precision —
GRAVITY-class relative astrometry, or anything where the separation reaches
tens of mas and the quoted uncertainty is below ~10 µas.

Turn `barycentric_lighttime` **on** for nearby, high-proper-motion systems, and
for any transit-timing or other absolute-timing data.

Turn both **off** for a large survey of distant unresolved systems at
Gaia-class precision, which is where the savings actually matter.

When in doubt, leave them on and measure whether it costs you anything — the
defaults are the accurate ones.

## What is not modelled at all

The flags above choose between two *source-side* models. Everything in this
section is outside both of them: `orbitsolve` does not compute it at any
setting, and no keyword turns it on. Several of these terms are *larger* than
corrections the flags do gate, so read this section before assuming a
sub-microarcsecond or sub-decimetre-per-second budget is met.

The organising fact is that **PlanetOrbits places the observer at the
solar-system barycentre**. `System` carries the barycentric direction, distance
and space motion of the target and nothing about where the telescope is; the
epochs you pass to `orbitsolve` are barycentric. Every term below is either a
consequence of the observer not actually being at the SSB, or a
special/general-relativistic term in the source frame.

### Observer-side geometry

- **Differential stellar aberration** — 9.9 µas at 100 mas separation, 99 µas
  at 1″. Depends on the observer's velocity, so it is removed by the data
  reduction or belongs in the instrument layer.

- **Annual–orbital parallax** (Kopeikin 1995). The Earth's ±1 AU excursion
  swings the line of sight to the target by the parallax angle over the year,
  which changes the *projection* of the orbit onto the sky. Barycentric
  corrections applied to timestamps and velocities do not absorb it, because it
  is a coupling between the observer's **position** and the target's own
  three-dimensional structure. Per 1 AU of observer displacement,

  ```
  Δθ ≈ 4.85 · z[AU] / d[pc]²  µas   ≡   4.85e-6 · z[mas] · ϖ[mas]  µas
  ```

  where `z` is the **line-of-sight** separation of the two bodies. Note what
  that is not: unlike every entry in the `observing_geometry` table, it does
  *not* scale with the sky separation ρ — a face-on orbit at 1″ has almost
  none of it, an edge-on orbit of the same size has all of it. Worked values:
  0.48 µas for a companion 10 AU deep at 10 pc, 3.9 µas at 20 AU and 5 pc,
  1.4 µas for a Barnard-like host at 1.8 pc with a 1 AU companion, and ~48 µas
  for a wide 1000 AU system at 10 pc. For pulsar timing the same coupling is
  the standard route to `i` and `Ω`; for astrometry it matters only for nearby,
  wide, inclined systems, and then it is not small.

- **The annual radial-velocity term** from the same coupling: the swing in
  viewing direction reprojects the target's transverse space velocity onto the
  line of sight, giving an annual signal of amplitude `v_t · ϖ` — 1.5 cm/s for
  30 km/s at 10 pc, but 24 cm/s for a Barnard-like 90 km/s at 1.8 pc, which is
  above the stability floor of current spectrographs.

### Relativistic terms in the radial velocity

`radvel` returns the coordinate line-of-sight velocity. It is not an
apparent spectroscopic velocity: the second-order Doppler (time-dilation)
term and gravitational redshift are both absent. Their combined size for
body 1 of a pair is

```
Δv = (1/c) · [ v₁²/2 + G·M₂/r ]
```

with a natural scale of `G·M₂/(a·c)` = **2.96 m/s · (M₂/M⊙) / (a/AU)**.

Which body you ask about changes the answer by three orders of magnitude, so
the two common uses of `radvel` have to be priced separately.

**Stellar reflex —** `radvel(sol, A, barycentre(sys))`. Here `M₂` is the
companion mass and `v₁` the star's reflex, so both terms are tiny, and most of
what remains is constant:

| term | size | observable? |
|---|---|---|
| star's own gravitational redshift | 636 m/s (solar) | constant — absorbed by the systemic offset γ |
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
already fits. The **variation** is not, and for a close-in eccentric companion
it is several m/s — comparable to the precision high-resolution
cross-correlation is now reaching on directly-detected companions. If you are
fitting relative RVs of an eccentric companion, this is the term on this page
most likely to matter to you.

Two further relativistic terms are also absent: the Shapiro delay (tens of µs
for a solar-mass companion — irrelevant to RV, relevant to edge-on timing) and
any post-Newtonian correction to the orbit itself, such as relativistic
periastron precession.

!!! note
    The Rømer-type light-travel terms *are* modelled — the common one by
    `barycentric_lighttime` and the per-body one by `observing_geometry`. It is
    only the Einstein (time-dilation plus redshift) and Shapiro terms of the
    usual timing triad that are missing.

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
