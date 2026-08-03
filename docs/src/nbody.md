# N-body integration

`KeplerianApprox` — the default — puts every hierarchy row on an independent
Keplerian orbit. That is exact for two bodies and is the classic approximation
for hierarchical systems, but it cannot express anything that arises from the
bodies actually pulling on each other: transit-timing variations, resonant
interactions, moon–planet–star coupling, secular exchange of eccentricity.

`AHL21` integrates the real N-body problem. It is the symplectic map of Agol,
Hernandez & Langford (2021), and it targets exactly the same trajectory output,
so **every observable works unchanged**.

!!! tip "Why there is no closed form"
    Two bodies have an exact analytic solution; three or more do not, so the
    trajectory has to be stepped forward numerically. *Orbital Mechanics &
    Astrodynamics* sets up the equations of motion and shows why in
    [Many-Body Problems](https://orbital-mechanics.space/the-n-body-problem/many-body-problems.html),
    building from
    [the two-body case](https://orbital-mechanics.space/the-n-body-problem/two-body-inertial-motion.html).

```@example nb
using PlanetOrbits
import PlanetOrbits as PO

A = PO.Body(mass=1.071, name=:A)
b = PO.Body(mass=1.32e-5, name=:b)
c = PO.Body(mass=2.42e-5, name=:c)

sys = PO.System((A, b, c), (
    PO.Orbit(b, about=A;      a=0.1153, e=0.022, i=π/2,      ω=0.4, Ω=0.0,  tp=60000.0),
    PO.Orbit(c, about=(A, b); a=0.1283, e=0.016, i=π/2+0.01, ω=1.3, Ω=0.05, tp=60002.0),
); plx=50.0)

epochs = collect(range(59980.0, 60060.0, length=12))
traj = orbitsolve(sys, epochs; method=AHL21(h=0.65, t0=60000.0))
raoff(traj[1], :b, :A)
```

Switching propagators costs one keyword. Compare against the Keplerian
approximation on the same system:

```@example nb
tk = orbitsolve(sys, epochs)
maximum(abs(raoff(tk[k], :c, :A) - raoff(traj[k], :c, :A)) for k in eachindex(epochs))
```

## The timestep

`AHL21` uses a **fixed** timestep `h` [days]. There is no adaptivity, and the
cost is one map evaluation per `h` of timespan covered.

Guidance is `h ≲ P_min/20` for the shortest period in the system; beyond that a
warning is emitted (once) and the error grows as `h²`. That warning comes from
the allocating `orbitsolve` only — `orbitsolve!` stays logging-free so it can
remain allocation-free.

This is the honest cost model, and it is worth internalising before putting a
tight moon in a model:

| N | cost per step |
|---|---|
| 2 | ≈ 0.5 µs |
| 3 | ≈ 1.4 µs |
| 4 | ≈ 2.9 µs |
| 8 | ≈ 13 µs |

Per-evaluation cost is roughly `(timespan / h) × step cost`. For the
Kepler-36-like system above — 12 epochs, 80 days, `h = 0.65` d — that is about
199 µs, against 2.5 µs for `KeplerianApprox`. Compact systems are the expensive
regime; wide-orbit imaging systems with periods of years need proportionally
few steps.

Both propagators are allocation-free.

## The osculating epoch `t0`

Under `AHL21` the orbital elements are **osculating elements at `t0`**, not
constants of the motion. Elements only set the initial conditions; after that
the integration takes over.

`t0` defaults to the frame's `ref_epoch` for systems with a full absolute
frame, and must be given explicitly otherwise:

```@example nb
try
    orbitsolve(sys, epochs; method=AHL21(h=0.65))
catch err
    println(sprint(showerror, err))
end
```

Epochs before `t0` are reached by marching the velocity-reversed state forward,
and output epochs are produced by a throwaway partial step from the last saved
state — so the main integration stays on the fixed-step symplectic backbone
rather than taking odd fractional steps.

!!! warning "Posteriors are not element-for-element comparable"
    Because the same element values mean different things under the two
    propagators — constants of the motion versus osculating at `t0` — fitting
    with `KeplerianApprox` and refining with `AHL21` does **not** give directly
    comparable posteriors. Re-evaluate, do not translate.

## Convention independence

Under `AHL21` the orbits only build initial conditions, so any invertible set
describing the same configuration integrates identically. Unlike
`KeplerianApprox`, the Jacobi-vs-astrocentric choice is pure bookkeeping here.

## What is exact

For a two-body system the map is exact and the two propagators agree to
roundoff. Energy shows the bounded `h²` oscillation characteristic of a
symplectic integrator — reached within the first few periods and never
exceeded — rather than a secular drift, and angular momentum is conserved to
roundoff. Zero-mass test particles are well-defined: they feel the massive
bodies without disturbing them.

That energy and angular momentum are conserved at all is the property being
tested here; both are derived in
[Constants of Orbital Motion](https://orbital-mechanics.space/constants-of-orbital-motion/energy-is-conserved-in-orbital-motion.html).
A non-symplectic integrator typically drifts in energy secularly, which over a
long baseline shows up as a spurious change in period.

## Differentiability

ForwardDiff `Dual`s flow through the integrator itself — construction, every
sub-step, and the observables. NbodyGradient's analytic Jacobian engine is
deliberately *not* merged; forward-mode duals cover the parameter counts these
models actually have, and the upstream authors flag bugs in its `dqdt` path.

```@example nb
import ForwardDiff
f = m -> begin
    s = PO.System((PO.Body(mass=m[1], name=:A), PO.Body(mass=m[2], name=:b), PO.Body(mass=m[3], name=:c)),
        (PO.Orbit(PO.Body(mass=m[2], name=:b), about=PO.Body(mass=m[1], name=:A);
                  a=0.1153, e=0.022, i=π/2, ω=0.4, Ω=0.0, tp=60000.0),
         PO.Orbit(PO.Body(mass=m[3], name=:c),
                  about=(PO.Body(mass=m[1], name=:A), PO.Body(mass=m[2], name=:b));
                  a=0.1283, e=0.016, i=π/2+0.01, ω=1.3, Ω=0.05, tp=60002.0)); plx=50.0)
    t = orbitsolve(s, epochs; method=AHL21(h=0.65, t0=60000.0))
    sum(raoff(x, :c, :A) for x in t)
end
ForwardDiff.gradient(f, [1.071, 1.32e-5, 2.42e-5])
```

## Citing

Please cite Agol, Hernandez & Langford (2021), MNRAS 507, 1582
([arXiv:2106.02188](https://arxiv.org/abs/2106.02188)) when publishing results
computed with `AHL21`.
