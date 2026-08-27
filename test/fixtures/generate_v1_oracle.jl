# Generates test/fixtures/v1_oracle.jl — the *external oracle* table for
# test/oracles.jl. See the header of that file for what the oracle is and why
# it gates the release.
#
# In one sentence: PlanetOrbits v1 (0.11.x) is an independently written
# implementation of the Campbell-element surface v2 reimplemented, so it can
# referee v2's shared math to roundoff — provided v2 is asked for the *same*
# physical model, which means `observing_geometry=false` (see below).
#
# Run against a v1 checkout:
#     julia --project=<path-to-v1-checkout> generate_v1_oracle.jl
#
# The output file is checked in so CI needs no v1 install. Regenerate only if
# the case list changes — and if a regenerated value *moves*, that is a v1
# behaviour change and wants investigating, not blessing.
#
# ---------------------------------------------------
# Why this is not the same thing as fixtures/v1_reference.jl
#
# `v1_reference.jl` is a small (9-case) regression fixture compared at
# rtol 3e-3 / 3e-2, because it is read through v2's *default* observing
# geometry, which applies four corrections v1 never had. It cannot detect a
# 1e-6 error anywhere.
#
# This table is the opposite trade: a wide parameter sweep, read through
# `observing_geometry=false` — which is documented in src/observe.jl as
# selecting "exactly what v1 ... computed" — and therefore gated near
# roundoff. Two different jobs; neither replaces the other.
#
# ---------------------------------------------------
# Scope, and what is deliberately left out
#
#   * Elliptical orbits only. v1's hyperbolic branch set the velocity
#     semiamplitude to zero (see the "hyperbolic orbits" testset in
#     runtests.jl), so v1 is a *known-wrong* referee there and must not be
#     used as one.
#   * `plx` frames only, never the full absolute frame. v1's absolute-frame
#     path carried two different speeds of light (2.998e8 m/s and a hardcoded
#     2.99792e5 km/s — see the note on `c_light_ms` in src/constants.jl), so
#     it disagrees with v2 at 1e-2 relative on the light-travel-shifted
#     phase. That comparison is the loose `v1_reference.jl` fixture's job.
#   * `pmra`/`pmdec` ARE recorded, but the oracle test gates them loosely and
#     separately: v2's proper motions are the exact derivative of `raoff`,
#     including the depth-rate term v1 dropped. See test/oracles.jl.
# ---------------------------------------------------

using PlanetOrbits

@assert !isdefined(PlanetOrbits, :System) "this script must run against PlanetOrbits v1"
@assert pkgversion(PlanetOrbits) < v"1" "expected PlanetOrbits 0.11.x, got $(pkgversion(PlanetOrbits))"

# Phases of the orbital period at which each case is sampled, relative to tp.
# Includes exact periastron (0.0) and exact apoastron (0.5); reaches several
# periods either side of tp so that the mean-anomaly reduction (which the two
# versions implement differently) is exercised rather than assumed.
const PHASES = [-3.37, -0.5, -0.15, 0.0, 0.02, 0.12, 0.37, 0.5, 0.63, 0.88, 4.73]

# Element values cycled through the sweep. `e` and `i` are enumerated fully
# (they select the numerical branches that matter: the near-circular
# ω-degeneracy, the near-parabolic solver regime, face-on/edge-on/retrograde
# geometry); everything else is cycled mod its length so each combination
# still lands on a different a / ω / Ω / M / plx / mass ratio.
const ES = [0.0, 1e-8, 0.3, 0.7, 0.9, 0.99]
const IS = [0.0, 0.4, π/2, 2.5]           # 2.5 rad > π/2 ⇒ retrograde
const AS = [0.05, 1.0, 8.0, 200.0, 0.4]
const ΩS = [0.0, 1.1, -2.3, 4.9]           # note: this is ω (arg. periastron)
const ΩCAPS = [0.0, 2.2, 5.9]
const MS = [0.08, 1.0, 30.0, 1.1]
const PLXS = [2.0, 24.5, 300.0, 55.0]
const FRACS = [0.001, 0.05, 0.3, 0.5]      # Mp / Mtot for the reflex columns

casedefs = Any[]
let n = 0
    for e in ES, i in IS
        n += 1
        cyc(v) = v[mod1(n, length(v))]
        M = cyc(MS)
        push!(casedefs, (
            name = "sweep-$(lpad(n,2,'0'))",
            params = (a=cyc(AS), e=e, i=i, ω=cyc(ΩS), Ω=cyc(ΩCAPS),
                      tp=58849.0, M=M, plx=cyc(PLXS)),
            Mp = cyc(FRACS) * M,
            epochs = nothing,
        ))
    end
end

# Explicit corner cases, on top of the sweep.
append!(casedefs, Any[
    # tp ~16000 periods before the observations: the mean anomaly is reduced
    # from a huge argument, which is the one place the two implementations'
    # different range reductions can diverge catastrophically rather than
    # gracefully. (v2 reduces with `rem2pi`/`vrem2pi`; v1 does not.)
    (name = "distant-tp-short-P",
     params = (a=0.05, e=0.4, i=0.6, ω=1.1, Ω=2.2, tp=0.0, M=1.0, plx=50.0),
     Mp = 0.02, epochs = collect(range(59000.0, 60000.0, length=11))),
    # Same idea with tp *after* the epochs (negative mean anomaly), and an
    # eccentricity high enough that the solver's starter matters.
    (name = "distant-tp-negative",
     params = (a=0.12, e=0.85, i=1.9, ω=-0.7, Ω=3.3, tp=90000.0, M=0.5, plx=120.0),
     Mp = 0.1, epochs = collect(range(58000.0, 58400.0, length=11))),
    # Longest period in the set: the observations cover ~0.5% of one orbit,
    # so the elements are barely constrained and every term is small.
    (name = "long-period",
     params = (a=500.0, e=0.15, i=1.05, ω=2.4, Ω=0.3, tp=58849.0, M=2.0, plx=8.0),
     Mp = 0.4, epochs = collect(range(58000.0, 60000.0, length=11))),
    # Near-parabolic: e = 0.999 is inside the elliptical branch but at the
    # eccentricity where Markley's starter is worst-conditioned.
    (name = "near-parabolic",
     params = (a=30.0, e=0.999, i=0.9, ω=0.2, Ω=1.7, tp=58849.0, M=1.4, plx=15.0),
     Mp = 0.2, epochs = nothing),
    # Exactly face-on and exactly circular at once: z ≡ 0, ω is degenerate,
    # and several products in both implementations are identically zero.
    (name = "circular-face-on",
     params = (a=3.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=58849.0, M=1.0, plx=100.0),
     Mp = 0.5, epochs = nothing),
    # Exactly edge-on: the RV amplitude is maximal and decoff-vs-posz mixing
    # is at its extreme.
    (name = "edge-on-ecc",
     params = (a=0.3, e=0.92, i=π/2, ω=-1.4, Ω=5.1, tp=58849.0, M=0.25, plx=210.0),
     Mp = 0.05, epochs = nothing),
])

results = Any[]
for c in casedefs
    o = orbit(; c.params...)
    P = period(o)
    epochs = isnothing(c.epochs) ? c.params.tp .+ P .* PHASES : c.epochs
    sols = [orbitsolve(o, t) for t in epochs]
    Mp = c.Mp
    data = Dict{Symbol,Any}(
        :posx => [posx(s) for s in sols],
        :posy => [posy(s) for s in sols],
        :posz => [posz(s) for s in sols],
        :velx => [velx(s) for s in sols],
        :vely => [vely(s) for s in sols],
        :velz => [velz(s) for s in sols],
        :radvel => [radvel(s) for s in sols],
        :raoff => [raoff(s) for s in sols],
        :decoff => [decoff(s) for s in sols],
        :pmra => [pmra(s) for s in sols],
        :pmdec => [pmdec(s) for s in sols],
        :projectedseparation => [projectedseparation(s) for s in sols],
        :posangle => [posangle(s) for s in sols],
        # Reflex of the primary about the system barycentre, via v1's
        # companion-mass overloads.
        :raoff_reflex => [raoff(s, Mp) for s in sols],
        :decoff_reflex => [decoff(s, Mp) for s in sols],
        :pmra_reflex => [pmra(s, Mp) for s in sols],
        :pmdec_reflex => [pmdec(s, Mp) for s in sols],
        :radvel_reflex => [radvel(s, Mp) for s in sols],
    )
    push!(results, (; c.name, c.params, Mp, period=P, epochs,
        data=NamedTuple(sort(collect(data), by=first))))
    println("generated: ", c.name, "  P=", P)
end

out = joinpath(@__DIR__, "v1_oracle.jl")
open(out, "w") do io
    println(io, "# AUTO-GENERATED by generate_v1_oracle.jl against PlanetOrbits v1 — do not edit.")
    println(io, "# See test/oracles.jl. $(length(results)) cases × $(length(first(results).epochs)) epochs.")
    print(io, "const V1_ORACLE = ")
    show(io, results)
    println(io)
end
println("wrote ", out, " (", round(filesize(out)/1024, digits=1), " KiB, ", length(results), " cases)")
