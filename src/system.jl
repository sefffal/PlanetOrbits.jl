# ---------------------------------------------------
# System: bodies + hierarchy + frame
# ---------------------------------------------------

"""
    OrbitDomainError(msg)

The *parameter values* — not the system's structure — put an orbit outside the
domain where it is defined: a non-finite mass, `e == 1` exactly, a semi-major
axis that does not size the conic (`a ≤ 0` or `a` non-finite for an ellipse,
`a == 0` or non-finite for a hyperbola) or a period that is not positive and
finite, a radial or zero-separation Cartesian state.

This is deliberately a distinct type from the `ErrorException`s the structural
checks throw, and the distinction is the useful one for a sampler. A structural
problem (a malformed hierarchy, a body with no orbit) is a property of the
*model*: it is true of every draw, so it deserves to be loud and to stop the
run. A domain problem is a property of *one proposal*: a sampler exploring a
wide prior — a parallel-tempering scheme's hottest chains especially — will
reach parameters that overflow to `Inf` under the unconstrained-space
transform, and the right response is to score that proposal `-Inf` and move on,
not to warn.

So callers that build a system per likelihood evaluation should catch
`OrbitDomainError` quietly and let anything else through.
"""
struct OrbitDomainError <: Exception
    msg::String
end
Base.showerror(io::IO, e::OrbitDomainError) = print(io, e.msg)
@noinline _err_domain(msg::AbstractString) = throw(OrbitDomainError(msg))

# Kept `@noinline`, and out of line of the branch that calls it, so that
# interpolating the offending value costs the construction path nothing (the
# `System` builder is on an allocation-free gate).
@noinline function _err_badsma(a, hyperbolic::Bool)
    hyperbolic && _err_domain(
        "a hyperbolic orbit (e > 1) needs a finite, nonzero semi-major axis; " *
        "got a = $a. The conic has no size there, and the derived constants " *
        "(n = √(GM/|a|³), J) are not finite, so every observable would come " *
        "back `NaN`.")
    _err_domain(
        "an elliptical orbit (e < 1) needs a finite, positive semi-major axis; " *
        "got a = $a. The derived constants (n = √(GM/a³) and J = 2πa/P) are not " *
        "finite there, so every observable would come back `NaN`. (A negative " *
        "`a` is the *hyperbolic* convention and needs e > 1. If this came from a " *
        "prior, note that both ends are reachable under the unconstrained-space " *
        "transform: `a ~ Uniform(0, 100)` returns exactly 0, and a " *
        "bounded-below prior returns `Inf`.)")
end

# Finiteness of a whole element group, on primals, as one branch. Folded with
# `&&` over a tuple rather than `all(isfinite, ...)` so it stays a chain of
# compares on scalars the compiler already has in registers — no iteration
# protocol, and nothing that could allocate on the construction path.
#
# `isfinite` of a `Dual` already tests the value alone (unlike `<`, `iszero`
# and `sign`, which ForwardDiff orders lexicographically), so `_primal` here
# is belt-and-braces rather than load-bearing — but it costs nothing, and
# writing the whole guard family the same way is what keeps the rule
# "classify on the primal" checkable by reading rather than by remembering
# which predicates ForwardDiff happens to define narrowly.
@inline _allfinite() = true
@inline _allfinite(x, xs...) = isfinite(_primal(x)) && _allfinite(xs...)

@noinline _err_badelements(e, i, ω, Ω, tp) = _err_domain(
    "orbital elements must be finite; got e = $e, i = $i, ω = $ω, Ω = $Ω, " *
    "tp = $tp. A non-finite angle survives the reduction as `NaN` " *
    "(`rem2pi(Inf, RoundDown)` is `NaN`) and a non-finite `tp` poisons the " *
    "mean anomaly, so the row would build without complaint and return " *
    "`NaN` for every observable. Unbounded priors reach `±Inf` under the " *
    "unconstrained-space transform, so a hot parallel-tempering chain gets " *
    "here on its own; catching it as a domain error scores the proposal " *
    "`-Inf` instead of poisoning the log-density.")

@noinline _err_badconstants(a, e, M, n, J) = _err_domain(
    "these elements are inside the admissible set but their derived " *
    "constants are not representable: a = $a, e = $e, M = $M give mean " *
    "motion n = $n and velocity factor J = $J. This is the float-range " *
    "floor rather than a statement about the orbit — `a³` underflows to " *
    "zero below a ≈ 1e-107, so `n = √(GM/a³)` overflows even though the " *
    "conic itself is perfectly well defined. Every observable would come " *
    "back `NaN`, so the row is rejected here instead.")

# One hierarchy row: the Keplerian orbit of a binary's exterior barycentre
# about its interior barycentre, with the same derived constants v0.11's
# KepOrbit precomputed. The row's gravitating mass is the total mass of every
# body it binds (the Keplerian approximation of the relative orbit).
struct Row{T<:Number}
    a::T; e::T; i::T; ω::T; Ω::T; tp::T
    M::T          # total mass of all member bodies [M⊙]
    n::T          # mean motion [rad / julian year]
    # √(1−e²) for ellipses; −√(e²−1) for hyperbolae, which is the analytic
    # continuation that makes the shared state kernel produce the correct
    # sin ν once cos E / sin E are read as cosh H / sinh H.
    sqrt1me2::T
    cosi::T; sini::T; cosΩ::T; sinΩ::T; cosω::T; sinω::T
    ecosω::T; esinω::T; cosi_cosΩ::T; cosi_sinΩ::T
    J::T          # velocity semiamplitude factor [AU / julian year]
    hyperbolic::Bool
end

function Row(a, e, i, ω, Ω, tp, M)
    # Same invariants as v0.11's KepOrbit inner constructor
    M = max(M, zero(M))
    # The angle reductions below map a non-finite argument to `NaN` without
    # complaint (`rem2pi(Inf, RoundDown)` is `NaN`), and a `NaN` angle then
    # travels through `sincos` into every observable — the silent-`NaN`
    # failure this whole guard family exists to prevent. `tp` never reaches a
    # reduction at all, but enters the mean anomaly as `n·(t − tp)`, so it has
    # the same standing. `e` is checked here rather than in either conic
    # branch because it is what *selects* the branch: `NaN` classifies
    # elliptical (`NaN >= 1` is false) and `Inf` hyperbolic, and both then
    # produce `NaN` derived constants a branch-local check would have to
    # duplicate. All five are reachable from ordinary priors under the
    # unconstrained-space transform — see the note on `a` below.
    _allfinite(e, i, ω, Ω, tp) || _err_badelements(_primal(e), _primal(i),
        _primal(ω), _primal(Ω), _primal(tp))
    i = rem(i, oftype(i, π), RoundDown)
    Ω = rem2pi(Ω, RoundDown)
    # `_primal` throughout the classification, for the reason spelled out at
    # the `a` guard below: ForwardDiff orders `Dual`s lexicographically, so
    # `Dual(1.0, +∂) >= 1` is `true` while `Dual(1.0, −∂) >= 1` is `false`.
    # Classifying on the perturbation would let one gradient component take
    # the hyperbolic branch and another the elliptical one for the *same*
    # orbit — and, worse, `Dual(1.0, +∂) > 1` is also `true`, so the exactly-
    # parabolic guard below would not fire and the row would be built with
    # sqrt1me2 = −√(e²−1) = −0 and an infinite derivative.
    hyperbolic = _primal(e) >= 1
    # a < 0 is the convention for hyperbolae; there is no valid a > 0 reading,
    # so a positive value is taken as |a|. `show` prints the stored (negative)
    # value, so the normalization is visible rather than silent.
    hyperbolic && _primal(a) > 0 && (a = -a)
    sini, cosi = sincos(i)
    sinω, cosω = sincos(ω)
    sinΩ, cosΩ = sincos(Ω)
    if hyperbolic
        _primal(e) > 1 ||
            _err_domain("parabolic orbits (e == 1) are not supported: the " *
                        "elements are degenerate there (a → ∞). Perturb e " *
                        "slightly, or use Cartesian initial conditions.")
        # After the normalization above the only values left here are a < 0 —
        # unless a is 0, NaN or −Inf, none of which size a conic. Same reason
        # as the elliptical guard below.
        -Inf < _primal(a) < 0 || _err_badsma(_primal(a), true)
        # n = √(μ/|a|³) and J = √(μ/p), both per *julian* year to match the
        # elliptical branch (for which these are algebraically identical to
        # 2π/P_yr and (2πa/P_yr)/√(1−e²)). p = a(1−e²) > 0 for a < 0, e > 1.
        μ = GM_sun_au3_julianyr2 * M
        n = √(μ / (-a)^3)
        sqrt1me2 = -√(e^2 - 1)
        J = √(μ / (a * (1 - e^2)))
    else
        # `a` sizes the conic and every derived constant divides by it: at
        # a == 0 the mean motion 2π/P = √(GM/a³) is `Inf` and J is `NaN`, and
        # nothing downstream errors — the observables are simply `NaN`, so a
        # likelihood gets a `NaN` log-density instead of a rejected proposal.
        # That is reachable from an ordinary prior: `a ~ Uniform(0, 100)`
        # inverse-transforms to *exactly* `0.0` once the sampler proposes a
        # sufficiently negative unconstrained coordinate, and the prior density
        # there is finite, so nothing else rejects the draw. The other end is
        # the same story as the overflowing mass of e3ab8f0 — the inverse of a
        # bounded-below transform is `lower + exp(y)`, so `a` reaches `Inf`
        # too, where J = 2πa/P is 0·∞ = `NaN`. So the admissible set is the
        # open interval, written as a comparison rather than as `isfinite` so
        # that `NaN` fails it as well.
        #
        # `_primal` is load-bearing, not decoration. ForwardDiff orders `Dual`s
        # *lexicographically*: `0 < Dual(0.0, 1.0)` is `true`, because the
        # values tie and the partial breaks it. A gradient-based sampler is
        # exactly where a == 0 arrives, and comparing the `Dual` directly would
        # let the one call that matters through while the primal call ahead of
        # it was rejected. Same idiom as `_validate_masses`.
        0 < _primal(a) < Inf || _err_badsma(_primal(a), false)
        period_days = √(a^3 / M) * kepler_year_to_julian_day_conversion_factor
        period_yrs = period_days / year2day_julian
        n = 2π / period_yrs
        sqrt1me2 = √(1 - e^2)
        J = (2π * a / period_yrs) / sqrt1me2
    end
    # The guards above admit `a` on the open interval, but admissible `a` is
    # not the same as representable `n`: `a³` underflows to zero below
    # a ≈ 1e-107 (and `μ/(-a)³` overflows at the same place on the hyperbolic
    # branch), so `n = √(μ/a³)` is `Inf`, the mean anomaly `n·(t − tp)` is
    # `Inf`, and every angle reduction of it — exact Payne–Hanek included —
    # is `NaN`. The ellipse is perfectly well defined there in exact
    # arithmetic; it is the *derived constants* that leave float range, which
    # is why this is a separate check on them rather than a tighter interval
    # on `a`. (Note it is emphatically not a fix for merely *large* `n`: a
    # finite mean anomaly of 1e60 is reduced correctly and must be — see
    # `VREM2PI_MAX` in kepsolve-simd.jl.) `sqrt1me2` joins them because it is
    # the third quantity the state kernel divides and multiplies by.
    #
    # Gated on `isfinite(M)`, which is not redundant. A row whose *total* mass
    # overflows to `Inf` — two bodies of 1e308 M⊙, say — is deliberately
    # admissible: `_normalize_masses` rescales so the A⁻¹ barycentre weights
    # stay exact, and e3ab8f0 made that a tested contract on the grounds that
    # an overflowing sum is not a malformed hierarchy. Its `n` is `Inf` too, so
    # an ungated check here would reject it and overturn that decision as a
    # side effect. This guard is about the *float-range floor under finite
    # inputs*; a non-finite total mass is the mass contract's business.
    isfinite(_primal(M)) && !_allfinite(n, J, sqrt1me2) &&
        _err_badconstants(_primal(a), _primal(e), _primal(M),
                          _primal(n), _primal(J))
    T = typeof(a)
    return Row{T}(a, e, i, ω, Ω, tp, M, n, sqrt1me2,
        cosi, sini, cosΩ, sinΩ, cosω, sinω,
        e * cosω, e * sinω, cosi * cosΩ, cosi * sinΩ, J, hyperbolic)
end

"""
    System(bodies, orbits; plx=…, ra=…, dec=…, pmra=…, pmdec=…, rv=…, ref_epoch=…)
    System(orbits; …)

A hierarchical system. `bodies` is a tuple of `Body` values defining the body
order (and hence the indices of `bodies(sys)`); `orbits` is a tuple of
`Orbit`s, one per degree of freedom — exactly `length(bodies) - 1` of them.
Given only `orbits`, the bodies are collected in order of first appearance.

A single body with no orbits is legal: `System((A,), ())`. The lone body *is*
the system barycentre, so its motion is purely the frame's (position, proper
motion, parallax, and the observing-geometry corrections); all displacement
observables of the body relative to the barycentre are identically zero. This
is the natural model for fitting the absolute astrometry of an isolated star.

    A = Body(mass=1.1, name=:A)
    b = Body(mass=8mjup, name=:b)
    c = Body(mass=2mjup, name=:c)
    System((A, b, c), (Orbit(b, about=A;      a=2.5),
                       Orbit(c, about=(A, b); a=8.0)))

The topology is whatever the orbits say — Jacobi chains, astrocentric sets,
moons, 2+2 quadruples, and mixtures of these are all expressible, and the
convention is never inferred. See `Jacobi` and `Astrocentric` for the two
standard chains, and `show(sys)` for the convention actually resolved.

The keyword arguments define the system's frame, attached to the system
barycentre:

  - none        → physical-unit observables only,
  - `plx` [mas] → angular observables in mas,
  - all of them → full absolute frame (see `AbsoluteFrame`).

`System` values are cheap, immutable, `isbits` structures rebuilt every
sample. Use `bodies(sys)` for the persistent per-body references.
"""
struct System{NB,NR,T<:Number,F<:AbstractFrame,Names,FL<:NamedTuple,L}
    masses::SVector{NB,T}
    rows::NTuple{NR,Row{T}}
    specs::NTuple{NR,RowSpec{NB}}
    # A⁻¹ action: absolute barycentric state of body j = Σ_k Ainv[j,k] × (row-k
    # relative state). (The system-barycentre term is identically zero in the
    # propagation frame.) Rebuilt per sample from the masses — this is how AD
    # flows through them.
    Ainv::SMatrix{NB,NR,T,L}
    # Per-band fluxes as body-vectors, for photocentre weights.
    fluxes::FL
    frame::F
end

function System(bodyspec, orbits; kwargs...)
    bods = _astuple(bodyspec)
    orbs = _astuple(orbits)
    frame = _make_frame(; kwargs...)
    NB = length(bods)
    NR = length(orbs)
    names = _bodynames(bods)
    _allunique_tuple(names) || error("body names must be unique, got $names")

    specs = map(o -> RowSpec(names, o), orbs)
    _validate_topology(specs, names, Val(NB), Val(NR))

    # Promote every scalar to a common float type. Body fluxes are included:
    # a `flux_<band>` variable may be sampled, or derived from a sampled one,
    # so under a gradient-based sampler it arrives as a `ForwardDiff.Dual`. If
    # the promotion ignored fluxes, `_collect_fluxes` would then try to narrow
    # a Dual back to Float64 and the system would fail to build — which blocked
    # sampled contrasts for PhotometryObs, ImageObs and interferometry alike.
    T = promote_type(map(l -> typeof(l.mass), bods)...,
        map(o -> typeof(o.a), orbs)...,
        map(_flux_scalar_type, bods)...,
        _frame_scalar_type(frame))
    masses = SVector{NB,T}(map(l -> l.mass, bods))
    # Before anything reads them. A non-finite mass makes every barycentre
    # downstream `NaN`, and the failure would otherwise surface as the
    # *singular-hierarchy* error below — which blames the topology for what is
    # purely a value problem, and is exactly the misdiagnosis a sampler
    # exploring an unconstrained parameter space runs into.
    _validate_masses(masses, names)

    rows = ntuple(Val(NR)) do k
        o = orbs[k]
        # Normally the row's gravitating mass is the total mass of every body
        # it binds, read from the *system's* masses so a stale Body value in
        # an Orbit spec can never disagree. `M=` overrides it verbatim.
        M = o.Moverride ? T(o.M) :
            _maskmass(masses, specs[k].int) + _maskmass(masses, specs[k].ext)
        Row(T(o.a), T(o.e), T(o.i), T(o.ω), T(o.Ω), T(o.tp), M)
    end

    Ainv = _build_ainv(masses, specs)
    # Backstop for a redundancy the per-row structural rules cannot see — a
    # cycle spread over three or more rows, say — which shows up only as a
    # singular H. Reaching it now really does mean the topology: masses are
    # validated above, and `_build_H` is scaled so that no admissible set of
    # them can overflow a weight.
    all(isfinite, Ainv) || _err_singular(names, specs)

    fluxes = _collect_fluxes(bods, Val(NB), T)
    frameT = _convert_frame(T, frame)
    return System{NB,NR,T,typeof(frameT),names,typeof(fluxes),NB * NR}(
        masses, rows, specs, Ainv, fluxes, frameT)
end

"""
    reframe(sys, frame)
    reframe(sys; plx=…, ra=…, dec=…, pmra=…, pmdec=…, rv=…, ref_epoch=…)

`sys` with its frame replaced and **nothing else changed**. The keyword form
takes exactly the frame arguments `System` does, so
`reframe(sys; plx=p)` and `reframe(sys)` (→ `NoFrame`) are spelled the same
way as construction.

Masses, rows, row specs, `A⁻¹` and fluxes are all frame-independent by
construction — the frame enters only in `frame_pass!`/`observe_pass!`, never
in the propagation — so they are carried over verbatim rather than rebuilt.
That is what makes this exact rather than merely equivalent: no element is
recomputed, so there is no room for a recomputation to differ in the last bit.

The motivating use is a *two-stage* frame: the barycentric `AbsoluteFrame`
inputs are not known until the model's own bodies have been solved, but
body–barycentre kinematics are fully defined at `Parallax` level (they need a
distance, not a sky position). So a system can be built and solved at
`Parallax` level, that solution used to derive the absolute frame, and the
frame then installed here. Octofitter's anchored-frame parameterization is
exactly this: sample the *anchor source's* observed catalog quantities, and
subtract the model's own motion of the anchor body about the system
barycentre to get the barycentric ones.

The scalar type only ever widens: `promote_type(T, …)` of the current type and
the new frame's. A frame narrower than the system (a `Float64` frame onto a
`Dual` system) leaves the system `Dual`, where fresh construction would have
given `Float64` — same values, different type. Widening cannot lose a
derivative; narrowing would.

See also [`System`](@ref).
"""
function reframe(sys::System{NB,NR,T,F,Names}, frame::AbstractFrame) where {NB,NR,T,F,Names}
    T2 = promote_type(T, _frame_scalar_type(frame))
    frameT = _convert_frame(T2, frame)
    masses = SVector{NB,T2}(sys.masses)
    rows = map(r -> _convert_row(T2, r), sys.rows)
    Ainv = SMatrix{NB,NR,T2}(sys.Ainv)
    fluxes = map(v -> SVector{NB,T2}(v), sys.fluxes)
    return System{NB,NR,T2,typeof(frameT),Names,typeof(fluxes),NB * NR}(
        masses, rows, sys.specs, Ainv, fluxes, frameT)
end

reframe(sys::System; kwargs...) = reframe(sys, _make_frame(; kwargs...))

# Field-wise, not `Row(a, e, i, ω, Ω, tp, M)`: re-running the outer constructor
# would recompute `n`, `J` and the sincos block at the new precision, and the
# whole point of `reframe` is that nothing about the orbit is recomputed.
@inline _convert_row(::Type{T}, r::Row) where {T} = Row{T}(
    r.a, r.e, r.i, r.ω, r.Ω, r.tp, r.M, r.n, r.sqrt1me2,
    r.cosi, r.sini, r.cosΩ, r.sinΩ, r.cosω, r.sinω,
    r.ecosω, r.esinω, r.cosi_cosΩ, r.cosi_sinΩ, r.J, r.hyperbolic)

# Bodies inferred from the orbits, in order of first appearance
# (interior endpoint before exterior, matching the natural reading order).
System(orbits::Union{Tuple{Vararg{Orbit}},AbstractVector{<:Orbit}}; kwargs...) =
    System(_collect_bodies(_astuple(orbits)), orbits; kwargs...)
System(o::Orbit; kwargs...) = System((o,); kwargs...)

@inline _astuple(t::Tuple) = t
_astuple(v::AbstractVector) = Tuple(v)
@inline _astuple(x) = (x,)

# Written as type-stable tuple recursion, not a loop with a growing
# accumulator: every decision here is a function of the *types* (body names
# are type parameters), so it constant-folds away. A loop version infers as
# Tuple and its instability propagates into `System`'s type, costing ~30 kB
# per likelihood evaluation downstream.
@inline _collect_bodies(orbs::Tuple) = _cb_orbits(orbs, ())
@inline _cb_orbits(::Tuple{}, acc) = acc
@inline function _cb_orbits(orbs::Tuple, acc)
    o = first(orbs)
    acc = _cb_add(_specbodies(o.interior), acc)
    acc = _cb_add(_specbodies(o.exterior), acc)
    return _cb_orbits(Base.tail(orbs), acc)
end
@inline _cb_add(::Tuple{}, acc) = acc
@inline _cb_add(bs::Tuple, acc) = _cb_add(Base.tail(bs), _cb_push(first(bs), acc))
@inline _cb_push(b, acc) = _hasname(acc, _name(b)) ? acc : (acc..., b)
@inline _hasname(::Tuple{}, ::Symbol) = false
@inline _hasname(acc::Tuple, nm::Symbol) =
    _name(first(acc)) === nm || _hasname(Base.tail(acc), nm)

# ---------------------------------------------------
# A⁻¹ from the hierarchy matrix H
#
# H is NB×NB: one row per hierarchy relationship, giving that row's relative
# coordinate as a weighted combination of body positions, plus a final row
# for the system barycentre. Inverting it and dropping the barycentre column
# (identically zero in the propagation frame) gives A⁻¹.
#
# This is one path for *every* convention. The closed form the Jacobi case
# admits is only ~20 ns cheaper on a 22–31 µs workload, and having a second
# path is how you ship an astrocentric system silently evaluated with the
# Jacobi formula (a 0.4 AU error, not a rounding difference).
# ---------------------------------------------------

@inline function _build_ainv(masses::SVector{NB,T},
                             specs::NTuple{NR,RowSpec{NB}}) where {NB,NR,T}
    Hi = inv(_build_H(masses, specs))
    return hcat(ntuple(k -> Hi[:, k], Val(NR))...)
end

# NR == 0 (a single body): there are no relative coordinates, so A⁻¹ is the
# empty NB×0 matrix and the body's absolute barycentric state is identically
# zero — the lone body sits at the system barycentre. `hcat()` of zero columns
# cannot spell that, so it is written explicitly.
@inline _build_ainv(masses::SVector{NB,T}, ::Tuple{}) where {NB,T} =
    SMatrix{NB,0,T,0}()

@inline function _build_H(masses::SVector{NB,T},
                          specs::NTuple{NR,RowSpec{NB}}) where {NB,NR,T}
    m = _rescale_masses(masses)
    hrows = ntuple(Val(NR)) do k
        (_groupweights(m, specs[k].ext) .- _groupweights(m, specs[k].int))'
    end
    # The system-barycentre row is `_groupweights` over *every* body, so the
    # massless limit (a system of test particles → their geometric centre) is
    # the one convention, written once. A literal `masses ./ sum(masses)` here
    # instead gave 0/0 = NaN for that case, and the NaN then read as a singular
    # hierarchy.
    return vcat(hrows..., _groupweights(m, _allmask(Val(NB)))')
end

@inline _allmask(::Val{NB}) where {NB} = SVector{NB,Bool}(ntuple(_ -> true, Val(NB)))

# Every row of H is a *ratio* of masses, so it is homogeneous of degree zero:
# scaling all the masses together cannot change H. Doing that scaling
# explicitly is what makes the construction unconditionally well-conditioned —
# without it, `sum(masses)` overflows for masses ≳1e308/NB and every weight
# collapses to zero, which presents as a singular hierarchy for a perfectly
# ordinary two-body system.
#
# The scale is the *power of two* bracketing the largest mass, so the division
# is exact and the resulting weights are bit-identical to the unscaled ones for
# every input that did not over/underflow to begin with. It is read from the
# primal, so under AD it is a locally constant factor — which is the exact
# derivative too, degree-zero homogeneity again.
#
# Measured at +11 ns per `System` construction (29 → 40 ns for `_build_ainv`,
# NB=3), still allocation-free. Done unconditionally rather than behind an
# `isfinite(sum(masses))` guard: the same reasoning that rejected a separate
# closed-form Jacobi path for being ~20 ns cheaper on a 22–31 µs workload
# applies here, and one path is one path.
@inline function _rescale_masses(masses::SVector{NB,T}) where {NB,T}
    mx = reduce(max, ntuple(j -> abs(_primal(masses[j])), Val(NB)))
    # Anything the scaling cannot help — all-zero, non-finite (`max` propagates
    # the NaN), or a non-float element type — is passed through untouched, so
    # this can only ever widen the set of inputs that build.
    (mx isa AbstractFloat && !iszero(mx) && isfinite(mx)) || return masses
    return masses ./ ldexp(one(mx), exponent(mx))
end

# Normalized barycentre weights of one member set. A *massless* set has no
# mass-weighted barycentre (0/0); its limit is the members' geometric centre,
# which for a single test particle is simply its own position. Without this,
# zero-mass bodies — the `n_planets`-prior pattern — produce NaN states.
@inline function _groupweights(masses::SVector{NB,T},
                               mask::SVector{NB,Bool}) where {NB,T}
    w = SVector{NB,T}(ntuple(j -> mask[j] ? masses[j] : zero(T), Val(NB)))
    tot = sum(w)
    # Test the *primal* value: a differentiated mass of 0.0 is a Dual whose
    # value is zero but whose partials are not, and `iszero` on that is false
    # — which would send a test particle down the 0/0 path under AD only.
    iszero(_primal(tot)) || return w ./ tot
    # A single-member massless group weights to 1 regardless of its mass, so
    # zero partials here is exact for the case that matters (test particles).
    ind = SVector{NB,T}(ntuple(j -> mask[j] ? one(T) : zero(T), Val(NB)))
    return ind ./ sum(ind)
end

# Primal value beneath any nesting of ForwardDiff Duals.
@inline _primal(x::Real) = x
@inline _primal(d::Dual) = _primal(value(d))

@inline _maskmass(masses::SVector{NB,T}, mask::SVector{NB,Bool}) where {NB,T} =
    sum(SVector{NB,T}(ntuple(j -> mask[j] ? masses[j] : zero(T), Val(NB))))

# ---------------------------------------------------
# Topology validation
#
# Structural only — a function of the masks, not the masses, so it says the
# same thing for every sample. Messages name the offending row.
# ---------------------------------------------------

function _validate_topology(specs::NTuple{NR,RowSpec{NB}}, names,
                            ::Val{NB}, ::Val{NR}) where {NB,NR}
    NB >= 1 || _err_nobodies()
    NR == NB - 1 || _err_wrongcount(NB, NR)
    # A single body needs no orbit: it *is* the system barycentre, and its
    # motion is purely the frame's (position, proper motion, parallax). Every
    # per-row check below is vacuous at NR == 0, including the unused-body one
    # — whose premise ("position undetermined") is false here, since the count
    # invariant already pins the lone body to the barycentre.
    NR == 0 && return nothing
    for k in 1:NR
        s = specs[k]
        any(s.int) || _err_emptyside(k, "interior (`about=`)")
        any(s.ext) || _err_emptyside(k, "exterior")
        for j in 1:NB
            (s.int[j] && s.ext[j]) && _err_bothsides(k, names[j])
        end
    end
    for k in 1:NR, l in (k+1):NR
        (specs[k].int == specs[l].int && specs[k].ext == specs[l].ext) &&
            _err_duplicate(k, l, names, specs[k])
        # the same relationship written backwards is equally redundant
        (specs[k].int == specs[l].ext && specs[k].ext == specs[l].int) &&
            _err_reversed(k, l, names, specs[k])
    end
    for j in 1:NB
        any(k -> specs[k].int[j] || specs[k].ext[j], 1:NR) ||
            _err_unused(names[j])
    end
    return nothing
end

# ---------------------------------------------------
# Mass validation
#
# Value-dependent, unlike the topology checks, so it throws `OrbitDomainError`
# and a sampler can treat it as a rejected proposal rather than a model defect.
# ---------------------------------------------------

@inline function _validate_masses(masses::SVector{NB,T}, names) where {NB,T}
    # A loop rather than `all(isfinite, masses)` so the message can name the
    # body; it costs the same, unrolls, and only the failing branch is called.
    # Tested on the primal — the value is what makes the barycentre `NaN`, and
    # a `Dual` with finite value and `NaN` partials is a different complaint.
    # This is the case worth a check on every construction: an `Inf` mass is
    # precisely what an unconstrained-space proposal produces once the inverse
    # of a bounded-below transform (`lower + exp(y)`) overflows.
    @inbounds for j in 1:NB
        isfinite(_primal(masses[j])) || _err_nonfinitemass(names[j], _primal(masses[j]))
    end
    return nothing
end

@noinline _err_nonfinitemass(nm::Symbol, m) = _err_domain(
    "body :$nm has a non-finite mass ($m). Every barycentre in the system is " *
    "a mass-weighted average, so no body's position is defined. If this came " *
    "from a sampler, the proposal is simply out of bounds — score it -Inf. " *
    "If it came from a model, check for a parameter that overflows (a " *
    "log-scaled or bounded-below variable pushed far enough) or a division " *
    "by zero upstream of the mass.")

_setnames(names, mask) = Tuple(names[j] for j in eachindex(names) if mask[j])

@noinline _err_nobodies() = error(
    "a system needs at least one body")
@noinline _err_wrongcount(NB, NR) = error(
    "a system of $NB bodies needs exactly $(NB - 1) orbits to determine every " *
    "body's motion; got $NR. (Each orbit supplies one relative coordinate; the " *
    "remaining degree of freedom is the system barycentre.)")
@noinline _err_emptyside(k, side) = error(
    "orbit $k has an empty $side endpoint")
@noinline _err_bothsides(k, nm) = error(
    "orbit $k lists body :$nm on both sides — a body cannot orbit a reference " *
    "that includes itself")
@noinline _err_duplicate(k, l, names, s) = error(
    "orbits $k and $l are the same relationship ($(_setnames(names, s.ext)) " *
    "about $(_setnames(names, s.int))), so they do not independently determine " *
    "the system. Did you mean a different `about=` for one of them?")
@noinline _err_reversed(k, l, names, s) = error(
    "orbits $k and $l are the same relationship written in opposite directions " *
    "($(_setnames(names, s.ext)) about $(_setnames(names, s.int)), and its " *
    "reverse), so they do not independently determine the system.")
@noinline _err_unused(nm) = error(
    "body :$nm does not appear in any orbit, so its position is undetermined. " *
    "Give it an orbit, or drop it from the body list.")

@noinline _err_singular(names, specs) = error(
    "these orbits do not independently determine every body's position — the " *
    "hierarchy is circular or redundant. Rows, as (exterior about interior):\n" *
    join(("  $k: $(_setnames(names, s.ext)) about $(_setnames(names, s.int))"
          for (k, s) in enumerate(specs)), "\n"))

# ---------------------------------------------------

_frame_scalar_type(::NoFrame) = Bool  # neutral in promote_type

# Scalar type of a body's declared fluxes, for the same promotion. A body with
# no flux contributes `Bool`, which is likewise neutral.
@inline _flux_scalar_type(l) = _flux_scalar_type(l.flux)
@inline _flux_scalar_type(::NamedTuple{(),Tuple{}}) = Bool
@inline _flux_scalar_type(nt::NamedTuple) = promote_type(map(typeof, values(nt))...)
_frame_scalar_type(fr::Parallax) = typeof(fr.plx)
_frame_scalar_type(fr::AbsoluteFrame) = typeof(fr.ra)

_convert_frame(::Type{T}, ::NoFrame) where {T} = NoFrame()
_convert_frame(::Type{T}, fr::Parallax) where {T} = Parallax{T}(fr.plx, fr.cart2angle)
_convert_frame(::Type{T}, fr::AbsoluteFrame) where {T} = AbsoluteFrame{T}(
    fr.ra, fr.dec, fr.plx, fr.pmra, fr.pmdec, fr.rv, fr.ref_epoch,
    fr.distance1, fr.x1, fr.y1, fr.z1, fr.dx, fr.dy, fr.dz,
    SMatrix{3,3,T,9}(fr.M1))

# Default names: positional :body1, :body2, … for unnamed bodies.
@inline function _bodynames(bods::NTuple{NB,Any}) where {NB}
    ntuple(Val(NB)) do j
        nm = _name(bods[j])
        nm === :auto ? Symbol(:body, j) : nm
    end
end

# allocation-free allunique for small tuples (Base's builds a Set)
@inline _allunique_tuple(::Tuple{}) = true
@inline _allunique_tuple(t::Tuple) =
    !(first(t) in Base.tail(t)) && _allunique_tuple(Base.tail(t))

# Union of all bands across bodies; missing bands are zero flux.
function _collect_fluxes(bods, ::Val{NB}, ::Type{T}) where {NB,T}
    protos = merge(map(l -> map(_ -> zero(T), l.flux), bods)...)
    bands = keys(protos)
    vals = map(bands) do band
        SVector{NB,T}(ntuple(j -> T(get(bods[j].flux, band, zero(T))), Val(NB)))
    end
    return NamedTuple{bands}(vals)
end

nbodies(::System{NB}) where {NB} = NB
nrows(::System{NB,NR}) where {NB,NR} = NR
_names(::System{NB,NR,T,F,Names}) where {NB,NR,T,F,Names} = Names

"""
    bodies(sys)

NamedTuple of persistent `BodyRef`s for the bodies of `sys`, keyed by the
names given at construction (`Body(… , name=:b)`), in the order the bodies
were listed. These are the resolved form of what observables accept —
guaranteed cheap in hot loops, and the only handles to *unnamed* bodies.
(Named `Body` values and `Symbol`s resolve to them automatically.)

    (; A, b) = bodies(sys)
    raoff(sol, b, A)
"""
bodies(::System{NB,NR,T,F,Names}) where {NB,NR,T,F,Names} =
    NamedTuple{Names}(ntuple(BodyRef, Val(NB)))

"""
    barycentre(sys)
    barycentre(sys, members...)

The mass-weighted barycentre of the whole system, or of the subsystem
spanned by the given members, as a `WeightedPoint`. Members can be given as
`BodyRef`s, named `Body` values, or `Symbol`s, e.g.
`barycentre(sys, jup, gan)`.
"""
function barycentre(sys::System{NB,NR,T}) where {NB,NR,T}
    return WeightedPoint{NB,T}(sys.masses ./ sum(sys.masses), false)
end
function barycentre(sys::System{NB,NR,T}, members::Union{BodyRef,Body,Symbol}...) where {NB,NR,T}
    return WeightedPoint{NB,T}(_subweights(sys.masses, _membermask(sys, members)), false)
end

# Member set → fixed-width indicator mask, then weights, exactly as row
# membership is stored (see `RowSpec`). The obvious spelling — an index tuple
# plus `j in idxs` inside an `ntuple` — is the same shape §12 records as a
# constant-folding cliff, and it is: at NB=5 under Dual{12} it cost 37 kB per
# call, invisible at NB=2 and invisible for Float64 at any NB.
@inline _membermask(sys::System{NB}, members::Tuple) where {NB} =
    _idxmask(map(m -> _resolve(_names(sys), m).idx, members), Val(NB))

@inline _idxmask(idxs::Tuple, ::Val{NB}) where {NB} =
    SVector{NB,Bool}(ntuple(j -> _inidxs(idxs, j), Val(NB)))
@inline _inidxs(::Tuple{}, ::Int) = false
@inline _inidxs(t::Tuple, j::Int) = first(t) === j || _inidxs(Base.tail(t), j)

# Normalized weights of one member set over a body-vector. No zero-total
# fallback here: `barycentre` keeps its existing behaviour and `photocentre`
# guards the case with an error before calling.
@inline function _subweights(v::SVector{NB,T}, mask::SVector{NB,Bool}) where {NB,T}
    w = SVector{NB,T}(ntuple(j -> mask[j] ? v[j] : zero(T), Val(NB)))
    return w ./ sum(w)
end

"""
    fluxes(sys)
    fluxes(sys, band)

Per-band fluxes of `sys`'s bodies, as they were declared on the `Body`
values: a `NamedTuple` of `SVector{NB}`s keyed by band, or the `SVector{NB}`
for one band. Bodies that declared no flux in a band read zero.

This is the entry point for likelihoods that need to build their *own*
weights — a per-epoch resolution taper, a sampled membership indicator — and
hand the result to [`photocentre`](@ref):

    f = fluxes(sys, :G)
    wp = photocentre(f .* member)     # a WeightedPoint, normalized

Vector order matches `bodies(sys)`.
"""
@inline fluxes(sys::System) = sys.fluxes
@inline fluxes(sys::System, band::Symbol) = _bandfluxes(sys, band)

# Band selection, shared by `fluxes` and `photocentre`. `nothing` means "the
# only band", which is an error when there is more than one.
#
# The lookup goes through the *position* of the band in the (type-level) key
# tuple and indexes the values tuple, rather than `sys.fluxes[band]`. Those
# are the same answer, but `getfield(::NamedTuple, ::Symbol)` needs the
# symbol to be a compile-time constant to stay on the stack, and a keyword
# argument threaded through a varargs method is exactly where that
# constant-propagation gives up — silently, and only under Duals, where the
# escaping SVector is big enough to be heap-allocated (12.9 kB per call at
# NB=5, Dual{12}; 0 bytes for Float64 at any NB). The values tuple is
# homogeneous, so a runtime index into it is neither.
@inline function _bandfluxes(sys::System{NB,NR,T,F,Names,FL}, band::Symbol) where {NB,NR,T,F,Names,FL}
    k = _findname(_bandnames(FL), band, 1)
    k === 0 && _err_noband(band, _bandnames(FL))
    return @inbounds values(sys.fluxes)[k]
end
@inline function _bandfluxes(sys::System, ::Nothing)
    isempty(sys.fluxes) && _err_noflux()
    length(sys.fluxes) == 1 || _err_multiband(keys(sys.fluxes))
    return @inbounds values(sys.fluxes)[1]
end

@inline _bandnames(::Type{<:NamedTuple{Bands}}) where {Bands} = Bands

@noinline _err_noflux() = error(
    "no fluxes defined: give at least one body a flux, e.g. Body(…, flux=(G=1.0,))")
@noinline _err_multiband(bands) = error(
    "multiple bands defined ($(bands)); pass band=…")
@noinline _err_noband(band, bands) = error(
    isempty(bands) ?
    "no body in this system declares a flux, so band :$band is undefined. " *
    "Give the bodies fluxes, e.g. Body(…, flux=(; $band=1.0))" :
    "no band :$band in this system (its bands are $(bands))")
@noinline _err_zeroflux(band) = error(
    "total flux is zero" * (band === nothing ? "" : " in band :$band") *
    "; cannot form a photocentre. (A structural photocentre over bodies that " *
    "are all dark is not a point on the sky — if membership is draw-dependent, " *
    "build the WeightedPoint in the likelihood and guard it there.)")

"""
    photocentre(sys; band=nothing)
    photocentre(sys, members...; band=nothing)
    photocentre(weights::StaticVector)

The flux-weighted photocentre of the system — the point astrometric
instruments observe for blended sources — as a `WeightedPoint`. With more
than one band defined on the system's bodies, pass `band` to select one.

Given `members`, the photocentre is over that subset only: weights
`f_j / Σ_members f_k` for members and zero for everything else. Members can
be given as `BodyRef`s, named `Body` values, or `Symbol`s, exactly as for
[`barycentre`](@ref):

    photocentre(sys, Aa, Ab; band=:G)

A subset whose total flux is zero is an error: a *structural* membership
declaration over bodies that are all dark has no meaning. Membership that
varies per draw or per epoch is the likelihood's business — build the weight
vector there and pass it to the third method, which normalizes it:

    photocentre(fluxes(sys, :G) .* member)

`WeightedPoint`s are `isbits`, so constructing one per epoch is free.
"""
function photocentre(sys::System{NB,NR,T}; band::Union{Nothing,Symbol}=nothing) where {NB,NR,T}
    fl = _bandfluxes(sys, band)
    total = sum(fl)
    # Test the *primal*: a differentiated zero flux is a Dual whose value is
    # zero but whose partials are not, and `iszero` on that is false — which
    # would send the guard's own failure case down the 0/0 path under AD only.
    iszero(_primal(total)) && _err_zeroflux(band)
    return WeightedPoint{NB,T}(fl ./ total, true)
end

# NB the `@inline`: a method that is *both* varargs and keyword-accepting is
# not inlined by default, and the un-inlined call returns its `WeightedPoint`
# through the heap — 12.4 kB per call at NB=5 under Dual{12}, and 0 bytes for
# Float64, so nothing but a wide-Dual gate sees it. `barycentre(sys,
# members...)` next door is varargs *without* keywords and needs no such
# annotation, which is exactly why this was easy to miss.
@inline function photocentre(sys::System{NB,NR,T}, members::Union{BodyRef,Body,Symbol}...;
                             band::Union{Nothing,Symbol}=nothing) where {NB,NR,T}
    fl = _bandfluxes(sys, band)
    mask = _membermask(sys, members)
    iszero(_primal(_masksum(fl, mask))) && _err_zeroflux(band)
    return WeightedPoint{NB,T}(_subweights(fl, mask), true)
end

@inline _masksum(v::SVector{NB,T}, mask::SVector{NB,Bool}) where {NB,T} =
    sum(SVector{NB,T}(ntuple(j -> mask[j] ? v[j] : zero(T), Val(NB))))

# Normalizing constructor for likelihood-built weights (the tier-2/tier-3
# pattern of design/g23h-v2-port.md): effective fluxes in, WeightedPoint out.
function photocentre(w::StaticVector{NB,T}) where {NB,T<:Number}
    total = sum(w)
    iszero(_primal(total)) && _err_zeroweights()
    return WeightedPoint{NB,T}(SVector{NB,T}(w) ./ total, true)
end

@noinline _err_zeroweights() = error(
    "photocentre weights sum to zero; there is no such point. Guard the " *
    "all-dark / no-member case in the likelihood (e.g. fall back to the host " *
    "body) before building the WeightedPoint.")

# ---------------------------------------------------
# Display — the resolved convention, per row.
#
# Mixed conventions are legal, so `show` never claims a single one for the
# system unless every row agrees.
# ---------------------------------------------------

# Laminar (generalized Jacobi): every row merges two disjoint groups already
# formed, ending with a single group. Covers Jacobi chains, 2+2 quadruples,
# and moons spliced into a chain. Astrocentric sets are deliberately *not*
# laminar — {A} is consumed by the first row — and that is exactly the
# distinction that makes them a different model under `KeplerianApprox`.
#
# Display only: mixed conventions are legal, so this classifies, it never
# rejects. Allocating is fine here; it is never called from the solve path.
function _is_laminar(specs::NTuple{NR,RowSpec{NB}}) where {NR,NB}
    groups = SVector{NB,Bool}[SVector{NB,Bool}(ntuple(i -> i == j, Val(NB))) for j in 1:NB]
    for s in specs
        ki = findfirst(==(s.int), groups)
        ke = findfirst(==(s.ext), groups)
        (ki === nothing || ke === nothing || ki == ke) && return false
        merged = s.int .| s.ext
        deleteat!(groups, sort!([ki, ke]))
        push!(groups, merged)
    end
    return length(groups) == 1
end

function _system_convention(specs::NTuple{NR,RowSpec{NB}}) where {NR,NB}
    NR == 0 && return :trivial
    NR == 1 && return :twobody
    # Astrocentric first: a star of single-body rows all sharing one centre.
    if all(k -> count(specs[k].int) == 1 && specs[k].int == specs[1].int, 1:NR)
        return :astrocentric
    end
    _is_laminar(specs) && return :jacobi
    return :mixed
end

_g(x, sig=6) = string(round(float(x); sigdigits=sig))

function Base.show(io::IO, ::MIME"text/plain", sys::System{NB,NR,T}) where {NB,NR,T}
    names = _names(sys)
    conv = _system_convention(sys.specs)
    label = conv === :jacobi ? "Jacobi (hierarchical)" :
            conv === :astrocentric ? "astrocentric" :
            conv === :twobody ? "two-body" :
            conv === :trivial ? "single body (at the barycentre)" :
            conv === :mixed ? "mixed conventions" : "custom"
    println(io, "System{$NB bodies, $NR orbits, $T} — $label")
    print(io, "  frame: ")
    show(io, MIME"text/plain"(), sys.frame)
    println(io)
    println(io, "  bodies:")
    for j in 1:NB
        println(io, "    ", rpad(string(names[j]), 10), " mass = ", _g(sys.masses[j]), " M⊙")
    end
    if NR == 0
        println(io, "  orbits: (none — the body is the system barycentre; ",
            "its motion is the frame's)")
        return nothing
    end
    println(io, "  orbits:")
    for k in 1:NR
        s = sys.specs[k]
        r = sys.rows[k]
        Pd = _period(r)
        # The explicit reference (`:A` vs `barycentre(:A, :b)`) *is* the
        # convention, per row — no separate tag can say more, and a tag would
        # contradict the system label on the innermost row, where Jacobi and
        # astrocentric coincide.
        println(io, "    ", k, ": ", _fmtset(_setnames(names, s.ext)),
            " about ", _fmtset(_setnames(names, s.int)))
        # P in both days and years: `P=` is in DAYS, and imaging users think
        # in years, so a 365× slip yields a plausible wrong answer.
        println(io, "       a = ", rpad(_g(r.a), 12), "AU   P = ", rpad(_g(Pd), 12),
            "d (= ", _g(Pd / year2day_julian), " yr)   e = ", _g(r.e, 4))
        println(io, "       i = ", rpad(_g(r.i), 12), "rad  ω = ", rpad(_g(r.ω), 12),
            "rad  Ω = ", rpad(_g(r.Ω), 12), "rad  tp = ", _g(r.tp, 9), " MJD")
        println(io, "       M = ", _g(r.M, 8), " M⊙",
            _row_moverride(sys, k) ? "   (M= override — compatibility only)" : "")
    end
    return nothing
end

_fmtset(t::Tuple) = length(t) == 1 ? ":$(t[1])" :
    "barycentre(" * join((":$n" for n in t), ", ") * ")"

# The override flag is not carried on Row; recover it by comparing against
# the mass the topology implies.
function _row_moverride(sys::System, k::Int)
    s = sys.specs[k]
    implied = _maskmass(sys.masses, s.int) + _maskmass(sys.masses, s.ext)
    return !(sys.rows[k].M ≈ implied)
end

Base.show(io::IO, sys::System{NB,NR,T}) where {NB,NR,T} =
    print(io, "System{$NB bodies, $NR orbits, $T}")

# `System`, `Body` and `Orbit` are deliberately NOT exported: Octofitter owns
# those three names unqualified, and its versions build these (§5, §8.5).
# Import them explicitly here —
#     import PlanetOrbits as PO;  PO.Body(mass=…)
#     using PlanetOrbits: Body, Orbit, System
# Everything else in the export surface is clash-free; `bodies`, `barycentre`
# and `photocentre` are functions Octofitter extends rather than shadows.
#
# Not exporting is what does the functional work. A `public` declaration would
# also mark them as API, but it is parse-level and 1.11+, and the floor here is
# 1.10 — add the bare keyword when the floor moves for an unrelated reason.
#
# `WeightedPoint` and `fluxes` are exported deliberately: likelihoods build
# their own weight vectors (per-epoch resolution tapers, sampled membership)
# and need to name the type and read the fluxes without reaching into
# `sys.fluxes`. Both were checked against everything Octofitter re-exports
# alongside PlanetOrbits — unlike `Trajectory`, neither clashes.
# `reframe` is exported although `System` is not: the frame *types* stay
# unexported (Octofitter builds frames through the same keyword grammar), but
# the two-stage frame is a modelling pattern users write, not an internal.
export Jacobi, Astrocentric, bodies, barycentre, photocentre, fluxes,
       BodyRef, WeightedPoint, reframe
