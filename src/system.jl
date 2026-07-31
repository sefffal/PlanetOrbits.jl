# ---------------------------------------------------
# System: bodies + hierarchy + frame
# ---------------------------------------------------

# One hierarchy row: the Keplerian orbit of a binary's exterior barycentre
# about its interior barycentre, with the same derived constants v1's
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
    # Same invariants as v1's KepOrbit inner constructor
    M = max(M, zero(M))
    i = rem(i, oftype(i, π), RoundDown)
    Ω = rem2pi(Ω, RoundDown)
    hyperbolic = e >= 1
    # a < 0 is the convention for hyperbolae; there is no valid a > 0 reading,
    # so a positive value is taken as |a|. `show` prints the stored (negative)
    # value, so the normalization is visible rather than silent.
    hyperbolic && a > 0 && (a = -a)
    sini, cosi = sincos(i)
    sinω, cosω = sincos(ω)
    sinΩ, cosΩ = sincos(Ω)
    if hyperbolic
        e > 1 || error("parabolic orbits (e == 1) are not supported: the " *
                       "elements are degenerate there (a → ∞). Perturb e " *
                       "slightly, or use Cartesian initial conditions.")
        # n = √(μ / |a|³) with μ = 4π²M in these units
        n = 2π * √(M / (-a)^3)
        sqrt1me2 = -√(e^2 - 1)
        # Perifocal velocity prefactor √(μ/p), p = a(1−e²) > 0 for a < 0, e > 1.
        # For ellipses this is algebraically identical to (2πa/P)/√(1−e²).
        J = 2π * √(M / (a * (1 - e^2)))
    else
        period_days = √(a^3 / M) * kepler_year_to_julian_day_conversion_factor
        period_yrs = period_days / year2day_julian
        n = 2π / period_yrs
        sqrt1me2 = √(1 - e^2)
        J = (2π * a / period_yrs) / sqrt1me2
    end
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

    # Promote every scalar to a common float type.
    T = promote_type(map(l -> typeof(l.mass), bods)...,
        map(o -> typeof(o.a), orbs)...,
        _frame_scalar_type(frame))
    masses = SVector{NB,T}(map(l -> l.mass, bods))

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
    # Structural rules cannot catch every redundancy (e.g. two rows that are
    # each other's reverse), and those show up as a singular H.
    all(isfinite, Ainv) || _err_singular(names, specs)

    fluxes = _collect_fluxes(bods, Val(NB), T)
    frameT = _convert_frame(T, frame)
    return System{NB,NR,T,typeof(frameT),names,typeof(fluxes),NB * NR}(
        masses, rows, specs, Ainv, fluxes, frameT)
end

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

@inline function _build_H(masses::SVector{NB,T},
                          specs::NTuple{NR,RowSpec{NB}}) where {NB,NR,T}
    hrows = ntuple(Val(NR)) do k
        (_groupweights(masses, specs[k].ext) .- _groupweights(masses, specs[k].int))'
    end
    return vcat(hrows..., (masses ./ sum(masses))')
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
    NR == NB - 1 || _err_wrongcount(NB, NR)
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

_setnames(names, mask) = Tuple(names[j] for j in eachindex(names) if mask[j])

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
_frame_scalar_type(fr::Parallax) = typeof(fr.plx)
_frame_scalar_type(fr::AbsoluteFrame) = typeof(fr.ra)

_convert_frame(::Type{T}, ::NoFrame) where {T} = NoFrame()
_convert_frame(::Type{T}, fr::Parallax) where {T} = Parallax{T}(fr.plx, fr.cart2angle)
_convert_frame(::Type{T}, fr::AbsoluteFrame) where {T} = AbsoluteFrame{T}(
    fr.ra, fr.dec, fr.plx, fr.pmra, fr.pmdec, fr.rv, fr.ref_epoch,
    fr.distance1, fr.x1, fr.y1, fr.z1, fr.dx, fr.dy, fr.dz)

@inline _sum_masses(masses, idxs::Tuple{}) = zero(eltype(masses))
@inline _sum_masses(masses, idxs::Tuple) =
    masses[first(idxs)] + _sum_masses(masses, Base.tail(idxs))

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
    return WeightedPoint{NB,T}(sys.masses ./ sum(sys.masses))
end
function barycentre(sys::System{NB,NR,T}, members::Union{BodyRef,Body,Symbol}...) where {NB,NR,T}
    idxs = map(m -> _resolve(_names(sys), m).idx, members)
    total = _sum_masses(sys.masses, idxs)
    w = ntuple(Val(NB)) do j
        j in idxs ? sys.masses[j] / total : zero(T)
    end
    return WeightedPoint{NB,T}(SVector{NB,T}(w))
end

"""
    photocentre(sys; band=nothing)

The flux-weighted photocentre of the system as a `WeightedPoint` — the point
astrometric instruments observe for blended sources. With more than one band
defined on the system's bodies, pass `band` to select one.
"""
function photocentre(sys::System{NB,NR,T}; band::Union{Nothing,Symbol}=nothing) where {NB,NR,T}
    isempty(sys.fluxes) &&
        error("no fluxes defined: give at least one body a flux, e.g. Body(…, flux=(G=1.0,))")
    fl = if band === nothing
        length(sys.fluxes) == 1 ? only(values(sys.fluxes)) :
            error("multiple bands defined ($(keys(sys.fluxes))); pass band=…")
    else
        sys.fluxes[band]
    end
    total = sum(fl)
    iszero(total) && error("total flux in band is zero; cannot form a photocentre")
    return WeightedPoint{NB,T}(fl ./ total)
end

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
            conv === :mixed ? "mixed conventions" : "custom"
    println(io, "System{$NB bodies, $NR orbits, $T} — $label")
    print(io, "  frame: ")
    show(io, MIME"text/plain"(), sys.frame)
    println(io)
    println(io, "  bodies:")
    for j in 1:NB
        println(io, "    ", rpad(string(names[j]), 10), " mass = ", _g(sys.masses[j]), " M⊙")
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

export System, Body, Orbit, Jacobi, Astrocentric, bodies, barycentre, photocentre, BodyRef
