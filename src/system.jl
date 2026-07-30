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
    sqrt1me2::T
    cosi::T; sini::T; cosΩ::T; sinΩ::T; cosω::T; sinω::T
    ecosω::T; esinω::T; cosi_cosΩ::T; cosi_sinΩ::T
    J::T          # velocity semiamplitude factor [AU / julian year]
end

function Row(a, e, i, ω, Ω, tp, M)
    # Same invariants as v1's KepOrbit inner constructor
    M = max(M, zero(M))
    i = rem(i, oftype(i, π), RoundDown)
    Ω = rem2pi(Ω, RoundDown)
    if e >= 1
        error("hyperbolic orbits (e ≥ 1) are not yet supported in PlanetOrbits v2")
    end
    period_days = √(a^3 / M) * kepler_year_to_julian_day_conversion_factor
    period_yrs = period_days / year2day_julian
    n = 2π / period_yrs
    sqrt1me2 = √(1 - e^2)
    sini, cosi = sincos(i)
    sinω, cosω = sincos(ω)
    sinΩ, cosΩ = sincos(Ω)
    J = (2π * a / period_yrs) / sqrt1me2
    T = typeof(a)
    return Row{T}(a, e, i, ω, Ω, tp, M, n, sqrt1me2,
        cosi, sini, cosΩ, sinΩ, cosω, sinω,
        e * cosω, e * sinω, cosi * cosΩ, cosi * sinΩ, J)
end

"""
    System(root::Union{Body,Binary}; plx=…, ra=…, dec=…, pmra=…, pmdec=…, rv=…, ref_epoch=…)

A hierarchical system of bodies. `root` is a nested tree of `Binary`s over
`Body`s. The keyword arguments define the system's frame, attached to the
system barycentre:

  - none        → physical-unit observables only,
  - `plx` [mas] → angular observables in mas,
  - all of them → full absolute frame (see `AbsoluteFrame`).

`System` values are cheap, immutable, `isbits` structures rebuilt every
sample. Use `bodies(sys)` for the persistent per-body references.
"""
struct System{NB,NR,T<:Number,F<:AbstractFrame,Names,FL<:NamedTuple,L}
    masses::SVector{NB,T}
    rows::NTuple{NR,Row{T}}
    # A⁻¹ action: absolute barycentric state of body j = Σ_k Ainv[j,k] × (row-k
    # relative state). (The system-barycentre term is identically zero in the
    # propagation frame.) Rebuilt per sample from the masses — this is how AD
    # flows through them.
    Ainv::SMatrix{NB,NR,T,L}
    # Per-band fluxes as body-vectors, for photocentre weights.
    fluxes::FL
    frame::F
end

function System(root::Union{Body,Binary}; kwargs...)
    frame = _make_frame(; kwargs...)
    leaves = _leaves(root)
    NB = length(leaves)
    rows_raw = _rows(root, 0)
    NR = length(rows_raw)
    names = _bodynames(leaves)
    _allunique_tuple(names) || error("body names must be unique, got $names")

    # Promote every scalar to a common float type.
    T = promote_type(map(l -> typeof(l.mass), leaves)...,
        map(r -> typeof(r.node.a), rows_raw)...,
        _frame_scalar_type(frame))
    masses = SVector{NB,T}(map(l -> l.mass, leaves))

    rows = map(rows_raw) do r
        bn = r.node
        M = _sum_masses(masses, r.int) + _sum_masses(masses, r.ext)
        Row(T(bn.a), T(bn.e), T(bn.i), T(bn.ω), T(bn.Ω), T(bn.tp), M)
    end

    # Ainv[j,k]: coefficient of row k's relative state in body j's absolute
    # barycentric state: -M_ext/M_row for interior members, +M_int/M_row for
    # exterior members, 0 otherwise (Hamers & Portegies Zwart 2016).
    Ainv = SMatrix{NB,NR,T}(_ainv_entries(masses, rows_raw, Val(NB))...)

    fluxes = _collect_fluxes(leaves, Val(NB), T)
    frameT = _convert_frame(T, frame)
    return System{NB,NR,T,typeof(frameT),names,typeof(fluxes),NB * NR}(
        masses, rows, Ainv, fluxes, frameT)
end

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

function _ainv_entries(masses::SVector{NB,T}, rows_raw, ::Val{NB}) where {NB,T}
    ntuple(Val(NB * length(rows_raw))) do lin
        k, j = divrem(lin - 1, NB) .+ (1, 1)
        r = rows_raw[k]
        Mint = _sum_masses(masses, r.int)
        Mext = _sum_masses(masses, r.ext)
        Mrow = Mint + Mext
        if j in r.int
            -Mext / Mrow
        elseif j in r.ext
            Mint / Mrow
        else
            zero(T)
        end
    end
end

# Default names: positional :body1, :body2, … for unnamed bodies.
@inline function _bodynames(leaves::NTuple{NB,Any}) where {NB}
    ntuple(Val(NB)) do j
        nm = _name(leaves[j])
        nm === :auto ? Symbol(:body, j) : nm
    end
end

# allocation-free allunique for small tuples (Base's builds a Set)
@inline _allunique_tuple(::Tuple{}) = true
@inline _allunique_tuple(t::Tuple) =
    !(first(t) in Base.tail(t)) && _allunique_tuple(Base.tail(t))

# Union of all bands across bodies; missing bands are zero flux.
function _collect_fluxes(leaves, ::Val{NB}, ::Type{T}) where {NB,T}
    protos = merge(map(l -> map(_ -> zero(T), l.flux), leaves)...)
    bands = keys(protos)
    vals = map(bands) do band
        SVector{NB,T}(ntuple(j -> T(get(leaves[j].flux, band, zero(T))), Val(NB)))
    end
    return NamedTuple{bands}(vals)
end

nbodies(::System{NB}) where {NB} = NB
nrows(::System{NB,NR}) where {NB,NR} = NR

"""
    bodies(sys)

NamedTuple of persistent `BodyRef`s for the bodies of `sys`, keyed by the
names given at construction (`Body(… , name=:b)`), in depth-first tree order.
These references are what observables accept:

    (; A, b) = bodies(sys)
    raoff(sol, b, A)
"""
bodies(::System{NB,NR,T,F,Names}) where {NB,NR,T,F,Names} =
    NamedTuple{Names}(ntuple(BodyRef, Val(NB)))

"""
    barycentre(sys)
    barycentre(sys, refs...)

The mass-weighted barycentre of the whole system, or of the subsystem
spanned by the given `BodyRef`s, as a `WeightedPoint`.
"""
function barycentre(sys::System{NB,NR,T}) where {NB,NR,T}
    return WeightedPoint{NB,T}(sys.masses ./ sum(sys.masses))
end
function barycentre(sys::System{NB,NR,T}, refs::BodyRef...) where {NB,NR,T}
    total = _sum_masses(sys.masses, map(r -> r.idx, refs))
    w = ntuple(Val(NB)) do j
        j in map(r -> r.idx, refs) ? sys.masses[j] / total : zero(T)
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

export System, Body, Binary, bodies, barycentre, photocentre, BodyRef
