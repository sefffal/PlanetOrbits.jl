# ---------------------------------------------------
# Bodies, Binaries, and references
# ---------------------------------------------------

"""
    Body(; mass, flux=(;), name=:auto)

A point mass participating in a `System`. `mass` is in solar masses (see the
`msun`, `mjup`, `mearth` constants). `flux` is an optional NamedTuple of
per-band fluxes (arbitrary but consistent units within a band), used to
compute photocentres.

`Body` values are construction-time inputs, typically rebuilt every MCMC
sample. They are *not* handles into a system: after building a `System`, use
`bodies(sys)` to obtain the persistent references that observables accept.
"""
struct Body{T<:Number,F<:NamedTuple,Name}
    mass::T
    flux::F
end
function Body(; mass, flux::NamedTuple=(;), name::Symbol=:auto)
    mass = float(mass)
    # The name is a *type* parameter: body identity is static per model (like
    # the topology itself), which keeps `System`'s type — and therefore the
    # whole construct→solve→query path — inferable. Pass `name` as a literal.
    return Body{typeof(mass),typeof(flux),name}(mass, flux)
end

_name(::Body{T,F,Name}) where {T,F,Name} = Name::Symbol

"""
    Binary(interior, exterior; a, e=0, i=0, ω=0, Ω=0, tp=0)

A Keplerian orbit binding `exterior` (a `Body` or another `Binary`) to
`interior` (likewise). The elements describe the orbit of the barycentre of
`exterior` **about the barycentre of `interior`** (Jacobi convention), so
argument order is meaningful.

Elements: semi-major axis `a` [AU], eccentricity `e`, inclination `i` [rad],
argument of periapsis `ω` [rad], longitude of ascending node `Ω` [rad], and
epoch of periastron passage `tp` [MJD].
"""
struct Binary{I,O,T<:Number}
    interior::I
    exterior::O
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    tp::T
end
function Binary(interior::Union{Body,Binary}, exterior::Union{Body,Binary};
                a, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0)
    a, e, i, ω, Ω, tp = promote(float(a), e, i, ω, Ω, tp)
    return Binary{typeof(interior),typeof(exterior),typeof(a)}(interior, exterior, a, e, i, ω, Ω, tp)
end

"""
    BodyRef

Persistent, `isbits` reference to a body of a `System` (an index into its
body list). Obtained from `bodies(sys)`; valid across samples for any system
sharing the same topology. This is what observables accept as a target or
reference point.
"""
struct BodyRef
    idx::Int
end

"""
    WeightedPoint

A normalized weighted combination of body states — the generalization of a
single body used to represent barycentres (mass weights) and photocentres
(flux weights). Construct via `barycentre` or `photocentre`. Valid only for
the sample whose masses/fluxes produced it.
"""
struct WeightedPoint{NB,T<:Number}
    w::SVector{NB,T}
end

const AbstractRef = Union{BodyRef,WeightedPoint}

# ---------------------------------------------------
# Tree flattening (type-stable tuple recursion; runs per sample)
# ---------------------------------------------------

# Leaves in depth-first order: interior subtree first, then exterior.
# This ordering defines the body indices used by BodyRef.
@inline _leaves(b::Body) = (b,)
@inline _leaves(bn::Binary) = (_leaves(bn.interior)..., _leaves(bn.exterior)...)

# Statically-known leaf counts (pure functions of the tree *type*, so they
# constant-fold and the Val-based ntuples below stay allocation-free).
@inline _nleaves(::Body) = 1
@inline _nleaves(bn::Binary) = _nleaves(bn.interior) + _nleaves(bn.exterior)

# Rows in post-order (children before parents). Each row records the Binary
# node and the index ranges of its interior/exterior member bodies. For a
# star + N-planet Jacobi chain this yields rows innermost-first.
@inline _rows(::Body, offset::Int) = ()
@inline function _rows(bn::Binary, offset::Int)
    nint = _nleaves(bn.interior)
    rows_int = _rows(bn.interior, offset)
    rows_ext = _rows(bn.exterior, offset + nint)
    row = (
        node = bn,
        int = ntuple(k -> offset + k, Val(_nleaves(bn.interior))),
        ext = ntuple(k -> offset + nint + k, Val(_nleaves(bn.exterior))),
    )
    return (rows_int..., rows_ext..., row)
end
