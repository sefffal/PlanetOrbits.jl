# ---------------------------------------------------
# Bodies, Orbits, and references
# ---------------------------------------------------

"""
    Body(; mass, flux=(;), name=:auto)

A point mass participating in a `System`. `mass` is in solar masses (see the
`msun`, `mjup`, `mearth` constants). `flux` is an optional NamedTuple of
per-band fluxes (arbitrary but consistent units within a band), used to
compute photocentres.

`Body` values are construction-time inputs, typically rebuilt every MCMC
sample. A *named* `Body` doubles as a handle: observables resolve it by name
against the solved system (`raoff(sol, b, A)` with the same variables used
to build it), reading only its name — a value left over from a previous
sample resolves identically. Unnamed bodies (and hot loops) use the
persistent references from `bodies(sys)` instead.
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
    Orbit(body, about=other; a, e=0, i=0, ω=0, Ω=0, tp=0)

A Keplerian orbit binding `body` (a `Body` or another `Orbit`) to `about`
(likewise). The elements describe the orbit of the barycentre of `body`
**about the barycentre of `about`** (Jacobi convention):

    Orbit(jup, about=sun; a=5.2, e=0.0489, …)          # Jupiter orbits the Sun
    Orbit(Orbit(gan, about=jup; …), about=sun; a=5.2)  # …with a moon on Jupiter

Elements: semi-major axis `a` [AU], eccentricity `e`, inclination `i` [rad],
argument of periapsis `ω` [rad], longitude of ascending node `Ω` [rad], and
epoch of periastron passage `tp` [MJD].
"""
struct Orbit{I,O,T<:Number}
    interior::I
    exterior::O
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    tp::T
end
function Orbit(body; about, a, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0)
    _check_member(body)
    _check_member(about)
    a, e, i, ω, Ω, tp = promote(float(a), e, i, ω, Ω, tp)
    return Orbit{typeof(about),typeof(body),typeof(a)}(about, body, a, e, i, ω, Ω, tp)
end

"""
    BodyRef

Persistent, `isbits` reference to a body of a `System` (an index into its
body list). Obtained from `bodies(sys)`; valid across samples for any system
sharing the same topology. This is the resolved form observables work with —
interactively, named `Body` values and `Symbol`s are accepted anywhere a
`BodyRef` is and resolve to one by name.
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

@inline _check_member(::Union{Body,Orbit}) = nothing
@noinline _check_member(::AbstractRef) = error(
    "got a reference into an existing System (as returned by `bodies`/`barycentre`) " *
    "where a `Body` value is required. `Orbit` builds new systems from construction-time " *
    "values, e.g. `Orbit(Body(mass=1mjup, name=:b), about=Body(mass=1.0, name=:A); a=…)`.")
@noinline _check_member(@nospecialize(x)) = error(
    "`Orbit` members must be `Body` or `Orbit` values, got a value of type $(typeof(x))")

# ---------------------------------------------------
# Name-based reference resolution
#
# Observables (and `barycentre`) accept, anywhere a ref is expected:
#   - a `BodyRef`/`WeightedPoint` (passed through untouched),
#   - a named `Body` value — resolved by its name type-parameter, so only
#     identity is read and "stale" values from a previous sample are fine,
#   - a `Symbol`.
# `names` is the system's name table (a type-level NTuple of Symbols), so
# resolution against literal names constant-folds to a plain BodyRef.
# ---------------------------------------------------

@inline _resolve(::Tuple, r::AbstractRef) = r
@inline function _resolve(names::Tuple, b::Body)
    nm = _name(b)
    nm === :auto && _resolve_err_auto()
    return _resolve(names, nm)
end
@inline function _resolve(names::Tuple, nm::Symbol)
    k = _findname(names, nm, 1)
    k === 0 && _resolve_err_missing(nm, names)
    return BodyRef(k)
end

@inline _findname(::Tuple{}, ::Symbol, ::Int) = 0
@inline _findname(names::Tuple, nm::Symbol, k::Int) =
    first(names) === nm ? k : _findname(Base.tail(names), nm, k + 1)

@noinline _resolve_err_auto() = error(
    "this Body was constructed without a name, so it cannot be resolved against a " *
    "system. Pass `name=…` at construction (`Body(mass=…, name=:b)`), or use the " *
    "references from `bodies(sys)`.")
@noinline _resolve_err_missing(nm::Symbol, names) = error(
    "no body named :$nm in this system (its bodies are named $names)")

# ---------------------------------------------------
# Tree flattening (type-stable tuple recursion; runs per sample)
# ---------------------------------------------------

# Leaves in depth-first order: interior subtree first, then exterior.
# This ordering defines the body indices used by BodyRef.
@inline _leaves(b::Body) = (b,)
@inline _leaves(bn::Orbit) = (_leaves(bn.interior)..., _leaves(bn.exterior)...)

# Statically-known leaf counts (pure functions of the tree *type*, so they
# constant-fold and the Val-based ntuples below stay allocation-free).
@inline _nleaves(::Body) = 1
@inline _nleaves(bn::Orbit) = _nleaves(bn.interior) + _nleaves(bn.exterior)

# Rows in post-order (children before parents). Each row records the Orbit
# node and the index ranges of its interior/exterior member bodies. For a
# star + N-planet Jacobi chain this yields rows innermost-first.
@inline _rows(::Body, offset::Int) = ()
@inline function _rows(bn::Orbit, offset::Int)
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
