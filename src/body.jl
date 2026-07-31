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

# ---------------------------------------------------
# Endpoint specs: what an orbit binds
#
# An orbit is (exterior spec, about = interior spec). Both take the same
# grammar: a single `Body`, or a tuple of `Body`s meaning *the barycentre of
# that set*. This is what makes convention explicit —
#
#   Orbit(b, about=A)       astrocentric: b about the star alone
#   Orbit(c, about=(A, b))  Jacobi: c about the A+b barycentre
#   Orbit((Ba,Bb), about=(Aa,Ab))   2+2 quadruple: set exteriors
#
# — and there is deliberately no default and no inference.
#
# The whole-system barycentre and photocentres are *observable* reference
# points, not orbit endpoints: an orbit about the system barycentre includes
# the orbiting body in its own reference, and a photocentre is not a
# dynamical quantity at all. Both remain available to observables via
# `barycentre(sys)` / `photocentre(sys)`.
# ---------------------------------------------------

const BodySpec = Union{Body,Tuple{Vararg{Body}}}

@inline _specmass(b::Body) = b.mass
@inline _specmass(t::Tuple{Vararg{Body}}) = _summasses(t)
@inline _summasses(::Tuple{}) = 0.0
@inline _summasses(t::Tuple) = first(t).mass + _summasses(Base.tail(t))

@inline _specnames(b::Body) = (_name(b),)
@inline _specnames(t::Tuple{Vararg{Body}}) = map(_name, t)

@inline _specbodies(b::Body) = (b,)
@inline _specbodies(t::Tuple{Vararg{Body}}) = t

@inline _check_spec(::BodySpec) = nothing
@noinline _check_spec(::AbstractRef) = error(
    "got a reference into an existing System (as returned by `bodies`/`barycentre`) " *
    "where a `Body` value is required. `Orbit` builds new systems from construction-time " *
    "values, e.g. `Orbit(Body(mass=1mjup, name=:b), about=Body(mass=1.0, name=:A); a=…)`.")
@noinline _check_spec(@nospecialize(x)) = error(
    "an orbit endpoint must be a `Body`, or a tuple of `Body`s meaning their " *
    "barycentre (e.g. `about=(A, b)`); got a value of type $(typeof(x))")

# ---------------------------------------------------
# Orbit
# ---------------------------------------------------

"""
    Orbit(exterior; about, <size>, e=0, i=0, ω=0, Ω=0, tp=0, M=nothing)

One Keplerian relationship in a system: the orbit of `exterior` **about**
`about`. Both endpoints take the same grammar — a `Body`, or a tuple of
`Body`s denoting their barycentre:

    Orbit(b, about=A;      a=5.2, e=0.05)   # astrocentric
    Orbit(c, about=(A, b); a=9.6)           # Jacobi
    Orbit((Ba, Bb), about=(Aa, Ab); a=120)  # 2+2 quadruple

The convention is always explicit: there is no default `about=` and no
inference from semi-major axis. See `Jacobi` and `Astrocentric` for the two
standard chains.

# Size (supply exactly one)
  - `a` — semi-major axis [AU]
  - `P` — period [**days**], converted via the row's gravitating mass

# Other elements
Eccentricity `e`, inclination `i` [rad], argument of periapsis `ω` [rad],
longitude of ascending node `Ω` [rad], epoch of periastron passage `tp` [MJD].

# `M` (compatibility escape hatch)
The row's gravitating mass is normally the total mass of every body the row
binds, taken from the `Body` values themselves. Passing `M` [M⊙] overrides
it. This is **physically inconsistent** — it decouples Kepler's third law
from the masses that drive the reflex amplitudes — and exists only to
reproduce published fits bit-for-bit and to match orbitize!/RadVel
conventions. It is a compatibility switch, not a modelling choice.
"""
struct Orbit{E,I,T<:Number}
    exterior::E
    interior::I
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    tp::T
    M::T          # gravitating mass of this row [M⊙]
    Moverride::Bool
end

function Orbit(exterior; about, a=nothing, P=nothing,
               e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0, M=nothing)
    _check_spec(exterior)
    _check_spec(about)
    # The row's gravitating mass is known here: `about` carries Body values,
    # so P→a can be done at construction (§4.3) rather than deferred.
    Mrow = M === nothing ? _specmass(exterior) + _specmass(about) : float(M)
    a_ = _size_from(a, P, Mrow)
    a_, e, i, ω, Ω, tp, Mrow = promote(float(a_), e, i, ω, Ω, tp, Mrow)
    return Orbit{typeof(exterior),typeof(about),typeof(a_)}(
        exterior, about, a_, e, i, ω, Ω, tp, Mrow, M !== nothing)
end

# --- size group: exactly one of `a` | `P` -----------------------------------
@inline _size_from(a, ::Nothing, M) = a
@inline _size_from(::Nothing, P, M) = _a_from_P(P, M)
@noinline _size_from(::Nothing, ::Nothing, M) = error(
    "supply exactly one of `a` [AU] or `P` [days] to `Orbit`; got neither")
@noinline _size_from(a, P, M) = error(
    "supply exactly one of `a` [AU] or `P` [days] to `Orbit`; got both " *
    "(a=$a, P=$P)")

# Kepler's third law, inverting Row's period_days = √(a³/M)·kepler_year_factor.
# NB: P is in DAYS, matching `period(sys)` so the two round-trip. Imaging
# users who think in years will get a plausible-looking 365× error, so `show`
# prints the period in both units.
@inline function _a_from_P(P, M)
    M > 0 || error("cannot convert P→a for a row of zero gravitating mass; " *
                   "give the bodies masses or pass `a` directly")
    Pyr = float(P) / kepler_year_to_julian_day_conversion_factor
    return cbrt(Pyr^2 * M)
end

# ---------------------------------------------------
# Convention constructors
#
# Both expand to explicit `about=` specs — greppable, and `show(sys)` prints
# the convention it resolved for each row.
# ---------------------------------------------------

"""
    Jacobi(inner, b => (; a=…, e=…), c => (; a=…), …)

Build a Jacobi chain: each body orbits the barycentre of everything interior
to it. Equivalent to spelling out

    Orbit(b, about=inner;       …)
    Orbit(c, about=(inner…, b); …)

Returns a tuple of `Orbit`s for `System`. `inner` is a `Body` or a tuple of
`Body`s (e.g. a tight binary at the centre of a circumbinary chain).
"""
function Jacobi(inner::BodySpec, pairs::Pair...)
    interiors = _accumulate_interiors(_specbodies(inner), pairs)
    return ntuple(length(pairs)) do k
        Orbit(first(pairs[k]); about=interiors[k], last(pairs[k])...)
    end
end

"""
    Astrocentric(centre, b => (; a=…), c => (; a=…), …)

Build an astrocentric set: every body orbits `centre` directly. Equivalent to
`Orbit(b, about=centre; …)`, `Orbit(c, about=centre; …)`, ….

Note this is a materially different model from `Jacobi` under
`KeplerianApprox` (not a relabelling): the rows *are* the approximation.
Under `AHL21` the two agree, since rows only set initial conditions.
"""
function Astrocentric(centre::BodySpec, pairs::Pair...)
    return ntuple(length(pairs)) do k
        Orbit(first(pairs[k]); about=centre, last(pairs[k])...)
    end
end

# interiors[k] = every body interior to the k-th exterior, as a flat tuple
@inline function _accumulate_interiors(inner::Tuple, pairs::Tuple)
    ntuple(length(pairs)) do k
        _flatten_upto(inner, pairs, k)
    end
end
@inline _flatten_upto(inner::Tuple, pairs::Tuple, k::Int) =
    k == 1 ? inner : (_flatten_upto(inner, pairs, k - 1)...,
                      _specbodies(first(pairs[k-1]))...)

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
# Row membership masks
#
# A row's interior/exterior member sets are stored as fixed-width indicator
# masks rather than variable-length index tuples. This is load-bearing for
# performance, not a style choice: with index tuples the rows tuple is
# *heterogeneous* ((1,), (1,2), (1,2,3), …), which defeats constant folding
# in the A⁻¹ build and sends it to the heap — measured at 10.2 µs and 26 kB
# for a 5-body system, versus 0.24 µs and 0 bytes with masks.
# ---------------------------------------------------

struct RowSpec{NB}
    int::SVector{NB,Bool}
    ext::SVector{NB,Bool}
end

@inline function RowSpec(names::NTuple{NB,Symbol}, o::Orbit) where {NB}
    return RowSpec{NB}(_mask(names, o.interior, Val(NB)),
                       _mask(names, o.exterior, Val(NB)))
end

@inline function _mask(names::NTuple{NB,Symbol}, spec, ::Val{NB}) where {NB}
    ns = _specnames(spec)
    return SVector{NB,Bool}(ntuple(j -> _innames(ns, names[j]), Val(NB)))
end

@inline _innames(::Tuple{}, ::Symbol) = false
@inline _innames(ns::Tuple, nm::Symbol) =
    first(ns) === nm || _innames(Base.tail(ns), nm)
