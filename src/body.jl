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

# Elements
Supply exactly one alternative from each group. Orientation is always `i`
[rad] and `Ω` [rad]; every angle is in radians and every epoch is MJD.

| group | alternatives | notes |
|---|---|---|
| size | `a` [AU] \\| `P` [**days**] | `P` uses the row's gravitating mass |
| shape | (`e`, `ω`) \\| (`secosω`, `sesinω`) \\| (`ecosω`, `esinω`) | default `e=0`, `ω=0` |
| phase | `tp` \\| `M0` + `epoch` \\| `θ` + `epoch` | default `tp=0` |

`secosω` = √e·cosω and `ecosω` = e·cosω sample the eccentricity disc rather
than the half-plane, removing the ω degeneracy at e → 0. `M0` is the mean
anomaly at `epoch` [rad]; `θ` is the sky-plane position angle at `epoch`
[rad]. `τ` is deliberately not accepted — it needs hidden period and
reference-epoch state and has no clean meaning under N-body integration.

**`P` is in days**, matching `period(sys)` so the two round-trip. Users who
think in years get a plausible-looking 365× error rather than a crash, so
`show` prints the period in both units.

# Cartesian initial conditions
Alternatively give the full relative state, which determines every element
and so replaces *all* of the groups above:

    Orbit(b, about=A; x=1.2, y=0.3, z=-0.1,       # [AU]
                      vx=-0.9, vy=2.1, vz=0.2,    # [AU / julian yr]
                      epoch=59000.0)              # [MJD]

Positions and velocities are of `exterior` relative to `about`, in the same
frame and units as `posx`/`velx`, so a state read back out of a solution
reconstructs the same orbit. Unbound states are fine — `a` comes out negative
and `e > 1` — which makes this the natural way to specify hyperbolic orbits.

# Eccentricity
`e > 1` is supported (unbound); `a` is negative by convention there and a
positive value is taken as |a|. `e == 1` exactly is rejected: the elements are
degenerate for parabolae. Use Cartesian initial conditions instead.

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

function Orbit(exterior; about,
               # size — exactly one
               a=nothing, P=nothing,
               # shape — at most one group
               e=nothing, ω=nothing,
               secosω=nothing, sesinω=nothing,
               ecosω=nothing, esinω=nothing,
               # orientation
               i=nothing, Ω=nothing,
               # phase — at most one, `M0`/`θ` need `epoch`
               tp=nothing, M0=nothing, θ=nothing, epoch=nothing,
               # joint: a Cartesian relative state replaces every group above
               x=nothing, y=nothing, z=nothing,
               vx=nothing, vy=nothing, vz=nothing,
               # compatibility
               M=nothing)
    _check_spec(exterior)
    _check_spec(about)
    # The row's gravitating mass is known here: `about` carries Body values,
    # so P→a — and every other reparametrization — can be done at
    # construction (§4.3) rather than deferred.
    Mrow = M === nothing ? _specmass(exterior) + _specmass(about) : float(M)

    el = if any(!isnothing, (x, y, z, vx, vy, vz))
        all(!isnothing, (x, y, z, vx, vy, vz)) || error(
            "Cartesian initial conditions need all six of `x`, `y`, `z` [AU] and " *
            "`vx`, `vy`, `vz` [AU/julian yr]")
        any(!isnothing, (a, P, e, ω, secosω, sesinω, ecosω, esinω, i, Ω, tp, M0, θ)) && error(
            "Cartesian initial conditions determine every orbital element, so they " *
            "cannot be combined with `a`/`P`, `e`/`ω`, `i`/`Ω`, or a phase keyword")
        epoch === nothing && error(
            "Cartesian initial conditions need `epoch=` [MJD]: a state is a state *at a time*")
        _elements_from_state(promote(float(x), y, z, vx, vy, vz)..., float(epoch), Mrow)
    else
        e_, ω_ = _shape_from(e, ω, secosω, sesinω, ecosω, esinω)
        i_ = i === nothing ? zero(e_) : float(i)
        Ω_ = Ω === nothing ? zero(e_) : float(Ω)
        a_ = _size_from(a, P, Mrow)
        tp_ = _phase_from(tp, M0, θ, epoch, a_, e_, i_, ω_, Ω_, Mrow)
        (a_, e_, i_, ω_, Ω_, tp_)
    end

    a_, e_, i_, ω_, Ω_, tp_, Mrow = promote(map(float, el)..., Mrow)
    return Orbit{typeof(exterior),typeof(about),typeof(a_)}(
        exterior, about, a_, e_, i_, ω_, Ω_, tp_, Mrow, M !== nothing)
end

# --- size group: exactly one of `a` | `P` -----------------------------------
@inline _size_from(a, ::Nothing, M) = float(a)
@inline _size_from(::Nothing, P, M) = _a_from_P(P, M)
@noinline _size_from(::Nothing, ::Nothing, M) = error(
    "supply exactly one of `a` [AU] or `P` [days] to `Orbit`; got neither")
@noinline _size_from(a, P, M) = error(
    "supply exactly one of `a` [AU] or `P` [days] to `Orbit`; got both " *
    "(a=$a, P=$P)")

# --- shape group: (e, ω) | (√e·cosω, √e·sinω) | (e·cosω, e·sinω) ------------
# The latter two sample the eccentricity disc rather than the half-plane, so
# they have no ω degeneracy at e → 0 — the usual choice for MCMC.
function _shape_from(e, ω, secosω, sesinω, ecosω, esinω)
    g1 = !isnothing(e) || !isnothing(ω)
    g2 = !isnothing(secosω) || !isnothing(sesinω)
    g3 = !isnothing(ecosω) || !isnothing(esinω)
    g1 + g2 + g3 > 1 && _err_shape()
    if g2
        (!isnothing(secosω) && !isnothing(sesinω)) ||
            error("`secosω` and `sesinω` must be given together")
        sc, ss = promote(float(secosω), float(sesinω))
        return sc^2 + ss^2, atan(ss, sc)
    elseif g3
        (!isnothing(ecosω) && !isnothing(esinω)) ||
            error("`ecosω` and `esinω` must be given together")
        ec, es = promote(float(ecosω), float(esinω))
        return hypot(ec, es), atan(es, ec)
    end
    ee = e === nothing ? 0.0 : float(e)
    return promote(ee, ω === nothing ? zero(ee) : float(ω))
end
@noinline _err_shape() = error(
    "supply at most one eccentricity parametrization: (`e`, `ω`), " *
    "(`secosω`, `sesinω`), or (`ecosω`, `esinω`)")

# --- phase group: tp | M0(+epoch) | θ(+epoch) -------------------------------
# `τ` is deliberately absent: it needs hidden period/reference-epoch state and
# has no clean meaning under N-body integration (§5).
function _phase_from(tp, M0, θ, epoch, a, e, i, ω, Ω, M)
    n = !isnothing(tp) + !isnothing(M0) + !isnothing(θ)
    n > 1 && _err_phase()
    tp !== nothing && return float(tp)
    n == 0 && return zero(float(a))
    epoch === nothing && error(
        "the `$(M0 !== nothing ? "M0" : "θ")` phase parametrization is measured at " *
        "an epoch: pass `epoch=` [MJD]")
    MA = M0 !== nothing ? float(M0) : _MA_from_θ(float(θ), e, i, ω, Ω)
    return float(epoch) - MA / _meanmotion(a, e, M) * year2day_julian
end
@noinline _err_phase() = error(
    "supply at most one phase parametrization: `tp`, `M0` (+`epoch`), or `θ` (+`epoch`)")

@inline function _meanmotion(a, e, M)
    μ = GM_sun_au3_julianyr2 * M
    e < 1 && return √(μ / a^3)
    aa = a > 0 ? -a : a
    return √(μ / (-aa)^3)
end

# Sky-plane position angle θ at an epoch → mean anomaly.
#
# Needs neither `a` nor the masses: the Thiele-Innes constants are built from
# the orientation angles alone, and deprojecting the *direction* (cos θ, sin θ)
# cancels the radius factor. (Octofitter's θ_at_epoch_to_tperi takes `a`/`M`
# only to reach the mean motion, and its `P` branch is broken — it reads
# `system.M` and `b.P`, neither in scope, and scales by the kepler-year factor
# where it should divide by its 2/3 power. That path is unreachable here: the
# constructor always has both a and the row mass.)
function _MA_from_θ(θ, e, i, ω, Ω)
    e < 1 || error("the `θ` phase parametrization is only defined for e < 1")
    sinΩ, cosΩ = sincos(Ω)
    sinω, cosω = sincos(ω)
    cosi = cos(i)
    A = cosΩ * cosω - sinΩ * sinω * cosi
    B = sinΩ * cosω + cosΩ * sinω * cosi
    F = -cosΩ * sinω - sinΩ * cosω * cosi
    G = -sinΩ * sinω + cosΩ * cosω * cosi
    det = A * G - B * F
    sθ, cθ = sincos(θ)
    xr = (G * cθ - F * sθ) / det          # T \ [cos θ; sin θ], T = [A F; B G]
    yr = (-B * cθ + A * sθ) / det
    ν = atan(yr, xr)
    sν, cν = sincos(ν)
    E = atan(√(1 - e^2) * sν, e + cν)
    return E - e * sin(E)
end

# --- joint group: Cartesian relative state → osculating elements ------------
#
# Conic-agnostic: `a` falls out of vis-viva (negative for unbound states) and
# `e` from the eccentricity vector, so hyperbolic initial conditions need no
# special case. The perifocal basis is read off the forward kernel in
# `_states_from_E`, which places
#     P̂ = (sinΩ, cosΩ, 0),  Q̂ = (cosi·cosΩ, −cosi·sinΩ, sini)
# so the orbit normal is Ŵ = P̂ × Q̂ = (sini·cosΩ, −sini·sinΩ, −cosi) — note
# the sign on the z component, which is where a textbook formula would differ.
function _elements_from_state(x, y, z, vx, vy, vz, epoch, M)
    T = typeof(x)
    μ = T(GM_sun_au3_julianyr2) * M        # AU³ / (julian yr)², matching Row
    r = SVector(x, y, z)
    v = SVector(vx, vy, vz)
    rn = √(r ⋅ r)
    rn > 0 || error("Cartesian initial conditions have zero separation")
    h = r × v
    hn = √(h ⋅ h)
    hn > 0 || error("Cartesian initial conditions are radial (zero angular " *
                    "momentum); the orbital plane is undefined")
    ĥ = h ./ hn
    inc = acos(clamp(-ĥ[3], -one(T), one(T)))
    sini = hypot(ĥ[1], ĥ[2])
    Ω = sini < eps(T)^(one(T) / 2) ? zero(T) : atan(-ĥ[2], ĥ[1])
    sΩ, cΩ = sincos(Ω)
    P̂ = SVector(sΩ, cΩ, zero(T))
    Q̂ = ĥ × P̂
    evec = (v × h) ./ μ .- r ./ rn
    e = √(evec ⋅ evec)
    a = inv(2 / rn - (v ⋅ v) / μ)
    u = atan(r ⋅ Q̂, r ⋅ P̂)                 # argument of latitude, ν + ω
    ω = e < eps(T)^(one(T) / 2) ? zero(T) : atan(evec ⋅ Q̂, evec ⋅ P̂)
    ν = u - ω
    sν, cν = sincos(ν)
    MA = if e < 1
        E = atan(√(1 - e^2) * sν, e + cν)
        E - e * sin(E)
    else
        H = asinh(√(e^2 - 1) * sν / (1 + e * cν))
        e * sinh(H) - H
    end
    tp = epoch - MA / _meanmotion(a, e, M) * year2day_julian
    return (a, e, inc, ω, Ω, tp)
end

# Kepler's third law, inverting Row's period_days = √(a³/M)·kepler_year_factor.
# NB: P is in DAYS, matching `period(sys)` so the two round-trip. Imaging
# users who think in years will get a plausible-looking 365× error, so `show`
# prints the period in both units.
#
# Unbound orbits have no period, so `P=` is elliptical-only — give `a` (negative,
# by convention) instead.
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
