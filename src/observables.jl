# ---------------------------------------------------
# Observables
#
# Every observable is a function of a *pair of references* — a target and a
# reference point, each a `BodyRef` or `WeightedPoint` (barycentre,
# photocentre) — evaluated on a per-epoch solution `sol = traj[k]`. Named
# `Body` values and `Symbol`s are accepted in either position and resolve
# by name (see `_resolve` in body.jl).
#
# Angular quantities [mas] require the system to have (at least) a parallax;
# this is enforced by dispatch on the trajectory's frame mode.
# ---------------------------------------------------

# State component reads: direct column load for a BodyRef, small dot product
# for a WeightedPoint.
for (fn, col) in ((:_sx, :x), (:_sy, :y), (:_sz, :z),
                  (:_svx, :vx), (:_svy, :vy), (:_svz, :vz))
    @eval @inline function $fn(sol::TrajectorySolution, ref::BodyRef)
        @inbounds sol.traj.$col[sol.k, ref.idx]
    end
    @eval @inline function $fn(sol::TrajectorySolution, p::WeightedPoint{NB}) where {NB}
        m = sol.traj.$col
        k = sol.k
        acc = zero(eltype(m))
        @inbounds for j in 1:NB
            acc = muladd(p.w[j], m[k, j], acc)
        end
        return acc
    end
end

# Solution aliases by frame mode. NB: written with the frame-mode parameter
# in an invariant position — a `<:Trajectory{<:Any,FM}` pattern would put FM
# only in another typevar's bound, where Julia's dispatcher ignores it.
const _SolWithMode{FM} = TrajectorySolution{Trajectory{T,FM,Names,E,V,M}} where
    {T<:Number,Names,E<:AbstractVector{<:Real},V<:AbstractVector{T},M<:AbstractMatrix{T}}
const _AngularSol = Union{_SolWithMode{ModeParallax},_SolWithMode{ModeAbsolute}}
const _AbsSol = _SolWithMode{ModeAbsolute}

# Frame guards: mas-valued observables need a distance. The AU -> mas factor
# is per *body*, not per epoch — a body at line-of-sight depth z subtends its
# offset over d+z (see observe.jl).
@inline _cart2angle(sol::_AngularSol, j::Int) =
    @inbounds sol.traj.cart2angle[sol.k, j]
@noinline _cart2angle(::TrajectorySolution, ::Int) =
    error("this system has no parallax: angular observables (mas) are unavailable. " *
          "Construct the System with plx=… (or a full absolute frame), or use the " *
          "physical-unit observables posx/posy/posz/velx/vely/velz.")

# Gnomonic (tangent-plane) coordinates ξ, η [mas] of a reference about the
# barycentre's apparent direction, and their exact time derivatives [mas/yr].
# `raoff`/`pmra` are differences of these, so the pair is consistent by
# construction. For a weighted point the weights apply to the *angular*
# coordinates — a photocentre is the flux-weighted mean of apparent
# positions, not of linear offsets.
for (fn, dfn, cp, cv) in ((:_ax, :_apmx, :x, :vx), (:_ay, :_apmy, :y, :vy))
    @eval @inline function $fn(sol::TrajectorySolution, ref::BodyRef)
        @inbounds sol.traj.$cp[sol.k, ref.idx] * _cart2angle(sol, ref.idx)
    end
    @eval @inline function $fn(sol::TrajectorySolution, p::WeightedPoint{NB}) where {NB}
        k = sol.k
        acc = zero(eltype(sol.traj.$cp))
        @inbounds for j in 1:NB
            acc = muladd(p.w[j] * sol.traj.$cp[k, j], _cart2angle(sol, j), acc)
        end
        return acc
    end
    # d/dt [ q / (d + z) ] = q̇/(d+z) − q(ḋ + ż)/(d+z)²
    @eval @inline function $dfn(sol::TrajectorySolution, ref::BodyRef)
        k = sol.k
        j = ref.idx
        c = _cart2angle(sol, j)
        @inbounds c * (sol.traj.$cv[k, j] -
                       sol.traj.$cp[k, j] * c * (sol.traj.bvz[k] + sol.traj.vz[k, j]) / rad2mas)
    end
    @eval @inline function $dfn(sol::TrajectorySolution, p::WeightedPoint{NB}) where {NB}
        k = sol.k
        acc = zero(eltype(sol.traj.$cp))
        @inbounds for j in 1:NB
            c = _cart2angle(sol, j)
            acc = muladd(p.w[j], c * (sol.traj.$cv[k, j] -
                    sol.traj.$cp[k, j] * c * (sol.traj.bvz[k] + sol.traj.vz[k, j]) / rad2mas),
                acc)
        end
        return acc
    end
end

# ---------------------------------------------------
# Physical-unit pairwise observables [AU, AU/julian year, m/s]
# ---------------------------------------------------

"""
    posx(sol, target, reference)

Position offset [AU] of `target` relative to `reference` along the
right-ascension (east) direction. `target`/`reference` are `BodyRef`s or
`WeightedPoint`s — or named `Body` values / `Symbol`s, resolved by name.

Resolved on the local East/North/line-of-sight triad of the system
barycentre's **apparent** direction at the observation epoch (not at
`ref_epoch`), with each body taken at its own light-travel-retarded time.
See `observe.jl`.
"""
@inline posx(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _sx(sol, t) - _sx(sol, r)
@inline posy(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _sy(sol, t) - _sy(sol, r)
@inline posz(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _sz(sol, t) - _sz(sol, r)
@inline velx(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svx(sol, t) - _svx(sol, r)
@inline vely(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svy(sol, t) - _svy(sol, r)
@inline velz(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svz(sol, t) - _svz(sol, r)

"""
    radvel(sol, target, reference)

Radial velocity [m/s] of `target` relative to `reference` along the line of
sight (positive receding). E.g. `radvel(sol, b, A)` for a relative RV, or
`radvel(sol, A, barycentre(sys))` for the stellar reflex.
"""
@inline radvel(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    velz(sol, t, r) * au2m / year2sec_julian

# ---------------------------------------------------
# Angular pairwise observables [mas, mas/julian year]
# ---------------------------------------------------

"""
    raoff(sol, target, reference)

Right-ascension offset [mas] of `target` relative to `reference`: the
difference of their gnomonic (tangent-plane) coordinates about the system
barycentre's *apparent* direction at the observation epoch, with each body
taken at its own light-travel-retarded time.

Note this is **not** `posx * cart2angle` for a single shared scale factor:
the two references are divided by their own `d + z`, which differs from the
shared-scale answer by ρ² in radians (≈ 4.85·ρ[″]² µas).
"""
@inline raoff(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    _ax(sol, t) - _ax(sol, r)

"""
    decoff(sol, target, reference)

Declination offset [mas] of `target` relative to `reference`. See `raoff`.
"""
@inline decoff(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    _ay(sol, t) - _ay(sol, r)

"""
    pmra(sol, target, reference)

Instantaneous relative proper motion [mas/julian year] of `target` with
respect to `reference` in right ascension — the exact time derivative of
`raoff`, so it carries the same per-body depth scaling.
"""
@inline pmra(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    _apmx(sol, t) - _apmx(sol, r)

@inline pmdec(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    _apmy(sol, t) - _apmy(sol, r)

"""
    projectedseparation(sol, target, reference)

Projected separation [mas] between `target` and `reference`.
"""
@inline function projectedseparation(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef)
    x = raoff(sol, t, r)
    y = decoff(sol, t, r)
    return sqrt(x^2 + y^2)
end

"""
    posangle(sol, target, reference)

Position angle [rad] of `target` about `reference`, measured from north
through east.

Nearly — but no longer exactly — parallax-free. For a body-vs-body pair the
per-body depth factor `1/(d + z)` is common to both components of the
separation and cancels, so the distance drops out as it always did. When the
two references sit at *different* line-of-sight depths (a body versus a
barycentre or photocentre) it does not cancel exactly, and the position angle
acquires a dependence on distance at the ~1e-6 rad level. On systems built
without `plx` the physical-unit fallback is used and no parallax is required.
"""
@inline function posangle(sol::_AngularSol, t::AbstractRef, r::AbstractRef)
    return atan(raoff(sol, t, r), decoff(sol, t, r))
end

@inline function posangle(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef)
    x = posx(sol, t, r)
    y = posy(sol, t, r)
    return atan(x, y)
end

# ---------------------------------------------------
# Propagated-frame accessors (AbsoluteFrame systems only)
#
# These return the *frame's* (i.e. the system barycentre's) apparent
# quantities at the epoch, after rigorous 3D space-motion propagation.
# Absolute quantities of a body compose as frame + pairwise-vs-barycentre:
#   pmra_abs(star) = frame_pmra(sol) + pmra(sol, A, barycentre(sys))
# ---------------------------------------------------

"""
    frame_ra(sol), frame_dec(sol)         [deg]
    frame_pmra(sol), frame_pmdec(sol)     [mas/yr]
    frame_rv(sol)                         [m/s]

Apparent position, proper motion, and radial velocity of the system
barycentre frame at this epoch, from rigorous 3D space-motion propagation of
the `AbsoluteFrame` catalog values. Compose with pairwise observables versus
`barycentre(sys)` to obtain absolute quantities of a body.
"""
frame_ra(sol::_AbsSol) = @inbounds sol.traj.ra2[sol.k]
frame_dec(sol::_AbsSol) = @inbounds sol.traj.dec2[sol.k]
frame_pmra(sol::_AbsSol) = @inbounds sol.traj.pmra2[sol.k]
frame_pmdec(sol::_AbsSol) = @inbounds sol.traj.pmdec2[sol.k]
frame_rv(sol::_AbsSol) = @inbounds sol.traj.rv2[sol.k]
export frame_ra, frame_dec, frame_pmra, frame_pmdec, frame_rv

# ---------------------------------------------------
# Name resolution: accept named `Body` values and `Symbol`s in either
# position. Both resolve type-stably through the trajectory's name table
# (constant-folding to a BodyRef when the names are compile-time constants),
# then dispatch to the AbstractRef methods above, which remain the hot-loop
# fast path.
# ---------------------------------------------------

const RefLike = Union{AbstractRef,Body,Symbol}

for fn in (:posx, :posy, :posz, :velx, :vely, :velz, :radvel,
           :raoff, :decoff, :pmra, :pmdec, :projectedseparation, :posangle)
    @eval @inline $fn(sol::TrajectorySolution, t::RefLike, r::RefLike) =
        $fn(sol, _resolve(_names(sol), t), _resolve(_names(sol), r))
end

# ---------------------------------------------------
# One-argument defaults: trivial two-body systems only
# ---------------------------------------------------

@inline function _default_pair(sol::TrajectorySolution)
    size(sol.traj.x, 2) == 2 || error(
        "one-argument observables are only defined for two-body systems; " *
        "pass an explicit target and reference — your Body values or their names, " *
        "e.g. raoff(sol, b, A) or raoff(sol, :b, :A)")
    return BodyRef(2), BodyRef(1)
end

for fn in (:posx, :posy, :posz, :velx, :vely, :velz, :radvel,
           :raoff, :decoff, :pmra, :pmdec, :projectedseparation, :posangle)
    @eval @inline function $fn(sol::TrajectorySolution)
        t, r = _default_pair(sol)
        return $fn(sol, t, r)
    end
end

export posx, posy, posz, velx, vely, velz, radvel,
       raoff, decoff, pmra, pmdec, projectedseparation, posangle
