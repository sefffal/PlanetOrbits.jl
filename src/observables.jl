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
    # The frame direction is the origin of the epoch triad by construction.
    @eval @inline $fn(sol::TrajectorySolution, ::FrameDirection) = zero(eltype(sol.traj.$col))
end

# Solution aliases by frame type. NB: written with the frame parameter in an
# invariant position — a `<:Trajectory{<:Any,FR}` pattern would put FR only in
# another typevar's bound, where Julia's dispatcher ignores it.
const _SolWithFrame{FR} = TrajectorySolution{Trajectory{T,FR,Names,E,V,M,FV}} where
    {T<:Number,Names,E<:AbstractVector{<:Real},V<:AbstractVector{T},
     M<:AbstractMatrix{T},FV<:AbstractVector{FR}}
const _AngularSol = Union{_SolWithFrame{<:Parallax},_SolWithFrame{<:AbsoluteFrame}}
const _AbsSol = _SolWithFrame{<:AbsoluteFrame}

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
    # The frame direction is the origin of the tangent plane and, by
    # definition, does not move within it: the frame's own drift is
    # `frame_pmra`/`frame_pmdec`, not a pairwise observable.
    @eval @inline $fn(sol::TrajectorySolution, ::FrameDirection) = zero(eltype(sol.traj.$cp))
    @eval @inline $dfn(sol::TrajectorySolution, ::FrameDirection) = zero(eltype(sol.traj.$cp))
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

"""
    posy(sol, target, reference)

Position offset [AU] of `target` relative to `reference` along the
declination (north) direction, on the same triad as [`posx`](@ref).
"""
@inline posy(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _sy(sol, t) - _sy(sol, r)

"""
    posz(sol, target, reference)

Position offset [AU] of `target` relative to `reference` along the line of
sight, positive away from the observer, on the same triad as [`posx`](@ref).
"""
@inline posz(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _sz(sol, t) - _sz(sol, r)

"""
    velx(sol, target, reference)

Velocity [AU / julian year] of `target` relative to `reference` along the
right-ascension (east) direction, on the same triad as [`posx`](@ref).
"""
@inline velx(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svx(sol, t) - _svx(sol, r)

"""
    vely(sol, target, reference)

Velocity [AU / julian year] of `target` relative to `reference` along the
declination (north) direction, on the same triad as [`posx`](@ref).
"""
@inline vely(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svy(sol, t) - _svy(sol, r)

"""
    velz(sol, target, reference)

Line-of-sight velocity [AU / julian year] of `target` relative to
`reference` — the **kinematic** quantity, for dynamics.

With the observing-geometry pass on, this already carries the projection onto
each body's own apparent direction (see `observe.jl`); what it does *not*
carry is the relativistic Einstein term. For the quantity a spectrograph
reports, in m/s, use [`radvel`](@ref). The distinction is kinematic
vs. spectroscopic, not coordinate vs. apparent.
"""
@inline velz(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) = _svz(sol, t) - _svz(sol, r)

# Einstein term [AU / julian year] of a reference point: a direct column read
# for a body, and zero for a barycentre — nothing emits from a barycentre, so
# `radvel(sol, A, barycentre(sys))` keeps the star's own term in full. A
# photocentre blends the light of its members, and so blends their terms.
@inline _ein(sol::TrajectorySolution, ref::BodyRef) =
    @inbounds sol.traj.ein[sol.k, ref.idx]
@inline function _ein(sol::TrajectorySolution, p::WeightedPoint{NB}) where {NB}
    m = sol.traj.ein
    acc = zero(eltype(m))
    p.emits || return acc
    k = sol.k
    @inbounds for j in 1:NB
        acc = muladd(p.w[j], m[k, j], acc)
    end
    return acc
end
@inline _ein(sol::TrajectorySolution, ::FrameDirection) = zero(eltype(sol.traj.ein))

"""
    radvel(sol, target, reference)

Radial velocity [m/s] of `target` relative to `reference` along the line of
sight, positive receding — the **spectroscopic** quantity, i.e. what a
spectrograph reports. E.g. `radvel(sol, b, A)` for a relative RV, or
`radvel(sol, A, barycentre(sys))` for the stellar reflex.

Two pieces: the kinematic projected line-of-sight velocity ([`velz`](@ref),
in physical units) and the **Einstein term** — the second-order Doppler and
gravitational-redshift difference between the two references,

```
Ein_i = ( ½|v_tot,i|² + Σ_{j≠i} G·mⱼ / r_ij ) / c
```

with `v_tot` the body's total barycentric velocity (orbital, plus the frame's
space velocity when the system has an `AbsoluteFrame`). Nothing emits from a
barycentre, so its Einstein term is zero and the stellar-reflex case carries
the star's own in full.

There is no keyword to decline this. The orbit-varying part depends on the
sampled orbit (e, masses, r(t)), so no reduction pipeline can have removed
it, and the constant part is absorbed by the instrument offset either way.
Its size, and which of the two uses of `radvel` it matters for, are tabulated
on the "Precision opt-outs" page — briefly, sub-cm/s for a stellar reflex
with a planetary companion, but several m/s of *variation* for the relative
RV of a close-in eccentric one.

Use [`velz`](@ref) instead when you want the kinematic velocity — for
dynamics, or to compare against a coordinate-velocity reference.

!!! note
    Masses therefore enter radial-velocity predictions, including their
    gradients. That is new in v2; see the migration guide.
"""
@inline radvel(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    (velz(sol, t, r) + (_ein(sol, t) - _ein(sol, r))) * au2m / year2sec_julian

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

"""
    pmdec(sol, target, reference)

Instantaneous relative proper motion [mas/julian year] of `target` with
respect to `reference` in declination — the exact time derivative of
[`decoff`](@ref). See [`pmra`](@ref).
"""
@inline pmdec(sol::TrajectorySolution, t::AbstractRef, r::AbstractRef) =
    _apmy(sol, t) - _apmy(sol, r)

# ---------------------------------------------------
# Observer-aware angular observables
#
# The observer is a property of each *observation*, not of the solve: one
# trajectory serves every likelihood in a model, and one model legitimately
# contains several observers (Hipparcos, Gaia at L2, a ground spectrograph).
# So the observer enters as an argument at read time and `orbitsolve` gains
# nothing — no ephemerides, no timescales, and no ambient observer state that
# could change an answer without appearing at the call site.
#
# The geometry is exact, not a series: each body's apparent direction is
# taken from the observer's actual position, so the annual–orbital (Kopeikin
# 1995) coupling and the exact per-body parallax factors fall out together.
# The tangent point is unchanged — the barycentre's apparent direction from
# the SSB — so these compose with the zero-argument forms.
# ---------------------------------------------------

# Observer position [AU, ICRS] resolved into this epoch's triad, which is the
# frame the stored body states are already in. R takes reference-triad
# components to epoch-triad components, and M1 takes reference-triad
# components to ICRS, so ICRS -> epoch triad is R * M1'.
@inline function _observer_offset(sol::_AbsSol, obs_pos)
    tr = sol.traj
    k = sol.k
    @inbounds r11 = tr.R11[k]
    isnan(_primal(r11)) && _err_observer_geometry()
    u = frame(sol).M1' * SVector(obs_pos[1], obs_pos[2], obs_pos[3])
    @inbounds begin
        ox = r11 * u[1] + tr.R12[k] * u[2] + tr.R13[k] * u[3]
        oy = tr.R21[k] * u[1] + tr.R22[k] * u[2] + tr.R23[k] * u[3]
        oz = tr.R31[k] * u[1] + tr.R32[k] * u[2] + tr.R33[k] * u[3]
    end
    return (ox, oy, oz)
end

@noinline _err_observer_geometry() = error(
    "observer-aware observables need a trajectory solved with " *
    "`observing_geometry=true`: the skipped path stores no per-epoch viewing " *
    "triad, and a µas-level observer coupling computed on top of the cheap " *
    "geometry would not be coherent anyway. Re-solve without " *
    "`observing_geometry=false`.")

@noinline _err_observer_frame() = error(
    "observer-aware observables need an `AbsoluteFrame`: an observer position " *
    "is given in ICRS, so placing it relative to the target requires the " *
    "target's ICRS direction. Build the `System` with ra, dec, plx, pmra, " *
    "pmdec, rv and ref_epoch.")

# Gnomonic coordinates as seen from `o`, per reference. Identical in form to
# `_ax`/`_ay` with the observer subtracted from both the transverse offset
# and the depth; at o = 0 they are the same expression.
for (fn, cp, co) in ((:_ax_obs, :x, 1), (:_ay_obs, :y, 2))
    @eval @inline function $fn(sol::TrajectorySolution, ref::BodyRef, o)
        tr = sol.traj
        k = sol.k
        j = ref.idx
        @inbounds (tr.$cp[k, j] - o[$co]) * rad2mas / (tr.d_au[k] + tr.z[k, j] - o[3])
    end
    @eval @inline function $fn(sol::TrajectorySolution, p::WeightedPoint{NB}, o) where {NB}
        tr = sol.traj
        k = sol.k
        acc = zero(eltype(tr.x))
        @inbounds for j in 1:NB
            acc = muladd(p.w[j],
                (tr.$cp[k, j] - o[$co]) * rad2mas / (tr.d_au[k] + tr.z[k, j] - o[3]), acc)
        end
        return acc
    end
    # A *direction* is unmoved by displacing the observer — it has no
    # parallax of its own. That asymmetry is the whole point of this
    # reference: against it, the target keeps its full parallax factor
    # (exactly, with no series expansion in ϖ) instead of having it cancelled
    # by a reference at the same distance.
    @eval @inline $fn(sol::TrajectorySolution, ::FrameDirection, o) = zero(eltype(sol.traj.x))
end

"""
    raoff(sol, target, reference, obs_pos)
    decoff(sol, target, reference, obs_pos)

Right-ascension / declination offset [mas] of `target` relative to
`reference` **as seen from `obs_pos`**, the observer's barycentric position
in ICRS Cartesian coordinates [AU] at this epoch (e.g. the Earth's, Gaia's at
L2, or `(0, 0, 0)` for the solar-system barycentre).

This is the seam for annual–orbital parallax (Kopeikin 1995) and exact
per-body parallax factors: the apparent direction of each body is computed
from the observer's actual position by the same exact geometry the
zero-argument forms use from the SSB, so the full coupling falls out with no
series expansion.

Which reference you name is what decides relative versus absolute, and there
is no flag. A body and a barycentre both sit at the system's distance, so
displacing the observer shifts target and reference together: the first-order
parallax cancels and only the *differential* (Kopeikin) part survives — that
is relative astrometry. [`framedirection`](@ref) is a direction, not a place,
and has no parallax of its own, so against it the target keeps its parallax
factor in full — that is absolute astrometry. Same code path, no special
case. (An absolute-astrometry likelihood that instead supplies its own
parallax term, as Gaia's published `parallax_factor_al` does, uses the
zero-argument forms against a barycentre and never passes `obs_pos` at all.)

Conventions: `obs_pos` is ICRS, in AU, and the epochs passed to `orbitsolve`
are barycentric (BJD\\_TDB-like MJD). Requires an `AbsoluteFrame` — an ICRS
observer position is meaningless without the target's ICRS direction — and a
trajectory solved with `observing_geometry=true`.

Ephemerides are deliberately not PlanetOrbits' business: the caller supplies
the position. A likelihood that needs one and has no ephemeris source should
say so when it is constructed, not degrade silently here.

    raoff(sol, b, A, earth_pos_au)              # relative: differential part only
    raoff(sol, A, framedirection, earth_pos_au) # absolute: the full parallax ellipse
    raoff(sol, A, barycentre(sys), (0,0,0)) == raoff(sol, A, barycentre(sys))
"""
@inline function raoff(sol::_AbsSol, t::AbstractRef, r::AbstractRef, obs_pos)
    o = _observer_offset(sol, obs_pos)
    return _ax_obs(sol, t, o) - _ax_obs(sol, r, o)
end

@inline function decoff(sol::_AbsSol, t::AbstractRef, r::AbstractRef, obs_pos)
    o = _observer_offset(sol, obs_pos)
    return _ay_obs(sol, t, o) - _ay_obs(sol, r, o)
end

"""
    projectedseparation(sol, target, reference, obs_pos)
    posangle(sol, target, reference, obs_pos)

Separation [mas] and position angle [rad] as seen from `obs_pos`. See the
four-argument [`raoff`](@ref).
"""
@inline function projectedseparation(sol::_AbsSol, t::AbstractRef, r::AbstractRef, obs_pos)
    o = _observer_offset(sol, obs_pos)
    x = _ax_obs(sol, t, o) - _ax_obs(sol, r, o)
    y = _ay_obs(sol, t, o) - _ay_obs(sol, r, o)
    return sqrt(x^2 + y^2)
end

@inline function posangle(sol::_AbsSol, t::AbstractRef, r::AbstractRef, obs_pos)
    o = _observer_offset(sol, obs_pos)
    return atan(_ax_obs(sol, t, o) - _ax_obs(sol, r, o),
                _ay_obs(sol, t, o) - _ay_obs(sol, r, o))
end

# Frames without an ICRS direction cannot place an ICRS observer at all.
for fn in (:raoff, :decoff, :projectedseparation, :posangle)
    @eval @noinline $fn(::TrajectorySolution, ::AbstractRef, ::AbstractRef, ::Any) =
        _err_observer_frame()
end

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
    frame_ra(sol)     [deg]

Apparent right ascension of the system barycentre frame at this epoch, from
rigorous 3D space-motion propagation of the `AbsoluteFrame` catalog values.

This and its four siblings — [`frame_dec`](@ref), [`frame_pmra`](@ref),
[`frame_pmdec`](@ref), [`frame_rv`](@ref) — describe the *frame*, not any one
body. Absolute quantities of a body compose as frame plus the pairwise
observable taken against `barycentre(sys)`, e.g.

    frame_pmra(sol) + pmra(sol, A, barycentre(sys))   # the star's absolute pmra

Requires a system built with a full absolute frame (`ra`, `dec`, `plx`,
`pmra`, `pmdec`, `rv`, `ref_epoch`).

`frame_ra` and `frame_dec` are computed on demand from the solved emission
epoch and the trajectory's frame rather than stored per epoch: they are the
only frame quantities requiring a transcendental, and nothing inside the
solver consumes them. Reading both costs ~33 ns/epoch; not storing them takes
`frame_pass!` from 76 to 40 ns/epoch for every model that never asks.
"""
@inline frame_ra(sol::_AbsSol) =
    _compensate_position(frame(sol), @inbounds sol.traj.t_em[sol.k]).ra2

"""
    frame_dec(sol)    [deg]

Apparent declination of the system barycentre frame at this epoch. Computed
on demand; see [`frame_ra`](@ref).
"""
@inline frame_dec(sol::_AbsSol) =
    _compensate_position(frame(sol), @inbounds sol.traj.t_em[sol.k]).dec2

"""
    frame_pmra(sol)   [mas / julian yr]

Apparent proper motion of the system barycentre frame in right ascension at
this epoch (already including the cos δ factor). Propagated, not frozen at
the catalog value, so perspective acceleration is present. See
[`frame_ra`](@ref).
"""
frame_pmra(sol::_AbsSol) = @inbounds sol.traj.pmra2[sol.k]

"""
    frame_pmdec(sol)  [mas / julian yr]

Apparent proper motion of the system barycentre frame in declination at this
epoch. See [`frame_ra`](@ref).
"""
frame_pmdec(sol::_AbsSol) = @inbounds sol.traj.pmdec2[sol.k]

"""
    frame_rv(sol)     [m/s]

Radial velocity of the system barycentre frame at this epoch, positive
receding — the propagated frame quantity, not a body's reflex. See
[`frame_ra`](@ref).
"""
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

# ... and for the observer-aware forms.
for fn in (:raoff, :decoff, :projectedseparation, :posangle)
    @eval @inline $fn(sol::TrajectorySolution, t::RefLike, r::RefLike, obs_pos) =
        $fn(sol, _resolve(_names(sol), t), _resolve(_names(sol), r), obs_pos)
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
