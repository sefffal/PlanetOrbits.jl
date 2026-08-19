# ---------------------------------------------------
# System-level frame information
#
# Three levels, mirroring what v1 expressed as KepOrbit / Visual / AbsoluteVisual:
#   NoFrame       — physical-unit observables only (AU, AU/yr)
#   Parallax      — adds plx [mas]: angular observables in mas
#   AbsoluteFrame — full barycentre frame: ra/dec/plx/pm/rv at ref_epoch,
#                   with rigorous 3D space-motion propagation and
#                   light-travel-time compensation.
# ---------------------------------------------------

# A `Trajectory` carries its frame's *type* as a parameter, so observables
# dispatch on what frame information is available directly on `NoFrame` /
# `Parallax` / `AbsoluteFrame` — there is no parallel hierarchy of mode tags.
abstract type AbstractFrame end

struct NoFrame <: AbstractFrame end

struct Parallax{T<:Number} <: AbstractFrame
    plx::T          # mas
    cart2angle::T   # mas per AU at the system distance (== plx numerically)
end
function Parallax(plx)
    plx = float(plx)
    # angle [mas] = x [AU] / dist [pc] * 1000 = x * plx — exact by the
    # definitions of AU, pc, and arcsecond.
    return Parallax{typeof(plx)}(plx, plx)
end

"""
    AbsoluteFrame(; ra, dec, plx, pmra, pmdec, rv, ref_epoch)

Absolute barycentric frame of a `System`: ICRS position `ra`, `dec` [deg],
parallax `plx` [mas], proper motion `pmra` (μα*), `pmdec` [mas/yr], radial
velocity `rv` [m/s], all valid at `ref_epoch` [MJD].

The epoch-independent half of the rigorous 3D space-motion propagation
(the unit-sphere position and space-velocity vector) is precomputed at
construction, so per-epoch compensation is shared by every body and
observable — it is computed once per system per epoch, never per planet.
"""
struct AbsoluteFrame{T<:Number} <: AbstractFrame
    ra::T        # deg
    dec::T       # deg
    plx::T       # mas
    pmra::T      # mas/yr
    pmdec::T     # mas/yr
    rv::T        # m/s
    ref_epoch::T # MJD
    # hoisted propagation state
    distance1::T                # pc
    x1::T; y1::T; z1::T         # pc
    # Space velocity [pc / julian year]. In a frame built from catalog values
    # this is the *catalog-convention* velocity (apparent rates propagated as
    # coordinate rates — the light-time-free standard model); the light-time
    # solve swaps in the de-Dopplered true velocity via `_effective_frame`.
    dx::T; dy::T; dz::T
    # Local (East, North, line-of-sight) triad at the *reference* direction,
    # as columns in ICRS. Orbital elements are expressed in this triad, so it
    # is the frame the propagator's body states come out in; `observe_pass!`
    # rotates them into the triad of the epoch's apparent direction.
    M1::SMatrix{3,3,T,9}
end

function AbsoluteFrame(; ra, dec, plx, pmra, pmdec, rv, ref_epoch)
    ra, dec, plx, pmra, pmdec, rv, ref_epoch =
        promote(float(ra), dec, plx, pmra, pmdec, rv, ref_epoch)
    T = typeof(ra)
    # `_primal`: ForwardDiff's `iszero` compares the whole `Dual`
    # lexicographically, so `iszero(Dual(0.0, ∂))` is false and a
    # differentiated zero parallax would walk straight into `1000/plx`.
    if iszero(_primal(plx))
        error("starting parallax is zero -- can't propagate barycentric motion")
    end
    rv_kms = rv / 1000
    distance1 = 1000 / plx
    sin_ra1, cos_ra1 = sincosd(ra)
    sin_dec1, cos_dec1 = sincosd(dec)
    # Pole guard: a classification, so it is asked of the value. The clamped
    # replacement keeps the `Dual`'s sign but drops its partials, which is the
    # intent — `1/cos(dec)` has no meaningful derivative at the pole either.
    if abs(_primal(cos_dec1)) < 1e-15
        cos_dec1 = oftype(cos_dec1, copysign(1e-15, _primal(cos_dec1)))
    end
    x1 = cos_ra1 * cos_dec1 * distance1
    y1 = sin_ra1 * cos_dec1 * distance1
    z1 = sin_dec1 * distance1
    M1 = _triad(x1, y1, z1, distance1)
    dx, dy, dz = _frame_velocity(sin_ra1, cos_ra1, sin_dec1, cos_dec1,
                                 distance1, pmra, pmdec, rv_kms)
    return AbsoluteFrame{T}(ra, dec, plx, pmra, pmdec, rv, ref_epoch,
        distance1, x1, y1, z1, dx, dy, dz, M1)
end

# Cartesian space velocity [pc / julian year] from tangential rates and a
# radial velocity, at the (sin, cos) of the reference direction. Shared by the
# constructor (catalog-convention velocity) and `_dedoppler` (true velocity),
# so the two differ by exactly the Doppler factor applied to `pmra`/`pmdec`
# and by nothing else.
@inline function _frame_velocity(sin_ra1, cos_ra1, sin_dec1, cos_dec1,
                                 distance1, pmra, pmdec, rv_kms)
    dra1 = deg2rad(pmra / 1000 / 60 / 60) / cos_dec1
    ddec1 = deg2rad(pmdec / 1000 / 60 / 60)
    ddist1 = rv_kms * one_over_pc2km_sec2yr
    dx = -sin_ra1 * cos_dec1 * distance1 * dra1 - cos_ra1 * sin_dec1 * distance1 * ddec1 + cos_ra1 * cos_dec1 * ddist1
    dy = cos_ra1 * cos_dec1 * distance1 * dra1 - sin_ra1 * sin_dec1 * distance1 * ddec1 + sin_ra1 * cos_dec1 * ddist1
    dz = cos_dec1 * distance1 * ddec1 + sin_dec1 * ddist1
    return dx, dy, dz
end

"""
The frame with its space velocity converted from catalog convention to the
*true* (coordinate) velocity, for the barycentric light-travel solve.

Catalog proper motions (Hipparcos, Gaia) are **apparent** quantities: the
derivative of the light-time-affected direction with respect to the time of
light arrival, `d/dt_obs` (Butkevich & Lindegren 2014, Sect. 1; Gaia DR2/DR3
documentation Sect. 3.3.3). They already carry the factor relating apparent to
true tangential motion, `μ_app = μ_true / (1 + v_r/c)`. Building a worldline
directly from them and *then* solving the light time would count that factor
twice — a spurious `μ·v_r/c` slope (3.8 mas/yr at Barnard's star kinematics)
in every position readout. So the light-time path de-Dopplers first:
`μ_true = μ_app · (1 + v_r/c)` (B&L Eqs. 11–13; the same observed→inertial
step ERFA's `starpv` applies to catalog input). The radial rate is unchanged —
the spectroscopic/apparent distinction there is second order and far below
every target of this package.

The light-time-free path (`barycentric_lighttime=false`) must NOT use this:
propagating the catalog values linearly as they stand *is* the standard model
those catalogs were reduced with (ESA 1997 Sect. 1.5.5; Lindegren et al.
2012/2021), and data reduced with it must be propagated with it (B&L
Sects. 5.5, 6.1). `_effective_frame` selects per solve.
"""
function _dedoppler(fr::AbsoluteFrame{T}) where {T}
    doppler = 1 + fr.rv / c_light_ms
    sin_ra1, cos_ra1 = sincosd(fr.ra)
    sin_dec1, cos_dec1 = sincosd(fr.dec)
    # Same pole guard as the constructor, and asked of the value for the same
    # reason (743d2e1): a `Dual` comparison would classify on the partials.
    if abs(_primal(cos_dec1)) < 1e-15
        cos_dec1 = oftype(cos_dec1, copysign(1e-15, _primal(cos_dec1)))
    end
    dx, dy, dz = _frame_velocity(sin_ra1, cos_ra1, sin_dec1, cos_dec1,
                                 fr.distance1, fr.pmra * doppler,
                                 fr.pmdec * doppler, fr.rv / 1000)
    return AbsoluteFrame{T}(fr.ra, fr.dec, fr.plx, fr.pmra, fr.pmdec, fr.rv,
        fr.ref_epoch, fr.distance1, fr.x1, fr.y1, fr.z1, dx, dy, dz, fr.M1)
end

# The frame a solve actually propagates: the de-Dopplered (true-velocity)
# frame when the barycentric light-travel solve is on, the catalog-convention
# frame when it is off. `frame_pass!` stores the result in `traj.frame`, so
# every downstream consumer — the stored pm/rv columns, the on-demand
# `frame_ra`/`frame_dec`, and the observing-geometry passes — reads one
# consistent worldline.
@inline _effective_frame(fr::AbsoluteFrame, lighttime::Bool) =
    lighttime ? _dedoppler(fr) : fr
@inline _effective_frame(fr::AbstractFrame, ::Bool) = fr

"""
Local orthonormal (East, North, line-of-sight) triad at the ICRS direction of
`(x, y, z)`, returned as columns. Built arithmetically from the Cartesian
position — no trigonometry, and no round trip through (ra, dec).

    ê_r = r̂                      (away from the observer, +receding)
    ê_α = (-y, x, 0)/√(x²+y²)    (East, ∂/∂α)
    ê_δ = ê_r × ê_α              (North, ∂/∂δ)

(ê_α, ê_δ, ê_r) is right-handed, matching the package's (x=RA, y=Dec, z=LOS)
state convention.
"""
@inline function _triad(x, y, z, d)
    invd = inv(d)
    er = SVector(x * invd, y * invd, z * invd)
    # Plain sqrt, matching `_frame_geometry_pass!`: the components are parsecs,
    # so `hypot`'s overflow/underflow guarding is unreachable, and this is the
    # better AD path. `_triad` is construction-time, so this is about keeping
    # the reference bit-identical to the production pass, not about speed.
    ρ = sqrt(x * x + y * y)
    # On the pole the East direction is degenerate; any choice is as good as
    # any other and this one keeps the triad orthonormal and finite.
    ρ = ρ < eps(oneunit(ρ)) ? oneunit(ρ) * eps(oneunit(ρ)) : ρ
    invρ = inv(ρ)
    ea = SVector(-y * invρ, x * invρ, zero(x))
    ed = er × ea
    return hcat(ea, ed, er)
end

"""
Rigorous 3D space-motion propagation of the barycentre to `t_em_days`, in
parsecs. The whole frame block reduces to this plus algebra, which is why
`ra2`/`dec2` are recomputed on demand rather than stored (see `frame_ra`).
"""
@inline function _propagate_pc(fr::AbsoluteFrame, t_em_days)
    Δt_jyear = (t_em_days - fr.ref_epoch) / year2day_julian
    return (fr.x1 + fr.dx * Δt_jyear,
            fr.y1 + fr.dy * Δt_jyear,
            fr.z1 + fr.dz * Δt_jyear)
end

# At exactly `ref_epoch` the propagation is a no-op and v1 nudged off it;
# preserved so the fixtures stay bit-comparable. Applied by `compensate` and
# by the on-demand `frame_ra`/`frame_dec`, but *not* by `_geometry`, which
# must stay bit-identical to `_frame_geometry_pass!`.
@inline function _nudge_ref(fr::AbsoluteFrame, t_em_days)
    # On the primals: `==` between `Dual`s also compares partials, so an epoch
    # that *is* the reference epoch would skip the nudge on the gradient path
    # while taking it on the value path. Identical for Float64 arguments, so
    # the v1 fixtures stay bit-comparable.
    return _primal(fr.ref_epoch) == _primal(t_em_days) ?
        t_em_days + eps(float(_primal(t_em_days))) : t_em_days
end

# Guarded barycentre distance [pc] at `t_em_days`, with the propagated
# components. Plain sqrt rather than `hypot`: the components are parsecs, so
# the overflow/underflow range `hypot` guards is unreachable (the zero case is
# handled explicitly), and it is a better AD path.
@inline function _propagate_dist(fr::AbsoluteFrame, t_em_days)
    x2, y2, z2 = _propagate_pc(fr, t_em_days)
    distance2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    # `_primal`, for the reason given at the parallax guard: a differentiated
    # zero distance is a `Dual` that `iszero` calls nonzero — and `sqrt` at
    # zero has an infinite derivative, so it is exactly the case that needs
    # catching.
    if iszero(_primal(distance2))
        x2 = y2 = z2 = zero(x2)
        # Assigned, not incremented: `sqrt` at exactly zero has an infinite
        # derivative, so `distance2` arrives here as `Dual(0.0, ±Inf)` (or
        # `NaN`) and `+= eps` would keep those partials while the components
        # around it have just been zeroed. Same Float64 value as before.
        distance2 = oftype(distance2, eps(one(_primal(distance2))))
    end
    return x2, y2, z2, distance2
end

"""
The epoch at which light emitted at `t_em_days` is *received*, the only thing
the light-travel fixed point in `frame_pass!` consumes.

Light-travel time is taken *relative to the reference epoch*: the constant d/c
is degenerate with `tp` and the linear part with the period, but the curvature
— driven by the perspective acceleration v_tan²/d — is not. Taking the norm of
a linearly-moving 3D vector reproduces that curvature exactly; propagating
(ra, dec, plx) separately would not.

v1 (and v2 before `42ed5b7`) had this subtraction the other way round, which
inverts the sign of the whole barycentric light-travel correction: it made a
receding system's apparent period *shorter* than its true period rather than
longer. See the sign test in `test/runtests.jl`.
"""
@inline function _received_epoch(fr::AbsoluteFrame, t_em_days)
    _, _, _, distance2 = _propagate_dist(fr, t_em_days)
    return t_em_days + (distance2 - fr.distance1) * pc2sec_light * sec2day
end

"""
Proper motion and radial velocity of the propagated frame at `t_em_days`.

Split out of `compensate` because it needs **no transcendental function at
all** — `ra2`/`dec2` are the only ones that do, and they are computed on
demand by `frame_ra`/`frame_dec`. The split has to be structural rather than
left to the compiler: `atand` and `asind` carry `throw` paths, so LLVM cannot
dead-code-eliminate them even when their results are unused (the same trap as
the domain-error branches in `kepsolve-simd.jl`). Measured: 70.5 → 39.6
ns/epoch on `frame_pass!`.

The trigonometric identities that remove a `hypot`, a `cos` and a duplicated
`sqrt` are noted inline; they agree with a literal transcription of v1's
`compensate_star_3d_motion` to ≤1.3e-15 relative, gated in
`test/observing-geometry.jl`.
"""
@inline function _compensate_kinematics(fr::AbsoluteFrame, t_em_days)
    t_em_days = _nudge_ref(fr, t_em_days)
    x2, y2, z2, distance2 = _propagate_dist(fr, t_em_days)
    invd = inv(distance2)
    # sin δ and cos δ algebraically. δ ∈ [-90°, 90°] so cos δ ≥ 0, hence
    # cos δ = √(1 − sin²δ) — and that single sqrt is simultaneously the
    # `cosd(dec2)` factor in μα* and the √(1 − z²/d²) denominator in ∂δ/∂t,
    # which were previously a `cos` and a `sqrt` computed independently.
    sindec2 = clamp(z2 * invd, -one(z2), one(z2))
    cosdec2 = sqrt(1 - sindec2 * sindec2)
    ddist2 = (x2 * fr.dx + y2 * fr.dy + z2 * fr.dz) * invd
    dra2 = (-y2 * fr.dx + x2 * fr.dy) / (x2^2 + y2^2)
    ddec2 = (-z2 * ddist2 * invd + fr.dz) * invd / cosdec2
    pmra2 = dra2 * rad2as_206265 * 1000 * cosdec2
    pmdec2 = ddec2 * 1000 * rad2as_206265
    rv2 = ddist2 * pc2km / year2sec_julian * 1e3   # m/s
    return (; pmra2, pmdec2, rv2)
end

"""
Apparent (ra, dec) [deg] of the propagated frame at `t_em_days` — the two
quantities needing a transcendental. Computed on demand by `frame_ra` /
`frame_dec` rather than stored per epoch.
"""
@inline function _compensate_position(fr::AbsoluteFrame, t_em_days)
    t_em_days = _nudge_ref(fr, t_em_days)
    x2, y2, z2, distance2 = _propagate_dist(fr, t_em_days)
    ra2 = (atand(y2, x2) + 360) % 360
    dec2 = asind(clamp(z2 / distance2, -one(z2), one(z2)))
    return (; ra2, dec2)
end

"""
Per-epoch half of the 3D space-motion compensation, complete. Algebraically
identical to PlanetOrbits v1's `compensate_star_3d_motion`, with the same
constants (verified to ≤3e-14), minus the setup hoisted into `AbsoluteFrame`.

Production splits this three ways — `_received_epoch` for the light-travel
fixed point, `_compensate_kinematics` for the stored pm/rv columns, and
`_compensate_position` for on-demand (ra, dec) — so that no path pays for
quantities it does not use. This assembled form is kept as the readable
reference the test suite gates those three against.
"""
@inline function compensate(fr::AbsoluteFrame, t_em_days)
    t_em_days = _nudge_ref(fr, t_em_days)
    x2, y2, z2, distance2 = _propagate_dist(fr, t_em_days)
    kin = _compensate_kinematics(fr, t_em_days)
    pos = _compensate_position(fr, t_em_days)
    return (; pos.ra2, pos.dec2, kin.pmra2, kin.pmdec2, kin.rv2,
            epoch2a_days=_received_epoch(fr, t_em_days), x2, y2, z2, distance2)
end

# ---------------------------------------------------
# Per-epoch observing geometry, consumed by `observe_pass!`.
#
# Returns (R, d_au, V):
#   R      — rotation taking body-state components from the reference triad
#            (in which the elements are expressed) to the triad of the
#            barycentre's *apparent* direction at this epoch.
#   d_au   — barycentre distance [AU] at this epoch, for the AU -> mas scale.
#   V      — barycentre space velocity [AU / julian year] resolved in that
#            same triad. Its transverse components are what makes two bodies
#            at different sky positions have different radial velocities.
# ---------------------------------------------------

@inline function _geometry(fr::AbsoluteFrame, t_em, ::Type{T}) where {T}
    x2, y2, z2 = _propagate_pc(fr, t_em)
    # sqrt, not hypot: see `_triad`. This is the reference `_frame_geometry_pass!`
    # is gated against, so it must use the same arithmetic the pass does.
    d2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    M2 = _triad(x2, y2, z2, d2)
    R = M2' * fr.M1
    V = (M2' * SVector(fr.dx, fr.dy, fr.dz)) * pc2au   # AU / julian year
    return (R, d2 * pc2au, V)
end

# Frame construction from System keyword arguments.
function _make_frame(; plx=nothing, ra=nothing, dec=nothing, pmra=nothing,
                       pmdec=nothing, rv=nothing, ref_epoch=nothing)
    absolute_args = (ra, dec, pmra, pmdec, rv, ref_epoch)
    if any(!isnothing, absolute_args)
        if any(isnothing, absolute_args) || isnothing(plx)
            error("an absolute frame requires all of ra, dec, plx, pmra, pmdec, rv, and ref_epoch")
        end
        return AbsoluteFrame(; ra, dec, plx, pmra, pmdec, rv, ref_epoch)
    elseif !isnothing(plx)
        return Parallax(plx)
    else
        return NoFrame()
    end
end

distance(fr::Parallax) = 1000 / fr.plx        # pc
distance(fr::AbsoluteFrame) = fr.distance1    # pc
