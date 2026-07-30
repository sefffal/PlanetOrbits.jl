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

abstract type AbstractFrame end

# Frame *modes*: plain singleton tags carried in a Trajectory's type so
# observables can dispatch on what frame information is available.
abstract type FrameMode end
struct ModeNone <: FrameMode end
struct ModeParallax <: FrameMode end
struct ModeAbsolute <: FrameMode end

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
    dx::T; dy::T; dz::T         # pc / julian year
end

function AbsoluteFrame(; ra, dec, plx, pmra, pmdec, rv, ref_epoch)
    ra, dec, plx, pmra, pmdec, rv, ref_epoch =
        promote(float(ra), dec, plx, pmra, pmdec, rv, ref_epoch)
    T = typeof(ra)
    if iszero(plx)
        error("starting parallax is zero -- can't propagate barycentric motion")
    end
    rv_kms = rv / 1000
    distance1 = 1000 / plx
    sin_ra1, cos_ra1 = sincosd(ra)
    sin_dec1, cos_dec1 = sincosd(dec)
    if abs(cos_dec1) < 1e-15
        cos_dec1 = copysign(1e-15, cos_dec1)
    end
    dra1 = deg2rad(pmra / 1000 / 60 / 60) / cos_dec1
    ddec1 = deg2rad(pmdec / 1000 / 60 / 60)
    ddist1 = rv_kms * one_over_pc2km_sec2yr
    x1 = cos_ra1 * cos_dec1 * distance1
    y1 = sin_ra1 * cos_dec1 * distance1
    z1 = sin_dec1 * distance1
    dx = -sin_ra1 * cos_dec1 * distance1 * dra1 - cos_ra1 * sin_dec1 * distance1 * ddec1 + cos_ra1 * cos_dec1 * ddist1
    dy = cos_ra1 * cos_dec1 * distance1 * dra1 - sin_ra1 * sin_dec1 * distance1 * ddec1 + sin_ra1 * cos_dec1 * ddist1
    dz = cos_dec1 * distance1 * ddec1 + sin_dec1 * ddist1
    return AbsoluteFrame{T}(ra, dec, plx, pmra, pmdec, rv, ref_epoch,
        distance1, x1, y1, z1, dx, dy, dz)
end

"""
Per-epoch half of the 3D space-motion compensation. Identical formulas and
constants to PlanetOrbits v1's `compensate_star_3d_motion` (verified to
≤3e-14), minus the setup hoisted into `AbsoluteFrame`, returning only the
fields observables actually consume.
"""
@inline function compensate(fr::AbsoluteFrame, t_em_days)
    if fr.ref_epoch == t_em_days
        t_em_days += eps(float(t_em_days))
    end
    Δt_jyear = (t_em_days - fr.ref_epoch) / year2day_julian
    x2 = fr.x1 + fr.dx * Δt_jyear
    y2 = fr.y1 + fr.dy * Δt_jyear
    z2 = fr.z1 + fr.dz * Δt_jyear
    distance2 = hypot(x2, y2, z2)
    if iszero(distance2)
        x2 = y2 = z2 = zero(x2)
        distance2 += eps(one(distance2))
    end
    ra2 = (atand(y2, x2) + 360) % 360
    arg = z2 / distance2
    arg = clamp(arg, -one(arg), one(arg))
    dec2 = asind(arg)
    ddist2 = (x2 * fr.dx + y2 * fr.dy + z2 * fr.dz) / distance2
    dra2 = (-y2 * fr.dx + x2 * fr.dy) / (x2^2 + y2^2)
    ddec2 = (-z2 * ddist2 / distance2 + fr.dz) / (distance2 * sqrt(1 - z2^2 / distance2^2))
    pmra2 = dra2 * rad2as_206265 * 1000 * cosd(dec2)
    pmdec2 = ddec2 * 1000 * rad2as_206265
    rv2 = ddist2 * pc2km / year2sec_julian * 1e3   # m/s
    delta_time = (distance2 - fr.distance1) * c_light_pc / 2.99792e5  # s
    epoch2a_days = t_em_days - delta_time * sec2day
    cart2angle = rad2as * 1e3 / (distance2 * pc2au)  # AU -> mas at current distance
    return (; cart2angle, ra2, dec2, pmra2, pmdec2, rv2, epoch2a_days)
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
