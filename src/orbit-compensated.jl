"""
    Compensated{OrbitType}(..., ref_epoch=, ra=, dec=, plx=, rv=, pmra=, pmdec=)

This wraps another orbit object to add parallax, proper motion, and
RV fields, at a given reference epoch. 

Like a Visual{OrbitType} this allows for calculating projected quantities,
eg. separation in milliarcseconds.

What this type additionally does is correct for the star's 3D motion through
space (RV and proper motion) and differential light travel-time compared to a
reference epoch when calculating various quantities. 
This becomes necessary when computing eg. RVs over a long time period.

ra        : degrees
dec       : degrees
plx       : mas
pmra      : mas/yr
pmdec     : mas/yr
rv        : km/s
ref_epoch : years

TODO: account for viewing angle differences and differential light travel
time between a planet and its host.
"""
struct CompensatedOrbit{T<:Number,O<:AbstractOrbit} <: AbstractOrbit{T}
    parent::O
    ref_epoch::T
    ra::T
    dec::T
    plx::T
    rv::T
    pmra::T
    pmdec::T
    dist::T
end
# TODO: distance vs time
distance(elem::CompensatedOrbit, t::Number) = elem.dist*au2pc

"""
    Compensated{OrbitType}(..., ref_epoch=, ra=, dec=, plx=, rv=, pmra=, pmdec=)

This wraps another orbit object to add parallax, proper motion, and
RV fields, at a given reference epoch. 

Like a Visual{OrbitType} this allows for calculating projected quantities,
eg. separation in milliarcseconds.

What this type additionally does is correct for the star's 3D motion through
space (RV and proper motion) and differential light travel-time compared to a
reference epoch when calculating various quantities. 
This becomes necessary when computing eg. RVs over a long time period.

ra        : degrees
dec       : degrees
parallax  : mas
pmra      : mas/yr
pmdec     : mas/yr
rv        : km/s
ref_epoch : years

TODO: account for viewing angle differences and differential light travel
time between a planet and its host.
"""
const Compensated{OrbitType} = CompensatedOrbit{T,OrbitType}  where T

function Compensated{OrbitType}(;ref_epoch, ra, dec, plx, rv, pmra, pmdec, kargs...,) where {OrbitType}
    dist = 1000/plx * pc2au # distance [AU]
    parent = OrbitType(;kargs...)
    T = _parent_num_type(parent)
    return CompensatedOrbit{T,OrbitType{T}}(parent, ref_epoch,  ra, dec, plx, rv, pmra, pmdec, dist)
end
function Compensated(parent::AbstractOrbit, ref_epoch,  ra, dec, plx, rv, pmra, pmdec,)
    dist = 1000/plx * pc2au # distance [AU]
    T = _parent_num_type(parent)
    return CompensatedOrbit{T,typeof(parent)}(parent, ref_epoch, ra, dec,  plx, rv, pmra, pmdec, dist)
end

export Compensated

struct OrbitSolutionCompensated{TEl<:AbstractOrbit,TSol<:AbstractOrbitSolution,T<:Number,TComp<:NamedTuple} <: AbstractOrbitSolution
    elem::TEl
    sol::TSol
    t::T
    compensated::TComp
end

"""
This function calculates how to account for stellar 3D motion 
when comparing measurements across epochs (epoch1 vs epoch2).

Typically `epoch1` is your reference epoch, `epoch2` is your measurement
epoch, and the remaining parameters are parameters you are hoping to fit.
You use this function to calculate their compensated values, and compare
these to data at `epoch2`.

Will also calculates light travel time, returning updated epochs
(epoch2a) due to change in distance between epoch1 and epoch2.
epoch2 will be when the light was detected, epoch2a will be the
"emitted" time accounting for the different positions between epoch1
and epoch 2.

Original Author: Eric Nielsen
"""
function compensate_star_3d_motion(elem::CompensatedOrbit,epoch2::Number #= years =#)
    ra1 = elem.ra             # degrees
    dec1 = elem.dec           # degrees
    parallax1 = elem.plx      # mas
    pmra1 = elem.pmra         # mas/yr
    pmdec1 = elem.pmdec       # mas/yr
    rv1 = elem.rv             # km/s
    epoch1 = elem.ref_epoch   # years

    T = promote_type(
        typeof(ra1),
        typeof(dec1),
        typeof(parallax1),
        typeof(pmra1),
        typeof(pmdec1),
        typeof(rv1),
        typeof(epoch1),
        typeof(epoch2)
    )

    mydtor = convert(T, π / 180)
    my206265 = convert(T, 180 / π * 60 * 60)
    sec2year = convert(T, 365.25 * 24 * 60 * 60)
    pc2km = convert(T, 3.08567758149137e13)

    distance1 = convert(T, 1000 / parallax1)
    
    # convert RV to pc/year, convert delta RA and delta Dec to radians/year
    dra1 = pmra1 / 1000 / my206265 / cos(dec1 * mydtor)
    ddec1 = pmdec1 / 1000 /my206265
    ddist1 = rv1 / pc2km * sec2year

    # convert first epoch to x,y,z and dx,dy,dz
    sin_ra1, cos_ra1 = sincos(ra1*mydtor)
    sin_dec1, cos_dec1 = sincos(dec1*mydtor)

    x₁ = cos_ra1 * cos_dec1 * distance1
    y₁ = sin_ra1 * cos_dec1 * distance1
    z₁ = sin_dec1 * distance1

    # Excellent.  Now dx,dy,dz,which are constants

    dx = -1 * sin_ra1 * cos_dec1 * distance1 * dra1 -
        cos_ra1 * sin_dec1 * distance1 * ddec1 +
        cos_ra1 * cos_dec1 * ddist1 

    dy = 1  * cos_ra1 * cos_dec1 * distance1 * dra1 -
        sin_ra1 * sin_dec1 * distance1 * ddec1 +
        sin_ra1 * cos_dec1 * ddist1

    dz = 1 * cos_dec1 * distance1 * ddec1 + sin_dec1 * ddist1

    # Now the simple part:

    x₂ = x₁ + dx * (epoch2-epoch1)
    y₂ = y₁ + dy * (epoch2-epoch1)
    z₂ = z₁ + dz * (epoch2-epoch1)

    # And done.  Now we just need to go backward.

    distance2 = sqrt(x₂^2 + y₂^2 + z₂^2)

    parallax2 = 1000/distance2

    ra2 = ((atan(y₂,x₂)/mydtor + 360) % 360)
    dec2 = asin(z₂ / distance2) / mydtor

    ddist2 = 1 / sqrt(x₂^2 + y₂^2 + z₂^2) * (x₂ * dx + y₂ * dy + z₂ * dz)
    dra2 = 1 / (x₂^2 + y₂^2) * (-1 * y₂ * dx + x₂ * dy)
    ddec2 = 1 / (distance2 * sqrt(1 - z₂^2 / distance2^2)) * (-1 * z₂ * ddist2 / distance2 + dz)

    pmra2 = dra2  * my206265 * 1000 * cos(dec2 * mydtor)
    pmdec2 = ddec2 * 1000 * my206265
    rv2 = ddist2 * pc2km / sec2year


    # dra1 = pmra1 / 1000d0 * distance1
    # ddec1 = pmdec1 / 1000d0 * distance1
    # ddist1 = rv1 * 3.24078e-14 * 3.15e7

    # stop

    # light travel time

    delta_time = (distance2 - distance1) * 3.085677e13 / 2.99792e5 # in seconds
    epoch2a = epoch2 - delta_time/3.154e7

    distance2_pc = distance2 * pc2au

    return (;
        distance2_pc,
        parallax2,
        ra2,
        dec2,
        ddist2,
        dra2,
        ddec2,
        pmra2,
        pmdec2,
        rv2,
        delta_time,
        epoch2a,
    )
end

# We have to override the generic `orbitsolve` for this case, as we have to adjust
# for light travel time here.
function orbitsolve(elem::CompensatedOrbit, t, method::AbstractSolver=Auto())
    
    # Epoch of periastron passage
    tₚ = periastron(elem)

    if t isa Integer
        t = float(t)
    end

    compensated = compensate_star_3d_motion(elem, t)

    # Mean anomaly
    MA = meanmotion(elem)/oftype(t, year2day) * (compensated.epoch2a - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem), method)
    
    # Calculate true anomaly
    ν = _trueanom_from_eccanom(elem, EA)

    return orbitsolve_ν(elem, ν, EA, t, compensated)
end

function orbitsolve_ν(elem::CompensatedOrbit, ν, EA, t, compensated::NamedTuple; kwargs...)
    # TODO: asking for a solution at a given ν is no longer well-defined,
    # as it will vary over time and not repeat over each orbital period.
    sol = orbitsolve_ν(elem.parent, ν, EA, compensated.epoch2a; kwargs...)
    return OrbitSolutionCompensated(elem, sol, t, compensated)
end
# The solution time is the time we asked for, not the true time accounting for light travel.
soltime(os::OrbitSolutionCompensated) = os.t

# Forward these functions to the underlying orbit object
solution_fun_list = (
    :trueanom,
    :eccanom,
    :meananom,
    :posx,
    :posy,
    :posz,
    :posangle,
    :velx,
    :vely,
    :velz,
)
for fun in solution_fun_list
    # TODO-1: several of these need to handle the varying parallax correctly
    # TODO-2: several more need to account for chaning viewing angle and planet light-travel time.
    @eval function ($fun)(os::OrbitSolutionCompensated, args...)
        return ($fun)(os.sol, args...)
    end
end
orbit_fun_list = (
    :eccentricity,
    :periastron,
    :period,
    :inclination,
    :semimajoraxis,
    :totalmass,
    :meanmotion,
    :semiamplitude,
    :_trueanom_from_eccanom,
)
for fun in orbit_fun_list
    @eval function ($fun)(elem::CompensatedOrbit, args...)
        return ($fun)(elem.parent, args...)
    end
end

function radvel(os::OrbitSolutionCompensated, args...)
    # Adjust RV to account for star's 3D motion through space.
    # We add the difference between the RV at the reference epoch
    # and the RV at the measurement epoch
    return radvel(os.sol, args...) + (os.elem.rv - os.compensated.rv2)
end
# TODO
function raoff(o::OrbitSolutionCompensated)
    xcart = posx(o) # [AU]
    cart2angle = rad2as*oftype(xcart, 1e3)/o.compensated.distance2_pc
    xang = xcart*cart2angle # [mas]
    return xang
end
# TODO
function decoff(o::OrbitSolutionCompensated)
    ycart = posy(o) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.compensated.distance2_pc
    yang = ycart*cart2angle # [mas]
    return yang
end
# TODO
function pmra(o::OrbitSolutionCompensated)
    ẋcart = o.elem.parent.J*(o.elem.parent.cosi_cosΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) - o.elem.parent.sinΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.compensated.distance2_pc
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang
end
# TODO
function pmdec(o::OrbitSolutionCompensated)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.compensated.distance2_pc
    ẏang = ẏcart*cart2angle # [mas/year]
    return ẏang
end
# TODO
function accra(o::OrbitSolutionCompensated)
    if eccentricity(o.elem) >= 1
        @warn "acceleration not tested for ecc >= 1 yet. Results are likely wrong."
    end
    ẍcart = -o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_cosΩ*o.sol.sinν_ω + o.elem.parent.sinΩ*o.sol.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ẍcart, 1e3)/o.compensated.distance2_pc
    ẍang = ẍcart*cart2angle # [mas/year^2] 
    return ẍang
end
# TODO
function accdec(o::OrbitSolutionCompensated)
    if eccentricity(o.elem) >= 1
        @warn "acceleration not tested for ecc >= 1 yet. Results are likely wrong."
    end
    ÿcart = o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_sinΩ*o.sol.sinν_ω - o.elem.parent.cosΩ*o.sol.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ÿcart, 1e3)/o.compensated.distance2_pc
    ÿang = ÿcart*cart2angle # [mas/year^2] 
    return ÿang
end



# Pretty printing
function Base.show(io::IO, mime::MIME"text/plain", elem::Compensated)
    show(io, mime, elem.parent)
    print(io, """\
    plx [mas] = $(round(elem.plx, sigdigits=3)) 
    distance    [pc  ] : $(round(distance(elem), digits=1)) 
    ──────────────────────────
    """)
end
