"""
    AbsoluteVisual{OrbitType}(..., ref_epoch=, ra=, dec=, plx=, rv=, pmra=, pmdec=)

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
rv        : m/s
ref_epoch : days

TODO: account for viewing angle differences and differential light travel
time between a planet and its host.
"""
struct AbsoluteVisualOrbit{T<:Number,O<:AbstractOrbit} <: AbstractOrbit{T}
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
distance(elem::AbsoluteVisualOrbit, t::Number) = elem.dist*au2pc

"""
    AbsoluteVisual{OrbitType}(..., ref_epoch=, ra=, dec=, plx=, rv=, pmra=, pmdec=)

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
rv        : m/s
ref_epoch : years

TODO: account for viewing angle differences and differential light travel
time between a planet and its host.
"""
const AbsoluteVisual{OrbitType} = AbsoluteVisualOrbit{T,OrbitType}  where T

function AbsoluteVisual{OrbitType}(;ref_epoch, ra, dec, plx, rv, pmra, pmdec, kargs...,) where {OrbitType}
    dist = 1000/plx * pc2au # distance [AU]
    parent = OrbitType(;kargs...)
    T = _parent_num_type(parent)
    T = promote_type(T, typeof(ref_epoch), typeof(ra), typeof(dec), typeof(plx), typeof(rv), typeof(pmra), typeof(pmdec))
    return AbsoluteVisualOrbit{T,OrbitType{T}}(parent, ref_epoch,  ra, dec, plx, rv, pmra, pmdec, dist)
end
function AbsoluteVisual(parent::AbstractOrbit, ref_epoch,  ra, dec, plx, rv, pmra, pmdec,)
    dist = 1000/plx * pc2au # distance [AU]
    T = _parent_num_type(parent)
    # TODO: we could have a conversion error here if the parent orbit uses a more restrictive number type and we cant convert these new properties to match
    return AbsoluteVisualOrbit{T,typeof(parent)}(parent, ref_epoch, ra, dec,  plx, rv, pmra, pmdec, dist)
end

export AbsoluteVisual

struct OrbitSolutionAbsoluteVisual{TEl<:AbstractOrbit,TSol<:AbstractOrbitSolution,T<:Number,TComp<:NamedTuple} <: AbstractOrbitSolution
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
function compensate_star_3d_motion(elem::AbsoluteVisualOrbit,epoch2_days::Number)
    ra1 = elem.ra             # degrees
    dec1 = elem.dec           # degrees
    parallax1 = elem.plx      # mas
    pmra1 = elem.pmra         # mas/yr
    pmdec1 = elem.pmdec       # mas/yr
    rv1 = elem.rv/1000        # m/s -> km/s
    epoch1_days = elem.ref_epoch # MJD

    # Be very careful now about assuming there are 365.245 days per year.
    # We want to use that value when using a number like proper motion in mas/yr.
    # We don't want to make that assumption when using delta years as a time unit,
    # since they are each either 365 or 366 days.

    # epoch1 = epoch1_days*days_per_average_year
    # epoch2 = epoch2_days*days_per_average_year

    # Guard against same epoch 
    # TODO: could just return arguments with appropriate units
    if epoch1_days == epoch2_days
        epoch2_days += eps(epoch2_days)
    end

    T = promote_type(
        typeof(ra1),
        typeof(dec1),
        typeof(parallax1),
        typeof(pmra1),
        typeof(pmdec1),
        typeof(rv1),
        typeof(epoch1_days),
        typeof(epoch2_days)
    )

    mydtor = convert(T, π / 180)
    my206265 = convert(T, 180 / π * 60 * 60)
    sec2year = convert(T, 365.25 * 24 * 60 * 60)
    pc2km = convert(T, 3.08567758149137e13)

    distance1 = convert(T, 1000 / parallax1)
    
    # convert RV to pc/year, convert delta RA and delta Dec to radians/year
    # These are differential quantities originally expressed per average-length-year.
    # We want them in units per day, which always have the same length
    dra1 = pmra1 / 1000 / my206265 / cosd(dec1)
    ddec1 = pmdec1 / 1000 /my206265
    ddist1 = rv1 / pc2km * sec2year

    # convert first epoch to x,y,z and dx,dy,dz
    sin_ra1, cos_ra1 = sincosd(ra1)
    sin_dec1, cos_dec1 = sincosd(dec1)

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

    
    # be careful here with units:
    delta_time_jyear = (epoch2_days - epoch1_days)/year2day_julian

    x₂ = x₁ + dx * delta_time_jyear#(epoch2-epoch1)
    y₂ = y₁ + dy * delta_time_jyear#(epoch2-epoch1)
    z₂ = z₁ + dz * delta_time_jyear#(epoch2-epoch1)

    # And done.  Now we just need to go backward.

    distance2 = sqrt(x₂^2 + y₂^2 + z₂^2)
    if distance2 == 0
        x₂=y₂=z₂=zero(x₂)
        distance2 += eps(distance2)
    end

    parallax2 = 1000/distance2

    ra2 = ((atand(y₂,x₂) + 360) % 360)
    arg = z₂ / distance2
    if 1.0 < arg < 1.0 + sqrt(eps(1.0))
        arg = one(arg)
    end
    dec2 = asind(arg)

    ddist2 = 1 / sqrt(x₂^2 + y₂^2 + z₂^2) * (x₂ * dx + y₂ * dy + z₂ * dz)
    dra2 = 1 / (x₂^2 + y₂^2) * (-1 * y₂ * dx + x₂ * dy)
    ddec2 = 1 / (distance2 * sqrt(1 - z₂^2 / distance2^2)) * (-1 * z₂ * ddist2 / distance2 + dz)

    pmra2 = dra2  * my206265 * 1000 * cosd(dec2)
    pmdec2 = ddec2 * 1000 * my206265
    rv2 = ddist2 * pc2km / sec2year

    # light travel time
    delta_time = (distance2 - distance1) * 3.085677e13 / 2.99792e5 # in seconds
    # epoch2a = epoch2 - delta_time/3.154e7
    epoch2a_days = epoch2_days - delta_time*sec2day

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
        rv2=rv2*1000,
        delta_time,
        epoch1_days,
        # epoch2,
        # epoch2a,
        epoch2a_days,
        x₁,
        y₁,
        z₁,
        x₂,
        y₂,
        z₂,
    )
end

# We have to override the generic `orbitsolve` for this case, as we have to adjust
# for light travel time here.
function orbitsolve(elem::AbsoluteVisualOrbit, t, method::AbstractSolver=Auto())
    
    # Epoch of periastron passage
    tₚ = periastron(elem)

    if t isa Integer
        t = float(t)
    end

    compensated = compensate_star_3d_motion(elem, t)

    # Mean anomaly
    MA = meanmotion(elem)/oftype(t, year2day_julian) * (compensated.epoch2a_days - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem), method)
    
    # Calculate true anomaly
    ν = _trueanom_from_eccanom(elem, EA)

    return orbitsolve_ν(elem, ν, EA, t, compensated)
end

function orbitsolve_meananom(elem::AbsoluteVisualOrbit, MA)
    
    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem))
    
    # Calculate true anomaly
    ν = 2*atan(elem.parent.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν, EA)
end

function orbitsolve_ν(
    elem::AbsoluteVisualOrbit,
    ν,
    # TODO: EA_from_ν is no longer accurate here, since the light travel
    # time can vary vs time, we can't determine time from nu directly.
    EA=EA_from_ν(elem.parent, ν),
    t=_time_from_EA(elem, EA),
    compensated::NamedTuple=compensate_star_3d_motion(elem,t);
    kwargs...
)
    # TODO: asking for a solution at a given ν is no longer well-defined,
    # as it will vary over time and not repeat over each orbital period.
    sol = orbitsolve_ν(elem.parent, ν, EA, compensated.epoch2a_days; kwargs...)
    return OrbitSolutionAbsoluteVisual(elem, sol, t, compensated)
end
# The solution time is the time we asked for, not the true time accounting for light travel.
soltime(os::OrbitSolutionAbsoluteVisual) = os.t

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
    @eval function ($fun)(os::OrbitSolutionAbsoluteVisual, args...)
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
    @eval function ($fun)(elem::AbsoluteVisualOrbit, args...)
        return ($fun)(elem.parent, args...)
    end
end

function radvel(os::OrbitSolutionAbsoluteVisual)
    # Adjust RV to account for star's 3D motion through space.
    # We add the difference between the RV at the reference epoch
    # and the RV at the measurement epoch
    return radvel(os.sol) + (os.compensated.rv2 - os.elem.rv)
end
function raoff(o::OrbitSolutionAbsoluteVisual)
    xcart = posx(o) # [AU]
    cart2angle = rad2as*oftype(xcart, 1e3)/o.compensated.distance2_pc
    xang = xcart*cart2angle # [mas]
    return xang 
end
function decoff(o::OrbitSolutionAbsoluteVisual)
    ycart = posy(o) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.compensated.distance2_pc
    yang = ycart*cart2angle # [mas]
    return yang
end
function pmra(o::OrbitSolutionAbsoluteVisual)
    ẋcart = o.elem.parent.J*(o.elem.parent.cosi_cosΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) - o.elem.parent.sinΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.compensated.distance2_pc
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang + (o.compensated.pmra2 - o.elem.pmra)
end
function pmdec(o::OrbitSolutionAbsoluteVisual)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.compensated.distance2_pc
    ẏang = ẏcart*cart2angle # [mas/year]
    return ẏang + (o.compensated.pmdec2 - o.elem.pmdec)
end

# The non-keplerian deviation due to system's 3D motion must be applied additively
# to both the
function radvel(o::OrbitSolutionAbsoluteVisual, M_planet)
    quantity = radvel(o.sol)
    M_tot = totalmass(o.elem)
    return -M_planet/M_tot*quantity + (o.compensated.rv2 - o.elem.rv)
end

function pmra(o::OrbitSolutionAbsoluteVisual, M_planet)
    ẋcart = o.elem.parent.J*(o.elem.parent.cosi_cosΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) - o.elem.parent.sinΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.compensated.distance2_pc
    quantity = ẋang = ẋcart*cart2angle # [mas/year]
    M_tot = totalmass(o.elem)
    return -M_planet/M_tot*quantity + (o.compensated.pmra2 - o.elem.pmra)
end
function pmdec(o::OrbitSolutionAbsoluteVisual, M_planet)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.compensated.distance2_pc
    quantity = ẏang = ẏcart*cart2angle # [mas/year]
    M_tot = totalmass(o.elem)
    return -M_planet/M_tot*quantity + (o.compensated.pmdec2 - o.elem.pmdec)
end

function accra(o::OrbitSolutionAbsoluteVisual)
    throw(NotImplementedException())
    # if eccentricity(o.elem) >= 1
    #     @warn "acceleration not tested for ecc >= 1 yet. Results are likely wrong."
    # end
    # ẍcart = -o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_cosΩ*o.sol.sinν_ω + o.elem.parent.sinΩ*o.sol.cosν_ω) # [AU/year^2]
    # cart2angle = rad2as*oftype(ẍcart, 1e3)/o.compensated.distance2_pc
    # ẍang = ẍcart*cart2angle # [mas/year^2] 
    # return ẍang
end
function accdec(o::OrbitSolutionAbsoluteVisual)
    # throw(NotImplementedException())
    # if eccentricity(o.elem) >= 1
    #     @warn "acceleration not tested for ecc >= 1 yet. Results are likely wrong."
    # end
    # ÿcart = o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_sinΩ*o.sol.sinν_ω - o.elem.parent.cosΩ*o.sol.cosν_ω) # [AU/year^2]
    # cart2angle = rad2as*oftype(ÿcart, 1e3)/o.compensated.distance2_pc
    # ÿang = ÿcart*cart2angle # [mas/year^2] 
    # return ÿang
end

"""
    PlanetOrbits.ra(orbit, t)

Get the instantaneous position of a companion in degrees of RA and Dec. 
For the relative position, see `raoff`.
"""
function ra(o::OrbitSolutionAbsoluteVisual, M_planet)
    # Already solved at correct epoch accoutning for light travel time
    # difference wrt. reference epoch.
    kep_offset_mas = raoff(o, M_planet)
    total = o.compensated.ra2 + kep_offset_mas/60/60/1000
    return total
end
"""
    PlanetOrbits.dec(orbit, t)

Get the instantaneous position of a companion in degrees of RA and Dec. 
For the relative position, see `decoff`.
"""
function dec(o::OrbitSolutionAbsoluteVisual, M_planet)
    # Already solved at correct epoch accoutning for light travel time
    # difference wrt. reference epoch.
    kep_offset_mas = decoff(o, M_planet)
    total = o.compensated.dec2 + kep_offset_mas/60/60/1000
    return total
end



# Pretty printing
function Base.show(io::IO, mime::MIME"text/plain", elem::AbsoluteVisual)
    show(io, mime, elem.parent)
    print(io, """\
    AbsoluteVisual
    ──────────────────────────
    reference epoch [days] = $(round(elem.ref_epoch, digits=1)) 
    plx [mas]      = $(round(elem.plx, digits=3)) 
    ra [°]         = $(round(elem.ra, digits=3)) 
    dec [°]        = $(round(elem.dec, digits=3)) 
    pmra [mas/yr]  = $(round(elem.pmra, digits=3)) 
    pmdec [mas/yr] = $(round(elem.pmdec, digits=3))
    rv [m/s]       = $(round(elem.rv, digits=3))
    ──────────────────────────
    """)
end
