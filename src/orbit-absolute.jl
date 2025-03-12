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

const c_light_ms = 2.998e+8
const c_light_pc = 3.085677e13


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
function compensate_star_3d_motion(elem::AbsoluteVisualOrbit, t_em_days::Number)
    ra1 = elem.ra             # degrees
    dec1 = elem.dec           # degrees
    parallax1 = elem.plx      # mas
    pmra1 = elem.pmra         # mas/yr
    pmdec1 = elem.pmdec       # mas/yr
    rv1 = elem.rv/1000        # m/s -> km/s
    epoch1_days = elem.ref_epoch # MJD

    # Guard against same epoch 
    # TODO: could just return arguments with appropriate units
    if epoch1_days == t_em_days
        t_em_days += eps(float(t_em_days))
    end

    T = promote_type(
        typeof(ra1),
        typeof(dec1),
        typeof(parallax1),
        typeof(pmra1),
        typeof(pmdec1),
        typeof(rv1),
        typeof(epoch1_days),
        typeof(t_em_days)
    )

    my206265 = convert(T, 180 / π * 60 * 60)
    sec2year = convert(T, 365.25 * 24 * 60 * 60) # Julia years
    pc2km = convert(T, 3.08567758149137e13)
    one_over_pc2km_sec2yr = 1.022712165045694034700736065713114217745793404987068055763763987835564887975633e-06

    distance1 = convert(T, 1000 / parallax1)

    if parallax1 == 0
        error("assertion error: starting parallax is zero -- can't propagate barycentric motion")
    end
    
    # convert RV to pc/year, convert delta RA and delta Dec to radians/year
    # These are differential quantities originally expressed per average-length-year.
    # We want them in units per day, which always have the same length

    # convert first epoch to x,y,z and dx,dy,dz
    sin_ra1, cos_ra1 = sincosd(ra1)
    sin_dec1, cos_dec1 = sincosd(dec1)
    if abs(cos_dec1) < 1e-15
        cos_dec1 = copysign(1e-15, cos_dec1)
    end

    dra1 = deg2rad(pmra1 / 1000 / 60 / 60) / cos_dec1
    ddec1 = deg2rad(pmdec1 / 1000 / 60 / 60)
    ddist1 = rv1 * one_over_pc2km_sec2yr#/ pc2km * sec2year

    x₁ = cos_ra1 * cos_dec1 * distance1
    y₁ = sin_ra1 * cos_dec1 * distance1
    z₁ = sin_dec1 * distance1

    # Now propagate through space linearly
    # We use compensated summation
    # dx = sum(@SVector [-1 * sin_ra1 * cos_dec1 * distance1 * dra1,
    #     -cos_ra1 * sin_dec1 * distance1 * ddec1,
    #     cos_ra1 * cos_dec1 * ddist1
    # ])
    dx = (-1 * sin_ra1 * cos_dec1 * distance1 * dra1 
        -cos_ra1 * sin_dec1 * distance1 * ddec1
        +cos_ra1 * cos_dec1 * ddist1
    )

    # dy = sum(@SVector[1  * cos_ra1 * cos_dec1 * distance1 * dra1,
    #     -sin_ra1 * sin_dec1 * distance1 * ddec1,
    #     sin_ra1 * cos_dec1 * ddist1,
    # ])
    dy = (1  * cos_ra1 * cos_dec1 * distance1 * dra1 
        -sin_ra1 * sin_dec1 * distance1 * ddec1
        +sin_ra1 * cos_dec1 * ddist1
    )

    dz = 1 * cos_dec1 * distance1 * ddec1 + sin_dec1 * ddist1

    
    # be careful here with units:
    # This is the change in time between requested epoch and reference epoch
    delta_time_jyear = (t_em_days - epoch1_days)/year2day_julian

    x₂ = x₁ + dx * delta_time_jyear#(epoch2-epoch1)
    y₂ = y₁ + dy * delta_time_jyear#(epoch2-epoch1)
    z₂ = z₁ + dz * delta_time_jyear#(epoch2-epoch1)

    # And done.  Now we just need to go backward.
    distance2 = hypot(x₂, y₂, z₂)
    if distance2 == 0
        x₂=y₂=z₂=zero(x₂)
        distance2 += eps(distance2)
    end

    parallax2 = 1000/distance2

    ra2 = ((atand(y₂,x₂) + 360) % 360)
    arg = z₂ / distance2
    arg = clamp(arg, -1.0, 1.0)
    dec2 = asind(arg)

    ddist2 = 1 / sqrt(x₂^2 + y₂^2 + z₂^2) * (x₂ * dx + y₂ * dy + z₂ * dz)
    dra2 = 1 / (x₂^2 + y₂^2) * (-1 * y₂ * dx + x₂ * dy)
    ddec2 = 1 / (distance2 * sqrt(1 - z₂^2 / distance2^2)) * (-1 * z₂ * ddist2 / distance2 + dz)

    pmra2 = dra2  * my206265 * 1000 * cosd(dec2)
    pmdec2 = ddec2 * 1000 * my206265
    rv2 = ddist2 * pc2km / sec2year

    # This is the light travel delta between the reference position and the solved position
    delta_time = (distance2 - distance1) * c_light_pc / 2.99792e5 # in seconds
    epoch2a_days = t_em_days - delta_time*sec2day

    distance2_pc = distance2

    return (;
        distance2_au=distance2_pc*pc2au,
        distance2_pc=distance2_pc,
        parallax1=parallax1,
        parallax2=parallax2,
        ra1=ra1,
        dec1=dec1,
        ra2=ra2,
        dec2=dec2,
        ddist2=ddist2,
        dra2=dra2,
        ddec2=ddec2,
        pmra2=pmra2,
        pmdec2=pmdec2,
        rv2=rv2*1e3,
        delta_time=delta_time,
        epoch1_days=epoch1_days,
        # epoch2,
        # epoch2a,
        t_em_days=t_em_days,
        epoch2a_days=epoch2a_days,
        x₁=x₁,
        y₁=y₁,
        z₁=z₁,
        x₂=x₂,
        y₂=y₂,
        z₂=z₂,
        dx=dx,
        dy=dy,
        dz=dz
    )
end


function _calculate_orbit_solution(elem, t_em, compensated, method)
    tₚ = periastron(elem)
    MA = meanmotion(elem)/year2day_julian * (t_em - tₚ)
    EA = kepler_solver(MA, eccentricity(elem), method)
    ν = _trueanom_from_eccanom(elem, EA)
    return orbitsolve_ν(elem, ν, EA, t_em, compensated)
end

#=

We want to calculate the emission time for an object at some distance for us to 
receive light at a given observation time. 

at t0, it is at distance 0
at t1, it is at distance z1
the time it takes for us to receive the signal from z1 is Δ
so at the emission time t1, we are going to receive it at observation time t1 + Δ

So when we're considering an observation time t1, to first order the emission time 
was t1 - Δ. There is a slight correction needed from the fact that at t1 - Δ, the
star wasn't as far away, meaning that Δ is an over correction. 

We need to solve implicitly for 
t_em = t_obs + Δt(distance_a_t_em)

How does our code fit into this?

t_obs + Δt(t_em) = compensate_star_3d_motion(elem, t_em).compensated.epoch_2a

This propagates the straight-line motion, and tells us how the light travel time has changed
either its gone positive (further away) or gone negative( closer to us) at the new epoch 
measured away from the reference epoch.
=#

# Calculate rigorous 3D motion propagation of star (though, not impacts of planets).
# Account for light travel time by calculating the retarded time iteratively.
function orbitsolve(elem::AbsoluteVisualOrbit{T}, t_obs::Number, method::AbstractSolver=Auto(); ) where T

    # t_obs is fixed--we know when we observed!
    # we adjust t_em to find the distance at which the barycentre was to emit light
    # that is then seen at t_obs

    # adjust t_em such that t_obs′ is approximately equal to t_obs within tolerance

    # First-order light travel time guess
    light_travel_time = elem.rv*(t_obs - elem.ref_epoch)*60*60*24/c_light_ms  # in seconds
    t_em = t_obs + light_travel_time * sec2day  # convert to days for MJD

    # Now iterate using our full 3D propagation code. 
    # On typical nearby stars, 1 extra iteration gets it correct
    # to better than an hour, 2 extra iterations get it correct
    # to within 10ms.
    local compensated
    for _ in 1:2
        compensated = compensate_star_3d_motion(elem, t_em)
        t_obs′ = compensated.epoch2a_days
        t_em += (t_obs - t_obs′)
    end
    out =  _calculate_orbit_solution(elem, t_em, compensated, method)
    return out

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
    return @inline OrbitSolutionAbsoluteVisual(elem, sol, t, compensated)
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
    cart2angle = rad2as*oftype(xcart, 1e3)/o.compensated.distance2_au
    xang = xcart*cart2angle # [mas]
    return xang 
end
function decoff(o::OrbitSolutionAbsoluteVisual)
    ycart = posy(o) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.compensated.distance2_au
    yang = ycart*cart2angle # [mas]
    return yang
end
function pmra(o::OrbitSolutionAbsoluteVisual)
    ẋcart = o.elem.parent.J*(o.elem.parent.cosi_cosΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) - o.elem.parent.sinΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.compensated.distance2_au
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang + (o.compensated.pmra2 - o.elem.pmra)
end
function pmdec(o::OrbitSolutionAbsoluteVisual)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.compensated.distance2_au
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
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.compensated.distance2_au
    quantity = ẋang = ẋcart*cart2angle # [mas/year]
    M_tot = totalmass(o.elem)
    return -M_planet/M_tot*quantity + (o.compensated.pmra2 - o.elem.pmra)
end
function pmdec(o::OrbitSolutionAbsoluteVisual, M_planet)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.compensated.distance2_au
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
    # cart2angle = rad2as*oftype(ẍcart, 1e3)/o.compensated.distance2_au
    # ẍang = ẍcart*cart2angle # [mas/year^2] 
    # return ẍang
end
function accdec(o::OrbitSolutionAbsoluteVisual)
    # throw(NotImplementedException())
    # if eccentricity(o.elem) >= 1
    #     @warn "acceleration not tested for ecc >= 1 yet. Results are likely wrong."
    # end
    # ÿcart = o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_sinΩ*o.sol.sinν_ω - o.elem.parent.cosΩ*o.sol.cosν_ω) # [AU/year^2]
    # cart2angle = rad2as*oftype(ÿcart, 1e3)/o.compensated.distance2_au
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
