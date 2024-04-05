"""
# PlanetOrbits
A package for calculating orbits in the context of direct imaging,
astrometry, and radial velocity.
"""

module PlanetOrbits

# ---------------------------------------------------
# Imports
# ---------------------------------------------------

using LinearAlgebra
using StaticArrays

# ---------------------------------------------------
# Constants
# ---------------------------------------------------

# radians <-> milliarcseconds
const rad2mas = 2.06264806e8
const mas2rad = 4.848136816903219e-9

# radians <-> arcseconds 
const rad2as = 206265
const as2rad = 4.848132257047972e-6

# parsecs <-> astronomical units
const pc2au = 206265
const au2pc = 4.848132257047972e-6

# astronomical units <-> metres
const au2m = 1.495978707e11
const m2au = 6.684587122268445e-12

# years <-> days
const year2day = 365.2422
const day2year = 2.737909255830788e-3

# years <-> seconds
const year2sec = 3.1556926e7
const sec2year = 3.168876461541279e-8

# days <-> seconds
const day2sec = 86400
const sec2day = 1.1574074074074073e-5

# jupiter masses <-> solar masses
const mjup2msol = 0.0009543

# ---------------------------------------------------
# Type Hierarchy
# ---------------------------------------------------

"""
    AbstractOrbit

Represents a orbit. Contains all the information to
calculate the location of a planet at a given time,
true anomaly, eccentric anomaly, or mean anomaly.
Different concrete implementations of AbstractOrbit
contain varying amounts of information.

Basic information about the orbit can be queried using
functions like `period(orbit)`.

Orbits can be solved using functions like `orbitsolve(orb)`.

See: `RadialVelocityOrbit`, `KepOrbit`, `VisualOrbit`
"""
abstract type AbstractOrbit{T} end
export AbstractOrbit

"""
    AbstractOrbitSolution

Represents the solution of an orbit. Contains all
the information of an AbstractOrbit, plus information
necessary to uniquely locate a planet.

The solution can be queried using a variety of functions
such as `radvel(solution)`.

The API for creating orbit solutions it not considered public
as the fields may change between minor versions. Instead,
create solutions only through the public `orbitsolve` and
`orbitsolve_...` functions.
"""
abstract type AbstractOrbitSolution end
export AbstractOrbitSolution

# Return the orbit solution type for a given orbit type.
# Register for each orbit type.
function _solution_type end
_solution_type(o::Any) = _solution_type(typeof(o))

# ---------------------------------------------------
# System Properties
# ---------------------------------------------------

"""
    period(orbit)

Period of an orbit [days].
"""
function period end
export period



"""
    totalmass(orbit)

Total mass of the system in solar masses
"""
function totalmass end
export totalmass

"""
    distance(orbit)

Distance to the system [pc].
"""
function distance end
export distance

"""
    meanmotion(orbit)

Mean motion [rad/year].
"""
function meanmotion end
export meanmotion


"""
    eccentricity(orbit)

Eccentricity of an orbit, between 0 and 1.
"""
function eccentricity end
export eccentricity


"""
    inclination(orbit)

Inclination of an orbit, if available [rad].
"""
function inclination end
export inclination

"""
    semimajoraxis(orbit)

Semi-major axis of an orbit, if available [au].
"""
function semimajoraxis end
export semimajoraxis


"""
   periastron(elements)

Compute the MJD of periastron passage most recently after the reference epoch tref specified in the orbit.
N.B. mjd of 58849 = 2020-01-01
"""
function periastron end
export periastron

"""
    semiamplitude(orbit)

Radial velocity semiamplitude [m/s].
"""
function semiamplitude end
export semiamplitude


# ---------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ---------------------------------------------------

"""
    orbitsolve(orbit, t, method=Auto())

Given an orbit object and a time `t` in days, get the position and
velocity of the secondary body (e.g. planet around a star).

This will output a struct that is a subtype of `AbstractOrbitSolution` which
we can then query with `raoff`, `decoff`, `radvel`, etc.
    
You can also calculate those quanitities individually (see their docstrings) 
but if you need more than one, it is most efficient to save the orbit solution
once.

Note: these calculations use the small angle approximation, so are only accurate when 
the star is much further way from the observer than the secondary is from the primary.

See also: `orbitsolve_ν`,  `orbitsolve_meananom`,  `orbitsolve_eccanom`, `projectedseparation`, `raoff`, `decoff`, `radvel`, `propmotionanom`.
"""
function orbitsolve end




export orbitsolve, orbitsolve_ν, orbitsolve_meananom, orbitsolve_eccanom

# ---------------------------------------------------
# Orbital Position and Motion
# ---------------------------------------------------

"""
    raoff(orbit, t)

Get the offset [mas] from the primary body in Right Ascension
at the time `t` [days].

    raoff(o)

Get the offset [mas] from the primary body in Right Ascension 
from an instance of `AbstractOrbitSolution`.
"""
function raoff end
export raoff

"""
    decoff(orbit, t)

Get the offset [mas] from the primary body in Declination
at the time `t` [days].

    decoff(orbit, t)

Get the offset [mas] from the primary body in Declination
from an instance of `AbstractOrbitSolution`.
"""
function decoff end
export decoff


"""
    posx(orbit, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posx(orbit, t)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posx end

"""
    posy(orbit, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posy(o)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posy end

"""
    posz(orbit, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posz(o)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posz end

"""
    posangle(orbit, t)

Calculate the position angle [rad] of the secondary about its primary
from our perspective at the time `t` [days].

    posangle(o)

Calculate the position angle [rad] of the secondary about its primary
from our perspective from an instance of `AbstractOrbitSolution`.

    posangle(elem, t, M_planet)

Calculate the position angle [rad] of the secondary about its primary
from our perspective at the time `t` [days].
In this case only, the value of M_planet can be arbitrary.

    posangle(o, M_planet)

Calculate the position angle [rad] of the **primary** 
from our perspective from an instance of `AbstractOrbitSolution`.
In this case only, the value of M_planet can be arbitrary.
"""
function posangle(o::AbstractOrbitSolution)
    x = posx(o)
    y = posy(o)
    return atan(x, y) # Note: the order of these arguments is *correct* in our conventions
end
export posangle

"""
    projectedseparation(orbit, t)

Calculate the projected separation [mas] of the secondary from its
primary at the time `t` [days].

    projectedseparation(o)

Calculate the projected separation [mas] of the secondary from its
primary from an instance of `AbstractOrbitSolution`.
"""
function projectedseparation(o::AbstractOrbitSolution)
    x = raoff(o)
    y = decoff(o)
    return sqrt(x^2 + y^2)
end
export projectedseparation

"""
    pmra(orbit, t)

Get the instantaneous proper motion anomaly [mas/year] in right-ascension of
the *secondary* at the time `t` [days].

    pmra(o)

Get the instantaneous proper motion anomaly [mas/year] in right-ascension of
the *secondary* from an instance of `AbstractOrbitSolution`.

    pmra(elem, t, M_planet)

Get the instantaneous proper motion anomaly [mas/year] in right-ascension of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.


    pmra(o, M_planet)

Same as above, but from an orbit solution.
"""
function pmra end


"""
    pmdec(orbit, t)

Get the instantaneous proper motion anomaly [mas/year] in declination of
the *secondary* at the time `t` [days].

    pmdec(o)

Get the instantaneous proper motion anomaly [mas/year] in declination of
the *secondary* from an instance of `AbstractOrbitSolution`.

    pmdec(elem, t, M_planet)

Get the instantaneous proper motion anomaly [mas/year] in declination of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.


    pmdec(o, M_planet)

Same as above, but from an orbit solution.
"""
function pmdec end

export pmra, pmdec

"""
    radvel(orbit, t)

Get the radial velocity [m/s] of the *secondary* along the
line of sight at the time `t` [days].

    radvel(o)

Get the radial velocity [m/s] of the *secondary* along the
line of sight from an instance of `AbstractOrbitSolution`.

    radvel(elem, t, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

    radvel(o, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight from an `AbstractOrbitSolution`. The units of `M_planet`
and `elem.M` must match.
"""
function radvel end


export radvel

"""
    accra(orbit, t)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *secondary* at the time `t` [days].

    accra(o)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *secondary* from an instance of `AbstractOrbitSolution`.

    accra(elem, t, M_planet)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

    accra(o)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *primary* from an instance of `AbstractOrbitSolution`. The units of
`M_planet` and `elem.M` must match.
"""
function accra end


"""
    accdec(orbit, t)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *secondary* at the time `t` [days].

    accdec(o)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *secondary* from an instance of `AbstractOrbitSolution`.

    accdec(elem, t, M_planet)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

    accdec(o)

Get the instantaneous acceleration [mas/year^2] in the right-ascension direction of
the *primary* from an instance of `AbstractOrbitSolution`. The units of
`M_planet` and `elem.M` must match.
"""
function accdec end
export accra, accdec


"""
    trueanom(orbit, t)

Get the true anomaly [radians] of the *secondary*
at the time `t` [days].

    trueanom(o)

Get the true anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.
"""
trueanom(os::AbstractOrbitSolution) = os.ν
trueanom(os::AbstractOrbitSolution, mass::Number) = trueanom(os) # Same for primary and secondary
"""
    eccanom(orbit, t)

Get the eccentric anomaly [radians] of the *secondary*
at the time `t` [days].

    eccanom(o)

Get the eccentric anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.

Note that for hyperbolic orbits, eccentric anomaly is not defined and the hyperbolic anomaly is returned instead.
"""
eccanom(os::AbstractOrbitSolution) = os.EA
eccanom(os::AbstractOrbitSolution, mass::Number) = eccanom(os) # Same for primary and secondary
"""
    meananom(orbit, t)

Get the mean anomaly [radians] of the *secondary*
at the time `t` [days].

    meananom(o)

Get the mean anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.
"""
function meananom(os::AbstractOrbitSolution)
    if os.elem.e < 1
        return eccanom(os) - os.elem.e * sin(eccanom(os))
    else
        return os.elem.e * sinh(eccanom(os)) - eccanom(os)
    end
end
meananom(os::AbstractOrbitSolution, mass::Number) = meananom(os) # Same for primary and secondary
export trueanom, eccanom, meananom

"""
    periapsis(orbit)

Return the periapsis of an orbit in AU.

Keywords: periastron, perihelion, perigee
"""
function periapsis(o::AbstractOrbit)
    if eccentricity(o) < 1
        semimajoraxis(o)*(1 - eccentricity(o))
    else
        -semimajoraxis(o)*(eccentricity(o) - 1)
    end
end

"""
    apoapsis(orbit)

Return the apoapsis of an orbit in AU.

Keywords: apoastron, apohelion, apogee
"""
function apoapsis(o::AbstractOrbit)
    if eccentricity(o) < 1
         semimajoraxis(o)*(1 + eccentricity(o))
    else
        -semimajoraxis(o)*(1 + eccentricity(o))
    end
end

"""
    semiminoraxis(orbit)

Return the semi-minor axis of an orbit in AU.
"""
function semiminoraxis(o::AbstractOrbit)
    if eccentricity(o) < 1
        semimajoraxis(o)*sqrt(1-eccentricity(o)^2)
    else
        semimajoraxis(o)*sqrt(eccentricity(o)^2 - 1)
    end
end

export periapsis, apoapsis, semiminoraxis

# Internal function used by each orbit type to map mean anomaly to true anomaly
function _trueanom_from_eccanom end


# Define iterate and length = 1 so that we can broadcast over elements.
Base.length(::AbstractOrbit) = 1
Base.iterate(elem::AbstractOrbit) = (elem, nothing)
Base.iterate(::AbstractOrbit, ::Nothing) = nothing




# ---------------------------------------------------
# Kepler Equation Solvers
# ---------------------------------------------------
abstract type AbstractSolver end

"""
    PlanetOrbits.Auto()

Automatic choice of Kepler solver algorithm.
Currently defaults to PlanetOrbits.Markley()
"""
struct Auto <: AbstractSolver end

include("kepsolve-goat.jl")
include("kepsolve-markley.jl")

"""
    PlanetOrbits.RootsMethod(method::Roots.PlanetOrbits.Roots.AbstractUnivariateZeroMethod, kwargs...)

Wraps a root finding method from Roots.jl. Requires Roots to be loaded first.
You can also pass keyword arguments that will be forwarded to Roots to control
the tolerance.

Examples:
```julia
method = PlanetOrbits.RootsMethod(Roots.Newton())
method = PlanetOrbits.RootsMethod(Roots.Thukral5B())
method = PlanetOrbits.RootsMethod(Roots.Bisection())
method = PlanetOrbits.RootsMethod(Roots.A42())
method = PlanetOrbits.RootsMethod(Roots.Newton(), rtol=1e-3, atol=1e-3)
```
"""
struct RootsMethod{M,K} <: AbstractSolver
    method::M
    kwargs::K
end
RootsMethod(method; kwargs...) = RootsMethod(method, kwargs)
include("kepsolve-roots.jl")

# Fallback kepler solver function.
# If algorithm is unspecified, select the best one here.
kepler_solver(MA, e) = kepler_solver(MA, e, Auto())
function kepler_solver(MA, e, ::Auto)
    if e < 1
        kepler_solver(MA, e, Markley())
    else
        # Halley() converged slightly faster than Newton() for hyperbolic orbits
        kepler_solver(MA, e, RootsMethod(Roots.Halley()))
        # kepler_solver(MA, e, RootsMethod(Roots.Bisection()))
    end
end



function orbitsolve(elem::AbstractOrbit, t, method::AbstractSolver=Auto())
    
    # Epoch of periastron passage
    tₚ = periastron(elem)

    if t isa Integer
        t = float(t)
    end
    # Mean anomaly
    MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem), method)
    
    # Calculate true anomaly
    ν = _trueanom_from_eccanom(elem, EA)

    return orbitsolve_ν(elem, ν, EA, t) # optimization: Don't have to recalculate EA and t.
end


"""
    orbitsolve_ν(elem, ν, EA)

Solve an orbit from a given true anomaly [rad].
See `orbitsolve` for the same function accepting a given time.
Can optionally pass eccentric anomaly (EA) if already computed.
"""
function orbitsolve_ν end

"""
    orbitsolve_meananom(elements, MA)

Same as `orbitsolve`, but solves orbit for a given mean anomaly instead of time.
"""
function orbitsolve_meananom(elem::AbstractOrbit, MA)
    
    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem))
    
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν, EA)
end

"""
    orbitsolve_eccanom(elements, EA)

Same as `orbitsolve`, but solves orbit for a given eccentric anomaly instead of time.
"""
function orbitsolve_eccanom(elem::AbstractOrbit, EA)
        
    # Calculate true anomaly
    ν = _trueanom_from_eccanom(elem, EA)
    return orbitsolve_ν(elem, ν)
end

function radvel(o::AbstractOrbitSolution)
    żcart = o.elem.K*(o.cosν_ω + o.elem.ecosω) # [m/s]
    return żcart
end

function _time_from_EA(sol::AbstractOrbitSolution, EA;)
    elem = sol.elem

    if eccentricity(elem) < 1
        # Epoch of periastron passage
        tₚ = periastron(elem)

        MA = EA - eccentricity(elem) * sin(EA) 

        # Mean anomaly    
        t = MA/meanmotion(elem)*oftype(EA, year2day) + tₚ
        
    else
        # Epoch of periastron passage
        tₚ = periastron(elem)
        MA = -EA + eccentricity(elem)*sinh(EA)
        t = MA/meanmotion(elem)*oftype(EA, year2day) + tₚ

    end

    return t
end

# Given an eccentric anomaly, calculate *a* time at which the body 
# would be at that location.
function _time_from_EA(elem::AbstractOrbit, EA;)

    if eccentricity(elem) < 1
        # Epoch of periastron passage
        tₚ = periastron(elem)

        MA = EA - eccentricity(elem) * sin(EA) 

        # Mean anomaly    
        t = MA/meanmotion(elem)*oftype(EA, year2day) + tₚ
        
    else
        # Epoch of periastron passage
        tₚ = periastron(elem)
        MA = -EA + eccentricity(elem)*sinh(EA)
        t = MA/meanmotion(elem)*oftype(EA, year2day) + tₚ

    end

    return t



    # # ---- Worked math for elliptical case --- 
    # # ν/2 = atan(elem.ν_fact*tan(EA/2))
    # # tan(ν/2) = elem.ν_fact*tan(EA/2)
    # # tan(ν/2)/elem.ν_fact = tan(EA/2)
    # # atan(tan(ν/2)/elem.ν_fact) = (EA/2)
    # # atan(tan(ν/2)/elem.ν_fact)*2 = EA
    # # EA = atan(tan(ν/2)/elem.ν_fact)*2

    # # Compute eccentric anomaly
    # MA = EA - elem.e * sin(EA)

    # # Epoch of periastron passage
    # tₚ = periastron(elem)
    
    
    # # MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)
    # # MA /  meanmotion(elem) * year2day = (t - tₚ)
    # # MA /  meanmotion(elem) * year2day + tₚ = t
    # t = MA /  meanmotion(elem) * year2day + tₚ - tref

end


include("orbit-keplerian.jl")
include("orbit-visual.jl")
include("orbit-compensated.jl")
include("orbit-thiele-innes.jl")
include("orbit-radvel.jl")
include("orbit-cartesian.jl")

function orbitsolve_meananom(elem::VisualOrbit, MA)
    
    # Compute eccentric anomaly
    EA = kepler_solver(MA, eccentricity(elem))
    
    # Calculate true anomaly
    ν = 2*atan(elem.parent.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν, EA)
end


"""
Get the position in the x direction in astronomical units.
"""
function posx(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    xcart = o.r*(o.cosν_ω*o.elem.sinΩ + o.sinν_ω*o.elem.cosi*o.elem.cosΩ) # [AU]
    return xcart
end
"""
Get the position in the y direction in astronomical units.
"""
function posy(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    ycart = o.r*(o.cosν_ω*o.elem.cosΩ - o.sinν_ω*o.elem.cosi*o.elem.sinΩ) # [AU]
    return ycart
end
"""
Get the position in the z direction in astronomical units.
"""
function posz(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    zcart = o.r*(o.sinν_ω*o.elem.sini) # [AU]
    return zcart
end
export posx, posy, posz

"""
Get the velocity in the x direction in astronomical units / year.
"""
function velx(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    ẋcart = o.elem.J*(o.elem.cosi_cosΩ*(o.cosν_ω + o.elem.ecosω) - o.elem.sinΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    return ẋcart
end
"""
Get the velocity in the y direction in astronomical units / year.
"""
function vely(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    ẏcart = -o.elem.J*(o.elem.cosi_sinΩ*(o.cosν_ω + o.elem.ecosω) + o.elem.cosΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    return ẏcart
end
"""
Get the velocity in the z direction in astronomical units / year.
"""
function velz(o::Union{OrbitSolutionKep, OrbitSolutionCartesian})
    żcart = radvel(o) * m2au * year2sec
    return żcart
end
export velx, vely, velz


"""
    orbit(...)

Construct an orbit from the provided keyword arguments. Will automatically select
a subclass of AbstractOrbit based on the information provided. This is a convenience
function that is not type stable and should not be used in performance sensitive
contexts. Instead, call one of the concrete constructors `KepOrbit`, `VisualOrbit`,
or `RadialVelocityOrbit` directly.
This function logs the kind of elements created so that it's easy to select the correct
constructor.

Required arguments:
- a: semi-major axis [AU]
- M: mass of primary [M⊙]

Optional arguments:
- tp: epoch of periastron passage, default=0
- e: eccentricity, default=0
- ω: argument of periapsis [rad], default=0
- i: inclination [rad]
- Ω: longitude of ascending node [rad]
- plx: parallax [mas]; defines the distance to the primary
"""
function orbit(;kwargs...)
    T = supportedorbit(kwargs)
    if !haskey(kwargs, :e)
        kwargs = (;kwargs...,e=0,ω=0)
    end
    if !haskey(kwargs, :tp)
        kwargs = (;kwargs...,tp=0)
    end
    return T(;kwargs...)
end
# Function to return what orbit type is supported based on precence
# or absence of properties
function supportedorbit(kwargs)
    OrbitType = 
        if haskey(kwargs, :x) && haskey(kwargs, :vx)
            CartesianOrbit
        elseif haskey(kwargs, :A)
            ThieleInnesOrbit
        elseif haskey(kwargs, :i)
            KepOrbit
        else
            RadialVelocityOrbit
        end
    if haskey(kwargs, :rv)
        return Compensated{OrbitType}
    elseif haskey(kwargs, :plx) && !(OrbitType==ThieleInnesOrbit)
        return Visual{OrbitType}
    else
        return OrbitType
    end
end
export orbit


# Define fallbacks for all accessor functions.
# If the user calls f(elems, t, args...) we compute the
# AbstractOrbitSolution for them.
fun_list = (
    :trueanom,
    :eccanom,
    :meananom,
    :posx,
    :posy,
    :posz,
    :raoff,
    :decoff,
    :posangle,
    :projectedseparation,
    :propmotionanom,
    :velx,
    :vely,
    :velz,
    :radvel,
    :pmra,
    :pmdec,
    :accra,
    :accdec,
    :acceleration,
)
for fun in fun_list
    @eval function ($fun)(orbit::AbstractOrbit, t::Real, args...)
        return ($fun)(orbitsolve(orbit, t), args...)
    end
end


# Define versions that compute the quantity on of the primary instead of 
# the secondary
mass_fun_list = (
    :posx,
    :posy,
    :posz,
    :radvel,
    :raoff,
    :decoff,
    :pmra,
    :pmdec,
    :accra,
    :accdec,
    :propmotionanom,
    :acceleration,
)
for fun in mass_fun_list
    @eval function ($fun)(o::AbstractOrbitSolution, M_planet)
        quantity = ($fun)(o)
        M_tot = totalmass(o.elem)
        return -M_planet/M_tot*quantity
    end
end
function projectedseparation(o::AbstractOrbitSolution, M_planet)
    quantity = projectedseparation(o)
    M_tot = totalmass(o.elem)
    return M_planet/M_tot*quantity
end
function posangle(o::AbstractOrbitSolution, M_planet)
    x = posx(o,M_planet)
    y = posy(o,M_planet)
    return atan(x, y) # Note: the order of these arguments is *correct* in our conventions
end



# ---------------------------------------------------
# Addional & Optional Features
# ---------------------------------------------------
include("recipes-plots.jl")
include("time.jl")
include("chain-rules.jl")
# include("transformation.jl")


include("precompile.jl")

end # module

# ---------------------------------------------------
