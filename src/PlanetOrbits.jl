"""
# PlanetOrbits
A package for calculating orbits in the context of direct imaging,
astrometry, and radial velocity.
"""

module PlanetOrbits

# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

using LinearAlgebra
using StaticArrays

# ----------------------------------------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------------------------------
# KeplerianElements
# ----------------------------------------------------------------------------------------------------------------------

abstract type AbstractOrbit end
export AbstractOrbit


# ----------------------------------------------------------------------------------------------------------------------
# Orbit Solutions
# ----------------------------------------------------------------------------------------------------------------------
abstract type AbstractOrbitSolution end


# ----------------------------------------------------------------------------------------------------------------------
# System Properties
# ----------------------------------------------------------------------------------------------------------------------

"""
    period(elem)

Period of an orbit [days].
"""
function period end
export period

"""
    distance(elem)

Distance to the system [pc].
"""
function distance end
export distance

"""
    meanmotion(elem)

Mean motion [rad/year].
"""
function meanmotion end
export meanmotion

"""
   periastron(elements, tref=58849)

Compute the MJD of periastron passage most recently after the reference epoch tref.
N.B. mjd of 58849 = 2020-01-01
"""
function periastron end
export periastron

"""
    semiamplitude(elem)

Radial velocity semiamplitude [m/s].
"""
function semiamplitude end
export semiamplitude

# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------

# TODO: update docs
"""
    orbitsolve(elements, t)

Given a set of orbital elements with a time `t` in days, get the position and
velocity of the secondary body (e.g. planet around a star).

This will output a struct that is a subtype of `AbstractOrbitSolution` which
we can then query with `raoff`, `decoff`, `radvel`, etc.
    
You can also calculate those quanitities individually (see their docstrings) 
but if you need more than one, it is most efficient to save the orbit solution
once.

Note: these calculations use the small angle approximation, so are only accurate when 
the star is much further way from the observer than the secondary is from the primary.

See also: `orbitsolve_??`, `projectedseparation`, `raoff`, `decoff`, `radvel`, `propmotionanom`.
"""
function orbitsolve end




export orbitsolve

# ----------------------------------------------------------------------------------------------------------------------
# Orbital Position and Motion
# ----------------------------------------------------------------------------------------------------------------------

"""
    raoff(elem, t)

Get the offset [mas] from the primary body in Right Ascension
at the time `t` [days].

    raoff(o)

Get the offset [mas] from the primary body in Right Ascension 
from an instance of `AbstractOrbitSolution`.
"""
function raoff end
export raoff

"""
    decoff(elem, t)

Get the offset [mas] from the primary body in Declination
at the time `t` [days].

    decoff(elem, t)

Get the offset [mas] from the primary body in Declination
from an instance of `AbstractOrbitSolution`.
"""
function decoff end
export decoff

"""
    posangle(elem, t)

Calculate the position angle [rad] of the secondary about its primary
from our perspective at the time `t` [days].

    posangle(o)

Calculate the position angle [rad] of the secondary about its primary
from our perspective from an instance of `AbstractOrbitSolution`.
"""
function posangle(o::AbstractOrbitSolution)
    x = raoff(o)
    y = decoff(o)
    return atan(x, y) # Note: the order of these arguments is *correct* in our conventions
end
export posangle

"""
    projectedseparation(elem, t)

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
    propmotionanom(elem, t)

Get the instantaneous proper motion anomaly [mas/year] of
the *secondary* at the time `t` [days].

    propmotionanom(o)

Get the instantaneous proper motion anomaly [mas/year] of
the *secondary* from an instance of `AbstractOrbitSolution`.
"""
function propmotionanom(o::AbstractOrbitSolution)
    ????_planet = SVector(pmra(o), pmdec(o))
    return ????_planet
end

function pmra end

function pmdec end

export pmra, pmdec

"""
    propmotionanom(elem, t, M_planet)

Get the instantaneous proper motion anomaly [mas/year] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.
"""
function propmotionanom(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    ????_planet = propmotionanom(o)
    ????_star = -(M_planet/(M_star + M_planet))*????_planet
    return ????_star
end
export propmotionanom

"""
    radvel(elem, t)

Get the radial velocity [m/s] of the *secondary* along the
line of sight at the time `t` [days].

    radvel(o)

Get the radial velocity [m/s] of the *secondary* along the
line of sight from an instance of `AbstractOrbitSolution`.
"""
function radvel end

"""
    radvel(elem, t, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

radvel(elem, t, M_planet)

Get the radial velocity [m/s] of the *primary* along the
line of sight from an `AbstractOrbitSolution`. The units of `M_planet`
and `elem.M` must match.
"""
function radvel(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    v_planet = radvel(o)
    v_star = -(M_planet/(M_star + M_planet))*v_planet 
    return v_star
end
export radvel

"""
    acceleration(elem, t)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* at the time `t` [days].

    acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* from an instance of `AbstractOrbitSolution`.
"""
function acceleration(o::AbstractOrbitSolution)
    acc_planet = SVector(accra(o), accdec(o))
    return acc_planet
end
function accra end
function accdec end
export accra, accdec

"""
    acceleration(elem, t, M_planet)

Get the instantaneous acceleration [mas/year^2] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *primary* from an instance of `AbstractOrbitSolution`. The units of
`M_planet` and `elem.M` must match.
"""
function acceleration(o::AbstractOrbitSolution, M_planet)
    M_star = o.elem.M
    acc_planet = acceleration(o)
    acc_star = -(M_planet/(M_star + M_planet))*acc_planet
    return acc_star
end
export acceleration


include("keplerian.jl")
include("radvel.jl")


function orbitsolve(elem::Union{KeplerianElements,RadialVelocityElements}, t; tref=58849)
    
    # Epoch of periastron passage
    t??? = periastron(elem, tref)

    if t isa Integer
        t = float(t)
    end
    # Mean anomaly    
    MA = meanmotion(elem)/oftype(t, year2day) * (t - t???)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ?? = 2*atan(elem.??_fact*tan(EA/2))

    return orbitsolve_??(elem, ??; EA)
end

function orbitsolve_meananom(elem::Union{KeplerianElements,RadialVelocityElements}, MA)
    
    # Compute eccentric anomaly
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ?? = 2*atan(elem.??_fact*tan(EA/2))

    return orbitsolve_??(elem, ??; EA)
end

function orbitsolve_eccanom(elem::Union{KeplerianElements,RadialVelocityElements}, EA)
        
    # Calculate true anomaly
    ?? = 2*atan(elem.??_fact*tan(EA/2))

    return orbitsolve_??(elem, ??)
end

function radvel(o::Union{OrbitSolutionKeplerian,OrbitSolutionRadialVelocity})
    z??cart = o.elem.K*(o.cos??_?? + o.elem.ecos??) # [m/s]
    return z??cart
end

# Define fallbacks for all accessor functions.
# If the user calls f(elems, t, args...) we compute the
# AbstractOrbitSolution for them.
fun_list = (
    :raoff,
    :decoff,
    :posangle,
    :projectedseparation,
    :propmotionanom,
    :radvel,
    :pmra,
    :pmdec,
    :accra,
    :accdec,
    :acceleration,
)
for fun in fun_list
    @eval function ($fun)(elem::AbstractOrbit, t::Real, args...)
        return ($fun)(orbitsolve(elem, t), args...)
    end
end

# ----------------------------------------------------------------------------------------------------------------------
# Kepler Equation Solvers
# ----------------------------------------------------------------------------------------------------------------------
include("kepsolve_markley.jl")

# Fallback kepler solver function.
# If algorithm is unspecified, select the best one here.
function kepler_solver(MA, e)
    kepler_solver(MA, e, Markley())
end


# ----------------------------------------------------------------------------------------------------------------------
# Addional & Optional Features
# ----------------------------------------------------------------------------------------------------------------------
include("diff-rules.jl")
include("transformation.jl")
include("time.jl")
include("plots-recipes.jl")

using Requires
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("makie-recipes.jl")

    # Small patch to allow symbolic tracing through the kepler solver.
    # Mean anomaly must still be in the range [0,2??] for the solution
    # to be valid.
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" include("symbolics.jl")
end

end # module

# ----------------------------------------------------------------------------------------------------------------------