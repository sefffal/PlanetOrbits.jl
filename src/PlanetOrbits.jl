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
abstract type AbstractOrbit end
export AbstractOrbit

"""
    AbstractOrbitSolution

Represents the solution of an orbit. Contains all
the information of an AbstractOrbit, plus information
necessary to uniquely locate a planet.

The solution can be queried using a variety of functions
such as `radvel(solution)`.
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

# ---------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ---------------------------------------------------

"""
    orbitsolve(elements, t, method=Auto())

Given a set of orbital elements with a time `t` in days, get the position and
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
    posx(elem, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posx(elem, t)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posx end

"""
    posy(elem, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posy(o)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posy end

"""
    posz(elem, t)

Get the offset [AU] from the primary body at the time `t` [days].

    posz(o)

Same as above, but from an instance of `AbstractOrbitSolution`.
"""
function posz end

"""
    posangle(elem, t)

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

    propmotionanom(elem, t, M_planet)

Get the instantaneous proper motion anomaly [mas/year] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.


    propmotionanom(o, M_planet)

Same as above, but from an orbit solution.
"""
function propmotionanom(o::AbstractOrbitSolution)
    Δμ_planet = SVector(pmra(o), pmdec(o))
    return Δμ_planet
end

function pmra end

function pmdec end

export pmra, pmdec

export propmotionanom

"""
    radvel(elem, t)

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
    acceleration(elem, t)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* at the time `t` [days].

    acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *secondary* from an instance of `AbstractOrbitSolution`.

    acceleration(elem, t, M_planet)

Get the instantaneous acceleration [mas/year^2] of 
the *primary* in at the time `t` [days]. The units of `M_planet`
and `elem.M` must match.

    acceleration(o)

Get the instantaneous acceleration [mas/year^2] of
the *primary* from an instance of `AbstractOrbitSolution`. The units of
`M_planet` and `elem.M` must match.
"""
function acceleration(o::AbstractOrbitSolution)
    acc_planet = SVector(accra(o), accdec(o))
    return acc_planet
end
function accra end
function accdec end
export accra, accdec
export acceleration


"""
    trueanom(orbit, t)

Get the true anomaly [radians] of the *secondary*
at the time `t` [days].

    trueanom(o)

Get the true anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.
"""
trueanom(os::AbstractOrbitSolution) = os.ν
"""
    eccanom(orbit, t)

Get the eccentric anomaly [radians] of the *secondary*
at the time `t` [days].

    eccanom(o)

Get the eccentric anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.
"""
eccanom(os::AbstractOrbitSolution) = os.EA
"""
    meananom(orbit, t)

Get the mean anomaly [radians] of the *secondary*
at the time `t` [days].

    meananom(o)

Get the mean anomaly [radians] of the *secondary*
from an instance of `AbstractOrbitSolution`.
"""
meananom(os::AbstractOrbitSolution) = eccanom(os) - os.elem.e * sin(eccanom(os))
export trueanom, eccanom, meananom


# Define iterate and length = 1 so that we can broadcast over elements.
Base.length(::AbstractOrbit) = 1
Base.iterate(elem::AbstractOrbit) = (elem, nothing)
Base.iterate(::AbstractOrbit, ::Nothing) = nothing


include("kep-elements.jl")
include("visual-elements.jl")

function posx(o::Union{OrbitSolutionKep, OrbitSolutionVisual})
    xcart = o.r*(o.cosν_ω*o.elem.sinΩ + o.sinν_ω*o.elem.cosi*o.elem.cosΩ) # [AU]
    return xcart
end
function posy(o::Union{OrbitSolutionKep, OrbitSolutionVisual})
    ycart = o.r*(o.cosν_ω*o.elem.cosΩ - o.sinν_ω*o.elem.cosi*o.elem.sinΩ) # [AU]
    return ycart
end
function posz(o::Union{OrbitSolutionKep, OrbitSolutionVisual})
    zcart = o.r*(o.cosν_ω*o.elem.sini) # [AU]
    return zcart
end


include("radvel-elements.jl")


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
- τ: epoch of periastron passage at MJD=0, default=0
- e: eccentricity, default=0
- ω: argument of periapsis [rad], default=0
- i: inclination [rad]
- Ω: longitude of ascending node [rad]
- plx: parallax [mas]; defines the distance to the primary
"""
function orbit(;a,M,e=0.0, ω=0.0, τ=0.0, kwargs...)
    T = supportedorbit(kwargs)
    if T == VisualOrbit
        @info "Selected VisualOrbit(;a, e, i, ω, Ω, τ, M, plx)"
        T(;a, e, ω, τ, M, kwargs...)
    elseif T == KepOrbit
        @info "Selected KepOrbit(;a, e, i, ω, Ω, τ, M)"
        T(;a, e, ω, τ, M, kwargs...)
    else
        @info "Selected RadialVelocityOrbit(;a, e, ω, τ, M)"
        T(;a, e, ω, τ, M)
    end
end
# Function to return what orbit type is supported based on precence
# or absence of properties
function supportedorbit(kwargs)
    if haskey(kwargs, :plx)
        VisualOrbit
    elseif haskey(kwargs, :i)
        KepOrbit
    else
        RadialVelocityOrbit
    end
end
export orbit



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

include("kepsolve_goat.jl")
include("kepsolve_markley.jl")

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
# include("kepsolve_roots.jl") # only pull this one in if the user has Roots loaded.

# Fallback kepler solver function.
# If algorithm is unspecified, select the best one here.
kepler_solver(MA, e) = kepler_solver(MA, e, Auto())
function kepler_solver(MA, e, ::Auto)
    kepler_solver(MA, e, Markley())
end



function orbitsolve(elem::AbstractOrbit, t, method::AbstractSolver=Auto(); tref=58849)
    
    # Epoch of periastron passage
    tₚ = periastron(elem, tref)

    if t isa Integer
        t = float(t)
    end
    # Mean anomaly    
    MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, elem.e, method)
    
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*ourtan(EA/2))

    return orbitsolve_ν(elem, ν; EA)
end


"""
    orbitsolve_ν(elem, ν; EA)

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
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν; EA)
end

"""
    orbitsolve_eccanom(elements, EA)

Same as `orbitsolve`, but solves orbit for a given eccentric anomaly instead of time.
"""
function orbitsolve_eccanom(elem::AbstractOrbit, EA)
        
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν)
end

function radvel(o::AbstractOrbitSolution)
    żcart = o.elem.K*(o.cosν_ω + o.elem.ecosω) # [m/s]
    return żcart
end

# Define fallbacks for all accessor functions.
# If the user calls f(elems, t, args...) we compute the
# AbstractOrbitSolution for them.
fun_list = (
    :trueanom,
    :eccanom,
    :posx,
    :posy,
    :posz,
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
        M_star = o.elem.M
        return -(M_planet/(M_star + M_planet))*quantity
    end
end
function projectedseparation(o::AbstractOrbitSolution, M_planet)
    quantity = projectedseparation(o)
    M_star = o.elem.M
    return (M_planet/(M_star + M_planet))*quantity
end
function posangle(o::AbstractOrbitSolution, M_planet)
    x = raoff(o,M_planet)
    y = decoff(o,M_planet)
    return atan(x, y) # Note: the order of these arguments is *correct* in our conventions
end



# ---------------------------------------------------
# Addional & Optional Features
# ---------------------------------------------------
include("diff-rules.jl")
include("transformation.jl")
include("time.jl")
include("plots-recipes.jl")

using Requires
function __init__()

    @require Roots="f2b01f46-fcfa-551c-844a-d8ac1e96c665" include("kepsolve_roots.jl")

    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("makie-recipes.jl")

    # Small patch to allow symbolic tracing through the kepler solver.
    # Mean anomaly must still be in the range [0,2π] for the solution
    # to be valid.
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" include("symbolics.jl")
end

end # module

# ---------------------------------------------------