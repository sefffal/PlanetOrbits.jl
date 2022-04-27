"""
# DirectOrbits
A package for calculating orbits in the context of direct imaging,
astrometry, and radial velocity.
"""

module DirectOrbits

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

This will output an `AbstractOrbitSolution` struct with the following properties:
 - `x`: δ right ascension [mas]
 - `y`: δ declination [mas]
 - `ẋ`: right ascension proper motion anomaly [mas/year]
 - `ẏ`: declination proper motion anomaly [mas/year]
 - `ż`: radial velocity of the *secondary* [m/s]
 - `ẍ`: right ascension acceleration [mas/year^2]
 - `ÿ`: declination acceleration [mas/year^2]

You can access the properties by name `.x`. There are helper functions to
calculate each of these properties individually, but if you need more than
one it is most efficient to calculate them in one go.

`radvel` can optionally accept the mass of the primary to calculate the impact
of the secondary body on radial velocity of the primary, instead of the radial
velocity of the secondary body itself.

Note: these calculations use the small angle approximation, so are only accurate when 
the star is much further way from the observer than the secondary is from the primary.

See also: `orbitsolve_ν`, `projectedseparation`, `raoff`, `decoff`, `radvel`, `propmotionanom`.
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
    Δμ_planet = SVector(pmra(o), pmdec(o))
    return Δμ_planet
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
    Δμ_planet = propmotionanom(o)
    Δμ_star = -(M_planet/(M_star + M_planet))*Δμ_planet
    return Δμ_star
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
    tₚ = periastron(elem, tref)

    if t isa Integer
        t = float(t)
    end
    # Mean anomaly    
    MA = meanmotion(elem)/oftype(t, year2day) * (t - tₚ)

    # Compute eccentric anomaly
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ν = 2*atan(elem.ν_fact*tan(EA/2))

    return orbitsolve_ν(elem, ν)
end

function radvel(o::Union{OrbitSolutionKeplerian,OrbitSolutionRadialVelocity})
    żcart = o.elem.K*(o.cosν_ω + o.elem.ecosω) # [m/s]
    return żcart
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
# Kepler Equation Solver
# ----------------------------------------------------------------------------------------------------------------------

# The following function is taken directly from AstroLib.jl
# We  remove one invariant check we handle elsewhere and also
# force inlining for about a 5% speedup.
# We also supply analytic gradients for use in autodiff packages.
@inline function kepler_solver(_M::Real, e::Real)
    # We already handle this invariant
    # @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    M = rem2pi(_M, RoundNearest)
    T = float(promote_type(typeof(M), typeof(e)))
    if iszero(M) || iszero(e)
        return T(M)
    end
    pi2 = abs2(T(pi))
    # equation (20)
    α = (3 * pi2 + 8 * (pi2 - pi * abs(M)) / (5 * (1 + e)))/(pi2 - 6)
    # equation (5)
    d = 3 * (1 - e) + α * e
    # equation (9)
    q = 2 * α * d * (1 - e) - M * M
    # equation (10)
    r = 3 * α * d * (d - 1 + e) * M + M * M * M
    # equation (14)
    w = cbrt(abs2(abs(r) + sqrt(q * q * q + r * r)))
    # equation (15)
    E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M)/d
    # equation (26) & equation (27)
    f2, f3 = e .* sincos(E1)
    # equation (21)
    f0 = E1 - f2 - M
    # equation (25)
    f1 = 1 - f3
    # equation (22)
    δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
    # equation (23)
    δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
    # equations (24) and (28)
    δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
    return E1 + δ5 # equation 29
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
    # Mean anomaly must still be in the range [0,2π] for the solution
    # to be valid.
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" include("symbolics.jl")
end

end # module

# ----------------------------------------------------------------------------------------------------------------------