"""
# DirectOrbits
A package for calculating orbits in the context of direct imaging.

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

const G = 6.6743e-11
const mas2rad = 4.8481368e-9
const rad2as = 206265
const pc2au = 206265
const au2m = 1.495978707e11
const year2days = 365.2422
const day2secs = 86400
const sun2kg = 1.989e30

# ----------------------------------------------------------------------------------------------------------------------
# KeplerianElements
# ----------------------------------------------------------------------------------------------------------------------

abstract type AbstractElements end

"""
    KeplerianElements(
        a=1.0, # semi-major axis [AU]
        i=π/2, # inclination [radians]
        e=0.1, # eccentricity
        τ=π/2, # epoch of periastron passage at MJD=0
        M=1.0, # mass of primary [M⊙]
        ω=π/2, # argument of periapsis [radians]
        Ω=π/2, # longitude of the ascending node [radians]
        plx=10.1, # parallax [milliarcseconds]; defines the distance to the object
    )

Represents one object's Keplerian elements. Values can be specified
by keyword argument or named tuple for convenience.

See also `KeplerianElementsDeg` for a convenience constructor accepting
units of degrees instead of radians.
"""
struct KeplerianElements{T<:Number} <: AbstractElements

    # Orbital properties
    a::T
    i::T
    e::T
    τ::T
    M::T
    ω::T
    Ω::T
    plx::T

    # Cached constants for these elements.
    dist::T
    T::T
    n::T
    ν_fact::T
    cos_Ω::T
    sin_Ω::T
    cos_i::T
    sin_i::T
    cos_ω::T
    K::T

    # Inner constructor to inforce invariants and pre-calculate a few
    # constants for these elements.
    function KeplerianElements(a, i, e, τ, M, ω, Ω, plx)

        # Ensure validity of parameters
        if a <= 0.0
            @warn "Invalid semi-major axis" a maxlog=50
        end
        if !(0 <= e < 1)
            @warn "Eccentricity out of range" e maxlog=50
        end
        if M < 0.0
            @warn "Invalid primary mass (<0.001 Msun)" M maxlog=50
        end

        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        M = max(M, zero(M))
        plx = max(plx, zero(plx))
        τ = mod(τ, one(τ))

        # Pre-calculate some factors that will be re-used when calculating orbitsolve at any time
        dist = 1/(plx/1000) * pc2au # distance [AU]
        period = √(a^3/M) * year2days # period [days]
        n = 2π/√(a^3/M) # mean motion
        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(i), typeof(e), typeof(τ),
            typeof(M), typeof(ω), typeof(Ω), typeof(plx),
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        # Pre-calculate geometric factors
        sin_Ω, cos_Ω = sincos(Ω)
        sin_i, cos_i = sincos(i)
        cos_ω = cos(ω)

        # Radial velocity semiamplitude [m/s]
        K = cbrt((2π*G)/(period*day2secs)) * cbrt(M*sun2kg)*sin(i) * (1 - e^2)^(-1//2) 

        new{T}(
            # Passed parameters that define the elements
            a, i, e, τ, M, ω, Ω, plx,
            # Cached calcuations
            dist, period, n, ν_fact,
            # Geometric factors
            cos_Ω, sin_Ω, cos_i, sin_i, cos_ω, K
        )
    end
end

# Allow arguments to be specified by keyword
KeplerianElements(;a, i, e, τ, M, ω, Ω, plx) = KeplerianElements(a, i, e, τ, M, ω, Ω, plx)
# And by a named tuple without splatting
KeplerianElements(nt) = KeplerianElements(nt.a, nt.i, nt.e, nt.τ, nt.M, nt.ω, nt.Ω, nt.plx)

export KeplerianElements

"""
    astuple(elements)

Return the parameters of a KeplerianElements value as a tuple.
"""
function astuple(elem::KeplerianElements)
    return (;elem.a,elem.i,elem.e,elem.τ,elem.M,elem.ω,elem.Ω,elem.plx)
end
export astuple

"""
    KeplerianElementsDeg(a, i, e, τ, M, ω, Ω, plx)

A convenience function for constructing KeplerianElements where
`i`, `ω`, and `Ω` are provided in units of degrees instead of radians.
"""
KeplerianElementsDeg(a, i, e, τ, M, ω, Ω, plx) = KeplerianElements(a, deg2rad(i), e, τ, M, deg2rad(ω), deg2rad(Ω), plx)
KeplerianElementsDeg(;a, i, e, τ, M, ω, Ω, plx) = KeplerianElementsDeg(a, i, e, τ, M, ω, Ω, plx)
export KeplerianElementsDeg

# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::KeplerianElements) = print(
    io, """
        $(typeof(elem))
        ─────────────────────────
        a   [au ] = $(round(elem.a,sigdigits=3)) 
        i   [°  ] = $(round(rad2deg(elem.i),sigdigits=3))
        e         = $(round(elem.e,sigdigits=3))
        τ         = $(round(elem.τ,sigdigits=3))
        M   [M⊙ ] = $(round(elem.M,sigdigits=3)) 
        ω   [°  ] = $(round(rad2deg(elem.ω),sigdigits=3))
        Ω   [°  ] = $(round(rad2deg(elem.Ω),sigdigits=3))
        plx [mas] = $(round(elem.plx,sigdigits=3)) 
        ──────────────────────────
        period      [yrs ] : $(round(period(elem)/year2days,digits=1)) 
        distance    [pc  ] : $(round(distance(elem),digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)),sigdigits=3)) 
        ──────────────────────────
        """
)

Base.show(io::IO, elem::KeplerianElements) = print(io,
    "KeplerianElements($(round(elem.a,sigdigits=3)), $(round(elem.i,sigdigits=3)), $(round(elem.e,sigdigits=3)), "*
    "$(round(elem.τ,sigdigits=3)), $(round(elem.M,sigdigits=3)), $(round(elem.ω,sigdigits=3)), "*
    "$(round(elem.Ω,sigdigits=3)), $(round(elem.plx,sigdigits=3)))"
)

# Pretty printing in notebooks as HTML
Base.show(io::IO, ::MIME"text/html", elem::KeplerianElements) = print(
    io, """
        <table style="font-family:monospace; text-align: right">
        <tr><th colspan=3 style="font-family:sans-serif; text-align: left">$(typeof(elem))</th></tr>
        <tr><td rowspan=8>Input</td><td>a   [au] =</td> <td>$(round(elem.a,sigdigits=3))</td></tr>
        <tr><td>i   [°] = </td><td>$(round(rad2deg(elem.i),sigdigits=3))</td></tr>
        <tr><td>e         = </td><td>$(round(elem.e,sigdigits=3))</td></tr>
        <tr><td>τ         = </td><td>$(round(elem.τ,sigdigits=3))</td></tr>
        <tr><td>M   [M⊙] = </td><td>$(round(elem.M,sigdigits=3)) </td></tr>
        <tr><td>ω   [°] = </td><td>$(round(rad2deg(elem.ω),sigdigits=3))</td></tr>
        <tr><td>Ω   [°] = </td><td>$(round(rad2deg(elem.Ω),sigdigits=3))</td></tr>
        <tr><td>plx [mas] = </td><td>$(round(elem.plx,sigdigits=3)) </td></tr>
        <tr><td rowspan=3>Computed</td><td>period      [yrs] : </td><td>$(round(period(elem)/DirectOrbits.year2days,digits=1)) </td></tr>
        <tr><td>distance    [pc] : </td><td>$(round(distance(elem),digits=1)) </td></tr>
        <tr><td>mean motion [°/yr] : </td><td>$(round(rad2deg(DirectOrbits.meanmotion(elem)),sigdigits=3)) </td></tr>
        </table>
        """
)

# Define iterate and length = 1 so that we can broadcast over elements.
Base.length(::AbstractElements) = 1
Base.iterate(elem::AbstractElements) = (elem, nothing)
Base.iterate(::AbstractElements, ::Nothing) = nothing

# ----------------------------------------------------------------------------------------------------------------------
# OrbitSolution
# ----------------------------------------------------------------------------------------------------------------------

"""
    OrbitSolution(
        x, # δ right ascension [milliarcseconds]
        y, # δ declination [milliarcseconds]
        ẋ, # right ascension proper motion anomaly [milliarcseconds/year]
        ẏ, # declination proper motion anomaly [milliarcseconds/year]
        ż, # radial velocity of the *secondary* [m/s]
    )

Represents the secondary's position on the sky in terms of offset from
the primary, its velocity on the sky, and its radial velocity.
"""
struct OrbitSolution{T<:Number}
    x::T
    y::T
    ẋ::T
    ẏ::T
    ż::T
end

# ----------------------------------------------------------------------------------------------------------------------
# System Properties
# ----------------------------------------------------------------------------------------------------------------------

"""
    period(elem)

Period of an orbit in days.
"""
period(elem::KeplerianElements) = elem.T
export period

"""
    distance(elem)

Distance to the system in parsecs.
"""
distance(elem::KeplerianElements) = elem.dist/pc2au
export distance

"""
    meanmotion(elem)

Mean motion, radians per year.
"""
meanmotion(elem::KeplerianElements) = elem.n
export meanmotion 

"""
   periastron(elements, tref=58849)

Compute the MJD of periastron passage most recently after the reference epoch tref.
N.B. mjd of 58849 = 2020-01-01
"""
function periastron(elem::AbstractElements, tref=58849)
    tₚ = elem.τ*period(elem) + tref
    return tₚ
end
export periastron

# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------

function orbitsolve_ν(elem::KeplerianElements, ν)
    # Semi-latus rectum    
    p = elem.a*(1 - elem.e^2) 
    r = p/(1 + elem.e*cos(ν))
    
    # Project back into Cartesian coordinates [AU]
    sin_ω_ν, cos_ω_ν = sincos(elem.ω+ν)
    xₐᵤ = r*(elem.cos_Ω*cos_ω_ν - elem.sin_Ω*sin_ω_ν*elem.cos_i)
    yₐᵤ = r*(elem.sin_Ω*cos_ω_ν + elem.cos_Ω*sin_ω_ν*elem.cos_i)

    # Specific angular momentum 
    h = sqrt(elem.M*p)
    
    # Separation from primary [mas]
    dist⁻¹ = 1/elem.dist
    # Note: use the small angle approximation since arctangent is relatively slow.
    xᵣ = xₐᵤ*dist⁻¹ # atan(xₐᵤ, elem.dist)
    yᵣ = yₐᵤ*dist⁻¹ # atan(yₐᵤ, elem.dist)

    xₘₐₛ = xᵣ * rad2as*oftype(xᵣ,1e3)
    yₘₐₛ = yᵣ * rad2as*oftype(yᵣ,1e3)

    # Factor out common sub-expressions
    A = h * elem.e / (r*p) * sin(ν)
    h_r = h / r

    # TODO: figure out units from first principles.
    ẋₐᵤ = xₐᵤ * A - h_r*(elem.cos_Ω*sin_ω_ν + elem.sin_Ω*cos_ω_ν*elem.cos_i)
    ẏₐᵤ = yₐᵤ * A - h_r*(elem.sin_Ω*sin_ω_ν - elem.cos_Ω*cos_ω_ν*elem.cos_i)

    # Note: use the small angle approximation since arctangent is relatively slow.
    ẋᵣ = ẋₐᵤ*dist⁻¹ # atan(ẋₐᵤ, elem.dist)
    ẏᵣ = ẏₐᵤ*dist⁻¹ # atan(ẏₐᵤ, elem.dist)

    # TODO: investigate source of 2pi factor
    ẋₘₐₛₐ = ẋᵣ * rad2as*oftype(ẋᵣ,1e3) * 2π
    ẏₘₐₛₐ = ẏᵣ * rad2as*oftype(ẏᵣ,1e3) * 2π
    
    # Radial velocity
    żₘₛ = elem.K*(cos_ω_ν + elem.e*elem.cos_ω)

    return OrbitSolution(xₘₐₛ, yₘₐₛ, ẋₘₐₛₐ, ẏₘₐₛₐ, żₘₛ)
end

"""
    orbitsolve(elements, t)

Given a set of orbital elements with a time `t` in days, get the position and
velocity of the secondary body (e.g. planet around a star).

This will output an `OrbitSolution` struct with the following properties:
 - `x`: δ right ascension [milliarcseconds]
 - `y`: δ declination [milliarcseconds]
 - `ẋ`: right ascension proper motion anomaly [milliarcseconds/year]
 - `ẏ`: declination proper motion anomaly [milliarcseconds/year]
 - `ż`: radial velocity of the *secondary* [m/s]

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
@inline function orbitsolve(elem::KeplerianElements{T}, t; tref=58849) where T
    T2 = promote_type(T, typeof(t))
    
    # Epoch of periastron passage
    tₚ = periastron(elem, tref)

    # Mean anomaly    
    MA = meanmotion(elem)/convert(T2, year2days) * (t - tₚ)

    # if !isfinite(MA)
    #     MA = zero(typeof(MA))
    #     @warn "non-finite mean anomaly" maxlog=50
    # end 

    # Compute eccentric anomaly
    EA = _kepler_solver_inline(MA, elem.e)

    # if !isfinite(EA)
    #     EA = MA
    #     @warn "non-finite eccentric anomaly" elem.e maxlog=50
    # end
    
    # Calculate true anomaly
    ν = convert(T2,2)*atan(elem.ν_fact*tan(EA/convert(T2,2)))

    # if !isfinite(ν)
    #     ν = zero(typeof(ν))
    #     @warn "non-finite true anomaly" maxlog=50
    # end

    return orbitsolve_ν(elem, ν)
end
export orbitsolve

# ----------------------------------------------------------------------------------------------------------------------
# Orbital Motion and Separation
# ----------------------------------------------------------------------------------------------------------------------

"""
    raoff(elements, t)

Get the offset from the primary body in Right Ascention in
milliarcseconds at some time `t` in days.
"""
function raoff(elements::AbstractElements, t)
    return orbitsolve(elements, t).x
end
export raoff

"""
    decoff(elements, t)

Get the offset from the primary body in Declination in
milliarcseconds at some time `t` in days.
"""
function decoff(elements::AbstractElements, t)
    return orbitsolve(elements, t).y
end
export decoff

"""
    propmotionanom(elements, t)

Calculate the instantenous proper motion anomaly of a secondary.
"""
function propmotionanom(elements::AbstractElements, t)
    o = orbitsolve(elements, t)
    Δμ_planet = -SVector(o.ẋ, o.ẏ) # milliarcseconds per year
    return Δμ_planet
end
export propmotionanom

"""
    propmotionanom(elements, t, M_planet)

Calculate the instantenous proper motion anomaly on a primary due 
to an orbiting secondary.
"""
function propmotionanom(elements::AbstractElements, t, M_planet)
    M_star = elements.M
    o = orbitsolve(elements, t)
    Δμ_planet = -SVector(o.ẋ, o.ẏ) # milliarcseconds per year
    Δμ_star = Δμ_planet * M_planet / (M_star + M_planet)
    return Δμ_star
end
export propmotionanom

"""
    posangle(elements, t)

Calculate the position angle in radians of a secondary about its primary
from our perspective.
"""
function posangle(elements::AbstractElements, t)
    o = orbitsolve(elements, t)
    return atan(o.y, o.x)
end
export posangle

"""
    projectedseparation(elements, t)

Projected separation in milliarcseconds from the primary at time t (days).
"""
function projectedseparation(elements::AbstractElements, t)
    o = orbitsolve(elements, t)
    return sqrt(o.x^2 + o.y^2)
end
export projectedseparation

"""
    radvel(elements, t)

Get the radial velocity of the *secondary* along the line of sight
at the time `t` in days, in units of m/s.
"""
function radvel(elements::AbstractElements, t)
    return orbitsolve(elements, t).ż
end
export radvel

"""
    radvel(elements, t, M_planet)

Get the radial velocity of the *primary* along the line of sight
at the time `t` in days, in units of m/s.
The mass of the star and planet must have consistent units.
"""
function radvel(elements::AbstractElements, t, M_planet)
    M_star = elements.M
    v_planet = orbitsolve(elements, t).ż
    v_star = -(M_planet/M_star)*v_planet 
    return v_star
end
export radvel

# ----------------------------------------------------------------------------------------------------------------------
# Kepler Equation Solver
# ----------------------------------------------------------------------------------------------------------------------

# The following function is taken directly from AstroLib.jl
# We  remove one invariant check we handle elsewhere and also
# force inlining for about a 5% speedup.
# We also supply analytic gradients for use in autodiff packages.
@inline function _kepler_solver_inline(_M::Real, e::Real)
    # We already handle this invariant
    # @assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    M = rem2pi_safe(_M)
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

# Using implicit differentiation, I found that the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing those here, upstream automatic differentiation libraries will be able
# to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule _kepler_solver_inline(MA, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# We try to support symbolic manipulation using Symbolics.jl, but it's
# not reasonable to use `remp2pi` on a symbolic variable.
# We therefore have a special fallback method for that case. We 
# define it when both packages get loaded by the user using Requires.
@inline rem2pi_safe(x) = rem2pi(x, RoundNearest)
# Define a scale rule to allow autodiff to diff through rem2pi
@scalar_rule rem2pi_safe(x) x

# ----------------------------------------------------------------------------------------------------------------------
# Module Requirements
# ----------------------------------------------------------------------------------------------------------------------

include("transformation.jl")
include("time.jl")
include("plots-recipes.jl")

using Requires
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("makie-recipes.jl")

    # Small patch to allow symbolic tracing through the kepler solver.
    # Mean anomaly must still be in the range [0,2π] for the solution
    # to be valid.
    @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
        @inline rem2pi_safe(x::Symbolics.Num) = x
    end
end

end # module

# ----------------------------------------------------------------------------------------------------------------------