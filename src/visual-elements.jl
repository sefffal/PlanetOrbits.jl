
"""
VisualElements(
    a, # semi-major axis [AU]
    e, # eccentricity
    i, # inclination [rad]
    ω, # argument of periapsis [rad]
    Ω, # longitude of ascending node [rad]
    τ, # epoch of periastron passage at MJD=0
    M, # mass of primary [M⊙]
    plx, # parallax [mas]; defines the distance to the primary
)

Represents the Keplerian elements of a secondary body orbiting a primary.
Values can be specified by keyword argument or named tuple for convenience.

See also `VisualElementsDeg` for a convenience constructor accepting
units of degrees instead of radians for `i`, `ω`, and `Ω`.
"""
struct VisualElements{T<:Number} <: AbstractOrbit

    # Orbital properties
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    τ::T
    M::T
    plx::T

    # Physical constants
    dist::T
    T::T
    n::T
    ν_fact::T
    p::T

    # Geometric factors
    cosi::T
    sini::T
    cosΩ::T
    sinΩ::T
    ecosω::T
    esinω::T
    cosi_cosΩ::T
    cosi_sinΩ::T

    # Semiamplitudes
    J::T
    K::T
    A::T

    # Inner constructor to enforce invariants and pre-calculate
    # constants from the orbital elements
    function VisualElements(a, e, i, ω, Ω, τ, M, plx)

        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        τ = mod(τ, one(τ))
        M = max(M, zero(M))
        plx = max(plx, zero(plx))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        dist = 1000/plx * pc2au # distance [AU]
        rootacubeoverm = √(a^3/M)
        periodyrs = rootacubeoverm
        period = periodyrs * year2day # period [days]
        n = 2π/√(a^3/M) # mean motion
        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor
        oneminusesq = (1 - e^2)
        p = a*oneminusesq # semi-latus rectum [AU]

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(e), typeof(i), typeof(ω),
            typeof(Ω), typeof(τ), typeof(M), typeof(plx),
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        # Geometric factors involving rotation angles
        sini, cosi = sincos(i)
        sinω, cosω = sincos(ω)
        sinΩ, cosΩ = sincos(Ω)
        ecosω = e*cosω
        esinω = e*sinω
        cosi_cosΩ = cosi*cosΩ
        cosi_sinΩ = cosi*sinΩ

        # Velocity and acceleration semiamplitudes
        J = ((2π*a)/periodyrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
        K = J*au2m*sec2year*sini # radial velocity semiamplitude [m/s]
        A = ((4π^2 * a)/periodyrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]

        new{T}(
            # Passed parameters that define the elements
            a, e, i, ω, Ω, τ, M, plx,
            # Cached calcuations
            dist, period, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
    end
end

# Allow arguments to be specified by keyword
VisualElements(;a, e, i, ω, Ω, τ, M, plx) = VisualElements(a, e, i, ω, Ω, τ, M, plx)
# Allow arguments to be specified by named tuple
VisualElements(nt) = VisualElements(nt.a, nt.e, nt.i, nt.ω, nt.Ω, nt.τ, nt.M, nt.plx)
export VisualElements

"""
VisualElementsDeg(a, e, i, ω, Ω, τ, M, plx)

A convenience function for constructing VisualElements where
`i`, `ω`, and `Ω` are provided in units of degrees instead of radians.
"""
VisualElementsDeg(a, e, i, ω, Ω, τ, M, plx) = VisualElements(a, e, deg2rad(i), deg2rad(ω), deg2rad(Ω), τ, M, plx)
VisualElementsDeg(;a, e, i, ω, Ω, τ, M, plx) = VisualElementsDeg(a, e, i, ω, Ω, τ, M, plx)
VisualElementsDeg(nt) = VisualElementsDeg(nt.a, nt.e, nt.i, nt.ω, nt.Ω, nt.τ, nt.M, nt.plx)
export VisualElementsDeg

"""
astuple(elements)

Return the parameters of a VisualElements value as a tuple.
"""
function astuple(elem::VisualElements)
return (;elem.a, elem.e, elem.i, elem.ω, elem.Ω, elem.τ, elem.M, elem.plx)
end
export astuple

# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::VisualElements) = print(
io, """
    $(typeof(elem))
    ─────────────────────────
    a   [au ] = $(round(elem.a, sigdigits=3))
    e         = $(round(elem.e, sigdigits=3))
    i   [°  ] = $(round(rad2deg(elem.i), sigdigits=3))
    ω   [°  ] = $(round(rad2deg(elem.ω), sigdigits=3))
    Ω   [°  ] = $(round(rad2deg(elem.Ω), sigdigits=3))
    τ         = $(round(elem.τ, sigdigits=3))
    M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
    plx [mas] = $(round(elem.plx, sigdigits=3)) 
    ──────────────────────────
    period      [yrs ] : $(round(period(elem)*day2year, digits=1)) 
    distance    [pc  ] : $(round(distance(elem), digits=1)) 
    mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
    ──────────────────────────
    """
)

Base.show(io::IO, elem::VisualElements) = print(io,
"VisualElements($(round(elem.a, sigdigits=3)), $(round(elem.e, sigdigits=3)), $(round(elem.i, sigdigits=3)), "*
"$(round(elem.ω, sigdigits=3)), $(round(elem.Ω, sigdigits=3)), $(round(elem.τ, sigdigits=3)), "*
"$(round(elem.M, sigdigits=3)), $(round(elem.plx, sigdigits=3)))"
)

# Pretty printing in notebooks as HTML
Base.show(io::IO, ::MIME"text/html", elem::VisualElements) = print(
io, """
    <table style="font-family:monospace; text-align: right">
    <tr><th colspan=3 style="font-family:sans-serif; text-align: left">$(typeof(elem))</th></tr>
    <tr><td rowspan=8>Input</td><td>a   [au] =</td> <td>$(round(elem.a, sigdigits=3))</td></tr>
    <tr><td>e         = </td><td>$(round(elem.e, sigdigits=3))</td></tr>
    <tr><td>i   [°] = </td><td>$(round(rad2deg(elem.i), sigdigits=3))</td></tr>
    <tr><td>ω   [°] = </td><td>$(round(rad2deg(elem.ω), sigdigits=3))</td></tr>
    <tr><td>Ω   [°] = </td><td>$(round(rad2deg(elem.Ω), sigdigits=3))</td></tr>
    <tr><td>τ         = </td><td>$(round(elem.τ, sigdigits=3))</td></tr>
    <tr><td>M   [M⊙] = </td><td>$(round(elem.M, sigdigits=3)) </td></tr>
    <tr><td>plx [mas] = </td><td>$(round(elem.plx, sigdigits=3)) </td></tr>
    <tr><td rowspan=3>Computed</td><td>period      [yrs] : </td><td>$(round(period(elem)*PlanetOrbits.day2year, digits=1)) </td></tr>
    <tr><td>distance    [pc] : </td><td>$(round(distance(elem), digits=1)) </td></tr>
    <tr><td>mean motion [°/yr] : </td><td>$(round(rad2deg(PlanetOrbits.meanmotion(elem)), sigdigits=3)) </td></tr>
    </table>
    """
)

# Define iterate and length = 1 so that we can broadcast over elements.
Base.length(::AbstractOrbit) = 1
Base.iterate(elem::AbstractOrbit) = (elem, nothing)
Base.iterate(::AbstractOrbit, ::Nothing) = nothing




"""
Represents a `VisualElements` evaluated to some position.
"""
struct OrbitSolutionVisual{T<:Number,TEl<:VisualElements} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    EA::T
    sinν_ω::T
    cosν_ω::T
    ecosν::T
    r::T
    function OrbitSolutionVisual(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r)
        promoted = promote(ν, EA, sinν_ω, cosν_ω, ecosν, r)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end
export OrbitSolutionVisual

# Printing
# Base.show(io::IO, os::AbstractOrbitSolution) = print(io,
#     "AbstractOrbitSolution(x = $(round(os.x, sigdigits=3)), y = $(round(os.y, sigdigits=3)), "*
#     "ẋ = $(round(os.ẋ, sigdigits=3)), ẏ = $(round(os.ẏ, sigdigits=3)), ż = $(round(os.ż, sigdigits=3)), "*
#     "ẍ = $(round(os.ẍ, sigdigits=3)), ÿ = $(round(os.ÿ, sigdigits=3)))"
# )

period(elem::VisualElements) = elem.T
distance(elem::VisualElements) = elem.dist*au2pc
meanmotion(elem::VisualElements) = elem.n
function periastron(elem::VisualElements, tref=58849)
    tₚ = elem.τ*period(elem) + tref
    return tₚ
end
semiamplitude(elem::VisualElements) = elem.K

# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------
function orbitsolve_ν(elem::VisualElements, ν; EA=2atan(tan(ν/2)/elem.ν_fact))
    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    r = elem.p/(1 + ecosν)
    return OrbitSolutionVisual(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r)
end

function posx(o::OrbitSolutionVisual)
    xcart = o.r*(o.cosν_ω*o.elem.sinΩ + o.sinν_ω*o.elem.cosi*o.elem.cosΩ) # [AU]
    return xcart
end

function posy(o::OrbitSolutionVisual)
    ycart = o.r*(o.cosν_ω*o.elem.cosΩ - o.sinν_ω*o.elem.cosi*o.elem.sinΩ) # [AU]
    return ycart
end


function posz(o::OrbitSolutionVisual)
    zcart = o.r*(o.cosν_ω*o.elem.sini) # [AU]
    return zcart
end


function raoff(o::OrbitSolutionVisual)
    xcart = posx(o) # [AU]
    cart2angle = rad2as*oftype(xcart, 1e3)/o.elem.dist
    xang = xcart*cart2angle # [mas]
    return xang
end

function decoff(o::OrbitSolutionVisual)
    ycart = posy(o) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.elem.dist
    yang = ycart*cart2angle # [mas]
    return yang
end

function pmra(o::OrbitSolutionVisual)
    ẋcart = o.elem.J*(o.elem.cosi_cosΩ*(o.cosν_ω + o.elem.ecosω) - o.elem.sinΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.elem.dist
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang
end

function pmdec(o::OrbitSolutionVisual)
    ẏcart = -o.elem.J*(o.elem.cosi_sinΩ*(o.cosν_ω + o.elem.ecosω) + o.elem.cosΩ*(o.sinν_ω + o.elem.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.elem.dist
    ẏang = ẏcart*cart2angle # [mas/year]
    return ẏang
end


function accra(o::AbstractOrbitSolution)
    ẍcart = -o.elem.A*(1 + o.ecosν)^2 * (o.elem.cosi_cosΩ*o.sinν_ω + o.elem.sinΩ*o.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ẍcart, 1e3)/o.elem.dist
    ẍang = ẍcart*cart2angle # [mas/year^2] 
    return ẍang
end
function accdec(o::AbstractOrbitSolution)
    ÿcart = o.elem.A*(1 + o.ecosν)^2 * (o.elem.cosi_sinΩ*o.sinν_ω - o.elem.cosΩ*o.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ÿcart, 1e3)/o.elem.dist
    ÿang = ÿcart*cart2angle # [mas/year^2] 
    return ÿang
end
