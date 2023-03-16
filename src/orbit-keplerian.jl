
"""
    KepOrbit(
        a, # semi-major axis [AU]
        e, # eccentricity
        i, # inclination [rad]
        ω, # argument of periapsis [rad]
        Ω, # longitude of ascending node [rad]
        τ, # epoch of periastron passage at MJD=0
        M, # mass of primary [M⊙]
    )

Represents the Keplerian elements of a secondary body orbiting a primary.
Use the traditional Campbell parameterization.
Values can be specified by keyword argument or named tuple for convenience.
"""
struct KepOrbit{T<:Number} <: AbstractOrbit{T}

    # Orbital properties
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    τ::T
    M::T

    # Reference epoch
    tref::T

    # Physical constants
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
    function KepOrbit(a, e, i, ω, Ω, τ, M, tref=58849)

        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        τ = mod(τ, one(τ))
        M = max(M, zero(M))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
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
            typeof(Ω), typeof(τ), typeof(M), typeof(tref)
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
            a, e, i, ω, Ω, τ, M, tref,
            # Cached calcuations
            period, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
    end
end

# Allow arguments to be specified by keyword
KepOrbit(;a, e, i, ω, Ω, τ, M, tref=58849, kwargs...) = KepOrbit(a, e, i, ω, Ω, τ, M, tref)
export KepOrbit


"""
    astuple(elements)

Return the parameters of a KepOrbit value as a tuple.
"""
function astuple(elem::KepOrbit)
return (;elem.a, elem.e, elem.i, elem.ω, elem.Ω, elem.τ, elem.M)
end
export astuple

# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::KepOrbit) = print(
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
    period      [yrs ] : $(round(period(elem)*day2year, digits=1)) 
    mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
    ──────────────────────────
    """
)

Base.show(io::IO, elem::KepOrbit) = print(io,
"KepOrbit($(round(elem.a, sigdigits=3)), $(round(elem.e, sigdigits=3)), $(round(elem.i, sigdigits=3)), "*
"$(round(elem.ω, sigdigits=3)), $(round(elem.Ω, sigdigits=3)), $(round(elem.τ, sigdigits=3)), "*
"$(round(elem.M, sigdigits=3)))"
)


"""
Represents a `KepOrbit` evaluated to some position.
"""
struct OrbitSolutionKep{T<:Number,TEl<:KepOrbit} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    EA::T
    sinν_ω::T
    cosν_ω::T
    ecosν::T
    r::T
    t::T
    function OrbitSolutionKep(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
        promoted = promote(ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end

_solution_type(::Type{KepOrbit}) = OrbitSolutionKep


period(elem::KepOrbit) = elem.T
meanmotion(elem::KepOrbit) = elem.n
eccentricity(o::KepOrbit) = o.e
hostmass(o::KepOrbit) = o.M
_trueanom_from_eccanom(o::KepOrbit, EA) =2*atan(o.ν_fact*tan(EA/2))
function periastron(elem::KepOrbit)
    tₚ = elem.τ*period(elem) + elem.tref
    return tₚ
end
semiamplitude(elem::KepOrbit) = elem.K


# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------
function orbitsolve_ν(elem::KepOrbit, ν, EA=2atan(tan(ν/2)/elem.ν_fact), t=_time_from_EA(elem, EA))
    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    r = elem.p/(1 + ecosν)
    return OrbitSolutionKep(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
end
soltime(os::OrbitSolutionKep) = os.t
