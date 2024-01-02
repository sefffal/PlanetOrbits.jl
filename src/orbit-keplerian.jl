
"""
    KepOrbit(
        a, # semi-major axis [AU]
        e, # eccentricity
        i, # inclination [rad]
        ω, # argument of periapsis [rad]
        Ω, # longitude of ascending node [rad]
        tp, # epoch of periastron passage at MJD=0
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
    tp::T
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
    function KepOrbit(a, e, i, ω, Ω, tp, M, tref=58849)

        # Enforce invariants on user parameters
        # a = max(a, zero(a))
        # e = max(zero(e), min(e, one(e)))
        M = max(M, zero(M))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        if e < 1
            rootacubeoverm = √(a^3/M)
            periodyrs = rootacubeoverm
            period = periodyrs * year2day # period [days]
            n = 2π/√(a^3/M) # mean motion
        else
            period = Inf
            n = √(M/-a^3) # mean motion
        end

        if e < 1
            ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor
        else
            ν_fact = √((1 + e)/(e - 1)) # true anomaly prefactor
        end
        oneminusesq = (1 - e^2)
        p = a*oneminusesq # semi-latus rectum [AU]

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(e), typeof(i), typeof(ω),
            typeof(Ω), typeof(tp), typeof(M), typeof(tref)
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        # Geometric factors involving rotation angles
        if e < 1
            sini, cosi = sincos(i)
            sinω, cosω = sincos(ω)
            sinΩ, cosΩ = sincos(Ω)
        else
            # Cached angular quantities become hyperbolic 
            # for hyperbolic orbits.
            sini, cosi = sinh(i), cosh(i)
            sinω, cosω = sinh(ω), cosh(ω)
            sinΩ, cosΩ = sinh(Ω), cosh(Ω)
        end
        ecosω = e*cosω
        esinω = e*sinω
        cosi_cosΩ = cosi*cosΩ
        cosi_sinΩ = cosi*sinΩ

        if e < 1

            # Velocity and acceleration semiamplitudes
            J = ((2π*a)/periodyrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
            K = J*au2m*sec2year*sini # radial velocity semiamplitude [m/s]
            A = ((4π^2 * a)/periodyrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]
        else
            @warn "velocity and acceleration not implemented for ecc >= 1 yet"
            J = K = A = 0.0
        end
        new{T}(
            # Passed parameters that define the elements
            a, e, i, ω, Ω, tp, M, tref,
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
KepOrbit(;a, e, i, ω, Ω, tp, M, tref=58849, kwargs...) = KepOrbit(a, e, i, ω, Ω, tp, M, tref)
export KepOrbit


"""
    astuple(elements)

Return the parameters of a KepOrbit value as a tuple.
"""
function astuple(elem::KepOrbit)
    return (;elem.a, elem.e, elem.i, elem.ω, elem.Ω, elem.tp, elem.M)
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
    tp  [day] = $(round(elem.tp, sigdigits=3))
    M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
    period      [yrs ] : $(round(period(elem)*day2year, digits=1)) 
    mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
    ──────────────────────────
    """
)

Base.show(io::IO, elem::KepOrbit) = print(io,
"KepOrbit($(round(elem.a, sigdigits=3)), $(round(elem.e, sigdigits=3)), $(round(elem.i, sigdigits=3)), "*
"$(round(elem.ω, sigdigits=3)), $(round(elem.Ω, sigdigits=3)), $(round(elem.tp, sigdigits=3)), "*
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
totalmass(o::KepOrbit) = o.M
inclination(o::KepOrbit) = o.i
semimajoraxis(o::KepOrbit) = o.a
function _trueanom_from_eccanom(o::KepOrbit, EA)
    if o.e < 1
        ν = 2*atan(o.ν_fact*tan(EA/2))
    else
        # true anomaly prefactor changed in constructor if hyperbolic
        ν = 2*atan(o.ν_fact*tanh(EA/2))
    end
    return ν
end
function periastron(elem::KepOrbit)
    return elem.tp
    # if eccentricity(elem) < 1
    #     tₚ = elem.τ*period(elem) + elem.tref
    # else
    #     @warn "TODO: periastron for hyperbolic"
    #     tₚ = elem.τ + elem.tref
    # end
    # return tₚ
end
semiamplitude(elem::KepOrbit) = elem.K


# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------
function EA_from_ν(elem::KepOrbit, ν)
    if elem.e < 1
        EA = 2atan(tan(ν/2)/elem.ν_fact)
    else
        println("TODO: confirm this h is in the right location")
        EA = 2atanh(tan(ν/2)/elem.ν_fact)
    end
    return EA
end
function orbitsolve_ν(elem::KepOrbit, ν, EA=EA_from_ν(elem, ν), t=_time_from_EA(elem, EA))

    @show EA t ν

    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    @show ecosν
    r = elem.p/(1 + ecosν)
    @show r


    return OrbitSolutionKep(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
end
soltime(os::OrbitSolutionKep) = os.t
