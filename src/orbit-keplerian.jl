
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
    function KepOrbit(a, e, i, ω, Ω, tp, M)

        # Enforce invariants on user parameters
        # a = max(a, zero(a))
        # e = max(zero(e), min(e, one(e)))
        M = max(M, zero(M))
        i = rem(i, π, RoundDown)
        Ω = rem2pi(Ω, RoundDown)

        if e >= 1 && a > 0
            @warn "Negative semi-major is required for hyperbolic (e>1) orbits. Flipping sign (maxlog=1)." maxlog=1
            a = -a
        end

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        if e < 1
            period_days = √(a^3/M)*kepler_year_to_julian_day_conversion_factor
            period_yrs = period_days/year2day_julian
            n = 2π/period_yrs # mean motion
        else
            period_days = Inf
            # TODO: Need to confirm where this 2pi is coming from 
            # mean motion
            n = 2π * √(M/-a^3)*kepler_year_to_julian_day_conversion_factor/year2day_julian
            # n = √(M/-a^3) # mean motion
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
            typeof(Ω), typeof(tp), typeof(M),
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

        if e < 1
            # Velocity and acceleration semiamplitudes
            J = ((2π*a)/period_yrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
            K = J*au2m*sec2year_julian*sini # radial velocity semiamplitude [m/s]
            A = ((4π^2 * a)/period_yrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]
        else
            J = -((2π*a)/√(M/-a^3)) / √(-oneminusesq) # horizontal velocity semiamplitude [AU/year]
            K = J*au2m*sec2year_julian*sini # radial velocity semiamplitude [m/s]
            # TODO: acceleration not verified for ecc >= 1 yet. Results will be silently wrong.
            A = ((4π^2 * a)/(M/-a^3)) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]
        end
        new{T}(
            # Passed parameters that define the elements
            a, e, i, ω, Ω, tp, M,
            # Cached calcuations
            period_days, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
    end
end

# Allow arguments to be specified by keyword
KepOrbit(;a, e, i, ω, Ω, tp, M, kwargs...) = KepOrbit(a, e, i, ω, Ω, tp, M)
export KepOrbit

Base.convert(::Type{KepOrbit{T2}}, orb::KepOrbit{T1}) where {T1<:Number,T2<:Number} = KepOrbit(
    T2(orb.a),
    T2(orb.e),
    T2(orb.i),
    T2(orb.ω),
    T2(orb.Ω),
    T2(orb.tp),
    T2(orb.M)
)

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
    e         = $(round(elem.e, sigdigits=8))
    i   [°  ] = $(round(rad2deg(elem.i), sigdigits=3))
    ω   [°  ] = $(round(rad2deg(elem.ω), sigdigits=3))
    Ω   [°  ] = $(round(rad2deg(elem.Ω), sigdigits=3))
    tp  [day] = $(round(elem.tp, sigdigits=3))
    M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
    period      [yrs ] : $(round(period(elem)*day2year_julian, digits=1)) 
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
    function OrbitSolutionKep(elem::KepOrbit{T}, ν, EA, sinν_ω, cosν_ω, ecosν, r, t) where T
        promoted = promote(ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
        return new{T,KepOrbit{T}}(elem, promoted...)
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
        β = eccentricity(o) / (1+sqrt(1-eccentricity(o)^2))
        sea, cea = sincos(EA)
        ν = EA + 2atan(β*sea/(1-β*cea))
    else
        # true anomaly prefactor changed in constructor if hyperbolic
        ν = 2*atan(o.ν_fact*tanh(EA/2))
    end
    return ν
end
function periastron(elem::KepOrbit)
    return elem.tp
end
semiamplitude(elem::KepOrbit) = elem.K


# ----------------------------------------------------------------------------------------------------------------------
# Solve Orbit in Cartesian Coordinates
# ----------------------------------------------------------------------------------------------------------------------
function EA_from_ν(elem::KepOrbit, ν)
    if elem.e < 1
        EA = 2atan(tan(ν/2)/elem.ν_fact)
    else
        arg = tan(ν/2)/elem.ν_fact
        if 1 < arg < 1+sqrt(eps(arg))
            arg = oftype(arg, 1)
        end
        EA = 2atanh(arg)
    end
    return EA
end
function orbitsolve_ν(elem::KepOrbit, ν, EA=EA_from_ν(elem, ν), t=_time_from_EA(elem, EA))

    # @show EA t ν

    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    # @show ecosν
    r = elem.p/(1 + ecosν)
    # @show r


    return OrbitSolutionKep(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
end
soltime(os::OrbitSolutionKep) = os.t
