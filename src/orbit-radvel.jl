"""
    RadialVelocityOrbit(a, e, ω, tp, M)

Represents an orbit of a planet with only the information
retrievable from radial velocity measurements.
That is, without inclination, longitude of ascending node,
or distance to the system.
"""
struct RadialVelocityOrbit{T<:Number} <: AbstractOrbit{T}
    a::T
    e::T
    ω::T
    tp::T
    M::T

    # Physical constants
    T::T
    n::T
    ν_fact::T

    # Geometric factors
    ecosω::T
        
    # Semiamplitude
    K::T
    
    # Inner constructor to enforce invariants and pre-calculate
    # constants from the orbital elements
    function RadialVelocityOrbit(a, e, ω, tp, M, tref=58849)
        # Enforce invariants on user parameters
        M = max(M, zero(M))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        period_days = √(a^3/M)*kepler_year_to_julian_day_conversion_factor
        period_yrs = period_days/year2day_julian
        n = 2π/period_yrs # mean motion
        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(e), typeof(ω),
            typeof(tp), typeof(M), typeof(tref)
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

         # Geometric factors involving rotation angles
        ecosω = e*cos(ω)

        # Velocity and acceleration semiamplitudes
        J = ((2π*a)/period_yrs) / √(1 - e^2) # horizontal velocity semiamplitude [AU/year]
        K = J*au2m*sec2year_julian # radial velocity semiamplitude [m/s]
        new{T}(
            # Passed parameters that define the elements
            a, e, ω, tp, M,
            # Cached calcuations
            period, n, ν_fact,
            # Geometric factors
            ecosω,
            # Semiamplitude
            K
        )
    end
end
# Allow arguments to be specified by keyword
RadialVelocityOrbit(;a, e, ω, tp, M, kwargs...) = RadialVelocityOrbit(a, e, ω, tp, M)
export RadialVelocityOrbit

period(elem::RadialVelocityOrbit) = elem.T
meanmotion(elem::RadialVelocityOrbit) = elem.n
eccentricity(o::RadialVelocityOrbit) = o.e
totalmass(o::RadialVelocityOrbit) = o.M
semimajoraxis(o::RadialVelocityOrbit) = o.a
_trueanom_from_eccanom(o::RadialVelocityOrbit, EA) =2*atan(o.ν_fact*tan(EA/2))
periastron(o::RadialVelocityOrbit) = o.tp
semiamplitude(elem::RadialVelocityOrbit) = elem.K

# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::RadialVelocityOrbit) = print(
io, """
    $(typeof(elem))
    ─────────────────────────
    a   [au ] = $(round(elem.a, sigdigits=3))
    e         = $(round(elem.e, sigdigits=8))
    ω   [°  ] = $(round(rad2deg(elem.ω), sigdigits=3))
    tp         = $(round(elem.tp, sigdigits=3))
    M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
    ──────────────────────────
    period        [yrs ] : $(round(period(elem)*day2year_julian, digits=1)) 
    mean motion   [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
    semiamplitude₂ [m/s] : $(round(semiamplitude(elem), digits=1)) 
    ──────────────────────────
    """
)


"""
    orbitsolve_ν(elem, ν)

Solve a keplerian orbit from a given true anomaly [rad].
See orbitsolve for the same function accepting a given time.
"""
function orbitsolve_ν(elem::RadialVelocityOrbit, ν, EA=2atan(tan(ν/2)/elem.ν_fact), t=_time_from_EA(elem, EA))
    cosν_ω = cos(elem.ω + ν)
    return OrbitSolutionRadialVelocity(elem, ν, EA, cosν_ω, t)
end


"""
Represents a `RadialVelocityOrbit` evaluated to some position.
"""
struct OrbitSolutionRadialVelocity{T<:Number,TEl<:RadialVelocityOrbit} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    EA::T
    cosν_ω::T
    t::T
    function OrbitSolutionRadialVelocity(elem, ν, EA, cosν_ω, t)
        promoted = promote(ν, EA, cosν_ω, t)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end
export OrbitSolutionRadialVelocity
soltime(os::OrbitSolutionRadialVelocity) = os.t
_solution_type(::Type{RadialVelocityOrbit}) = OrbitSolutionRadialVelocity
