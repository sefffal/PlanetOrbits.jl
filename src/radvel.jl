
struct RadialVelocityElements{T<:Number} <: AbstractOrbit
    a::T
    e::T
    ω::T
    τ::T
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
    function RadialVelocityElements(a, e, ω, τ, M)
        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        τ = mod(τ, one(τ))
        M = max(M, zero(M))

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        periodyrs = √(a^3/M)
        period = periodyrs * year2day # period [days]
        n = 2π/√(a^3/M) # mean motion
        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor

        # Get type of parameters
        T = promote_type(
            typeof(a), typeof(e), typeof(ω),
            typeof(τ), typeof(M)
        )

        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

         # Geometric factors involving rotation angles
        ecosω = e*cos(ω)

        # Velocity and acceleration semiamplitudes
        J = ((2π*a)/periodyrs) / √(1 - e^2) # horizontal velocity semiamplitude [AU/year]
        K = J*au2m*sec2year # radial velocity semiamplitude [m/s]
        new{T}(
            # Passed parameters that define the elements
            a, e, ω, τ, M, 
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
RadialVelocityElements(;a, e, ω, τ, M) = RadialVelocityElements(a, e, ω, τ, M)
# Allow arguments to be specified by named tuple
RadialVelocityElements(nt) = RadialVelocityElements(nt.a, nt.e, nt.ω, nt.τ, nt.M)
export RadialVelocityElements


period(elem::RadialVelocityElements) = elem.T
distance(elem::RadialVelocityElements) = elem.dist*au2pc
meanmotion(elem::RadialVelocityElements) = elem.n
function periastron(elem::RadialVelocityElements, tref=58849)
    tₚ = elem.τ*period(elem) + tref
    return tₚ
end
semiamplitude(elem::RadialVelocityElements) = elem.K


"""
    orbitsolve_ν(elem, ν)

Solve a keplerian orbit from a given true anomaly [rad].
See orbitsolve for the same function accepting a given time.
"""
function orbitsolve_ν(elem::RadialVelocityElements, ν)
    cosν_ω = cos(elem.ω + ν)
    return OrbitSolutionRadialVelocity(elem, ν, cosν_ω)
end

struct OrbitSolutionRadialVelocity{T<:Number,TEl<:RadialVelocityElements} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    cosν_ω::T
    function OrbitSolutionRadialVelocity(elem, ν, cosν_ω,)
        promoted = promote(ν, cosν_ω)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end
export OrbitSolutionRadialVelocity