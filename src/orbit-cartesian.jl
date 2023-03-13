


struct CartesianOrbit{T<:Number} <: AbstractOrbit
    # Note: these position and velocity values are in *barycentric* coordinates
    x::T    # AU (increasing to the left)
    y::T    # AU (increasing upwards)
    z::T    # AU (increasing away)
    vx::T   # AU/yr
    vy::T   # AU/yr
    vz::T   # AU/yr
    M::T    # Host mass (solar masses)
    tref::T


    # Orbital properties
    a::T
    e::T
    i::T
    ω::T
    Ω::T
    τ::T

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
    function CartesianOrbit(x,y,z,vx,vy,vz,M,tref=58849)
        x,y,z,vx,vy,vz,M,tref=promote(x,y,z,vx,vy,vz,M,tref)

        # We convert to a keplerian orbit (campbell parameters)
        # so that we can propagate easily.
        # Future work could adapt `orbitsolve` to use an integrator
        # and therefore support multibody physics

        # Compute the specific angular momentum and check for a degenerate orbit
        # r⃗ = @SVector([ x, y, z]) * au2m
        # v⃗ = @SVector([vx,vy,vz]) * au2m / year2sec # au/yr -> m/s
        r⃗ = @SVector([ x, y, z]) * au2m
        v⃗ = @SVector([vx,vy,vz]) * au2m / year2sec # au/yr -> m/s
        h⃗ = r⃗ × v⃗
        h = sqrt(h⃗[1]^2 + h⃗[2]^2 + h⃗[3]^2)

        @show r⃗ v⃗


        # Compute the radius, r, and velocity, v
        r = sqrt(r⃗[1]^2 + r⃗[2]^2 + r⃗[3]^2)
        v² = v⃗[1]^2 + v⃗[2]^2 + v⃗[3]^2

        @show r sqrt(v²)

        # μ = (M + mass*PlanetOrbits.mjup2msol) # * G, where G=1

        # Assume planet mass is negligible
        # μ = M # * G, where G=1 [L^3 T^-2] AU, yr
        
        G = 6.67e-11
        M_kg = M*1.989e+30
        μ = M_kg * G # 1 [L^3 T^-2] AU, yr

        # Compute the specific energy, E, and verify elliptical motion
        # E = v²/2 - μ/r # [L^2 T^-2]
        E = v²/2 - μ/r # [L^2 T^-2]

        # Compute semi-major axis, a
        a_m = -μ/2E
        a = a_m / au2m

        # Compute eccentricity, e
        e² = 1 - h^2/(a_m*μ)

        e  = sqrt(e²)

        # Compute inclination, i
        i = pi- acos(h⃗[3]/h)

        # Compute right ascension of the ascending node, Ω
        # Ω = atan(h⃗[1],-h⃗[2])-π
        Ω0 = atan(h⃗[1],-h⃗[2])
        Ω = rem2pi(-(atan(h⃗[1],-h⃗[2])-π/2),RoundDown)

        # Compute argument of latitude, ω
        ω_plus_ν = atan(z/sin(i), x*cos(Ω0)+y*sin(Ω0))

        # Compute true anomaly, ν 
        oneminusesq = (1 - e^2)
        p = a*oneminusesq
        ν = atan(sqrt(p/μ), p-r)

        # Compute argument of periapse, ω
        ω = (ω_plus_ν - ν ) + π

        # Compute eccentric anomaly, EA
        ν_fact = √((1 + e)/(1 - e)) 
        EA = 2atan(tan(ν/2)/ν_fact)

        # Compute mean motion
        a³ = a^3
        n = sqrt(μ/a)

        # Compute the time of periapse passage, T
        T = tref - 1/n * (EA - e*sin(EA))

        # TODO:
        τ = 0.0

        # Geometric factors involving rotation angles
        sini, cosi = sincos(i)
        sinω, cosω = sincos(ω)
        sinΩ, cosΩ = sincos(Ω)
        ecosω = e*cosω
        esinω = e*sinω
        cosi_cosΩ = cosi*cosΩ
        cosi_sinΩ = cosi*sinΩ

        # Velocity and acceleration semiamplitudes
        periodyrs = √(a³/M)
        period = periodyrs * year2day # period [days]
        J = ((2π*a)/periodyrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
        K = J*au2m*sec2year*sini # radial velocity semiamplitude [m/s]
        A = ((4π^2 * a)/periodyrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]

        return new{typeof(M)}(
            # Passed parameters that define the elements
            x,y,z,vx,vy,vz,M,tref,
            # Converted campbell elements
            a, e, i, ω, Ω, τ,
            # Cached calcuations
            period, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
    end
end


# struct VisualCartesianOrbit{T<:Number,O<:CartesianOrbit{T}} <: AbstractOrbit
#     cartorbit::O
#     plx::T
# end


"""
Represents a `KepOrbit` evaluated to some position.
"""
struct OrbitSolutionCartesian{T<:Number,TEl<:CartesianOrbit} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    EA::T
    sinν_ω::T
    cosν_ω::T
    ecosν::T
    r::T
    t::T
    function OrbitSolutionCartesian(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
        promoted = promote(ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end

export CartesianOrbit
export OrbitSolutionCartesian


# Solve orbit to a new cartesian position given true anomaly
function orbitsolve_ν(elem::CartesianOrbit, ν, EA=2atan(tan(ν/2)/elem.ν_fact), t=_time_from_EA(elem, EA))
    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e*cos(ν)
    r = elem.p/(1 + ecosν)
    return OrbitSolutionCartesian(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
end

# TODO: we can accelerate this since we already know some parameters
"""
Convert an existing orbit object to a CartesianOrbit. 
"""
function CartesianOrbit(o::AbstractOrbit, tref=58849)
    s  = orbitsolve(o,tref)

    x = PlanetOrbits.posx(s)
    y = PlanetOrbits.posy(s)
    z = PlanetOrbits.posz(s)
    vx = PlanetOrbits.xvel(s)
    vy = PlanetOrbits.yvel(s)
    vz = PlanetOrbits.zvel(s)
    return CartesianOrbit(
        x,
        y,
        z,
        vx,
        vy,
        vz,
        o.M,
        tref
    )
end


#=
o = orbit(a=1.0,i=0,ω=π/2,e=0.5,Ω=0,M=1,plx=100.,τ=0.0)

=#
