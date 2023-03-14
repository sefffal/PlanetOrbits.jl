


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
    function CartesianOrbit(x, y, z, vx, vy, vz, M, tref)
        if M isa Integer
            M = float(M)
        end
        x, y, z, vx, vy, vz, M, tref = promote(x, y, z, vx, vy, vz, M, tref)

        r⃗ = @SVector([x, y, z])
        v⃗ = @SVector([vx, vy, vz]) ./ 2π # ??
        h⃗ = r⃗ × v⃗
        h = norm(h⃗)

        # Compute the radius, r, and velocity, v
        r = norm(r⃗)

        tmp⃗ = v⃗ × h⃗
        e⃗ = tmp⃗ / M - r⃗ / r
        e = norm(e⃗)

        if e > 1
            error(lazy"Unbound orbit: e is greater than 1 (e=$e)")
        end

        a = (h⃗ ⋅ h⃗) / (M * (1 - e^2))

        i⃗ = @SVector([1.0, 0.0, 0.0])
        j⃗ = @SVector([0.0, 1.0, 0.0])
        k⃗ = @SVector([0.0, 0.0, 1.0])

        i = pi - acos((k⃗ ⋅ h⃗) / h)

        n⃗ = k⃗ × h⃗
        n = norm(n⃗)
        # Ω = acos((i⃗ ⋅ n⃗) / n)
        Ω = asin((i⃗ ⋅ n⃗) / n)
        if n⃗ ⋅ j⃗ < 0
            Ω = 2π - Ω
        end


        ω = acos((n⃗ ⋅ e⃗) / (n * e))
        if e⃗ ⋅ k⃗ < 0
            ω = 2π - ω
        end

        arg = (e⃗ ⋅ r⃗) / (e * r)
        if 1 < arg < 1+3eps()
            arg = one(arg)
        end
        θ = acos(arg)
        if θ < 0
            θ = 2π - θ
        end

        EA = acos((e + cos(θ)) / (1 + e * cos(θ)))
        if π < θ < 2π
            EA = 2π - EA
        end
        MA = EA - e * sin(EA)


        tₚ = MA / n * PlanetOrbits.year2day

        a³ = a^3
        oneminusesq = (1 - e^2)
        ν_fact = √((1 + e) / (1 - e))
        p = a * oneminusesq
        n = 2π / √(a^3 / M) # mean motion

        periodyrs = √(a³ / M)
        period = periodyrs * year2day # period [days]

        τ = rem(tₚ / period, 1, RoundDown)

        # Geometric factors involving rotation angles
        sini, cosi = sincos(i)
        sinω, cosω = sincos(ω)
        sinΩ, cosΩ = sincos(Ω)
        ecosω = e * cosω
        esinω = e * sinω
        cosi_cosΩ = cosi * cosΩ
        cosi_sinΩ = cosi * sinΩ

        # Velocity and acceleration semiamplitudes
        J = ((2π * a) / periodyrs) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
        K = J * au2m * sec2year * sini # radial velocity semiamplitude [m/s]
        A = ((4π^2 * a) / periodyrs^2) / oneminusesq^2 # horizontal acceleration semiamplitude [AU/year^2]



        orbit = new{typeof(M)}(
            # Passed parameters that define the elements
            x, y, z, vx, vy, vz, M, tref,
            # Converted campbell elements
            a, e, i, ω, Ω, τ,
            # Cached calcuations
            period, n, ν_fact, p,
            # Geometric factors
            cosi, sini, cosΩ, sinΩ, ecosω, esinω, cosi_cosΩ, cosi_sinΩ,
            # Semiamplitudes
            J, K, A
        )
        # Test
        # os = orbitsolve(orbit, tref)
        # @show posx(os), x
        # @show posy(os), y
        # @show posz(os), z

        return orbit
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
function orbitsolve_ν(elem::CartesianOrbit, ν, EA=2atan(tan(ν / 2) / elem.ν_fact), t=_time_from_EA(elem, EA))
    sinν_ω, cosν_ω = sincos(elem.ω + ν)
    ecosν = elem.e * cos(ν)
    r = elem.p / (1 + ecosν)
    return OrbitSolutionCartesian(elem, ν, EA, sinν_ω, cosν_ω, ecosν, r, t)
end

# TODO: we can accelerate this since we already know some parameters
"""
Convert an existing orbit object to a CartesianOrbit. 
"""
function CartesianOrbit(os::AbstractOrbitSolution)
    x = PlanetOrbits.posx(os)
    y = PlanetOrbits.posy(os)
    z = PlanetOrbits.posz(os)
    vx = PlanetOrbits.xvel(os)
    vy = PlanetOrbits.yvel(os)
    vz = PlanetOrbits.zvel(os)
    return CartesianOrbit(
        x,
        y,
        z,
        vx,
        vy,
        vz,
        os.elem.M,
        os.t
        # os.elem.tref
    )
end


#=
o = orbit(a=1.0,i=0,ω=π/2,e=0.5,Ω=0,M=1,plx=100.,τ=0.0)

##
x  = -1.05
y  = 3.782338790704024e-16
z  = -2.5048146051777413e-17
vx = -3.490253699036788e-16 * 2pi
vy = -0.9464377445249709 * 2pi
vz = 0.09496052074620637 * 2pi
oc = CartesianOrbit(x,y,z,vx,vy,vz,1,0)

sc = orbitsolve(oc, 0)
x  = PlanetOrbits.posx(sc)
y  = PlanetOrbits.posy(sc)
z  = PlanetOrbits.posz(sc)
vx = PlanetOrbits.xvel(sc)
vy = PlanetOrbits.yvel(sc)
vz = PlanetOrbits.zvel(sc)
oc2= CartesianOrbit(x,y,z,vx,vy,vz,1,0)
plot(orbitsolve(oc,0));plot!(orbitsolve(oc2,0))
# i appears to be going backwards


o = orbit(
    a = 1,
    i = π/4,
    Ω = 0.001,
    ω = π/4,
    e = 0.5,
    τ = 0.0,
    M = 1,
    tref=0
)
oc3 = CartesianOrbit(orbitsolve(o,0))
oc4 = CartesianOrbit(orbitsolve(oc3,0))
oc5 = CartesianOrbit(orbitsolve(oc4,0))

plot(orbitsolve(o,0),label="o", lw=2, ls=:dash, color=1)
plot!(orbitsolve(oc3,0), label="oc3", color=2)
plot!(orbitsolve(oc4,0), label="oc4", color=3)
plot!(orbitsolve(oc5,0), label="oc5", color=4)

I think there are three things left:
* coordinate system is flipped
* velocity is scaled by 2pi for all orbits (huh? Need to check radial velocity math for all orbits)
* something about the dates / times stamping of CartesianOrbit is not making sense.
=#
