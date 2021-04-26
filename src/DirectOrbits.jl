"""
# DirectOrbits
A package for calculating orbits in the context of direct imaging.

"""
module DirectOrbits

using LinearAlgebra
# using CoordinateTransformations
using StaticArrays
using Roots # For solving for eccentric anomaly
# import Dates
# import Base.inv

using AstroLib: kepler_solver


const mas2rad = 4.8481368E-9
const rad2as = 206265
const pc2au = 206265
const au2m = 1.495978707e11
const year2days = 365.2422


abstract type AbstractElements end

"""
    Orbit(
        a=1.0, # semi-major axis, AU
        i=π/2, # inclination, radians
        e=0.1, # eccentricity
        τ=π/2, # fraction of elements past periastron passage at MJD=0,
        μ=1.0, # graviational parameter, solar masses
        ω=π/2, # argument of periapsis
        Ω=π/2, # longitude of the ascending node
        plx=10.1, # paralax in milliarcseconds. Defines the distance to the object
    )

Represents one object's Keplerian elementsal elements. Values can be specified
by keyword argument for convinience, or kep2cart for efficiency.

See also `KeplerianElementsDeg` for a convinience constructor accepting
units of degrees instead of radians.
"""
struct KeplerianElements{T<:Number} <: AbstractElements

    # Orbital properties
    a::T
    i::T
    e::T
    τ::T
    μ::T
    ω::T
    Ω::T
    plx::T

    # Cached constants for these elements
    dist::T
    T::T
    n::T
    ν_fact::T
    cos_Ω::T
    sin_Ω::T
    cos_i::T
    sin_i::T

    # Inner constructor to inforce invariants and pre-calculate a few
    # constants for these elements.
    function KeplerianElements(a, i, e, τ, μ, ω, Ω, plx)


        # Enforce invariants on user parameters
        a = max(a, zero(a))
        e = max(zero(e), min(e, one(e)))
        μ = max(μ, zero(μ))
        plx = max(plx, zero(plx))
        # Pre-calculate some factors that will be re-used when calculating kep2cart at any time
        # Distance in AU
        dist = 1/(plx/1000) * pc2au
        # Compute period (days)
        period = √(a^3/μ) * year2days
        # Mean motion
        n = 2π/√(a^3/μ)
        # Factor in calculating the true anomaly
        ν_fact = √((1+e)/(1-e))

        τ = mod(τ, one(τ))

        T = promote_type(
            typeof(a),
            typeof(i),
            typeof(e),
            typeof(τ),
            typeof(μ),
            typeof(ω),
            typeof(Ω),
            typeof(plx),
        )
        # The user might pass in integers, but it makes no sense to do these
        # calculations on integers. Assume they mean to use floats.
        if T <: Integer
            T = promote_type(T, Float64)
        end

        sin_Ω, cos_Ω = sincos(Ω)
        sin_i, cos_i = sincos(i)
        new{T}(
            # Passed parameters that define the elements
            a,
            i,
            e,
            τ,
            μ,
            ω,
            Ω,
            plx,
            # Cached calcuations
            dist,            
            period,
            n,
            ν_fact,
            # Geometric factors
            cos_Ω,
            sin_Ω,
            cos_i,
            sin_i,
        )
    end
end
# Allow arguments to be specified by keyword.
KeplerianElements(;a, i, e, τ, μ, ω, Ω, plx) = KeplerianElements(a, i, e, τ, μ, ω, Ω, plx)
# And by a named tuple without splatting
KeplerianElements(nt::NamedTuple) = KeplerianElements(nt.a, nt.i, nt.e, nt.τ, nt.μ, nt.ω, nt.Ω, nt.plx)
export KeplerianElements

"""
    astuple(elements)

Return the parameters of a KeplerianElements value as a tuple.
"""
function astuple(elem::KeplerianElements)
    return (;elem.a,elem.i,elem.e,elem.τ,elem.μ,elem.ω,elem.Ω,elem.plx)
end

"""
    KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx)

A convinience function for constructing KeplerianElements where
`i`, `ω`, and `Ω` are provided in units of degrees instead of radians.
"""
KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx) = KeplerianElements(a, deg2rad(i), e, τ, μ, deg2rad(ω), deg2rad(Ω), plx)
KeplerianElementsDeg(;a, i, e, τ, μ, ω, Ω, plx) = KeplerianElementsDeg(a, i, e, τ, μ, ω, Ω, plx)
export KeplerianElementsDeg

function Orbit(args...; kwargs...)
    @warn "Orbit is deprecated in favour of KeplerianElements"
    return KeplerianElements(args...; kwrags...)
end
export Orbit

# Better printing
Base.show(io::IO, ::MIME"text/plain", elem::KeplerianElements) = print(
    io, """
        $(typeof(elem))
        ─────────────────────────
        a   [au ] = $(round(elem.a,sigdigits=3)) 
        i   [°  ] = $(round(rad2deg(elem.i),sigdigits=3))
        e         = $(round(elem.e,sigdigits=3))
        τ         = $(round(elem.τ,sigdigits=3))
        μ   [M⊙ ] = $(round(elem.μ,sigdigits=3)) 
        ω   [°  ] = $(round(rad2deg(elem.ω),sigdigits=3))
        Ω   [°  ] = $(round(rad2deg(elem.Ω),sigdigits=3))
        plx [mas] = $(round(elem.plx,sigdigits=3)) 
        ──────────────────────────
        period      [yrs ] : $(round(period(elem)/year2days,digits=1)) 
        distance    [pc  ] : $(round(distance(elem),digits=1)) 
        mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)),sigdigits=3)) 
        ──────────────────────────
        """)
Base.show(io::IO, elem::KeplerianElements) = print(io,
    "KeplerianElements($(round(elem.a,sigdigits=3)), $(round(elem.i,sigdigits=3)), $(round(elem.e,sigdigits=3)), "*
    "$(round(elem.τ,sigdigits=3)), $(round(elem.μ,sigdigits=3)), $(round(elem.ω,sigdigits=3)), "*
    "$(round(elem.Ω,sigdigits=3)), $(round(elem.plx,sigdigits=3)))"
)

Base.show(io::IO, ::MIME"text/html", elem::KeplerianElements) = print(
    io, """
        <table style="font-family:monospace; text-align: right">
        <tr><th colspan=3 style="font-family:sans-serif; text-align: left">$(typeof(elem))</th></tr>
        <tr><td rowspan=8>Input</td><td>a   [au] =</td> <td>$(round(elem.a,sigdigits=3))</td></tr>
        <tr><td>i   [°] = </td><td>$(round(rad2deg(elem.i),sigdigits=3))</td></tr>
        <tr><td>e         = </td><td>$(round(elem.e,sigdigits=3))</td></tr>
        <tr><td>τ         = </td><td>$(round(elem.τ,sigdigits=3))</td></tr>
        <tr><td>μ   [M⊙] = </td><td>$(round(elem.μ,sigdigits=3)) </td></tr>
        <tr><td>ω   [°] = </td><td>$(round(rad2deg(elem.ω),sigdigits=3))</td></tr>
        <tr><td>Ω   [°] = </td><td>$(round(rad2deg(elem.Ω),sigdigits=3))</td></tr>
        <tr><td>plx [mas] = </td><td>$(round(elem.plx,sigdigits=3)) </td></tr>
        <tr><td rowspan=3>Computed</td><td>period      [yrs] : </td><td>$(round(period(elem)/DirectOrbits.year2days,digits=1)) </td></tr>
        <tr><td>distance    [pc] : </td><td>$(round(distance(elem),digits=1)) </td></tr>
        <tr><td>mean motion [°/yr] : </td><td>$(round(rad2deg(DirectOrbits.meanmotion(elem)),sigdigits=3)) </td></tr>
        </table>
        """)

import Base: length, iterate
length(::AbstractElements) = 1
iterate(elem::AbstractElements) = (elem, nothing)
iterate(::AbstractElements, ::Nothing) = nothing


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

"""
    kep2cart(elements, t)

Given an set of elementsal elements with a time `t` in days to get
a projected displacement x, y, and z in milliarcseconds.
X is increasing to the West, Y increasing to the North, and Z 
away from the observer.

See also: `projectedseparation`, `raoff`, `decoff`, and `losoff`.

In pathalogical cases solving for eccentric anomaly might fail.
This is very unlikely for any reasonable elements with e ≤ 1, but if using
this routine as part of an image distortion step (via e.g. CoordinateTransformations)
than this can occur near the origin. A warning will be generated
and the function will use the mean anomaly in place of the eccentric anomaly.
"""
function kep2cart(elem::KeplerianElements{T}, t; tref=58849) where T
    T2 = promote_type(T, typeof(t))
    

    # Compute mean anomaly
    
    # elem.τ = (tₚ - tᵣ)/period(elem)
    # elem.τ*period(elem) = (tₚ - tᵣ)
    tₚ = elem.τ*period(elem) + tref
    MA = meanmotion(elem)/convert(T2, year2days) * (t - tₚ)

    # τ = (tₚ - tᵣ)/period(elem)  


    MA = rem2pi(MA, RoundDown)

    # EA = eccentric_anomaly(elem.e, MA)
    # EA = eccentric_anomaly_goat(elem.e, MA)
    EA = kepler_solver(MA, elem.e)
    
    # Calculate true anomaly
    ν = convert(T2,2)*atan(elem.ν_fact*tan(EA/convert(T2,2)))

    return kep2cart_ν(elem, ν)
end
export kep2cart

# Kep2cart, only it directly accepts the true anomaly.
# This is used as part of `kep2cart` but also independently
# when drawing orbits.
function kep2cart_ν(elem::KeplerianElements{T}, ν::T) where T
    
    # Semi-latus of rectum    
    p = elem.a*(1-elem.e^2) 
    r = p/(1+elem.e*cos(ν))

    
    # Project back into Cartesian coordinates (AU).
    sin_ω_ν, cos_ω_ν = sincos(elem.ω+ν)
    x_au = r*(elem.sin_Ω*cos_ω_ν + elem.cos_Ω*sin_ω_ν*elem.cos_i)
    y_au = r*(elem.cos_Ω*cos_ω_ν - elem.sin_Ω*sin_ω_ν*elem.cos_i)
    z_au = r*(sin(elem.i)*sin_ω_ν)

    # Radial velocity
    h = sqrt(elem.μ*elem.a*(1-elem.e^2)) # Specific angular momentum
    rv_au = z_au * h * elem.e / (r * p) + h/r * sin(elem.i) * cos_ω_ν

    x_rad = atan(x_au, elem.dist)
    y_rad = atan(y_au, elem.dist)
    z_rad = atan(z_au, elem.dist)

    x_mas = x_rad * rad2as*1e3
    y_mas = y_rad * rad2as*1e3
    z_mas = z_rad * rad2as*1e3

    # rv_kms⁻¹ = rv_au#*au2m*1e-3/year2days
    # Currently au/year?
    # rv_kms⁻¹ = rv_au*au2m*1e-3#/4.904847694504482e6
    rv_kms⁻¹ = rv_au*au2m*1e-3/4.84814e6
    #/year2days/24/60/60

    # coords_AU = SVector(x,y,z)
    # # coords_AU = MVector(x,y,z)
    # # coords_AU = [x,y,z]
    # dist_proj_rad = atan.(coords_AU, elem.dist)
    # dist_proj_mas = dist_proj_rad .* convert(eltype(dist_proj_rad),rad2as*1e3) # radians -> mas

    # return dist_proj_mas
    # return (;x=x_mas, y=y_mas, z=z_mas, rv=rv_kms⁻¹)
    return SVector(x_mas, y_mas, z_mas, rv_kms⁻¹)
end


"""
    eccentric_anomaly(elem, MA)

From an elements and mean anomaly, calculate the eccentric anomaly
numerically (Kepler's equation).

In pathalogical cases solving for eccentric anomaly might fail.
This is unlikely for any reasonable elements, but if using this routine
as part of an image distortion step (via e.g. CoordinateTransformations)
than this can occur near the origin. A warning will be generated
and the function will return (0,0,0). Specifying `throw_ea=true`
turns that warning into an error.
"""
function eccentric_anomaly(e, MA)
    throw_ea = false

    # Numerically solve EA = MA + e * sin(EA) for EA, given MA and e.


    # Fast path for perfectly circular orbits
    if e == 0
        return MA
    end


    # if e ≥ 1
    #     if throw_ea
    #         error("Parabolic and hyperbolic orbits are not yet supported (e≥1)")
    #     else
    #         @warn "Parabolic and hyperbolic orbits are not yet supported (e≥1, e=$e)" maxlog=5
    #         return MA
    #     end
    # end

    # @show e MA


    # Solve for eccentric anomaly
    # The let-block is to capture the bindings of e and M1 directly (performance)
    f= let e = e, MA=MA
        @inline f(EA) = EA - MA - e*sin(EA)
    end

    # After experimentation, Roots finds the root the 
    # fastest / with least allocations using zeroth-order 
    # methods without derivatives. This was surprising.
    # Therefore, we only pass the ojective and no
    # derivatives even though they are trivial to provide.

    # For pathalogical cases, this may not converge.
    # In that case, throw a warning and send the point to the origin


    # For cases very close to one, use a method based on the bisection
    # method immediately
    if isapprox(e, 1, rtol=1e-3)
        try
            # This is a modification of the bisection method. It should be 
            # very very robust.
            EA = find_zero(f, (MA-1, MA+1), FalsePosition(), maxevals=200)
        catch err
            if typeof(err) <: InterruptException
                rethrow(err)
            end
            @warn "Solving for eccentric anomaly near 1 failed. Pass `throw_ea=true` to turn this into an error." e exception=err maxlog=5
            return MA
        end
    end

    # In general, on the other hand:
    local EA
    try
        # Our initial start point EA₀ begins at the mean anomaly.
        # This is a common prescription for fast convergence,
        # though there are more elaborate ways to get better values.
        EA₀ = MA
        # Begin the initial conditions differntly for highly eccentric orbits,
        # another common prescription.
        if e > 0.8
            EA₀ = oftype(MA, π)
        end
        # In benchmarking, the most consistently fast method for solving this root
        # is actually not Newton's method, but the default zeroth order method.
        # We bail out very quickly though if it is not converging (see below)
        EA = find_zero(f, EA₀, maxevals=150)
    catch err
        if typeof(err) <: InterruptException
            rethrow(err)
        end
        # If it fails to converge in some pathalogical case,
        # try a different root finding algorithm.
        # This is a modification of the bisection method. It should be 
        # very very robust.
        # TODO: there are precriptions on how to choose the initial 
        # upper and lower bounds that should be implemented here.
        try
            # EA = find_zero(f, (-2π, 2π), FalsePosition(), maxevals=100)
            EA = find_zero(f, (MA-1, MA+1), FalsePosition(), maxevals=100)
        catch err
            if typeof(err) <: InterruptException
                rethrow(err)
            end
            @warn "Solving for eccentric anomaly failed twice. Pass `throw_ea=true` to turn this into an error." e exception=err maxlog=5
            return MA
        end
    end
end


# Implementation from https://arxiv.org/abs/2103.15829 and https://github.com/oliverphilcox/Keplers-Goat-Herd
function eccentric_anomaly_goat(e, 𝓁)

    # This function implements the 🐐 GOAT algorithm for 
    # solving Kepler's equation. It is approximately
    # 4x faster than the other methods implemented
    # here.

    if isapprox(e, 0)
        return 𝓁
    end

    if isapprox(rem(𝓁,π), 0)
        return 𝓁
    end

    N_it = 15
    N_points = N_it-2
    N_fft = (N_it)*2

    radius = e / 2

    # Generate e^{ikx} sampling points and precompute real and imaginary parts
    # Keep these on the stack inside an MVector -> no allocations
    exp2R = @MVector zeros(typeof(e), N_points)
    exp2I = @MVector zeros(typeof(e), N_points)
    exp4R = @MVector zeros(typeof(e), N_points)
    exp4I = @MVector zeros(typeof(e), N_points)
    coshI = @MVector zeros(typeof(e), N_points)
    sinhI = @MVector zeros(typeof(e), N_points)
    ecosR = @MVector zeros(typeof(e), N_points)
    esinR = @MVector zeros(typeof(e), N_points)
    @inbounds for j in 1:N_points
        freq = 2π*j/N_fft
        cf = cos(freq)
        sf = sin(freq)
        exp2R[j] = cf
        exp2I[j] = sf
        exp4R[j] = cf*cf-sf*sf
        exp4I[j] = 2.0*cf*sf
        coshI[j] = cosh(radius*exp2I[j])
        sinhI[j] = sinh(radius*exp2I[j])
        ecosR[j] = e*cos(radius*exp2R[j])
        esinR[j] = e*sin(radius*exp2R[j])
    end

    esinRadius = e*sin(radius)
    ecosRadius = e*cos(radius)


    # Define contour center for each ell and precompute sin(center), cos(center)
    if 𝓁 < π
        center = 𝓁 + e/2
    else
        center = 𝓁 - e/2
    end
    sinC = sin(center)
    cosC = cos(center)
    output = center

    # Accumulate Fourier coefficients
    # NB: we halve the range by symmetry, absorbing factor of 2 into ratio

    #######
    # Separate out j = 0 piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center + radius

    # Compute e*sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius+cosC*esinRadius # sin(zR)

    # Compute f(z(x)) in real and imaginary parts (fxI = 0)
    fxR = zR - tmpsin - 𝓁

    # Add to array, with factor of 1/2 since an edge
    ft_gx2 = 0.5/fxR
    ft_gx1 = 0.5/fxR

    #######
    # Compute for j = 1 to N_points
    @inbounds @simd for j in 1:N_points

        # Compute z in real and imaginary parts
        zR = center + radius*exp2R[j]
        zI = radius*exp2I[j]

        # Compute f(z(x)) in real and imaginary parts
        # can use precomputed cosh / sinh / cos / sin for this!
        tmpcosh = coshI[j] # cosh(zI)
        tmpsinh = sinhI[j] # sinh(zI)
        tmpsin = sinC*ecosR[j]+cosC*esinR[j] # e sin(zR)
        tmpcos = cosC*ecosR[j]-sinC*esinR[j] # e cos(zR)

        fxR = zR - tmpsin*tmpcosh-𝓁
        fxI = zI - tmpcos*tmpsinh

        # Compute 1/f(z) and append to array
        ftmp = fxR*fxR+fxI*fxI
        fxR /= ftmp
        fxI /= ftmp

        ft_gx2 += (exp4R[j]*fxR+exp4I[j]*fxI)
        ft_gx1 += (exp2R[j]*fxR+exp2I[j]*fxI)
    end

    #######
    # Separate out j = N_it piece, which is simpler

    # Compute z in real and imaginary parts (zI = 0 here)
    zR = center - radius

    # Compute sin(zR) from precomputed quantities
    tmpsin = sinC*ecosRadius-cosC*esinRadius # sin(zR)

    # Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
    fxR = zR - tmpsin-𝓁

    # Add to sum, with 1/2 factor for edges
    ft_gx2 += 0.5/fxR
    ft_gx1 += -0.5/fxR

    #######
    # Compute E(ell)
    output += radius*ft_gx2/ft_gx1;
    
    return output
end


# Using implicit differentiation, the derivatives of eccentric anomaly
# have closed form solutions once the primal value is known. 
# By providing thoesehere, upstream  automatic differentiation libraries
# will be able to efficiently diff through Kepler's equation.
using ChainRulesCore
@scalar_rule eccentric_anomaly(e, MA) @setup(u = 1 - e*cos(Ω)) (sin(Ω) / u, 1 / u)
@scalar_rule eccentric_anomaly_goat(e, MA) @setup(u = 1 - e*cos(Ω)) (sin(Ω) / u, 1 / u)


"""
    raoff(elements, t)

Get the offset from the central body in Right Ascention in
milliarcseconds at some time `t` in days.
"""
function raoff(elements::AbstractElements, t)
    return kep2cart(elements, t)[1]
end
export raoff

"""
    decoff(elements, t)

Get the offset from the central body in Declination in
milliarcseconds at some time `t` in days.
"""
function decoff(elements::AbstractElements, t)
    return kep2cart(elements, t)[2]
end
export decoff

"""
    losoff(elements, t)

Get the offset from the central body in the line of sight towards
the system at time `t` in days, also in milliarcseconds. Of course, we can't observe this
displacement, but we use the same units for consistency.
"""
function losoff(elements::AbstractElements, t)
    return kep2cart(elements, t)[3]
end
export losoff


function radvel(elements::AbstractElements, t)
    return kep2cart(elements, t)[4]
end
export radvel

"""
    projectedseparation(elements, t)

Projected separation in mas from the central body at time t (days).
"""
function projectedseparation(elements::AbstractElements, t)
    x,y,z = kep2cart(elements,t)
    return sqrt(x^2 + y^2 + z^2)
end
export projectedseparation

# TODO: take steps of equal projected distance instead of equal time.

using RecipesBase
@recipe function f(elem::AbstractElements)
    νs = range(-π, π, length=180)
    coords = kep2cart_ν.(elem, νs)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true

    return xs, ys
end

@recipe function f(elems::AbstractArray{<:AbstractElements})

    # Step through true anomaly instead of time.
    # This produces far nicer looking plots, especially if
    # the orbits in question vary significantly in period
    # or are eccentric
    νs = range(-π, π, length=180)
    coords = kep2cart_ν.(elems, νs')

    xs = [c[1] for c in coords]'
    ys = [c[2] for c in coords]'

    # Treat as one long series interrupted by NaN
    xs = reduce(vcat, [[x; NaN] for x in eachcol(xs)])
    ys = reduce(vcat, [[y; NaN] for y in eachcol(ys)])

    # We almost always want to see spatial coordinates with equal step sizes
    aspect_ratio --> 1
    # And we almost always want to reverse the RA coordinate to match how we
    # see it in the sky.
    xflip --> true
    xguide --> "ΔRA - mas"
    yguide --> "ΔDEC - mas"

    seriesalpha --> 30/length(elems)


    return xs, ys
end

include("Fitting.jl")
include("Transformation.jl")

end # module
