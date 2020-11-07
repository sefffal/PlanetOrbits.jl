"""
Module for calculating orbits and converting into cartiesian coordinates.

# Available Functions
- Orbit(...)
- τ2t0(τ,P)
- period(orbit)
- xyz(orbit, t)
- plot(orbit)
"""
module DirectOrbits

using DirectImages


# Required packages
using LinearAlgebra
using CoordinateTransformations
using StaticArrays
using Roots # For solving for eccentric anomaly
import Dates
import Base.inv
using Statistics: mean
using RecipesBase


"""
Convert from (fraction of orbit past periastron past MJD=0)
to the time of periastron passage in MJD.
P in years.
τ [0,1]
"""
function τ2t0(τ,P)
    # if tau == 0, then T0 = 0 MJD
    P*365.25 - τ*(P*365.25)# - 10.0*365#9*365.
end

## Function that takes orbital parameters and time, and maps them to an image location
struct Orbit{T<:Real}
    a::T
    i::T
    e::T
    τ::T
    μ::T
    ω::T
    Ω::T
    plx::T
end
Orbit(;a, i, e, τ, μ, ω, Ω, plx) = Orbit(a, i, e, τ, μ, ω, Ω, plx)
export Orbit

import Base: length, iterate
length(::Orbit) = 1
iterate(orb::Orbit) = (orb, nothing)
iterate(::Orbit, ::Nothing) = nothing

"""
Period of orbit in years.
"""
period(orb::Orbit) = sqrt(orb.a^3/orb.μ)

"""
    xyz(orbit, t)
Call an Orbit object with a time (as a modified Julian date) to get
a projected position x,y,and z in milliarcseconds.
"""
function xyz(orb::Orbit, t)
    mas2rad = 4.8481368E-9
    rad2mas = 206264806.2471
    pc2au = 206265e3
    au2m = 1.495978707e11
    
    dist = 1/(orb.plx*1e-3) * pc2au #* au2m

    # Compute period (years)
    T = sqrt(orb.a^3/orb.μ)
    n = 2π/T
    # Compute mean anomaly
    T0 = τ2t0(orb.τ, T)
    MA = n*(t-T0)/365.25
    MA = rem2pi(MA, RoundUp)
    # @show dist/pc2au orb.a T
    
    # Solve for eccentric anomaly
    # The let-block is to capture the bindings of e and M1 directly (performance)
    f, fp = let e = orb.e, MA=MA
        @inline f(EA) = EA - MA - e*sin(EA)
        @inline fp(EA) = 1 - cos(EA)
        f,fp
    end
    # For pathalogical cases, this may not converge.
    # In that case, throw a warning and send the point to the origin (not much else we can do)
    local EA
    try
        # Newton's method
        # EA = Roots.newton(f, fp, π, MA)
        # Bisection method.
        # This seems both slight faster and more robust.
        # EA = find_zero(f, (-π, π))
        EA = find_zero(f, MA)
    catch err
        if typeof(err) <: InterruptException
            rethrow(err)
        end
        @warn "Solving for eccentric anomaly failed" exception=err maxlog=5
        return SVector(0.,0.,0.0)
    end
    # @show EA
    # EA = MA
    # Calculate true anomaly
    ν = 2atan(√((1+orb.e)/(1-orb.e))*tan(EA/2))
    # New orbital radius.
    # This is the semi-major axis, modified by the eccentricity
    r = orb.a*(1-orb.e*cos(EA))
    # h = sqrt(orb.μ*orb.a*(1-orb.e^2))
    # x = r*(cos(orb.Ω)*cos(orb.ω + ν) - sin(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    # y = r*(sin(orb.Ω)*cos(orb.ω + ν) - cos(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    # # Project back into Cartesian coordinates (meters).
    x = r*(cos(orb.Ω)*cos(orb.ω+ν) - sin(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    y = r*(sin(orb.Ω)*cos(orb.ω+ν) + cos(orb.Ω)*sin(orb.ω+ν)*cos(orb.i))
    # Now convert back to projected separation in pixels
    x = atan(x,dist)*rad2mas # radians -> mas
    y = atan(y,dist)*rad2mas # radians -> mas
    z = atan(r*(sin(orb.i)*sin(orb.ω+ν)),dist)*rad2mas # radians -> mas
    return SVector(y,x,z)
end



@recipe function f(::Type{T}, orb::T) where{T<:Orbit}
    ts = range(0, 365.25period(orb), length=200)
    coords = xyz.(orb, ts)
    xs = [c[1] for c in coords].*1e3
    ys = [c[2] for c in coords].*1e3

    xguide --> "x - mas"
    yguide --> "y - mas"
    framestyle --> :box
    fontfamily --> "Times"
    minorticks --> true
    # return (xs, ys)
    # return (xs..., ys...)
    return xs, ys
end

export orbitplot

using Requires
@require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
    function orbitplot(orb::Orbit, args...; kwargs...)
        ts = range(0, 365.25period(orb), length=200)
        coords = xyz.(orb, ts)
        xs = [c[1] for c in coords].*1e3
        ys = [c[2] for c in coords].*1e3
        # zs = [c[3] for c in coords]
        return Plots.plot(
            xs,
            ys,
            args...;
            aspectratio=1,
            xlabel="x - mas",
            ylabel="y - mas",
            framestyle=:box,
            fontfamily = "Times",
            minorticks = true,
            kwargs...
        )
    end
end


end # module
