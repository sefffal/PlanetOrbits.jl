"""
    ThieleInnesOrbit(e, τ, M, plx, A, B, F, G)

Represents a visual orbit of a planet using Thiele-Innes
orbital elements. Convertable to and from a VisualOrbit.
This parameterization does not have the issue that traditional
angular parameters have where the argument of periapsis and
longitude of ascending node become undefined for circular and face
on orbits respectively.
"""
struct ThieleInnesOrbit{T<:Number} <: AbstractOrbit

    # Orbital properties
    e::T
    τ::T
    M::T
    plx::T
    A::T
    B::T
    F::T
    G::T

    # Constants
    C::T
    H::T
    T::T
    n::T
    ν_fact::T

    # Inner constructor to enforce invariants and pre-calculate
    # constants from the orbital elements
    function ThieleInnesOrbit(e, τ, M, plx, A, B, F, G)
        e, τ, M, plx, A, B, F, G = promote(e, τ, M, plx, A, B, F, G)
        T = typeof(e)

        # TODO: confirm these following lines are necessary for us to calculate
        # the period and mean motion (that's all we really need)
        u = (A^2 + B^2 + F^2 + G^2)/2
        v = A*G - B * F
        α = sqrt(u + sqrt((u+v)*(u-v)))
        a = α/plx

        # Calculate coefficients C & H.
        # TODO: Is this necessary?
        ω_p_Ω = atan((B-F),(A+G))
        ω_m_Ω = atan((B+F),(G-A))
        ω = (ω_p_Ω+ω_m_Ω)/2+π/2 # There is a convention difference we account for with this phase shift. We want the ω of the planet not the primary.
        Ω = (ω_p_Ω-ω_m_Ω)/2
        if Ω < 0
            ω += π
            Ω += π
        end
        s,c = sincos(ω-Ω)
        d₁ = abs((A+G)*c)
        d₂ = abs((F-B)*s)
        s2,c2 = sincos(ω+Ω)
        if d₁ >= d₂
            i = 2atan(sqrt(abs((A-G)*c2)/d₁))
        else
            i = 2atan(sqrt(abs((B+F)*s2)/d₂))
        end
        C = a*sin(ω)*sin(i)
        H = a*cos(ω)*sin(i)

        periodyrs = √(a^3/M)
        period = periodyrs * year2day # period [days]
        n = 2π/√(a^3/M) # mean motion

        ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor

        new{T}(e, τ, M, plx, A, B, F, G, C, H, period, n, ν_fact)
    end
end
ThieleInnesOrbit(;e, τ, M, plx, A, B, F, G) = ThieleInnesOrbit(e, τ, M, plx, A, B, F, G)
ThieleInnesOrbit(nt) = ThieleInnesOrbit(nt.e, nt.τ,nt. M, nt.plx, nt.A, nt.B, nt.F, nt.G)

export ThieleInnesOrbit

"""
Represents a `ThieleInnesOrbit` evaluated to some position.
"""
struct OrbitSolutionThieleInnes{T<:Number,TEl<:ThieleInnesOrbit} <: AbstractOrbitSolution
    elem::TEl
    ν::T
    EA::T
    x::T
    y::T
    ẋ::T
    ẏ::T
    t::T
    
    function OrbitSolutionThieleInnes(elem, ν, EA, x, y,ẋ, ẏ,t)
        promoted = promote(ν, EA, x, y,ẋ,ẏ, t)
        return new{eltype(promoted),typeof(elem)}(elem, promoted...)
    end
end
export OrbitSolutionThieleInnes

function orbitsolve_ν(elem::ThieleInnesOrbit, ν, EA=2atan(tan(ν/2)/elem.ν_fact), t=_time_from_EA(elem, EA))
    # https://arxiv.org/ftp/arxiv/papers/1008/1008.3416.pdf
    sea, cea = sincos(EA)
    x = cea - elem.e
    y = sea * sqrt(1 - elem.e^2)
    μ = elem.n
    ẋ = -μ*sin(EA)/(1-elem.e*cos(EA))
    ẏ = √(1-elem.e^2)*μ*cos(EA)/(1-elem.e*cos(EA))
    return OrbitSolutionThieleInnes(elem, ν, EA, x, y, ẋ, ẏ, t)
end


function raoff(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.B + sol.y*sol.elem.G
end

function decoff(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.A + sol.y*sol.elem.F
end

# Radial velocity not currently right. Z position is correct.
function radvel(sol::OrbitSolutionThieleInnes)
    (sol.ẋ*sol.elem.C + sol.ẏ*sol.elem.H)*au2m*sec2year
end


function posx(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.B/sol.elem.plx + sol.y*sol.elem.G/sol.elem.plx
end
function posy(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.A/sol.elem.plx + sol.y*sol.elem.F/sol.elem.plx
end
function posz(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.C + sol.y*sol.elem.H
end

function pmra(sol::OrbitSolutionThieleInnes)
    sol.ẋ*sol.elem.B + sol.ẏ*sol.elem.G
end

function pmdec(sol::OrbitSolutionThieleInnes)
    sol.ẋ*sol.elem.A + sol.ẏ*sol.elem.F
end


function ThieleInnesOrbit(orbit::VisualOrbit)
    α = orbit.a*orbit.plx 

    A = α*( orbit.cosΩ*cos(orbit.ω)-orbit.sinΩ*sin(orbit.ω)*orbit.cosi)
    B = α*( orbit.sinΩ*cos(orbit.ω)+orbit.cosΩ*sin(orbit.ω)*orbit.cosi)
    F = α*(-orbit.cosΩ*sin(orbit.ω)-orbit.sinΩ*cos(orbit.ω)*orbit.cosi)
    G = α*(-orbit.sinΩ*sin(orbit.ω)+orbit.cosΩ*cos(orbit.ω)*orbit.cosi)

    ThieleInnesOrbit(orbit.e, orbit.τ, orbit.M, orbit.plx, A, B, F, G)
end



function VisualOrbit(orbit::ThieleInnesOrbit)

    # TODO: There is something incorrect here in this conversion to do with omega

    u = (orbit.A^2 + orbit.B^2 + orbit.F^2 + orbit.G^2)/2
    v = orbit.A*orbit.G - orbit.B * orbit.F
    α = sqrt(u + sqrt((u+v)*(u-v)))
    ω_p_Ω = atan((orbit.B-orbit.F),(orbit.A+orbit.G))
    ω_m_Ω = atan((orbit.B+orbit.F),(orbit.G-orbit.A))
    ω = (ω_p_Ω+ω_m_Ω)/2+π/2 # There is a convention difference we account for with this phase shift. We want the ω of the planet not the primary.
    Ω = (ω_p_Ω-ω_m_Ω)/2
    if Ω < 0
        ω += π
        Ω += π
    end
    s,c = sincos(ω-Ω)
    d₁ = abs((orbit.A+orbit.G)*c)
    d₂ = abs((orbit.F-orbit.B)*s)
    s2,c2 = sincos(ω+Ω)
    if d₁ >= d₂
        i = 2atan(sqrt(abs((orbit.A-orbit.G)*c2)/d₁))
    else
        i = 2atan(sqrt(abs((orbit.B+orbit.F)*s2)/d₂))
    end
    a = α/orbit.plx

    return VisualOrbit(a, orbit.e, i, ω, Ω, orbit.τ, orbit.M, orbit.plx)
end