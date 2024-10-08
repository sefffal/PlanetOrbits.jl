"""
    ThieleInnesOrbit(e, tp, M, plx, A, B, F, G)

Represents a visual orbit of a planet using Thiele-Innes
orbital elements. Convertable to and from a VisualOrbit.
This parameterization does not have the issue that traditional
angular parameters have where the argument of periapsis and
longitude of ascending node become undefined for circular and face
on orbits respectively.

!!! warning
    There is a remaining bug in this implementation for pi <= Ω < 2pi
"""
struct ThieleInnesOrbit{T<:Number} <: AbstractOrbit{T}

    # Orbital properties
    e::T
    tp::T
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
    function ThieleInnesOrbit(e, tp, M, plx, A, B, F, G)
        e, tp, M, plx, A, B, F, G = promote(e, tp, M, plx, A, B, F, G)
        T = typeof(e)

        # TODO: confirm these following lines are necessary for us to calculate
        # the period and mean motion (that's all we really need)
        u = (A^2 + B^2 + F^2 + G^2)/2
        v = A*G - B * F
        α = sqrt(u + sqrt((u+v)*(u-v)))
        a = α/plx
        if e > 1
            a = -a
            @warn "Support for hyperbolic Thiele-Innes orbits is not yet complete. Results will be silently wrong."
        end

        ω_p_Ω = atan((B-F),(A+G))
        ω_m_Ω = atan((B+F),(G-A))
        
        
        # sign of ω_p_Ω: sin(ω_p_Ω) same sign as (B-F)
        # sign of ω_m_Ω: sin(ω_m_Ω) same sign as (-B-F)
        if sign(sin(ω_p_Ω)) != sign(B-F)
            ω_p_Ω = ω_p_Ω + pi
        end; 
        if sign(sin(ω_m_Ω)) != sign(-B-F)
            ω_m_Ω = ω_m_Ω + pi
        end;
        ω_p_Ω = rem2pi(ω_p_Ω, RoundDown)
        ω_m_Ω = rem2pi(ω_m_Ω, RoundDown)
        ω = (ω_p_Ω + ω_m_Ω)/2
        Ω = (ω_p_Ω - ω_m_Ω)/2
        # TODO: we have an error with the following orbit orbit(M=1, e=0.4, a=1, i=2, Ω=4, ω=1, tp=0, plx=10.0);
        # Error seems to happen when Ω + ω > 2pi
        # The issue is that ω ends up off by pi.
        if Ω < 0 
            ω += π
            Ω += π
        end
        ω = rem2pi(ω, RoundDown)
        
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

        # Pre-calculate factors to be re-used by orbitsolve
        # Physical constants of system and orbit
        if e < 1
            period_days = √(a^3/M)*kepler_year_to_julian_day_conversion_factor
            period_yrs = period_days/year2day_julian
            n = 2π/period_yrs # mean motio
        else
            period = Inf
            # TODO: Need to confirm where this 2pi is coming from 
            n = 2π * √(M/-a^3)*kepler_year_to_julian_day_conversion_factor/year2day_julian
            # n = √(M/-a^3) # mean motion
        end

        if e < 1
            ν_fact = √((1 + e)/(1 - e)) # true anomaly prefactor
        else
            ν_fact = √((1 + e)/(e - 1)) # true anomaly prefactor
        end


        new{T}(e, tp, M, plx, A, B, F, G, C, H, period_days, n, ν_fact)
    end
end
ThieleInnesOrbit(;e, tp, M, plx, A, B, F, G, kwargs...) = ThieleInnesOrbit(e, tp, M, plx, A, B, F, G)

export ThieleInnesOrbit

period(elem::ThieleInnesOrbit) = elem.T
meanmotion(elem::ThieleInnesOrbit) = elem.n
eccentricity(o::ThieleInnesOrbit) = o.e
totalmass(o::ThieleInnesOrbit) = o.M
function semimajoraxis(o::ThieleInnesOrbit)
    (;A,B,F,G,plx) = o
    u = (A^2 + B^2 + F^2 + G^2)/2
    v = A*G - B * F
    α = sqrt(u + sqrt((u+v)*(u-v)))
    a = α/plx
    return a
end
function inclination(o::ThieleInnesOrbit)
    # TODO: test

    ω_p_Ω = atan((o.B-o.F),(o.A+o.G))
    ω_m_Ω = atan((o.B+o.F),(o.G-o.A))
    
    
    # sign of ω_p_Ω: sin(ω_p_Ω) same sign as (B-F)
    # sign of ω_m_Ω: sin(ω_m_Ω) same sign as (-B-F)
    if sign(sin(ω_p_Ω)) != sign(o.B-o.F)
        ω_p_Ω = ω_p_Ω + pi
    end; 
    if sign(sin(ω_m_Ω)) != sign(-o.B-o.F)
        ω_m_Ω = ω_m_Ω + pi
    end;
    ω_p_Ω = rem2pi(ω_p_Ω, RoundDown)
    ω_m_Ω = rem2pi(ω_m_Ω, RoundDown)
    ω = (ω_p_Ω + ω_m_Ω)/2
    Ω = (ω_p_Ω - ω_m_Ω)/2
    Ω, ω
    if Ω < 0
        ω += π
        Ω += π
    end
    
    s,c = sincos(ω-Ω)
    d₁ = abs((o.A+o.G)*c)
    d₂ = abs((o.F-o.B)*s)
    s2,c2 = sincos(ω+Ω)
    if d₁ >= d₂
        i = 2atan(sqrt(abs((o.A-o.G)*c2)/d₁))
    else
        i = 2atan(sqrt(abs((o.B+o.F)*s2)/d₂))
    end
    return i
end
function _trueanom_from_eccanom(o::ThieleInnesOrbit, EA)
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
periastron(o::ThieleInnesOrbit) = o.tp
function semiamplitude(o::ThieleInnesOrbit)
    # TODO: test implementation
    oneminusesq = (1 - eccentricity(o)^2)
    a = semimajoraxis(o)
    sini = sin(inclination(o))
    J = ((2π*a)/period(o)*day2year_julian) / √oneminusesq # horizontal velocity semiamplitude [AU/year]
    K = J*au2m*sec2year_julian*sini # radial velocity semiamplitude [m/s]
    return K
end
distance(o::ThieleInnesOrbit) = 1000/o.plx * pc2au

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

function orbitsolve_ν(elem::ThieleInnesOrbit, ν, EA=2atan(tan(ν/2)/elem.ν_fact), t=_time_from_EA(elem, EA))
    # https://arxiv.org/ftp/arxiv/papers/1008/1008.3416.pdf
    sea, cea = sincos(EA)
    x = cea - elem.e
    if elem.e < 1
        y = sea * sqrt(1 - elem.e^2)
        ẏ = √(1-elem.e^2)*elem.n*cos(EA)/(1-elem.e*cos(EA))
    else
        # TODO: this is just a guess
        y = sea * sqrt(elem.e^2 - 1)
        ẏ = √(elem.e^2-1)*elem.n*cos(EA)/(elem.e*cos(EA)-1)
    end
    ẋ = -elem.n*sin(EA)/(1-elem.e*cos(EA))
    return OrbitSolutionThieleInnes(elem, ν, EA, x, y, ẋ, ẏ, t)
end
soltime(os::OrbitSolutionThieleInnes) = os.t


function raoff(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.B + sol.y*sol.elem.G
end

function decoff(sol::OrbitSolutionThieleInnes)
    sol.x*sol.elem.A + sol.y*sol.elem.F
end

# Radial velocity not currently right. Z position is correct.
function radvel(sol::OrbitSolutionThieleInnes)
    (sol.ẋ*sol.elem.C + sol.ẏ*sol.elem.H)*au2m*sec2year_julian
end


function posx(sol::OrbitSolutionThieleInnes)
    raoff(sol)/sol.elem.plx
end
function posy(sol::OrbitSolutionThieleInnes)
    decoff(sol)/sol.elem.plx
end
function posz(sol::OrbitSolutionThieleInnes)
    (sol.x*sol.elem.C + sol.y*sol.elem.H)
end

function pmra(sol::OrbitSolutionThieleInnes)
    sol.ẋ*sol.elem.B + sol.ẏ*sol.elem.G
end

function pmdec(sol::OrbitSolutionThieleInnes)
    sol.ẋ*sol.elem.A + sol.ẏ*sol.elem.F
end


function ThieleInnesOrbit(orbit::Visual{KepOrbit{T1},T2}) where {T1,T2}
    a = semimajoraxis(orbit)
    α = a*orbit.plx 
    elem = orbit.parent
    A = α*( elem.cosΩ*cos(elem.ω)-elem.sinΩ*sin(elem.ω)*elem.cosi)
    B = α*( elem.sinΩ*cos(elem.ω)+elem.cosΩ*sin(elem.ω)*elem.cosi)
    F = α*(-elem.cosΩ*sin(elem.ω)-elem.sinΩ*cos(elem.ω)*elem.cosi)
    G = α*(-elem.sinΩ*sin(elem.ω)+elem.cosΩ*cos(elem.ω)*elem.cosi)

    ThieleInnesOrbit(
        eccentricity(orbit),
        periastron(orbit),
        totalmass(orbit),
        1000/distance(orbit),
        A, B, F, G
    )
end



function Visual{KepOrbit}(o::ThieleInnesOrbit)
    u = (o.A^2 + o.B^2 + o.F^2 + o.G^2)/2
    v = o.A*o.G - o.B * o.F
    α = sqrt(u + sqrt((u+v)*(u-v)))
    
   
    ω_p_Ω = atan((o.B-o.F),(o.A+o.G))
    ω_m_Ω = atan((o.B+o.F),(o.G-o.A))
    
    
    # sign of ω_p_Ω: sin(ω_p_Ω) same sign as (B-F)
    # sign of ω_m_Ω: sin(ω_m_Ω) same sign as (-B-F)
    if sign(sin(ω_p_Ω)) != sign(o.B-o.F)
        ω_p_Ω = ω_p_Ω + pi
    end; 
    if sign(sin(ω_m_Ω)) != sign(-o.B-o.F)
        ω_m_Ω = ω_m_Ω + pi
    end;
    ω_p_Ω = rem2pi(ω_p_Ω, RoundDown)
    ω_m_Ω = rem2pi(ω_m_Ω, RoundDown)
    ω = (ω_p_Ω + ω_m_Ω)/2
    Ω = (ω_p_Ω - ω_m_Ω)/2
    if Ω < 0
        ω += π
        Ω += π
    end
    
    
    s,c = sincos(ω-Ω)
    d₁ = abs((o.A+o.G)*c)
    d₂ = abs((o.F-o.B)*s)
    s2,c2 = sincos(ω+Ω)
    if d₁ >= d₂
        i = 2atan(sqrt(abs((o.A-o.G)*c2)/d₁))
    else
        i = 2atan(sqrt(abs((o.B+o.F)*s2)/d₂))
    end

    a = α/o.plx
    return Visual{KepOrbit}(;a, o.e, i, ω, Ω, o.tp, o.M, o.plx)
end


# Pretty printing
Base.show(io::IO, ::MIME"text/plain", elem::ThieleInnesOrbit) = print(
io, """
    $(typeof(elem))
    ─────────────────────────
    A   [mas] = $(round(elem.A, sigdigits=3))
    B   [mas] = $(round(elem.B, sigdigits=3))
    F   [mas] = $(round(elem.F, sigdigits=3))
    G   [mas] = $(round(elem.G, sigdigits=3))
    e         = $(round(elem.e, sigdigits=8))
    tp         = $(round(elem.tp, sigdigits=3))
    M   [M⊙ ] = $(round(elem.M, sigdigits=3)) 
    period      [yrs ] : $(round(period(elem)*day2year_julian, digits=1)) 
    mean motion [°/yr] : $(round(rad2deg(meanmotion(elem)), sigdigits=3)) 
    ──────────────────────────
    """
)
