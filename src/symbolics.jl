
Symbolics.@register kepler_solver(EA, e)

# Derivative w.r.t. the first argument
function Symbolics.derivative(::typeof(kepler_solver), args::NTuple{2,Any}, ::Val{1})
    MA, e = args
    EA = kepler_solver(MA, e)
    u = 1 - e*cos(EA)
    return 1 / u            
end
# Derivative w.r.t. the second argument
function Symbolics.derivative(::typeof(kepler_solver), args::NTuple{2,Any}, ::Val{2})
    MA, e = args
    EA = kepler_solver(MA, e)
    u = 1 - e*cos(EA)
    return sin(EA) / u    
end

Base.show(io::IO, ::MIME"text/plain", elem::KeplerianElements) = print(
    io, """
        $(typeof(elem))
        ─────────────────────────
        a   [au ] = $(elem.a) 
        i   [rad] = $(elem.i)
        e         = $(elem.e)
        τ         = $(elem.τ)
        M   [M⊙ ] = $(elem.M) 
        ω   [rad] = $(elem.ω)
        Ω   [rad] = $(elem.Ω)
        plx [mas] = $(elem.plx) 
        ──────────────────────────
        period      [yrs ]   : $(period(elem)) 
        distance    [pc  ]   : $(distance(elem)) 
        mean motion [rad/yr] : $(rad2deg(meanmotion(elem))) 
        ──────────────────────────
        """)