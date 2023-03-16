"""
    Visual{OrbitType}(..., plx=...)

This wraps another orbit to add the parallax distance field `plx`,
thus allowing projected quantities to be calculated.
It forwards everything else to the parent orbit.

For example, the KepOrbit type supports calculating x and y positions in AU.
A Visual{KepOrbit} additionally supports calculating projected
right ascension and declination offsets.

Note: the `ThieleInnesOrbit` type does not need to be wrapped in `Visual`
as it the Thiele-Innes constants are already expressed in milliarcseconds and
thus it always requires a `plx` value.
"""
struct VisualOrbit{T<:Number,O<:AbstractOrbit} <: AbstractOrbit{T}
    parent::O
    plx::T
    dist::T
end
distance(elem::VisualOrbit) = elem.dist*au2pc

"""
    Visual{OrbitType}(..., plx=...)

This wraps another orbit to add the parallax distance field `plx`,
thus allowing projected quantities to be calculated.
It forwards everything else to the parent orbit.

For example, the KepOrbit type supports calculating x and y positions in AU.
A Visual{KepOrbit} additionally supports calculating projected
right ascension and declination offsets.

!!! note
    The `ThieleInnesOrbit` type does not need to be wrapped in `Visual`
    as it the Thiele-Innes constants are already expressed in milliarcseconds and
    thus it always requires a `plx` value.
"""
const Visual{OrbitType} = VisualOrbit{T,OrbitType}  where T

_parent_num_type(orbit::AbstractOrbit{T}) where T = T
function Visual{OrbitType}(;plx, args...,) where {OrbitType}
    dist = 1000/plx * pc2au # distance [AU]
    parent = OrbitType(;args...)
    T = _parent_num_type(parent)
    return VisualOrbit{T,OrbitType}(parent, plx, dist)
end
function Visual(parent::AbstractOrbit, plx,)
    dist = 1000/plx * pc2au # distance [AU]
    T = _parent_num_type(parent)
    return VisualOrbit{T,typeof(parent)}(parent, plx, dist)
end

export Visual

struct OrbitSolutionVisual{TEl<:AbstractOrbit,TSol<:AbstractOrbitSolution} <: AbstractOrbitSolution
    elem::TEl
    sol::TSol
end

function orbitsolve_ν(elem::VisualOrbit, ν, args...; kwargs...)
    sol = orbitsolve_ν(elem.parent, ν, args...; kwargs...)
    return OrbitSolutionVisual(elem, sol)
end
soltime(os::OrbitSolutionVisual) = soltime(os.sol)


# Forward these functions to the underlying orbit object
solution_fun_list = (
    :trueanom,
    :eccanom,
    :meananom,
    :posx,
    :posy,
    :posz,
    :posangle,
    :velx,
    :vely,
    :velz,
    :radvel,
)
for fun in solution_fun_list
    @eval function ($fun)(os::OrbitSolutionVisual, args...)
        return ($fun)(os.sol, args...)
    end
end
orbit_fun_list = (
    :eccentricity,
    :periastron,
    :period,
    :hostmass,
    :meanmotion,
    :semiamplitude,
    :_trueanom_from_eccanom,
)
for fun in orbit_fun_list
    @eval function ($fun)(elem::VisualOrbit, args...)
        return ($fun)(elem.parent, args...)
    end
end

function raoff(o::OrbitSolutionVisual)
    xcart = posx(o) # [AU]
    cart2angle = rad2as*oftype(xcart, 1e3)/o.elem.dist
    xang = xcart*cart2angle # [mas]
    return xang
end

function decoff(o::OrbitSolutionVisual)
    ycart = posy(o) # [AU]
    cart2angle = rad2as*oftype(ycart, 1e3)/o.elem.dist
    yang = ycart*cart2angle # [mas]
    return yang
end

function pmra(o::OrbitSolutionVisual)
    ẋcart = o.elem.parent.J*(o.elem.parent.cosi_cosΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) - o.elem.parent.sinΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẋcart, 1e3)/o.elem.dist
    ẋang = ẋcart*cart2angle # [mas/year]
    return ẋang
end

function pmdec(o::OrbitSolutionVisual)
    ẏcart = -o.elem.parent.J*(o.elem.parent.cosi_sinΩ*(o.sol.cosν_ω + o.elem.parent.ecosω) + o.elem.parent.cosΩ*(o.sol.sinν_ω + o.elem.parent.esinω)) # [AU/year]
    cart2angle = rad2as*oftype(ẏcart, 1e3)/o.elem.dist
    ẏang = ẏcart*cart2angle # [mas/year]
    return ẏang
end

function accra(o::OrbitSolutionVisual{<:Visual{KepOrbit}})
    ẍcart = -o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_cosΩ*o.sol.sinν_ω + o.elem.parent.sinΩ*o.sol.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ẍcart, 1e3)/o.elem.dist
    ẍang = ẍcart*cart2angle # [mas/year^2] 
    return ẍang
end
function accdec(o::OrbitSolutionVisual{<:Visual{KepOrbit}})
    ÿcart = o.elem.parent.A*(1 + o.sol.ecosν)^2 * (o.elem.parent.cosi_sinΩ*o.sol.sinν_ω - o.elem.parent.cosΩ*o.sol.cosν_ω) # [AU/year^2]
    cart2angle = rad2as*oftype(ÿcart, 1e3)/o.elem.dist
    ÿang = ÿcart*cart2angle # [mas/year^2] 
    return ÿang
end



# Pretty printing
function Base.show(io::IO, mime::MIME"text/plain", elem::Visual)
    show(io, mime, elem.parent)
    print(io, """\
    plx [mas] = $(round(elem.plx, sigdigits=3)) 
    distance    [pc  ] : $(round(distance(elem), digits=1)) 
    ──────────────────────────
    """)
end
