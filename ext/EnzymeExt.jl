"""
This extension package adds rule support for Enzyme.
There are a few useful rules we can define for the sake of performance
but the most important thing is a rule for the Kepler solver.

"""
module EnzymeExt

using PlanetOrbits
using Enzyme: EnzymeRules

# @scalar_rule kepler_solver(M, e) @setup(u = 1 - e*cos(Ω)) (1 / u,sin(Ω) / u)

# These could be added to improve performance. The units would have to be adapted.
# @scalar_rule raoff(o::OrbitSolution) pmra(o)
# @scalar_rule decoff(o::OrbitSolution) pmdec(o)
# @scalar_rule propmotionanom(o::OrbitSolution) acceleration(o)
# @scalar_rule orbitsolve_ν(elem::VisualOrbit{T}, t; tref=58849)

# EnzymeRules.reverse(::EnzymeRules.Config, func::EnzymeRules.Annotation{typeof(PlanetOrbits.kepler_solver)}, dret::Active, tape, args::EnzymeRules.Annotation...) = 1

# function EnzymeRules.reverse(
#     # config::EnzymeRules.ConfigWidth{1},
#     config,
#     # func::typeof(PlanetOrbits.kepler_solver)
#     func,
#     args...
#     # ::Type{<:Active},
#     # ::Any,
#     # ma::Active,
#     # e::Active
# )
#     @show config func ma e
#     if needs_primal(config)
#         return AugmentedReturn(func.val(x.val),nothing, nothing)
#     else
#         return AugmentedReturn(nothing, nothing, nothing)
#     end
# end

# function EnzymeRules.reverse(
#     conf::Config,
#     func::Annotation{typeof(f)},
#     ann::Type{<:Annotation},
#     tape,
#     args::Annotation...
# )
#     @show conf func ann tape args
# end

# function forward(::Const{typeof(PlanetOrbits.kepler_solver)}, ::Type{<:DuplicatedNoNeed}, x::Duplicated)
#     return 10+2*x.val*x.dval
# end

# function forward(::Const{typeof(PlanetOrbits.kepler_solver)}, ::Type{<:BatchDuplicatedNoNeed}, x::BatchDuplicated{T, N}) where {T, N}
#     return NTuple{N, T}(1000+2*x.val*dv for dv in x.dval)
# end

# function forward(func::Const{typeof(PlanetOrbits.kepler_solver)}, ::Type{<:Duplicated}, x::Duplicated)
#     return Duplicated(func.val(x.val), 100+2*x.val*x.dval)
# end

# function forward(func::Const{typeof(PlanetOrbits.kepler_solver)}, ::Type{<:BatchDuplicated}, x::BatchDuplicated{T, N}) where {T,N}
#     return BatchDuplicated(func.val(x.val), NTuple{N, T}(10000+2*x.val*dv for dv in x.dval))
# end

# function forward(::Const{Core.typeof(f_ip)}, ::Type{<:Const}, x::Duplicated)
#     ld = x.val[1]
#     x.val[1] *= ld
#     x.dval[1] *= 2 * ld + 10
#     return nothing
# end

end
