# ---------------------------------------------------
# Analytic Kepler-solve derivatives for Enzyme.
#
# The same implicit-function rule ForwardDiff gets from the `Dual` methods in
# `kepsolve.jl` and ChainRules consumers get from the `@scalar_rule`: with E
# defined by E − e sin E = M,
#
#     dE = (dM + sin(E) de) / (1 − e cos E).
#
# Without this, Enzyme differentiates *through* Markley's closed-form
# approximation sequence. That is not merely slower — the derivative of an
# approximation is a worse approximation than the value it came from, and the
# `rem2pi` at the top of the solver is an additional hazard.
#
# Elliptical branch only, matching `kepsolve.jl`: hyperbolic solves fall
# through to Enzyme's generic handling of the iterative path.
#
# NOTE (2026-08-01): Enzyme cannot yet differentiate a whole `orbitsolve!`
# — `solve_row_simd!` raises `EnzymeNoDerivativeError` and AHL21 segfaults —
# so this rule is exercised only at the kernel level today. See §10.4.2 of
# design/planetorbits-v2-nbody-migration.md for the measurements behind
# staying ForwardDiff-first.
# ---------------------------------------------------
module PlanetOrbitsEnzymeCoreExt

using PlanetOrbits: PlanetOrbits, kepler_solver, Markley, Goat, RootsMethod, vcbrt
using EnzymeCore
using EnzymeCore.EnzymeRules

const EllipticSolver = Union{Markley,Goat,RootsMethod}

# ---------------------------------------------------
# `vcbrt`: the one thing in the SIMD Kepler path Enzyme cannot handle.
#
# Its seed is an exponent bit-hack (`reinterpret(Float64, u ÷ 3 + magic)`),
# which lowers to an LLVM `bitcast` and raises `EnzymeNoDerivativeError` —
# blocking `solve_row_simd!`, and with it the whole SIMD value path, which is
# ~2× faster than the scalar fallback. Everything else in the kernel
# (`vrem2pi`, `vsincos`, `vsincos_small`) differentiates fine; measured.
#
# The seed's derivative is irrelevant: four Newton steps follow, so the value
# converges to cbrt(x) whatever the seed, and the derivative is simply
# d/dx x^(1/3) = 1/(3 y²). Same principle as the Kepler rule below —
# differentiate the solution, not the solver.
# ---------------------------------------------------
@inline _dcbrt(x, y) = inv(3 * y * y)

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{1},
                             func::Const{typeof(vcbrt)}, ::Type{RT},
                             x::Annotation{<:Real}) where {RT}
    y = vcbrt(x.val)
    RT <: Const && return y
    d = x isa Const ? zero(x.val) : x.dval * _dcbrt(x.val, y)
    return RT <: DuplicatedNoNeed ? d : Duplicated(y, d)
end

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{W},
                             func::Const{typeof(vcbrt)}, ::Type{RT},
                             x::Annotation{<:Real}) where {W,RT}
    y = vcbrt(x.val)
    RT <: Const && return y
    g = _dcbrt(x.val, y)
    d = ntuple(i -> x isa Const ? zero(x.val) : x.dval[i] * g, Val(W))
    return RT <: BatchDuplicatedNoNeed ? d : BatchDuplicated(y, d)
end

function EnzymeRules.augmented_primal(config::EnzymeRules.RevConfig,
                                      func::Const{typeof(vcbrt)}, ::Type{<:Active},
                                      x::Annotation{<:Real})
    y = vcbrt(x.val)
    return EnzymeRules.AugmentedReturn(
        EnzymeRules.needs_primal(config) ? y : nothing, nothing, _dcbrt(x.val, y))
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{1},
                             func::Const{typeof(vcbrt)}, dret::Active, tape,
                             x::Annotation{<:Real})
    return (x isa Const ? nothing : dret.val * tape,)
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{W},
                             func::Const{typeof(vcbrt)}, dret::Active, tape,
                             x::Annotation{<:Real}) where {W}
    return (x isa Const ? nothing : ntuple(i -> dret.val[i] * tape, Val(W)),)
end

# ∂E/∂M and ∂E/∂e at the solved E. Applied unconditionally, matching the
# `@scalar_rule` in `kepsolve.jl`: these solver types *are* the elliptical
# branch (`Auto` dispatches e ≥ 1 to `HyperbolicHalley`, which has no rule
# here), so an `e < 1` guard would only add a return-type union that Enzyme
# rejects — forward rules must return a concrete shadow type.
@inline function _kepler_partials(M, e, method)
    E = kepler_solver(M, e, method)
    s, c = sincos(E)
    invu = inv(1 - e * c)
    return E, invu, s * invu
end

# Shadow `i` of an argument, or zero when it carries no derivative. Enzyme
# calls forward rules in batch mode (`Enzyme.gradient(Forward, …)` does),
# where `.dval` is an NTuple rather than a scalar — a rule written for width 1
# fails with a `*(::Float64, ::Tuple)` MethodError.
@inline _sh(x::Const, ::Int) = zero(x.val)
@inline _sh(x::Annotation, i::Int) = x.dval[i]
@inline _sh1(x::Const) = zero(x.val)
@inline _sh1(x::Annotation) = x.dval

# The width is taken from the config *type* so it is a compile-time constant:
# reading it from the value leaves `ntuple` returning `Tuple{Any,Any}`, which
# Enzyme rejects as an incorrect shadow return type.
function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{1},
                             func::Const{typeof(kepler_solver)},
                             ::Type{RT},
                             M::Annotation{<:Real},
                             e::Annotation{<:Real},
                             method::Const{<:EllipticSolver}) where {RT}
    E, dEdM, dEde = _kepler_partials(M.val, e.val, method.val)
    RT <: Const && return E
    d = dEdM * _sh1(M) + dEde * _sh1(e)
    return RT <: DuplicatedNoNeed ? d : Duplicated(E, d)
end

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{W},
                             func::Const{typeof(kepler_solver)},
                             ::Type{RT},
                             M::Annotation{<:Real},
                             e::Annotation{<:Real},
                             method::Const{<:EllipticSolver}) where {W,RT}
    E, dEdM, dEde = _kepler_partials(M.val, e.val, method.val)
    RT <: Const && return E
    d = ntuple(i -> dEdM * _sh(M, i) + dEde * _sh(e, i), Val(W))
    return RT <: BatchDuplicatedNoNeed ? d : BatchDuplicated(E, d)
end

function EnzymeRules.augmented_primal(config::EnzymeRules.RevConfig,
                                      func::Const{typeof(kepler_solver)},
                                      ::Type{<:Active},
                                      M::Annotation{<:Real},
                                      e::Annotation{<:Real},
                                      method::Const{<:EllipticSolver})
    E, dEdM, dEde = _kepler_partials(M.val, e.val, method.val)
    primal = EnzymeRules.needs_primal(config) ? E : nothing
    # Tape only the two scalar partials: cheaper than re-deriving E on the
    # reverse pass, and it keeps the reverse rule free of the solver.
    return EnzymeRules.AugmentedReturn(primal, nothing, (dEdM, dEde))
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{1},
                             func::Const{typeof(kepler_solver)},
                             dret::Active,
                             tape,
                             M::Annotation{<:Real},
                             e::Annotation{<:Real},
                             method::Const{<:EllipticSolver})
    dEdM, dEde = tape
    dM = M isa Const ? nothing : dret.val * dEdM
    de = e isa Const ? nothing : dret.val * dEde
    return (dM, de, nothing)
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{W},
                             func::Const{typeof(kepler_solver)},
                             dret::Active,
                             tape,
                             M::Annotation{<:Real},
                             e::Annotation{<:Real},
                             method::Const{<:EllipticSolver}) where {W}
    dEdM, dEde = tape
    dM = M isa Const ? nothing : ntuple(i -> dret.val[i] * dEdM, Val(W))
    de = e isa Const ? nothing : ntuple(i -> dret.val[i] * dEde, Val(W))
    return (dM, de, nothing)
end

end # module
