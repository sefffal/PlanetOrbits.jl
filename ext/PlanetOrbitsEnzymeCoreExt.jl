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
# NOTE (2026-08-01, corrected 2026-08-04): Enzyme cannot yet differentiate a
# whole `orbitsolve!` — `solve_row_simd!` raises `EnzymeNoDerivativeError` — so
# this rule is exercised only at the kernel level today. The earlier claim here
# that **AHL21 segfaults is false**: §10.6.6 measured Enzyme 0.13.198 running
# `ahl21_step` cleanly in both modes, with and without runtime activity. See
# §10.4.2 and §10.6.6 of design/planetorbits-v2-nbody-migration.md for the
# measurements behind staying ForwardDiff-first.
#
# The universal-variable γ of the AHL21 Kepler drift has its own rule below.
# ---------------------------------------------------
module PlanetOrbitsEnzymeCoreExt

using PlanetOrbits: PlanetOrbits, kepler_solver, Markley, RootsMethod, _pow23
using EnzymeCore
using EnzymeCore.EnzymeRules
using ForwardDiff: ForwardDiff, Dual, Partials

const EllipticSolver = Union{Markley,RootsMethod}

# ---------------------------------------------------
# `_pow23` (x^(2/3)): the one thing in the SIMD Kepler path Enzyme cannot
# handle.
#
# Its seed is an exponent bit-hack (`reinterpret(Float64, magic - u ÷ 3)`),
# which lowers to an LLVM `bitcast` and raises `EnzymeNoDerivativeError` —
# blocking `solve_row_simd!`, and with it the whole SIMD value path, which is
# ~2× faster than the scalar fallback. Everything else in the kernel
# (`vrem2pi`, `vsincos`, `vsincos_small`) differentiates fine; measured.
#
# The seed's derivative is irrelevant: the cubic iterations that follow
# converge to x^(2/3) whatever the seed, and the derivative is simply
# d/dx x^(2/3) = (2/3)·x^(-1/3) = 2y/(3x). Same principle as the Kepler rule
# below — differentiate the solution, not the solver.
# ---------------------------------------------------
@inline _dpow23(x, y) = 2 * y / (3 * x)

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{1},
                             func::Const{typeof(_pow23)}, ::Type{RT},
                             x::Annotation{<:Real}) where {RT}
    y = _pow23(x.val)
    RT <: Const && return y
    d = x isa Const ? zero(x.val) : x.dval * _dpow23(x.val, y)
    return RT <: DuplicatedNoNeed ? d : Duplicated(y, d)
end

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{W},
                             func::Const{typeof(_pow23)}, ::Type{RT},
                             x::Annotation{<:Real}) where {W,RT}
    y = _pow23(x.val)
    RT <: Const && return y
    g = _dpow23(x.val, y)
    d = ntuple(i -> x isa Const ? zero(x.val) : x.dval[i] * g, Val(W))
    return RT <: BatchDuplicatedNoNeed ? d : BatchDuplicated(y, d)
end

function EnzymeRules.augmented_primal(config::EnzymeRules.RevConfig,
                                      func::Const{typeof(_pow23)}, ::Type{<:Active},
                                      x::Annotation{<:Real})
    y = _pow23(x.val)
    return EnzymeRules.AugmentedReturn(
        EnzymeRules.needs_primal(config) ? y : nothing, nothing, _dpow23(x.val, y))
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{1},
                             func::Const{typeof(_pow23)}, dret::Active, tape,
                             x::Annotation{<:Real})
    return (x isa Const ? nothing : dret.val * tape,)
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{W},
                             func::Const{typeof(_pow23)}, dret::Active, tape,
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

# ---------------------------------------------------
# `_solve_gamma`: the universal-variable anomaly γ of the AHL21 Kepler drift.
#
# γ is defined implicitly by F(γ; p) = 0 (Kepler's equation in universal
# variables), so ∂γ/∂pᵢ = −(∂F/∂pᵢ)/(∂F/∂γ) — and ∂F/∂γ is already returned by
# `_gamma_F_dF` as `dF`, so the rule costs one `_gamma_F_dF` evaluation instead
# of a Newton loop.
#
# Without it Enzyme differentiates the loop itself. Measured at k = 2.5–3.3
# iterations (§10.6.8) that is a ~2.6–2.9× arithmetic penalty for batched
# forward at W = 21 — but the bigger objection is structural: the loop has a
# **dynamic trip count**, which forces a dynamically-sized cache inside the
# differentiated region (§10.6.9). This rule deletes that loop from the region.
#
# Rather than transcribe eight partials of F by hand, the rule seeds
# `_gamma_F_dF` with ForwardDiff `Dual`s and reads `partials(−F/dF)` — which is
# exactly what `_solve_gamma`'s own `Dual` method does in `nbody-kernels.jl`.
# The two rules therefore agree by *construction* rather than by transcription,
# including the `value(F) ≈ 0` term both drop (~1 ulp at the converged root).
# `perf/nbody-enzyme.jl`'s (A) case gates that agreement.
#
# `r0inv` is an argument of `_solve_gamma` but never enters `F` — it appears
# only in the initial guess, which the converged root does not depend on — so
# ∂γ/∂r0inv ≡ 0.
# ---------------------------------------------------
struct _GammaTag end

# F's parameters, in `_gamma_F_dF` order (i.e. `_solve_gamma`'s minus `r0inv`).
@inline _gamma_p(gm, r0, beta0, signb, sqb, zeta, eta, h) =
    promote(gm, r0, beta0, signb, sqb, zeta, eta, h)

# Seed the eight parameters with W-wide tangents and read off dγ. One
# `_gamma_F_dF` evaluation covers all W directions at once.
@inline function _gamma_tangents(gamma, p::NTuple{8,T},
                                 dp::NTuple{8,NTuple{W,T}}) where {T,W}
    d = ntuple(k -> Dual{_GammaTag}(p[k], Partials{W,T}(dp[k])), Val(8))
    F, dF = PlanetOrbits._gamma_F_dF(gamma, d...)
    delta = -F / dF
    return ntuple(i -> ForwardDiff.partials(delta, i), Val(W))
end

@inline _gamma_primal(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h) =
    PlanetOrbits._solve_gamma(gm.val, r0.val, r0inv.val, beta0.val, signb.val,
                              sqb.val, zeta.val, eta.val, h.val)

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{1},
                             func::Const{typeof(PlanetOrbits._solve_gamma)},
                             ::Type{RT},
                             gm::Annotation{<:Real}, r0::Annotation{<:Real},
                             r0inv::Annotation{<:Real}, beta0::Annotation{<:Real},
                             signb::Annotation{<:Real}, sqb::Annotation{<:Real},
                             zeta::Annotation{<:Real}, eta::Annotation{<:Real},
                             h::Annotation{<:Real}) where {RT}
    gamma = _gamma_primal(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
    RT <: Const && return gamma
    p = _gamma_p(gm.val, r0.val, beta0.val, signb.val, sqb.val, zeta.val, eta.val, h.val)
    s = _gamma_p(_sh1(gm), _sh1(r0), _sh1(beta0), _sh1(signb),
                 _sh1(sqb), _sh1(zeta), _sh1(eta), _sh1(h))
    d = _gamma_tangents(gamma, p, ntuple(k -> (s[k],), Val(8)))[1]
    return RT <: DuplicatedNoNeed ? d : Duplicated(gamma, d)
end

function EnzymeRules.forward(config::EnzymeRules.FwdConfigWidth{W},
                             func::Const{typeof(PlanetOrbits._solve_gamma)},
                             ::Type{RT},
                             gm::Annotation{<:Real}, r0::Annotation{<:Real},
                             r0inv::Annotation{<:Real}, beta0::Annotation{<:Real},
                             signb::Annotation{<:Real}, sqb::Annotation{<:Real},
                             zeta::Annotation{<:Real}, eta::Annotation{<:Real},
                             h::Annotation{<:Real}) where {W,RT}
    gamma = _gamma_primal(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
    RT <: Const && return gamma
    p = _gamma_p(gm.val, r0.val, beta0.val, signb.val, sqb.val, zeta.val, eta.val, h.val)
    args = (gm, r0, beta0, signb, sqb, zeta, eta, h)
    dp = ntuple(k -> ntuple(i -> oftype(p[k], _sh(args[k], i)), Val(W)), Val(8))
    d = _gamma_tangents(gamma, p, dp)
    return RT <: BatchDuplicatedNoNeed ? d : BatchDuplicated(gamma, d)
end

function EnzymeRules.augmented_primal(config::EnzymeRules.RevConfig,
                                      func::Const{typeof(PlanetOrbits._solve_gamma)},
                                      ::Type{<:Active},
                                      gm::Annotation{<:Real}, r0::Annotation{<:Real},
                                      r0inv::Annotation{<:Real}, beta0::Annotation{<:Real},
                                      signb::Annotation{<:Real}, sqb::Annotation{<:Real},
                                      zeta::Annotation{<:Real}, eta::Annotation{<:Real},
                                      h::Annotation{<:Real})
    gamma = _gamma_primal(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
    p = _gamma_p(gm.val, r0.val, beta0.val, signb.val, sqb.val, zeta.val, eta.val, h.val)
    T = typeof(p[1])
    # Tape the eight partials, not the solve: the reverse pass is then a
    # contraction, and never re-enters the solver. One 8-wide `_gamma_F_dF`.
    g = _gamma_tangents(gamma, p, ntuple(k -> ntuple(l -> T(l == k), Val(8)), Val(8)))
    primal = EnzymeRules.needs_primal(config) ? gamma : nothing
    return EnzymeRules.AugmentedReturn(primal, nothing, g)
end

@inline _gadj(x::Const, gi, λ) = nothing
@inline _gadj(x::Annotation, gi, λ) = λ * gi
@inline _gadjW(x::Const, gi, λ, ::Val{W}) where {W} = nothing
@inline _gadjW(x::Annotation, gi, λ, ::Val{W}) where {W} =
    ntuple(i -> λ[i] * gi, Val(W))

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{1},
                             func::Const{typeof(PlanetOrbits._solve_gamma)},
                             dret::Active, tape,
                             gm::Annotation{<:Real}, r0::Annotation{<:Real},
                             r0inv::Annotation{<:Real}, beta0::Annotation{<:Real},
                             signb::Annotation{<:Real}, sqb::Annotation{<:Real},
                             zeta::Annotation{<:Real}, eta::Annotation{<:Real},
                             h::Annotation{<:Real})
    g = tape
    λ = dret.val
    z = zero(eltype(g))
    return (_gadj(gm, g[1], λ), _gadj(r0, g[2], λ), _gadj(r0inv, z, λ),
            _gadj(beta0, g[3], λ), _gadj(signb, g[4], λ), _gadj(sqb, g[5], λ),
            _gadj(zeta, g[6], λ), _gadj(eta, g[7], λ), _gadj(h, g[8], λ))
end

function EnzymeRules.reverse(config::EnzymeRules.RevConfigWidth{W},
                             func::Const{typeof(PlanetOrbits._solve_gamma)},
                             dret::Active, tape,
                             gm::Annotation{<:Real}, r0::Annotation{<:Real},
                             r0inv::Annotation{<:Real}, beta0::Annotation{<:Real},
                             signb::Annotation{<:Real}, sqb::Annotation{<:Real},
                             zeta::Annotation{<:Real}, eta::Annotation{<:Real},
                             h::Annotation{<:Real}) where {W}
    g = tape
    λ = dret.val
    z = zero(eltype(g))
    V = Val(W)
    return (_gadjW(gm, g[1], λ, V), _gadjW(r0, g[2], λ, V), _gadjW(r0inv, z, λ, V),
            _gadjW(beta0, g[3], λ, V), _gadjW(signb, g[4], λ, V), _gadjW(sqb, g[5], λ, V),
            _gadjW(zeta, g[6], λ, V), _gadjW(eta, g[7], λ, V), _gadjW(h, g[8], λ, V))
end

end # module
