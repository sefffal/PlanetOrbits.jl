# Gate for the `_solve_gamma` EnzymeRules rule (§10.6.8/§10.6.9).
#
#   julia --project=<env-with-Enzyme> perf/nbody-gamma-rule.jl
#
# Correctness alone cannot tell you the rule fired: differentiating *through*
# the Newton loop also gives the right answer (that is exactly what Enzyme was
# doing before, at 2.8e-16). So this checks three separate things.
#
#  1. **Registration.** `EnzymeRules.has_frule_from_sig` / `has_rrule_from_sig`
#     are the predicates Enzyme itself consults, so asking them about the real
#     call signature catches the most likely failure — a rule whose argument
#     types don't quite match and is therefore silently never used.
#  2. **Agreement.** Enzyme (forward, batched forward, reverse) against
#     ForwardDiff, which reaches `_solve_gamma`'s own `Dual` method. Both apply
#     the same implicit rule and drop the same `value(F) ≈ 0` term, so they
#     should agree far below the BigFloat-truth error — that tight agreement is
#     itself evidence the same code path is being used.
#  3. **Truth.** Both against BigFloat(256) central differences, since agreeing
#     with each other proves nothing if the shared rule is wrong.

using StaticArrays
using ForwardDiff
using BenchmarkTools
using Printf
using PlanetOrbits: PlanetOrbits, _solve_gamma, _delxv_gamma, G_au3_day2
using Enzyme
using Enzyme.EnzymeRules: has_frule_from_sig, has_rrule_from_sig

setprecision(BigFloat, 256)

# `_solve_gamma`'s nine arguments as a function of a pair state, so the probe
# exercises the rule at arguments the propagator actually produces. Mirrors
# `_delxv_gamma`'s prologue.
function gamma_of(u::AbstractVector{T}) where {T}
    x0 = SVector{3,T}(u[1], u[2], u[3])
    v0 = SVector{3,T}(u[4], u[5], u[6])
    gm = u[7]
    h = u[8]
    rtmp = x0 - h * v0
    r0 = sqrt(sum(abs2, rtmp))
    r0inv = inv(r0)
    beta0 = 2 * gm * r0inv - sum(abs2, v0)
    signb = sign(beta0)
    sqb = sqrt(signb * beta0)
    zeta = gm - r0 * beta0
    eta = sum(rtmp .* v0)
    return _solve_gamma(gm, r0, r0inv, beta0, signb, sqb, zeta, eta, h)
end

const U = [0.09, 0.02, 0.001, -0.5, 2.3, 0.05,
           G_au3_day2 * (1.071 + 1.36e-5), 0.69090]

# BigFloat central differences on the *whole* map, per-group scaled so an
# exactly-zero component never sets the step.
function fd_reference(u)
    b = big.(u)
    scales = (maximum(abs, b[1:3]), maximum(abs, b[1:3]), maximum(abs, b[1:3]),
              maximum(abs, b[4:6]), maximum(abs, b[4:6]), maximum(abs, b[4:6]),
              abs(b[7]), abs(b[8]))
    map(1:8) do i
        δ = scales[i] * big(2.0)^-40
        up = copy(b); up[i] += δ
        dn = copy(b); dn[i] -= δ
        (gamma_of(up) - gamma_of(dn)) / (2δ)
    end
end

function main()
    # (1) Registration — against the concrete signature Enzyme will see.
    sig = Tuple{typeof(_solve_gamma), ntuple(_ -> Float64, 9)...}
    @printf("rule registered:  forward %s   reverse %s\n",
            has_frule_from_sig(sig) ? "yes" : "*** NO ***",
            has_rrule_from_sig(sig) ? "yes" : "*** NO ***")

    ref = fd_reference(U)
    scale = maximum(abs, ref)
    rel(g) = Float64(maximum(abs, BigFloat.(collect(g)) .- ref) / scale)

    gfd = ForwardDiff.gradient(gamma_of, U)
    gef = collect(Enzyme.gradient(Forward, gamma_of, U)[1])
    ger = Enzyme.gradient(Reverse, gamma_of, U)[1]

    @printf("\n%-28s %12s %12s\n", "", "vs BigFloat", "vs ForwardDiff")
    @printf("%-28s %12.2e %12s\n", "ForwardDiff (Dual rule)", rel(gfd), "—")
    @printf("%-28s %12.2e %12.2e\n", "Enzyme forward", rel(gef),
            maximum(abs, gef .- gfd) / maximum(abs, gfd))
    @printf("%-28s %12.2e %12.2e\n", "Enzyme reverse", rel(ger),
            maximum(abs, ger .- gfd) / maximum(abs, gfd))

    println("\ntimings (γ solve + prologue, 8 parameters)")
    for (label, f) in (("ForwardDiff", () -> ForwardDiff.gradient(gamma_of, U)),
                       ("Enzyme forward", () -> Enzyme.gradient(Forward, gamma_of, U)),
                       ("Enzyme reverse", () -> Enzyme.gradient(Reverse, gamma_of, U)))
        b = @benchmark $f() samples=200 evals=1
        @printf("  %-22s %9.3f µs\n", label, minimum(b).time * 1e-3)
    end
end

main()
