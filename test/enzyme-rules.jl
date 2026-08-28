# ---------------------------------------------------
# EnzymeRules extension gate.
#
# Skipped unless Enzyme is loadable, so it costs CI nothing by default —
# Enzyme is a heavy dependency and (as of 2026-08-01) cannot differentiate a
# whole `orbitsolve!` anyway; see §10.4.2 of the design doc. Run it with an
# environment that has Enzyme to exercise `ext/PlanetOrbitsEnzymeCoreExt.jl`.
#
# Referee is BigFloat, via the implicit relation E − e sin E = M solved to
# 256-bit precision, not the Float64 solver.
# ---------------------------------------------------
if isnothing(Base.identify_package("Enzyme"))
    @info "skipping EnzymeRules gate: Enzyme not available in this environment"
else
    @eval using Enzyme
    @testset "EnzymeRules: analytic Kepler derivatives" begin
        setprecision(BigFloat, 256)
        @test !isnothing(Base.get_extension(PlanetOrbits, :PlanetOrbitsEnzymeCoreExt))

        # ∂E/∂M, ∂E/∂e from a high-precision Newton solve of Kepler's equation.
        function ref_partials(M, e)
            Mb = big(M); eb = big(e)
            E = Mb
            for _ in 1:200
                E -= (E - eb*sin(E) - Mb) / (1 - eb*cos(E))
            end
            u = 1 - eb*cos(E)
            return (1/u, sin(E)/u)
        end

        f(x) = PlanetOrbits.kepler_solver(x[1], x[2], PlanetOrbits.Markley())
        for M in (0.05, 0.7, 2.0, 3.1, -1.2), e in (0.0, 0.1, 0.5, 0.9, 0.99)
            gf = collect(Enzyme.gradient(Forward, f, [M, e])[1])
            gr = Enzyme.gradient(Reverse, f, [M, e])[1]
            t = ref_partials(M, e)
            # Compared on an absolute scale: ∂E/∂e = sin(E)/u passes through
            # zero at M = 0 and M = π, where a relative test is meaningless.
            for (a, r) in ((gf[1], t[1]), (gr[1], t[1]), (gf[2], t[2]), (gr[2], t[2]))
                @test isapprox(big(a), r; rtol=1e-12, atol=1e-13)
            end
            @test gf ≈ gr                       # forward and reverse agree
        end

        # `_pow23`'s exponent bit-hack seed lowers to an LLVM `bitcast`, which
        # raises EnzymeNoDerivativeError and blocks `solve_row_simd!` — and
        # with it the entire SIMD value path. The rule replaces it with
        # d/dx x^(2/3) = 2/(3 cbrt(x)). Everything else in the SIMD kernel
        # (`vrem2pi`, `vsincos`, `vsincos_small`) differentiates natively.
        @testset "_pow23 unblocks the SIMD Kepler kernel" begin
            for x in (1e-8, 0.5, 1.0, 8.0, 1e3, 1e8)
                g = Enzyme.autodiff(Reverse, y -> PlanetOrbits._pow23(y), Active, Active(x))[1][1]
                @test isapprox(g, 2/(3*cbrt(x)); rtol=1e-10)
            end
            # the kernel that was blocked now differentiates end to end
            for (M, e) in ((0.7, 0.3), (2.0, 0.7), (-1.2, 0.05))
                g = Enzyme.autodiff(Reverse,
                        m -> PlanetOrbits.markley_sincosE(m, e)[1], Active, Active(M))[1][1]
                ref = ForwardDiff.derivative(
                        m -> sin(PlanetOrbits.kepler_solver(m, e, PlanetOrbits.Markley())), M)
                @test isapprox(g, ref; rtol=1e-12)
            end
        end

        # The rule's reason for existing: `kepler_solver` short-circuits on
        # `iszero(e)` and returns M, so differentiating *through* it yields
        # ∂E/∂e = 0 at e = 0, where the true value is sin(M). Circular-orbit
        # fits sit exactly here.
        for M in (0.05, 0.7, 2.0, -1.2)
            g = Enzyme.gradient(Reverse, f, [M, 0.0])[1]
            @test isapprox(g[2], sin(M); rtol=1e-13)
            @test !isapprox(g[2], 0.0; atol=1e-8)   # the bug this rule removes
        end

        # ---------------------------------------------------
        # `_solve_gamma`: the universal-variable anomaly of the AHL21 Kepler
        # drift. Without the rule Enzyme differentiates the Newton loop — still
        # *correct*, so accuracy alone cannot detect a rule that silently fails
        # to apply. The registration check is what catches that: it asks the
        # same predicates Enzyme itself consults, against the concrete
        # signature, so an argument-type mismatch fails here rather than
        # quietly reverting to the slow path.
        # ---------------------------------------------------
        @testset "_solve_gamma implicit rule" begin
            sig = Tuple{typeof(PlanetOrbits._solve_gamma), ntuple(_ -> Float64, 9)...}
            @test Enzyme.EnzymeRules.has_frule_from_sig(sig)
            @test Enzyme.EnzymeRules.has_rrule_from_sig(sig)

            # `_delxv_gamma`'s prologue, so the rule is exercised at arguments
            # the propagator actually produces. `drift_first` selects whether
            # the removed drift precedes the Kepler step.
            function gamma_of(u::AbstractVector{T}, drift_first::Bool) where {T}
                x0 = @view u[1:3]
                v0 = @view u[4:6]
                gm = u[7]; h = u[8]
                rtmp = x0 .- h .* v0
                r0 = drift_first ? sqrt(sum(abs2, rtmp)) : sqrt(sum(abs2, x0))
                r0inv = inv(r0)
                beta0 = 2 * gm * r0inv - sum(abs2, v0)
                signb = sign(beta0)
                sqb = sqrt(signb * beta0)
                zeta = gm - r0 * beta0
                eta = drift_first ? sum(rtmp .* v0) : sum(x0 .* v0)
                return PlanetOrbits._solve_gamma(gm, r0, r0inv, beta0, signb,
                                                 sqb, zeta, eta, h)
            end

            # Elliptic (β > 0) and hyperbolic (β < 0), both drift orders.
            states = ([0.09, 0.02, 0.001, -0.5, 2.3, 0.05, 0.0002959, 0.6909],
                      [0.13, -0.07, 0.02, 1.1, -0.9, 0.2, 0.0002959, 0.2303],
                      [0.09, 0.02, 0.001, -8.0, 12.0, 0.5, 0.0002959, 0.05])
            for u in states, df in (true, false)
                g(x) = gamma_of(x, df)
                gfd = ForwardDiff.gradient(g, u)           # `_solve_gamma`'s Dual rule
                gef = collect(Enzyme.gradient(Forward, g, u)[1])
                ger = Enzyme.gradient(Reverse, g, u)[1]
                s = maximum(abs, gfd)
                @test maximum(abs, gef .- gfd) / s < 1e-12
                @test maximum(abs, ger .- gfd) / s < 1e-12
            end

            # Both rules apply the same implicit relation, so agreeing with
            # each other would not catch a shared error. BigFloat central
            # differences on the whole map are the independent referee; one
            # state is enough, since what could be wrong is the argument
            # mapping, not the state.
            u = states[1]
            gfd = ForwardDiff.gradient(x -> gamma_of(x, true), u)
            s = maximum(abs, gfd)
            for i in eachindex(u)
                δ = abs(big(u[i])) * big(2.0)^-40
                up = big.(u); up[i] += δ
                dn = big.(u); dn[i] -= δ
                ref = (gamma_of(up, true) - gamma_of(dn, true)) / (2δ)
                @test abs(big(gfd[i]) - ref) / s < 1e-11
            end
        end
    end
end
