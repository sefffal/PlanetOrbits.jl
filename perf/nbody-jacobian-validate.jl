# Validate and cost the analytic pair Jacobian (src/nbody-jacobian.jl).
#
# Run:
#   julia --project=perf perf/nbody-jacobian-validate.jl
#
# Three questions, in order:
#   1. Is the port correct?  Analytic J vs BigFloat(256) ForwardDiff of the
#      same primal. A transcription slip shows up as a single bad column.
#   2. Does upstream's separate `jac_mass` actually recover the digits the
#      naive ∂/∂m_i combination loses? (nbody-mass-derivative.jl says the
#      naive form is down to ~11 digits at γ = 0.01.)
#   3. What does the assembly cost, against the budget the chain-rule probe
#      left for it?

using StaticArrays
using ForwardDiff
using BenchmarkTools
using Printf
using PlanetOrbits: _delxv_gamma, _solve_gamma, G_au3_day2

include(joinpath(dirname(@__DIR__), "src", "nbody-jacobian.jl"))

setprecision(BigFloat, 256)

function pair_kernel(u::AbstractVector{T}, drift_first::Bool) where {T}
    x0 = SVector{3,T}(u[1], u[2], u[3])
    v0 = SVector{3,T}(u[4], u[5], u[6])
    dx, dv = _delxv_gamma(x0, v0, u[7], u[8], drift_first)
    return SVector{6,T}(dx[1], dx[2], dx[3], dv[1], dv[2], dv[3])
end

function achieved_gamma(x0, v0, gm, h, drift_first)
    rtmp = x0 - h * v0
    r0 = drift_first ? sqrt(sum(abs2, rtmp)) : sqrt(sum(abs2, x0))
    beta0 = 2gm / r0 - sum(abs2, v0)
    signb = sign(beta0); sqb = sqrt(signb * beta0)
    zeta = gm - r0 * beta0
    eta = drift_first ? sum(rtmp .* v0) : sum(x0 .* v0)
    return _solve_gamma(gm, r0, inv(r0), beta0, signb, sqb, zeta, eta, h)
end

function pair_case(γ_target, ecc, mratio, drift_first)
    mi, mj = 1.0, mratio
    gm = G_au3_day2 * (mi + mj)
    a = 0.1
    r = a * abs(1 - ecc)
    x0 = SVector(r, 0.0, 0.0)
    vmag = sqrt(gm * (1 + ecc) / r)
    v0 = SVector(0.0, vmag * 0.9, vmag * 0.1)
    beta0 = 2gm / r - sum(abs2, v0)
    lo, hi = 0.0, γ_target * r / sqrt(abs(beta0))
    for _ in 1:200
        achieved_gamma(x0, v0, gm, hi, drift_first) >= γ_target && break
        lo = hi; hi *= 1.5
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        achieved_gamma(x0, v0, gm, mid, drift_first) < γ_target ? (lo = mid) : (hi = mid)
    end
    return (mi = mi, mj = mj, x0 = x0, v0 = v0, gm = gm, h = 0.5 * (lo + hi))
end

const CASES = [(γ, e, df) for γ in (0.01, 0.1, 0.5, 1.0, 2.0, 3.0)
                          for e in (0.0, 0.6, 0.95, 1.4)
                          for df in (true, false)]

function check_jacobian()
    println("### 1. Analytic pair Jacobian vs BigFloat(256), relative to max|J|")
    @printf("%6s %6s %6s | %10s %10s | %10s\n",
            "γ", "e", "drift1", "analytic", "ForwardDiff", "primal ===")
    println("-"^62)
    worst = 0.0
    for (γt, ecc, df) in CASES
        c = pair_case(γt, ecc, 1e-3, df)
        u64 = [c.x0..., c.v0..., c.gm, c.h]
        Jbig = ForwardDiff.jacobian(u -> pair_kernel(u, df), BigFloat.(u64))
        Jfd = ForwardDiff.jacobian(u -> pair_kernel(u, df), u64)
        dx, dv, Jan, _ = _delxv_gamma_jac(c.x0, c.v0, c.gm, c.h, df; G=G_au3_day2)
        pdx, pdv = _delxv_gamma(c.x0, c.v0, c.gm, c.h, df)
        scale = maximum(abs, Jbig)
        ean = Float64(maximum(abs, BigFloat.(Jan) .- Jbig) / scale)
        efd = Float64(maximum(abs, BigFloat.(Jfd) .- Jbig) / scale)
        worst = max(worst, ean)
        @printf("%6.2f %6.2f %6s | %10.2e %10.2e | %10s\n", γt, ecc, df, ean, efd,
                (dx === pdx && dv === pdv) ? "yes" : "NO")
    end
    @printf("\nworst analytic error over the grid: %.2e\n", worst)
end

function check_mass()
    println("\n### 2. ∂(m_j/(m_i+m_j)·Δx)/∂m_i — naive AD vs upstream's closed form")
    println("`cancel` = operand/result; relerr relative to the result.")
    @printf("%6s %6s %6s %8s | %8s | %11s %11s\n",
            "γ", "e", "drift1", "m_j/m_i", "cancel", "AD relerr", "closed relerr")
    println("-"^76)
    for mratio in (1e-3, 1e-8), (γt, ecc, df) in
            [(γ, e, d) for γ in (0.01, 0.1, 0.5, 2.0) for e in (0.0, 0.6) for d in (true, false)]
        c = pair_case(γt, ecc, mratio, df)
        # Reference and the naive AD combination, both on the same expression.
        incr(mi_, T) = begin
            gm = T(G_au3_day2) * (mi_ + T(c.mj))
            dx, dv = _delxv_gamma(SVector{3,T}(c.x0), SVector{3,T}(c.v0), gm, T(c.h), df)
            w = T(c.mj) / (mi_ + T(c.mj))
            SVector{6}(w*dx[1], w*dx[2], w*dx[3], w*dv[1], w*dv[2], w*dv[3])
        end
        dbig = ForwardDiff.derivative(m -> incr(m, BigFloat), big(c.mi))
        d64 = ForwardDiff.derivative(m -> incr(m, Float64), c.mi)
        _, _, _, jm = _delxv_gamma_jac(c.x0, c.v0, c.gm, c.h, df; G=G_au3_day2)
        dclosed = jm .* c.mj                      # upstream's jac_ij[:,7]
        res = maximum(abs, dbig)
        dx, dv = _delxv_gamma(c.x0, c.v0, c.gm, c.h, df)
        opd = c.mj / (c.mi + c.mj)^2 * maximum(abs, SVector(dx..., dv...))
        @printf("%6.2f %6.2f %6s %8.0e | %8.1e | %11.2e %11.2e\n",
                γt, ecc, df, mratio, opd / Float64(res),
                Float64(maximum(abs, BigFloat.(d64) .- dbig) / res),
                Float64(maximum(abs, BigFloat.(dclosed) .- dbig) / res))
    end
end

function check_cost()
    println("\n### 3. Assembly cost: analytic 6×8 Jacobian vs the primal it comes from")
    c = pair_case(0.5, 0.6, 1e-3, true)
    for df in (true, false)
        bp = @benchmark _delxv_gamma($(c.x0), $(c.v0), $(c.gm), $(c.h), $df) samples=2000 evals=10
        ba = @benchmark _delxv_gamma_jac($(c.x0), $(c.v0), $(c.gm), $(c.h), $df;
                                         G=$G_au3_day2) samples=2000 evals=10
        tp = minimum(bp).time * 1e-3
        ta = minimum(ba).time * 1e-3
        @printf("drift_first=%-5s  primal %7.4f µs   primal+analytic J %7.4f µs   ratio %5.2f×\n",
                df, tp, ta, ta / tp)
    end
end

check_jacobian()
check_mass()
check_cost()
