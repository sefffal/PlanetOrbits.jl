# Gaia-DR4-shaped workload benchmark — the §11 headline performance gate of
# design/planetorbits-v2-nbody-migration.md: 150 epochs of absolute
# astrometry, 1–3 companions, per-sample value and ForwardDiff-gradient
# evaluations, single core.
#
# Setup (once):
#   julia --project=perf -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
# Run:
#   julia --project=perf perf/gaia-workload.jl
#
# For the per-pass breakdown of the gradient path, see `dual-passes.jl`.
#
# Reference numbers (AMD EPYC Milan / AVX2, Julia 1.12, 2026-07-30 prototype,
# design/prototype-soa/NOTES.md): v0.11 path 70/141/210 µs for Np=1/2/3;
# prototype SIMD 24/28/32 µs; gradient 481→99 µs (Np=1), 2547→497 µs (Np=3).

include("workload.jl")

using BenchmarkTools
using ForwardDiff
using Printf

println("Gaia-DR4-shaped workload: $NEP epochs, value + gradient evals, 1 core")
println("time per eval [µs]  (ns per epoch·companion in parens)")
@printf("%-4s  %12s  %12s  %14s  %6s\n", "Np", "scalar", "simd", "gradient(nθ)", "allocs")
for NP in 1:3
    θ = example_theta(NP)
    np = Val(NP)
    sys = build(θ, np)
    traj = Trajectory(sys, EPOCHS)
    m_scal = KeplerianApprox(simd=false)
    m_simd = KeplerianApprox(simd=true)
    evalonce(θ, np, m_simd, traj)
    allocs = @allocated evalonce(θ, np, m_simd, traj)
    t_scal = @belapsed evalonce($θ, $np, $m_scal, $traj)
    t_simd = @belapsed evalonce($θ, $np, $m_simd, $traj)
    f = θ -> evalonce(θ, np, m_simd)
    g = θ -> ForwardDiff.gradient(f, θ)
    g(θ)
    t_grad = @belapsed $g($θ)
    @printf("%-4d  %9.1f µs  %9.1f µs  %9.1f µs(%2d)  %6d\n",
        NP, 1e6 * t_scal, 1e6 * t_simd, 1e6 * t_grad, length(θ),
        allocs)
    @printf("      %12s  (%8.1f)  \n", "", 1e9 * t_simd / (NEP * NP))
end
