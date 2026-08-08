# The evidence behind design §10.4.7: why the observing pass is *not* short of
# SIMD under Duals, and what actually moves the gradient instead.
#
# Run:
#   julia --project=perf perf/dual-layout.jl
#
# Three measurements:
#   1. SoA vs AoS on `_rotate_pass!`'s arithmetic — the structure-of-arrays Dual
#      layout §10.4.1 proposed, prototyped directly. It wins only at N=1.
#   2. Vector-instruction counts in the emitted IR — ForwardDiff's Partials are
#      an NTuple and LLVM already packs them, so there is no lost SIMD.
#   3. Gradient cost vs ForwardDiff chunk width — the lever that *does* pay.

include("workload.jl")

using BenchmarkTools, ForwardDiff, Printf, InteractiveUtils
using ForwardDiff: Dual, Partials
using PlanetOrbits: _rotate_pass!, _observe_pass!

const METHOD = KeplerianApprox(simd=true)
const NB = 2

# Seed θ the way a chunked ForwardDiff.gradient does. Partial *values* never
# affect timing (nothing branches on them).
seed(θ::Vector{Float64}, ::Val{N}) where {N} =
    [Dual{Nothing}(θ[i], Partials(ntuple(j -> Float64(j == mod1(i, N)), Val(N))))
     for i in eachindex(θ)]

# ---------------------------------------------------------------------------
# 1. SoA vs AoS
#
# `_rotate_pass!` is the cleanest test case in the observing pass: pure
# multiply-add, no divide, no sqrt, no branch — and it carries the pass's worst
# Dual penalty. If SoA cannot beat AoS here it cannot beat it anywhere.
# ---------------------------------------------------------------------------

function rotate_aos!(x, y, z, R, ::Val{NB}) where {NB}
    nk = size(x, 1)
    @inbounds for j in 1:NB
        @simd for k in 1:nk
            r11 = R[1][k]; r12 = R[2][k]; r13 = R[3][k]
            r21 = R[4][k]; r22 = R[5][k]; r23 = R[6][k]
            r31 = R[7][k]; r32 = R[8][k]; r33 = R[9][k]
            xk = x[k, j]; yk = y[k, j]; zk = z[k, j]
            x[k, j] = r11 * xk + r12 * yk + r13 * zk
            y[k, j] = r21 * xk + r22 * yk + r23 * zk
            z[k, j] = r31 * xk + r32 * yk + r33 * zk
        end
    end
    return x
end

# Value and partial streams as separate Float64 arrays; partial index is the
# slow axis so the inner loop over epochs stays contiguous.
function rotate_soa!(xv, xp, yv, yp, zv, zp, Rv, Rp, ::Val{NB}, ::Val{N}) where {NB,N}
    nk = size(xv, 1)
    @inbounds for j in 1:NB
        # Partials first: they read the old values, which the value stream
        # below overwrites. This is also where SoA loses — each partial has to
        # re-read both value operands, because the operation is bilinear.
        for p in 1:N
            @simd for k in 1:nk
                r11 = Rv[1][k]; r12 = Rv[2][k]; r13 = Rv[3][k]
                r21 = Rv[4][k]; r22 = Rv[5][k]; r23 = Rv[6][k]
                r31 = Rv[7][k]; r32 = Rv[8][k]; r33 = Rv[9][k]
                xk = xv[k, j]; yk = yv[k, j]; zk = zv[k, j]
                dx = xp[k, p, j]; dy = yp[k, p, j]; dz = zp[k, p, j]
                xp[k, p, j] = r11 * dx + r12 * dy + r13 * dz +
                              Rp[1][k, p] * xk + Rp[2][k, p] * yk + Rp[3][k, p] * zk
                yp[k, p, j] = r21 * dx + r22 * dy + r23 * dz +
                              Rp[4][k, p] * xk + Rp[5][k, p] * yk + Rp[6][k, p] * zk
                zp[k, p, j] = r31 * dx + r32 * dy + r33 * dz +
                              Rp[7][k, p] * xk + Rp[8][k, p] * yk + Rp[9][k, p] * zk
            end
        end
        @simd for k in 1:nk
            r11 = Rv[1][k]; r12 = Rv[2][k]; r13 = Rv[3][k]
            r21 = Rv[4][k]; r22 = Rv[5][k]; r23 = Rv[6][k]
            r31 = Rv[7][k]; r32 = Rv[8][k]; r33 = Rv[9][k]
            xk = xv[k, j]; yk = yv[k, j]; zk = zv[k, j]
            xv[k, j] = r11 * xk + r12 * yk + r13 * zk
            yv[k, j] = r21 * xk + r22 * yk + r23 * zk
            zv[k, j] = r31 * xk + r32 * yk + r33 * zk
        end
    end
    return xv
end

mkdual(::Val{N}) where {N} = Dual{Nothing}(rand(), Partials(ntuple(_ -> rand(), Val(N))))

function soa_case(::Val{N}) where {N}
    D = Dual{Nothing,Float64,N}
    x = D[mkdual(Val(N)) for _ in 1:NEP, _ in 1:NB]
    y = copy(x); z = copy(x)
    R = ntuple(_ -> D[mkdual(Val(N)) for _ in 1:NEP], 9)
    xv = [ForwardDiff.value(x[k, j]) for k in 1:NEP, j in 1:NB]
    yv = [ForwardDiff.value(y[k, j]) for k in 1:NEP, j in 1:NB]
    zv = [ForwardDiff.value(z[k, j]) for k in 1:NEP, j in 1:NB]
    xp = [ForwardDiff.partials(x[k, j], p) for k in 1:NEP, p in 1:N, j in 1:NB]
    yp = [ForwardDiff.partials(y[k, j], p) for k in 1:NEP, p in 1:N, j in 1:NB]
    zp = [ForwardDiff.partials(z[k, j], p) for k in 1:NEP, p in 1:N, j in 1:NB]
    Rv = ntuple(c -> [ForwardDiff.value(R[c][k]) for k in 1:NEP], 9)
    Rp = ntuple(c -> [ForwardDiff.partials(R[c][k], p) for k in 1:NEP, p in 1:N], 9)
    return (; x, y, z, R, xv, yv, zv, xp, yp, zp, Rv, Rp)
end

println("1. SoA vs AoS Dual layout, `_rotate_pass!` arithmetic")
println("   $NEP epochs, $NB bodies, µs\n")
@printf("   %-4s %9s %9s %9s | %9s %9s\n",
        "N", "AoS", "SoA", "Float64", "SoA gain", "AoS pen.")
for N in (1, 4, 8, 12)
    c = soa_case(Val(N))
    # The two must agree exactly before either timing means anything.
    xa = copy(c.x); ya = copy(c.y); za = copy(c.z)
    xv = copy(c.xv); yv = copy(c.yv); zv = copy(c.zv)
    xp = copy(c.xp); yp = copy(c.yp); zp = copy(c.zp)
    rotate_aos!(xa, ya, za, c.R, Val(NB))
    rotate_soa!(xv, xp, yv, yp, zv, zp, c.Rv, c.Rp, Val(NB), Val(N))
    ev = maximum(abs, [ForwardDiff.value(xa[k, j]) - xv[k, j] for k in 1:NEP, j in 1:NB])
    ep = maximum(abs, [ForwardDiff.partials(xa[k, j], p) - xp[k, p, j]
                       for k in 1:NEP, p in 1:N, j in 1:NB])
    @assert ev < 1e-12 && ep < 1e-12 "SoA disagrees with AoS: $ev $ep"

    t_aos = @belapsed rotate_aos!(x, y, z, $(c.R), $(Val(NB))) seconds = 1 setup =
        (x = copy($(c.x)); y = copy($(c.y)); z = copy($(c.z)))
    t_soa = @belapsed rotate_soa!(xv, xp, yv, yp, zv, zp, $(c.Rv), $(c.Rp),
                                  $(Val(NB)), $(Val(N))) seconds = 1 setup =
        (xv = copy($(c.xv)); yv = copy($(c.yv)); zv = copy($(c.zv));
         xp = copy($(c.xp)); yp = copy($(c.yp)); zp = copy($(c.zp)))
    Rf = ntuple(_ -> rand(NEP), 9)
    t_f64 = @belapsed rotate_aos!(x, y, z, $Rf, $(Val(NB))) seconds = 1 setup =
        (x = rand(NEP, NB); y = rand(NEP, NB); z = rand(NEP, NB))
    @printf("   %-4d %9.2f %9.2f %9.2f | %8.2fx %8.1fx\n",
            N, 1e6t_aos, 1e6t_soa, 1e6t_f64, t_aos / t_soa, t_aos / t_f64)
end

# ---------------------------------------------------------------------------
# 2. Does the Dual path vectorize at all?
# ---------------------------------------------------------------------------

function irstats(f, types)
    io = IOBuffer()
    code_llvm(io, f, types; debuginfo=:none, optimize=true)
    s = String(take!(io))
    n(re) = count(_ -> true, eachmatch(re, s))
    return (v2=n(r"<2 x double>"), v4=n(r"<4 x double>"), v8=n(r"<8 x double>"))
end

println("\n2. LLVM vector-instruction counts in the observing-pass inner loops")
if Base.JLOptions().check_bounds == 1
    println("   skipped: bounds checking is forced on, which blocks the vectorizer")
else
    θ0 = example_theta(1)
    @printf("   %-16s %-10s %7s %7s %7s\n", "kernel", "eltype", "<2xd>", "<4xd>", "<8xd>")
    for (label, N) in (("Float64", 0), ("Dual{1}", 1), ("Dual{4}", 4),
                       ("Dual{8}", 8), ("Dual{12}", 12))
        θ = N == 0 ? θ0 : seed(θ0, Val(N))
        sys = build(θ, Val(1))
        traj = Trajectory{eltype(θ)}(sys, EPOCHS)
        orbitsolve!(traj, sys; method=METHOD)
        for (kname, fn) in (("_rotate_pass!", _rotate_pass!),
                            ("_observe_pass!", _observe_pass!))
            st = irstats(fn, (typeof(traj), Val{NB}))
            @printf("   %-16s %-10s %7d %7d %7d\n", kname, label, st.v2, st.v4, st.v8)
        end
    end
end

# ---------------------------------------------------------------------------
# 3. Gradient cost vs chunk width
# ---------------------------------------------------------------------------

println("\n3. Gradient cost vs ForwardDiff chunk width ($NEP epochs, absolute frame)")
@printf("   %-4s %-4s %6s %8s %11s %9s\n", "Np", "nθ", "chunk", "chunks", "grad [µs]", "vs best")
for NP in 1:3
    θ = example_theta(NP)
    np = Val(NP)
    f = θ -> evalonce(θ, np, METHOD)
    P = length(θ)
    dflt = ForwardDiff.pickchunksize(P)
    widths = filter(<=(P), unique(sort([1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24,
                                        P, dflt, cld(P, 2), cld(P, 3)])))
    results = map(widths) do N
        cfg = ForwardDiff.GradientConfig(f, θ, ForwardDiff.Chunk{N}())
        ForwardDiff.gradient(f, θ, cfg)
        (N, cld(P, N), 1e6 * @belapsed(ForwardDiff.gradient($f, $θ, $cfg), seconds = 1))
    end
    best = minimum(r[3] for r in results)
    for (N, nc, t) in results
        mark = N == dflt ? "  <- ForwardDiff default" : (t == best ? "  <- best" : "")
        @printf("   %-4d %-4d %6d %8d %11.1f %8.2fx%s\n", NP, P, N, nc, t, t / best, mark)
    end
end
