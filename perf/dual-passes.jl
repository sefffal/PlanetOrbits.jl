# Per-pass cost of the *gradient* path — the decomposition behind §10.4.1 of
# design/planetorbits-v2-nbody-migration.md.
#
# Run:
#   julia --project=perf perf/dual-passes.jl
#
# Measures each pass of `orbitsolve!` at Float64 and at ForwardDiff chunk
# widths 1/4/8/12, on the Gaia-DR4-shaped workload (150 epochs, Np=1, absolute
# frame). The point is the *penalty* column: a pass that only pays for the
# extra partial arithmetic should cost ~(N+1)/2-ish more than Float64; a pass
# that also loses SIMD across epochs costs much more than that, and that excess
# is what is recoverable by restructuring.

include("workload.jl")

using BenchmarkTools
using ForwardDiff
using ForwardDiff: Dual, Partials
using Printf

using PlanetOrbits: frame_pass!, propagate!, observe_pass!,
                    _frame_geometry_pass!, _rotate_pass!, _accjerk_pass!,
                    _observe_pass!, c_au_per_julianyr, rad2mas,
                    solve_row_dual_simd!, solve_row!, solve_rows!,
                    DUAL_SIMD_BLOCK, Auto

const NP = 1
const NPV = Val(NP)
const METHOD = KeplerianApprox(simd=true)

# Seed θ the way a chunked ForwardDiff.gradient does: every entry is a
# Dual{_,Float64,N}, with one-hot partials cycling through the chunk slots.
# Partial *values* never affect timing (nothing branches on them).
function seed(θ::Vector{Float64}, ::Val{N}) where {N}
    return [Dual{Nothing}(θ[i], Partials(ntuple(j -> Float64(j == mod1(i, N)), Val(N))))
            for i in eachindex(θ)]
end

# --- ablation kernels, for attributing the observe-pass penalty -------------

# `_observe_pass!` with its one sqrt and two divides replaced by multiplies.
# Everything else — the loads, stores, retardation polynomial — is identical,
# so the difference is exactly what the scalar transcendental units cost.
function _observe_pass_nodiv!(traj, ::Val{NB}) where {NB}
    T = eltype(traj.x)
    invc = inv(T(c_au_per_julianyr))
    r2m = T(rad2mas)
    nk = length(traj.epochs)
    @inbounds for j in 1:NB
        @simd for k in 1:nk
            x = traj.x[k, j]; y = traj.y[k, j]; z = traj.z[k, j]
            vx = traj.vx[k, j]; vy = traj.vy[k, j]; vz = traj.vz[k, j]
            ax = traj.ax[k, j]; ay = traj.ay[k, j]; az = traj.az[k, j]
            jx = traj.jx[k, j]; jy = traj.jy[k, j]; jz = traj.jz[k, j]
            Vx = traj.bvx[k]; Vy = traj.bvy[k]; Vz = traj.bvz[k]
            d_au = traj.d_au[k]
            dt = -z * invc
            dt = -(z + ((vz + Vz) * dt + az * dt * dt / 2)) * invc
            dt2 = dt * dt
            px = x + (vx + Vx) * dt + ax * dt2 / 2 + jx * dt2 * dt / 6
            py = y + (vy + Vy) * dt + ay * dt2 / 2 + jy * dt2 * dt / 6
            pz = z + (vz + Vz) * dt + az * dt2 / 2 + jz * dt2 * dt / 6
            nvx = vx + ax * dt + jx * dt2 / 2
            nvy = vy + ay * dt + jy * dt2 / 2
            nvz = vz + az * dt + jz * dt2 / 2
            dzp = d_au + pz
            invr = px * px + py * py + dzp * dzp        # was inv(sqrt(·))
            traj.x[k, j] = px; traj.y[k, j] = py; traj.z[k, j] = pz
            traj.vx[k, j] = nvx; traj.vy[k, j] = nvy
            traj.vz[k, j] = ((Vx + nvx) * px + (Vy + nvy) * py +
                             (Vz + nvz) * dzp) * invr - Vz
            traj.cart2angle[k, j] = r2m * dzp                   # was r2m / dzp
        end
    end
    return traj
end

# --- driver ----------------------------------------------------------------

struct Case{T}
    label::String
    θ::Vector{T}
    sys::Any
    traj::Any
end

function makecase(label, θraw)
    θ = θraw
    sys = build(θ, NPV)
    traj = Trajectory{eltype(θ)}(sys, EPOCHS)
    orbitsolve!(traj, sys; method=METHOD)
    return Case(label, θ, sys, traj)
end

θ0 = example_theta(NP)
cases = [makecase("Float64", θ0),
         makecase("Dual{1}",  seed(θ0, Val(1))),
         makecase("Dual{4}",  seed(θ0, Val(4))),
         makecase("Dual{8}",  seed(θ0, Val(8))),
         makecase("Dual{12}", seed(θ0, Val(12)))]

const NB = 2   # Np=1 -> star + companion

# Each entry: (name, closure taking a Case -> elapsed seconds). `setup=` re-runs
# the upstream passes before every sample, because the observing pass mutates
# the state columns in place and is not idempotent.
function bench_passes(c::Case)
    sys = c.sys; traj = c.traj; fr = sys.frame
    t_frame = @belapsed frame_pass!($traj, $fr, true)
    t_prop  = @belapsed propagate!($traj, $sys, $METHOD)
    t_obs   = @belapsed(observe_pass!($traj, $sys),
                        setup = (propagate!($traj, $sys, $METHOD)), evals = 1)
    t_read  = @belapsed readout($traj, $sys, $(eltype(c.θ)))
    t_total = @belapsed orbitsolve!($traj, $sys; method=$METHOD)
    return (; t_frame, t_prop, t_obs, t_read, t_total)
end

function bench_subpasses(c::Case)
    sys = c.sys; traj = c.traj; fr = sys.frame
    pre() = (propagate!(traj, sys, METHOD); nothing)
    t_geom = @belapsed(_frame_geometry_pass!($traj, $fr), setup = ($pre()), evals = 1)
    t_rot  = @belapsed(_rotate_pass!($traj, $(Val(NB))), setup = ($pre()), evals = 1)
    t_acc  = @belapsed(_accjerk_pass!($traj, $sys, $(Val(NB))), setup = ($pre()), evals = 1)
    t_obs  = @belapsed(_observe_pass!($traj, $(Val(NB))), setup = ($pre()), evals = 1)
    t_nod  = @belapsed(_observe_pass_nodiv!($traj, $(Val(NB))), setup = ($pre()), evals = 1)
    return (; t_geom, t_rot, t_acc, t_obs, t_nod)
end

μ(t) = 1e6 * t

println("Gaia-DR4-shaped gradient decomposition: $NEP epochs, Np=$NP, absolute frame")
println("(design §10.4.1)\n")

r = [bench_passes(c) for c in cases]
@printf("%-14s %8s %8s %8s %8s %8s | %8s\n",
        "pass", "Float64", "Dual{1}", "Dual{4}", "Dual{8}", "Dual{12}", "penalty")
for (name, f) in (("frame_pass!", :t_frame), ("propagate!", :t_prop),
                  ("observe_pass!", :t_obs), ("readout", :t_read),
                  ("total", :t_total))
    v = [getfield(x, f) for x in r]
    @printf("%-14s %8.2f %8.2f %8.2f %8.2f %8.2f | %7.1fx\n",
            name, μ(v[1]), μ(v[2]), μ(v[3]), μ(v[4]), μ(v[5]), v[5] / v[1])
end

println("\nobserve_pass! sub-passes (µs):\n")
s = [bench_subpasses(c) for c in cases]
@printf("%-24s %8s %8s %8s %8s %8s | %8s\n",
        "sub-pass", "Float64", "Dual{1}", "Dual{4}", "Dual{8}", "Dual{12}", "penalty")
for (name, f) in (("A1 _frame_geometry", :t_geom), ("A2 _rotate", :t_rot),
                  ("B  _accjerk", :t_acc), ("C  _observe", :t_obs),
                  ("C' _observe (no sqrt/div)", :t_nod))
    v = [getfield(x, f) for x in s]
    @printf("%-24s %8.2f %8.2f %8.2f %8.2f %8.2f | %7.1fx\n",
            name, μ(v[1]), μ(v[2]), μ(v[3]), μ(v[4]), μ(v[5]), v[5] / v[1])
end

# --- straight-line block width for the Dual Kepler batch --------------------
#
# `solve_row_dual_simd!` unrolls B independent primal solves so LLVM's SLP pass
# can pack them; B trades lane occupancy against the ~20 live values each
# `markley_sincosE` carries. Sweep it against the scalar row solve.

println("\nDual Kepler batch: straight-line block width (one row, µs)")
println("(default DUAL_SIMD_BLOCK = $DUAL_SIMD_BLOCK)\n")
@printf("%-14s %8s %8s %8s %8s\n", "block", "Dual{1}", "Dual{4}", "Dual{8}", "Dual{12}")
let dualcases = cases[2:end]
    row_of(c) = c.sys.rows[1]
    times = Dict{Any,Vector{Float64}}()
    for B in (1, 2, 4, 8, 16)
        times[B] = [@belapsed(solve_row_dual_simd!($(c.traj), $(row_of(c)), 1, $(Val(B))))
                    for c in dualcases]
    end
    times[:scalar] = [@belapsed(solve_row!($(c.traj), $(row_of(c)), 1, $(Auto())))
                      for c in dualcases]
    v = times[:scalar]
    @printf("%-14s %8.2f %8.2f %8.2f %8.2f\n", "scalar", μ.(v)...)
    for B in (1, 2, 4, 8, 16)
        v = times[B]
        @printf("%-14s %8.2f %8.2f %8.2f %8.2f\n", "B=$B", μ.(v)...)
    end
end
