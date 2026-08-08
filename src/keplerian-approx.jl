# ---------------------------------------------------
# KeplerianApprox propagator
#
# Formalizes the physics PlanetOrbits v1 + Octofitter always used: each
# hierarchy row's Keplerian is solved independently, and absolute
# barycentric per-body states follow from the (mass-weighted) A⁻¹ transform.
# Three passes over the trajectory:
#   1. frame_pass!   — per-epoch frame compensation, shared by all bodies
#   2. solve_row!    — per-row Kepler solve, batched across epochs
#   3. combine!      — A⁻¹: relative row states → absolute body states
# ---------------------------------------------------

abstract type AbstractPropagator end

"""
    KeplerianApprox(; solver=Auto(), simd=true)

Propagator in which every hierarchy row evolves on an independent Keplerian
orbit (exact for two bodies; the classic approximation for hierarchical
systems). `solver` selects the Kepler-equation algorithm (see
`PlanetOrbits.Markley`, `PlanetOrbits.Goat`, `PlanetOrbits.RootsMethod`).

With `simd=true` (the default), solves with the Markley/Auto solver batch
across epochs through branch-free vectorizable kernels (≈4× on AVX2; agrees
with the scalar solver to ≤4e-15). This covers both Float64 elements and
first-order ForwardDiff `Dual`s, which solve their *primal* roots through the
same kernel and attach partials analytically — so a gradient evaluation
carries a value bit-identical to a plain Float64 evaluation. Other solvers,
nested `Dual`s (Hessians) and hyperbolic orbits use the scalar path, which
applies the same implicit rule one epoch at a time.
"""
struct KeplerianApprox{S<:AbstractSolver} <: AbstractPropagator
    solver::S
    simd::Bool
end
KeplerianApprox(; solver::AbstractSolver=Auto(), simd::Bool=true) =
    KeplerianApprox{typeof(solver)}(solver, simd)

"""
    orbitsolve(sys::System, epochs; method=KeplerianApprox())

Solve `sys` at the sorted `epochs` [MJD], returning a `Trajectory`. Index it
(`traj[k]`) for per-epoch solutions to pass to observables.

Allocates the trajectory storage; for allocation-free hot loops see
`orbitsolve!`.
"""
function orbitsolve(sys::System, epochs::AbstractVector;
                    method::AbstractPropagator=KeplerianApprox(),
                    observing_geometry::Bool=true,
                    barycentric_lighttime::Bool=true)
    _check_method(sys, method)
    traj = Trajectory(sys, epochs)
    return orbitsolve!(traj, sys; method, observing_geometry, barycentric_lighttime)
end

# Propagator-specific sanity checks (e.g. AHL21's h vs P_min guidance) run in
# the allocating convenience entry point; `orbitsolve!` stays logging-free.
_check_method(::System, ::AbstractPropagator) = nothing

"""
    orbitsolve(sys::System, t::Real; method=KeplerianApprox())

Single-epoch convenience: solve at one epoch [MJD] and return the solution
directly. Allocates; batch epochs into a vector for performance.
"""
function orbitsolve(sys::System, t::Real; method::AbstractPropagator=KeplerianApprox(),
                    observing_geometry::Bool=true,
                    barycentric_lighttime::Bool=true)
    traj = orbitsolve(sys, SVector(float(t)); method, observing_geometry,
                      barycentric_lighttime)
    return traj[1]
end

# Epochs per task below which a chunk is not worth its ~µs spawn cost: the
# per-epoch work across all passes is O(100 ns), so 512 epochs ≈ tens of µs
# per chunk and the orchestration stays a few percent.
const MIN_EPOCHS_PER_TASK = 512

"""
    orbitsolve!(traj::Trajectory, sys::System; method=KeplerianApprox())

Fill a caller-allocated `Trajectory` with per-body barycentric states of
`sys` at `traj.epochs`. Performs no allocation itself (unless `threads > 1`):
with caller-provided (e.g. bump-allocated) column storage the whole
construct → solve → query path is allocation-free.

The two precision opt-outs are independent, and gate *different* corrections —
see the "Precision opt-outs" page in the manual before setting either.

- `observing_geometry=false` skips the observing-geometry pass, whose terms all
  scale with the system's angular extent ρ. See `observe_pass!`.
- `barycentric_lighttime=false` skips the barycentric light-travel solve, a
  whole-system *timing* correction that scales with proximity and proper
  motion, not with ρ. See `frame_pass!`.

`threads=n` (default 1) splits the epochs into up to `n` contiguous chunks
solved on concurrent tasks. Every pass is epoch-local, each chunk writes a
disjoint epoch range of the same storage, and epochs keep their identity —
so the result is identical to the serial solve, bit for bit. Only the
`KeplerianApprox` propagator supports this (`AHL21` marches through time
sequentially); with any other method, or too few epochs for the task
overhead to amortize (`$MIN_EPOCHS_PER_TASK` per task), the solve silently
runs serial.
"""
function orbitsolve!(traj::Trajectory, sys::System{NB,NR};
                     method::AbstractPropagator=KeplerianApprox(),
                     observing_geometry::Bool=true,
                     barycentric_lighttime::Bool=true,
                     threads::Integer=1) where {NB,NR}
    size(traj.x, 2) == NB || throw(DimensionMismatch(
        "trajectory body storage has $(size(traj.x,2)) columns but the system has $NB bodies"))
    nchunks = _solve_chunks(method, nepochs(traj), threads)
    if nchunks > 1
        # The chunk views carry private frame holders (see `_epochview`), so
        # the shared frame column is written here, once.
        @inbounds traj.frame[1] = sys.frame
        # Chunk boundaries are aligned to the Dual solve's straight-line block
        # (see `DUAL_SIMD_BLOCK`): every chunk then tiles into blocks exactly
        # as the serial solve does, and only the final chunk carries the same
        # tail epochs the serial solve would. Misaligned boundaries would put
        # different epochs through the block vs. tail code, whose roundings
        # differ in the last ulp — and "identical to serial, bit for bit" is
        # the contract that makes `threads=` safe to flip on anywhere.
        len = DUAL_SIMD_BLOCK * cld(cld(nepochs(traj), nchunks), DUAL_SIMD_BLOCK)
        @sync for c in 1:nchunks
            lo = (c - 1) * len + 1
            hi = min(c * len, nepochs(traj))
            lo <= hi || continue
            sub = _epochview(traj, lo:hi)
            # Positional, not the keyword front door: a keyword call from a
            # task closure goes through the kwarg sorter on every spawn.
            Threads.@spawn _solve_serial!(sub, sys, method, observing_geometry,
                                          barycentric_lighttime)
        end
        return traj
    end
    return _solve_serial!(traj, sys, method, observing_geometry, barycentric_lighttime)
end

function _solve_serial!(traj::Trajectory, sys::System, method::AbstractPropagator,
                        observing_geometry::Bool, barycentric_lighttime::Bool)
    frame_pass!(traj, sys.frame, barycentric_lighttime)
    propagate!(traj, sys, method)
    if observing_geometry
        observe_pass!(traj, sys)
    else
        observe_skip!(traj, sys)
    end
    return traj
end

# Threading splits over epochs, so it is only available to propagators whose
# work is epoch-local. AHL21 is a time march and must see the epochs in order.
_solve_chunks(::KeplerianApprox, nep::Int, threads::Integer) =
    max(1, min(Int(threads), Threads.nthreads(), fld(nep, MIN_EPOCHS_PER_TASK)))
_solve_chunks(::AbstractPropagator, nep::Int, threads::Integer) = 1

# Propagator seam: the propagator owns everything between the two
# propagator-independent passes. Each method fills the per-body absolute
# barycentric state columns of the trajectory at traj.t_em, in the reference
# triad; `observe_pass!` then converts those into what is observed.
function propagate!(traj::Trajectory, sys::System, method::KeplerianApprox)
    solve_rows!(traj, sys, method)
    combine!(traj, sys)
    return traj
end

# ---------------------------------------------------
# Pass 1: frame compensation (once per system per epoch)
# ---------------------------------------------------

"""
    frame_pass!(traj, frame, barycentric_lighttime=true)

Fill the per-epoch frame columns: the emission epoch `t_em` and the
barycentre's propagated proper motion and radial velocity.

`barycentric_lighttime=false` takes `t_em = t_obs`, skipping the light-travel
solve. It does **not** skip the space-motion propagation — `pmra2`, `pmdec2`
and `rv2` are still propagated to each epoch, because that (and specifically
its perspective-acceleration curvature) *is* the absolute-astrometry signal.
The frame-level equivalent of "do not propagate at all" is not a solve-time
flag; it is choosing a `Parallax` frame when the `System` is built.
"""
function frame_pass!(traj::Trajectory, fr::NoFrame, ::Bool=true)
    @inbounds traj.frame[1] = fr
    copyto!(traj.t_em, traj.epochs)
    return traj
end

# No barycentric light-travel correction exists without a distance and a
# space velocity, so the flag is a no-op for both cheap frames.
function frame_pass!(traj::Trajectory, fr::Parallax, ::Bool=true)
    @inbounds traj.frame[1] = fr
    copyto!(traj.t_em, traj.epochs)
    return traj
end

function frame_pass!(traj::Trajectory, fr::AbsoluteFrame,
                     barycentric_lighttime::Bool=true)
    # Recorded here, not at construction: the hot loop rebuilds `sys` from θ
    # every sample while reusing trajectory buffers, so anything captured
    # earlier would be a previous sample's frame. See `Trajectory`.
    @inbounds traj.frame[1] = fr
    # Branch outside the loop so each body stays branch-free.
    if barycentric_lighttime
        _frame_pass_kernel!(traj, fr, true)
    else
        _frame_pass_kernel!(traj, fr, false)
    end
    return traj
end

@inline function _frame_pass_kernel!(traj::Trajectory, fr::AbsoluteFrame,
                                     lighttime::Bool)
    epochs = traj.epochs
    @inbounds for k in eachindex(epochs)
        tobs = epochs[k]
        t_em = tobs
        if lighttime
            # Solve t_em = t_obs − (d(t_em) − d_ref)/c: seed with the linear
            # frame-RV estimate, then two fixed-point steps (the map contracts
            # by v_r/c ~ 1e-4 per step, and the seed is already good to the
            # perspective-acceleration term).
            ltt = fr.rv * (tobs - fr.ref_epoch) * 60 * 60 * 24 / c_light_ms
            t_em = tobs - ltt * sec2day
            t_em += tobs - _received_epoch(fr, _nudge_ref(fr, t_em))
            t_em += tobs - _received_epoch(fr, _nudge_ref(fr, t_em))
        end
        # `ra2`/`dec2` are deliberately *not* stored: they are pure leaf
        # outputs (only `frame_ra`/`frame_dec` read them), and the only two
        # quantities here needing a transcendental — everything else is
        # algebra. They cannot merely be left unused, because `atand`/`asind`
        # carry `throw` paths that block dead-code elimination, so the split
        # is structural. 113.6 -> 39.6 ns/epoch across this and the
        # simplification in `_compensate_kinematics`. See design §10.3.
        kin = _compensate_kinematics(fr, t_em)
        traj.t_em[k] = t_em
        traj.pmra2[k] = kin.pmra2
        traj.pmdec2[k] = kin.pmdec2
        traj.rv2[k] = kin.rv2
    end
    return traj
end

# ---------------------------------------------------
# Pass 2: per-row Kepler solve, batched across epochs
# ---------------------------------------------------

function solve_rows!(traj::Trajectory, sys::System{NB,NR}, method::KeplerianApprox) where {NB,NR}
    for j in 1:NR
        row = sys.rows[j]
        if _use_simd(method, traj, row)
            solve_row_simd!(traj, row, j)
        elseif _use_dual_simd(method, traj, row)
            solve_row_dual_simd!(traj, row, j)
        else
            solve_row!(traj, row, j, method.solver)
        end
    end
    return traj
end

# The SIMD batch path applies only to Float64 storage/elements with the
# Markley (or Auto, which selects Markley for e < 1) solver; everything else
# — other solvers, nested Duals, hyperbolae — is compile-time routed to the
# scalar path, where `_anomaly_sincos` still applies the implicit rule.
@inline _use_simd(::KeplerianApprox, ::Trajectory, ::Row) = false
@inline _use_simd(method::KeplerianApprox{<:Union{Auto,Markley}},
                  ::Trajectory{Float64}, row::Row{Float64}) =
    method.simd && row.e < 1

# First-order ForwardDiff Duals over Float64 solve their primal roots through
# the same batch kernel (see `solve_row_dual_simd!`). The row may be plain
# Float64 — differentiating only the frame — but the trajectory must be Dual,
# since that is what carries `t_em` and the state columns.
@inline _use_dual_simd(::KeplerianApprox, ::Trajectory, ::Row) = false
@inline _use_dual_simd(method::KeplerianApprox{<:Union{Auto,Markley}},
                       ::Trajectory{Dual{Tg,Float64,N}},
                       row::Row{<:Union{Float64,Dual{Tg,Float64,N}}}) where {Tg,N} =
    method.simd && row.e < 1

# Solve the row's Kepler equation at mean anomaly `MA` and return the pair
# the state kernel consumes: (sin E, cos E) for ellipses, (sinh H, cosh H)
# for hyperbolae. `_states_from_E` is then identical for both conics — with
# a < 0 and sqrt1me2 = −√(e²−1), its algebra is the analytic continuation of
# the elliptical case.
@inline function _anomaly_sincos(row::Row, MA, solver::AbstractSolver)
    if row.hyperbolic
        H = kepler_solver(MA, row.e, HyperbolicHalley())
        return sinh(H), cosh(H)
    end
    E = kepler_solver(MA, row.e, solver)
    return sincos(E)
end

# Dual mean anomaly: solve on primal values and attach the partials with
# `_dual_sincosE` directly, instead of routing through `kepler_solver`'s Dual
# methods and then calling `sincos` on the Dual root they return. The two are
# bit-identical — the same primal solve, the same implicit rule — but the
# generic path evaluates sincos twice, once inside the rule and once on the
# root, and only the first is needed. Hyperbolic rows keep the generic path.
@inline function _anomaly_sincos(row::Row, MA::Dual{Tg,V,N},
                                 solver::AbstractSolver) where {Tg,V,N}
    if row.hyperbolic
        H = kepler_solver(MA, row.e, HyperbolicHalley())
        return sinh(H), cosh(H)
    end
    e = convert(Dual{Tg,V,N}, row.e)
    E = kepler_solver(value(MA), value(e), solver)
    sE, cE = sincos(E)
    return _dual_sincosE(MA, e, sE, cE)
end

# From sincos(E) to position/velocity, all algebraic: sin/cos of the true
# anomaly are derived from sincos(E) directly, so exactly one transcendental
# evaluation happens per solve.
@inline function _states_from_E(row::Row, sE, cE)
    temp = 1 - row.e * cE            # = 1 - e cos E = r/a
    r = row.a * temp
    invtemp = inv(temp)
    cosν = (cE - row.e) * invtemp
    sinν = row.sqrt1me2 * sE * invtemp
    cosν_ω = cosν * row.cosω - sinν * row.sinω
    sinν_ω = sinν * row.cosω + cosν * row.sinω
    x = r * (cosν_ω * row.sinΩ + sinν_ω * row.cosi_cosΩ)
    y = r * (cosν_ω * row.cosΩ - sinν_ω * row.cosi_sinΩ)
    z = r * (sinν_ω * row.sini)
    vfac1 = cosν_ω + row.ecosω
    vfac2 = sinν_ω + row.esinω
    vx = row.J * (row.cosi_cosΩ * vfac1 - row.sinΩ * vfac2)
    vy = -row.J * (row.cosi_sinΩ * vfac1 + row.cosΩ * vfac2)
    vz = row.J * row.sini * vfac1
    return x, y, z, vx, vy, vz
end

function solve_row!(traj::Trajectory, row::Row, j::Int, solver::AbstractSolver)
    n_per_day = row.n / year2day_julian
    t_em = traj.t_em
    rx = traj.rx; ry = traj.ry; rz = traj.rz
    rvx = traj.rvx; rvy = traj.rvy; rvz = traj.rvz
    @inbounds for k in eachindex(t_em)
        MA = n_per_day * (t_em[k] - row.tp)
        sE, cE = _anomaly_sincos(row, MA, solver)
        x, y, z, vx, vy, vz = _states_from_E(row, sE, cE)
        rx[k, j] = x; ry[k, j] = y; rz[k, j] = z
        rvx[k, j] = vx; rvy[k, j] = vy; rvz[k, j] = vz
    end
    return traj
end

# ---------------------------------------------------
# Pass 3: A⁻¹ combine — relative row states → absolute body states
# ---------------------------------------------------

function combine!(traj::Trajectory, sys::System{NB,NR}) where {NB,NR}
    Ainv = sys.Ainv
    _combine_component!(traj.x, traj.rx, Ainv)
    _combine_component!(traj.y, traj.ry, Ainv)
    _combine_component!(traj.z, traj.rz, Ainv)
    _combine_component!(traj.vx, traj.rvx, Ainv)
    _combine_component!(traj.vy, traj.rvy, Ainv)
    _combine_component!(traj.vz, traj.rvz, Ainv)
    return traj
end

function _combine_component!(dst::AbstractMatrix, rel::AbstractMatrix,
                             Ainv::SMatrix{NB,NR}) where {NB,NR}
    @inbounds for k in axes(dst, 1)
        # NB and NR are compile-time constants: these loops unroll.
        for j in 1:NB
            acc = zero(eltype(dst))
            for r in 1:NR
                acc = muladd(Ainv[j, r], rel[k, r], acc)
            end
            dst[k, j] = acc
        end
    end
    return dst
end

export orbitsolve, orbitsolve!, KeplerianApprox, Trajectory
