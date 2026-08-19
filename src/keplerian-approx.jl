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
                     threads::Integer=1,
                     row_cache=nothing) where {NB,NR}
    size(traj.x, 2) == NB || throw(DimensionMismatch(
        "trajectory body storage has $(size(traj.x,2)) columns but the system has $NB bodies"))
    nchunks = _solve_chunks(method, nepochs(traj), threads)
    if nchunks > 1
        # The chunk views carry private frame holders (see `_epochview`), so
        # the shared frame column is written here, once — the same effective
        # frame every chunk's `frame_pass!` will compute for itself.
        @inbounds traj.frame[1] = _effective_frame(sys.frame, barycentric_lighttime)
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
            # `row_cache` is deliberately NOT forwarded: the chunks would
            # share (and race on) one cache, and their view columns are fresh
            # objects every call so it could never hit anyway.
            Threads.@spawn _solve_serial!(sub, sys, method, observing_geometry,
                                          barycentric_lighttime, nothing)
        end
        return traj
    end
    return _solve_serial!(traj, sys, method, observing_geometry,
                          barycentric_lighttime, row_cache)
end

function _solve_serial!(traj::Trajectory, sys::System, method::AbstractPropagator,
                        observing_geometry::Bool, barycentric_lighttime::Bool,
                        row_cache=nothing)
    frame_pass!(traj, sys.frame, barycentric_lighttime)
    propagate!(traj, sys, method, row_cache)
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
#
# `row_cache` is an optional solver-state reuse hint (see the caching block
# below); propagators that cannot honor it ignore it, so the 4-arg form is
# always safe to call.
propagate!(traj::Trajectory, sys::System, method::AbstractPropagator, row_cache) =
    propagate!(traj, sys, method)
function propagate!(traj::Trajectory, sys::System, method::KeplerianApprox,
                    row_cache=nothing)
    solve_rows!(traj, sys, method, row_cache)
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

The two settings are two *conventions*, both anchored to the catalog values at
`ref_epoch`, and the frame readouts are internally consistent under each:

- **off** — the light-time-free standard model the astrometric catalogs are
  reduced with (ESA 1997 Sect. 1.5.5; Lindegren et al. 2012/2021): catalog
  values propagated linearly in 3D, readouts are that worldline's direction
  and rates at `t_obs`. Use this for catalog-convention absolute astrometry
  (Butkevich & Lindegren 2014, Sects. 5.5, 6.1).
- **on** — the rigorous apparent path: the catalog (apparent) proper motions
  are first de-Dopplered to the true space velocity (`_dedoppler`), the
  emission epoch is solved along that worldline, and the position *and* rate
  readouts are apparent quantities, `d/dt_obs`. Both conventions reproduce
  the catalog values exactly at `ref_epoch`; away from it they differ only by
  the genuine (second-order) light-time terms.
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
    # The light-time solve runs on the de-Dopplered *true* worldline; the
    # light-time-free path propagates the catalog values as they stand. See
    # `_dedoppler` for why using the catalog velocity under the light-time
    # solve would double-count the μ·v_r/c Doppler factor.
    fr = _effective_frame(fr, barycentric_lighttime)
    # Recorded here, not at construction: the hot loop rebuilds `sys` from θ
    # every sample while reusing trajectory buffers, so anything captured
    # earlier would be a previous sample's frame. See `Trajectory`. What is
    # recorded is the *effective* frame, so the on-demand `frame_ra`/
    # `frame_dec` reads propagate the same worldline the solve did.
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
        if lighttime
            # `_compensate_kinematics` returns coordinate rates d/dt_em along
            # the worldline. The *angular* readouts pair with the observation
            # epoch, so convert them to apparent rates:
            #     d/dt_obs = (d/dt_em) / (dt_obs/dt_em),  dt_obs/dt_em = 1 + ḋ/c
            # (B&L 2014 Eqs. 10–13). At `ref_epoch` this exactly undoes the
            # `_dedoppler` factor, reproducing the catalog proper motions
            # identically; away from it the two conventions differ only by the
            # genuine light-time terms. Without it the PM readouts would be
            # coordinate rates while the position readouts slew at apparent
            # rates — a μ·v_r/c inter-readout inconsistency (3.8 mas/yr at
            # Barnard's star).
            #
            # `rv2` deliberately does NOT take this factor. The angular
            # convention is forced on us: catalog proper motions *are*
            # apparent, and they must pair with positions read at t_obs. The
            # radial one is not, and taking it would answer the wrong
            # question. A spectrograph measures the Doppler shift of light
            # emitted at t_em, which is set by the source's velocity at that
            # event — the coordinate rate ḋ(t_em) — not by the rate at which
            # the light-time-affected distance changes against arrival time.
            # (Lindegren & Dravins 2003 draw the same astrometric-vs-
            # spectroscopic distinction.) Gaia itself publishes apparent
            # proper motions beside a spectroscopic radial velocity, so this
            # mixed convention is the data's own — and it keeps
            # `frame_rv(ref_epoch) == rv` exactly under both settings, which
            # is what every consumer composing `frame_rv` against a catalog
            # radial velocity relies on.
            s = inv(1 + kin.rv2 / c_light_ms)
            traj.pmra2[k] = kin.pmra2 * s
            traj.pmdec2[k] = kin.pmdec2 * s
            traj.rv2[k] = kin.rv2
        else
            traj.pmra2[k] = kin.pmra2
            traj.pmdec2[k] = kin.pmdec2
            traj.rv2[k] = kin.rv2
        end
    end
    return traj
end

# ---------------------------------------------------
# Pass 2: per-row Kepler solve, batched across epochs
# ---------------------------------------------------

function solve_rows!(traj::Trajectory, sys::System{NB,NR}, method::KeplerianApprox,
                     row_cache::Nothing=nothing) where {NB,NR}
    tspan = _solve_tspan(traj, method)
    for j in 1:NR
        _solve_one_row!(traj, sys.rows[j], j, method, tspan)
    end
    return traj
end

@inline function _solve_one_row!(traj::Trajectory, row::Row, j::Int,
                                 method::KeplerianApprox,
                                 tspan::Tuple{Float64,Float64})
    inrange = _ma_in_kernel_range(row, tspan)
    if _use_simd(method, traj, row, inrange)
        solve_row_simd!(traj, row, j)
    elseif _use_dual_simd(method, traj, row, inrange)
        solve_row_dual_simd!(traj, row, j)
    else
        solve_row!(traj, row, j, method.solver)
    end
    return traj
end

# ---------------------------------------------------
# Per-row solve caching (primal Float64 path only)
#
# A coordinate-wise sampler — Pigeons' SliceSampler above all — evaluates the
# log density many times in a row with *all but one* parameter bit-identical.
# Each hierarchy row's Kepler solve reads only that row's elements and `t_em`,
# so a row whose inputs are unchanged since the previous solve into the same
# storage would recompute exactly the states already sitting in its columns.
# At ~90 epochs the Markley batch is the single largest term in the primal
# evaluation, so skipping it when a nuisance parameter (jitter, an offset, a
# proper motion, plx) — or another row's element — is the thing being varied
# is worth 15–40% depending on how many rows the model has.
#
# Correctness rests on three things, all checked per call:
#
#   1. *Storage identity.* The cache remembers the `rx` column and `epochs`
#      objects and compares with `===`. For heap `Array`s that is object
#      identity, immune to address recycling. For Bumper-backed columns
#      (`UnsafeArray`, an isbits struct) egal compares pointer+dims — and
#      there address reuse is exactly the design: Octofitter's generated
#      evaluation carves the same columns at the same offsets out of the same
#      task-local slab every evaluation, so pointer+dims equality means "the
#      bytes my previous solve wrote are still there". Two *different* live
#      models never collide because their epoch vectors are distinct objects.
#
#   2. *Row inputs.* The key holds the row's seven defining elements. `==` on
#      the tuple makes NaN a guaranteed miss, and −0.0 == 0.0 folds two
#      parameterizations of the same orbit together, both conservative-or-exact.
#
#   3. *`t_em` identity.* `t_em` is a deterministic function of (epochs,
#      frame, barycentric_lighttime). The key carries the frame's seven
#      catalog scalars plus three samples of `t_em` itself (first, middle,
#      last). With the frame scalars pinned, the on/off lighttime difference
#      is (d(t_em) − d_ref)/c whose zero set along a linear space motion has
#      at most two roots — so three samples cannot all coincide and a flag
#      flip is a guaranteed miss without threading the flag through
#      `propagate!`'s seam.
#
# The skip is bit-exact, not approximate: a hit re-uses values the identical
# computation already wrote. The one behavior change is contractual: a caller
# who manually scribbles into `traj.rx` between two solves of an *identical*
# system no longer gets the scribbles overwritten — which is why the cache is
# strictly opt-in, per call: pass `row_cache = PlanetOrbits._rowcache(traj)`
# to `orbitsolve!`. The default (`nothing`) keeps re-fill semantics and keeps
# the solve *statically* allocation-free — `_rowcache` touches task-local
# storage, whose first use on a task allocates its IdDict, and an opt-in at
# the call boundary is the only place that cost can live without failing the
# AllocCheck gate on the default path. (Octofitter opts in from its generated
# likelihood, which is where the sampler's evaluation loop lives.)
#
# NB for inert-row elision (indicator-off companions): a row skipped for
# *that* reason must write `_NO_ROW_KEY` into its cache slot, not a valid key
# — its columns then hold zeros, not the key's states.
# ---------------------------------------------------

const _ROW_KEY_LEN = 18
const _NO_ROW_KEY = ntuple(_ -> NaN, _ROW_KEY_LEN)

# The method changes the bits a solve writes (SIMD and scalar Markley agree to
# 4e-15, not exactly; other solvers more so), so it is part of the key: one
# slot encoding the solver type and the simd flag.
@inline _rowcache_methodkey(method::KeplerianApprox) =
    Float64(hash(typeof(method.solver)) & 0x000fffff) + (method.simd ? 0.5 : 0.0)

# Parametric on the trajectory's column and epoch types so the identity
# compares below are unboxed: an `Any` field `===` an isbits `UnsafeArray`
# (the Bumper-backed column type) would box the operand on every likelihood
# evaluation — 32 bytes per call on the path whose contract is zero.
mutable struct _RowSolveCache{M,E}
    rx::Union{M,Nothing}      # column object of the trajectory last solved
    epochs::Union{E,Nothing}  # its epochs object
    keys::Vector{NTuple{_ROW_KEY_LEN,Float64}}
    hits::Int                 # per-task counters, for tests and diagnostics
    misses::Int
    _RowSolveCache{M,E}() where {M,E} =
        new{M,E}(nothing, nothing, NTuple{_ROW_KEY_LEN,Float64}[], 0, 0)
end

# One cache per task (and per column type), alongside Bumper's own task-local
# buffer: same storage discipline, same thread-safety argument. The type
# object is the key (interned, so `===` across calls). First touch on a task
# allocates the task's TLS dict — call this OUTSIDE any statically-checked
# zero-alloc region and pass the result in as `row_cache`.
@inline _rowcache(::Type{M}, ::Type{E}) where {M,E} =
    get!(() -> _RowSolveCache{M,E}(), task_local_storage(),
        _RowSolveCache{M,E})::_RowSolveCache{M,E}
@inline _rowcache(traj::Trajectory) =
    _rowcache(typeof(traj.rx), typeof(traj.epochs))

@inline _rowcache_framekey(::NoFrame) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# plx does not enter the row solve (t_em ≡ epochs on this path).
@inline _rowcache_framekey(::Parallax) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _rowcache_framekey(fr::AbsoluteFrame) =
    (fr.ra, fr.dec, fr.plx, fr.pmra, fr.pmdec, fr.rv, fr.ref_epoch)

function solve_rows!(traj::Trajectory, sys::System{NB,NR},
                     method::KeplerianApprox,
                     cache::_RowSolveCache{M,E}) where {NB,NR,M,E}
    nep = nepochs(traj)
    # The cache only applies to the primal Float64 path with matching column
    # types; anything else — Duals above all — quietly takes the plain loop.
    # (`isa` on types the compiler knows, so the branch folds.)
    if !(traj isa Trajectory{Float64} && sys isa System{NB,NR,Float64} &&
         traj.rx isa M && traj.epochs isa E && NR > 0 && nep > 0)
        return solve_rows!(traj, sys, method)
    end
    if !(cache.rx === traj.rx && cache.epochs === traj.epochs &&
         length(cache.keys) == NR)
        cache.rx = traj.rx
        cache.epochs = traj.epochs
        resize!(cache.keys, NR)
        fill!(cache.keys, _NO_ROW_KEY)
    end
    fk = _rowcache_framekey(sys.frame)
    tk = @inbounds (traj.t_em[1], traj.t_em[(nep + 1) >> 1], traj.t_em[nep])
    mk = _rowcache_methodkey(method)
    # Whether a row falls back from the batch kernel to the scalar solver is a
    # function of the row's elements and the epochs — both already pinned by
    # the key — so it needs no slot of its own, on exactly the argument §3
    # makes for `t_em` itself.
    tspan = _solve_tspan(traj, method)
    for j in 1:NR
        row = sys.rows[j]
        key = (row.a, row.e, row.i, row.ω, row.Ω, row.tp, row.M, fk..., tk..., mk)
        if @inbounds(cache.keys[j]) == key
            cache.hits += 1
            continue
        end
        cache.misses += 1
        _solve_one_row!(traj, row, j, method, tspan)
        @inbounds cache.keys[j] = key
    end
    return traj
end

# The SIMD batch path applies only to Float64 storage/elements with the
# Markley (or Auto, which selects Markley for e < 1) solver; everything else
# — other solvers, nested Duals, hyperbolae — is compile-time routed to the
# scalar path, where `_anomaly_sincos` still applies the implicit rule.
#
# Two runtime conditions ride along with the type-level ones:
#
#   * `!row.hyperbolic` rather than a fresh `row.e < 1`. The conic was
#     classified once, on the primal, in `Row`; re-deriving it here would let
#     a `Dual` eccentricity be classified lexicographically and route a row
#     the constructor called hyperbolic into the elliptical kernel.
#
#   * `inrange`, the mean-anomaly bound described at `VREM2PI_MAX`. The batch
#     kernel's angle reduction is a branch-free Cody–Waite step, valid only
#     while the quotient M/2π is an exactly-representable integer; the scalar
#     path's `Base.rem2pi` is Payne–Hanek and exact for every finite argument.
#     Routing the out-of-range rows to the scalar path is what keeps the two
#     paths answering the same question — and it is a *row*-level decision, so
#     the hot loop stays branch-free and the serial/threaded bit-identity
#     contract (which needs the choice to be a property of the row, not of a
#     chunk) is preserved.
@inline _use_simd(::KeplerianApprox, ::Trajectory, ::Row, inrange::Bool) = false
@inline _use_simd(method::KeplerianApprox{<:Union{Auto,Markley}},
                  ::Trajectory{Float64}, row::Row{Float64}, inrange::Bool) =
    method.simd && !row.hyperbolic && inrange

# First-order ForwardDiff Duals over Float64 solve their primal roots through
# the same batch kernel (see `solve_row_dual_simd!`). The row may be plain
# Float64 — differentiating only the frame — but the trajectory must be Dual,
# since that is what carries `t_em` and the state columns.
@inline _use_dual_simd(::KeplerianApprox, ::Trajectory, ::Row, inrange::Bool) = false
@inline _use_dual_simd(method::KeplerianApprox{<:Union{Auto,Markley}},
                       ::Trajectory{Dual{Tg,Float64,N}},
                       row::Row{<:Union{Float64,Dual{Tg,Float64,N}}},
                       inrange::Bool) where {Tg,N} =
    method.simd && !row.hyperbolic && inrange

# Three-argument forms: "will this row take the batch path?", answered from
# the trajectory alone. `_solve_one_row!` passes the range flag explicitly
# because it hoists the epoch reduction out of the row loop; everything else
# (tests, and anyone reasoning about routing at the REPL) wants the question
# without having to reproduce that hoist.
#
# PRECONDITION: `frame_pass!` has run. `t_em` is `undef` on a freshly built
# `Trajectory` — it is written by the frame pass, which `_solve_serial!` always
# runs before `solve_rows!`, so the pipeline satisfies this by construction.
# Asking earlier reads uninitialized memory; the answer degrades to "scalar
# path" rather than to anything unsafe (see `_ma_in_kernel_range`'s `≤`), but
# it is not the answer the solve will actually give.
@inline _use_simd(method::KeplerianApprox, traj::Trajectory, row::Row) =
    _use_simd(method, traj, row, _ma_in_kernel_range(row, _solve_tspan(traj, method)))
@inline _use_dual_simd(method::KeplerianApprox, traj::Trajectory, row::Row) =
    _use_dual_simd(method, traj, row, _ma_in_kernel_range(row, _solve_tspan(traj, method)))

"""
    PlanetOrbits._ma_in_kernel_range(row, tspan) -> Bool

Whether every mean anomaly this row will be evaluated at lies inside the
batch kernel's validated reduction range (`VREM2PI_MAX`).

`tspan` is `(tmin, tmax)` over the trajectory's emission epochs, so the bound
is `|n_per_day| · max(|tmin − tp|, |tmax − tp|)` — O(1) per row, off one
O(nepochs) reduction per solve, rather than a per-epoch test inside the loop.

Written as `≤` so that a `NaN` bound (a `NaN` epoch, or a `NaN` `tp` that
somehow reached a hand-built `Row`) answers `false` and takes the scalar path,
which propagates the `NaN` honestly instead of feeding it to a kernel whose
contract excludes it.
"""
@inline function _ma_in_kernel_range(row::Row, tspan::Tuple{Float64,Float64})
    tmin, tmax = tspan
    tp = _primal(row.tp)
    n_per_day = _primal(row.n) / year2day_julian
    maxdt = max(abs(tmin - tp), abs(tmax - tp))
    return abs(n_per_day) * maxdt <= VREM2PI_MAX
end

# Primal extrema of the emission epochs. One pass per *solve*, hoisted out of
# the row loop rather than repeated per row; `_primal` because the Dual path's
# `t_em` carries partials and the reduction range is a statement about values.
#
# Measured at n = 2000 epochs: 1.27 µs against 36 µs for `solve_rows!` and
# 219 µs for the whole `orbitsolve!` — 3.5% of the row solve, 0.6% of a call.
# `min`/`max` rather than a vectorizable `ifelse` pair on purpose: they
# propagate `NaN`, which is what makes a `NaN` epoch route to the scalar path
# instead of into a kernel whose contract excludes it.
@inline function _t_em_span(traj::Trajectory)
    t_em = traj.t_em
    isempty(t_em) && return (0.0, 0.0)
    tmin = tmax = _primal(@inbounds t_em[1])
    @inbounds for k in eachindex(t_em)
        t = _primal(t_em[k])
        tmin = min(tmin, t)
        tmax = max(tmax, t)
    end
    return (tmin, tmax)
end

# …and not even that when the batch kernel is switched off: the span exists
# only to decide whether the kernel's reduction domain is respected. `NaN`
# rather than a sentinel that reads as "in range" — every routing predicate
# then answers `false` for the same reason a `NaN` epoch does, so turning
# `simd` off cannot accidentally turn the kernel back on.
@inline _solve_tspan(traj::Trajectory, method::KeplerianApprox) =
    method.simd ? _t_em_span(traj) : (NaN, NaN)

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
