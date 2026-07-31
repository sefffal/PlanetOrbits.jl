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

With `simd=true` (the default), Float64 solves with the Markley/Auto solver
batch across epochs through branch-free vectorizable kernels (≈4× on AVX2;
agrees with the scalar solver to ≤4e-15). Other element types (e.g.
ForwardDiff Duals) and solvers always use the scalar path.
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
                    observing_geometry::Bool=true)
    _check_method(sys, method)
    traj = Trajectory(sys, epochs)
    return orbitsolve!(traj, sys; method, observing_geometry)
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
                    observing_geometry::Bool=true)
    traj = orbitsolve(sys, SVector(float(t)); method, observing_geometry)
    return traj[1]
end

"""
    orbitsolve!(traj::Trajectory, sys::System; method=KeplerianApprox())

Fill a caller-allocated `Trajectory` with per-body barycentric states of
`sys` at `traj.epochs`. Performs no allocation itself: with caller-provided
(e.g. bump-allocated) column storage the whole construct → solve → query
path is allocation-free.

`observing_geometry=false` skips the observing-geometry pass — see
`observe_pass!` and the note on `orbitsolve`.
"""
function orbitsolve!(traj::Trajectory, sys::System{NB,NR};
                     method::AbstractPropagator=KeplerianApprox(),
                     observing_geometry::Bool=true) where {NB,NR}
    size(traj.x, 2) == NB || throw(DimensionMismatch(
        "trajectory body storage has $(size(traj.x,2)) columns but the system has $NB bodies"))
    frame_pass!(traj, sys.frame)
    propagate!(traj, sys, method)
    if observing_geometry
        observe_pass!(traj, sys)
    else
        observe_skip!(traj, sys)
    end
    return traj
end

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

function frame_pass!(traj::Trajectory, fr::NoFrame)
    @inbounds traj.frame[1] = fr
    copyto!(traj.t_em, traj.epochs)
    return traj
end

function frame_pass!(traj::Trajectory, fr::Parallax)
    @inbounds traj.frame[1] = fr
    copyto!(traj.t_em, traj.epochs)
    return traj
end

function frame_pass!(traj::Trajectory, fr::AbsoluteFrame)
    epochs = traj.epochs
    # Recorded here, not at construction: the hot loop rebuilds `sys` from θ
    # every sample while reusing trajectory buffers, so anything captured
    # earlier would be a previous sample's frame. See `Trajectory`.
    @inbounds traj.frame[1] = fr
    @inbounds for k in eachindex(epochs)
        tobs = epochs[k]
        # Solve t_em = t_obs − (d(t_em) − d_ref)/c: seed with the linear
        # frame-RV estimate, then two fixed-point steps (the map contracts by
        # v_r/c ~ 1e-4 per step, and the seed is already good to the
        # perspective-acceleration term).
        ltt = fr.rv * (tobs - fr.ref_epoch) * 60 * 60 * 24 / c_light_ms
        t_em = tobs - ltt * sec2day
        t_em += tobs - _received_epoch(fr, _nudge_ref(fr, t_em))
        t_em += tobs - _received_epoch(fr, _nudge_ref(fr, t_em))
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
        else
            solve_row!(traj, row, j, method.solver)
        end
    end
    return traj
end

# The SIMD batch path applies only to Float64 storage/elements with the
# Markley (or Auto, which selects Markley for e < 1) solver; everything else
# — Duals, other solvers — is compile-time routed to the scalar path.
@inline _use_simd(::KeplerianApprox, ::Trajectory, ::Row) = false
@inline _use_simd(method::KeplerianApprox{<:Union{Auto,Markley}},
                  ::Trajectory{Float64}, row::Row{Float64}) =
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
