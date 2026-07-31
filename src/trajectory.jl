# ---------------------------------------------------
# Trajectory: struct-of-arrays solution storage
#
# Epoch-fastest layout: every column is contiguous over epochs, so the
# per-row Kepler solve batches (and later SIMD-vectorizes) across epochs,
# and per-chunk epoch ranges are contiguous for threading.
#
# The frame *type* FM is carried in the trajectory's type so observables can
# dispatch on what frame information is available; the per-epoch frame data
# lives in the columns. The system's body-name table `Names` is carried too,
# so observables can resolve `Body` values and `Symbol`s by name at zero
# cost (see `_resolve` in body.jl).
# ---------------------------------------------------

struct Trajectory{T<:Number,FR<:AbstractFrame,Names,E<:AbstractVector{<:Real},
                  V<:AbstractVector{T},M<:AbstractMatrix{T},
                  FV<:AbstractVector{FR}}
    epochs::E    # observation epochs [MJD], sorted
    # The frame this trajectory was last solved against, as a one-element
    # column. It is written by `frame_pass!` rather than captured at
    # construction, because `orbitsolve!(traj, sys)` takes both and the hot
    # loop rebuilds `sys` from θ every sample while reusing trajectory
    # buffers — a construction-time capture would silently be the *previous*
    # sample's frame. Observables that need the frame (`frame_ra`,
    # `frame_dec`) read it from here, so they cannot go stale.
    frame::FV
    # per-epoch frame columns (shared by every body and observable)
    t_em::V          # light-travel-time–corrected emission epoch [MJD]
    pmra2::V; pmdec2::V  # propagated frame proper motion [mas/yr]
    rv2::V           # propagated frame radial velocity [m/s]
    # Observing geometry, per epoch: barycentre distance and space velocity
    # resolved in that epoch's triad. Hoisted out of the per-epoch loop so the
    # per-body passes can run with epochs innermost (see observe.jl).
    d_au::V          # barycentre distance [AU]
    bvx::V; bvy::V; bvz::V   # barycentre velocity [AU / julian year]
    # Reference-triad -> epoch-triad rotation, unpacked into scalar columns so
    # the rotation pass can run with epochs innermost and vectorize. Gated
    # against the readable `_geometry` reference in the test suite.
    R11::V; R12::V; R13::V
    R21::V; R22::V; R23::V
    R31::V; R32::V; R33::V
    # AU -> mas conversion, **per body**: rad2mas / (d_barycentre + z_body).
    # Bodies at different depths along the line of sight do not share a scale
    # factor; the difference is ρ² in radians (≈ 4.85·ρ[″]² µas), which is
    # 43 µas for a companion at 3″ and 0.05 µas at 100 mas.
    cart2angle::M
    # per-body absolute barycentric states [AU, AU/julian year], expressed in
    # the local triad of the barycentre's *apparent* direction at the epoch,
    # each body taken at its own light-travel-retarded time.
    x::M; y::M; z::M; vx::M; vy::M; vz::M
    # per-body acceleration and jerk, scratch for the observing pass: needed
    # for every body before any body is retarded, so it cannot be fused away
    ax::M; ay::M; az::M; jx::M; jy::M; jz::M
    # per-row Jacobi-relative scratch
    rx::M; ry::M; rz::M; rvx::M; rvy::M; rvz::M
end

"""
    Trajectory(sys, epochs)
    Trajectory{T}(sys, epochs)

Allocate storage for solving `sys` at `epochs` (sorted, MJD). The element
type defaults to the system's scalar type; pass `T` explicitly when solving
with a different element type (e.g. ForwardDiff Duals).

For allocation-free hot loops, the columns can instead be caller-provided
(e.g. Bumper-allocated) using the full inner constructor; see `orbitsolve!`.
"""
Trajectory(sys::System{NB,NR,T}, epochs::AbstractVector) where {NB,NR,T} =
    Trajectory{T}(sys, epochs)
function Trajectory{T}(sys::System{NB,NR}, epochs::AbstractVector) where {T,NB,NR}
    nep = length(epochs)
    vk() = Vector{T}(undef, nep)
    mkb() = Matrix{T}(undef, nep, NB)
    mkr() = Matrix{T}(undef, nep, NR)
    FR = typeof(sys.frame)
    return Trajectory{T,FR,_names(sys),typeof(epochs),Vector{T},Matrix{T},Vector{FR}}(
        epochs,
        Vector{FR}(undef, 1),
        vk(), vk(), vk(), vk(), vk(), vk(), vk(), vk(),
        vk(), vk(), vk(), vk(), vk(), vk(), vk(), vk(), vk(),
        mkb(),
        mkb(), mkb(), mkb(), mkb(), mkb(), mkb(),
        mkb(), mkb(), mkb(), mkb(), mkb(), mkb(),
        mkr(), mkr(), mkr(), mkr(), mkr(), mkr())
end

nepochs(traj::Trajectory) = length(traj.epochs)
_names(::Trajectory{<:Any,<:Any,Names}) where {Names} = Names

"""
    PlanetOrbits.frame(traj)
    PlanetOrbits.frame(sol)

The frame the trajectory was last solved against, written by `frame_pass!`.
Observables derive on-demand frame quantities from it; see `frame_ra`.
"""
@inline frame(traj::Trajectory) = @inbounds traj.frame[1]

# ---------------------------------------------------
# Per-epoch solution view
# ---------------------------------------------------

"""
    sol = traj[k]

Immutable zero-cost view of a `Trajectory` at epoch index `k`. Observables
(`raoff(sol, b, A)`, `radvel(sol, A, barycentre(sys))`, …) read from it.
"""
struct TrajectorySolution{TR<:Trajectory}
    traj::TR
    k::Int
end

Base.getindex(traj::Trajectory, k::Integer) = TrajectorySolution(traj, Int(k))
_names(sol::TrajectorySolution) = _names(sol.traj)
@inline frame(sol::TrajectorySolution) = frame(sol.traj)
Base.length(traj::Trajectory) = nepochs(traj)
Base.firstindex(traj::Trajectory) = 1
Base.lastindex(traj::Trajectory) = nepochs(traj)
Base.eachindex(traj::Trajectory) = Base.OneTo(nepochs(traj))
Base.iterate(traj::Trajectory, k=1) = k > nepochs(traj) ? nothing : (traj[k], k + 1)
Base.eltype(::Type{TR}) where {TR<:Trajectory} = TrajectorySolution{TR}

"""
    soltime(sol)

The *observation* epoch [MJD] this solution corresponds to — identically the
value passed into `orbitsolve`, preserving the epoch-indexing contract
(`soltime(traj[k]) === epochs[k]`).
"""
soltime(sol::TrajectorySolution) = sol.traj.epochs[sol.k]
export soltime
