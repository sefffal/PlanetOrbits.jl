# ---------------------------------------------------
# Trajectory: struct-of-arrays solution storage
#
# Epoch-fastest layout: every column is contiguous over epochs, so the
# per-row Kepler solve batches (and later SIMD-vectorizes) across epochs,
# and per-chunk epoch ranges are contiguous for threading.
#
# The frame *type* FM is carried in the trajectory's type so observables can
# dispatch on what frame information is available; the per-epoch frame data
# lives in the columns.
# ---------------------------------------------------

struct Trajectory{T<:Number,FM<:FrameMode,E<:AbstractVector{<:Real},
                  V<:AbstractVector{T},M<:AbstractMatrix{T}}
    epochs::E    # observation epochs [MJD], sorted
    # per-epoch frame columns (shared by every body and observable)
    t_em::V          # light-travel-time–corrected emission epoch [MJD]
    cart2angle::V    # AU -> mas at the epoch's distance
    ra2::V; dec2::V  # propagated frame position [deg]
    pmra2::V; pmdec2::V  # propagated frame proper motion [mas/yr]
    rv2::V           # propagated frame radial velocity [m/s]
    # per-body absolute barycentric states [AU, AU/julian year]
    x::M; y::M; z::M; vx::M; vy::M; vz::M
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
    return Trajectory{T,_frame_mode(sys.frame),typeof(epochs),Vector{T},Matrix{T}}(
        epochs,
        vk(), vk(), vk(), vk(), vk(), vk(), vk(),
        mkb(), mkb(), mkb(), mkb(), mkb(), mkb(),
        mkr(), mkr(), mkr(), mkr(), mkr(), mkr())
end

_frame_mode(::NoFrame) = ModeNone
_frame_mode(::Parallax) = ModeParallax
_frame_mode(::AbsoluteFrame) = ModeAbsolute

nepochs(traj::Trajectory) = length(traj.epochs)

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
