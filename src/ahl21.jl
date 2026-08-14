# ---------------------------------------------------
# AHL21 propagator: symplectic N-body integration targeting the same
# Trajectory output as KeplerianApprox, so every observable works unchanged.
#
# Initial conditions: each hierarchy row's Keplerian solved at the osculating
# epoch t0 (elements are *osculating at t0* under this propagator), combined
# to absolute barycentric states via A⁻¹ — i.e. exactly the KeplerianApprox
# state at t0. The A-matrix guarantees zero total momentum, so the system
# barycentre stays at the origin and the frame layer composes identically.
#
# Output at the requested epochs: the map marches with the fixed step h;
# each epoch's state is produced by a *throwaway* partial step from the
# saved pre-step state (the NbodyGradient `findtransit!` pattern), never by
# odd fractional steps mid-sequence — the main march stays symplectic.
# Epochs before t0 are reached by marching the velocity-reversed state
# forward (the map is time-reversal equivariant).
# ---------------------------------------------------

# G in AU³ day⁻² M⊙⁻¹, derived from the same Kepler-year constant Row uses
# for its mean motion — the two propagators share one gravity by construction.
const G_au3_day2 = 4π^2 / kepler_year_to_julian_day_conversion_factor^2

"""
    AHL21(; h, t0=nothing)

N-body symplectic propagator of Agol, Hernandez & Langford (2021), MNRAS
507, 1582 (arXiv:2106.02188) — please cite this paper when publishing
results computed with it. Pass as `orbitsolve(sys, epochs; method=AHL21(h=…))`.

- `h`: fixed timestep [days]. Guidance: `h ≲ P_min/20` for the shortest
  period in the system; a warning is emitted (once) beyond that. Tight moons
  make likelihood evaluations expensive — the cost is one map evaluation per
  `h` of timespan covered, there is no adaptivity.
- `t0`: osculating epoch [MJD] — the orbital elements are interpreted as
  osculating elements at `t0`. Defaults to the frame `ref_epoch` for systems
  with an absolute frame; must be given explicitly otherwise.

Unlike `KeplerianApprox`, every body gravitates: hierarchical approximation
error disappears, and phenomena it cannot express (TTVs, resonant
interactions, moon-planet-star coupling) are captured. For a two-body
system the map is exact and the two propagators agree to roundoff.
"""
struct AHL21{TT<:Union{Nothing,Float64}} <: AbstractPropagator
    h::Float64
    t0::TT
end
function AHL21(; h::Real, t0::Union{Nothing,Real}=nothing)
    h > 0 || error("AHL21 timestep must be positive; got h=$h (epochs before t0 are handled automatically)")
    return AHL21{t0 === nothing ? Nothing : Float64}(Float64(h), t0 === nothing ? nothing : Float64(t0))
end

@inline _ahl21_t0(method::AHL21{Float64}, ::System) = method.t0
@inline _ahl21_t0(::AHL21{Nothing}, sys::System) = _default_t0(sys.frame)
@inline _default_t0(fr::AbsoluteFrame) = fr.ref_epoch
@noinline _default_t0(::AbstractFrame) = error(
    "AHL21 needs an osculating epoch for the orbital elements: pass AHL21(h=…, t0=…) " *
    "[MJD]. (Systems constructed with a full absolute frame default to its ref_epoch.)")

# h-guidance check: runs in the allocating `orbitsolve` entry points only —
# `orbitsolve!` is the allocation-free hot path and must not touch logging.
function _check_method(sys::System, method::AHL21)
    if method.h > _min_period_days(sys) / 20
        @warn "AHL21 timestep h exceeds P_min/20 for this system; integration error " *
              "grows as h² — reduce h or expect degraded accuracy." h = method.h maxlog = 1
    end
    return nothing
end

function propagate!(traj::Trajectory, sys::System{NB,NR}, method::AHL21) where {NB,NR}
    T = eltype(traj.x)
    t0 = T(_ahl21_t0(method, sys))
    h = method.h
    Gm = SVector{NB,T}(G_au3_day2 .* sys.masses)
    st0 = _initial_state(T, sys, t0)
    # The forward/backward split index is a property of the *values*: with a
    # Dual `t_em`, `searchsortedfirst`'s comparisons are ForwardDiff's
    # lexicographic ones, so an epoch sitting exactly on `t0` would be placed
    # in a different half depending on the sign of a partial — and the two
    # halves integrate in opposite time directions. `by=_primal` makes the
    # partition identical to the Float64 run's, which is the same guarantee
    # the step count below needs.
    kf = searchsortedfirst(traj.t_em, t0; by=_primal)
    _march!(traj, st0, Gm, h, t0, kf:nepochs(traj), 1)
    _march!(traj, _reverse_v(st0), Gm, h, t0, (kf-1):-1:1, -1)
    return traj
end

@inline _min_period_days(sys::System) =
    2π / maximum(r -> r.n, sys.rows) * year2day_julian

# The KeplerianApprox state at t0: per-row Kepler solve + A⁻¹ combine,
# with velocities converted from Trajectory units [AU/julian yr] to the
# integrator's [AU/day].
function _initial_state(::Type{T}, sys::System{NB,NR}, t0) where {T,NB,NR}
    rel = ntuple(Val(NR)) do j
        row = sys.rows[j]
        MA = row.n * day2year_julian * (t0 - row.tp)
        sE, cE = _anomaly_sincos(row, MA, Markley())
        NTuple{6,T}(_states_from_E(row, sE, cE))
    end
    xmat = SMatrix{3,NB,T}(_combine_ic(sys.Ainv, rel, 0, Val(3NB)))
    vmat = SMatrix{3,NB,T}(map(vc -> vc * T(day2year_julian),
                               _combine_ic(sys.Ainv, rel, 3, Val(3NB))))
    return NBodyState(xmat, vmat)
end

@inline function _combine_ic(Ainv::SMatrix{NB,NR}, rel, comp0::Int, ::Val{L}) where {NB,NR,L}
    ntuple(Val(L)) do lin
        c = (lin - 1) % 3 + 1 + comp0
        jb = (lin - 1) ÷ 3 + 1
        acc = zero(eltype(Ainv))
        for r in 1:NR
            acc = muladd(Ainv[jb, r], rel[r][c], acc)
        end
        acc
    end
end

# March over the epochs `ks` in pseudo-time τ = dir·(t_em − t0) ≥ 0,
# ascending along `ks`. The full-step chain is the symplectic backbone;
# outputs are throwaway partial steps.
function _march!(traj::Trajectory, st::NBodyState, Gm::SVector, h::Float64,
                 t0, ks, dir::Int)
    t_em = traj.t_em
    istep = 0
    @inbounds for k in ks
        τ = dir * (t_em[k] - t0)
        # How many full symplectic steps to take, and whether a partial step
        # is needed at all, are decided on the primal. Under Duals these are
        # comparisons ForwardDiff orders lexicographically, so an epoch landing
        # exactly on a step boundary — the common case, since `t0` is usually
        # the frame's `ref_epoch` and epochs are often regular — would take a
        # different *number of steps* depending on the sign of a partial. The
        # value carried through a gradient evaluation would then stop being
        # bit-identical to the plain Float64 one, and the derivative would be
        # of a different discretization than the value it accompanies.
        τp = _primal(τ)
        while τp - istep * h >= h
            st = ahl21_step(st, Gm, h)
            istep += 1
        end
        δ = τ - istep * h
        out = _primal(δ) > 0 ? ahl21_step(st, Gm, δ) : st
        _write_state!(traj, k, out, dir)
    end
    return nothing
end

@inline function _write_state!(traj::Trajectory, k::Int, st::NBodyState{N}, dir::Int) where {N}
    v2yr = dir * year2day_julian   # AU/day -> AU/julian yr, undoing reversal
    @inbounds for j in 1:N
        traj.x[k, j] = st.x[1, j] + st.xerr[1, j]
        traj.y[k, j] = st.x[2, j] + st.xerr[2, j]
        traj.z[k, j] = st.x[3, j] + st.xerr[3, j]
        traj.vx[k, j] = v2yr * (st.v[1, j] + st.verr[1, j])
        traj.vy[k, j] = v2yr * (st.v[2, j] + st.verr[2, j])
        traj.vz[k, j] = v2yr * (st.v[3, j] + st.verr[3, j])
    end
    return nothing
end

export AHL21
