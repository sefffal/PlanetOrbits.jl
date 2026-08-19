# ---------------------------------------------------
# Observing geometry — the pass that turns barycentric states into what an
# observer at the solar-system barycentre actually sees.
#
# `propagate!` leaves per-body states in the *reference* triad (the local
# East/North/line-of-sight frame in which the orbital elements are expressed)
# at the common, barycentre-retarded emission epoch. Four things separate that
# from an observation, and all four are applied here — once, to the stored
# states, so that every observable downstream is unchanged:
#
#   1. Viewing-direction change. The barycentre's own 3D motion swings our
#      line of sight; the tangent frame at the epoch is rotated relative to
#      the one at `ref_epoch`. The error from ignoring it is ρ·μT — 19 µas
#      over 40 yr for a companion at 0.5″ around a 200 mas/yr star, 1.1 mas
#      for a Barnard's-star-like host. Needs a full absolute frame.
#
#   2. Differential (per-body) light-travel time. Bodies at different depths
#      z along the line of sight are seen at different retarded times, t−z/c.
#      The positional consequence is ≈ 21–24 µas for a wide range of systems,
#      and the classic Rømer delay for timing. Needs no frame information at
#      all — z is in AU natively — so this one applies to every system.
#
#   3. Per-body angular scale. A body at depth z subtends its offset over
#      d+z, not d. Ignoring it costs ρ² in radians (≈ 4.85·ρ[″]² µas). Needs
#      a parallax.
#
#   4. Line-of-sight projection. Two bodies at different sky positions are
#      seen along *different* unit vectors, so a velocity common to both —
#      above all the barycentre's transverse space velocity — projects
#      differently onto each. This is the same geometry as the convergent
#      point of a moving cluster, and for a resolved pair it contributes
#      ≈ 0.023·μ[″/yr]·a[AU] m/s to the relative radial velocity — 24 cm/s
#      for a Barnard's-star-like host, and *independent of distance*. Unlike
#      1–3 it is not secular: it is there at every epoch, and because it
#      tracks the sky-plane separation it appears at the orbital period,
#      phase-locked to the astrometric orbit.
#
# Each correction is applied exactly when the frame carries the information to
# define it. A system whose bodies sit at inconsistent times or in
# inconsistent frames is not a physical model, so none of this is a runtime
# option.
#
# ---------------------------------------------------
# Loop structure
#
# The trajectory is epoch-fastest (§3.4), so every column is contiguous over
# epochs and per-epoch work is independent. The pass is therefore organised
# **epochs innermost**, which is what lets the divides and square roots that
# dominate it vectorize:
#
#   Pass A1 per epoch,  — frame geometry: barycentre distance, space velocity
#           vectorized     in the epoch triad, and the reference->epoch
#                          rotation, unpacked into scalar columns.
#   Pass A2 per body,   — apply that rotation to every body's state.
#           epochs
#           innermost
#   Pass B  per pair,   — barycentric acceleration and jerk of every body,
#           epochs        accumulated into scratch columns, and the potential
#           innermost     half of the Einstein term alongside them. This is the
#                         expensive one: it grows as NB² and was 51–64% of the
#                         pass when written epoch-outer.
#   Pass C  per body,   — retardation, line-of-sight projection and angular
#           epochs        scale, fused into a single traversal so each body's
#           innermost     state column is read and written exactly once.
#
# Pass B needs its output for *all* bodies before Pass C runs (retarding body
# i in place would perturb body j's acceleration), hence the scratch columns
# rather than a fused single pass.
# ---------------------------------------------------

"""
    observe_pass!(traj, sys)

Third and final pass of `orbitsolve!`. See the notes at the top of this file.
"""
# The four corrections above are deliberately gated by *one* flag, even
# though their costs differ by an order of magnitude (the NB² accel/jerk pair
# loop dominates; the depth scaling is a divide). Two reasons, so this is not
# re-litigated: they share the ρ-scaling that makes the opt-out decidable
# from a single property of the data, and splitting them would give a test
# matrix — and a documentation surface — nobody can keep honest. The Einstein
# term deliberately sits *outside* the flag: it is filled on both paths,
# because a flag about precision may not change what `radvel` means.
#
# Note the ρ-scaling argument is about the *astrometric* corrections. The
# line-of-sight projection (#4) does not vanish for an unresolved system in
# radial velocity — see the "Precision opt-outs" page.
# Reads the trajectory's *effective* frame (written by `frame_pass!`), not
# `sys.frame`: on the light-time path the two differ by the de-Doppler factor
# on the space velocity, and the observing geometry must see the same
# worldline the frame pass propagated.
observe_pass!(traj::Trajectory, sys::System) = observe_pass!(traj, sys, frame(traj))

# --- no frame: differential light-travel time only -------------------------

function observe_pass!(traj::Trajectory, sys::System{NB}, ::NoFrame) where {NB}
    T = eltype(traj.x)
    fill!(traj.d_au, T(NaN))
    fill!(traj.bvx, zero(T)); fill!(traj.bvy, zero(T)); fill!(traj.bvz, zero(T))
    fill!(traj.cart2angle, T(NaN))
    _accjerk_pass!(traj, sys, Val(NB))
    _retard_pass!(traj, Val(NB))
    return traj
end

# --- parallax: + angular scale and line-of-sight projection ----------------

function observe_pass!(traj::Trajectory, sys::System{NB}, fr::Parallax) where {NB}
    T = eltype(traj.x)
    fill!(traj.d_au, T(1000 / fr.plx * pc2au))
    # No space motion is known, so only the bodies' own velocities project.
    fill!(traj.bvx, zero(T)); fill!(traj.bvy, zero(T)); fill!(traj.bvz, zero(T))
    _accjerk_pass!(traj, sys, Val(NB))
    _observe_pass!(traj, Val(NB))
    return traj
end

# --- absolute frame: + viewing-direction change ----------------------------

function observe_pass!(traj::Trajectory, sys::System{NB}, fr::AbsoluteFrame) where {NB}
    T = eltype(traj.x)
    _frame_geometry_pass!(traj, fr)
    _rotate_pass!(traj, Val(NB))
    _accjerk_pass!(traj, sys, Val(NB))
    _observe_pass!(traj, Val(NB))
    return traj
end

# ---------------------------------------------------
# Pass A1: per-epoch observing geometry.
#
# Algebraically identical to `_geometry`, which stays as the readable
# single-epoch reference (and is what the test suite checks these columns
# against); written out in scalars here so the loop vectorizes. The local
# triad is built arithmetically from the Cartesian barycentre position — no
# trigonometry, and no round trip through (ra, dec).
# ---------------------------------------------------

function _frame_geometry_pass!(traj::Trajectory, fr::AbsoluteFrame)
    T = eltype(traj.x)
    M1 = fr.M1
    invyr = inv(T(year2day_julian))
    p2a = T(pc2au)
    @inbounds @simd for k in eachindex(traj.epochs)
        Δ = (traj.t_em[k] - fr.ref_epoch) * invyr
        x2 = fr.x1 + fr.dx * Δ
        y2 = fr.y1 + fr.dy * Δ
        z2 = fr.z1 + fr.dz * Δ
        d2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2)
        invd = inv(d2)
        erx = x2 * invd; ery = y2 * invd; erz = z2 * invd
        ρ = sqrt(x2 * x2 + y2 * y2)
        # On the pole East is degenerate; any finite choice keeps the triad
        # orthonormal. Selected branchlessly so the loop still vectorizes.
        invρ = iszero(_primal(ρ)) ? zero(ρ) : inv(ρ)
        eax = -y2 * invρ; eay = x2 * invρ          # eaz ≡ 0
        edx = -erz * eay; edy = erz * eax; edz = erx * eay - ery * eax
        # R = M2ᵀ M1, with M2's columns the (East, North, line-of-sight) triad
        traj.R11[k] = eax * M1[1, 1] + eay * M1[2, 1]
        traj.R12[k] = eax * M1[1, 2] + eay * M1[2, 2]
        traj.R13[k] = eax * M1[1, 3] + eay * M1[2, 3]
        traj.R21[k] = edx * M1[1, 1] + edy * M1[2, 1] + edz * M1[3, 1]
        traj.R22[k] = edx * M1[1, 2] + edy * M1[2, 2] + edz * M1[3, 2]
        traj.R23[k] = edx * M1[1, 3] + edy * M1[2, 3] + edz * M1[3, 3]
        traj.R31[k] = erx * M1[1, 1] + ery * M1[2, 1] + erz * M1[3, 1]
        traj.R32[k] = erx * M1[1, 2] + ery * M1[2, 2] + erz * M1[3, 2]
        traj.R33[k] = erx * M1[1, 3] + ery * M1[2, 3] + erz * M1[3, 3]
        traj.d_au[k] = d2 * p2a
        traj.bvx[k] = (eax * fr.dx + eay * fr.dy) * p2a
        traj.bvy[k] = (edx * fr.dx + edy * fr.dy + edz * fr.dz) * p2a
        traj.bvz[k] = (erx * fr.dx + ery * fr.dy + erz * fr.dz) * p2a
    end
    return traj
end

# ---------------------------------------------------
# Pass A2: rotate positions and velocities from the reference triad into the
# triad of the barycentre's apparent direction at each epoch.
# ---------------------------------------------------

function _rotate_pass!(traj::Trajectory, ::Val{NB}) where {NB}
    nk = length(traj.epochs)
    @inbounds for j in 1:NB
        @simd for k in 1:nk
            r11 = traj.R11[k]; r12 = traj.R12[k]; r13 = traj.R13[k]
            r21 = traj.R21[k]; r22 = traj.R22[k]; r23 = traj.R23[k]
            r31 = traj.R31[k]; r32 = traj.R32[k]; r33 = traj.R33[k]
            x = traj.x[k, j]; y = traj.y[k, j]; z = traj.z[k, j]
            vx = traj.vx[k, j]; vy = traj.vy[k, j]; vz = traj.vz[k, j]
            traj.x[k, j] = r11 * x + r12 * y + r13 * z
            traj.y[k, j] = r21 * x + r22 * y + r23 * z
            traj.z[k, j] = r31 * x + r32 * y + r33 * z
            traj.vx[k, j] = r11 * vx + r12 * vy + r13 * vz
            traj.vy[k, j] = r21 * vx + r22 * vy + r23 * vz
            traj.vz[k, j] = r31 * vx + r32 * vy + r33 * vz
        end
    end
    return traj
end

# ---------------------------------------------------
# Pass B: barycentric acceleration and jerk, by direct pair sum.
#
# The jerk carries the *velocity* retardation to second order. With velocity
# corrected only to first order the leading omission is ½ȧ(z/c)², which peaks
# at ~0.08 cm/s for a few-day-period orbit with a substantial reflex — inside
# one order of magnitude of the < 1 cm/s target this package is built for.
# With it, the same configurations agree with an exactly-retarded reference to
# ~1e-6 m/s.
#
# Both are accumulated in ONE traversal of each pair: the separation, r², its
# square root and the reciprocal are shared, so once the acceleration has been
# paid for the jerk costs only multiply-adds — no extra divide or sqrt, which
# is what actually sets the cost here.
#
# The accelerations are the direct N-body sum over the system's masses. That
# is exactly AHL21's own force; under `KeplerianApprox` it differs by the same
# approximation that propagator already makes, which is far below what a
# correction of relative size v/c ~ 1e-4 requires.
# ---------------------------------------------------

function _accjerk_pass!(traj::Trajectory, sys::System{NB}, ::Val{NB}) where {NB}
    T = eltype(traj.x)
    G = T(GM_sun_au3_julianyr2)
    nk = length(traj.epochs)
    fill!(traj.ax, zero(T)); fill!(traj.ay, zero(T)); fill!(traj.az, zero(T))
    fill!(traj.jx, zero(T)); fill!(traj.jy, zero(T)); fill!(traj.jz, zero(T))
    fill!(traj.ein, zero(T))
    @inbounds for i in 1:NB, j in 1:NB
        i == j && continue
        gm = G * sys.masses[j]
        # Test the *primal*: `iszero` on a Dual whose value is 0 but whose
        # partials are not returns false, so a differentiated zero mass would
        # otherwise reach the 0/0 below with a finite primal and NaN
        # gradients (§12).
        iszero(_primal(gm)) && continue
        @simd for k in 1:nk
            dx = traj.x[k, j] - traj.x[k, i]
            dy = traj.y[k, j] - traj.y[k, i]
            dz = traj.z[k, j] - traj.z[k, i]
            r2 = dx * dx + dy * dy + dz * dz
            # Branchless coincident-body guard: `inv(0)` is Inf, but the
            # select replaces it with 0 before it can reach a 0·Inf = NaN.
            # A `continue` here would defeat vectorization.
            invr2 = iszero(_primal(r2)) ? zero(r2) : inv(r2)
            invr = sqrt(invr2)
            f = gm * invr2 * invr
            # The Einstein term's potential half, free-riding on the `invr`
            # the force already needed: one multiply-add per pair-epoch.
            traj.ein[k, i] += gm * invr
            dvx = traj.vx[k, j] - traj.vx[k, i]
            dvy = traj.vy[k, j] - traj.vy[k, i]
            dvz = traj.vz[k, j] - traj.vz[k, i]
            rv3 = 3 * (dx * dvx + dy * dvy + dz * dvz) * invr2
            traj.ax[k, i] += f * dx
            traj.ay[k, i] += f * dy
            traj.az[k, i] += f * dz
            traj.jx[k, i] += f * (dvx - rv3 * dx)
            traj.jy[k, i] += f * (dvy - rv3 * dy)
            traj.jz[k, i] += f * (dvz - rv3 * dz)
        end
    end
    return traj
end

# ---------------------------------------------------
# Pass C: retardation, line-of-sight projection, angular scale.
#
# Retardation is by local Taylor extrapolation rather than re-propagation: the
# residual against an exactly-solved retarded time is ≤1e-3 µas and ≤0.1 mm/s
# across every regime the package targets. That matters structurally as well
# as for cost — re-propagating would force AHL21 into NB separate off-grid
# partial steps per epoch, and would make `t_em` a per-body rather than a
# per-epoch quantity.
#
# The position shift uses the body's *total* velocity: over the retardation
# interval the whole system translates by V·dt as well as the body orbiting by
# v·dt. Positions are measured from the barycentre at the common epoch t_em —
# the anchor the gnomonic projection uses — so the bulk term does not cancel,
# and it dominates whenever V_tan > v_orb, i.e. for most nearby fast-moving
# hosts. Velocities stay barycentre-relative, and V is constant, so only the
# orbital acceleration and jerk enter there.
#
# `velz` is then replaced by the velocity projected onto the body's *own*
# apparent direction, still measured relative to the barycentre so `radvel`
# stays a relative observable and `frame_rv` is not double-counted:
#
#     vz_obs = (V + v)·n̂  −  V·ê_r ,   n̂ = (x, y, d+z)/‖(x, y, d+z)‖
#
# Exact — no expansion in ρ. Reduces to vz identically for a body on the line
# of sight, and the V·ê_r terms cancel in any pairwise difference. Consequence
# for the physical-unit observables: `velz` is the *observed* line-of-sight
# velocity, so it is no longer the third Cartesian component of the same
# vector as `velx`/`vely`. That is what a spectrograph measures, and it is
# what keeps `radvel` honest.
# ---------------------------------------------------

function _observe_pass!(traj::Trajectory, ::Val{NB}) where {NB}
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
            # Retardation, with one refinement of dt: the depth itself changes
            # as the body is retarded, and for a fast-moving host that
            # second-order term is ~0.03 µas — small, but within 3× of the
            # 0.1 µas relative-astrometry target.
            dt = -z * invc
            dt = -(z + ((vz + Vz) * dt + az * dt * dt / 2)) * invc
            dt2 = dt * dt
            px = x + (vx + Vx) * dt + ax * dt2 / 2 + jx * dt2 * dt / 6
            py = y + (vy + Vy) * dt + ay * dt2 / 2 + jy * dt2 * dt / 6
            pz = z + (vz + Vz) * dt + az * dt2 / 2 + jz * dt2 * dt / 6
            nvx = vx + ax * dt + jx * dt2 / 2
            nvy = vy + ay * dt + jy * dt2 / 2
            nvz = vz + az * dt + jz * dt2 / 2
            # Line-of-sight projection and per-body angular scale.
            dzp = d_au + pz
            invr = inv(sqrt(px * px + py * py + dzp * dzp))
            traj.x[k, j] = px; traj.y[k, j] = py; traj.z[k, j] = pz
            # Einstein term: the potential half is already accumulated in
            # `ein` by `_accjerk_pass!`; add ½|v_tot|² and scale to a
            # velocity. `v_tot` is the body's *total* barycentric velocity,
            # orbital plus the frame's space velocity, which is what carries
            # the v_sys·v_orb/c cross term (~3 m/s for an SB system).
            tvx = nvx + Vx; tvy = nvy + Vy; tvz = nvz + Vz
            traj.ein[k, j] = (traj.ein[k, j] +
                              (tvx * tvx + tvy * tvy + tvz * tvz) / 2) * invc
            traj.vx[k, j] = nvx; traj.vy[k, j] = nvy
            traj.vz[k, j] = ((Vx + nvx) * px + (Vy + nvy) * py +
                             (Vz + nvz) * dzp) * invr - Vz
            traj.cart2angle[k, j] = r2m / dzp
        end
    end
    return traj
end

# Frameless variant: retardation only. There is no distance, so neither the
# angular scale nor the line-of-sight projection is defined, and `V` is zero.
function _retard_pass!(traj::Trajectory, ::Val{NB}) where {NB}
    T = eltype(traj.x)
    invc = inv(T(c_au_per_julianyr))
    nk = length(traj.epochs)
    @inbounds for j in 1:NB
        @simd for k in 1:nk
            x = traj.x[k, j]; y = traj.y[k, j]; z = traj.z[k, j]
            vx = traj.vx[k, j]; vy = traj.vy[k, j]; vz = traj.vz[k, j]
            ax = traj.ax[k, j]; ay = traj.ay[k, j]; az = traj.az[k, j]
            jx = traj.jx[k, j]; jy = traj.jy[k, j]; jz = traj.jz[k, j]
            dt = -z * invc
            dt = -(z + (vz * dt + az * dt * dt / 2)) * invc
            dt2 = dt * dt
            traj.x[k, j] = x + vx * dt + ax * dt2 / 2 + jx * dt2 * dt / 6
            traj.y[k, j] = y + vy * dt + ay * dt2 / 2 + jy * dt2 * dt / 6
            traj.z[k, j] = z + vz * dt + az * dt2 / 2 + jz * dt2 * dt / 6
            nvx = vx + ax * dt + jx * dt2 / 2
            nvy = vy + ay * dt + jy * dt2 / 2
            nvz = vz + az * dt + jz * dt2 / 2
            traj.vx[k, j] = nvx
            traj.vy[k, j] = nvy
            traj.vz[k, j] = nvz
            # No frame, so no space velocity: `v_tot` is the orbital velocity
            # alone. See `_observe_pass!`.
            traj.ein[k, j] = (traj.ein[k, j] +
                              (nvx * nvx + nvy * nvy + nvz * nvz) / 2) * invc
        end
    end
    return traj
end

# ---------------------------------------------------
# The opt-out: `orbitsolve(...; observing_geometry=false)`.
#
# Not a physics switch with a "no correction" branch — it selects the *cheap*
# geometry, which is exactly what v1 and PlanetOrbits v2 stage 1 computed: one
# shared AU→mas scale per epoch at the barycentre's distance, no viewing-
# direction rotation, no per-body retardation, no line-of-sight projection.
#
# It exists because all four corrections scale with the angular excursion
# actually observed, ρ:
#
#     rotation      ≈ 4.85e-3 · ρ[mas] · μ[″/yr] · T[yr]   µas
#     per-body LTT  ≈ 0.099   · ρ[mas] · √(M/a[AU])        µas
#     depth scale   ≈ 4.85e-6 · ρ[mas]²                    µas
#
# and for *absolute* astrometry ρ is the photocentre reflex, not the relative
# orbit. A Jupiter analogue at 10 pc gives ρ ≈ 0.475 mas, i.e. 0.005 / 0.021 /
# 1e-6 µas against a 30–100 µas per-epoch precision — four orders of magnitude
# of headroom, for ~2.5x on the gradient path. Fitting a few million such
# systems, that is worth declining; fitting one GRAVITY+/CRIRES+ target at
# 2 pc, where ρ is the full separation, it very much is not.
#
# **Defaults to on.** Skipping is an assertion about the precision of the data,
# which PlanetOrbits cannot see — only the caller can make it. It is a keyword
# on `orbitsolve`/`orbitsolve!` rather than a field of `System` deliberately:
# it is a statement about the *observations*, not about the physical system,
# and `System` is rebuilt from parameters every MCMC sample.
# ---------------------------------------------------

# Same effective-frame read as `observe_pass!`.
observe_skip!(traj::Trajectory, sys::System) = observe_skip!(traj, sys, frame(traj))

function observe_skip!(traj::Trajectory, sys::System{NB}, ::NoFrame) where {NB}
    T = eltype(traj.x)
    fill!(traj.d_au, T(NaN))
    fill!(traj.bvx, zero(T)); fill!(traj.bvy, zero(T)); fill!(traj.bvz, zero(T))
    fill!(traj.cart2angle, T(NaN))
    _ein_skip!(traj, sys, zero(T), zero(T), zero(T), Val(NB))
    return traj
end

function observe_skip!(traj::Trajectory, sys::System{NB}, fr::Parallax) where {NB}
    T = eltype(traj.x)
    d_au = T(1000 / fr.plx * pc2au)
    fill!(traj.d_au, d_au)
    fill!(traj.bvx, zero(T)); fill!(traj.bvy, zero(T)); fill!(traj.bvz, zero(T))
    fill!(traj.cart2angle, T(rad2mas) / d_au)
    _ein_skip!(traj, sys, zero(T), zero(T), zero(T), Val(NB))
    return traj
end

function observe_skip!(traj::Trajectory, sys::System{NB}, fr::AbsoluteFrame) where {NB}
    T = eltype(traj.x)
    r2m = T(rad2mas)
    p2a = T(pc2au)
    invyr = inv(T(year2day_julian))
    # No epoch triad is computed on this path, and the observer-aware
    # observables need one. Poison the first rotation column rather than
    # leave the caller's (bump-allocated, uninitialized) storage to be read
    # as though it were a rotation; `_observer_offset` tests it and throws.
    fill!(traj.R11, T(NaN))
    fill!(traj.bvx, zero(T)); fill!(traj.bvy, zero(T)); fill!(traj.bvz, zero(T))
    @inbounds @simd for k in eachindex(traj.epochs)
        Δ = (traj.t_em[k] - fr.ref_epoch) * invyr
        x2 = fr.x1 + fr.dx * Δ
        y2 = fr.y1 + fr.dy * Δ
        z2 = fr.z1 + fr.dz * Δ
        traj.d_au[k] = sqrt(x2 * x2 + y2 * y2 + z2 * z2) * p2a
    end
    @inbounds for j in 1:NB
        @simd for k in eachindex(traj.epochs)
            traj.cart2angle[k, j] = r2m / traj.d_au[k]
        end
    end
    # Skip mode leaves the states in the *reference* triad and `bv*` zeroed,
    # so the space velocity is resolved once from the frame directly rather
    # than read per epoch — it is constant in that triad.
    V = (fr.M1' * SVector(fr.dx, fr.dy, fr.dz)) * p2a
    _ein_skip!(traj, sys, T(V[1]), T(V[2]), T(V[3]), Val(NB))
    return traj
end

# ---------------------------------------------------
# The Einstein term on the skipped path.
#
# `radvel` is the spectroscopic velocity at *every* setting of
# `observing_geometry`: the flag chooses the precision of the geometry, and
# is not allowed to change what an observable means. So the skipped path pays
# its own pair loop here — the potential is not otherwise computed, because
# `_accjerk_pass!` never ran.
#
# Cost is NB² multiply-adds plus one sqrt per pair-epoch, taken
# unconditionally. There is deliberately no third opt-out flag for it: the
# combinatorics of the existing two are already the limit of what can be
# tested and documented honestly.
# ---------------------------------------------------

function _ein_skip!(traj::Trajectory, sys::System{NB}, Vx, Vy, Vz,
                    ::Val{NB}) where {NB}
    T = eltype(traj.x)
    G = T(GM_sun_au3_julianyr2)
    invc = inv(T(c_au_per_julianyr))
    nk = length(traj.epochs)
    fill!(traj.ein, zero(T))
    @inbounds for i in 1:NB, j in 1:NB
        i == j && continue
        gm = G * sys.masses[j]
        # Same zero-mass-dual and coincident-body guards as `_accjerk_pass!`.
        iszero(_primal(gm)) && continue
        @simd for k in 1:nk
            dx = traj.x[k, j] - traj.x[k, i]
            dy = traj.y[k, j] - traj.y[k, i]
            dz = traj.z[k, j] - traj.z[k, i]
            r2 = dx * dx + dy * dy + dz * dz
            invr2 = iszero(_primal(r2)) ? zero(r2) : inv(r2)
            traj.ein[k, i] += gm * sqrt(invr2)
        end
    end
    @inbounds for j in 1:NB
        @simd for k in 1:nk
            tvx = traj.vx[k, j] + Vx
            tvy = traj.vy[k, j] + Vy
            tvz = traj.vz[k, j] + Vz
            traj.ein[k, j] = (traj.ein[k, j] +
                              (tvx * tvx + tvy * tvy + tvz * tvz) / 2) * invc
        end
    end
    return traj
end
