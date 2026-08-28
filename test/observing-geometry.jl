# Observing geometry (src/observe.jl): the four corrections that separate a
# barycentric state from what an observer at the solar-system barycentre
# actually sees.
#
# The load-bearing test here is the brute-force reference: it rebuilds the
# entire observing chain from full 3D ICRS vectors, with an *exactly* solved
# per-body retarded time rather than the production Taylor step, and with the
# local triads built from (ra, dec) by trigonometry rather than from Cartesian
# components. It shares only the Kepler solve with production.

using LinearAlgebra
using PlanetOrbits: Body, Orbit, System, Trajectory, frame_pass!, propagate!,
                    KeplerianApprox, _geometry, _triad, observe_pass!,
                    pc2sec_light, au2pc, pc2au, pc2m, rad2mas, au2m,
                    year2day_julian, year2sec_julian, c_au_per_julianyr,
                    rad2as_206265, pc2km, sec2day, c_light_ms,
                    GM_sun_au3_julianyr2

const AS2RAD = π / 180 / 3600
const MAS2RAD = AS2RAD / 1000

# Raw reference-triad states of the *unretarded, unrotated* system at time t.
function _rawstate(sysnf, t)
    traj = Trajectory(sysnf, [t])
    frame_pass!(traj, sysnf.frame)
    propagate!(traj, sysnf, KeplerianApprox())
    n = size(traj.x, 2)
    r = [SVector(traj.x[1, j], traj.y[1, j], traj.z[1, j]) for j in 1:n]
    v = [SVector(traj.vx[1, j], traj.vy[1, j], traj.vz[1, j]) for j in 1:n]
    return r, v
end

# Independent reimplementation of the whole observing chain.
function brute_force(sys, sysnf, t_obs)
    fr = sys.frame
    sa, ca = sincosd(fr.ra)
    sd, cd = sincosd(fr.dec)
    ea = SVector(-sa, ca, 0.0)
    ed = SVector(-sd * ca, -sd * sa, cd)
    er = SVector(cd * ca, cd * sa, sd)
    M1 = hcat(ea, ed, er)
    d1 = 1000 / fr.plx                                   # pc
    B1 = er * d1
    # Catalog proper motions are apparent rates; the true worldline the
    # light-time solve runs on carries the de-Dopplered tangential velocity
    # μ·(1 + v_r/c) — B&L 2014 Eqs. 11–13, same as production's `_dedoppler`
    # (gated independently against the rigorous oracle in lighttime-bl2014.jl).
    dop = 1 + fr.rv / c_light_ms
    V = (fr.pmra * dop / 1000 * AS2RAD * d1) * ea +      # pc / julian yr
        (fr.pmdec * dop / 1000 * AS2RAD * d1) * ed +
        (fr.rv * year2sec_julian / pc2m) * er
    Bt(t) = B1 + V * ((t - fr.ref_epoch) / year2day_julian)

    # Barycentre emission epoch. Constant d/c dropped (degenerate with tp),
    # matching production's convention.
    t_em = t_obs
    for _ in 1:80
        t_em = t_obs - (norm(Bt(t_em)) - d1) * pc2sec_light / 86400
    end
    B2 = Bt(t_em)
    d2 = norm(B2)
    er2 = B2 / d2
    ρ = hypot(B2[1], B2[2])
    ea2 = SVector(-B2[2] / ρ, B2[1] / ρ, 0.0)
    ed2 = er2 × ea2

    nb = length(_rawstate(sysnf, t_obs)[1])
    ξ = zeros(nb); η = zeros(nb); vr = zeros(nb)
    speeds = zeros(nb)
    for j in 1:nb
        # Per-body emission epoch, solved exactly — no Taylor expansion.
        t_j = t_em
        local P, Vel
        for _ in 1:80
            r, v = _rawstate(sysnf, t_j)
            P = Bt(t_j) + M1 * (r[j] * au2pc)
            Vel = V + M1 * (v[j] * au2pc)
            t_j = t_obs - (norm(P) - d1) * pc2sec_light / 86400
        end
        ξ[j] = (P ⋅ ea2) / (P ⋅ er2) / MAS2RAD
        η[j] = (P ⋅ ed2) / (P ⋅ er2) / MAS2RAD
        n̂ = P / norm(P)
        vr[j] = (Vel ⋅ n̂) * pc2m / year2sec_julian       # m/s
        speeds[j] = norm(Vel) * pc2m / year2sec_julian   # m/s, total
    end
    # Einstein term, added so `radvel` is checked whole rather than only its
    # kinematic half. The potential is evaluated at the *common* emission
    # epoch and from separations, matching production's Pass B (which runs
    # before retardation); the speed is the per-body retarded one. Both
    # choices are second-order differences, inside the tolerances below.
    r_em, _ = _rawstate(sysnf, t_em)
    gfac = (au2m / year2sec_julian)^2                    # (AU/yr)² -> (m/s)²
    for j in 1:nb
        Φ = 0.0
        for i in 1:nb
            i == j && continue
            Φ += GM_sun_au3_julianyr2 * sys.masses[i] / norm(r_em[j] - r_em[i])
        end
        vr[j] += (speeds[j]^2 / 2 + Φ * gfac) / c_light_ms
    end
    return (; ξ, η, vr)
end

function geometry_case(; d, mu_as, a, M, m2, ra=45.0, dec=20.0, rv=-30e3, T=40.0)
    plx = 1000 / d
    ref_epoch = 57388.0
    A = Body(mass=M - m2, name=:A)
    b = Body(mass=m2, name=:b)
    orb = Orbit(b, about=A; a=a, e=0.2, i=0.9, ω=1.1, Ω=2.2, tp=57000.0)
    kw = (; plx, ra, dec, pmra=mu_as * 1000 / sqrt(2), pmdec=mu_as * 1000 / sqrt(2),
          rv, ref_epoch)
    return (System((A, b), (orb,); kw...), System((A, b), (orb,)),
            ref_epoch + T * 365.25)
end

@testset "observing geometry" begin

    @testset "vs brute-force 3D reference" begin
        # (name, d[pc], μ["/yr], a[AU], M, m2, rv[m/s])
        cases = [
            ("Barnard-like",      1.83, 10.357, 1.0,  0.16, 1.5e-6,  -110e3),
            ("Solar analog 10pc", 10.0,  0.200, 5.0,  1.00, 9.5e-4,    20e3),
            ("beta Pic-like",     19.6,  0.050, 10.0, 1.75, 1.05e-2,   20e3),
            ("tau Ceti-like",     3.65,  1.922, 1.0,  0.78, 9.1e-6,   -17e3),
            ("GJ 876-like",       4.68,  1.000, 0.21, 0.37, 2.2e-3,   -1.6e3),
        ]
        for (nm, d, mu, a, M, m2, rv) in cases
            @testset "$nm" begin
                sys, sysnf, t_obs = geometry_case(; d, mu_as=mu, a, M, m2, rv)
                sol = orbitsolve(sys, [t_obs])[1]
                R = brute_force(sys, sysnf, t_obs)
                rA, rb = bodies(sys)
                # 1e-3 µas — a thousand times below the 0.1 µas relative
                # astrometry design target. What remains is the third-order
                # term of the (twice-refined) retardation.
                @test abs(raoff(sol, rb, rA) - (R.ξ[2] - R.ξ[1])) * 1e3 < 1e-3
                @test abs(decoff(sol, rb, rA) - (R.η[2] - R.η[1])) * 1e3 < 1e-3
                # 0.1 mm/s, i.e. ≥100x below the < 1 cm/s design target. The
                # residual is the second-order term of the *velocity*
                # retardation (v is corrected to first order in the
                # acceleration, so the jerk term survives); it is largest for
                # tight orbits with a big reflex, hence GJ 876 dominating.
                @test abs(radvel(sol, rb, rA) - (R.vr[2] - R.vr[1])) * 1e3 < 1e-1
            end
        end
    end

    @testset "barycentric light-travel sign" begin
        # A receding system's light path lengthens with time, so its apparent
        # period must be LONGER than its true period. v0.11 had this inverted.
        A = Body(mass=1.0, name=:A)
        b = Body(mass=1e-6, name=:b)
        orb = Orbit(b, about=A; a=0.1, e=0.0, i=π / 2, ω=0.0, Ω=0.0, tp=57388.0)
        base = (; plx=100.0, ra=45.0, dec=20.0, pmra=0.0, pmdec=0.0, ref_epoch=57388.0)
        ts = [57388.0, 67388.0]
        for (rv, longer) in ((+300e3, true), (-300e3, false))
            sys = System((A, b), (orb,); base..., rv)
            traj = Trajectory(sys, ts)
            frame_pass!(traj, sys.frame)
            ratio = (ts[2] - ts[1]) / (traj.t_em[2] - traj.t_em[1])  # P_obs / P_true
            @test (ratio > 1) == longer
            # magnitude: dt_em/dt_obs = 1/(1 + v_r/c)
            @test ratio ≈ 1 + rv / 2.99792458e8 rtol = 1e-4
        end
    end

    @testset "per-body angular scale is exact" begin
        # With a massless companion the host sits at the barycentre, so
        # raoff·(d + z) == posx·rad2mas identically: the depth factor is the
        # body's own, not a shared per-epoch scale.
        A = Body(mass=1.0, name=:A)
        b = Body(mass=0.0, name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=5.0, e=0.3, i=1.0, ω=1.1, Ω=2.2,
                                    tp=57000.0),); plx=100.0)
        d_au = 1000 / 100.0 * pc2au
        rA, rb = bodies(sys)
        for sol in orbitsolve(sys, [57500.0, 58000.0, 59000.0])
            @test raoff(sol, rb, rA) * (d_au + posz(sol, rb, rA)) ≈ posx(sol, rb, rA) * rad2mas rtol = 1e-13
            @test decoff(sol, rb, rA) * (d_au + posz(sol, rb, rA)) ≈ posy(sol, rb, rA) * rad2mas rtol = 1e-13
            # ...and the size of the departure from a shared scale is ρ² in
            # radians, i.e. ≈ 4.85·ρ[″]² µas.
            shared = posx(sol, rb, rA) * rad2mas / d_au
            ρ_as = projectedseparation(sol, rb, rA) / 1000
            @test abs(raoff(sol, rb, rA) - shared) * 1e3 < 6 * ρ_as^2
        end
    end

    @testset "radial velocity: differing lines of sight" begin
        # Two bodies at different sky positions are seen along different unit
        # vectors, so the barycentre's transverse velocity projects
        # differently onto each. Isolated by differencing against the same
        # system with zero proper motion, evaluated at ref_epoch so that no
        # time has elapsed (no rotation, no secular light-travel term).
        ref_epoch = 57388.0
        A = Body(mass=0.16, name=:A)
        b = Body(mass=1.5e-6, name=:b)
        orb = Orbit(b, about=A; a=1.0, e=0.2, i=0.9, ω=1.1, Ω=2.2, tp=57000.0)
        d = 1.83
        base = (; plx=1000 / d, ra=45.0, dec=20.0, rv=-110e3, ref_epoch)
        pm = 10.357 * 1000 / sqrt(2)
        sys_pm = System((A, b), (orb,); base..., pmra=pm, pmdec=pm)
        sys_0 = System((A, b), (orb,); base..., pmra=0.0, pmdec=0.0)
        rA, rb = bodies(sys_pm)
        sol_pm = orbitsolve(sys_pm, [ref_epoch])[1]
        sol_0 = orbitsolve(sys_0, [ref_epoch])[1]
        # The kinematic quantity: the Einstein term differs between the two
        # systems too (|v_tot|² contains the transverse space velocity, worth
        # a constant ~13 m/s here), and that is not the effect under test.
        Δ = kinrv(sol_pm, rb, rA) - kinrv(sol_0, rb, rA)
        # Closed form: V_tan · ρ⃗, with V_tan[m/s] = 4.74047·μ["/yr]·d[pc]·1e3
        # and ρ⃗ the sky-plane separation in radians.
        Vα = 4.74047 * (pm / 1000) * d * 1e3
        Vδ = Vα
        pred = Vα * raoff(sol_pm, rb, rA) * MAS2RAD + Vδ * decoff(sol_pm, rb, rA) * MAS2RAD
        @test Δ ≈ pred rtol = 0.02
        # ...and it is a real effect at the package's stated RV precision.
        @test abs(Δ) > 0.05   # m/s
    end

    @testset "viewing-direction rotation" begin
        sys, _, t_obs = geometry_case(; d=1.83, mu_as=10.357, a=1.0, M=0.16,
                                      m2=1.5e-6, rv=-110e3)
        fr = sys.frame
        R, d_au, V = _geometry(fr, t_obs, Float64)
        @test R' * R ≈ I atol = 1e-14
        @test det(R) ≈ 1 atol = 1e-14
        # The rotation angle is the barycentre's own angular displacement.
        θ = acos(clamp(R[3, 3], -1, 1))
        μ_total = hypot(fr.pmra, fr.pmdec) / 1000 * AS2RAD          # rad/yr
        # Only approximately μT: the proper motion is quoted at ref_epoch and
        # the distance changes by ~0.5% over the baseline, so the accumulated
        # angular displacement departs from the linear estimate at that level.
        @test θ ≈ μ_total * (t_obs - fr.ref_epoch) / year2day_julian rtol = 1e-2
        # At the reference epoch it is the identity.
        R0, _, _ = _geometry(fr, fr.ref_epoch, Float64)
        @test R0 ≈ I atol = 1e-14
    end

    @testset "observing_geometry=false is exactly the v0.11 / stage-1 geometry" begin
        # The opt-out selects the *cheap* geometry, not a different physics:
        # one shared AU->mas scale per epoch, no rotation, no retardation, no
        # line-of-sight projection. For the frameless and parallax-only
        # fixtures that is bit-for-bit what v0.11 computed, which is what proves
        # the observing pass is the only thing that changed for them.
        for c in V0_REFERENCE
            c.kind === :absvis && continue
            sys, _ = fixture_system(c)
            traj = orbitsolve(sys, c.epochs; observing_geometry=false)
            d = c.data
            for (k, sol) in enumerate(traj)
                @test posx(sol) ≈ d.posx[k] rtol = 1e-13 atol = 1e-13
                @test posz(sol) ≈ d.posz[k] rtol = 1e-13 atol = 1e-13
                @test velz(sol) ≈ d.velz[k] rtol = 1e-13 atol = 1e-13
                if c.kind === :visual
                    @test raoff(sol) ≈ d.raoff[k] rtol = 1e-13 atol = 1e-13
                    @test decoff(sol) ≈ d.decoff[k] rtol = 1e-13 atol = 1e-13
                end
            end
        end
        # For the absolute-frame fixtures a residual remains, and it is
        # entirely the barycentric light-travel *sign* fix — the one v0.11 bug
        # that is not part of the observing pass and so is not opted out of.
        for c in V0_REFERENCE
            c.kind === :absvis || continue
            sys, _ = fixture_system(c)
            traj = orbitsolve(sys, c.epochs; observing_geometry=false)
            dev = maximum(abs(posx(traj[k]) - c.data.posx[k]) /
                          max(abs(c.data.posx[k]), 1e-30) for k in eachindex(c.epochs))
            @test 1e-4 < dev < 1e-1
        end
    end

    @testset "observing_geometry=false is cheaper, allocation-free, and differs" begin
        sys, _, t_obs = geometry_case(; d=1.83, mu_as=10.357, a=1.0, M=0.16,
                                      m2=1.5e-6, rv=-110e3)
        eps = collect(range(57388.0, t_obs, length=64))
        rA, rb = bodies(sys)
        full = orbitsolve(sys, eps)
        cheap = orbitsolve(sys, eps; observing_geometry=false)
        # The corrections are worth ~1 mas here, so the two must not agree.
        Δ = maximum(abs(raoff(full[k], rb, rA) - raoff(cheap[k], rb, rA))
                    for k in eachindex(eps))
        @test Δ > 0.1                                  # mas
        traj = Trajectory(sys, eps)
        orbitsolve!(traj, sys; observing_geometry=false)
        if DYNAMIC_ALLOC_GATE
            @test (@allocated orbitsolve!(traj, sys; observing_geometry=false)) == 0
        else
            @test_skip (@allocated orbitsolve!(traj, sys; observing_geometry=false)) == 0
        end
    end

    @testset "vectorized geometry pass matches the _geometry reference" begin
        # `_frame_geometry_pass!` is `_geometry` written out in scalars so the
        # loop vectorizes. Nothing else keeps the two in step, so assert it.
        sys, _, _ = geometry_case(; d=1.83, mu_as=10.357, a=1.0, M=0.16,
                                  m2=1.5e-6, rv=-110e3)
        fr = sys.frame
        eps = collect(range(fr.ref_epoch - 8000, fr.ref_epoch + 8000, length=23))
        traj = Trajectory(sys, eps)
        frame_pass!(traj, fr)
        PlanetOrbits._frame_geometry_pass!(traj, fr)
        for k in eachindex(eps)
            R, d_au, V = _geometry(fr, traj.t_em[k], Float64)
            @test traj.d_au[k] ≈ d_au rtol = 1e-15
            @test traj.bvx[k] ≈ V[1] rtol = 1e-12
            @test traj.bvy[k] ≈ V[2] rtol = 1e-12
            @test traj.bvz[k] ≈ V[3] rtol = 1e-12
            Rcols = (traj.R11[k], traj.R12[k], traj.R13[k],
                     traj.R21[k], traj.R22[k], traj.R23[k],
                     traj.R31[k], traj.R32[k], traj.R33[k])
            @test all(isapprox(Rcols[3(r-1)+c], R[r, c]; atol=1e-14)
                      for r in 1:3, c in 1:3)
        end
    end

    @testset "barycentric_lighttime=false skips only the light-travel solve" begin
        # The two opt-outs gate *different* corrections. `observing_geometry`
        # gates terms scaling with the system's angular extent ρ;
        # `barycentric_lighttime` gates a whole-system timing correction
        # scaling with proximity and proper motion. This asserts the second one
        # does exactly that and nothing more — in particular that it does NOT
        # degrade into "stop propagating the frame", which would silently
        # delete the perspective acceleration.
        sys, _, _ = geometry_case(; d=1.83, mu_as=10.357, a=1.0, M=0.16,
                                  m2=1.5e-6, rv=-110e3)
        fr = sys.frame
        eps = collect(range(fr.ref_epoch - 3652.5, fr.ref_epoch + 3652.5, length=25))
        on = orbitsolve(sys, eps)
        off = orbitsolve(sys, eps; barycentric_lighttime=false)

        # t_em is exactly the observation epoch when the solve is skipped...
        @test all(off.t_em .== eps)
        # ...and is not, when it is not.
        @test maximum(abs, on.t_em .- eps) > 0.1        # days, this is a fast nearby case

        # The frame is still propagated: pm/rv vary across epochs and match a
        # direct evaluation at t_obs — they are NOT frozen at catalog values.
        for k in eachindex(eps)
            kin = PlanetOrbits._compensate_kinematics(fr, eps[k])
            @test frame_pmra(off[k]) == kin.pmra2
            @test frame_pmdec(off[k]) == kin.pmdec2
            @test frame_rv(off[k]) == kin.rv2
        end
        @test frame_pmra(off[1]) != frame_pmra(off[end])
        @test frame_rv(off[1]) != frame_rv(off[end])
        # Perspective acceleration survives: μ is not linear in t, so the
        # second difference is nonzero.
        d2 = frame_pmdec(off[1]) - 2frame_pmdec(off[13]) + frame_pmdec(off[25])
        @test abs(d2) > 0

        # It is a real correction, not a no-op: turning it off moves observables.
        rA, rb = bodies(sys)
        @test maximum(k -> abs(raoff(on[k], rb, rA) - raoff(off[k], rb, rA)),
                      eachindex(eps)) > 0

        # Orthogonal to observing_geometry: all four combinations are distinct
        # and each flag changes the same thing regardless of the other.
        combos = [orbitsolve(sys, eps; observing_geometry=g, barycentric_lighttime=l)
                  for g in (true, false), l in (true, false)]
        vals = [[raoff(t[k], rb, rA) for k in eachindex(eps)] for t in combos]
        for i in 1:4, j in (i+1):4
            @test vals[i] != vals[j]
        end

        # No-op on frames that have no barycentric light-travel time at all.
        mkcheap(; kw...) = let A = Body(mass=1.0, name=:A), b = Body(mass=0.0, name=:b)
            System((A, b), (Orbit(b, about=A; a=1.0, e=0.0, i=0.0, ω=0.0,
                                  Ω=0.0, tp=57000.0),); kw...)
        end
        for cheap in (mkcheap(), mkcheap(plx=25.0))
            a = orbitsolve(cheap, eps)
            b = orbitsolve(cheap, eps; barycentric_lighttime=false)
            @test a.t_em == b.t_em
            @test all(posx(a[k]) == posx(b[k]) for k in eachindex(eps))
        end

        # Allocation-free, and cheaper.
        traj = Trajectory(sys, eps)
        orbitsolve!(traj, sys; barycentric_lighttime=false)
        if DYNAMIC_ALLOC_GATE
            @test (@allocated orbitsolve!(traj, sys; barycentric_lighttime=false)) == 0
        else
            @test_skip (@allocated orbitsolve!(traj, sys; barycentric_lighttime=false)) == 0
        end
    end

    @testset "frame_ra/frame_dec are on demand and cannot go stale" begin
        # `ra2`/`dec2` are no longer stored; `frame_ra`/`frame_dec` recompute
        # from `t_em` and the trajectory's frame. That frame is written by
        # `frame_pass!`, not captured at construction — precisely so that
        # reusing trajectory buffers across a `sys` rebuild (the sampling hot
        # loop, and `perf/gaia-workload.jl`) cannot silently serve the previous
        # sample's frame. Assert both halves.
        mk(ra, dec, pmra) = let A = Body(mass=1.0, name=:A), b = Body(mass=0.0, name=:b)
            System((A, b), (Orbit(b, about=A; a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0,
                                  tp=57000.0),);
                   plx=50.0, ra=ra, dec=dec, pmra=pmra, pmdec=20.0, rv=1e4,
                   ref_epoch=57388.5)
        end
        eps = collect(range(57000.0, 61000.0, length=7))
        sys1 = mk(45.0, -30.0, 100.0)
        sys2 = mk(200.0, 55.0, -900.0)

        # Built against sys1, then solved against sys2 with the SAME buffers.
        traj = Trajectory(sys1, eps)
        orbitsolve!(traj, sys2)
        ref = orbitsolve(sys2, eps)
        for k in eachindex(eps)
            @test frame_ra(traj[k]) ≈ frame_ra(ref[k]) rtol = 1e-14
            @test frame_dec(traj[k]) ≈ frame_dec(ref[k]) rtol = 1e-14
        end
        # ...and it really is sys2's frame, not sys1's — the *effective*
        # (de-Dopplered) form of it, since the default solve has the
        # light-travel path on.
        @test !isapprox(frame_ra(traj[1]), 45.0; atol=1.0)
        @test PlanetOrbits.frame(traj) === PlanetOrbits._effective_frame(sys2.frame, true)

        # Re-solving against sys1 through the same buffers switches back.
        orbitsolve!(traj, sys1)
        ref1 = orbitsolve(sys1, eps)
        for k in eachindex(eps)
            @test frame_ra(traj[k]) ≈ frame_ra(ref1[k]) rtol = 1e-14
            @test frame_dec(traj[k]) ≈ frame_dec(ref1[k]) rtol = 1e-14
        end

        # Reading them is allocation-free (they are computed, not stored).
        sol = traj[3]
        frame_ra(sol); frame_dec(sol)
        if DYNAMIC_ALLOC_GATE
            @test (@allocated frame_ra(sol)) == 0
            @test (@allocated frame_dec(sol)) == 0
        else
            @test_skip (@allocated frame_ra(sol)) == 0
            @test_skip (@allocated frame_dec(sol)) == 0
        end
    end

    @testset "compensate: algebraic form matches the literal transcription" begin
        # `compensate` is the single largest term in an absolute-frame
        # `orbitsolve!`, so it uses three identities to drop a `hypot`, a `cos`
        # and a duplicated `sqrt`:
        #     hypot(x,y,z)  ==  √(x²+y²+z²)          (components are parsecs)
        #     cosd(dec2)    ==  √(1 − sin²δ)          (δ ∈ [−90°, 90°])
        #     √(1 − sin²δ)  ==  the ∂δ/∂t denominator factor √(1 − z²/d²)
        # Nothing else keeps the optimized form honest, so assert it against a
        # literal transcription of v0.11's `compensate_star_3d_motion`.
        function compensate_literal(fr, t_em_days)
            if fr.ref_epoch == t_em_days
                t_em_days += eps(float(t_em_days))
            end
            Δt_jyear = (t_em_days - fr.ref_epoch) / year2day_julian
            x2 = fr.x1 + fr.dx * Δt_jyear
            y2 = fr.y1 + fr.dy * Δt_jyear
            z2 = fr.z1 + fr.dz * Δt_jyear
            distance2 = hypot(x2, y2, z2)
            ra2 = (atand(y2, x2) + 360) % 360
            arg = clamp(z2 / distance2, -1.0, 1.0)
            dec2 = asind(arg)
            ddist2 = (x2 * fr.dx + y2 * fr.dy + z2 * fr.dz) / distance2
            dra2 = (-y2 * fr.dx + x2 * fr.dy) / (x2^2 + y2^2)
            ddec2 = (-z2 * ddist2 / distance2 + fr.dz) /
                    (distance2 * sqrt(1 - z2^2 / distance2^2))
            pmra2 = dra2 * rad2as_206265 * 1000 * cosd(dec2)
            pmdec2 = ddec2 * 1000 * rad2as_206265
            rv2 = ddist2 * pc2km / year2sec_julian * 1e3
            epoch2a_days = t_em_days +
                           (distance2 - fr.distance1) * pc2sec_light * sec2day
            return (; ra2, dec2, pmra2, pmdec2, rv2, epoch2a_days, distance2)
        end

        # Spread over the sky, including a fast, nearby, high-declination case
        # where cos δ is small and the ∂δ/∂t denominator is worst-conditioned,
        # and a southern retreating one.
        cases = (
            (ra=45.0, dec=-30.0, plx=24.5, pmra=100.0, pmdec=-50.0, rv=25e3),
            (ra=269.4, dec=4.7, plx=546.98, pmra=-802.8, pmdec=10362.5, rv=-110e3),
            (ra=10.0, dec=87.5, plx=120.0, pmra=-1500.0, pmdec=900.0, rv=60e3),
            (ra=350.0, dec=-72.0, plx=8.0, pmra=25.0, pmdec=-15.0, rv=-5e3),
            (ra=0.0, dec=0.0, plx=50.0, pmra=0.0, pmdec=0.0, rv=0.0),
        )
        for c in cases
            fr = PlanetOrbits.AbsoluteFrame(; c..., ref_epoch=57388.5)
            # ±100 yr, well beyond any real epoch span, and straddling the
            # reference epoch to exercise the `ref_epoch == t` guard.
            for t in (fr.ref_epoch, range(fr.ref_epoch - 36525, fr.ref_epoch + 36525,
                                          length=41)...)
                got = PlanetOrbits.compensate(fr, t)
                ref = compensate_literal(fr, t)
                # Both forms evaluate the same mathematics, so the only honest
                # gate is "agrees to a few ulps". Two conditioning regimes:
                #
                # ra2/dec2/rv2/distance2 are well conditioned — ≤ 21 ulps
                # (≤ 4e-15 relative) across every case here — so 1e-13 relative
                # is already ~30x of headroom over what either form can produce.
                for f in (:ra2, :dec2, :rv2, :distance2)
                    a = getfield(got, f)
                    b = getfield(ref, f)
                    @test a ≈ b rtol = 1e-13 atol = 1e-13
                end
                # pmra2/pmdec2 are not: ∂δ/∂t is
                #     (-z2*ddist2/distance2 + fr.dz) / (distance2*√(1 - z2²/distance2²))
                # and in the high-declination, fast, nearly-radially-moving case
                # (ra=10, dec=87.5, rv=60 km/s) both halves of that numerator
                # agree to ~3 digits, so last-ulp differences in `distance2` are
                # amplified ~1.5e3x. Measured worst case is 1489 ulps / 1.9e-13
                # relative on 1.10 (8.7e-14 on 1.12), which is cancellation, not
                # disagreement — 1e-12 keeps ~5x margin over it while staying
                # four orders tighter than any physically meaningful slop.
                for f in (:pmra2, :pmdec2)
                    a = getfield(got, f)
                    b = getfield(ref, f)
                    @test a ≈ b rtol = 1e-12 atol = 1e-13
                end
                # t_em feeds the Kepler solve directly, so hold it in absolute
                # days rather than relatively — but scaled to the resolution
                # that actually exists there. These are MJDs of magnitude ~5e4,
                # where eps() is already 7.3e-12 d: the previous fixed
                # atol=1e-13 was a hundredth of one ulp, i.e. demanding
                # bit-identity while claiming not to, and it held on 1.12 only
                # by rounding luck. Measured worst case is 3 ulps.
                @test isapprox(got.epoch2a_days, ref.epoch2a_days;
                               atol=8 * eps(ref.epoch2a_days))
            end
        end
    end

    @testset "allocation-free and inferred at NB ≥ 5, Float64 and Duals" begin
        # Per §12: allocation gates that only ever used a 2-body system missed
        # a cliff that appeared from NB = 3 up, and only under Duals.
        function mk(T, NB)
            A = Body(mass=T(1.0), name=:A)
            bs = ntuple(i -> Body(mass=T(1e-3), name=Symbol("b", i)), NB - 1)
            orbs = ntuple(i -> Orbit(bs[i], about=A; a=T(1.0 + i), e=T(0.1), i=T(0.5),
                                     ω=T(1.1), Ω=T(2.2), tp=T(57000.0)), NB - 1)
            System((A, bs...), orbs; plx=T(50.0), ra=T(45.0), dec=T(20.0),
                   pmra=T(200.0), pmdec=T(150.0), rv=T(-30e3), ref_epoch=T(57388.0))
        end
        eps = collect(range(57388.0, 58000.0, length=64))
        for NB in (2, 5), T in (Float64, ForwardDiff.Dual{Nothing,Float64,6})
            sys = mk(T, NB)
            traj = Trajectory{T}(sys, eps)
            orbitsolve!(traj, sys)                     # warm up
            if DYNAMIC_ALLOC_GATE
                @test (@allocated orbitsolve!(traj, sys)) == 0
            else
                @test_skip (@allocated orbitsolve!(traj, sys)) == 0
            end
            @test (@inferred observe_pass!(traj, sys)) isa Trajectory
        end
    end

    @testset "differential light-travel time is present and scales as z/c" begin
        # Face-on: z ≡ 0, so retardation is identically inert. Edge-on: the
        # apparent phase lag is z/c, so the position differs from the
        # unretarded state by v·z/c.
        A = Body(mass=1.0, name=:A)
        b = Body(mass=0.0, name=:b)
        for (i, retarded) in ((0.0, false), (π / 2, true))
            orb = Orbit(b, about=A; a=5.0, e=0.0, i=i, ω=0.0, Ω=0.0, tp=57000.0)
            sys = System((A, b), (orb,))
            sysraw = sys
            t = 57900.0
            sol = orbitsolve(sys, [t])[1]
            r, v = _rawstate(sysraw, t)
            rA, rb = bodies(sys)
            Δ = hypot(posx(sol, rb, rA) - (r[2][1] - r[1][1]),
                      posy(sol, rb, rA) - (r[2][2] - r[1][2]))
            if retarded
                zb = r[2][3] - r[1][3]
                vb = hypot(v[2][1] - v[1][1], v[2][2] - v[1][2])
                @test Δ ≈ abs(vb * zb / c_au_per_julianyr) rtol = 1e-3
                @test Δ > 0
            else
                @test Δ < 1e-15
            end
        end
    end
end

# ---------------------------------------------------
# The Einstein term in `radvel` (src/observe.jl, src/observables.jl).
#
# `radvel` is the spectroscopic velocity: the kinematic projection plus the
# second-order Doppler and gravitational-redshift difference between its two
# references. The magnitudes below are the ones tabulated on the "Precision
# opt-outs" manual page, so this testset is also what keeps that page honest.
# ---------------------------------------------------

# The Einstein term of a pair, isolated: `radvel` minus the kinematic
# quantity in the same units.
einstein(sol, t, r) =
    radvel(sol, t, r) - velz(sol, t, r) * PlanetOrbits.au2m / year2sec_julian

@testset "radvel Einstein term" begin
    A = Body(mass=1.0, name=:A)

    @testset "magnitudes match the documented table" begin
        # Circular hot Jupiter, 0.05 AU: ~89 m/s, and constant.
        b = Body(mass=1mjup, name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=0.05, e=0.0, i=0.7, ω=0.3,
                                    Ω=1.1, tp=57000.0),); plx=10.0)
        eps = collect(range(57000.0, 57030.0, length=41))
        traj = orbitsolve(sys, eps)
        e = [einstein(traj[k], :b, :A) for k in eachindex(eps)]
        @test all(x -> isapprox(x, 89.0; rtol=0.02), e)
        @test maximum(e) - minimum(e) < 1e-6           # constant to << 1 mm/s

        # 1 AU, e = 0.4 Jupiter: 2.7–8.4 m/s, varying by ~5.6 m/s.
        sys2 = System((A, b), (Orbit(b, about=A; a=1.0, e=0.4, i=0.7, ω=0.3,
                                     Ω=1.1, tp=57000.0),); plx=10.0)
        eps2 = collect(range(57000.0, 57000.0 + 365.25, length=201))
        traj2 = orbitsolve(sys2, eps2)
        e2 = [einstein(traj2[k], :b, :A) for k in eachindex(eps2)]
        @test minimum(e2) ≈ 2.7 rtol = 0.05
        @test maximum(e2) ≈ 8.4 rtol = 0.05
        @test maximum(e2) - minimum(e2) ≈ 5.6 rtol = 0.05

        # The same system's stellar reflex is three orders of magnitude
        # smaller and sub-cm/s in its variation — the asymmetry the manual
        # prices, and the reason the two uses of `radvel` are documented
        # separately.
        bary = barycentre(sys2)
        r2 = [einstein(traj2[k], :A, bary) for k in eachindex(eps2)]
        @test maximum(abs, r2) < 1e-2                        # < 1 cm/s
        @test (maximum(e2) - minimum(e2)) / (maximum(r2) - minimum(r2)) > 1e3
    end

    @testset "nothing emits from a barycentre" begin
        b = Body(mass=0.3, name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=1.0, e=0.1, i=0.7, ω=0.3,
                                    Ω=1.1, tp=57000.0),); plx=10.0)
        sol = orbitsolve(sys, [57100.0])[1]
        bary = barycentre(sys)
        rA, rb = bodies(sys)
        # A barycentre's Einstein term is identically zero, so the reflex
        # case carries the star's own term in full rather than differencing
        # it against a mass-weighted average.
        @test PlanetOrbits._ein(sol, bary) == 0
        @test PlanetOrbits._ein(sol, framedirection) == 0
        @test PlanetOrbits._ein(sol, rA) != 0
        @test einstein(sol, rA, bary) ≈ PlanetOrbits._ein(sol, rA) *
                                        PlanetOrbits.au2m / year2sec_julian rtol = 1e-8
        # A photocentre does emit: it is the flux-weighted blend.
        Af = Body(mass=1.0, name=:A, flux=(; G=1.0))
        bf = Body(mass=0.3, name=:b, flux=(; G=0.25))
        sysf = System((Af, bf), (Orbit(bf, about=Af; a=1.0, e=0.1, i=0.7, ω=0.3,
                                       Ω=1.1, tp=57000.0),); plx=10.0)
        solf = orbitsolve(sysf, [57100.0])[1]
        rAf, rbf = bodies(sysf)
        p = photocentre(sysf)
        @test PlanetOrbits._ein(solf, p) ≈
              (1.0 * PlanetOrbits._ein(solf, rAf) +
               0.25 * PlanetOrbits._ein(solf, rbf)) / 1.25 rtol = 1e-13
    end

    @testset "v_sys · v_orb / c cross term needs an absolute frame" begin
        # An SB-like pair: 30 km/s systemic, tens of km/s reflex. `v_tot` is
        # the body's total barycentric velocity, so the cross term is present
        # with an AbsoluteFrame and definitionally absent with a Parallax.
        b = Body(mass=0.8, name=:b)
        orb = Orbit(b, about=A; a=0.5, e=0.0, i=π / 2, ω=0.0, Ω=0.0, tp=57000.0)
        base = (; plx=10.0)
        abs_sys = System((A, b), (orb,); base..., ra=45.0, dec=20.0,
                         pmra=0.0, pmdec=0.0, rv=30e3, ref_epoch=57000.0)
        plx_sys = System((A, b), (orb,); base...)
        eps = collect(range(57000.0, 57000.0 + 130.0, length=81))
        ea = [einstein(orbitsolve(abs_sys, eps)[k], :A, barycentre(abs_sys))
              for k in eachindex(eps)]
        ep = [einstein(orbitsolve(plx_sys, eps)[k], :A, barycentre(plx_sys))
              for k in eachindex(eps)]
        # Constant offset v_sys²/2c aside, the *varying* part is the cross
        # term and it is at the m/s level.
        @test (maximum(ea) - minimum(ea)) - (maximum(ep) - minimum(ep)) > 1.0
        @test (maximum(ep) - minimum(ep)) < 0.5
    end

    @testset "the geometry flag may not change what radvel means" begin
        # `radvel` legitimately differs across `observing_geometry` — the
        # kinematic line-of-sight projection is one of the four corrections.
        # What may *not* differ is the Einstein term itself: if the skipped
        # path left it out, the flag would change the observable's meaning
        # rather than its precision. So gate the Einstein delta, and require
        # it to agree to the size of the geometry corrections themselves.
        b = Body(mass=1mjup, name=:b)
        orb = Orbit(b, about=A; a=1.0, e=0.4, i=0.7, ω=0.3, Ω=1.1, tp=57000.0)
        sys = System((A, b), (orb,); plx=546.5, ra=45.0, dec=20.0,
                     pmra=7323.0, pmdec=7323.0, rv=-110e3, ref_epoch=57388.0)
        eps = collect(range(57388.0, 57388.0 + 365.25, length=51))
        full = orbitsolve(sys, eps)
        cheap = orbitsolve(sys, eps; observing_geometry=false)
        de = maximum(abs(einstein(full[k], :b, :A) - einstein(cheap[k], :b, :A))
                     for k in eachindex(eps))
        dr = maximum(abs(radvel(full[k], :b, :A) - radvel(cheap[k], :b, :A))
                     for k in eachindex(eps))
        @test de < dr                       # the Einstein term is not the difference
        @test de < 1e-2                      # ...and is below a cm/s here
        @test dr > 1.0                       # while the kinematic part is m/s
        # It is genuinely filled on the cheap path, not zeroed.
        @test maximum(abs(einstein(cheap[k], :b, :A)) for k in eachindex(eps)) > 1.0
        # Every frame level fills it.
        for kw in ((;), (; plx=10.0))
            s = System((A, b), (orb,); kw...)
            for og in (true, false)
                so = orbitsolve(s, [57500.0]; observing_geometry=og)[1]
                @test abs(einstein(so, :b, :A)) > 1.0
            end
        end
    end

    @testset "masses reach the RV gradient" begin
        b = Body(mass=1mjup, name=:b)
        f = function (m)
            bb = Body(mass=m, name=:b)
            s = System((A, bb), (Orbit(bb, about=A; a=1.0, e=0.4, i=0.7, ω=0.3,
                                       Ω=1.1, tp=57000.0),); plx=10.0)
            radvel(orbitsolve(s, [57100.0])[1], :b, :A)
        end
        g = ForwardDiff.derivative(f, 1mjup)
        @test isfinite(g)
        @test g ≈ FiniteDiff.finite_difference_derivative(f, 1mjup) rtol = 1e-5
        # And it is the Einstein term that put it there for the *relative* RV
        # of a massless companion: the kinematic part of `radvel(b, A)` does
        # depend on m through the row's total mass too, so compare against
        # the derivative of the kinematic quantity, which must differ.
        fk = function (m)
            bb = Body(mass=m, name=:b)
            s = System((A, bb), (Orbit(bb, about=A; a=1.0, e=0.4, i=0.7, ω=0.3,
                                       Ω=1.1, tp=57000.0),); plx=10.0)
            sol = orbitsolve(s, [57100.0])[1]
            velz(sol, :b, :A) * PlanetOrbits.au2m / year2sec_julian
        end
        @test !isapprox(g, ForwardDiff.derivative(fk, 1mjup); rtol=1e-6)
    end
end

# ---------------------------------------------------
# Observer-aware observables (src/observables.jl).
# ---------------------------------------------------

@testset "observer-aware observables" begin
    ref_epoch = 57388.0
    A = Body(mass=1.0, name=:A)
    b = Body(mass=1e-6, name=:b)
    # No proper motion and read at `ref_epoch`, so the epoch triad is the
    # reference triad and an ICRS observer position built from M1 has known
    # components in it.
    function obscase(; d=10.0, a=10.0, i=π / 2)
        orb = Orbit(b, about=A; a, e=0.0, i, ω=0.0, Ω=0.0, tp=57000.0)
        System((A, b), (orb,); plx=1000 / d, ra=45.0, dec=20.0, pmra=0.0,
               pmdec=0.0, rv=0.0, ref_epoch)
    end

    @testset "an observer at the SSB is the zero-argument form" begin
        sys = obscase()
        sol = orbitsolve(sys, [ref_epoch])[1]
        bary = barycentre(sys)
        # Not bit-identical only because the shared `cart2angle` factor is
        # rounded once and reciprocated here — same expression, same value to
        # the last few ulp.
        same(x, y) = isapprox(x, y; rtol=1e-14, atol=1e-25)
        for (t, r) in ((:b, :A), (:A, bary), (:A, framedirection))
            @test same(raoff(sol, t, r, (0.0, 0.0, 0.0)), raoff(sol, t, r))
            @test same(decoff(sol, t, r, (0.0, 0.0, 0.0)), decoff(sol, t, r))
        end
        @test same(projectedseparation(sol, :b, :A, SVector(0.0, 0.0, 0.0)),
                   projectedseparation(sol, :b, :A))
        @test same(posangle(sol, :b, :A, SVector(0.0, 0.0, 0.0)),
                   posangle(sol, :b, :A))
    end

    @testset "against a direction, the full parallax ellipse" begin
        d = 10.0
        sys = obscase(; d)
        sol = orbitsolve(sys, [ref_epoch])[1]
        M1 = sys.frame.M1
        # 1 AU East and 1 AU North, in ICRS.
        east = M1 * SVector(1.0, 0.0, 0.0)
        north = M1 * SVector(0.0, 1.0, 0.0)
        @test raoff(sol, :A, framedirection, east) - raoff(sol, :A, framedirection) ≈
              -1000 / d rtol = 1e-9
        @test decoff(sol, :A, framedirection, east) - decoff(sol, :A, framedirection) ≈
              0 atol = 1e-12
        @test decoff(sol, :A, framedirection, north) - decoff(sol, :A, framedirection) ≈
              -1000 / d rtol = 1e-9
        # Against the barycentre — a point at the *same* distance — the same
        # displacement leaves only the differential part, five orders down.
        bary = barycentre(sys)
        @test abs(raoff(sol, :A, bary, east) - raoff(sol, :A, bary)) < 1e-4
    end

    @testset "Kopeikin: differential parallax survives a relative pair" begin
        # Δθ ≈ 4.85 · z[AU] / d[pc]² µas per AU of observer displacement,
        # with z the *line-of-sight* separation — so it is there for an
        # edge-on orbit and absent for a face-on one of the same size.
        for d in (5.0, 10.0)
            sys = obscase(; d, a=10.0, i=π / 2)
            sol = orbitsolve(sys, [ref_epoch])[1]
            east = sys.frame.M1 * SVector(1.0, 0.0, 0.0)
            Δ = (raoff(sol, :b, :A, east) - raoff(sol, :b, :A)) * 1000   # µas
            @test Δ ≈ 4.84814 * posz(sol, :b, :A) / d^2 rtol = 1e-3
            @test abs(Δ) > 0.05
            # First-order parallax has cancelled: what remains is smaller
            # than the parallax factor itself by ~z/d.
            @test abs(Δ) < 1e-3 * 1000 * (1000 / d)
        end
        face = obscase(; d=10.0, a=10.0, i=0.0)
        solf = orbitsolve(face, [ref_epoch])[1]
        eastf = face.frame.M1 * SVector(1.0, 0.0, 0.0)
        @test abs(raoff(solf, :b, :A, eastf) - raoff(solf, :b, :A)) < 1e-12
    end

    @testset "first-order parallax factors, exactly" begin
        # The migration path for a likelihood that today consumes SSB
        # observables plus hand-rolled parallax factors (Hipparcos IAD, HST
        # FGS): the observer-aware read must reproduce those factors at first
        # order, and differ from them only by the depth/Kopeikin terms the
        # factor math cannot express.
        #
        # No ephemeris is needed to check that — a circular Earth in the
        # ecliptic is enough, since the claim is about the projection, not
        # about where the Earth actually is.
        d = 30.0
        sys = obscase(; d, a=2.0, i=0.9)
        ε = deg2rad(23.4392911)
        # ICRS position of a 1 AU circular observer in the ecliptic plane.
        earth(λ) = SVector(cos(λ), sin(λ) * cos(ε), sin(λ) * sin(ε))
        eA = sys.frame.M1[:, 1]     # ê_α, ICRS
        eD = sys.frame.M1[:, 2]     # ê_δ
        sol = orbitsolve(sys, [ref_epoch])[1]
        for λ in range(0, 2π; length=13)[1:end-1]
            o = earth(λ)
            # The classical first-order factors: the observer's transverse
            # displacement over the distance, in mas.
            fα = -(o ⋅ eA) * (1000 / d)
            fδ = -(o ⋅ eD) * (1000 / d)
            Δα = raoff(sol, :A, framedirection, o) - raoff(sol, :A, framedirection)
            Δδ = decoff(sol, :A, framedirection, o) - decoff(sol, :A, framedirection)
            @test Δα ≈ fα rtol = 1e-6 atol = 1e-9
            @test Δδ ≈ fδ rtol = 1e-6 atol = 1e-9
            # ...and the residual is not zero: it is the second-order part,
            # which is exactly what the exact geometry buys.
            @test 0 < abs(Δα - fα) < 1e-4        # mas, i.e. sub-0.1 µas here
        end
    end

    @testset "error paths" begin
        sys = obscase()
        east = sys.frame.M1 * SVector(1.0, 0.0, 0.0)
        cheap = orbitsolve(sys, [ref_epoch]; observing_geometry=false)[1]
        @test_throws ErrorException raoff(cheap, :b, :A, east)
        @test_throws ErrorException decoff(cheap, :b, :A, east)
        # A Parallax (or frameless) system has no ICRS direction to place an
        # ICRS observer against.
        orb = Orbit(b, about=A; a=10.0, e=0.0, i=π / 2, ω=0.0, Ω=0.0, tp=57000.0)
        for kw in ((; plx=100.0), (;))
            s = System((A, b), (orb,); kw...)
            so = orbitsolve(s, [ref_epoch])[1]
            @test_throws ErrorException raoff(so, :b, :A, east)
        end
    end

    @testset "allocation-free and inferred" begin
        sys = obscase()
        traj = orbitsolve(sys, [ref_epoch, ref_epoch + 100])
        sol = traj[1]
        o = SVector(0.3, -0.8, 0.1)
        rA, rb = bodies(sys)
        @test (@inferred raoff(sol, rb, rA, o)) isa Float64
        @test (@inferred decoff(sol, rb, rA, o)) isa Float64
        @test (@inferred raoff(sol, rA, framedirection, o)) isa Float64
        f(traj, rb, rA, o) = sum(raoff(traj[k], rb, rA, o) for k in eachindex(traj))
        f(traj, rb, rA, o)
        if DYNAMIC_ALLOC_GATE
            @test (@allocated f(traj, rb, rA, o)) == 0
        else
            @test_skip (@allocated f(traj, rb, rA, o)) == 0
        end
    end
end
