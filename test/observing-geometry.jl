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
                    pc2sec_light, au2pc, pc2au, pc2m, rad2mas,
                    year2day_julian, year2sec_julian, c_au_per_julianyr

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
    V = (fr.pmra / 1000 * AS2RAD * d1) * ea +            # pc / julian yr
        (fr.pmdec / 1000 * AS2RAD * d1) * ed +
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
        # period must be LONGER than its true period. v1 had this inverted.
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
        Δ = radvel(sol_pm, rb, rA) - radvel(sol_0, rb, rA)
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

    @testset "observing_geometry=false is exactly the v1 / stage-1 geometry" begin
        # The opt-out selects the *cheap* geometry, not a different physics:
        # one shared AU->mas scale per epoch, no rotation, no retardation, no
        # line-of-sight projection. For the frameless and parallax-only
        # fixtures that is bit-for-bit what v1 computed, which is what proves
        # the observing pass is the only thing that changed for them.
        for c in V1_REFERENCE
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
        # entirely the barycentric light-travel *sign* fix — the one v1 bug
        # that is not part of the observing pass and so is not opted out of.
        for c in V1_REFERENCE
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
        @test (@allocated orbitsolve!(traj, sys; observing_geometry=false)) == 0
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
            @test (@allocated orbitsolve!(traj, sys)) == 0
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
