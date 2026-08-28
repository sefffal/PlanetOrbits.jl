# Barycentric light time against a rigorous Butkevich & Lindegren (2014)
# oracle.
#
# The load-bearing references here are built from full 3D ICRS vectors and
# share no propagation code with production:
#
#   * `bl_rigorous`  — the rigorous *apparent* place and proper motion: the
#     catalog (apparent, d/dt_obs) proper motions are de-Dopplered to the true
#     space velocity, μ_true = μ_app·(1 + v_r/c) (B&L Eqs. 11–13; ERFA
#     `starpv`'s observed→inertial step), the emission parameter τ(t) is
#     solved from τ + (d(τ) − d_ref)/c = t, the direction is read at τ, and
#     the rates are converted back to apparent, d/dt_obs = (d/dτ)/(1 + ḋ/c).
#     This is what `barycentric_lighttime=true` must produce.
#
#   * `bl_standard`  — the light-time-free standard model the astrometric
#     catalogs are reduced with (ESA 1997 Sect. 1.5.5; Lindegren et al.
#     2012/2021): the catalog values propagated linearly as they stand.
#     This is what `barycentric_lighttime=false` must produce.
#
# Both conventions reproduce the catalog values identically at `ref_epoch` and
# differ away from it only by the genuine second-order light-time terms.
#
# History this guards against (fixed together; each alone merely moves the
# error between readouts): building the light-time worldline from apparent
# proper motions *without* de-Dopplering double-counts μ·v_r/c in every
# position readout, while `_compensate_kinematics`' coordinate rates (d/dτ)
# accidentally cancel it in the PM readouts — an internal position-vs-PM
# inconsistency of μ·|v_r|/c (3.8 mas/yr at Barnard's star) that a fit can
# only resolve by inventing a companion. The internal-consistency testset at
# the bottom detects that whole bug class without reference to any oracle.

using LinearAlgebra
using PlanetOrbits: frame_ra, frame_dec, frame_pmra, frame_pmdec, frame_rv,
                    year2day_julian, year2sec_julian, pc2km, c_light_ms

@testset "barycentric light time: B&L 2014 oracle" begin
    c_kms = c_light_ms / 1000
    c_pc_yr = c_kms * year2sec_julian / pc2km
    rad2mas_ = 180 / π * 3600 * 1000

    bl_triad(ra_deg, dec_deg) = begin
        sra, cra = sincosd(ra_deg)
        sdec, cdec = sincosd(dec_deg)
        ([-sra, cra, 0.0], [-cra * sdec, -sra * sdec, cdec],
         [cra * cdec, sra * cdec, sdec])
    end

    function bl_worldlines(s)
        E0, N0, R0 = bl_triad(s.ra, s.dec)
        d0 = 1000 / s.plx
        vr_pcyr = s.rv / 1000 * year2sec_julian / pc2km
        β = s.rv / 1000 / c_kms
        vel(f) = f * d0 * s.pmra / rad2mas_ .* E0 .+
                 f * d0 * s.pmdec / rad2mas_ .* N0 .+ vr_pcyr .* R0
        return (b0=d0 .* R0, v_true=vel(1 + β), v_tilde=vel(1.0), d0=d0)
    end

    function bl_state(w, v, τ_yr)
        b = w.b0 .+ v .* τ_yr
        d = norm(b)
        ra = atand(b[2], b[1])
        ra < 0 && (ra += 360)
        dec = asind(clamp(b[3] / d, -1, 1))
        e, n, _ = bl_triad(ra, dec)
        (; ra, dec, pmra=dot(v, e) / d * rad2mas_,
           pmdec=dot(v, n) / d * rad2mas_, ddot=dot(b, v) / d)
    end

    function bl_tau(w, v, t_yr)
        τ = t_yr
        for _ in 1:8
            τ = t_yr - (norm(w.b0 .+ v .* τ) - w.d0) / c_pc_yr
        end
        return τ
    end

    # The apparent factor goes on the angular rates only; `rv` stays the
    # coordinate rate ḋ at the emission event, which is what a Doppler
    # measurement of that event corresponds to. See the convention note in
    # `_frame_pass_kernel!`.
    function bl_rigorous(w, t_yr)
        st = bl_state(w, w.v_true, bl_tau(w, w.v_true, t_yr))
        s = 1 / (1 + st.ddot / c_pc_yr)
        (; st.ra, st.dec, pmra=st.pmra * s, pmdec=st.pmdec * s,
           rv=st.ddot * pc2km / year2sec_julian * 1e3)
    end

    function bl_standard(w, t_yr)
        st = bl_state(w, w.v_tilde, t_yr)
        (; st.ra, st.dec, st.pmra, st.pmdec,
           rv=st.ddot * pc2km / year2sec_julian * 1e3)
    end

    # Gaia DR3 values. Barnard's star is B&L's own worked example and the
    # worst catalog case there is (nearest single star, fastest PM, |v_r| =
    # 110 km/s: μ·|v_r|/c = 3.83 mas/yr). Kapteyn's star is the fast *receding*
    # counterpart with its perspective acceleration far from parallel to μ;
    # ups And is a garden-variety planet host; the last case has v_r = 0, where
    # the de-Doppler factor is exactly the identity.
    stars = (
        (name="Barnard", ra=269.44850252543836, dec=4.739420051112412,
         plx=546.975939730948, pmra=-801.5509783684709,
         pmdec=10362.394206546573, rv=-110468.22357177734),
        (name="Kapteyn", ra=77.9599373502188, dec=-45.0438126993602,
         plx=254.19859326384577, pmra=6491.223339061598,
         pmdec=-5708.614150045243, rv=245053.13110351562),
        (name="ups And", ra=24.198321215249642, dec=41.4037617615712,
         plx=74.19396780666189, pmra=-171.89184082707104,
         pmdec=-381.8151588992299, rv=-28745.87059020996),
        (name="v_r = 0", ra=100.0, dec=-25.0, plx=200.0,
         pmra=1500.0, pmdec=-2500.0, rv=0.0),
    )
    ref_epoch = 57388.0   # J2016.0 (TCB), MJD
    ep_off = [-24.75, -12.3, -1.4, 0.0, 1.4, 10.0]

    mksys(s) = let A = Body(mass=0.3, name=:A), b = Body(mass=0.0, name=:b)
        System((A, b), (Orbit(b, about=A; a=2.0, e=0.0, i=0.0, ω=0.0, Ω=0.0,
                              tp=ref_epoch),);
               ra=s.ra, dec=s.dec, plx=s.plx, pmra=s.pmra, pmdec=s.pmdec,
               rv=s.rv, ref_epoch=ref_epoch)
    end

    for s in stars
        @testset "$(s.name)" begin
            w = bl_worldlines(s)
            sys = mksys(s)
            epochs = ref_epoch .+ year2day_julian .* ep_off

            # Both flag settings against their own oracle. The limits sit well
            # above the measured agreement (≤2e-7 mas, ≤2e-11 mas/yr,
            # ≤1e-10 m/s — the floor is double-precision cancellation in
            # differencing two ~270° angles, so it deserves room to move
            # between machines) and far below anything physical: the failure
            # this guards against is the μ·v_r/c convention error, five to
            # seven orders of magnitude larger at 94 mas / 3.8 mas/yr.
            for (lt, oracle) in ((true, bl_rigorous), (false, bl_standard))
                traj = orbitsolve(sys, epochs; barycentric_lighttime=lt)
                for (k, ep) in enumerate(epochs)
                    o = oracle(w, (ep - ref_epoch) / year2day_julian)
                    sol = traj[k]
                    @test abs(frame_ra(sol) - o.ra) * 3.6e6 * cosd(o.dec) < 1e-5
                    @test abs(frame_dec(sol) - o.dec) * 3.6e6 < 1e-5
                    @test abs(frame_pmra(sol) - o.pmra) < 1e-8
                    @test abs(frame_pmdec(sol) - o.pmdec) < 1e-8
                    @test abs(frame_rv(sol) - o.rv) < 1e-7
                end
            end

            # Catalog identity at ref_epoch, under BOTH conventions: the
            # anchor of the whole construction (B&L Sect. 1 — the catalog
            # values are apparent quantities at their own epoch).
            for lt in (true, false)
                sol = orbitsolve(sys, [ref_epoch]; barycentric_lighttime=lt)[1]
                @test abs(frame_ra(sol) - s.ra) * 3.6e6 * cosd(s.dec) < 1e-5
                @test abs(frame_dec(sol) - s.dec) * 3.6e6 < 1e-5
                @test abs(frame_pmra(sol) - s.pmra) < 1e-8
                @test abs(frame_pmdec(sol) - s.pmdec) < 1e-8
                # Including the radial one: `rv2` carries no apparent factor,
                # so this identity holds under BOTH conventions — which is
                # what lets a consumer compose `frame_rv` against a catalog
                # radial velocity without knowing which way the flag went.
                @test abs(frame_rv(sol) - s.rv) < 1e-6
            end

            # Internal consistency of the apparent path: the numeric
            # d/dt_obs of the position readouts must BE the PM readouts.
            # This is the oracle-free detector for the whole convention-error
            # bug class: with either half of the fix missing it fails at
            # μ·|v_r|/c (3.8 mas/yr at Barnard), against a finite-difference
            # truncation floor ~1e-5 mas/yr.
            h = 2.0  # days
            for dt in (-24.65, 0.17, 8.0)
                ep = ref_epoch + dt * year2day_julian
                tr = orbitsolve(sys, [ep - h, ep, ep + h];
                                barycentric_lighttime=true)
                num_pmra = (frame_ra(tr[3]) - frame_ra(tr[1])) *
                           cosd(frame_dec(tr[2])) * 3.6e6 / (2h / year2day_julian)
                num_pmdec = (frame_dec(tr[3]) - frame_dec(tr[1])) * 3.6e6 /
                            (2h / year2day_julian)
                @test abs(num_pmra - frame_pmra(tr[2])) < 1e-3
                @test abs(num_pmdec - frame_pmdec(tr[2])) < 1e-3
            end

            # The two conventions differ ONLY by the genuine second-order
            # light-time terms: ≤0.2 mas / ≤0.02 mas/yr over ±25 yr even at
            # these extreme kinematics (B&L Sect. 5.5's negligibility scale).
            # This bound is what lets a *catalog-convention* consumer treat
            # the flag as harmless either way — and for v_r ≠ 0 the
            # difference must also be genuinely present, so the apparent path
            # cannot silently degrade into the standard model.
            on = orbitsolve(sys, epochs; barycentric_lighttime=true)
            off = orbitsolve(sys, epochs; barycentric_lighttime=false)
            worst_pos = maximum(
                max(abs(frame_ra(on[k]) - frame_ra(off[k])) * 3.6e6 * cosd(s.dec),
                    abs(frame_dec(on[k]) - frame_dec(off[k])) * 3.6e6)
                for k in eachindex(epochs))
            worst_pm = maximum(
                max(abs(frame_pmra(on[k]) - frame_pmra(off[k])),
                    abs(frame_pmdec(on[k]) - frame_pmdec(off[k])))
                for k in eachindex(epochs))
            @test worst_pos < 0.2       # mas
            @test worst_pm < 0.02       # mas/yr
            if !iszero(s.rv) && hypot(s.pmra, s.pmdec) > 1000
                @test worst_pos > 1e-3
            end
        end
    end
end
