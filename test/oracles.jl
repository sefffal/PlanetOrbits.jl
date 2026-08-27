# ---------------------------------------------------
# External oracles — release gates.
#
# Everything else in test/ checks PlanetOrbits against *itself*: round trips,
# invariants, one code path against another, AD against finite differences.
# That is the right way to catch a regression, and it cannot catch a value
# that has been wrong since the first commit. These testsets exist for the
# other failure mode: each one compares v2 against a number PlanetOrbits did
# not produce.
#
# Three oracles, chosen for coverage per unit of trust required:
#
#   1. PlanetOrbits v1 (0.11.4) — an independently written implementation of
#      the same Campbell-element surface. Covers the *whole* pipeline
#      (elements -> Kepler solve -> Cartesian state -> observables) over a
#      30-case parameter sweep. Values are pre-generated into
#      fixtures/v1_oracle.jl so CI needs no v1 install; see
#      fixtures/generate_v1_oracle.jl to regenerate.
#
#   2. Kepler's equation in 256-bit BigFloat, solved by bisection on a bracket
#      where the residual is provably monotone, then Newton-polished. Covers
#      the solvers alone, but to a precision and over a range of (M, e) that
#      no cross-implementation check can reach. Runs live — it is cheap.
#
#   3. Closed-form and published values: IAU nominal constants, the definition
#      of the parsec, the textbook RV semi-amplitude, and the Thiele-Innes
#      constants. Each is written out here from its source rather than read
#      back out of the package.
#
# Deliberately NOT (re)done here:
#
#   * The AHL21 N-body propagator already has an external oracle — the
#     NbodyGradient v0.2.3 trajectories in fixtures/nbody_reference.jl — plus
#     two-body-vs-KeplerianApprox and symplectic energy/angular-momentum
#     gates. See test/nbody.jl. Adding a fourth N-body comparison here would
#     be duplication, not coverage.
#   * v1 as an oracle for hyperbolic orbits or for the absolute frame. v1 is
#     known-wrong on both (zero velocity semiamplitude for e > 1; two
#     different speeds of light in the light-travel term). See the header of
#     generate_v1_oracle.jl.
#
# Tolerances below are the measured worst case rounded up to a round number,
# with the margin stated. If one of these starts failing, the number that
# moved is the interesting thing — do not widen the tolerance.
# ---------------------------------------------------

include("fixtures/v1_oracle.jl")

using PlanetOrbits: au2m, year2sec_julian, year2day_julian, pc2au, rad2mas,
                    kepler_year_to_julian_day_conversion_factor,
                    kepler_solver, Markley, Auto, Goat, HyperbolicHalley,
                    RootsMethod, markley_sincosE, VREM2PI_MAX

# v1's `radvel` is the kinematic line-of-sight velocity; v2's adds the
# Einstein term (see the note by `kinrv` in runtests.jl). Same definition,
# repeated so this file reads standalone.
_orc_kinrv(sol, refs...) = velz(sol, refs...) * au2m / year2sec_julian

# ---------------------------------------------------
# Oracle 1 — PlanetOrbits v1
# ---------------------------------------------------
#
# Read through `observing_geometry=false`, which src/observe.jl documents as
# selecting "exactly what v1 ... computed": one shared AU->mas scale per epoch
# at the barycentre's distance, no viewing-direction rotation, no per-body
# light-travel retardation, no line-of-sight projection. Asked for the same
# physical model, two independent implementations must agree to roundoff, and
# they do.
#
# This is a much sharper instrument than the "v1 regression fixtures" testset
# in runtests.jl, which reads the *default* geometry and is therefore pinned
# at rtol 3e-3. Both are kept: that one gates the corrections' overall size,
# this one gates the shared math.

@testset "oracle: PlanetOrbits v1 cross-check" begin
    # Measured worst deviation across all 30 cases x 11 epochs, relative to
    # each case's own characteristic scale: 4.9e-14 (`near-parabolic`, e=0.999).
    # Gated at 1e-12, a 20x margin.
    RTOL = 1e-12
    # Proper motion is the exception and is handled separately below.
    pm_worst_inclined = 0.0

    for c in V1_ORACLE
        @testset "$(c.name)" begin
            p = c.params
            A = Body(mass=p.M - c.Mp, name=:A)
            b = Body(mass=c.Mp, name=:b)
            sys = System(Orbit(b, about=A; a=p.a, e=p.e, i=p.i, ω=p.ω, Ω=p.Ω, tp=p.tp);
                         plx=p.plx)
            refs = bodies(sys)
            bary = barycentre(sys)
            traj = orbitsolve(sys, c.epochs; observing_geometry=false)
            d = c.data

            # Kepler's third law, through two independent code paths.
            @test period(sys) ≈ c.period rtol = RTOL

            # Per-case scales, so that a quantity passing through zero (posz
            # on a face-on orbit, velx at a turning point) is compared against
            # the size of the thing it belongs to rather than against itself.
            s_pos = maximum(abs, d.posx)
            s_vel = maximum(abs, d.velx)
            s_rv = s_vel * au2m / year2sec_julian
            s_ang = maximum(abs, d.raoff)
            s_angr = maximum(abs, d.raoff_reflex)
            s_rvr = maximum(abs, d.radvel_reflex)

            ap(got, want, scale) = isapprox(got, want; rtol=RTOL, atol=RTOL * scale)

            for (k, sol) in enumerate(traj)
                # --- barycentric-relative state, AU and AU/julian yr --------
                @test ap(posx(sol, refs.b, refs.A), d.posx[k], s_pos)
                @test ap(posy(sol, refs.b, refs.A), d.posy[k], s_pos)
                @test ap(posz(sol, refs.b, refs.A), d.posz[k], s_pos)
                @test ap(velx(sol, refs.b, refs.A), d.velx[k], s_vel)
                @test ap(vely(sol, refs.b, refs.A), d.vely[k], s_vel)
                @test ap(velz(sol, refs.b, refs.A), d.velz[k], s_vel)
                @test ap(_orc_kinrv(sol, refs.b, refs.A), d.radvel[k], s_rv)

                # --- sky-plane observables, mas ------------------------------
                @test ap(raoff(sol, refs.b, refs.A), d.raoff[k], s_ang)
                @test ap(decoff(sol, refs.b, refs.A), d.decoff[k], s_ang)
                @test ap(projectedseparation(sol, refs.b, refs.A),
                         d.projectedseparation[k], s_ang)
                # Position angle wraps; compare the shorter arc.
                @test abs(rem(posangle(sol, refs.b, refs.A) - d.posangle[k],
                              2π, RoundNearest)) < RTOL

                # --- reflex of the primary about the barycentre --------------
                # Exercises the A^-1 mass split against v1's closed-form
                # companion-mass overloads.
                @test ap(raoff(sol, refs.A, bary), d.raoff_reflex[k], s_angr)
                @test ap(decoff(sol, refs.A, bary), d.decoff_reflex[k], s_angr)
                @test ap(_orc_kinrv(sol, refs.A, bary), d.radvel_reflex[k], s_rvr)
            end

            # --- proper motion -------------------------------------------
            #
            # v2's `pmra`/`pmdec` are the *exact* derivative of `raoff`/
            # `decoff`, which carries a term v1 dropped: the AU->mas scale is
            # itself a function of time through the body's line-of-sight
            # depth, so d/dt[q/(d+z)] has a second piece, -q*ż/(d+z)². The
            # residual is therefore a real correction, not a disagreement,
            # and it scales as rho/d — the angular extent over the distance.
            #
            # Two gates, so the correction is pinned from both sides:
            #
            #   * On the six exactly face-on cases (i = 0) z ≡ ż ≡ 0, the
            #     extra term vanishes identically, and v2 must reproduce v1
            #     to roundoff. Measured worst 8.0e-16; gated at 1e-12.
            #   * Everywhere else, the departure is bounded (measured worst
            #     2.5e-5 relative, at a=200 AU seen from 3.3 pc) but must be
            #     *present* — the assertion after the loop stops a deleted
            #     correction from passing as "agrees with v1".
            s_pm = max(maximum(abs, d.pmra), maximum(abs, d.pmdec))
            s_pmr = max(maximum(abs, d.pmra_reflex), maximum(abs, d.pmdec_reflex))
            pmdev = 0.0
            for (k, sol) in enumerate(traj)
                pmdev = max(pmdev,
                    abs(pmra(sol, refs.b, refs.A) - d.pmra[k]) / s_pm,
                    abs(pmdec(sol, refs.b, refs.A) - d.pmdec[k]) / s_pm,
                    abs(pmra(sol, refs.A, bary) - d.pmra_reflex[k]) / s_pmr,
                    abs(pmdec(sol, refs.A, bary) - d.pmdec_reflex[k]) / s_pmr)
            end
            if p.i == 0.0
                @test pmdev < RTOL
            else
                @test pmdev < 3e-4          # 12x the measured 2.5e-5
                pm_worst_inclined = max(pm_worst_inclined, pmdev)
            end
        end
    end

    # The depth-rate term is actually being applied somewhere in the sweep.
    # Measured 2.5e-5; asserted at 1e-6 so this fails loudly if the term is
    # ever removed and every case silently starts matching v1 exactly.
    @test pm_worst_inclined > 1e-6
end

# ---------------------------------------------------
# Oracle 2 — Kepler's equation in 256-bit BigFloat
# ---------------------------------------------------

# Elliptical branch. On M reduced to [0, 2π) the residual
# f(E) = E - e sin E - M is strictly increasing (f' = 1 - e cos E >= 1 - e > 0)
# on [0, 2π], so plain bisection cannot converge to the wrong root, whatever
# the eccentricity. 150 halvings take the bracket to ~2π·2^-150 ≈ 4e-45;
# Newton then polishes to the working precision. Nothing here shares a line
# of code, a starter, or a range reduction with the package.
function _big_kepler(M::Float64, e::Float64)
    Mb = BigFloat(M)                       # exact: M is a Float64
    eb = BigFloat(e)
    twoπ = 2 * BigFloat(π)
    Mr = Mb - twoπ * floor(Mb / twoπ)      # exact reduction at 256 bits
    lo = zero(BigFloat)
    hi = twoπ
    for _ in 1:150
        mid = (lo + hi) / 2
        if mid - eb * sin(mid) - Mr < 0
            lo = mid
        else
            hi = mid
        end
    end
    E = (lo + hi) / 2
    for _ in 1:8
        E -= (E - eb * sin(E) - Mr) / (1 - eb * cos(E))
    end
    return E
end

# Hyperbolic branch: e sinh H - H = M, also strictly increasing in H for
# e > 1. Bracket by doubling outward from ±1, then bisect and polish.
function _big_hyperbolic(M::Float64, e::Float64)
    Mb = BigFloat(M)
    eb = BigFloat(e)
    lo = -one(BigFloat)
    hi = one(BigFloat)
    while eb * sinh(lo) - lo - Mb > 0
        lo *= 2
    end
    while eb * sinh(hi) - hi - Mb < 0
        hi *= 2
    end
    for _ in 1:400
        mid = (lo + hi) / 2
        if eb * sinh(mid) - mid - Mb < 0
            lo = mid
        else
            hi = mid
        end
    end
    H = (lo + hi) / 2
    for _ in 1:10
        H -= (eb * sinh(H) - H - Mb) / (eb * cosh(H) - 1)
    end
    return H
end

@testset "oracle: Kepler's equation vs 256-bit BigFloat" begin
    setprecision(BigFloat, 256)

    # Eccentricities: exactly circular, the near-circular band where E ≈ M and
    # cancellation bites, the ordinary range, and the near-parabolic tail where
    # Markley's starter is worst-conditioned.
    ES = (0.0, 1e-12, 1e-8, 1e-4, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999,
          0.999999, 0.9999999)
    # Mean anomalies within one revolution: the endpoints and the midpoint of
    # the reduced interval are where a starter fitted on [-π, π] is worst, and
    # 2π ± 1e-9 straddles the wrap.
    MBASE = (0.0, 1e-12, 1e-6, 1e-3, 0.1, 1.0, 3.0, 3.1415, float(π), 3.15, 4.0,
             6.0, 2π - 1e-9, 2π, 2π + 1e-9, -1.0, -float(π), -3.0)
    # ...times whole revolutions, out to 1e20. `M = n(t - tp)` really does
    # reach these: an ordinary 1 AU orbit with tp 1e16 days away, or a very
    # small `a` that a hot tempering chain proposes. This is the range
    # reduction under test as much as the solver.
    MULT = (0.0, 1.0, 7.0, 1000.0, 1e6, 1e9, 1e12, 1e14, 1e15, 1e17, 1e20)

    # Compared as (sin E, cos E) rather than E, because that pair is what the
    # propagator actually consumes and it is insensitive to which revolution
    # a solver chooses to report. Both are bounded by 1, so the tolerance is
    # absolute.
    #
    # Measured worst |Δsin|,|Δcos| over the full 2376-point grid:
    #   e <= 0.99      3.7e-16      e <= 0.9999      8.4e-15
    #   e <= 0.999     2.1e-15      e <= 0.9999999   1.6e-13
    # The growth is the 1/sqrt(1-e) conditioning of the equation itself, not a
    # solver defect. Tiered so the ordinary range keeps a tight gate.
    keptol(e) = e <= 0.999 ? 5e-15 : e <= 0.99999 ? 5e-14 : 5e-13

    @testset "Markley (the default solver)" begin
        worst = 0.0
        for e in ES, m0 in MBASE, k in MULT
            M = m0 + 2π * k
            sE, cE = sincos(kepler_solver(M, e, Markley()))
            Eb = _big_kepler(M, e)
            err = max(abs(sE - Float64(sin(Eb))), abs(cE - Float64(cos(Eb))))
            err > keptol(e) && @error "Markley vs BigFloat" M e err
            worst = max(worst, err / keptol(e))
        end
        @test worst <= 1
    end

    @testset "Auto dispatches to a solver of the same quality" begin
        worst = 0.0
        for e in ES, m0 in MBASE, k in (0.0, 1.0, 1e9, 1e17)
            M = m0 + 2π * k
            sE, cE = sincos(kepler_solver(M, e, Auto()))
            Eb = _big_kepler(M, e)
            worst = max(worst, max(abs(sE - Float64(sin(Eb))),
                                   abs(cE - Float64(cos(Eb)))) / keptol(e))
        end
        @test worst <= 1
    end

    @testset "markley_sincosE (SIMD kernel) inside its documented domain" begin
        # `VREM2PI_MAX` is the package's own claim about where the branch-free
        # Cody-Waite reduction is a valid stand-in for `Base.rem2pi`; rows
        # beyond it are routed to the scalar solver by `_ma_in_kernel_range`.
        # Here that claim is checked against BigFloat rather than against the
        # scalar path, so a shared reduction bug could not hide.
        worst = 0.0
        n = 0
        for e in ES, m0 in MBASE, k in MULT
            M = m0 + 2π * k
            abs(M) <= VREM2PI_MAX || continue
            sE, cE = markley_sincosE(M, e)
            Eb = _big_kepler(M, e)
            worst = max(worst, max(abs(sE - Float64(sin(Eb))),
                                   abs(cE - Float64(cos(Eb)))) / keptol(e))
            n += 1
        end
        @test n > 1000              # the domain filter did not empty the grid
        @test worst <= 1
    end

    @testset "RootsMethod(A42)" begin
        # A bracketing method from an external package: a second opinion on
        # the same equation that shares nothing with Markley but the interface.
        worst = 0.0
        for e in ES, m0 in MBASE, k in (0.0, 1.0, 1000.0, 1e12)
            M = m0 + 2π * k
            E = kepler_solver(M, e, RootsMethod(PlanetOrbits.Roots.A42()))
            sE, cE = sincos(E)
            Eb = _big_kepler(M, e)
            worst = max(worst, max(abs(sE - Float64(sin(Eb))),
                                   abs(cE - Float64(cos(Eb)))) / keptol(e))
        end
        @test worst <= 1
    end

    @testset "HyperbolicHalley" begin
        # Measured worst relative |ΔH| 1.1e-14; gated at 1e-12 (90x margin).
        # runtests.jl already checks the *residual* of this solver; a residual
        # cannot tell a correct root from a correct root of the wrong equation,
        # which is what an independent solve of e sinh H - H = M adds.
        worst = 0.0
        for e in (1.0000001, 1.001, 1.01, 1.1, 1.5, 2.0, 5.0, 20.0, 100.0),
            M in (-1e6, -1e3, -10.0, -1.0, -1e-8, 0.0, 1e-8, 1.0, 10.0, 1e3, 1e6, 1e9)
            H = kepler_solver(M, e, HyperbolicHalley())
            Hb = Float64(_big_hyperbolic(M, e))
            worst = max(worst, abs(H - Hb) / max(abs(Hb), 1.0))
        end
        @test worst < 1e-12
    end

    @testset "Goat" begin
        # KNOWN BROKEN — inherited unchanged from v1, opt-in only, and
        # documented in its own docstring as "here for comparison purposes
        # only". It is not reached by `Auto()`, so nothing in the package uses
        # it unless a caller asks for it by name.
        #
        # Two independent failures, both reproducible without BigFloat:
        #
        #   1. It does not solve the equation to anything near machine
        #      precision. kepler_solver(1.0, 0.1, Goat()) = 1.0919981742869485
        #      against a true root of 1.0885977523978936 — a residual
        #      |E - e sin E - M| of 3.2e-3. The residual grows with e, reaching
        #      7.9e-2 at (M, e) = (3.0, 0.9999999).
        #   2. It returns NaN for small but nonzero e near a multiple of π —
        #      e.g. (M, e) = (3.1415, 1e-12) — because the `isapprox(e, 0)`
        #      short circuit has no absolute tolerance and so does not fire.
        #
        # Left as @test_broken rather than fixed or deleted: correcting a
        # numerical result is the maintainer's call, and the port would need
        # re-deriving against the paper (arXiv 2103.15829) rather than patched.
        # Flip these to @test if the solver is ever repaired.
        residual(M, e) = (E = kepler_solver(M, e, Goat()); abs(E - e * sin(E) - M))
        @test_broken residual(1.0, 0.1) < 1e-12
        @test_broken residual(3.0, 0.9) < 1e-12
        @test_broken isfinite(kepler_solver(3.1415, 1e-12, Goat()))
        # The failure is exactly as characterised above — if these stop
        # holding, Goat's behaviour changed and the notes above are stale.
        @test residual(1.0, 0.1) > 1e-4
        @test isnan(kepler_solver(3.1415, 1e-12, Goat()))
    end
end

# ---------------------------------------------------
# Oracle 3 — closed-form and published values
# ---------------------------------------------------

@testset "oracle: closed-form and published values" begin

    @testset "IAU nominal constants" begin
        # Published, exact-by-definition inputs, written out here rather than
        # read from the package:
        #   GM_sun   IAU 2015 Resolution B3, nominal solar mass parameter
        #   GM_jup   IAU 2015 Resolution B3, nominal jovian mass parameter
        #   GM_earth IAU 2015 Resolution B3, nominal terrestrial mass parameter
        #   au       IAU 2012 Resolution B2, exact
        GM_sun = 1.3271244e20        # m^3 / s^2
        GM_jup = 1.2668653e17
        GM_earth = 3.986004e14
        au = 149_597_870_700.0       # m
        day = 86_400.0               # s, exact

        # The "Kepler year": the period of a 1 AU orbit about 1 solar mass.
        # This is the constant inside Kepler's third law throughout the
        # package, and it differs from the Julian year by 1.9e-5 — an error
        # that looks right and quietly drifts, which is exactly why it is
        # worth pinning to a recomputation from first principles.
        P_kepler_days = sqrt(4π^2 / GM_sun * au^3) / day
        @test kepler_year_to_julian_day_conversion_factor ≈ P_kepler_days rtol = 1e-15

        # ...and Kepler's third law in the package must actually use it.
        earthlike = System((Body(mass=1.0, name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A);
                   a=1.0, e=0.0167, i=0.0, ω=0.0, Ω=0.0, tp=58849.0),))
        @test period(earthlike) ≈ P_kepler_days rtol = 1e-14
        # For orientation: the anomalistic year is 365.259636 d (Astronomical
        # Almanac). The 1 AU / 1 Msun period is not the same quantity — the
        # Earth-Moon system is not massless and Earth's semi-major axis is not
        # exactly 1 AU — so they agree only to 7.5e-6 relative (2.7e-3 days).
        # Gated at 1e-5: too loose to catch a precision bug, tight enough to
        # catch a period that is not a year.
        @test isapprox(period(earthlike), 365.259636; rtol=1e-5)

        # The two mass units are exact ratios of nominal mass parameters.
        @test mjup == GM_jup / GM_sun
        @test mearth == GM_earth / GM_sun
        @test msun == 1.0

        # The Julian year is 365.25 days exactly (IAU), and the parsec is
        # 648000/π AU exactly (IAU 2015 Resolution B2).
        @test year2day_julian == 365.25
        @test pc2au ≈ 648_000 / π rtol = 1e-16
    end

    @testset "the parsec, by definition" begin
        # A parsec is the distance at which 1 AU subtends 1 arcsecond. So a
        # companion 1 AU from its (massless-companion) primary, seen face-on
        # from plx = 1000 mas = 1 pc, is at 1000.000000 mas separation. This
        # gates the entire AU -> mas chain against the definition rather than
        # against another line of the same package.
        #
        # Face-on so that z ≡ 0 and the per-body depth scaling vanishes
        # identically: the identity then holds on the *default* observing
        # geometry too, not only with the opt-out.
        A = Body(mass=1.0, name=:A)
        b = Body(mass=0.0, name=:b)
        sys = System(Orbit(b, about=A; a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=58849.0);
                     plx=1000.0)
        for og in (true, false)
            sol = orbitsolve(sys, [58849.0]; observing_geometry=og)[1]
            @test projectedseparation(sol, :b, :A) ≈ 1000.0 rtol = 1e-14
            @test decoff(sol, :b, :A) ≈ 1000.0 rtol = 1e-14
            @test abs(raoff(sol, :b, :A)) < 1e-11
        end
        # Half the distance, twice the parallax, half the angle — linear, and
        # the same 1 AU physical separation.
        sys2 = System(Orbit(b, about=A; a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=58849.0);
                      plx=500.0)
        sol2 = orbitsolve(sys2, [58849.0])[1]
        @test projectedseparation(sol2, :b, :A) ≈ 500.0 rtol = 1e-14
    end

    @testset "RV semi-amplitude vs the closed form" begin
        # The textbook stellar reflex semi-amplitude
        #
        #     K1 = 2π a1 sin(i) / (P sqrt(1 - e^2)),    a1 = a m2/(m1+m2)
        #
        # and the closed-form RV curve v(ν) = -K1 [cos(ν+ω) + e cos ω], where ν
        # is the companion's true anomaly. (Sign: PlanetOrbits' line-of-sight
        # axis is positive *towards* the observer for the target-relative
        # reads used here, so the star's velocity is negative when the
        # companion is at ν + ω = 0. That convention is pinned by v1
        # agreement in oracle 1; what is tested here is the amplitude and its
        # phase dependence.)
        #
        # The two extremal epochs are located analytically from ν, through the
        # standard eccentric-anomaly relation and Kepler's equation written
        # out here — so the epoch, not just the amplitude, is refereed.
        #
        # Read with `observing_geometry=false` and through `velz` rather than
        # `radvel`: the closed form is Newtonian and knows nothing about the
        # Einstein term or the light-travel retardation, both of which v2
        # applies on the default path.
        for (m1, m2, a, e, i, ω) in ((1.0, 0.001, 5.0, 0.0, π / 2, 0.0),
                                     (1.1, 0.02, 2.0, 0.4, 1.0, 0.7),
                                     (0.5, 0.3, 0.6, 0.85, 1.3, -2.0),
                                     (2.2, 0.15, 12.0, 0.6, 0.35, 4.4))
            tp = 58849.0
            A = Body(mass=m1, name=:A)
            b = Body(mass=m2, name=:b)
            sys = System(Orbit(b, about=A; a=a, e=e, i=i, ω=ω, Ω=1.234, tp=tp))
            P = period(sys)                                # days
            a1 = a * m2 / (m1 + m2)                        # AU
            K = 2π * a1 * sin(i) / (P / year2day_julian * sqrt(1 - e^2)) *
                au2m / year2sec_julian                     # m/s

            # Epoch at which the companion's true anomaly is ν.
            function epoch_at_ν(ν)
                E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(ν / 2))
                return tp + (E - e * sin(E)) * P / 2π
            end
            t0 = epoch_at_ν(rem(-ω, 2π, RoundNearest))     # ν + ω = 0
            tπ = epoch_at_ν(rem(π - ω, 2π, RoundNearest))  # ν + ω = π
            ts = sort([t0, tπ])
            traj = orbitsolve(sys, ts; observing_geometry=false)
            bary = barycentre(sys)
            v = Dict(ts[k] => _orc_kinrv(traj[k], :A, bary) for k in 1:2)

            # Amplitude and phase dependence, both signed.
            @test v[t0] ≈ -K * (1 + e * cos(ω)) rtol = 1e-12
            @test v[tπ] ≈ K * (1 - e * cos(ω)) rtol = 1e-12
            # ...hence the half-range of the curve is K, sign conventions and
            # all cancelled out.
            @test abs(v[t0] - v[tπ]) / 2 ≈ K rtol = 1e-12
        end
    end

    @testset "Thiele-Innes constants vs the closed form" begin
        # A, B, F, G written out from the standard definitions (e.g. Wright &
        # Howard 2009 eq. 12-15; Gaia DR3 non-single-star documentation), and
        # the sky-plane offsets assembled from them independently of the
        # package's element -> state path:
        #
        #     X = cos E - e            Y = sqrt(1 - e^2) sin E
        #     Δδ (north) = A X + F Y   Δα* (east) = B X + G Y
        #
        # E comes from a plain Newton solve of Kepler's equation written here.
        for (a, e, i, ω, Ω) in ((5.0, 0.3, 0.7, 1.1, 2.2),
                                (0.8, 0.0, 0.0, 0.0, 0.0),
                                (20.0, 0.85, 2.4, 5.0, 5.5),
                                (2.0, 0.5, π / 2, 3.0, 0.1))
            tp = 59000.0
            A_ = a * (cos(ω) * cos(Ω) - sin(ω) * sin(Ω) * cos(i))
            B_ = a * (cos(ω) * sin(Ω) + sin(ω) * cos(Ω) * cos(i))
            F_ = -a * (sin(ω) * cos(Ω) + cos(ω) * sin(Ω) * cos(i))
            G_ = -a * (sin(ω) * sin(Ω) - cos(ω) * cos(Ω) * cos(i))

            star = Body(mass=1.1, name=:A)
            comp = Body(mass=0.0, name=:b)
            sys = System(Orbit(comp, about=star; a=a, e=e, i=i, ω=ω, Ω=Ω, tp=tp))
            P = period(sys)

            # The package's own constants must be those numbers.
            ti = thieleinnes(sys)
            @test ti.A ≈ A_ rtol = 1e-14 atol = 1e-14 * a
            @test ti.B ≈ B_ rtol = 1e-14 atol = 1e-14 * a
            @test ti.F ≈ F_ rtol = 1e-14 atol = 1e-14 * a
            @test ti.G ≈ G_ rtol = 1e-14 atol = 1e-14 * a

            # ...and the trajectory must be the one they describe. Measured
            # worst departure 5e-15 AU on a ~5 AU orbit; gated at 1e-12
            # relative to a.
            for frac in (0.0, 0.13, 0.4, 0.77, 0.96)
                MA = 2π * frac
                E = MA
                for _ in 1:200
                    E -= (E - e * sin(E) - MA) / (1 - e * cos(E))
                end
                X = cos(E) - e
                Y = sqrt(1 - e^2) * sin(E)
                sol = orbitsolve(sys, [tp + frac * P]; observing_geometry=false)[1]
                @test isapprox(posx(sol, :b, :A), B_ * X + G_ * Y;
                               rtol=1e-12, atol=1e-12 * a)
                @test isapprox(posy(sol, :b, :A), A_ * X + F_ * Y;
                               rtol=1e-12, atol=1e-12 * a)
            end
        end
    end

    @testset "vis-viva at 1 AU vs the published mean orbital speed" begin
        # A coarse but genuinely external units check: a circular 1 AU orbit
        # about 1 solar mass has speed 2π AU / (Kepler year) = 29.7847 km/s,
        # against Earth's published mean orbital speed of 29.78 km/s (NASA
        # Earth fact sheet; the Earth's own orbit is slightly eccentric and
        # the Sun is not exactly 1 nominal solar mass, hence the loose gate).
        # It cannot catch a 1e-12 error; it catches a factor of 365.25.
        A = Body(mass=1.0, name=:A)
        b = Body(mass=0.0, name=:b)
        sys = System(Orbit(b, about=A; a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=58849.0))
        sol = orbitsolve(sys, [58900.0]; observing_geometry=false)[1]
        speed_ms = hypot(velx(sol, :b, :A), vely(sol, :b, :A), velz(sol, :b, :A)) *
                   au2m / year2sec_julian
        @test isapprox(speed_ms / 1000, 29.78; rtol=1e-3)
    end
end
