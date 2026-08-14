# ---------------------------------------------------
# `reframe` and the Cartesian-state orbit constructor
#
# These two exist to serve a two-stage frame: build and solve at `Parallax`
# level, derive the absolute frame from that solution, install it. The tests
# below gate the two properties that makes rely on — that reframing is exact
# rather than approximate, and that a relative state round-trips through the
# elements.
# ---------------------------------------------------

@testset "reframe" begin

    # A system with enough structure that a silent reuse-vs-rebuild difference
    # would show: unequal masses (so A⁻¹ is not symmetric), a hierarchical
    # third body, and fluxes (which are promoted alongside the frame).
    function refsystem(; kwargs...)
        A = Body(mass=1.29, name=:A, flux=(G=1.0,))
        b = Body(mass=0.0069, name=:b, flux=(G=1e-4,))
        c = Body(mass=0.0019, name=:c)
        return System((A, b, c),
            (Orbit(b, about=A; a=0.059, e=0.012, i=0.8, ω=1.1, Ω=2.2, tp=58849.0),
             Orbit(c, about=(A, b); a=0.83, e=0.26, i=1.2, ω=0.4, Ω=5.1, tp=58200.0));
            kwargs...)
    end

    epochs = [57000.0, 58849.0, 59123.5, 60000.0]

    noframe   = (;)
    parallax  = (; plx=74.19)
    absolute  = (; plx=74.19, ra=24.199, dec=41.405, pmra=-172.57, pmdec=-381.03,
                   rv=-28.0e3, ref_epoch=57388.0)

    # Query the whole observable surface a frame can touch, so an unconverted
    # field cannot hide in a component nothing reads.
    function probe(sys, kind)
        traj = orbitsolve(sys, epochs)
        (; A, b, c) = bodies(sys)
        out = Float64[]
        for i in eachindex(epochs)
            sol = traj[i]
            append!(out, (posx(sol, b, A), posy(sol, b, A), posz(sol, b, A),
                          velx(sol, c, A), vely(sol, c, A), velz(sol, c, A),
                          radvel(sol, A, barycentre(sys))))
            kind === :noframe && continue
            append!(out, (raoff(sol, b, A), decoff(sol, c, barycentre(sys))))
            kind === :parallax && continue
            append!(out, (frame_ra(sol), frame_dec(sol),
                          frame_pmra(sol), frame_pmdec(sol), frame_rv(sol),
                          raoff(sol, A, framedirection), decoff(sol, A, framedirection)))
        end
        return out
    end

    @testset "matches fresh construction, $(name)" for (name, kw, kind) in
            (("NoFrame", noframe, :noframe),
             ("Parallax", parallax, :parallax),
             ("AbsoluteFrame", absolute, :absolute))
        fresh = refsystem(; kw...)
        # Every starting frame, reframed *to* this one, must land on `fresh`.
        for (from, fromkw) in (("NoFrame", noframe), ("Parallax", parallax),
                               ("AbsoluteFrame", absolute))
            got = reframe(refsystem(; fromkw...); kw...)
            @test typeof(got) === typeof(fresh)
            # Bit-identical, not `approx`: nothing is recomputed, so anything
            # short of equality means a field was rebuilt rather than reused.
            @test probe(got, kind) == probe(fresh, kind)
            @test got.masses == fresh.masses
            @test got.Ainv == fresh.Ainv
            @test got.rows == fresh.rows
            @test got.fluxes == fresh.fluxes
        end
    end

    @testset "frame object form" begin
        sys = refsystem(; parallax...)
        fr = PlanetOrbits._make_frame(; absolute...)
        @test probe(reframe(sys, fr), :absolute) == probe(refsystem(; absolute...), :absolute)
        @test reframe(sys, PlanetOrbits.NoFrame()).frame === PlanetOrbits.NoFrame()
    end

    @testset "round trip through every level" begin
        # Reframing is a replacement, not an accumulation: going out to an
        # absolute frame and back must return the original exactly.
        for kw in (noframe, parallax, absolute)
            sys = refsystem(; kw...)
            there_and_back = reframe(reframe(sys; absolute...); kw...)
            @test there_and_back.frame == sys.frame
            @test there_and_back.rows == sys.rows
        end
    end

    @testset "scalar type promotion" begin
        # A wider frame widens the whole system: this is the anchored-frame
        # case, where a Float64-mass system takes a Dual-valued frame.
        sys = refsystem(; parallax...)
        @test eltype(sys.masses) === Float64
        d = ForwardDiff.Dual{:tag}(74.19, 1.0)
        wide = reframe(sys; plx=d)
        @test eltype(wide.masses) <: ForwardDiff.Dual
        @test eltype(wide.Ainv) <: ForwardDiff.Dual
        @test eltype(values(wide.fluxes)[1]) <: ForwardDiff.Dual
        @test ForwardDiff.value(wide.frame.plx) == 74.19

        # ... and a narrower frame does not narrow it (documented behaviour:
        # widening cannot lose a derivative, narrowing would).
        back = reframe(wide; parallax...)
        @test eltype(back.masses) <: ForwardDiff.Dual
    end

    @testset "AD through reframe matches finite differences" begin
        # Differentiate an observable w.r.t. the *frame* quantities with the
        # frame installed by `reframe`, which is the path the anchored
        # parameterization takes.
        base = [74.19, -172.57, -381.03, -28.0e3]
        function f(p)
            sys = reframe(refsystem(; plx=p[1]);
                plx=p[1], ra=24.199, dec=41.405, pmra=p[2], pmdec=p[3],
                rv=p[4], ref_epoch=57388.0)
            traj = orbitsolve(sys, epochs)
            (; A, b) = bodies(sys)
            return sum(i -> raoff(traj[i], A, framedirection) +
                            decoff(traj[i], A, framedirection) +
                            raoff(traj[i], b, A), eachindex(epochs))
        end
        ad = ForwardDiff.gradient(f, base)
        fd = FiniteDiff.finite_difference_gradient(f, base)
        @test all(isfinite, ad)
        @test ad ≈ fd rtol = 1e-5

        # And w.r.t. masses, which must still flow through the carried-over
        # A⁻¹ rather than being frozen by the reframe.
        function g(m)
            A = Body(mass=m[1], name=:A)
            b = Body(mass=m[2], name=:b)
            sys = System((A, b), (Orbit(b, about=A; a=0.059, e=0.012, i=0.8,
                                        ω=1.1, Ω=2.2, tp=58849.0),); plx=74.19)
            traj = orbitsolve(reframe(sys; absolute...), epochs)
            return sum(i -> raoff(traj[i], A, framedirection), eachindex(epochs))
        end
        mbase = [1.29, 0.0069]
        @test ForwardDiff.gradient(g, mbase) ≈
              FiniteDiff.finite_difference_gradient(g, mbase) rtol = 1e-5
    end
end

@testset "Cartesian-state orbit constructor" begin

    # `Orbit(…; x, y, z, vx, vy, vz, epoch)` is the wide-companion
    # parameterization: sample a relative state, derive the elements.

    A = Body(mass=1.29, name=:A)
    B = Body(mass=0.2256, name=:B)
    Mrow = 1.29 + 0.2256

    # State of `ext` relative to `int` at `epoch`. Deliberately *raw* — the
    # propagator's states at the common emission epoch, with no `observe_pass!`.
    #
    # This is the distinction the constructor lives on: it takes a dynamical
    # state, and `posx`/`velx` off a full `orbitsolve` are the *apparent* one,
    # each body retarded to its own emission time. Measured on a 4 AU orbit the
    # two differ by 7.7e-6 of a period in `tp` — v/c, exactly as designed, and
    # far above any tolerance worth setting. Round-tripping through
    # `orbitsolve` would therefore be testing the size of the retardation, not
    # the constructor.
    function stateof(sys, epoch, ext, int)
        traj = Trajectory(sys, [float(epoch)])
        PlanetOrbits.frame_pass!(traj, sys.frame)
        PlanetOrbits.propagate!(traj, sys, KeplerianApprox())
        sol = traj[1]
        return (posx(sol, ext, int), posy(sol, ext, int), posz(sol, ext, int),
                velx(sol, ext, int), vely(sol, ext, int), velz(sol, ext, int))
    end

    @testset "elements → state → elements round trip" begin
        # A deterministic sweep rather than rand(): a failure has to be
        # reproducible from the printed case, and this covers the grid the
        # acceptance criterion asks for (e ∈ [0, 0.95], all angles, both z
        # signs — the latter via i above and below π/2).
        nfail = 0
        for a in (0.3, 4.0, 750.0), e in (0.0, 0.21, 0.63, 0.95),
            i in (0.2, 1.3, 1.9, 2.9), ω in (0.0, 2.1, 4.7),
            Ω in (0.3, 3.3, 5.9), M0 in (0.1, 2.5, 4.9)

            epoch = 57388.0
            o = Orbit(B, about=A; a, e, i, ω, Ω, M0, epoch)
            sys = System((A, B), (o,))
            st = stateof(sys, epoch, :B, :A)
            o2 = Orbit(B, about=A; x=st[1], y=st[2], z=st[3],
                       vx=st[4], vy=st[5], vz=st[6], epoch)

            # Angles compare modulo 2π: `Row` normalizes Ω into [0, 2π) while
            # the extraction's `atan` returns (-π, π], so the two spellings of
            # the same orientation differ by exactly 2π.
            angleq(x, y) = isapprox(rem2pi(x - y, RoundNearest), 0.0; atol=1e-6)
            # `tp` is periapsis *passage*, so it is only defined modulo the
            # period, and the extraction lands on a different branch than the
            # `M0` parametrization did: `atan` returns M ∈ (−π, π], so a model
            # given M0 = 4.9 comes back with M0 − 2π and a `tp` one whole
            # period later. Compare the phase itself, which is what either
            # spelling encodes.
            P = PlanetOrbits._period(System((A, B), (o,)).rows[1])
            MAof(orb) = 2π * (epoch - orb.tp) / P

            ok = isapprox(o2.a, o.a; rtol=1e-9) &&
                 isapprox(o2.e, o.e; rtol=1e-8, atol=1e-9) &&
                 angleq(o2.i, o.i) && angleq(o2.Ω, o.Ω)
            # ω and the phase are separately meaningful only away from the
            # circular limit: at e → 0 there is no periapsis to measure ω
            # from, and the extraction reports ω = 0 having folded the whole
            # angle into the phase. What survives the limit is their sum —
            # the argument of latitude, ν = M at e = 0 — so that is what a
            # circular case is held to.
            ok &= e > 1e-6 ? (angleq(o2.ω, o.ω) && angleq(MAof(o2), MAof(o))) :
                             angleq(o2.ω + MAof(o2), o.ω + MAof(o))
            ok || (nfail += 1; @info "round-trip mismatch" a e i ω Ω M0 got=(o2.a, o2.e, o2.i, o2.ω, o2.Ω, MAof(o2), MAof(o)))

            # The state is the invariant, and the one that subsumes every
            # branch choice above: rebuild from the extracted elements and
            # re-solve, at the construction epoch and away from it.
            sys2 = System((A, B), (o2,))
            for dt in (0.0, 0.037P, 0.61P)
                @test all(isapprox.(stateof(sys2, epoch + dt, :B, :A),
                                    stateof(sys, epoch + dt, :B, :A);
                                    rtol=1e-7, atol=1e-11))
            end
        end
        @test nfail == 0
    end

    @testset "solve at t_ref reproduces the state" begin
        # Off the construction epoch too: the elements must describe the same
        # orbit, not merely match at one instant.
        st = (12.0, -4.5, 3.1, -0.8, 1.4, 0.25)
        epoch = 59000.0
        o = Orbit(B, about=A; x=st[1], y=st[2], z=st[3],
                  vx=st[4], vy=st[5], vz=st[6], epoch)
        sys = System((A, B), (o,))
        @test all(isapprox.(stateof(sys, epoch, :B, :A), st; rtol=1e-10, atol=1e-12))

        # Vis-viva is the energy invariant the constructor inverts; check it
        # holds at an unrelated epoch, which pins `a` independently.
        st_later = stateof(sys, epoch + 900.0, :B, :A)
        μ = PlanetOrbits.GM_sun_au3_julianyr2 * Mrow
        r = hypot(st_later[1], st_later[2], st_later[3])
        v2 = st_later[4]^2 + st_later[5]^2 + st_later[6]^2
        @test 2 / r - v2 / μ ≈ 1 / o.a rtol = 1e-9
    end

    @testset "hyperbolic states are supported, not rejected" begin
        # Unbound initial conditions are a *designed* use (the docstring calls
        # them the natural way to specify a hyperbolic orbit), so this gates
        # against a well-meaning "guard" being added later.
        μ = PlanetOrbits.GM_sun_au3_julianyr2 * Mrow
        r = 30.0
        vesc = sqrt(2μ / r)
        o = Orbit(B, about=A; x=r, y=0.0, z=0.0,
                  vx=0.2 * vesc, vy=1.4 * vesc, vz=0.0, epoch=59000.0)
        @test o.e > 1
        @test o.a < 0
        sys = System((A, B), (o,))
        @test all(isfinite, stateof(sys, 59000.0, :B, :A))
    end

    @testset "near-parabolic states fail loudly" begin
        # `a = inv(2/r − v²/μ)` blows up as the specific energy passes through
        # zero, so a state within rounding of parabolic yields a huge-but-finite
        # `a` and a downstream answer that is meaningless rather than wrong-by-a-
        # little. Error, naming the energy.
        μ = PlanetOrbits.GM_sun_au3_julianyr2 * Mrow
        r = 30.0
        vesc = sqrt(2μ / r)
        for scale in (1.0, 1 + 1e-14, 1 - 1e-14)
            @test_throws ErrorException Orbit(B, about=A;
                x=r, y=0.0, z=0.0, vx=0.0, vy=scale * vesc, vz=0.0, epoch=59000.0)
        end
        # ... but a state only mildly eccentric on either side is fine.
        for scale in (0.9, 1.1)
            o = Orbit(B, about=A; x=r, y=0.0, z=0.0,
                      vx=0.0, vy=scale * vesc, vz=0.0, epoch=59000.0)
            @test isfinite(o.a) && isfinite(o.e)
        end
    end

    @testset "AD through the constructor matches finite differences" begin
        base = [12.0, -4.5, 3.1, -0.8, 1.4, 0.25]
        epochs2 = [58900.0, 59000.0, 59400.0]
        function f(s)
            o = Orbit(B, about=A; x=s[1], y=s[2], z=s[3],
                      vx=s[4], vy=s[5], vz=s[6], epoch=59000.0)
            sys = System((A, B), (o,); plx=74.19)
            traj = orbitsolve(sys, epochs2)
            return sum(i -> raoff(traj[i], :B, :A) + decoff(traj[i], :B, :A) +
                            radvel(traj[i], :A, barycentre(sys)), eachindex(epochs2))
        end
        ad = ForwardDiff.gradient(f, base)
        fd = FiniteDiff.finite_difference_gradient(f, base)
        @test all(isfinite, ad)
        @test ad ≈ fd rtol = 1e-5

        # Through the masses too — they set μ, so they enter the element
        # extraction itself, not just the propagation.
        function g(m)
            Aa = Body(mass=m[1], name=:A)
            Bb = Body(mass=m[2], name=:B)
            o = Orbit(Bb, about=Aa; x=12.0, y=-4.5, z=3.1,
                      vx=-0.8, vy=1.4, vz=0.25, epoch=59000.0)
            traj = orbitsolve(System((Aa, Bb), (o,); plx=74.19), epochs2)
            return sum(i -> raoff(traj[i], :B, :A), eachindex(epochs2))
        end
        mbase = [1.29, 0.2256]
        @test ForwardDiff.gradient(g, mbase) ≈
              FiniteDiff.finite_difference_gradient(g, mbase) rtol = 1e-5
    end
end
