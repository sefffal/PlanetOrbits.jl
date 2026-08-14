# ---------------------------------------------------
# Degenerate parameter values: classification on the primal, and the two
# float-range bands where a well-posed orbit used to come back `NaN`.
#
# The shared theme is that a *branch* — elliptical or hyperbolic, batch kernel
# or scalar solver, guard fired or not — must be a function of the parameter's
# value alone. ForwardDiff orders `Dual`s lexicographically, so `<`, `>=`,
# `==`, `iszero` and `sign` all consult the partials when the values tie; a
# predicate written on the `Dual` therefore answers a different question on the
# gradient path than on the value path, and the two disagree precisely at the
# degeneracies the predicates exist to catch.
#
# These are `NaN`-hunting tests, so they mostly assert the *absence* of `NaN`
# and the *presence* of an `OrbitDomainError` — which is the contract a sampler
# needs (see the `OrbitDomainError` docstring): a bad proposal scores `-Inf`,
# it does not poison the log-density.
# ---------------------------------------------------

const D1 = ForwardDiff.Dual{:degen}
_dual(x, ∂=0.0) = D1(float(x), float(∂))

# Raw barycentric states — no observing pass. The observing pass legitimately
# diverges for these orbits (see "superluminal" below), and what is under test
# here is the Kepler solve.
function _raw_rowstates(sys, epochs; method=KeplerianApprox())
    traj = Trajectory(sys, epochs)
    PlanetOrbits.frame_pass!(traj, sys.frame)
    PlanetOrbits.propagate!(traj, sys, method)
    return (copy(traj.rx), copy(traj.ry), copy(traj.rz),
            copy(traj.rvx), copy(traj.rvy), copy(traj.rvz))
end

# Inside its documented domain the kernel reduction must be indistinguishable
# from `Base.rem2pi` — same residue class, and inside [-π, π] where Markley's
# starter is fitted.
function vrem2pi_ok(x)
    v = PlanetOrbits.vrem2pi(x)
    return abs(v) <= π + 1e-12 && abs(v - rem2pi(x, RoundNearest)) <= 1e-15
end

_tiny_system(a; e=0.1, tp=0.0) = System(
    (Body(mass=1.0, name=:A), Body(mass=1e-9, name=:b)),
    (Orbit(Body(mass=1e-9, name=:b), about=Body(mass=1.0, name=:A);
           a=a, e=e, i=0.5, ω=0.1, Ω=0.2, tp=tp),))

@testset "classification predicates use the primal" begin
    # --- e == 1 exactly: the parabolic guard must fire for both signs -------
    #
    # This is the reported failure. `hyperbolic = e >= 1` sent
    # `Dual(1.0, +∂)` down the hyperbolic branch (lexicographically greater),
    # where `e > 1` was *also* true, so the parabolic guard never fired and the
    # row was built with sqrt1me2 = −√(e²−1) = −0.0 carrying an infinite
    # derivative. `Dual(1.0, −∂)` took the elliptical branch instead and
    # divided by √(1−e²) = 0. Same orbit, two different failures, neither of
    # them the error the value path raises.
    for ∂ in (1.0, -1.0, 0.0)
        e = _dual(1.0, ∂)
        @test_throws PlanetOrbits.OrbitDomainError PlanetOrbits.Row(
            _dual(1.0), e, _dual(0.5), _dual(0.1), _dual(0.2), _dual(0.0), _dual(1.0))
    end
    # …and through the public constructor, where the elements promote.
    @test_throws PlanetOrbits.OrbitDomainError System(
        (Body(mass=_dual(1.0), name=:A), Body(mass=_dual(1e-3), name=:b)),
        (Orbit(Body(mass=_dual(1e-3), name=:b), about=Body(mass=_dual(1.0), name=:A);
               a=_dual(-3.0), e=_dual(1.0, 1.0), i=_dual(0.5),
               ω=_dual(0.1), Ω=_dual(0.2), tp=_dual(0.0)),))

    # --- The conic a row lands on does not depend on the perturbation ------
    for ev in (0.8, 0.999999, 1.000001, 1.4), ∂ in (1.0, -1.0)
        a = ev > 1 ? -3.0 : 3.0
        r = PlanetOrbits.Row(_dual(a), _dual(ev, ∂), _dual(0.5), _dual(0.1),
                             _dual(0.2), _dual(0.0), _dual(1.0))
        rf = PlanetOrbits.Row(a, ev, 0.5, 0.1, 0.2, 0.0, 1.0)
        @test r.hyperbolic == rf.hyperbolic
        # …and the derived constants agree with the Float64 row's in value.
        @test ForwardDiff.value(r.n) == rf.n
        @test ForwardDiff.value(r.sqrt1me2) == rf.sqrt1me2
        @test ForwardDiff.value(r.J) == rf.J
    end

    # --- `a` normalization for hyperbolae is a classification too ----------
    # `hyperbolic && a > 0 && (a = -a)`: at a == 0 exactly, `Dual(0.0, +∂) > 0`
    # is true, so the sign flip fired and `-0.0` then read as "a < 0, fine".
    @test_throws PlanetOrbits.OrbitDomainError PlanetOrbits.Row(
        _dual(0.0, 1.0), _dual(1.4), _dual(0.5), _dual(0.1),
        _dual(0.2), _dual(0.0), _dual(1.0))

    # --- `_meanmotion` must classify the same way `Row` does ---------------
    # It runs one layer earlier, converting an `M0`/`θ` phase to `tp`, so a
    # disagreement would place `tp` using the wrong conic's mean motion.
    for ev in (0.8, 1.4), ∂ in (1.0, -1.0)
        a = ev > 1 ? -3.0 : 3.0
        @test PlanetOrbits._meanmotion(_dual(a), _dual(ev, ∂), _dual(1.0)) ==
              PlanetOrbits._meanmotion(_dual(a), _dual(ev, 0.0), _dual(1.0))
    end
    @test_throws PlanetOrbits.OrbitDomainError PlanetOrbits._MA_from_θ(
        _dual(0.3), _dual(1.0, 1.0), _dual(0.5), _dual(0.1), _dual(0.2))

    # --- `Auto` picks the solver on the value ------------------------------
    # `Dual(0.999…, +∂)` must still reach the elliptical solver, and the
    # implicit-rule methods must agree with `Auto` about which branch it is.
    for ev in (0.5, 0.99), ∂ in (1.0, -1.0)
        MA = _dual(0.7, 0.3)
        @test PlanetOrbits.kepler_solver(MA, _dual(ev, ∂), PlanetOrbits.Auto()) ===
              PlanetOrbits.kepler_solver(MA, _dual(ev, ∂), PlanetOrbits.Markley())
    end

    # --- `iszero` is lexicographic: a differentiated zero is not zero ------
    # Zero parallax reached `1000/plx` because `iszero(Dual(0.0, ∂))` is false.
    @test_throws "parallax is zero" PlanetOrbits.AbsoluteFrame(
        ra=_dual(45.0), dec=_dual(-30.0), plx=_dual(0.0, 1.0), pmra=_dual(100.0),
        pmdec=_dual(-50.0), rv=_dual(25e3), ref_epoch=_dual(57388.5))
end

@testset "non-finite elements are a domain error, not a NaN" begin
    # `rem2pi(Inf, RoundDown)` is `NaN` and `sincos(NaN)` is `(NaN, NaN)`, so
    # before this guard a non-finite angle built a row without complaint and
    # returned `NaN` for every observable. Unbounded priors reach ±Inf under
    # the unconstrained-space transform.
    base = (1.0, 0.1, 0.5, 0.1, 0.2, 0.0, 1.0)   # a, e, i, ω, Ω, tp, M
    for slot in 2:6, bad in (Inf, -Inf, NaN)
        args = collect(Float64, base)
        args[slot] = bad
        @test_throws PlanetOrbits.OrbitDomainError PlanetOrbits.Row(args...)
        # …and on the gradient path, where the guard's own comparison is the
        # thing at risk.
        @test_throws PlanetOrbits.OrbitDomainError PlanetOrbits.Row(
            map(x -> _dual(x, 1.0), args)...)
    end
    # The finite elements it is meant to let through still pass.
    @test PlanetOrbits.Row(base...) isa PlanetOrbits.Row
end

@testset "derived constants that leave float range are a domain error" begin
    # `a³` underflows to zero below a ≈ 1e-107, so `n = √(GM/a³)` is `Inf`,
    # the mean anomaly is `Inf`, and every reduction of it — Payne-Hanek
    # included — is `NaN`. The ellipse is well defined; the constants are not.
    for a in (1e-108, 1e-120, 1e-200, 5e-324)
        @test_throws PlanetOrbits.OrbitDomainError _tiny_system(a)
        @test_throws PlanetOrbits.OrbitDomainError _tiny_system(-a; e=1.4)
    end
    # The guard is on the constants, not a blanket ban on small `a`: one
    # decade above the underflow the row builds and solves finitely.
    for a in (1e-100, 1e-60, 1e-20)
        sys = _tiny_system(a)
        @test isfinite(sys.rows[1].n)
        @test all(isfinite, _raw_rowstates(sys, [58000.0, 59000.0])[1])
    end
    # It costs the construction path no allocation (branch on a value).
    _dc_gate(a) = _tiny_system(a)
    _dc_gate(8.0); _dc_gate(_dual(8.0, 1.0))
    if DYNAMIC_ALLOC_GATE
        @test (@allocated _dc_gate(8.0)) == 0
        @test (@allocated _dc_gate(_dual(8.0, 1.0))) == 0
    else
        @test_skip (@allocated _dc_gate(8.0)) == 0
    end
    errs = filter(!_ac_benign, AllocCheck.check_allocs(_dc_gate, (Float64,)))
    isempty(errs) || display(errs[1])
    @test isempty(errs)
end

@testset "angle reduction: the batch kernel's domain" begin
    # `vrem2pi` is one Cody-Waite step and is exact where it is used at all.
    for x in (0.0, 1.0, -1.0, 1e3, 1e8, 1e12, PlanetOrbits.VREM2PI_MAX,
              -PlanetOrbits.VREM2PI_MAX, nextfloat(-1.0), 2.0^47 + 0.5)
        @test vrem2pi_ok(x)
    end
    # Beyond the domain it is *not* a stand-in for `Base.rem2pi` — which is
    # the whole reason for the routing below, so pin the fact rather than
    # leaving it to be rediscovered.
    @test abs(PlanetOrbits.vrem2pi(1e20)) > π
    @test isnan(first(PlanetOrbits.markley_sincosE(1e63, 0.3)))
end

@testset "huge mean anomalies: SIMD and scalar agree, and neither is NaN" begin
    # Two independent ways to a mean anomaly past 2^48, both reachable from
    # ordinary priors: a short period, and a `tp` far from the epochs.
    #
    # Before the routing fix the batch kernel silently disagreed with the
    # scalar solver from |M| ≈ 2^52 and returned `NaN` past ~2^156, while
    # `Base.rem2pi` reduced the same argument exactly. A likelihood that
    # happened to be evaluated with `simd=true` therefore saw a different
    # orbit — or none at all — from one evaluated with `simd=false`.
    epochs = collect(range(58000.0, 60000.0, length=37))
    simd = KeplerianApprox(simd=true)
    scal = KeplerianApprox(simd=false)

    # Agreement is asserted two ways, and the difference between them is the
    # point. Inside the kernel's domain the two paths are genuinely different
    # code and agree to the documented ~4e-15; outside it, routing sends the
    # batch path to the scalar solver, so they agree *bit for bit*. A
    # regression that widened the domain would break the `==`, and one that
    # broke the reduction would break the tolerance.
    function agree(sys, eps)
        vs = _raw_rowstates(sys, eps; method=simd)
        vr = _raw_rowstates(sys, eps; method=scal)
        nan = any(c -> any(isnan, c), vr) || any(c -> any(isnan, c), vs)
        scale = maximum(c -> maximum(abs, c), vr)
        rel = maximum(map((s, r) -> maximum(abs, s .- r), vs, vr)) / scale
        traj = Trajectory(sys, eps)
        PlanetOrbits.frame_pass!(traj, sys.frame)
        fellback = !PlanetOrbits._use_simd(simd, traj, sys.rows[1])
        return (; nan, rel, exact=(vs == vr), fellback)
    end

    @testset "short period (a = 1e$(round(Int, log10(a))))" for a in
            (1e-3, 1e-8, 1e-11, 1e-15, 1e-20, 1e-40, 1e-60, 1e-100)
        r = agree(_tiny_system(a), epochs)
        @test !r.nan
        @test r.rel <= 4e-15
        @test r.exact == r.fellback
    end

    @testset "distant tp (tp = -1e$(round(Int, log10(-tp))))" for tp in
            (-1e6, -1e12, -1e16, -1e20, -1e40, -1e100, -1e300)
        r = agree(_tiny_system(1.0; tp=tp), epochs)   # an ordinary 1 AU orbit
        @test !r.nan
        @test r.rel <= 4e-15
        @test r.exact == r.fellback
    end

    # The routing predicate is what makes that true, so assert it directly:
    # in range → batch kernel, out of range → scalar path. (`frame_pass!`
    # first: `t_em` is what the mean anomaly is built from, and it is `undef`
    # until the frame pass writes it.)
    routed(sys, eps) = begin
        traj = Trajectory{eltype(sys.masses)}(sys, eps)
        PlanetOrbits.frame_pass!(traj, sys.frame)
        PlanetOrbits._use_simd(simd, traj, sys.rows[1])
    end
    @test routed(_tiny_system(1.0), epochs)
    @test !routed(_tiny_system(1.0; tp=-1e20), epochs)
    # A `NaN` epoch must route to the scalar path too, rather than reach a
    # kernel whose contract excludes it.
    @test !routed(_tiny_system(1.0), [58000.0, NaN])

    # The Dual batch path shares the kernel, so it shares the domain: a row
    # out of range must still differentiate, on the scalar fallback.
    #
    # No finite-difference cross-check here, deliberately. At tp = -1e20 the
    # mean anomaly is ~1.7e18 rad, so ∂M/∂a ~ 2.6e18 and *any* usable FD step
    # sweeps billions of complete orbits — FD returns an aliased number
    # (~1e6 against ForwardDiff's ~7e19) and would be measuring its own step
    # size, not the derivative. The meaningful check is that the two solver
    # routes give the same derivative, which is what the routing must preserve.
    f(method) = θ -> begin
        s = _tiny_system(θ[1]; tp=-1e20)
        traj = Trajectory{eltype(θ)}(s, epochs)
        PlanetOrbits.frame_pass!(traj, s.frame)
        PlanetOrbits.propagate!(traj, s, method)
        sum(traj.rx)
    end
    g = ForwardDiff.gradient(f(simd), [1.0])
    @test all(isfinite, g)
    @test g == ForwardDiff.gradient(f(scal), [1.0])
    let sysd = _tiny_system(_dual(1.0, 1.0); tp=_dual(-1e20))
        trajd = Trajectory{eltype(sysd.masses)}(sysd, epochs)
        PlanetOrbits.frame_pass!(trajd, sysd.frame)
        @test !PlanetOrbits._use_dual_simd(simd, trajd, sysd.rows[1])
    end
end

@testset "AHL21 step count is a function of the values" begin
    # `while τ - istep*h >= h` and `δ > 0` decide how many symplectic steps to
    # take. Under `Dual`s those are lexicographic comparisons, so an epoch
    # sitting exactly on a step boundary — the common case, since `t0` defaults
    # to `ref_epoch` and epoch grids are often regular — took a different
    # *number of steps* depending on the sign of a partial. The value carried
    # through a gradient evaluation would then no longer be bit-identical to
    # the plain Float64 one.
    h = 5.0
    t0 = 58000.0
    epochs = [t0, t0 + h, t0 + 2h, t0 + 3h, t0 - h, t0 - 2h]   # all on boundaries
    mk(m) = System((Body(mass=m, name=:A), Body(mass=oftype(m, 1e-3), name=:b)),
                   (Orbit(Body(mass=oftype(m, 1e-3), name=:b),
                          about=Body(mass=m, name=:A);
                          a=oftype(m, 1.0), e=oftype(m, 0.1), i=oftype(m, 0.5),
                          ω=oftype(m, 0.1), Ω=oftype(m, 0.2), tp=oftype(m, 57900.0)),);
                   plx=oftype(m, 24.5), ra=oftype(m, 45.0), dec=oftype(m, -30.0),
                   pmra=oftype(m, 100.0), pmdec=oftype(m, -50.0),
                   rv=oftype(m, 25e3), ref_epoch=oftype(m, t0))
    method = AHL21(h=h, t0=t0)
    solve(sys) = begin
        traj = Trajectory{eltype(sys.masses)}(sys, epochs)
        PlanetOrbits.frame_pass!(traj, sys.frame)
        PlanetOrbits.propagate!(traj, sys, method)
        copy(traj.x)
    end
    ref = solve(mk(1.0))
    for ∂ in (1.0, -1.0)
        got = solve(mk(_dual(1.0, ∂)))
        @test ForwardDiff.value.(got) == ref     # bit-identical, not ≈
        @test all(isfinite, ForwardDiff.partials.(got, 1))
    end
end
