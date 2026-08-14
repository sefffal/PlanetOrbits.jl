using PlanetOrbits
using PlanetOrbits: orbit, rvorbit   # v1-compat sugar; deliberately unexported
using PlanetOrbits: Body, Orbit, System   # not exported: Octofitter owns these names (§5)
using Test
import ForwardDiff
import FiniteDiff
using StaticArrays
using InteractiveUtils: code_llvm
import AllocCheck

include("fixtures/v1_reference.jl")

approx(a, b) = isapprox(a, b; rtol=1e-11, atol=1e-10)

# v1's `radvel` was the *kinematic* line-of-sight velocity. v2's is the
# spectroscopic one: it adds the Einstein (second-order Doppler plus
# gravitational-redshift) difference between the two references, which is a
# term v1 never had rather than a change of precision. The v1 fixtures
# therefore gate this — the kinematic half, in v1's units — and the Einstein
# term has its own testset ("radvel Einstein term").
kinrv(sol, refs...) = velz(sol, refs...) * PlanetOrbits.au2m / PlanetOrbits.year2sec_julian

# Solve *without* `observe_pass!`, i.e. raw barycentric states in the
# reference triad at the common emission epoch. Physical invariants of the
# Keplerian/N-body solution (vis-viva, energy, angular momentum, A⁻¹
# properties) are properties of those states; the observing pass deliberately
# perturbs them at v/c order by retarding each body to its own emission time.
function rawsolve(sys, epochs)
    ts = epochs isa AbstractVector ? epochs : [epochs]
    traj = Trajectory(sys, ts)
    PlanetOrbits.frame_pass!(traj, sys.frame)
    PlanetOrbits.propagate!(traj, sys, KeplerianApprox())
    # `observe_pass!` owns the `cart2angle`, `d_au` and barycentre-velocity
    # columns, and they are meaningless without it: an angular observable is
    # about an observation, which is exactly what a raw state is not. Poison
    # them rather than synthesising a shared-scale stand-in — a NaN
    # propagating out of `raoff` is a loud, correct failure, whereas a
    # plausible v1-style value would quietly test a mode that no longer
    # exists. Use posx/posy/posz and velx/vely/velz on raw states.
    fill!(traj.d_au, NaN)
    fill!(traj.bvx, NaN); fill!(traj.bvy, NaN); fill!(traj.bvz, NaN)
    fill!(traj.cart2angle, NaN)
    return epochs isa AbstractVector ? traj : traj[1]
end

# Build the v2 equivalent of a v1 fixture case. Returns (sys, refs).
# `Mp` splits the total mass between primary and secondary so that reflex
# quantities can be tested; the row's total mass — and hence the relative
# orbit — is identical either way.
function fixture_system(c; Mp=nothing)
    p = c.params
    m2 = isnothing(Mp) ? 0.0 : Mp
    A = Body(mass=p.M - m2, name=:A)
    b = Body(mass=m2, name=:b)
    framekw = (;)
    if haskey(p, :plx)
        framekw = (; plx=p.plx)
    end
    if haskey(p, :ra)
        framekw = (; plx=p.plx, ra=p.ra, dec=p.dec, pmra=p.pmra, pmdec=p.pmdec,
                     rv=p.rv, ref_epoch=p.ref_epoch)
    end
    sys = System(Orbit(b, about=A; a=p.a, e=p.e, i=p.i, ω=p.ω, Ω=p.Ω, tp=p.tp); framekw...)
    return sys, bodies(sys)
end

@testset "v1 regression fixtures" begin
    for c in V1_REFERENCE
        @testset "$(c.name)" begin
            # v2 applies four observing-geometry corrections v1 did not (see
            # src/observe.jl), plus the barycentric light-travel sign fix, so
            # these fixtures deliberately no longer agree bit-for-bit. The
            # tolerances below are the *measured* deviations rounded up; the
            # exactness gate for the new physics is "observing geometry vs
            # brute-force 3D reference" further down, not this testset.
            #
            # `kep-face-on` is the exception and stays an exact v1 gate: the
            # orbit is face-on, so z ≡ 0 at every epoch and all four
            # corrections vanish identically.
            tol = c.name == "kep-face-on" ? 1e-11 :
                  c.kind === :absvis ? 3e-2 : 3e-3
            ap(x, y) = isapprox(x, y; rtol=tol, atol=1e-8)
            sys, refs = fixture_system(c)
            @test ap(period(sys), c.period)
            traj = orbitsolve(sys, c.epochs)
            d = c.data
            # Guard against the loose tolerance masking a *deleted*
            # correction: wherever z ≢ 0 the departure from v1 must actually
            # be there, not merely small.
            if c.name != "kep-face-on"
                @test maximum(abs(posx(traj[k]) - d.posx[k]) for k in eachindex(c.epochs)) > 1e-12
            end
            for (k, sol) in enumerate(traj)
                @test ap(posx(sol), d.posx[k])
                @test ap(posy(sol), d.posy[k])
                @test ap(posz(sol), d.posz[k])
                @test ap(velx(sol), d.velx[k])
                @test ap(vely(sol), d.vely[k])
                @test ap(velz(sol), d.velz[k])
            end
            if c.kind == :kep
                for (k, sol) in enumerate(traj)
                    @test ap(kinrv(sol), d.radvel[k])
                end
            end
            if c.kind in (:visual, :absvis)
                for (k, sol) in enumerate(traj)
                    @test ap(raoff(sol), d.raoff[k])
                    @test ap(decoff(sol), d.decoff[k])
                    @test ap(projectedseparation(sol), d.projectedseparation[k])
                    @test ap(posangle(sol), d.posangle[k])
                end
            end
            if c.kind == :visual
                for (k, sol) in enumerate(traj)
                    @test ap(kinrv(sol), d.radvel[k])
                    @test ap(pmra(sol), d.pmra[k])
                    @test ap(pmdec(sol), d.pmdec[k])
                end
            end
            if c.kind == :absvis
                # v1 adds the propagated-frame drift to pm and rv (but not to
                # the position offsets); in v2 that composition is explicit.
                p = c.params
                for (k, sol) in enumerate(traj)
                    @test ap(kinrv(sol) + (frame_rv(sol) - p.rv), d.radvel[k])
                    @test ap(pmra(sol) + (frame_pmra(sol) - p.pmra), d.pmra[k])
                    @test ap(pmdec(sol) + (frame_pmdec(sol) - p.pmdec), d.pmdec[k])
                    @test ap(frame_ra(sol), d.comp_ra2[k])
                    @test ap(frame_dec(sol), d.comp_dec2[k])
                    @test ap(frame_pmra(sol), d.comp_pmra2[k])
                    @test ap(frame_pmdec(sol), d.comp_pmdec2[k])
                    @test ap(frame_rv(sol), d.comp_rv2[k])
                end
            end
            # Reflex: same case with the mass split between the two bodies.
            if !isnothing(c.Mp)
                sysr, refsr = fixture_system(c; Mp=c.Mp)
                bary = barycentre(sysr)
                trajr = orbitsolve(sysr, c.epochs)
                p = c.params
                for (k, sol) in enumerate(trajr)
                    # relative quantities are unchanged by the mass split
                    @test ap(raoff(sol, refsr.b, refsr.A), d.raoff[k])
                    @test ap(raoff(sol, refsr.A, bary), d.raoff_reflex[k])
                    @test ap(decoff(sol, refsr.A, bary), d.decoff_reflex[k])
                    if c.kind == :visual
                        @test ap(pmra(sol, refsr.A, bary), d.pmra_reflex[k])
                        @test ap(pmdec(sol, refsr.A, bary), d.pmdec_reflex[k])
                        @test ap(kinrv(sol, refsr.A, bary), d.radvel_reflex[k])
                    else
                        @test ap(pmra(sol, refsr.A, bary) + (frame_pmra(sol) - p.pmra), d.pmra_reflex[k])
                        @test ap(pmdec(sol, refsr.A, bary) + (frame_pmdec(sol) - p.pmdec), d.pmdec_reflex[k])
                        @test ap(kinrv(sol, refsr.A, bary) + (frame_rv(sol) - p.rv), d.radvel_reflex[k])
                    end
                end
            end
        end
    end
end

@testset "physical invariants" begin
    # 3-body Jacobi chain: star + two planets
    A = Body(mass=1.1, name=:A)
    b = Body(mass=8mjup, name=:b)
    c = Body(mass=2mjup, name=:c)
    sys = System((A, b, c), (
        Orbit(b, about=A;      a=2.5, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
        Orbit(c, about=(A, b); a=8.0, e=0.3, i=0.6, ω=0.4, Ω=2.0, tp=57000.0)); plx=50.0)
    traj = orbitsolve(sys, [58000.0, 59000.0, 60000.0])
    refs = bodies(sys)
    bary = barycentre(sys)
    for sol in traj
        # barycentre is the origin of the propagation frame: mass-weighted
        # position and momentum vanish
        for f in (posx, posy, posz, velx, vely, velz)
            @test abs(sum(sys.masses[j] * f(sol, BodyRef(j), bary) for j in 1:3) /
                      sum(sys.masses)) < 1e-12
        end
        # pairwise antisymmetry
        @test raoff(sol, refs.b, refs.A) ≈ -raoff(sol, refs.A, refs.b)
    end

    # circumbinary: planet around a tight binary
    B1 = Body(mass=0.6, name=:Aa)
    B2 = Body(mass=0.4, name=:Ab)
    p = Body(mass=1mjup, name=:b)
    csys = System((B1, B2, p), (
        Orbit(B2, about=B1;       a=0.2, e=0.05, i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
        Orbit(p,  about=(B1, B2); a=3.0, e=0.1,  i=0.4, ω=1.0, Ω=2.0, tp=58000.0)); plx=80.0)
    ctraj = orbitsolve(csys, [58900.0])
    crefs = bodies(csys)
    sol = ctraj[1]
    innerbary = barycentre(csys, crefs.Aa, crefs.Ab)
    # the two binary members are on opposite sides of their barycentre,
    # scaled by the mass ratio
    @test posx(sol, crefs.Aa, innerbary) ≈ -(0.4 / 0.6) * posx(sol, crefs.Ab, innerbary)
    # Planet's separation from the inner barycentre solves the outer row.
    # This is a property of the A⁻¹ combine, so it is asserted on the raw
    # propagated states: `observe_pass!` afterwards retards each body to its
    # own emission time, which deliberately breaks the equality at the v/c
    # level (the row scratch columns are never retarded).
    rawtraj = Trajectory(csys, [58900.0])
    PlanetOrbits.frame_pass!(rawtraj, csys.frame)
    PlanetOrbits.propagate!(rawtraj, csys, KeplerianApprox())
    rawsol = rawtraj[1]
    @test approx(
        hypot(posx(rawsol, crefs.b, innerbary), posy(rawsol, crefs.b, innerbary),
              posz(rawsol, crefs.b, innerbary)),
        hypot(rawsol.traj.rx[1, 2], rawsol.traj.ry[1, 2], rawsol.traj.rz[1, 2]))
    # ...and the retarded states depart from it by exactly the expected order.
    @test !approx(
        hypot(posx(sol, crefs.b, innerbary), posy(sol, crefs.b, innerbary),
              posz(sol, crefs.b, innerbary)),
        hypot(sol.traj.rx[1, 2], sol.traj.ry[1, 2], sol.traj.rz[1, 2]))

    # zero-mass companion degrades gracefully: star sits at the barycentre
    zsys = System(Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A);
        a=5.0, e=0.2, i=0.4, ω=1.0, Ω=0.5, tp=58849.0); plx=40.0)
    zrefs = bodies(zsys)
    ztraj = orbitsolve(zsys, [58900.0])
    @test abs(raoff(ztraj[1], zrefs.A, barycentre(zsys))) < 1e-13
end

@testset "photocentre" begin
    A = Body(mass=1.0, flux=(G=1.0,), name=:A)
    b = Body(mass=0.2, flux=(G=0.0,), name=:b)
    sys = System(Orbit(b, about=A; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    refs = bodies(sys)
    sol = orbitsolve(sys, [59000.0])[1]
    # dark companion: photocentre is the star
    pc = photocentre(sys)
    @test raoff(sol, pc, barycentre(sys)) ≈ raoff(sol, refs.A, barycentre(sys))
    # equal-brightness pair: photocentre at the midpoint
    A2 = Body(mass=1.0, flux=(G=1.0,), name=:A)
    b2 = Body(mass=0.2, flux=(G=1.0,), name=:b)
    sys2 = System(Orbit(b2, about=A2; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    refs2 = bodies(sys2)
    sol2 = orbitsolve(sys2, [59000.0])[1]
    pc2 = photocentre(sys2)
    mid = (raoff(sol2, refs2.A, barycentre(sys2)) + raoff(sol2, refs2.b, barycentre(sys2))) / 2
    @test raoff(sol2, pc2, barycentre(sys2)) ≈ mid
    # multi-band requires selection
    A3 = Body(mass=1.0, flux=(G=1.0, K=0.5), name=:A)
    b3 = Body(mass=0.2, flux=(K=0.4,), name=:b)
    sys3 = System(Orbit(b3, about=A3; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    @test_throws ErrorException photocentre(sys3)
    @test photocentre(sys3; band=:G).w[2] == 0
    @test photocentre(sys3; band=:K).w[2] ≈ 0.4 / 0.9
    # `fluxes` is the public read of what the bodies declared, in body order;
    # a band no body declares is zero, not missing.
    @test fluxes(sys3, :G) == SVector(1.0, 0.0)
    @test fluxes(sys3, :K) == SVector(0.5, 0.4)
    @test keys(fluxes(sys3)) === (:G, :K)
    @test fluxes(sys3).K === fluxes(sys3, :K)
    @test_throws "no band :H" fluxes(sys3, :H)
    # a system with no fluxes at all says so, and points at `Body(flux=…)`
    nofl = System(Orbit(Body(mass=0.2, name=:b), about=Body(mass=1.0, name=:A);
        a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    @test isempty(fluxes(nofl))
    @test_throws "no body in this system declares a flux" fluxes(nofl, :G)
    @test_throws "no fluxes defined" photocentre(nofl)
    # an all-dark band is a loud failure, not a NaN point
    dark = System(Orbit(Body(mass=0.2, flux=(G=0.0,), name=:b),
            about=Body(mass=1.0, flux=(G=0.0,), name=:A);
            a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    @test_throws "total flux is zero" photocentre(dark)
end

@testset "body fluxes participate in the type promotion" begin
    # A `flux_<band>` variable may be sampled, or derived from a sampled one, so
    # under a gradient-based sampler it arrives as a `ForwardDiff.Dual`. Before
    # this was fixed, `System`'s promotion covered masses, elements and the frame
    # scalar but not fluxes, so `_collect_fluxes` narrowed the Dual back to
    # Float64 and the model failed to build -- which blocked sampled contrasts
    # for photometry, images and interferometry alike.
    function offset_from_photocentre(x)
        A = Body(mass=1.0, flux=(; H=x[1]), name=:A)
        b = Body(mass=0.01, flux=(; H=x[2]), name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=3.0, e=0.1, i=0.5,
            ω=0.2, Ω=0.3, tp=55000.0),); plx=25.0)
        tr = orbitsolve(sys, [56000.0])
        return raoff(tr[1], BodyRef(2), photocentre(sys; band=:H))
    end

    x0 = [1.0, 0.25]
    @test offset_from_photocentre(x0) isa Float64

    # The system builds at all with a Dual flux, and carries it through.
    dualsys = let x = ForwardDiff.Dual.(x0, 1.0, 0.0)
        A = Body(mass=1.0, flux=(; H=x[1]), name=:A)
        b = Body(mass=0.01, flux=(; H=x[2]), name=:b)
        System((A, b), (Orbit(b, about=A; a=3.0, e=0.1, i=0.5,
            ω=0.2, Ω=0.3, tp=55000.0),); plx=25.0)
    end
    @test eltype(fluxes(dualsys, :H)) <: ForwardDiff.Dual

    # And the gradient is right, not merely finite. For two bodies the
    # photocentre offset is f_A·Δ/(f_A + f_b), so
    #   ∂/∂f_A =  f_b·Δ/(f_A + f_b)²,   ∂/∂f_b = −f_A·Δ/(f_A + f_b)²
    g = ForwardDiff.gradient(offset_from_photocentre, x0)
    fA, fb = x0
    Δ = offset_from_photocentre(x0) * (fA + fb) / fA
    @test g[1] ≈ fb * Δ / (fA + fb)^2 rtol = 1e-12
    @test g[2] ≈ -fA * Δ / (fA + fb)^2 rtol = 1e-12

    # Mixed declared/undeclared fluxes must promote too: a body with no flux
    # contributes a neutral type rather than pinning the system to Float64.
    mixed(x) = begin
        A = Body(mass=1.0, flux=(; H=x), name=:A)
        b = Body(mass=0.01, name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=3.0, e=0.1, i=0.5,
            ω=0.2, Ω=0.3, tp=55000.0),); plx=25.0)
        return sum(fluxes(sys, :H))
    end
    @test ForwardDiff.derivative(mixed, 2.0) ≈ 1.0
end

@testset "subset photocentre" begin
    # Two pairs, four different fluxes: everything below is checked against
    # the member fluxes by hand rather than against another code path.
    Aa = Body(mass=0.6, flux=(G=1.0, K=0.2), name=:Aa)
    Ab = Body(mass=0.4, flux=(G=0.25,), name=:Ab)
    Ba = Body(mass=0.8, flux=(G=0.8,), name=:Ba)
    Bb = Body(mass=0.7, flux=(G=0.5,), name=:Bb)
    sys = System((Aa, Ab, Ba, Bb), (
        Orbit(Ab, about=Aa; a=0.5, e=0.1, i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
        Orbit(Bb, about=Ba; a=0.6, e=0.2, i=0.4, ω=0.3, Ω=0.2, tp=58849.0),
        Orbit((Ba, Bb), about=(Aa, Ab); a=50.0, e=0.3, i=0.5, ω=0.4, Ω=0.3, tp=58000.0),
    ); plx=20.0)
    refs = bodies(sys)
    sol = orbitsolve(sys, 59000.0)

    # weights: members share f_j / Σ_members f_k, non-members are exactly zero
    wA = photocentre(sys, refs.Aa, refs.Ab; band=:G).w
    @test wA ≈ SVector(1.0 / 1.25, 0.25 / 1.25, 0.0, 0.0)
    @test wA[3] === 0.0 && wA[4] === 0.0
    @test sum(wA) ≈ 1
    wB = photocentre(sys, refs.Ba, refs.Bb; band=:G).w
    @test wB ≈ SVector(0.0, 0.0, 0.8 / 1.3, 0.5 / 1.3)
    # the whole-system photocentre is the all-members subset
    @test photocentre(sys, refs.Aa, refs.Ab, refs.Ba, refs.Bb; band=:G).w ≈
          photocentre(sys; band=:G).w

    # member resolution: BodyRef / named Body value / Symbol are interchangeable
    @test photocentre(sys, Aa, Ab; band=:G).w == wA
    @test photocentre(sys, :Aa, :Ab; band=:G).w == wA
    @test photocentre(sys, Body(mass=99.9, flux=(G=42.0,), name=:Aa), refs.Ab; band=:G).w == wA
    @test_throws "no body named :nope" photocentre(sys, :nope; band=:G)
    @test_throws ErrorException photocentre(sys, Body(mass=1.0, flux=(G=1.0,)); band=:G)

    # band selection errors, on the subset path too
    @test_throws ErrorException photocentre(sys, refs.Aa, refs.Ab)      # two bands
    @test_throws "no band :H" photocentre(sys, refs.Aa, refs.Ab; band=:H)
    # …and a subset that is dark in the selected band. Structural membership
    # over bodies that are all dark is not a point on the sky.
    @test_throws "total flux is zero" photocentre(sys, refs.Ab, refs.Ba; band=:K)

    # a single-member subset is observably the bare body
    @test photocentre(sys, refs.Aa; band=:G).w == SVector(1.0, 0.0, 0.0, 0.0)
    @test raoff(sol, photocentre(sys, refs.Aa; band=:G), refs.Ba) ===
          raoff(sol, refs.Aa, refs.Ba)
    # …even when that body is the only one lit in its band
    @test raoff(sol, photocentre(sys, refs.Aa, refs.Ab; band=:K), refs.Ba) ≈
          raoff(sol, refs.Aa, refs.Ba)

    # DEFINITIONAL CHECK: the subset photocentre's offset is the flux-weighted
    # mean of its members' offsets, against any reference.
    for R in (refs.Ba, barycentre(sys), photocentre(sys; band=:G))
        fa, fb = 1.0, 0.25
        @test raoff(sol, photocentre(sys, refs.Aa, refs.Ab; band=:G), R) ≈
              (fa * raoff(sol, refs.Aa, R) + fb * raoff(sol, refs.Ab, R)) / (fa + fb) rtol = 1e-13
        @test decoff(sol, photocentre(sys, refs.Aa, refs.Ab; band=:G), R) ≈
              (fa * decoff(sol, refs.Aa, R) + fb * decoff(sol, refs.Ab, R)) / (fa + fb) rtol = 1e-13
        @test radvel(sol, photocentre(sys, refs.Aa, refs.Ab; band=:G), R) ≈
              (fa * radvel(sol, refs.Aa, R) + fb * radvel(sol, refs.Ab, R)) / (fa + fb) rtol = 1e-13
    end

    # A hand-built WeightedPoint goes anywhere a ref goes, and `photocentre(w)`
    # is the normalizing constructor for the per-epoch/per-draw pattern.
    member = SVector(1.0, 1.0, 0.0, 0.0)
    wp = photocentre(fluxes(sys, :G) .* member)
    @test wp isa WeightedPoint
    @test wp.w ≈ wA
    @test raoff(sol, wp, refs.Ba) ≈ raoff(sol, photocentre(sys, refs.Aa, refs.Ab; band=:G), refs.Ba)
    @test photocentre(SVector(2.0, 2.0, 2.0, 2.0)).w == SVector(0.25, 0.25, 0.25, 0.25)
    @test_throws "sum to zero" photocentre(SVector(0.0, 0.0, 0.0, 0.0))
    # an explicitly-built WeightedPoint is accepted verbatim (no renormalizing)
    @test raoff(sol, WeightedPoint(SVector(1.0, 0.0, 0.0, 0.0); emits=true), refs.Ba) ===
          raoff(sol, refs.Aa, refs.Ba)
    # `emits` has no default: it is invisible in every observable but `radvel`,
    # so a hand-built point must state it rather than inherit a guess.
    @test_throws UndefKeywordError WeightedPoint(SVector(1.0, 0.0, 0.0, 0.0))
    @test PlanetOrbits._ein(sol, WeightedPoint(SVector(1.0, 0.0, 0.0, 0.0); emits=false)) == 0
    @test PlanetOrbits._ein(sol, WeightedPoint(SVector(1.0, 0.0, 0.0, 0.0); emits=true)) ==
          PlanetOrbits._ein(sol, refs.Aa)

    # inference and isbits-ness: this is meant to be built inside a scan loop
    @test isbits(wp)
    @inferred photocentre(sys, refs.Aa, refs.Ab; band=:G)
    @inferred photocentre(sys, :Aa, :Ab; band=:G)
    @inferred fluxes(sys, :G)
end

@testset "subset photocentre: 2+2 quadruple carries every level of motion" begin
    # Two tight pairs on a wide mutual orbit. Each catalog "source" is the
    # photocentre of one pair; the claim under test is that one dot product
    # over *absolute* states carries both the wide-orbit motion and the
    # intra-pair photocentric wobble, with no per-level bookkeeping.
    els_A = (; a=0.5, e=0.1, i=0.3, ω=0.2, Ω=0.1, tp=58849.0)
    els_B = (; a=0.6, e=0.2, i=0.4, ω=0.3, Ω=0.2, tp=58849.0)
    els_W = (; a=50.0, e=0.3, i=0.5, ω=0.4, Ω=0.3, tp=58000.0)
    epochs = collect(range(58000.0, 59800.0, length=64))

    function quad(fAb, fBb)
        Aa = Body(mass=0.6, flux=(G=1.0,), name=:Aa)
        Ab = Body(mass=0.4, flux=(G=fAb,), name=:Ab)
        Ba = Body(mass=0.8, flux=(G=1.0,), name=:Ba)
        Bb = Body(mass=0.7, flux=(G=fBb,), name=:Bb)
        return System((Aa, Ab, Ba, Bb), (
            Orbit(Ab, about=Aa; els_A...),
            Orbit(Bb, about=Ba; els_B...),
            Orbit((Ba, Bb), about=(Aa, Ab); els_W...)); plx=20.0)
    end

    # Dark secondaries: each source degrades to its pair's primary, and the
    # source-to-source track is dominated by the wide orbit.
    dark = quad(0.0, 0.0)
    lit = quad(0.6, 0.45)      # both secondaries luminous
    srcA(s) = photocentre(s, :Aa, :Ab; band=:G)
    srcB(s) = photocentre(s, :Ba, :Bb; band=:G)
    baryA(s) = barycentre(s, :Aa, :Ab)
    baryB(s) = barycentre(s, :Ba, :Bb)

    td = orbitsolve(dark, epochs)
    tl = orbitsolve(lit, epochs)
    for k in eachindex(epochs)
        # dark ⇒ photocentre *is* the primary
        @test raoff(td[k], srcA(dark), :Aa) ≈ 0 atol = 1e-12
        @test raoff(td[k], srcB(dark), :Ba) ≈ 0 atol = 1e-12
    end

    # (1) The wide-orbit motion is present in each source. Compare the
    # source's track against the pair barycentre's track: the difference is
    # the intra-pair wobble alone, and it is bounded by the pair's own size.
    wide = [raoff(tl[k], baryB(lit), baryA(lit)) for k in eachindex(epochs)]
    srcsep = [raoff(tl[k], srcB(lit), srcA(lit)) for k in eachindex(epochs)]
    @test maximum(abs, wide) > 200            # mas: the 50 AU orbit at 20 mas plx
    # the wide orbit dominates the source-to-source separation…
    @test maximum(abs, srcsep .- wide) < 0.2 * maximum(abs, wide)
    # …but the sources are NOT the barycentres: the wobble is really there
    @test maximum(abs, srcsep .- wide) > 1.0

    # (2) The intra-pair wobble is exactly the flux-weighted photocentre
    # offset of the pair, i.e. the source-vs-pair-barycentre excursion is a
    # function of the pair's own orbit only. Its amplitude follows the
    # standard |f/(1+f) − m/M| coefficient of the pair's relative separation.
    # NB the sign is load-bearing: with f/(1+f) < m₂/M the photocentre sits on
    # the *opposite* side of the barycentre from the secondary.
    fA = 0.6
    mA1, mA2 = 0.6, 0.4
    coeff = fA / (1 + fA) - mA2 / (mA1 + mA2)
    @test coeff < 0
    wob = [raoff(tl[k], srcA(lit), baryA(lit)) for k in eachindex(epochs)]
    rel = [raoff(tl[k], :Ab, :Aa) for k in eachindex(epochs)]
    @test maximum(abs, wob .- coeff .* rel) < 1e-9 * maximum(abs, rel)

    # (3) …and the same wobble is superposed on the wide motion, additively,
    # in the source-to-source signal — the decomposition is exact.
    wobB = [raoff(tl[k], srcB(lit), baryB(lit)) for k in eachindex(epochs)]
    @test maximum(abs, srcsep .- (wide .+ wobB .- wob)) < 1e-9 * maximum(abs, wide)

    # (4) Changing only a secondary's flux moves that source and nothing else.
    lit2 = quad(0.6, 0.9)
    t2 = orbitsolve(lit2, epochs)
    @test all(raoff(t2[k], srcA(lit2), baryA(lit2)) ≈ wob[k] for k in eachindex(epochs))
    @test any(!isapprox(raoff(t2[k], srcB(lit2), baryB(lit2)), wobB[k]) for k in eachindex(epochs))
end

@testset "error paths" begin
    A = Body(mass=1.0, name=:A)
    b = Body(mass=0.001, name=:b)
    # e > 1 is supported now (see the "hyperbolic orbits" testset); only the
    # degenerate parabolic case e == 1 is rejected
    @test_throws "parabolic" System(Orbit(b, about=A; a=1.0, e=1.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0))
    # angular observables need a parallax
    sys = System(Orbit(b, about=A; a=1.0, e=0.1, i=0.2, ω=0.0, Ω=0.0, tp=58849.0))
    sol = orbitsolve(sys, 58900.0)
    @test_throws ErrorException raoff(sol)
    @test_throws ErrorException projectedseparation(sol)
    @test posx(sol) isa Float64
    # …but `posangle` does NOT: it is a ratio of sky-plane coordinates, so the
    # distance scaling cancels. Must work without plx, and agree with the
    # same system given one.
    #
    # No longer *exactly*, though: with a parallax each body is divided by its
    # own d+z, and that factor is only common to the two sky components when
    # both references sit at the same depth. Here they do not (the host has a
    # reflex), so the agreement is to roundoff rather than bit-identical, and
    # for a body-vs-barycentre pair it departs at the ~1e-6 rad level.
    withplx = System(Orbit(b, about=A; a=1.0, e=0.1, i=0.2, ω=0.0, Ω=0.0, tp=58849.0); plx=24.5)
    @test posangle(sol, :b, :A) isa Float64
    @test posangle(sol, :b, :A) ≈ posangle(orbitsolve(withplx, 58900.0), :b, :A) rtol = 1e-12
    # partial absolute frame
    @test_throws ErrorException System(Orbit(b, about=A; a=1.0, e=0.1, i=0.2, ω=0.0, Ω=0.0, tp=0.0);
        plx=10.0, ra=45.0)
    # duplicate names
    @test_throws ErrorException System(Orbit(Body(mass=0.1, name=:x), about=Body(mass=1.0, name=:x);
        a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0))
    # one-argument observables on >2 bodies
    m = Body(mass=0.0001, name=:m)
    sys3 = System((A, b, m), (
        Orbit(b, about=A;      a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0),
        Orbit(m, about=(A, b); a=5.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0)); plx=10.0)
    sol3 = orbitsolve(sys3, 58900.0)
    @test_throws ErrorException raoff(sol3)
end

@testset "name-based reference resolution" begin
    A = Body(mass=1.0, name=:A)
    b = Body(mass=2mjup, name=:b)
    sys = System(Orbit(b, about=A; a=3.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=25.0)
    refs = bodies(sys)
    sol = orbitsolve(sys, 58900.0)
    expected = raoff(sol, refs.b, refs.A)
    # Body values, Symbols, and refs are interchangeable in observables
    @test raoff(sol, b, A) === expected
    @test raoff(sol, :b, :A) === expected
    @test raoff(sol, b, refs.A) === expected
    @test radvel(sol, A, barycentre(sys)) === radvel(sol, refs.A, barycentre(sys))
    # resolution reads only the name: a "stale" Body from another sample works
    @test raoff(sol, Body(mass=99.9, name=:b), A) === expected
    # …and in barycentre membership
    @test barycentre(sys, A, b).w == barycentre(sys, refs.A, refs.b).w
    @test barycentre(sys, :A, :b).w == barycentre(sys, refs.A, refs.b).w
    # unnamed bodies and unknown names error clearly
    @test_throws ErrorException raoff(sol, Body(mass=1.0), A)
    @test_throws ErrorException raoff(sol, :nope, :A)
    # resolution is type-stable
    @inferred raoff(sol, b, A)
    @inferred raoff(sol, :b, :A)
    # Orbit refuses references into an existing system (and other non-members)
    @test_throws ErrorException Orbit(refs.b, about=A; a=1.0)
    @test_throws ErrorException Orbit(b, about=barycentre(sys); a=1.0)
    @test_throws ErrorException Orbit(b, about=:A; a=1.0)
    # orbit() two-body sugar is opt-in, not exported
    @test !Base.isexported(PlanetOrbits, :orbit)
end

@testset "soltime contract & trajectory interface" begin
    o = orbit(M=1.1, a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0, plx=24.5)
    epochs = [58849.0, 59000.0, 59500.0]
    traj = orbitsolve(o, epochs)
    @test length(traj) == 3
    for k in eachindex(traj)
        @test soltime(traj[k]) === epochs[k]
    end
    @test collect(soltime(s) for s in traj) == epochs
    # …including under an absolute frame with large light-travel corrections
    o2 = orbit(M=1.26, a=9.5, e=0.1, i=0.5, ω=1.0, Ω=1.0, tp=58000.0,
        plx=21.9, rv=1000.0, ra=3.65, dec=-7.19, pmra=1e7, pmdec=-1e7, ref_epoch=57388.5)
    traj2 = orbitsolve(o2, epochs)
    for k in eachindex(traj2)
        @test soltime(traj2[k]) === epochs[k]
        @test traj2.t_em[k] != epochs[k]  # emission time genuinely differs
    end
end

@testset "SIMD batch path" begin
    # kernel agreement with the scalar Markley solver
    worst = 0.0
    for e in 0.0:0.05:0.95, M in range(-40.0, 40.0, length=401)
        sE, cE = PlanetOrbits.markley_sincosE(M, e)
        E = PlanetOrbits.kepler_solver(M, e, PlanetOrbits.Markley())
        s0, c0 = sincos(E)
        worst = max(worst, abs(sE - s0), abs(cE - c0))
    end
    @test worst ≤ 4e-15

    # full-trajectory agreement: SIMD vs scalar path
    A = Body(mass=1.1, name=:A); b = Body(mass=5mjup, name=:b); c = Body(mass=2mjup, name=:c)
    sys = System((A, b, c), (
        Orbit(b, about=A;      a=2.5, e=0.6, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
        Orbit(c, about=(A, b); a=8.0, e=0.1, i=0.6, ω=0.4, Ω=2.0, tp=57000.0));
        plx=24.5, ra=45.0, dec=-30.0, pmra=100.0, pmdec=-50.0, rv=25e3, ref_epoch=57388.5)
    epochs = collect(range(56000.0, 61000.0, length=307))
    t_simd = orbitsolve(sys, epochs; method=KeplerianApprox(simd=true))
    t_scal = orbitsolve(sys, epochs; method=KeplerianApprox(simd=false))
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test maximum(abs, getfield(t_simd, f) .- getfield(t_scal, f)) ≤ 1e-12
    end

    # the batch loop must actually vectorize (pack lanes): look for vector
    # ops on double in the emitted IR. Under `Pkg.test` bounds checking is
    # forced on (`--check-bounds=yes`), which voids @inbounds and blocks the
    # vectorizer — the check is only meaningful when bounds checks are not
    # forced (e.g. `include("runtests.jl")` or the perf/ harness).
    if Base.JLOptions().check_bounds == 1
        @info "skipping vectorization IR check: bounds checking is forced on"
    else
        io = IOBuffer()
        code_llvm(io, PlanetOrbits.solve_row_simd!,
            (typeof(t_simd), PlanetOrbits.Row{Float64}, Int); optimize=true)
        ir = String(take!(io))
        @test occursin(r"<\d+ x double>", ir)
    end
end

# --- Micro performance gates (§11 of the design doc) ---

function _build_workload_system(θ)
    A = Body(mass=θ[1], name=:A)
    b = Body(mass=θ[2], name=:b)
    return System(Orbit(b, about=A; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
        plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
end

function _eval_workload(θ, epochs, traj=nothing; method=KeplerianApprox())
    sys = _build_workload_system(θ)
    if traj === nothing
        traj = Trajectory{eltype(θ)}(sys, epochs)
    end
    orbitsolve!(traj, sys; method)
    refs = bodies(sys)
    bary = barycentre(sys)
    acc = zero(eltype(θ))
    for sol in traj
        acc += raoff(sol, refs.b, refs.A) + decoff(sol, refs.b, refs.A) +
               pmra(sol, refs.A, bary) + frame_pmra(sol) +
               radvel(sol, refs.A, bary) + frame_rv(sol)
    end
    return acc
end

const θ0 = [1.1, 5mjup, 8.0, 0.1, 0.5, 1.1, 2.2, 58849.0,
            24.5, 45.0, -30.0, 100.0, -50.0, 25e3, 57388.5]

@testset "allocation-free hot path" begin
    epochs = collect(range(58000.0, 60000.0, length=50))
    θ = SVector{15}(θ0)
    A = Body(mass=θ[1], name=:A); b = Body(mass=θ[2], name=:b)
    sys = System(Orbit(b, about=A; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
        plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
    # preallocated trajectory stands in for Bumper-owned buffers
    traj = Trajectory(sys, epochs)
    _eval_workload(θ, epochs, traj)  # warm up
    allocs = @allocated _eval_workload(θ, epochs, traj)
    @test allocs == 0
end

# Static counterpart to the @allocated smoke test above: AllocCheck proves
# allocation-freedom across *all* compiled paths, not just the one executed
# (ignore_throw=true excludes the deliberate error() guard branches).
# NB: solver=Markley() explicitly — with Auto() the compiled-but-unreachable
# hyperbolic fallback would drag allocating Roots.jl machinery into the check.
_ac_build(θ) = System(
    Orbit(Body(mass=θ[2], name=:b), about=Body(mass=θ[1], name=:A);
        a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
    plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
# Through `_solve_serial!`, not the `orbitsolve!` front door: the `threads > 1`
# branch reachable from `orbitsolve!` spawns tasks and builds epoch views,
# which allocate by design (and its docstring says so). The static
# allocation-freedom contract is on the serial solve every call runs when
# `threads == 1` — the default, and the only path the sampler hot loop takes.
_ac_solve!(traj, sys) = PlanetOrbits._solve_serial!(
    traj, sys, KeplerianApprox(solver=PlanetOrbits.Markley()), true, true)
_ac_solve_scalar!(traj, sys) = PlanetOrbits._solve_serial!(
    traj, sys, KeplerianApprox(solver=PlanetOrbits.Markley(), simd=false), true, true)
_ac_query(sol, b, A, w) = raoff(sol, b, A) + decoff(sol, b, A) +
    pmra(sol, A, w) + radvel(sol, A, w) + frame_rv(sol)

# Known-benign static findings: Base.rem2pi's `abs(x) < π` Irrational
# comparison compiles in a BigFloat conversion branch — including its
# setprecision ScopedValue/HAMT plumbing — that is dynamically dead for
# Float64 (verified: @allocated rem2pi(x, RoundNearest) == 0 even at
# x = 1e300). Everything else must be provably allocation-free.
function _ac_benign(err)
    s = sprint(show, err)
    return occursin("rem2pi", s) || occursin("ScopedValue", s) ||
           occursin("MPFR", s) || occursin("BigFloat", s)
end

@testset "static allocation-freedom (AllocCheck)" begin
    θ = SVector{15}(θ0)
    sys = _ac_build(θ)
    traj = Trajectory(sys, [58900.0, 59000.0])
    sol = traj[1]
    w = barycentre(sys)
    bB = Body(mass=0.1, name=:b); AB = Body(mass=1.0, name=:A)
    for (f, types) in (
        (_ac_build, (typeof(θ),)),
        (_ac_solve!, (typeof(traj), typeof(sys))),
        (_ac_solve_scalar!, (typeof(traj), typeof(sys))),
        (_ac_query, (typeof(sol), BodyRef, BodyRef, typeof(w))),
        # same query with Body-value keys: name resolution must fold away
        (_ac_query, (typeof(sol), typeof(bB), typeof(AB), typeof(w))),
    )
        errs = filter(!_ac_benign, AllocCheck.check_allocs(f, types))
        isempty(errs) || display(errs[1])
        @test isempty(errs)
    end
end

@testset "type stability" begin
    epochs = [58900.0, 59000.0]
    A = Body(mass=1.1, name=:A); b = Body(mass=5mjup, name=:b)
    root = Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
    sys = @inferred System(root; plx=24.5, ra=45.0, dec=-30.0, pmra=100.0,
        pmdec=-50.0, rv=25e3, ref_epoch=57388.5)
    @test isbits(sys)
    traj = @inferred orbitsolve(sys, epochs)
    refs = @inferred bodies(sys)
    sol = traj[1]
    @inferred raoff(sol, refs.b, refs.A)
    @inferred radvel(sol, refs.A, barycentre(sys))
    @inferred barycentre(sys)
end

@testset "ForwardDiff gradient" begin
    epochs = collect(range(58000.0, 60000.0, length=20))
    f(θ) = _eval_workload(θ, epochs)
    g_fd = ForwardDiff.gradient(f, θ0)
    g_ref = FiniteDiff.finite_difference_gradient(f, θ0)
    for i in eachindex(θ0)
        @test isapprox(g_fd[i], g_ref[i]; rtol=1e-5, atol=1e-6 * max(1.0, abs(g_fd[i])))
    end
    @test f(θ0) isa Float64
end

# Dual-valued rows solve their *primal* roots through the same
# `markley_sincosE` batch kernel the Float64 path uses, then attach partials
# with the implicit rule, instead of running the scalar solver once per epoch.
@testset "Dual SIMD batch path" begin
    D = ForwardDiff.Dual
    P = ForwardDiff.Partials
    simd = KeplerianApprox(simd=true)
    scal = KeplerianApprox(simd=false)

    # The rule must be bit-identical to the route it replaced — solving through
    # `kepler_solver`'s Dual methods and then calling `sincos` on the Dual root.
    # `===`, not `≈`: this is the same arithmetic, minus one redundant sincos.
    @testset "implicit rule is bit-identical to the two-step route" begin
        N = 8
        for e0 in (0.0, 0.01, 0.3, 0.7, 0.95, 0.99),
            M0 in (-40.0, -3.1, 0.05, 0.7, 2.0, 3.1, 40.0)

            MA = D{Nothing}(M0, P(ntuple(i -> 0.1i, N)))
            ed = D{Nothing}(e0, P(ntuple(i -> 0.03i, N)))
            sref, cref = sincos(PlanetOrbits.kepler_solver(MA, ed, PlanetOrbits.Markley()))
            E = PlanetOrbits.kepler_solver(M0, e0, PlanetOrbits.Markley())
            sE, cE = sincos(E)
            snew, cnew = PlanetOrbits._dual_sincosE(MA, ed, sE, cE)
            @test snew === sref
            @test cnew === cref
        end
    end

    # Routing: only first-order Duals over Float64, with Markley/Auto, e < 1.
    @testset "routing" begin
        epochs = collect(range(58000.0, 60000.0, length=20))
        θd = [D{Nothing}(θ0[i], P(ntuple(j -> Float64(j == i), 15))) for i in 1:15]
        sysd = _build_workload_system(θd)
        trajd = Trajectory{eltype(θd)}(sysd, epochs)
        sysf = _build_workload_system(θ0)
        trajf = Trajectory(sysf, epochs)
        @test PlanetOrbits._use_dual_simd(simd, trajd, sysd.rows[1])
        @test !PlanetOrbits._use_dual_simd(scal, trajd, sysd.rows[1])
        @test !PlanetOrbits._use_dual_simd(simd, trajf, sysf.rows[1])
        # Float64 rows with a differentiated frame are still batched.
        @test PlanetOrbits._use_dual_simd(simd, trajd, sysf.rows[1])
        # Nested Duals (Hessians) fall through to the scalar path, which keeps
        # the implicit rule via `_anomaly_sincos`.
        θdd = [D{Nothing}(θd[i], P(ntuple(j -> zero(eltype(θd)), 2))) for i in 1:15]
        sysdd = _build_workload_system(θdd)
        trajdd = Trajectory{eltype(θdd)}(sysdd, epochs)
        @test !PlanetOrbits._use_dual_simd(simd, trajdd, sysdd.rows[1])
        @test ForwardDiff.hessian(θ -> _eval_workload(θ, epochs), θ0) isa Matrix
    end

    @testset "gradient agrees with the scalar path" begin
        epochs = collect(range(58000.0, 60000.0, length=97))
        g_simd = ForwardDiff.gradient(θ -> _eval_workload(θ, epochs; method=simd), θ0)
        g_scal = ForwardDiff.gradient(θ -> _eval_workload(θ, epochs; method=scal), θ0)
        # Against the gradient *norm*: several components are near-zero by
        # construction (∂/∂ra is ~1e-55 — the readout is invariant to the
        # tangent-point longitude), so an element-wise relative comparison
        # against them measures nothing. See design §10.4.2.
        @test maximum(abs, g_simd .- g_scal) / maximum(abs, g_simd) ≤ 1e-14

        # The batch kernel is shared with the value path, so the value carried
        # through a gradient evaluation is now bit-identical to a plain Float64
        # evaluation rather than merely agreeing to the kernels' ~4e-15.
        θd = [D{Nothing}(θ0[i], P(ntuple(j -> Float64(j == i), 12))) for i in 1:15]
        @test ForwardDiff.value(_eval_workload(θd, epochs; method=simd)) ===
              _eval_workload(θ0, epochs; method=simd)
    end

    # Hyperbolic Dual rows are excluded from the batch (the kernel is bound-orbit
    # only) and must still differentiate correctly on the scalar path.
    @testset "hyperbolic Duals fall back correctly" begin
        epochs = collect(range(58000.0, 60000.0, length=40))
        θh = copy(θ0); θh[4] = 1.4     # e > 1
        sysh = _build_workload_system(θh)
        @test !PlanetOrbits._use_dual_simd(simd, Trajectory(sysh, epochs), sysh.rows[1])
        f(θ) = _eval_workload(θ, epochs)
        g_fd = ForwardDiff.gradient(f, θh)
        g_ref = FiniteDiff.finite_difference_gradient(f, θh)
        @test maximum(abs, g_fd .- g_ref) / maximum(abs, g_fd) ≤ 1e-6
    end

    @testset "allocation-free under Duals" begin
        epochs = collect(range(58000.0, 60000.0, length=50))
        θd = [D{Nothing}(θ0[i], P(ntuple(j -> Float64(j == i), 12))) for i in 1:15]
        sysd = _build_workload_system(θd)
        trajd = Trajectory{eltype(θd)}(sysd, epochs)
        orbitsolve!(trajd, sysd; method=simd)
        @test (@allocated orbitsolve!(trajd, sysd; method=simd)) == 0
    end
end

# ---------------------------------------------------------------------------
# Topology, conventions, and the generalized A⁻¹ (§4, §7 of the design doc)
# ---------------------------------------------------------------------------

# The definitional property of A⁻¹, checked without reference to any closed
# form: given arbitrary per-row relative states ρ, the absolute states
# r = A⁻¹ρ must reproduce every row's relative coordinate (exterior
# barycentre minus interior barycentre) and carry zero total momentum.
# A massless member set has no mass-weighted barycentre; its limit is the
# members' geometric centre.
function ainv_residual(sys)
    m = sys.masses
    NB, NR = length(m), length(sys.specs)
    ρ = [0.37k + 1.1 for k in 1:NR]
    r = sys.Ainv * ρ
    groupbary(mask) = begin
        tot = sum(m[j] * mask[j] for j in 1:NB)
        iszero(tot) ? sum(r[j] * mask[j] for j in 1:NB) / count(mask) :
                      sum(m[j] * r[j] * mask[j] for j in 1:NB) / tot
    end
    worst = 0.0
    for k in 1:NR
        s = sys.specs[k]
        worst = max(worst, abs((groupbary(s.ext) - groupbary(s.int)) - ρ[k]))
    end
    mom = abs(sum(m[j] * r[j] for j in 1:NB)) / sum(m)
    return max(worst, mom)
end

@testset "topology: A⁻¹ definitional property" begin
    A = Body(mass=1.1, name=:A); b = Body(mass=8mjup, name=:b)
    c = Body(mass=2mjup, name=:c); d = Body(mass=0.4, name=:d)
    els(x) = (; a=x, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
    Aa = Body(mass=1.1, name=:Aa); Ab = Body(mass=0.9, name=:Ab)
    Ba = Body(mass=0.8, name=:Ba); Bb = Body(mass=0.7, name=:Bb)
    cases = (
        ("two-body", System((A, b), (Orbit(b, about=A; els(3.0)...),))),
        ("Jacobi 3", System((A, b, c), (Orbit(b, about=A; els(2.5)...),
                                        Orbit(c, about=(A, b); els(8.0)...)))),
        ("astrocentric 3", System((A, b, c), (Orbit(b, about=A; els(2.5)...),
                                              Orbit(c, about=A; els(8.0)...)))),
        ("moon (mixed)", System((A, b, c), (Orbit(b, about=A; els(5.2)...),
                                            Orbit(c, about=b; els(0.02)...)))),
        ("Jacobi 4", System((A, b, c, d), (Orbit(b, about=A; els(2.5)...),
                                           Orbit(c, about=(A, b); els(8.0)...),
                                           Orbit(d, about=(A, b, c); els(20.0)...)))),
        ("astrocentric 4", System((A, b, c, d), (Orbit(b, about=A; els(2.5)...),
                                                 Orbit(c, about=A; els(8.0)...),
                                                 Orbit(d, about=A; els(20.0)...)))),
        ("2+2 quadruple", System((Aa, Ab, Ba, Bb), (
            Orbit(Ab, about=Aa; els(0.5)...),
            Orbit(Bb, about=Ba; els(0.6)...),
            Orbit((Ba, Bb), about=(Aa, Ab); els(50.0)...)))),
        # zero-mass members: the `n_planets`-prior pattern must keep working
        ("zero-mass chain", System(
            (Body(mass=1.0, name=:A), Body(mass=0.0, name=:b), Body(mass=0.0, name=:c)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A); els(3.0)...),
             Orbit(Body(mass=0.0, name=:c),
                   about=(Body(mass=1.0, name=:A), Body(mass=0.0, name=:b)); els(9.0)...)))),
    )
    for (nm, sys) in cases
        @testset "$nm" begin
            @test ainv_residual(sys) < 1e-13
            @test all(isfinite, sys.Ainv)
        end
    end
end

@testset "topology: conventions are distinct, not relabellings" begin
    A = Body(mass=1.1, name=:A); b = Body(mass=80mjup, name=:b); c = Body(mass=20mjup, name=:c)
    e1 = (; a=2.5, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
    e2 = (; a=8.0, e=0.3, i=0.6, ω=0.4, Ω=2.0, tp=57000.0)
    jac = System((A, b, c), (Orbit(b, about=A; e1...), Orbit(c, about=(A, b); e2...)); plx=50.0)
    ast = System((A, b, c), (Orbit(b, about=A; e1...), Orbit(c, about=A;      e2...)); plx=50.0)

    # Under KeplerianApprox the rows *are* the model, so the two must differ —
    # a "generalization" that silently collapsed them would be a real bug.
    @test !(jac.Ainv ≈ ast.Ainv)
    epochs = [58000.0, 59000.0, 60000.0]
    tj = orbitsolve(jac, epochs); ta = orbitsolve(ast, epochs)
    diffs = [abs(raoff(tj[k], :c, :A) - raoff(ta[k], :c, :A)) for k in 1:3]
    @test maximum(diffs) > 1.0           # mas — a real, observable difference
    # and the row masses differ exactly as the conventions prescribe
    @test jac.rows[2].M ≈ 1.1 + 80mjup + 20mjup
    @test ast.rows[2].M ≈ 1.1 + 20mjup

    # Both conventions are reported honestly
    @test PlanetOrbits._system_convention(jac.specs) === :jacobi
    @test PlanetOrbits._system_convention(ast.specs) === :astrocentric

    # A massless intermediate makes barycentre(A,b) ≡ A, so the two
    # conventions describe the same configuration and must agree *exactly*.
    b0 = Body(mass=0.0, name=:b)
    jac0 = System((A, b0, c), (Orbit(b0, about=A; e1...), Orbit(c, about=(A, b0); e2...)); plx=50.0)
    ast0 = System((A, b0, c), (Orbit(b0, about=A; e1...), Orbit(c, about=A;       e2...)); plx=50.0)
    @test jac0.Ainv ≈ ast0.Ainv
    for prop in (KeplerianApprox(), AHL21(h=1.0, t0=58000.0))
        t1 = orbitsolve(jac0, epochs; method=prop)
        t2 = orbitsolve(ast0, epochs; method=prop)
        @test maximum(abs, t1.x .- t2.x) < 1e-12
        @test maximum(abs, t1.vx .- t2.vx) < 1e-12
    end

    # Two-body: only one spelling exists, and A⁻¹ is the analytic reflex split
    two = System((A, b), (Orbit(b, about=A; e1...),))
    M = 1.1 + 80mjup
    @test two.Ainv[1, 1] ≈ -80mjup / M
    @test two.Ainv[2, 1] ≈ 1.1 / M
end

@testset "topology: moons and set exteriors" begin
    # A moon orbiting its host under KeplerianApprox — impossible in v1 and
    # in the nested v2 tree, since the host would have to appear twice.
    A = Body(mass=1.0, name=:A); b = Body(mass=10mjup, name=:b); m = Body(mass=1mearth, name=:m)
    sys = System((A, b, m), (
        Orbit(b, about=A; a=5.2, e=0.05, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
        Orbit(m, about=b; a=0.007, e=0.0, i=0.4, ω=0.0, Ω=1.0, tp=58849.0)); plx=50.0)
    sol = orbitsolve(sys, 58900.0)
    # e = 0, so the moon's 3D separation from its host is a — up to the
    # per-body light-travel retardation, which moves host and moon to
    # slightly different points on their orbits (a v/c-order effect).
    @test hypot(posx(sol, :m, :b), posy(sol, :m, :b), posz(sol, :m, :b)) ≈ 0.007 rtol = 1e-4
    # projected on the sky it is foreshortened by inclination: a·cos(i) … a
    @test 0.007 * cos(0.4) * 50.0 - 1e-9 ≤ projectedseparation(sol, :m, :b) ≤ 0.007 * 50.0 + 1e-9
    @test sys.rows[2].M ≈ 10mjup + 1mearth      # moon row: host + moon only
    @test ainv_residual(sys) < 1e-13

    # 2+2 quadruple via set exteriors: the wide row's endpoints really are
    # the two inner barycentres.
    Aa = Body(mass=1.1, name=:Aa); Ab = Body(mass=0.9, name=:Ab)
    Ba = Body(mass=0.8, name=:Ba); Bb = Body(mass=0.7, name=:Bb)
    q = System((Aa, Ab, Ba, Bb), (
        Orbit(Ab, about=Aa; a=0.5, e=0.1, i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
        Orbit(Bb, about=Ba; a=0.6, e=0.2, i=0.4, ω=0.3, Ω=0.2, tp=58849.0),
        Orbit((Ba, Bb), about=(Aa, Ab); a=50.0, e=0.3, i=0.5, ω=0.4, Ω=0.3, tp=58000.0)); plx=20.0)
    bA = barycentre(q, :Aa, :Ab); bB = barycentre(q, :Ba, :Bb)
    # Asserted on the raw propagated states: the row scratch columns are
    # never retarded, so comparing them against observed states would fail at
    # the v/c level by construction (see "physical invariants").
    qs = Trajectory(q, [58900.0])
    PlanetOrbits.frame_pass!(qs, q.frame)
    PlanetOrbits.propagate!(qs, q, KeplerianApprox())
    qs = qs[1]
    @test hypot(posx(qs, bB, bA), posy(qs, bB, bA), posz(qs, bB, bA)) ≈
          hypot(qs.traj.rx[1, 3], qs.traj.ry[1, 3], qs.traj.rz[1, 3])
    @test q.rows[3].M ≈ 3.5
    @test ainv_residual(q) < 1e-13
end

@testset "topology: validation errors name the offending row" begin
    A = Body(mass=1.1, name=:A); b = Body(mass=8mjup, name=:b)
    c = Body(mass=2mjup, name=:c); d = Body(mass=0.5, name=:d)
    @test_throws "needs exactly 2 orbits" System((A, b, c), (Orbit(b, about=A; a=1.0),))
    @test_throws "on both sides" System((A, b), (Orbit(b, about=(A, b); a=1.0),))
    @test_throws "orbits 1 and 2 are the same relationship" System(
        (A, b, c), (Orbit(b, about=A; a=1.0), Orbit(b, about=A; a=2.0)))
    @test_throws "opposite directions" System(
        (A, b, c), (Orbit(c, about=(A, b); a=1.0), Orbit((A, b), about=c; a=2.0)))
    @test_throws "does not appear in any orbit" System(
        (A, b, c, d), (Orbit(b, about=A; a=1.0), Orbit(c, about=A; a=2.0),
                       Orbit(c, about=(A, b); a=3.0)))
    @test_throws "must be unique" System((A, A), (Orbit(A, about=A; a=1.0),))
end

@testset "size group: a | P" begin
    A = Body(mass=1.1, name=:A); b = Body(mass=8mjup, name=:b)
    sysa = System((A, b), (Orbit(b, about=A; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),))
    # P is in DAYS and round-trips exactly through period(sys)
    sysp = System((A, b), (Orbit(b, about=A; P=period(sysa), e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),))
    @test semimajoraxis(sysp) ≈ 8.0 rtol = 1e-14
    @test period(sysp) ≈ period(sysa) rtol = 1e-14
    @test_throws "got neither" Orbit(b, about=A; e=0.1)
    @test_throws "got both" Orbit(b, about=A; a=1.0, P=365.0)
    # P uses the row's own gravitating mass, so the same P under different
    # conventions gives different a
    c = Body(mass=200mjup, name=:c)
    jr = Orbit(c, about=(A, b); P=5000.0)
    ar = Orbit(c, about=A; P=5000.0)
    @test jr.a > ar.a

    # M= override: labelled compatibility, changes the row mass verbatim
    ov = System((A, b), (Orbit(b, about=A; a=8.0, M=2.5, e=0.1),))
    @test ov.rows[1].M == 2.5
    @test occursin("M= override", sprint(show, MIME"text/plain"(), ov))
end

# Regression gate for the pre-existing A⁻¹ cliff: building a 5-body system
# used to cost 10.2 µs and 26 kB (a flat ntuple(Val(NB*NR)) over a
# heterogeneous rows tuple, which fell off the heap-allocation threshold at
# NB=5). Fixed-width masks brought it to ~0.24 µs and 0 bytes.
_build5(m) = System(
    (Body(mass=m[1], name=:A), Body(mass=m[2], name=:b), Body(mass=m[3], name=:c),
     Body(mass=m[4], name=:d), Body(mass=m[5], name=:e)),
    (Orbit(Body(mass=m[2], name=:b), about=Body(mass=m[1], name=:A); a=2.5),
     Orbit(Body(mass=m[3], name=:c),
           about=(Body(mass=m[1], name=:A), Body(mass=m[2], name=:b)); a=5.0),
     Orbit(Body(mass=m[4], name=:d),
           about=(Body(mass=m[1], name=:A), Body(mass=m[2], name=:b),
                  Body(mass=m[3], name=:c)); a=9.0),
     Orbit(Body(mass=m[5], name=:e),
           about=(Body(mass=m[1], name=:A), Body(mass=m[2], name=:b),
                  Body(mass=m[3], name=:c), Body(mass=m[4], name=:d)); a=17.0)))

@testset "many-body construction stays allocation-free" begin
    m5 = SVector(1.1, 5mjup, 3mjup, 2mjup, 1mjup)
    _build5(m5)
    @test (@allocated _build5(m5)) == 0
    @test ainv_residual(_build5(m5)) < 1e-13
    # …and under Duals, where the old code allocated from NB=3 upward
    md = SVector(ForwardDiff.Dual(1.1, 1.0, 0.0), ForwardDiff.Dual(5mjup, 0.0, 1.0),
                 ForwardDiff.Dual(3mjup, 0.0, 0.0), ForwardDiff.Dual(2mjup, 0.0, 0.0),
                 ForwardDiff.Dual(1mjup, 0.0, 0.0))
    _build5(md)
    @test (@allocated _build5(md)) == 0
end

# The subset-photocentre query path, end to end: construct → solve → query,
# at NB = 5 and under Duals. NB ≥ 5 is the point — a heterogeneous-tuple
# constant-folding failure is a cliff, not a slope, and every gate that used
# a 2-body system missed the last one (§12). Two of the five bodies are
# blended into one source, which is the shape a catalog likelihood has.
function _build_photo_sys(θ)
    A = Body(mass=θ[1], flux=(G=θ[6],), name=:A)
    b = Body(mass=θ[2], flux=(G=θ[7],), name=:b)
    c = Body(mass=θ[3], flux=(G=θ[8],), name=:c)
    d = Body(mass=θ[4], flux=(G=θ[9],), name=:d)
    e = Body(mass=θ[5], flux=(G=θ[10],), name=:e)
    return System((A, b, c, d, e), (
            Orbit(b, about=A; a=θ[11], e=θ[12], i=0.3, ω=0.2, Ω=0.1, tp=58849.0),
            Orbit(c, about=(A, b); a=θ[13], e=0.05, i=0.4, ω=0.3, Ω=0.2, tp=58849.0),
            Orbit(d, about=(A, b, c); a=θ[14], e=0.1, i=0.5, ω=0.4, Ω=0.3, tp=58000.0),
            Orbit(e, about=(A, b, c, d); a=θ[15], e=0.2, i=0.6, ω=0.5, Ω=0.4, tp=58000.0));
        plx=θ[16])
end

# solve → query only; `_photo_workload` adds the construction on top.
function _photo_query(sys, traj)
    orbitsolve!(traj, sys; method=KeplerianApprox(solver=PlanetOrbits.Markley()))
    src = photocentre(sys, :A, :b; band=:G)   # the blended source
    bary = barycentre(sys)
    acc = zero(eltype(sys.masses))
    for sol in traj
        acc += raoff(sol, src, bary) + decoff(sol, src, bary) + radvel(sol, src, bary)
    end
    return acc
end

function _photo_workload(θ, epochs, traj=nothing)
    sys = _build_photo_sys(θ)
    if traj === nothing
        traj = Trajectory{eltype(θ)}(sys, epochs)
    end
    return _photo_query(sys, traj)
end

const θ_photo = [1.1, 0.3, 5mjup, 3mjup, 2mjup,      # masses
                 1.0, 0.4, 0.02, 0.0, 0.0,           # G fluxes
                 2.5, 0.15, 6.0, 12.0, 25.0,         # a's and one e
                 24.5]                               # plx

@testset "subset photocentre query is allocation-free (NB=5, Float64 + Dual)" begin
    epochs = collect(range(58000.0, 60000.0, length=50))
    θ = SVector{16}(θ_photo)
    @test isfinite(_photo_workload(θ, epochs))   # warm up, and it must run
    traj = Trajectory(_build_photo_sys(θ), epochs)
    _photo_workload(θ, epochs, traj)
    @test (@allocated _photo_workload(θ, epochs, traj)) == 0

    D = ForwardDiff.Dual
    P = ForwardDiff.Partials
    θd = SVector{16}([D{Nothing}(θ_photo[i], P(ntuple(j -> Float64(j == i), 12))) for i in 1:16])
    sysd = _build_photo_sys(θd)
    trajd = Trajectory{eltype(θd)}(sysd, epochs)
    # Solve → query only under Duals. `System` construction itself allocates
    # 14–33 kB at NB = 5 once the chunk width reaches ~6 (0 bytes at Dual{2},
    # 0 for Float64 at any width, and unaffected by fluxes or by the frame) —
    # measured identical on the branch point, so it is a pre-existing hole
    # that the NB=5 gate above never saw because it uses Dual{2}. Not this
    # change's to fix; flagged separately.
    _photo_query(sysd, trajd)
    @test (@allocated _photo_query(sysd, trajd)) == 0
    # the gradient is right, not merely cheap
    f(x) = _photo_workload(x, epochs)
    g_fd = ForwardDiff.gradient(f, collect(θ_photo))
    g_ref = FiniteDiff.finite_difference_gradient(f, collect(θ_photo))
    @test maximum(abs, g_fd .- g_ref) / maximum(abs, g_fd) ≤ 1e-6
end

# Static counterpart, so the gate covers every compiled path rather than the
# one executed — including the Symbol member resolution folding to indices.
_ac_photo_query(sys, sol) = begin
    src = photocentre(sys, :A, :b; band=:G)
    raoff(sol, src, barycentre(sys)) + decoff(sol, src, photocentre(sys; band=:G))
end

@testset "subset photocentre static allocation-freedom (AllocCheck)" begin
    θ = SVector{16}(θ_photo)
    sys = _build_photo_sys(θ)
    traj = Trajectory(sys, [58900.0, 59000.0])
    θd = SVector{16}([ForwardDiff.Dual{Nothing}(θ_photo[i],
        ForwardDiff.Partials(ntuple(j -> Float64(j == i), 12))) for i in 1:16])
    sysd = _build_photo_sys(θd)
    trajd = Trajectory{eltype(θd)}(sysd, [58900.0, 59000.0])
    for (f, types) in (
        (_ac_photo_query, (typeof(sys), typeof(traj[1]))),
        (_ac_photo_query, (typeof(sysd), typeof(trajd[1]))),
        (_build_photo_sys, (typeof(θ),)),
    )
        errs = filter(!_ac_benign, AllocCheck.check_allocs(f, types))
        isempty(errs) || display(errs[1])
        @test isempty(errs)
    end
end

@testset "hyperbolic orbits (e > 1)" begin
    # Solver: residual of e·sinh(H) − H = M over a wide (M, e) grid
    worst = 0.0
    for e in (1.0000001, 1.001, 1.01, 1.1, 1.5, 2.0, 5.0, 20.0, 100.0),
        MA in (-1e6, -1e3, -10.0, -1.0, -1e-8, 0.0, 1e-8, 1.0, 10.0, 1e3, 1e6)
        H = PlanetOrbits.kepler_solver(MA, e, PlanetOrbits.HyperbolicHalley())
        worst = max(worst, abs(e * sinh(H) - H - MA) / max(abs(MA), 1.0))
    end
    @test worst < 1e-13
    @test (PlanetOrbits.kepler_solver(3.0, 2.0, PlanetOrbits.HyperbolicHalley());
           @allocated PlanetOrbits.kepler_solver(3.0, 2.0, PlanetOrbits.HyperbolicHalley())) == 0
    # e == 1 is degenerate in element space and must say so
    @test_throws "parabolic" System((Body(mass=1.0, name=:A), Body(mass=0.0, name=:b)),
        (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A); a=-5.0, e=1.0),))

    # Physics: vis-viva, and conservation of energy and angular momentum.
    # These are independent of the implementation — v1's hyperbolic branch
    # set the velocity semiamplitude to zero and would fail all three.
    for (a, e, Mtot) in ((-5.0, 1.2, 1.0), (-2.0, 3.0, 1.1), (-50.0, 1.05, 0.8))
        A = Body(mass=Mtot, name=:A); b = Body(mass=0.0, name=:b)
        sys = System((A, b), (Orbit(b, about=A; a=a, e=e, i=0.5, ω=1.1, Ω=2.2, tp=59000.0),))
        @test period(sys) == Inf
        @test semimajoraxis(sys) < 0
        # AU³ / julian-yr² — the unit velx/vely/velz are in. A bare 4π² here
        # is GM per *kepler* year and silently misses by 1.9e-5, which is
        # exactly the error this test failed to catch the first time.
        μ = PlanetOrbits.GM_sun_au3_julianyr2 * Mtot
        Es = Float64[]; Ls = Float64[]
        for t in range(59000.0 - 4000, 59000.0 + 4000, length=41)
            s = rawsolve(sys, t)
            x, y, z = posx(s, :b, :A), posy(s, :b, :A), posz(s, :b, :A)
            vx, vy, vz = velx(s, :b, :A), vely(s, :b, :A), velz(s, :b, :A)
            r = hypot(x, y, z); v2 = vx^2 + vy^2 + vz^2
            @test v2 ≈ μ * (2 / r - 1 / a) rtol = 1e-12      # vis-viva
            push!(Es, v2 / 2 - μ / r)
            push!(Ls, hypot(y * vz - z * vy, z * vx - x * vz, x * vy - y * vx))
        end
        @test (maximum(Es) - minimum(Es)) / abs(Es[1]) < 1e-12
        @test (maximum(Ls) - minimum(Ls)) / Ls[1] < 1e-12
    end

    # a > 0 with e > 1 has no valid reading; it is taken as |a|
    sp = System((Body(mass=1.0, name=:A), Body(mass=0.0, name=:b)),
        (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A); a=5.0, e=1.5),))
    @test semimajoraxis(sp) == -5.0

    # Hot path stays allocation-free, and AD agrees with finite differences
    A = Body(mass=1.1, name=:A); b = Body(mass=0.0, name=:b)
    hsys = System((A, b), (Orbit(b, about=A; a=-5.0, e=1.4, i=0.5, ω=1.1, Ω=2.2, tp=59000.0),); plx=25.0)
    ep = collect(range(58000.0, 60000.0, length=25))
    htraj = Trajectory(hsys, ep)
    orbitsolve!(htraj, hsys)
    @test (@allocated orbitsolve!(htraj, hsys)) == 0
    hf = θ -> begin
        s = System((Body(mass=θ[1], name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=θ[1], name=:A);
                   a=θ[2], e=θ[3], i=θ[4], ω=1.1, Ω=2.2, tp=59000.0 + θ[5]),); plx=25.0)
        sum(raoff(x, :b, :A) + radvel(x, :b, :A) for x in orbitsolve(s, ep))
    end
    θh = [1.1, -5.0, 1.4, 0.5, 0.0]
    gad = ForwardDiff.gradient(hf, θh)
    gfd = FiniteDiff.finite_difference_gradient(hf, θh)
    for j in eachindex(θh)
        @test isapprox(gad[j], gfd[j]; rtol=1e-5, atol=1e-6 * max(1.0, abs(gad[j])))
    end

    # Hyperbolic rows compose with the hierarchy: an unbound companion in a
    # 3-body system still satisfies the A⁻¹ definition.
    c = Body(mass=5mjup, name=:c)
    mix = System((A, c, Body(mass=1mjup, name=:d)), (
        Orbit(c, about=A; a=3.0, e=0.2, i=0.4, ω=1.0, Ω=0.5, tp=59000.0),
        Orbit(Body(mass=1mjup, name=:d), about=(A, c); a=-20.0, e=1.6, i=0.7, ω=0.3, Ω=1.2, tp=59000.0)))
    @test ainv_residual(mix) < 1e-13
    @test period(mix, 2) == Inf
    @test all(isfinite, orbitsolve(mix, [58500.0, 59000.0, 59500.0]).x)
end

@testset "parametrization groups" begin
    mk(; kw...) = System(
        (Body(mass=1.1, name=:A), Body(mass=0.0, name=:b)),
        (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.1, name=:A); kw...),); plx=25.0)

    @testset "Cartesian round-trip" begin
        # elements → state at an epoch → elements must reconstruct the orbit.
        # This is what pins the perifocal basis convention; a textbook sign
        # error in the orbit normal shows up here and nowhere else.
        for (nm, kw) in (
            ("circular-ish", (; a=5.0, e=0.01, i=0.5, ω=1.1, Ω=2.2, tp=59000.0)),
            ("eccentric", (; a=5.0, e=0.7, i=0.5, ω=1.1, Ω=2.2, tp=59000.0)),
            ("very eccentric", (; a=5.0, e=0.95, i=0.5, ω=1.1, Ω=2.2, tp=59000.0)),
            ("i > π/2", (; a=5.0, e=0.3, i=2.6, ω=1.1, Ω=2.2, tp=59000.0)),
            ("near edge-on", (; a=5.0, e=0.3, i=π / 2 - 1e-6, ω=0.4, Ω=1.0, tp=59000.0)),
            ("hyperbolic", (; a=-5.0, e=1.4, i=0.5, ω=1.1, Ω=2.2, tp=59000.0)),
            ("hyperbolic wide", (; a=-50.0, e=1.05, i=0.9, ω=2.9, Ω=0.3, tp=59000.0)))
            @testset "$nm" begin
                s0 = mk(; kw...)
                ep = 59123.0
                # Raw states: Cartesian *initial conditions* are dynamical
                # state at a coordinate time, not an observation, so the
                # round-trip must not go through the observing pass.
                sol = rawsolve(s0, ep)
                s1 = mk(; x=posx(sol, :b, :A), y=posy(sol, :b, :A), z=posz(sol, :b, :A),
                    vx=velx(sol, :b, :A), vy=vely(sol, :b, :A), vz=velz(sol, :b, :A),
                    epoch=ep)
                r0, r1 = s0.rows[1], s1.rows[1]
                @test r1.a ≈ r0.a rtol = 1e-12
                @test r1.e ≈ r0.e atol = 1e-12
                @test abs(rem(r1.i - r0.i, 2π, RoundNearest)) < 1e-12
                @test abs(rem(r1.ω - r0.ω, 2π, RoundNearest)) < 1e-11
                @test abs(rem(r1.Ω - r0.Ω, 2π, RoundNearest)) < 1e-12
                # the actual requirement: identical trajectories
                ts = collect(range(58500.0, 59600.0, length=17))
                t0 = rawsolve(s0, ts); t1 = rawsolve(s1, ts)
                @test maximum(abs, t0.x .- t1.x) < 1e-10
                @test maximum(abs, t0.vx .- t1.vx) < 1e-10
            end
        end
    end

    @testset "shape groups" begin
        for (e, ω) in ((0.0, 0.0), (0.3, 1.1), (0.85, 4.0), (0.5, -2.0))
            base = mk(; a=5.0, e=e, ω=ω, i=0.5, Ω=2.2, tp=59000.0)
            se = mk(; a=5.0, secosω=√e * cos(ω), sesinω=√e * sin(ω), i=0.5, Ω=2.2, tp=59000.0)
            ee = mk(; a=5.0, ecosω=e * cos(ω), esinω=e * sin(ω), i=0.5, Ω=2.2, tp=59000.0)
            ts = [58800.0, 59000.0, 59400.0]
            @test maximum(abs, orbitsolve(base, ts).x .- orbitsolve(se, ts).x) < 1e-14
            @test maximum(abs, orbitsolve(base, ts).x .- orbitsolve(ee, ts).x) < 1e-14
        end
    end

    @testset "phase groups" begin
        for (e, ω, i, Ω) in ((0.0, 0.0, 0.5, 2.2), (0.4, 1.1, 0.5, 2.2), (0.8, 3.0, 1.2, 0.4))
            base = mk(; a=5.0, e=e, ω=ω, i=i, Ω=Ω, tp=59000.0)
            ep = 59211.0
            n_per_day = base.rows[1].n / PlanetOrbits.year2day_julian
            viaM0 = mk(; a=5.0, e=e, ω=ω, i=i, Ω=Ω, M0=n_per_day * (ep - 59000.0), epoch=ep)
            # θ is the sky-plane position angle at `epoch`; recovering tp from it
            # needs no mass and no `a` (the radius factor cancels)
            # `θ` is the *orbital* sky-plane position angle at `epoch` — a
            # dynamical quantity, so it is read off the raw state rather than
            # the observed (retarded) one, and from the physical-unit
            # components rather than the angular ones.
            rs = rawsolve(base, ep)
            θraw = atan(posx(rs, :b, :A), posy(rs, :b, :A))
            viaθ = mk(; a=5.0, e=e, ω=ω, i=i, Ω=Ω, θ=θraw, epoch=ep)
            @test viaM0.rows[1].tp ≈ base.rows[1].tp atol = 1e-8
            @test viaθ.rows[1].tp ≈ base.rows[1].tp atol = 1e-8
            ts = [58800.0, 59000.0, 59400.0]
            @test maximum(abs, orbitsolve(base, ts).x .- orbitsolve(viaM0, ts).x) < 1e-12
            @test maximum(abs, orbitsolve(base, ts).x .- orbitsolve(viaθ, ts).x) < 1e-12
        end
    end

    @testset "group validation" begin
        A = Body(mass=1.1, name=:A); b = Body(mass=0.0, name=:b)
        @test_throws "at most one eccentricity" Orbit(b, about=A; a=1.0, e=0.1, secosω=0.2, sesinω=0.1)
        @test_throws "must be given together" Orbit(b, about=A; a=1.0, secosω=0.2)
        @test_throws "at most one phase" Orbit(b, about=A; a=1.0, tp=5.0, M0=1.0, epoch=59000.0)
        @test_throws "pass `epoch=`" Orbit(b, about=A; a=1.0, M0=1.0)
        @test_throws "pass `epoch=`" Orbit(b, about=A; a=1.0, θ=1.0)
        @test_throws "all six of" Orbit(b, about=A; x=1.0, y=2.0)
        @test_throws "determine every orbital element" Orbit(b, about=A; a=1.0,
            x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2, vz=0.3, epoch=59000.0)
        @test_throws "need `epoch=`" Orbit(b, about=A; x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2, vz=0.3)
        @test_throws "radial" Orbit(b, about=A; x=1.0, y=0.0, z=0.0, vx=0.5, vy=0.0, vz=0.0, epoch=59000.0)
    end

    @testset "Cartesian hot path" begin
        buildc(θ) = System((Body(mass=θ[1], name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=θ[1], name=:A);
                   x=θ[2], y=θ[3], z=θ[4], vx=θ[5], vy=θ[6], vz=θ[7], epoch=59000.0),); plx=25.0)
        sθ = SVector(1.1, 3.0, 1.0, 0.5, -1.5, 2.0, 0.3)
        buildc(sθ)
        @test (@allocated buildc(sθ)) == 0
        ep = collect(range(58800.0, 59200.0, length=15))
        f = θ -> sum(raoff(s, :b, :A) + radvel(s, :b, :A) for s in orbitsolve(buildc(θ), ep))
        θc = [1.1, 3.0, 1.0, 0.5, -1.5, 2.0, 0.3]
        gad = ForwardDiff.gradient(f, θc)
        gfd = FiniteDiff.finite_difference_gradient(f, θc)
        for j in eachindex(θc)
            @test isapprox(gad[j], gfd[j]; rtol=1e-5, atol=1e-6 * max(1.0, abs(gad[j])))
        end
    end
end

@testset "Thiele-Innes, RV-only, and export policy" begin
    # Thiele-Innes round-trips against the Campbell elements, including the
    # regions where v1's ThieleInnesOrbit was documented as wrong: Ω ≥ π and
    # ω + Ω > 2π.
    for (a, i, ω, Ω) in ((5.0, 0.5, 1.1, 2.2), (5.0, 0.5, 1.1, 4.0),
                         (5.0, 0.5, 5.0, 5.5), (2.0, 2.6, 3.0, 3.5),
                         (1.0, 0.01, 0.2, 0.3), (8.0, π / 2 - 1e-8, 6.0, 0.1))
        ti = thieleinnes(System(
            (Body(mass=1.1, name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.1, name=:A);
                   a=a, e=0.3, i=i, ω=ω, Ω=Ω, tp=59000.0),)))
        got = ThieleInnes(; ti...)
        @test got.a ≈ a rtol = 1e-12
        @test got.i ≈ i atol = 1e-10
        # Recovery is exact only up to the ±180° node ambiguity: (ω, Ω) and
        # (ω+π, Ω+π) give identical constants, so accept either branch. The
        # returned one always has Ω ∈ [0, π).
        @test 0 <= got.Ω < π + 1e-12
        dω = abs(rem(got.ω - ω, 2π, RoundNearest))
        dΩ = abs(rem(got.Ω - Ω, 2π, RoundNearest))
        @test (dω < 1e-10 && dΩ < 1e-10) ||
              (abs(dω - π) < 1e-10 && abs(dΩ - π) < 1e-10)
        # …and reconstructs the same sky-plane trajectory through the
        # constructor either way (the two nodes differ only along the line of
        # sight).
        #
        # Compared on raw states: the two branches have opposite-sign z, so
        # per-body light-travel retardation moves them to different points on
        # the orbit and the sky-plane tracks separate at v/c (~1e-4 AU here).
        # The node degeneracy is therefore only *almost* exact once light
        # travel time is modelled — far below any real measurement, but not
        # zero.
        base = System((Body(mass=1.1, name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.1, name=:A);
                   a=a, e=0.3, i=i, ω=ω, Ω=Ω, tp=59000.0),); plx=25.0)
        viaTI = System((Body(mass=1.1, name=:A), Body(mass=0.0, name=:b)),
            (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.1, name=:A);
                   got..., e=0.3, tp=59000.0),); plx=25.0)
        ts = [58800.0, 59000.0, 59400.0]
        @test maximum(abs, rawsolve(base, ts).x .- rawsolve(viaTI, ts).x) < 1e-10
        @test maximum(abs, orbitsolve(base, ts).x .- orbitsolve(viaTI, ts).x) < 1e-3
    end
    # mas ↔ AU via plx
    timas = thieleinnes(System(
        (Body(mass=1.0, name=:A), Body(mass=0.0, name=:b)),
        (Orbit(Body(mass=0.0, name=:b), about=Body(mass=1.0, name=:A);
               a=4.0, e=0.1, i=0.6, ω=1.0, Ω=2.0, tp=59000.0),)); plx=25.0)
    @test ThieleInnes(; timas..., plx=25.0).a ≈ 4.0 rtol = 1e-12

    # RV-only convention: no parallax, i = π/2, Ω = 0
    rv = rvorbit(M=1.1, msini=2mjup, P=400.0, e=0.2, ω=1.0, tp=59000.0)
    @test inclination(rv) ≈ π / 2
    @test periastron(rv) == 59000.0
    @test period(rv) ≈ 400.0 rtol = 1e-12
    @test rv.masses[2] ≈ 2mjup
    sol = orbitsolve(rv, 59100.0)
    @test isfinite(radvel(sol, :A, barycentre(rv)))
    @test_throws ErrorException raoff(sol)      # no parallax ⇒ no angular obs

    # §5 export policy: Octofitter owns these three names unqualified
    for nm in (:System, :Body, :Orbit)
        @test !Base.isexported(PlanetOrbits, nm)
        @test isdefined(PlanetOrbits, nm)
    end
    for nm in (:bodies, :barycentre, :photocentre, :fluxes, :Jacobi, :Astrocentric,
               :ThieleInnes, :thieleinnes, :orbitsolve, :period,
               :BodyRef, :WeightedPoint)
        @test Base.isexported(PlanetOrbits, nm)
    end
    @test !Base.isexported(PlanetOrbits, :orbit)
    @test !Base.isexported(PlanetOrbits, :rvorbit)
end

include("observing-geometry.jl")
include("reframe.jl")
include("nbody.jl")
include("row-cache.jl")
include("enzyme-rules.jl")
include("plotting.jl")
