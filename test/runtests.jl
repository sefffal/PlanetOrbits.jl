using PlanetOrbits
using PlanetOrbits: orbit   # v1-compat two-body sugar; deliberately unexported
using Test
import ForwardDiff
import FiniteDiff
using StaticArrays
using InteractiveUtils: code_llvm
import AllocCheck

include("fixtures/v1_reference.jl")

approx(a, b) = isapprox(a, b; rtol=1e-11, atol=1e-10)

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
            sys, refs = fixture_system(c)
            @test approx(period(sys), c.period)
            traj = orbitsolve(sys, c.epochs)
            d = c.data
            for (k, sol) in enumerate(traj)
                @test approx(posx(sol), d.posx[k])
                @test approx(posy(sol), d.posy[k])
                @test approx(posz(sol), d.posz[k])
                @test approx(velx(sol), d.velx[k])
                @test approx(vely(sol), d.vely[k])
                @test approx(velz(sol), d.velz[k])
            end
            if c.kind == :kep
                for (k, sol) in enumerate(traj)
                    @test approx(radvel(sol), d.radvel[k])
                end
            end
            if c.kind in (:visual, :absvis)
                for (k, sol) in enumerate(traj)
                    @test approx(raoff(sol), d.raoff[k])
                    @test approx(decoff(sol), d.decoff[k])
                    @test approx(projectedseparation(sol), d.projectedseparation[k])
                    @test approx(posangle(sol), d.posangle[k])
                end
            end
            if c.kind == :visual
                for (k, sol) in enumerate(traj)
                    @test approx(radvel(sol), d.radvel[k])
                    @test approx(pmra(sol), d.pmra[k])
                    @test approx(pmdec(sol), d.pmdec[k])
                end
            end
            if c.kind == :absvis
                # v1 adds the propagated-frame drift to pm and rv (but not to
                # the position offsets); in v2 that composition is explicit.
                p = c.params
                for (k, sol) in enumerate(traj)
                    @test approx(radvel(sol) + (frame_rv(sol) - p.rv), d.radvel[k])
                    @test approx(pmra(sol) + (frame_pmra(sol) - p.pmra), d.pmra[k])
                    @test approx(pmdec(sol) + (frame_pmdec(sol) - p.pmdec), d.pmdec[k])
                    @test approx(frame_ra(sol), d.comp_ra2[k])
                    @test approx(frame_dec(sol), d.comp_dec2[k])
                    @test approx(frame_pmra(sol), d.comp_pmra2[k])
                    @test approx(frame_pmdec(sol), d.comp_pmdec2[k])
                    @test approx(frame_rv(sol), d.comp_rv2[k])
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
                    @test approx(raoff(sol, refsr.b, refsr.A), d.raoff[k])
                    @test approx(raoff(sol, refsr.A, bary), d.raoff_reflex[k])
                    @test approx(decoff(sol, refsr.A, bary), d.decoff_reflex[k])
                    if c.kind == :visual
                        @test approx(pmra(sol, refsr.A, bary), d.pmra_reflex[k])
                        @test approx(pmdec(sol, refsr.A, bary), d.pmdec_reflex[k])
                        @test approx(radvel(sol, refsr.A, bary), d.radvel_reflex[k])
                    else
                        @test approx(pmra(sol, refsr.A, bary) + (frame_pmra(sol) - p.pmra), d.pmra_reflex[k])
                        @test approx(pmdec(sol, refsr.A, bary) + (frame_pmdec(sol) - p.pmdec), d.pmdec_reflex[k])
                        @test approx(radvel(sol, refsr.A, bary) + (frame_rv(sol) - p.rv), d.radvel_reflex[k])
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
    # planet's separation from the inner barycentre solves the outer row
    @test approx(
        hypot(posx(sol, crefs.b, innerbary), posy(sol, crefs.b, innerbary), posz(sol, crefs.b, innerbary)),
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
    @test posx(sol) isa Float64
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

function _eval_workload(θ, epochs, traj=nothing)
    A = Body(mass=θ[1], name=:A)
    b = Body(mass=θ[2], name=:b)
    sys = System(Orbit(b, about=A; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
        plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
    if traj === nothing
        traj = Trajectory{eltype(θ)}(sys, epochs)
    end
    orbitsolve!(traj, sys)
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
_ac_solve!(traj, sys) = orbitsolve!(traj, sys; method=KeplerianApprox(solver=PlanetOrbits.Markley()))
_ac_solve_scalar!(traj, sys) =
    orbitsolve!(traj, sys; method=KeplerianApprox(solver=PlanetOrbits.Markley(), simd=false))
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
    # e = 0, so the moon's 3D separation from its host is exactly a
    @test hypot(posx(sol, :m, :b), posy(sol, :m, :b), posz(sol, :m, :b)) ≈ 0.007 rtol = 1e-12
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
    qs = orbitsolve(q, 58900.0)
    bA = barycentre(q, :Aa, :Ab); bB = barycentre(q, :Ba, :Bb)
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
        μ = 4π^2 * Mtot                       # AU³ / kepler-yr²
        Es = Float64[]; Ls = Float64[]
        for t in range(59000.0 - 4000, 59000.0 + 4000, length=41)
            s = orbitsolve(sys, t)
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

include("nbody.jl")
