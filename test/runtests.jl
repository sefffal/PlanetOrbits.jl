using PlanetOrbits
using Test
import ForwardDiff
import FiniteDiff
using StaticArrays
using InteractiveUtils: code_llvm

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
    sys = System(Binary(A, b; a=p.a, e=p.e, i=p.i, ω=p.ω, Ω=p.Ω, tp=p.tp); framekw...)
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
    inner = Binary(A, b; a=2.5, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
    sys = System(Binary(inner, c; a=8.0, e=0.3, i=0.6, ω=0.4, Ω=2.0, tp=57000.0); plx=50.0)
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
    tight = Binary(B1, B2; a=0.2, e=0.05, i=0.3, ω=0.2, Ω=0.1, tp=58849.0)
    csys = System(Binary(tight, p; a=3.0, e=0.1, i=0.4, ω=1.0, Ω=2.0, tp=58000.0); plx=80.0)
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
    zsys = System(Binary(Body(mass=1.0, name=:A), Body(mass=0.0, name=:b);
        a=5.0, e=0.2, i=0.4, ω=1.0, Ω=0.5, tp=58849.0); plx=40.0)
    zrefs = bodies(zsys)
    ztraj = orbitsolve(zsys, [58900.0])
    @test abs(raoff(ztraj[1], zrefs.A, barycentre(zsys))) < 1e-13
end

@testset "photocentre" begin
    A = Body(mass=1.0, flux=(G=1.0,), name=:A)
    b = Body(mass=0.2, flux=(G=0.0,), name=:b)
    sys = System(Binary(A, b; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    refs = bodies(sys)
    sol = orbitsolve(sys, [59000.0])[1]
    # dark companion: photocentre is the star
    pc = photocentre(sys)
    @test raoff(sol, pc, barycentre(sys)) ≈ raoff(sol, refs.A, barycentre(sys))
    # equal-brightness pair: photocentre at the midpoint
    A2 = Body(mass=1.0, flux=(G=1.0,), name=:A)
    b2 = Body(mass=0.2, flux=(G=1.0,), name=:b)
    sys2 = System(Binary(A2, b2; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    refs2 = bodies(sys2)
    sol2 = orbitsolve(sys2, [59000.0])[1]
    pc2 = photocentre(sys2)
    mid = (raoff(sol2, refs2.A, barycentre(sys2)) + raoff(sol2, refs2.b, barycentre(sys2))) / 2
    @test raoff(sol2, pc2, barycentre(sys2)) ≈ mid
    # multi-band requires selection
    A3 = Body(mass=1.0, flux=(G=1.0, K=0.5), name=:A)
    b3 = Body(mass=0.2, flux=(K=0.4,), name=:b)
    sys3 = System(Binary(A3, b3; a=4.0, e=0.1, i=0.5, ω=1.0, Ω=2.0, tp=58849.0); plx=30.0)
    @test_throws ErrorException photocentre(sys3)
    @test photocentre(sys3; band=:G).w[2] == 0
    @test photocentre(sys3; band=:K).w[2] ≈ 0.4 / 0.9
end

@testset "error paths" begin
    A = Body(mass=1.0, name=:A)
    b = Body(mass=0.001, name=:b)
    @test_throws ErrorException System(Binary(A, b; a=1.0, e=1.5, i=0.0, ω=0.0, Ω=0.0, tp=0.0))
    # angular observables need a parallax
    sys = System(Binary(A, b; a=1.0, e=0.1, i=0.2, ω=0.0, Ω=0.0, tp=58849.0))
    sol = orbitsolve(sys, 58900.0)
    @test_throws ErrorException raoff(sol)
    @test posx(sol) isa Float64
    # partial absolute frame
    @test_throws ErrorException System(Binary(A, b; a=1.0, e=0.1, i=0.2, ω=0.0, Ω=0.0, tp=0.0);
        plx=10.0, ra=45.0)
    # duplicate names
    @test_throws ErrorException System(Binary(Body(mass=1.0, name=:x), Body(mass=0.1, name=:x);
        a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0))
    # one-argument observables on >2 bodies
    m = Body(mass=0.0001, name=:m)
    sys3 = System(Binary(Binary(A, b; a=1.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0), m;
        a=5.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=0.0); plx=10.0)
    sol3 = orbitsolve(sys3, 58900.0)
    @test_throws ErrorException raoff(sol3)
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
    inner = Binary(A, b; a=2.5, e=0.6, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
    sys = System(Binary(inner, c; a=8.0, e=0.1, i=0.6, ω=0.4, Ω=2.0, tp=57000.0);
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
    sys = System(Binary(A, b; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
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
    sys = System(Binary(A, b; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
        plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
    # preallocated trajectory stands in for Bumper-owned buffers
    traj = Trajectory(sys, epochs)
    _eval_workload(θ, epochs, traj)  # warm up
    allocs = @allocated _eval_workload(θ, epochs, traj)
    @test allocs == 0
end

@testset "type stability" begin
    epochs = [58900.0, 59000.0]
    A = Body(mass=1.1, name=:A); b = Body(mass=5mjup, name=:b)
    root = Binary(A, b; a=8.0, e=0.1, i=0.5, ω=1.1, Ω=2.2, tp=58849.0)
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
