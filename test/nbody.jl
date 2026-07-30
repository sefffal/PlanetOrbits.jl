# ---------------------------------------------------
# AHL21 propagator gates (§11 "N-body sanity" of the design doc).
# Included from runtests.jl (uses _ac_benign and the θ0 workload defined there).
# ---------------------------------------------------

using PlanetOrbits: NBodyState, ahl21_step, G_au3_day2, _initial_state, _reverse_v

include("fixtures/nbody_reference.jl")

# Total energy and angular momentum of an integrator state [M⊙, AU, days].
function _nbody_energy(st::NBodyState{N}, m) where {N}
    E = 0.0
    for i in 1:N
        E += 0.5 * m[i] * (st.v[1,i]^2 + st.v[2,i]^2 + st.v[3,i]^2)
        for j in i+1:N
            r = sqrt((st.x[1,i]-st.x[1,j])^2 + (st.x[2,i]-st.x[2,j])^2 + (st.x[3,i]-st.x[3,j])^2)
            E -= G_au3_day2 * m[i] * m[j] / r
        end
    end
    return E
end
function _nbody_angmom(st::NBodyState{N}, m) where {N}
    L = zeros(3)
    for i in 1:N
        L[1] += m[i] * (st.x[2,i]*st.v[3,i] - st.x[3,i]*st.v[2,i])
        L[2] += m[i] * (st.x[3,i]*st.v[1,i] - st.x[1,i]*st.v[3,i])
        L[3] += m[i] * (st.x[1,i]*st.v[2,i] - st.x[2,i]*st.v[1,i])
    end
    return L
end

@testset "AHL21 kernel equivalence vs NbodyGradient v0.2.3" begin
    # Same Cartesian ICs, masses, G, and step counts as upstream — the merged
    # kernels must reproduce upstream trajectories to near machine precision.
    # (Measured deviation ≤ 8e-15 after 1000 steps; gate at 20× margin.)
    for case in NBODY_CASES
        N = length(case.m)
        Gm = SVector{N}(NBODY_GNEWT .* case.m)
        st = NBodyState(SMatrix{3,N}(case.x0), SMatrix{3,N}(case.v0))
        nstep = 0
        @testset "$(case.name)" begin
            for (si, target) in enumerate(NBODY_SNAP_STEPS)
                while nstep < target
                    st = ahl21_step(st, Gm, case.h)
                    nstep += 1
                end
                @test maximum(abs, Matrix(st.x) .- case.snaps[si].x) < 2e-13
                @test maximum(abs, Matrix(st.v) .- case.snaps[si].v) < 2e-13
            end
        end
    end
end

@testset "AHL21 two-body exactness vs KeplerianApprox" begin
    # For one pair the corrector vanishes and the map is exact two-body
    # evolution: the propagators must agree to roundoff at ANY h, including
    # epochs before t0 (backward march) and off-grid epochs (partial steps).
    A = Body(mass=1.1, name=:A)
    b = Body(mass=5mjup, name=:b)
    sys = System(Orbit(b, about=A; a=3.0, e=0.45, i=0.6, ω=1.1, Ω=2.2, tp=58849.0); plx=24.5)
    P = PlanetOrbits.period(sys, 1)
    epochs = collect(range(58849.0 - 2.3P, 58849.0 + 3.7P, length=41))
    tk = orbitsolve(sys, epochs)
    ta = orbitsolve(sys, epochs; method=AHL21(h=P/23, t0=58849.0))
    for k in eachindex(epochs), f in (raoff, decoff, radvel)
        @test isapprox(f(tk[k], :b, :A), f(ta[k], :b, :A); rtol=1e-9, atol=1e-8)
    end

    # Absolute-frame system: t0 defaults to ref_epoch; frame composition
    # (LTT, space motion) is shared, so observables agree end-to-end.
    θ = θ0
    Aa = Body(mass=θ[1], name=:A); ba = Body(mass=θ[2], name=:b)
    sysa = System(Orbit(ba, about=Aa; a=θ[3], e=θ[4], i=θ[5], ω=θ[6], Ω=θ[7], tp=θ[8]);
        plx=θ[9], ra=θ[10], dec=θ[11], pmra=θ[12], pmdec=θ[13], rv=θ[14], ref_epoch=θ[15])
    Pa = PlanetOrbits.period(sysa, 1)
    epochs2 = collect(range(58000.0, 60000.0, length=17))
    tka = orbitsolve(sysa, epochs2)
    taa = orbitsolve(sysa, epochs2; method=AHL21(h=Pa/40))   # t0 = ref_epoch
    w = barycentre(sysa)
    for k in eachindex(epochs2)
        @test isapprox(raoff(tka[k], :b, :A), raoff(taa[k], :b, :A); rtol=1e-9, atol=1e-8)
        @test isapprox(radvel(tka[k], :A, w), radvel(taa[k], :A, w); rtol=1e-9, atol=1e-8)
        @test frame_pmra(tka[k]) == frame_pmra(taa[k])
    end
end

@testset "AHL21 energy & angular momentum conservation" begin
    # Kepler-36-like three-body over ~50 inner periods at h ≈ P/20.
    A = Body(mass=1.071, name=:A)
    b = Body(mass=1.32e-5, name=:b)
    c = Body(mass=2.42e-5, name=:c)
    sys = System(
        Orbit(c, about=Orbit(b, about=A; a=0.1153, e=0.022, i=π/2, ω=0.4, Ω=0.0, tp=60000.0);
              a=0.1283, e=0.016, i=π/2 + 0.01, ω=1.3, Ω=0.05, tp=60002.0); plx=50.0)
    st = _initial_state(Float64, sys, 60000.0)
    m = sys.masses
    Gm = SVector{3}(G_au3_day2 .* m)
    E0 = _nbody_energy(st, m)
    L0 = _nbody_angmom(st, m)
    h = 0.7
    # Symplectic signature: the energy error is a bounded h² oscillation
    # (measured 4.8e-8 at h ≈ P/20 for this system, reached within the first
    # periods and never exceeded), NOT a secular drift.
    worst_first = 0.0
    worst_all = 0.0
    for s in 1:1000
        st = ahl21_step(st, Gm, h)
        rel = abs(_nbody_energy(st, m) - E0) / abs(E0)
        s <= 500 && (worst_first = max(worst_first, rel))
        worst_all = max(worst_all, rel)
    end
    L1 = _nbody_angmom(st, m)
    @test worst_all < 1e-6                          # bounded oscillation
    @test worst_all <= 1.5 * worst_first + 1e-14    # no secular growth
    @test maximum(abs, L1 .- L0) / sqrt(sum(abs2, L0)) < 1e-12
end

@testset "AHL21 zero-mass test particles" begin
    # A massless body feels the massive pair but does not disturb it: the
    # massive pair must remain exactly two-body (matching KeplerianApprox),
    # and the test particle's states must be finite.
    A = Body(mass=1.0, name=:A)
    b = Body(mass=10mjup, name=:b)
    tp = Body(mass=0.0, name=:tp)
    inner = Orbit(b, about=A; a=1.0, e=0.2, i=0.3, ω=0.5, Ω=1.0, tp=59000.0)
    sys = System(Orbit(tp, about=inner; a=6.0, e=0.05, i=0.35, ω=2.0, Ω=1.1, tp=58000.0); plx=30.0)
    P = PlanetOrbits.period(sys, 1)
    epochs = collect(range(59000.0, 59000.0 + 3P, length=21))
    tk = orbitsolve(sys, epochs)
    ta = orbitsolve(sys, epochs; method=AHL21(h=P/25, t0=59000.0))
    for k in eachindex(epochs)
        @test isapprox(raoff(tk[k], :b, :A), raoff(ta[k], :b, :A); rtol=1e-8, atol=1e-7)
        @test isfinite(raoff(ta[k], :tp, :A))
        @test isfinite(radvel(ta[k], :tp, :A))
    end
    # Two test particles (gm == 0 pair) must not NaN either.
    tp2 = Body(mass=0.0, name=:tp2)
    sys2 = System(Orbit(tp2, about=Orbit(tp, about=A; a=2.0, tp=59000.0); a=5.0, tp=59000.0); plx=30.0)
    ta2 = orbitsolve(sys2, epochs; method=AHL21(h=10.0, t0=59000.0))
    @test all(isfinite, ta2.x)
end

@testset "AHL21 t0 handling & error paths" begin
    A = Body(mass=1.0, name=:A); b = Body(mass=1mjup, name=:b)
    sys = System(Orbit(b, about=A; a=1.0, tp=59000.0); plx=30.0)
    @test_throws "osculating epoch" orbitsolve(sys, [59000.0]; method=AHL21(h=5.0))
    @test_throws "must be positive" AHL21(h=-1.0)
    @test_throws "must be positive" AHL21(h=0.0)
    # explicit t0 on a parallax-only system works
    sol = orbitsolve(sys, [59100.0]; method=AHL21(h=5.0, t0=59000.0))[1]
    @test isfinite(raoff(sol, :b, :A))
end

# AHL21 analog of the KeplerianApprox micro workload: 3-body system, mixed
# forward/backward epochs, full observable set. Used for gradient +
# allocation gates.
function _eval_workload_nbody(θ, epochs, traj=nothing)
    A = Body(mass=θ[1], name=:A)
    b = Body(mass=θ[2], name=:b)
    c = Body(mass=θ[3], name=:c)
    # tp parameters are offsets from a base epoch: finite-difference reference
    # gradients need steps ≪ the 14-day inner period, which a raw ~59000 MJD
    # parameter's relative step would violate.
    sys = System(
        Orbit(c, about=Orbit(b, about=A; a=θ[4], e=θ[5], i=θ[6], ω=θ[7], Ω=θ[8], tp=59000.0 + θ[9]);
              a=θ[10], e=θ[11], i=θ[12], ω=θ[13], Ω=θ[14], tp=59000.0 + θ[15]);
        plx=θ[16], ra=45.0, dec=-30.0, pmra=100.0, pmdec=-50.0, rv=25e3, ref_epoch=59005.0)
    if traj === nothing
        traj = Trajectory{eltype(θ)}(sys, epochs)
    end
    orbitsolve!(traj, sys; method=AHL21(h=0.65, t0=59005.0))
    refs = bodies(sys)
    bary = barycentre(sys)
    acc = zero(eltype(θ))
    for sol in traj
        acc += raoff(sol, refs.b, refs.A) + decoff(sol, refs.c, refs.A) +
               pmra(sol, refs.A, bary) + radvel(sol, refs.A, bary) + frame_rv(sol)
    end
    return acc
end

const θ_nb = [1.071, 1.32e-5, 2.42e-5,
              0.1153, 0.022, π/2, 0.4, 0.0, 0.0,
              0.1283, 0.016, π/2 + 0.01, 1.3, 0.05, 2.0,
              50.0]
const epochs_nb = collect(range(58980.0, 59060.0, length=12))  # straddles t0

@testset "AHL21 ForwardDiff gradient" begin
    f(θ) = _eval_workload_nbody(θ, epochs_nb)
    g_fd = ForwardDiff.gradient(f, θ_nb)
    g_ref = FiniteDiff.finite_difference_gradient(f, θ_nb)
    for i in eachindex(θ_nb)
        @test isapprox(g_fd[i], g_ref[i]; rtol=1e-4, atol=1e-4 * max(1.0, abs(g_fd[i])))
    end
    @test f(θ_nb) isa Float64
end

@testset "AHL21 type stability & allocation-free hot path" begin
    θ = SVector{16}(θ_nb)
    A = Body(mass=θ[1], name=:A); b = Body(mass=θ[2], name=:b); c = Body(mass=θ[3], name=:c)
    sys = System(
        Orbit(c, about=Orbit(b, about=A; a=θ[4], e=θ[5], i=θ[6], ω=θ[7], Ω=θ[8], tp=θ[9]);
              a=θ[10], e=θ[11], i=θ[12], ω=θ[13], Ω=θ[14], tp=θ[15]);
        plx=θ[16], ra=45.0, dec=-30.0, pmra=100.0, pmdec=-50.0, rv=25e3, ref_epoch=59005.0)
    traj = @inferred orbitsolve(sys, epochs_nb; method=AHL21(h=0.65, t0=59005.0))
    @inferred raoff(traj[1], :c, :A)
    # Allocation gates: forced bounds checking (`Pkg.test` passes
    # --check-bounds=yes) inserts throw branches that capture the MMatrix
    # state buffers, forcing them to heap — the gates are only meaningful
    # when bounds checks are not forced (same caveat as the SIMD IR check;
    # exercised via `include("runtests.jl")` and the perf/ harness, where
    # the AHL21 hot path is verified to be 0-alloc).
    if Base.JLOptions().check_bounds == 1
        @info "skipping AHL21 allocation gates: bounds checking is forced on"
    else
        trajbuf = Trajectory(sys, epochs_nb)
        _eval_workload_nbody(θ, epochs_nb, trajbuf)  # warm up
        allocs = @allocated _eval_workload_nbody(θ, epochs_nb, trajbuf)
        @test allocs == 0
    end
end

_ac_solve_ahl21!(traj, sys) = orbitsolve!(traj, sys; method=AHL21(h=0.65, t0=59005.0))

@testset "AHL21 static allocation-freedom (AllocCheck)" begin
    if Base.JLOptions().check_bounds == 1
        @info "skipping AHL21 AllocCheck gate: bounds checking is forced on"
    else
        θ = SVector{16}(θ_nb)
        A = Body(mass=θ[1], name=:A); b = Body(mass=θ[2], name=:b); c = Body(mass=θ[3], name=:c)
        sys = System(
            Orbit(c, about=Orbit(b, about=A; a=θ[4], e=θ[5], i=θ[6], ω=θ[7], Ω=θ[8], tp=θ[9]);
                  a=θ[10], e=θ[11], i=θ[12], ω=θ[13], Ω=θ[14], tp=θ[15]);
            plx=θ[16], ra=45.0, dec=-30.0, pmra=100.0, pmdec=-50.0, rv=25e3, ref_epoch=59005.0)
        traj = Trajectory(sys, epochs_nb)
        errs = filter(!_ac_benign, AllocCheck.check_allocs(_ac_solve_ahl21!, (typeof(traj), typeof(sys))))
        isempty(errs) || display(errs[1])
        @test isempty(errs)
    end
end
