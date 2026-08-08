# Backend-free plotting support: the resolver table, epoch grids, and phase.
# The Makie extension itself is exercised by docs/logo.jl and Octofitter's
# plotting tests; nothing here needs a plotting package.

using PlanetOrbits: orbit

@testset "plotinfo / plotlabel" begin
    @test plotlabel(radvel) == "radial velocity [m/s]"
    @test plotlabel(raoff) == "Δα* [mas]"
    @test plotinfo(raoff).flip
    @test !plotinfo(decoff).flip
    @test plotinfo(posangle).wrap == 2π
    @test plotinfo(projectedseparation).wrap === nothing
end

@testset "orbit_track_epochs" begin
    o = orbit(a=1.2, e=0.7, i=0.4, ω=1.0, Ω=2.0, tp=59000.0, plx=100.0, M=1.0)
    P = period(o)
    ts = orbit_track_epochs(o; n=100)
    @test length(ts) == 100
    @test ts[1] == 59000.0                       # anchored at periastron
    @test ts[end] - ts[1] ≈ P rtol = 1e-12       # closed track
    # EA-uniform: dense at periastron, sparse at apastron
    dts = diff(ts)
    @test dts[1] < dts[49] && dts[end] < dts[49]
    # the shifted cycle contains tstart
    ts2 = orbit_track_epochs(o; n=100, tstart=59000.0 + 7.3P)
    @test ts2[1] ≈ 59000.0 + 7P rtol = 1e-12

    # circular orbit: uniform everywhere
    oc = orbit(a=1.0, e=0.0, tp=0.0, M=1.0)
    dtc = diff(orbit_track_epochs(oc; n=50))
    @test all(x -> isapprox(x, dtc[1], rtol=1e-8), dtc)

    # hyperbolic: refused with a pointer
    oh = orbit(a=-2.0, e=1.4, tp=59000.0, M=1.0)
    @test_throws ErrorException orbit_track_epochs(oh)
end

@testset "plot_epochs" begin
    o = orbit(a=1.0, e=0.1, tp=59000.0, plx=50.0, M=1.0)
    ts = plot_epochs(o, 59000.0, 61000.0)
    @test issorted(ts) && allunique(ts)
    @test 200 <= length(ts) <= 1000
    @test ts[1] == 59000.0 && ts[end] == 61000.0
    @test_throws ErrorException plot_epochs(o, 100.0, 100.0)

    # a short-period inner orbit densifies the grid; the cap holds
    A = PlanetOrbits.Body(mass=1.0, name=:A)
    b = PlanetOrbits.Body(mass=0.001, name=:b)
    c = PlanetOrbits.Body(mass=0.001, name=:c)
    sys = PlanetOrbits.System((A, b, c), (
        PlanetOrbits.Orbit(b, about=A; P=10.0, e=0.1, tp=0.0),
        PlanetOrbits.Orbit(c, about=(A, b); P=10000.0, e=0.1, tp=0.0),
    ))
    ts2 = plot_epochs(sys, 59000.0, 60000.0)
    @test length(ts2) == 1000                      # clamped
end

@testset "orbit_phase" begin
    o = orbit(a=1.0, e=0.3, tp=59000.0, M=1.0)
    P = period(o)
    @test orbit_phase(o, 59000.0) == 0.0
    @test orbit_phase(o, 59000.0 + P / 4) ≈ π / 2 atol = 1e-8
    @test orbit_phase(o, 59000.0 - P / 4) ≈ 3π / 2 atol = 1e-8   # wraps
    oh = orbit(a=-2.0, e=1.4, tp=59000.0, M=1.0)
    @test isnan(orbit_phase(oh, 59000.0))
end

@testset "trajectory broadcasting" begin
    o = orbit(a=1.0, e=0.3, tp=59000.0, plx=100.0, M=1.0)
    traj = orbitsolve(o, [59000.0, 59100.0, 59200.0])
    xs = raoff.(traj)
    @test xs isa Vector && length(xs) == 3
    @test xs == [raoff(traj[k]) for k in 1:3]
    ys = radvel.(traj, :b, :A)
    @test length(ys) == 3
end

@testset "orbit() phase-group passthrough" begin
    # M0+epoch through the two-body sugar must match the equivalent tp.
    o1 = orbit(a=1.0, e=0.2, M0=π / 2, epoch=59000.0, plx=100.0, M=1.0)
    P = period(o1)
    tp_expect = 59000.0 - (π / 2) / (2π) * P
    o2 = orbit(a=1.0, e=0.2, tp=tp_expect, plx=100.0, M=1.0)
    s1 = orbitsolve(o1, 59123.4)
    s2 = orbitsolve(o2, 59123.4)
    @test raoff(s1) ≈ raoff(s2) rtol = 1e-12
    @test decoff(s1) ≈ decoff(s2) rtol = 1e-12
end
