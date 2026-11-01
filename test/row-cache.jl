# ---------------------------------------------------
# Per-row solve caching (see the block comment in keplerian-approx.jl) and
# the branch-free `vlog` kernel.
#
# The cache's contract is that a hit is *bit-exact*: it restores, from its
# own snapshot, values the identical computation already produced. So every
# test here compares against a fresh, uncached solve for equality, not
# approximate agreement. The cache is opt-in per call (`row_cache=`), so the
# reference solves simply don't pass one.
# ---------------------------------------------------

@testset "vlog kernel" begin
    # <1 ulp against Base.log across the σ² ranges a likelihood feeds it
    # (mas² scales through wildly mis-scaled jitters), plus exact powers of
    # two and values straddling the [√2/2, √2) reduction boundary.
    worst = 0.0
    for x in vcat(exp10.(range(-12, 8, length=4001)),
                  [2.0^k for k in -40:40],
                  nextfloat.(fill(sqrt(2) / 2, 3), 0:2))
        worst = max(worst, abs(PlanetOrbits.vlog(x) - log(x)) / max(abs(log(x)), 1e-300))
    end
    @test worst ≤ 1e-15
    # generic fallback: Duals go through Base.log, so derivatives are exact
    d = ForwardDiff.derivative(PlanetOrbits.vlog, 3.7)
    @test d == 1 / 3.7
end

@testset "row solve cache" begin
    # Heap trajectories in this testset: Matrix/Vector columns.
    cache = PlanetOrbits._rowcache(Matrix{Float64}, Vector{Float64})
    epochs = collect(range(57000.0, 60000.0, length=93))

    buildsys(; mb=5mjup, a=2.5, framekw...) = System(
        (Body(mass=1.1, name=:A), Body(mass=mb, name=:b), Body(mass=2mjup, name=:c)),
        (Orbit(Body(mass=mb, name=:b), about=Body(mass=1.1, name=:A);
               a=a, e=0.3, i=0.5, ω=1.1, Ω=2.2, tp=58849.0),
         Orbit(Body(mass=2mjup, name=:c), about=Body(mass=1.1, name=:A);
               a=8.0, e=0.1, i=0.6, ω=0.4, Ω=2.0, tp=57000.0));
        framekw...)

    # Reference values: plain solves, no cache involved.
    sys0 = buildsys(plx=24.5)
    ref = orbitsolve(sys0, epochs)

    # One trajectory buffer re-solved repeatedly, as the sampler hot loop does.
    traj = Trajectory(sys0, epochs)

    # First solve into fresh storage: all rows miss.
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj, sys0; row_cache=cache)
    @test cache.misses - m0 == 2 && cache.hits - h0 == 0
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traj, f) == getfield(ref, f)
    end

    # Identical re-solve: all rows hit, bit-identical states.
    h0 = cache.hits
    orbitsolve!(traj, sys0; row_cache=cache)
    @test cache.hits - h0 == 2
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traj, f) == getfield(ref, f)
    end

    # The slice-sampler pattern: change ONE row's element; only that row
    # re-solves, and the result equals a fresh uncached solve.
    sys1 = buildsys(plx=24.5, a=2.6)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj, sys1; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (1, 1)
    ref1 = orbitsolve(sys1, epochs)
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traj, f) == getfield(ref1, f)
    end

    # A companion-mass change touches that row's M (and the A⁻¹ combine),
    # but row c is untouched: one hit, one miss.
    sys2 = buildsys(plx=24.5, a=2.6, mb=6mjup)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj, sys2; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (1, 1)

    # plx is not a row input on a Parallax frame: full hit, and the observe
    # pass still sees the new parallax.
    sys3 = buildsys(plx=30.0, a=2.6, mb=6mjup)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj, sys3; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (2, 0)
    ref3 = orbitsolve(sys3, epochs)
    @test traj.cart2angle == ref3.cart2angle

    # AbsoluteFrame: t_em depends on the frame, so a frame change misses; a
    # lighttime-flag flip (with rv ≠ 0) misses through the t_em samples.
    fkw = (; plx=24.5, ra=45.0, dec=-30.0, pmra=100.0, pmdec=-50.0,
           rv=25e3, ref_epoch=57388.5)
    sysa = buildsys(; a=2.6, mb=6mjup, fkw...)
    traja = Trajectory(sysa, epochs)
    orbitsolve!(traja, sysa; row_cache=cache)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traja, sysa; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (2, 0)
    sysb = buildsys(; a=2.6, mb=6mjup, fkw..., rv=26e3)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traja, sysb; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (0, 2)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traja, sysb; row_cache=cache, barycentric_lighttime=false)
    @test cache.misses - m0 == 2
    refb = orbitsolve(sysb, epochs; barycentric_lighttime=false)
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traja, f) == getfield(refb, f)
    end

    # The method is part of the key: flipping simd= on the same storage and
    # elements must not serve the other kernel's bits.
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traja, sysb; row_cache=cache, barycentric_lighttime=false,
                method=KeplerianApprox(simd=false))
    @test cache.misses - m0 == 2

    # A different trajectory at identical elements and epochs hits: the
    # states come from the cache's snapshot, not from the storage, so they
    # are correct in storage that has never been solved into. (Seed with the
    # default method and flags first: the solve above used other ones.)
    orbitsolve!(traja, sysb; row_cache=cache)
    traj2 = Trajectory(sysb, epochs)
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj2, sysb; row_cache=cache)
    @test (cache.hits - h0, cache.misses - m0) == (2, 0)
    refb2 = orbitsolve(sysb, epochs)
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traj2, f) == getfield(refb2, f)
    end

    # Storage scribbled between two identical solves: a hit must not serve
    # what is sitting in the columns.
    for f in (:rx, :ry, :rz, :rvx, :rvy, :rvz, :x, :y, :z, :vx, :vy, :vz)
        fill!(getfield(traj2, f), NaN)
    end
    h0 = cache.hits
    orbitsolve!(traj2, sysb; row_cache=cache)
    @test cache.hits - h0 == 2
    for f in (:x, :y, :z, :vx, :vy, :vz)
        @test getfield(traj2, f) == getfield(refb2, f)
    end

    # A Dual-valued solve handed the (Float64-typed) cache takes the plain
    # loop via the type guard: correct values, counters untouched.
    h0, m0 = cache.hits, cache.misses
    g = ForwardDiff.derivative(2.6) do a
        s = buildsys(plx=24.5, a=a, mb=6mjup)
        t = Trajectory{typeof(a)}(s, epochs)
        orbitsolve!(t, s; row_cache=cache)
        t.x[45, 2]
    end
    @test isfinite(g)
    @test (cache.hits, cache.misses) == (h0, m0)

    # No row_cache argument: strict re-fill semantics, counters untouched.
    h0, m0 = cache.hits, cache.misses
    orbitsolve!(traj2, sysb)
    @test (cache.hits, cache.misses) == (h0, m0)

    # Warm steady state allocates nothing (runtime check; the static
    # AllocCheck gate certifies the default, cache-free path).
    solve_warm(t, s, c) = orbitsolve!(t, s; row_cache=c)
    solve_warm(traj2, sysb, cache)
    if DYNAMIC_ALLOC_GATE
        @test @allocated(solve_warm(traj2, sysb, cache)) == 0
    else
        @test_skip @allocated(solve_warm(traj2, sysb, cache)) == 0
    end
end

# Octofitter #162. Octofitter carves each evaluation's trajectory out of one
# task-local bump arena, so successive trajectories have the *same* column
# objects (an `UnsafeArray` is egal on pointer+dims) — but the bytes behind
# them do not survive between evaluations: a gradient evaluation, another
# model, or any other `@alloc` user of the arena writes over them. The cache
# used to read a pointer match as "my previous solve's states are still
# there" and skip the solve, returning leftovers as orbital states — a
# silently wrong ℓπ, and a NaN only when the leftovers happened to be NaN.
# A view-backed arena reproduces that storage model without Bumper.
@testset "row solve cache: arena-backed storage (Octofitter #162)" begin
    epochs = collect(range(56000.0, 58500.0, length=6))
    sys = System(
        (Body(mass=1.0, name=:A), Body(mass=1e-3, name=:b)),
        (Orbit(Body(mass=1e-3, name=:b), about=Body(mass=1.0, name=:A);
               a=3.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=55000.0),);
        plx=25.0)
    ref = orbitsolve(sys, epochs)

    arena = zeros(100_000)
    offset = Ref(0)
    function arenalloc(::Type{S}, dims::Integer...) where {S}
        S === Float64 || return Array{S}(undef, dims...)
        n = prod(dims)
        v = view(arena, offset[]+1:offset[]+n)
        offset[] += n
        return length(dims) == 1 ? v : reshape(v, dims)
    end
    # One "evaluation": reset the bump pointer, carve, solve.
    function evaluate(s)
        offset[] = 0
        traj = Trajectory(arenalloc, Float64, s, epochs)
        orbitsolve!(traj, s; row_cache=PlanetOrbits._rowcache(traj))
        return traj
    end

    t1 = evaluate(sys)
    cache = PlanetOrbits._rowcache(t1)
    for poison in (NaN, 0.0, 1e3, -7.5)
        fill!(arena, poison)   # someone else used the arena in between
        h0 = cache.hits
        t = evaluate(sys)
        @test t.rx === t1.rx   # same storage object: the #162 precondition
        @test cache.hits - h0 == 1
        for f in (:x, :y, :z, :vx, :vy, :vz, :cart2angle)
            @test getfield(t, f) == getfield(ref, f)
        end
    end
end

@testset "_POISON_TRAJECTORIES debug switch" begin
    epochs = collect(range(57000.0, 58000.0, length=5))
    sys = System(Orbit(Body(mass=1e-3, name=:b), about=Body(mass=1.0, name=:A);
                       a=3.0, e=0.1, i=0.2, ω=0.3, Ω=0.4, tp=57000.0); plx=25.0)
    old = PlanetOrbits._POISON_TRAJECTORIES[]
    try
        PlanetOrbits._POISON_TRAJECTORIES[] = true
        t = Trajectory(sys, epochs)
        @test all(isnan, t.rx) && all(isnan, t.x) && all(isnan, t.R11)
        @test t.epochs == epochs
        orbitsolve!(t, sys)
        @test t.x == orbitsolve(sys, epochs).x
        # Dual columns are poisoned too.
        td = Trajectory{ForwardDiff.Dual{Nothing,Float64,1}}(sys, epochs)
        @test all(x -> isnan(ForwardDiff.value(x)), td.x)
    finally
        PlanetOrbits._POISON_TRAJECTORIES[] = old
    end
end
