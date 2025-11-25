using Test
using PlanetOrbits

@testset "Vectorized orbitsolve API" begin
    # Create a test orbit
    testorbit = orbit(
        a = 1.0,
        e = 0.5,
        i = 45.0,
        ω = 30.0,
        Ω = 60.0,
        tp = 0.0,
        M = 1.0,
        plx = 100.0
    )

    # Test epochs
    epochs = [0.0, 100.0, 200.0, 365.24, 500.0]

    @testset "Generic fallback (Auto solver)" begin
        # Compute using generic fallback
        solutions_vec = orbitsolve(testorbit, epochs)

        # Compute individually
        solutions_individual = [orbitsolve(testorbit, t) for t in epochs]

        # Compare results
        @test length(solutions_vec) == length(epochs)
        for i in eachindex(solutions_vec)
            @test solutions_vec[i].t ≈ solutions_individual[i].t
            @test solutions_vec[i].ν ≈ solutions_individual[i].ν rtol=1e-10
            @test solutions_vec[i].EA ≈ solutions_individual[i].EA rtol=1e-10
            @test posx(solutions_vec[i]) ≈ posx(solutions_individual[i]) rtol=1e-10
            @test posy(solutions_vec[i]) ≈ posy(solutions_individual[i]) rtol=1e-10
            @test posz(solutions_vec[i]) ≈ posz(solutions_individual[i]) rtol=1e-10
        end
    end

    @testset "Markley vectorized solver" begin
        # Compute using vectorized Markley solver
        solutions_vec = orbitsolve(testorbit, epochs, Markley())

        # Compute individually with Markley
        solutions_individual = [orbitsolve(testorbit, t, Markley()) for t in epochs]

        # Compare results
        @test length(solutions_vec) == length(epochs)
        for i in eachindex(solutions_vec)
            @test solutions_vec[i].t ≈ solutions_individual[i].t
            @test solutions_vec[i].ν ≈ solutions_individual[i].ν rtol=1e-10
            @test solutions_vec[i].EA ≈ solutions_individual[i].EA rtol=1e-10
            @test posx(solutions_vec[i]) ≈ posx(solutions_individual[i]) rtol=1e-10
            @test posy(solutions_vec[i]) ≈ posy(solutions_individual[i]) rtol=1e-10
            @test posz(solutions_vec[i]) ≈ posz(solutions_individual[i]) rtol=1e-10
        end
    end

    @testset "Consistency between Auto and Markley" begin
        # For elliptical orbits, Auto should use Markley
        solutions_auto = orbitsolve(testorbit, epochs)
        solutions_markley = orbitsolve(testorbit, epochs, Markley())

        @test length(solutions_auto) == length(solutions_markley)
        for i in eachindex(solutions_auto)
            @test solutions_auto[i].t ≈ solutions_markley[i].t
            @test solutions_auto[i].ν ≈ solutions_markley[i].ν rtol=1e-10
            @test solutions_auto[i].EA ≈ solutions_markley[i].EA rtol=1e-10
        end
    end

    @testset "Different orbit types" begin
        # Test with RadialVelocityOrbit
        rvorbit = RadialVelocityOrbit(
            a = 1.0,
            e = 0.3,
            ω = 45.0,
            tp = 0.0,
            M = 1.0
        )

        solutions_vec = orbitsolve(rvorbit, epochs, Markley())
        solutions_individual = [orbitsolve(rvorbit, t, Markley()) for t in epochs]

        @test length(solutions_vec) == length(epochs)
        for i in eachindex(solutions_vec)
            @test solutions_vec[i].ν ≈ solutions_individual[i].ν rtol=1e-10
            @test radvel(solutions_vec[i]) ≈ radvel(solutions_individual[i]) rtol=1e-10
        end
    end

    @testset "Large number of epochs" begin
        # Test with many epochs to verify vectorization works correctly
        many_epochs = range(0, 1000, length=1000)

        solutions_vec = orbitsolve(testorbit, collect(many_epochs), Markley())

        @test length(solutions_vec) == length(many_epochs)
        # Spot check a few solutions
        for i in [1, 250, 500, 750, 1000]
            sol_individual = orbitsolve(testorbit, many_epochs[i], Markley())
            @test solutions_vec[i].ν ≈ sol_individual.ν rtol=1e-10
            @test solutions_vec[i].EA ≈ sol_individual.EA rtol=1e-10
        end
    end
end
