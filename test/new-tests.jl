using Test
using DirectOrbits
using ForwardDiff

@testset begin
    
    # Create an idealized orbit like the Earth's at 1pc distance.
    circular_face_on_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.0,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )

    # Basics: check locations.
    os = DirectOrbits.orbitsolve_ν(circular_face_on_1AU_1Msun_1pc, 0)
    # All zero angles should point at celestial north pole
    @test decoff(os) ≈ 1000
    @test raoff(os) ≈ 0
    # Planet is at maximum vertical separation so instantenous PM should be zero
    @test os.ẏ ≈ 0
    # Travelling CCW (CW in plane of the sky)
    @test sign(os.ẋ) == -1

end


@testset "Derivatives" begin
    
    # Create an idealized orbit like the Earth's at 1pc distance.
    circular_face_on_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0,
        i = 0.0,
        e = 0.0,
        τ = 0.0,
        M = 1.0,
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas <-> 1pc
    )

    # Basics: check locations.
    os = DirectOrbits.orbitsolve(circular_face_on_1AU_1Msun_1pc, 100.0)
    os_d = ForwardDiff.derivative(
        t->DirectOrbits.orbitsolve(circular_face_on_1AU_1Msun_1pc, t),
        100.0
    )
    @test os.ẋ ≈ os_d.x
    @test os.ẋ ≈ os_d.x

end
