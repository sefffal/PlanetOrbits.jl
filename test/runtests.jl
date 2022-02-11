using Test
using DirectOrbits
using ForwardDiff
using FiniteDiff

# 10 steps per day for one year
one_year_range = 0.0:0.1:365.24

# Relative tolerance for certain tests
rtol=1e-6
# Absolute tolerance for certain tests
atol=1e-6


@testset "Positions" begin
    
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
    @test pmdec(os) ≈ 0
    # Travelling CCW (CW in plane of the sky)
    @test sign(pmra(os)) == +1

end



# Begin with basic consistency tests
@testset "Circular" begin
    
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

    # Period of the orbit is 1 year by definition
    @test period(circular_face_on_1AU_1Msun_1pc) == 1.0*DirectOrbits.year2day

    # 1000mas==1as of paralax means 1pc, by definition
    @test distance(circular_face_on_1AU_1Msun_1pc) == 1

    # Mean motion is one revolution per year
    @test DirectOrbits.meanmotion(circular_face_on_1AU_1Msun_1pc) == 2π
    
    # Face on orbits should have their x & y positions add to 1AU separation
    x = raoff.(circular_face_on_1AU_1Msun_1pc, one_year_range)
    y = decoff.(circular_face_on_1AU_1Msun_1pc, one_year_range)
    sep = sqrt.(x.^2 .+ y.^2)
    @test all(≈(1000.0), sep)
    # And should have no radial velocity
    @test all(radvel.(circular_face_on_1AU_1Msun_1pc, one_year_range) .≈ 0)
    
    # We should return to our initial position after exactly one period
    @test orbitsolve(circular_face_on_1AU_1Msun_1pc, period(circular_face_on_1AU_1Msun_1pc)) ≈ orbitsolve(circular_face_on_1AU_1Msun_1pc, 0.0)

    # We should be halfway around the orbit after half a period (in this case)
    @test -orbitsolve(circular_face_on_1AU_1Msun_1pc, period(circular_face_on_1AU_1Msun_1pc)/2) ≈ orbitsolve(circular_face_on_1AU_1Msun_1pc, 0.0)

end


@testset "Inclined" begin
    # Now add some inclination
    circular_inclined_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0, # AU
        i = deg2rad(90.0),
        e = 0.0,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    @test period(circular_inclined_1AU_1Msun_1pc) == 1.0*DirectOrbits.year2day
    @test distance(circular_inclined_1AU_1Msun_1pc) == 1
    @test DirectOrbits.meanmotion(circular_inclined_1AU_1Msun_1pc) == 2π
    @test -orbitsolve(circular_inclined_1AU_1Msun_1pc, period(circular_inclined_1AU_1Msun_1pc)/2) ≈ orbitsolve(circular_inclined_1AU_1Msun_1pc, 0.0)

    ys = decoff.(circular_inclined_1AU_1Msun_1pc, one_year_range)
    @test maximum(ys) ≈ 1000 rtol=rtol
    @test minimum(ys) ≈ -1000 rtol=rtol

    xs = raoff.(circular_inclined_1AU_1Msun_1pc, one_year_range)
    @test all(≈(0.0, atol=atol), xs)

    # Now adjust Ω
    inclined_rot_90deg = KeplerianElements(
        a = 1.0, # AU
        i = deg2rad(90.0),
        e = 0.0,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = deg2rad(90.0),
        plx = 1000.0, # 1000 mas == 1pc
    )
    
    xs = raoff.(inclined_rot_90deg, one_year_range)
    ys = decoff.(inclined_rot_90deg, one_year_range)
    @test all(≈(0.0, atol=atol), ys)
    @test maximum(xs) ≈ 1000 rtol=rtol
    @test minimum(xs) ≈ -1000 rtol=rtol

    # Intermediate inclination
    inclined_45deg_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0, # AU
        i = deg2rad(45.0),
        e = 0.0,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    
    ys = decoff.(inclined_45deg_1AU_1Msun_1pc, one_year_range)
    @test maximum(ys) ≈ 1000 rtol=rtol
    @test minimum(ys) ≈ -1000 rtol=rtol

    xs = raoff.(inclined_45deg_1AU_1Msun_1pc, one_year_range)
    @test maximum(xs) ≈ 1000*√2/2 rtol=rtol
    @test minimum(xs) ≈ -1000*√2/2 rtol=rtol

end



@testset "Eccentricity" begin
    # Basic eccentric orbit
    eccentric_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.5,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(eccentric_1AU_1Msun_1pc, one_year_range)
    ys = decoff.(eccentric_1AU_1Msun_1pc, one_year_range)
    ps = projectedseparation.(eccentric_1AU_1Msun_1pc, one_year_range)

    @test period(eccentric_1AU_1Msun_1pc) == 1.0*DirectOrbits.year2day
    @test distance(eccentric_1AU_1Msun_1pc) == 1
    
    # Mean motion should be the same
    @test DirectOrbits.meanmotion(eccentric_1AU_1Msun_1pc) == 2π

    # The separation should now be varying
    # By definition of eccentricity 0.5, 1AU and 1PC
    @test maximum(ps) ≈ 1500 rtol=rtol
    @test minimum(ps) ≈ 500 rtol=rtol

    # When argument of periapsis and periastron are both zero, periastron should be in the East, apoastron in the West
    @test maximum(ys) ≈ 500 rtol=rtol
    @test minimum(ys) ≈ -1500 rtol=rtol


    # Rotate Ω
    ecc_rot_Ω = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.5,
        τ = 90.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = deg2rad(90),
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_Ω, one_year_range)
    ys = decoff.(ecc_rot_Ω, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test minimum(xs) ≈ -1500 rtol=rtol
    @test maximum(xs) ≈ 500 rtol=rtol

    # Rotate τ
    ecc_rot_ω = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.5,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = deg2rad(90.0),
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_ω, one_year_range)
    ys = decoff.(ecc_rot_ω, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test minimum(xs) ≈ -1500 rtol=rtol
    @test maximum(xs) ≈ 500 rtol=rtol

    # Rotate Ω & τ
    ecc_rot_Ωτ = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.5,
        τ = 90.0,
        M = 1.0, # M_sun
        ω = deg2rad(-90),
        Ω = deg2rad(90),
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_Ωτ, one_year_range)
    ys = decoff.(ecc_rot_Ωτ, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test maximum(ys) ≈ 500 rtol=rtol
    @test minimum(ys) ≈ -1500 rtol=rtol

    
    # Highly eccentric 
    ecc09 = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 0.9,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc09, one_year_range)
    ys = decoff.(ecc09, one_year_range)
    ps = projectedseparation.(ecc09, one_year_range)
    # Loosen the tolerance on these
    @test maximum(ps) ≈ 1900 rtol=1e-4
    @test minimum(ps) ≈ 100 rtol=1e-4

    # Extremely eccentric 
    ecc09 = KeplerianElements(
        a = 1.0, # AU
        i = 0.0,
        e = 1-1e-3,
        τ = 0.0,
        M = 1.0, # M_sun
        ω = 0.0,
        Ω = 0.0,
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc09, one_year_range)
    ys = decoff.(ecc09, one_year_range)
    ps = projectedseparation.(ecc09, one_year_range)
    @test maximum(ps) ≈ 1999 rtol=1e-4
    # Loosen the tolerance on these even more (periastron flies by very quickly)
    @test minimum(ps) ≈ 1 rtol=1e1

    
end;  

##

@testset "Chain rules" begin
    # These tests are broken at MA===0, e>0

    # First test analytic chain rules
    k1(MA) =e->DirectOrbits.kepler_solver(MA, e)
    k2(e) = MA->DirectOrbits.kepler_solver(MA, e)
    
    for e in 0:0.1:0.9
        for MA in 0.001:0.1:2π
            @test FiniteDiff.finite_difference_derivative(k2(e), MA) ≈ ForwardDiff.derivative(k2(e), MA) rtol=rtol
        end
    end

    for e = 0.001:0.1:0.9
        for MA in 0.001:0.1:2π
            @test FiniteDiff.finite_difference_derivative(k1(MA), e) ≈ ForwardDiff.derivative(k1(MA), e) rtol=rtol
        end
    end
    
end

##
@testset "PMA & Accel." begin


    # Check analytic derivative properties against ForwardDiff over a big range of orbits
    for t in 0.:35:356.,
        a in 0.1:0.2:3,
        e in 0:0.1:0.9,
        i in deg2rad.([-45, 0, 45, 90, ]),
        ω in deg2rad.([-45, 0, 45, 90, ]),
        Ω in deg2rad.([-45, 0, 45, 90, ])

        elems = KeplerianElements(;
            a,
            i = 0.0,
            e,
            τ = 0.0,
            M = 1.0,
            ω = 0.0,
            Ω = 0.0,
            plx = 1000.0, # 1000 mas <-> 1pc
        )

        @test pmra(elems, 100.0) ≈ ForwardDiff.derivative(
            t->raoff(elems, t),
            100.0
        )*DirectOrbits.year2day

        @test pmdec(elems, 100.0) ≈ ForwardDiff.derivative(
            t->decoff(elems, t),
            100.0
        )*DirectOrbits.year2day

        @test acceleration(elems, 100.0)[1] ≈ ForwardDiff.derivative(
            t->pmra(elems, t),
            100.0
        )*DirectOrbits.year2day

        @test acceleration(elems, 100.0)[2] ≈ ForwardDiff.derivative(
            t->pmdec(elems, t),
            100.0
        )*DirectOrbits.year2day

        
    end

end




# @testset "Motion" begin
#     circ = KeplerianElements(
#         a = 1.0,
#         i = 0.0,
#         e = 0.0,
#         τ = 0.0,
#         M = 1.0,
#         ω = 0.0,
#         Ω = 0.0,
#         plx = 1000.0,
#     )

#     @test kep2cart(circ, 0.0, tref=0.)[1:2] ≈ [0., 1000.0] rtol=rtol
#     @test kep2cart(circ, period(circ)/4, tref=0.)[1:2] ≈ [1000., 0.0] rtol=rtol
#     @test kep2cart(circ, period(circ)/2, tref=0.)[1:2] ≈ [0.0, -1000.0] rtol=rtol
#     @test kep2cart(circ, period(circ)*3/4, tref=0.)[1:2] ≈ [-1000.0, 0.0] rtol=rtol
#     @test kep2cart(circ, period(circ), tref=0.)[1:2] ≈ [0., 1000.0]  rtol=rtol

#     ecc_rot_ω = KeplerianElements(
#         a = 1.0, # AU
#         i = 0.0,
#         e = 0.5,
#         τ = 0.0,
#         M = 1.0, # M_sun
#         ω = deg2rad(90.0),
#         Ω = 0.0,
#         plx = 1000.0, # 1000 mas == 1pc
#     )

#     @test kep2cart(ecc_rot_ω, 0.0, tref=0.)[1:2] ≈ [500.0, 0.0] rtol=rtol
#     @test kep2cart(ecc_rot_ω, period(ecc_rot_ω)/2, tref=0.)[1:2] ≈ [-1500., 0.0] rtol=rtol


#     circt2 = KeplerianElements(
#         a = 1.0,
#         i = 0.0,
#         e = 0.0,
#         τ = 0.5,
#         M = 1.0,
#         ω = 0.0,
#         Ω = 0.0,
#         plx = 1000.0,
#     )

#     @test kep2cart(circt2, 0.0, tref=0.)[1:2] ≈ [0., -1000.0] rtol=rtol
#     @test kep2cart(circt2, period(circt2)/4, tref=0.)[1:2] ≈ [-1000., 0.0] rtol=rtol
#     @test kep2cart(circt2, period(circt2)/2, tref=0.)[1:2] ≈ [0.0, 1000.0] rtol=rtol
#     @test kep2cart(circt2, period(circt2)*3/4, tref=0.)[1:2] ≈ [1000.0, 0.0] rtol=rtol
#     @test kep2cart(circt2, period(circt2), tref=0.)[1:2] ≈ [0., -1000.0]  rtol=rtol
   
# end


# Next step is integrating examples from the literature

# ##

# @testset "Transformations" begin
    
#     # Test that our orbital transformations code produces the same results as our forward direct code, tested above.

#     # for i in 0:0.1:2π
#     # for i=0.0, a=0.5:1.0:10, e=0.0, ω=0.00, Ω=0:1:2π, plx=1000.0, M=1.0, τ=0.8, t₀=0:10:300
#     for i=0.0, a=1.0, e=0.0, ω=0.01, Ω=0.01, plx=1000.0, M=1.0, τ=0.8, t₀=0.
#         el = KeplerianElements(;a, i, e, ω, Ω, plx, M, τ)
#         pos1 = kep2cart(el, t₀)
#         # for dt in range(0, stop=period(el), length=4)
#         for dt=10.0
#             pos2 = kep2cart(el, t₀+dt)[1:2]
#             ot = OrbitalTransformation(;i, e, ω, Ω, plx, M, platescale=1., dt=dt)
#             pos2′ = ot(pos1)
#             # @show pos2 pos2′
#             @test pos2′ ≈ pos2 rtol=1e-2
#         end
#     end
# end