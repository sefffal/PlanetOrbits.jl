# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

using Test
using PlanetOrbits
using ForwardDiff
using FiniteDiff

# ----------------------------------------------------------------------------------------------------------------------
# Constants and Helper Functions
# ----------------------------------------------------------------------------------------------------------------------

# 10 steps per day for one year
one_year_range = 0.0:0.1:365.24
# Relative tolerance for certain tests
rtol = 1e-6
# Absolute tolerance for certain tests
atol = 1e-6

# ----------------------------------------------------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------------------------------------------------

## Test relationships between inverse constants
@testset "Constants" begin
    @test PlanetOrbits.mas2rad == 1/PlanetOrbits.rad2mas
    @test PlanetOrbits.as2rad == 1/PlanetOrbits.rad2as
    @test PlanetOrbits.au2pc == 1/PlanetOrbits.pc2au
    @test PlanetOrbits.m2au == 1/PlanetOrbits.au2m
    @test PlanetOrbits.day2year == 1/PlanetOrbits.year2day
    @test PlanetOrbits.sec2year == 1/PlanetOrbits.year2sec
    @test PlanetOrbits.sec2day == 1/PlanetOrbits.day2sec
end



## Idealized face-on Earth with circular orbit at 1 pc 
@testset "Earth, i = 0, e = 0, d = 1 pc" begin
    idealearth = orbit(
        a = 1.0,
        e = 0.0,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        tp = 0.0,
        M = 1.0,
        plx = 1000.0
    )

    # Test basic orbit properties
    @test period(idealearth) ≈ PlanetOrbits.year2day
    @test distance(idealearth) ≈ 1.0
    @test meanmotion(idealearth) ≈ 2π
    @test periastron(idealearth) ≈ 58849
    @test semiamplitude(idealearth) ≈ 0.0

    # Orbit solutions at quarters of the orbit
    oq1 = PlanetOrbits.orbitsolve_ν(idealearth, 0.0)
    oq2 = PlanetOrbits.orbitsolve_ν(idealearth, π/2)
    oq3 = PlanetOrbits.orbitsolve_ν(idealearth, π)
    oq4 = PlanetOrbits.orbitsolve_ν(idealearth, 3π/2)

    # Test orbit values at first quarter
    @test raoff(oq1) ≈ 0.0 atol=atol
    @test decoff(oq1) ≈ 1000.0 rtol=rtol
    @test posangle(oq1) ≈ 0.0 atol=atol
    @test projectedseparation(oq1) ≈ 1000.0 rtol=rtol
    
    @test sign(pmra(oq1)) == +1
    @test pmdec(oq1) ≈ 0.0 atol=atol
    @test radvel(oq1) ≈ 0.0 atol=atol

    @test accra(oq1) ≈ 0.0 atol=atol
    @test sign(accdec(oq1)) == -1

    # Test orbit values at second quarter
    @test raoff(oq2) ≈ 1000.0 rtol=rtol
    @test decoff(oq2) ≈ 0.0 atol=atol
    @test posangle(oq2) ≈ π/2 rtol=rtol
    @test projectedseparation(oq2) ≈ 1000.0 rtol=rtol

    @test pmra(oq2) ≈ 0.0 atol=atol
    @test sign(pmdec(oq2)) == -1
    @test radvel(oq2) ≈ 0.0 atol=atol

    @test sign(accra(oq2)) == -1
    @test accdec(oq2) ≈ 0.0 atol=atol

    # Test orbit values at third quarter
    @test raoff(oq3) ≈ 0.0 atol=atol
    @test decoff(oq3) ≈ -1000.0 rtol=rtol
    @test posangle(oq3) ≈ π rtol=rtol
    @test projectedseparation(oq3) ≈ 1000.0 rtol=rtol

    @test sign(pmra(oq3)) == -1
    @test pmdec(oq3) ≈ 0.0 atol=atol
    @test radvel(oq3) ≈ 0.0 atol=atol

    @test accra(oq3) ≈ 0.0 atol=atol 
    @test sign(accdec(oq3)) == +1

    # Test orbit values at fourth quarter
    @test raoff(oq4) ≈ -1000.0 rtol=rtol
    @test decoff(oq4) ≈ 0.0 atol=atol
    @test posangle(oq4) ≈ -π/2 rtol=rtol
    @test projectedseparation(oq4) ≈ 1000.0 rtol=rtol
    
    @test pmra(oq4) ≈ 0.0 atol=atol
    @test sign(pmdec(oq4)) == +1
    @test radvel(oq4) ≈ 0.0 atol=atol

    @test sign(accra(oq4)) == +1
    @test accdec(oq4) ≈ 0.0 atol=atol

    # Compare velocities and accelerations
    @test pmra(oq1) ≈ -pmra(oq3) rtol=rtol
    @test pmdec(oq2) ≈ -pmdec(oq4) rtol=rtol
    @test accdec(oq1) ≈ -accdec(oq3) rtol=rtol
    @test accra(oq2) ≈ -accra(oq4) rtol=rtol
end

## Idealized edge-on Earth with circular orbit at 1 pc 
@testset "Earth, i = 90, e = 0, d = 1 pc" begin
    idealearth = orbit(
        a = 1.0,
        e = 0.0,
        i = π/2,
        ω = 0.0,
        Ω = 0.0,
        tp = 0.0,
        M = 1.0,
        plx = 1000.0
    )

    # Test basic orbit properties
    @test period(idealearth) == PlanetOrbits.year2day
    @test distance(idealearth) == 1.0
    @test meanmotion(idealearth) == 2π
    @test periastron(idealearth) == 58849
    @test semiamplitude(idealearth) ≈ 29785.89 rtol=1e-3

    # Orbit solutions at quarters of the orbit
    oq1 = PlanetOrbits.orbitsolve_ν(idealearth, 0.0)
    oq2 = PlanetOrbits.orbitsolve_ν(idealearth, π/2)
    oq3 = PlanetOrbits.orbitsolve_ν(idealearth, π)
    oq4 = PlanetOrbits.orbitsolve_ν(idealearth, 3π/2)

    # Test orbit values at first quarter
    @test raoff(oq1) ≈ 0.0 atol=atol
    @test decoff(oq1) ≈ 1000.0 rtol=rtol
    @test projectedseparation(oq1) ≈ 1000.0 rtol=rtol
    
    @test pmra(oq1) ≈ 0.0 atol=atol
    @test pmdec(oq1) ≈ 0.0 atol=atol
    @test radvel(oq1) ≈ 29785.89 rtol=1e-3

    @test accra(oq1) ≈ 0.0 atol=atol
    @test sign(accdec(oq1)) == -1

    # Test orbit values at second quarter
    @test raoff(oq2) ≈ 0.0 atol=atol
    @test decoff(oq2) ≈ 0.0 atol=atol
    @test projectedseparation(oq2) ≈ 0.0 atol=atol

    @test pmra(oq2) ≈ 0.0 atol=atol
    @test sign(pmdec(oq2)) == -1
    @test radvel(oq2) ≈ 0.0 atol=atol

    @test accra(oq2) ≈ 0.0 atol=atol
    @test accdec(oq2) ≈ 0.0 atol=atol

    # Test orbit values at third quarter
    @test raoff(oq3) ≈ 0.0 atol=atol
    @test decoff(oq3) ≈ -1000.0 rtol=rtol
    @test projectedseparation(oq3) ≈ 1000.0 rtol=rtol

    @test pmra(oq3) ≈ 0.0 atol=atol
    @test pmdec(oq3) ≈ 0.0 atol=atol
    @test radvel(oq3) ≈ -29785.89 rtol=1e-3

    @test accra(oq3) ≈ 0.0 atol=atol 
    @test sign(accdec(oq3)) == +1

    # Test orbit values at fourth quarter
    @test raoff(oq4) ≈ 0.0 atol=atol
    @test decoff(oq4) ≈ 0.0 atol=atol
    @test projectedseparation(oq4) ≈ 0.0 atol=atol
    
    @test pmra(oq4) ≈ 0.0 atol=atol
    @test sign(pmdec(oq4)) == +1
    @test radvel(oq4) ≈ 0.0 atol=atol

    @test sign(accra(oq4)) == +1
    @test accdec(oq4) ≈ 0.0 atol=atol

    # Compare velocities and accelerations
    @test pmdec(oq2) ≈ -pmdec(oq4) rtol=rtol
    @test accdec(oq1) ≈ -accdec(oq3) rtol=rtol
end

## Test varying eccentricity
@testset "Eccentricity" begin
    # Basic eccentric orbit
    eccentric_1AU_1Msun_1pc = orbit(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(eccentric_1AU_1Msun_1pc, one_year_range)
    ys = decoff.(eccentric_1AU_1Msun_1pc, one_year_range)
    ps = projectedseparation.(eccentric_1AU_1Msun_1pc, one_year_range)

    @test period(eccentric_1AU_1Msun_1pc) == 1.0*PlanetOrbits.year2day
    @test distance(eccentric_1AU_1Msun_1pc) == 1
    
    # Mean motion should be the same
    @test PlanetOrbits.meanmotion(eccentric_1AU_1Msun_1pc) == 2π

    # The separation should now be varying
    # By definition of eccentricity 0.5, 1AU and 1PC
    @test maximum(ps) ≈ 1500 rtol=rtol
    @test minimum(ps) ≈ 500 rtol=rtol

    # When argument of periapsis and periastron are both zero, periastron should be in the East, apoastron in the West
    @test maximum(ys) ≈ 500 rtol=rtol
    @test minimum(ys) ≈ -1500 rtol=rtol

    # Rotate Ω
    ecc_rot_Ω = orbit(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = 0.0,
        Ω = deg2rad(90),
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_Ω, one_year_range)
    ys = decoff.(ecc_rot_Ω, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test minimum(xs) ≈ -1500 rtol=rtol
    @test maximum(xs) ≈ 500 rtol=rtol

    # Rotate τ
    ecc_rot_ω = orbit(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = deg2rad(90.0),
        Ω = 0.0,
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_ω, one_year_range)
    ys = decoff.(ecc_rot_ω, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test minimum(xs) ≈ -1500 rtol=rtol
    @test maximum(xs) ≈ 500 rtol=rtol

    # Rotate Ω & τ
    ecc_rot_Ωτ = orbit(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = deg2rad(-90),
        Ω = deg2rad(90),
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc_rot_Ωτ, one_year_range)
    ys = decoff.(ecc_rot_Ωτ, one_year_range)
    # Recall, East is left in the sky.
    # We have rotated  90 degrees CCW.
    @test maximum(ys) ≈ 500 rtol=rtol
    @test minimum(ys) ≈ -1500 rtol=rtol

    # Highly eccentric 
    ecc09 = orbit(
        a = 1.0, # AU
        e = 0.9,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc09, one_year_range)
    ys = decoff.(ecc09, one_year_range)
    ps = projectedseparation.(ecc09, one_year_range)
    # Loosen the tolerance on these
    @test maximum(ps) ≈ 1900 rtol=1e-4
    @test minimum(ps) ≈ 100 rtol=1e-4

    # Extremely eccentric 
    ecc09 = orbit(
        a = 1.0, # AU
        e = 1-1e-3,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        tp = 0.0,
        M = 1.0, # M_sun
        plx = 1000.0, # 1000 mas == 1pc
    )
    xs = raoff.(ecc09, one_year_range)
    ys = decoff.(ecc09, one_year_range)
    ps = projectedseparation.(ecc09, one_year_range)
    @test maximum(ps) ≈ 1999 rtol=1e-4
    # Loosen the tolerance on these even more (periastron flies by very quickly)
    @test minimum(ps) ≈ 1 rtol=1e1
end 

## Test chain rules
@testset "Chain Rules" begin
    # These tests are broken at MA===0, e>0

    # First test analytic chain rules
    k1(MA) = e->PlanetOrbits.kepler_solver(MA, e)
    k2(e) = MA->PlanetOrbits.kepler_solver(MA, e)
    
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

## Test analytic derivatives match numeric derivatives
@testset "PMA & Accel." begin

    # Check analytic derivative properties against ForwardDiff over a big range of orbits
    for t in 0.:35:356.,
        a in 0.1:0.2:3,
        e in 0:0.1:0.9,
        i in deg2rad.([-45, 0, 45, 90, ]),
        ω in deg2rad.([-45, 0, 45, 90, ]),
        Ω in deg2rad.([-45, 0, 45, 90, ])

        elems = orbit(;
            a,
            e,
            i = 0.0,
            ω = 0.0,
            Ω = 0.0,
            tp = 0.0,
            M = 1.0,
            plx = 1000.0, # 1000 mas <-> 1pc
        )

        @test pmra(elems, 100.0) ≈ ForwardDiff.derivative(
            t->raoff(elems, t),
            100.0
        )*PlanetOrbits.year2day

        @test pmdec(elems, 100.0) ≈ ForwardDiff.derivative(
            t->decoff(elems, t),
            100.0
        )*PlanetOrbits.year2day

        @test accra(elems, 100.0) ≈ ForwardDiff.derivative(
            t->pmra(elems, t),
            100.0
        )*PlanetOrbits.year2day

        @test accdec(elems, 100.0) ≈ ForwardDiff.derivative(
            t->pmdec(elems, t),
            100.0
        )*PlanetOrbits.year2day    
    end
end



@testset "Orbit selection" begin
    @test typeof(orbit(;a=1.0, e=0.0, ω=0.0, tp=0.0, M=1.0)) <: RadialVelocityOrbit
    @test typeof(orbit(;a=1.0, e=0.0, ω=0.0, tp=0.0, M=1.0, i=0.1, Ω=0.0)) <: KepOrbit
    @test typeof(orbit(;a=1.0, e=0.0, ω=0.0, tp=0.0, M=1.0, i=0.1, Ω=0.0, plx=100.0).parent) <: KepOrbit
    @test typeof(orbit(;A=100.0, B=100.0, F=100.0, G=-100.0, e=0.5, tp=0.0, M=1.0, plx=100.0)) <: ThieleInnesOrbit
end

# ----------------------------------------------------------------------------------------------------------------------
