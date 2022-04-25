# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

using Test
using DirectOrbits
using ForwardDiff
using FiniteDiff
import Distributions: Uniform

# ----------------------------------------------------------------------------------------------------------------------
# Constants and Helper Functions
# ----------------------------------------------------------------------------------------------------------------------

# 10 steps per day for one year
one_year_range = 0.0:0.1:365.24
# Relative tolerance for certain tests
rtol = 1e-6
# Absolute tolerance for certain tests
atol = 1e-6

# Randomly generate orbit parameters for tests where fixed values are unimportant
function randomparamsrad()
    a = rand(Uniform(0.1, 50))
    e = rand() 
    i = rand(Uniform(0.0, 2π))
    ω = rand(Uniform(0.0, 2π))
    Ω = rand(Uniform(0.0, 2π))
    τ = rand()
    M = rand(Uniform(0.1, 3.0))
    plx = rand(Uniform(1.0, 1000.0))
    params = (a, e, i, ω, Ω, τ, M, plx)
    return params
end

function randomparamsdeg()
    a = rand(Uniform(0.1, 50))
    e = rand() 
    i = rand(Uniform(0.0, 360.0))
    ω = rand(Uniform(0.0, 360.0))
    Ω = rand(Uniform(0.0, 360.0))
    τ = rand()
    M = rand(Uniform(0.1, 3.0))
    plx = rand(Uniform(1.0, 1000.0))
    params = (a, e, i, ω, Ω, τ, M, plx)
    return params
end

# Randomly generate orbit values for tests where fixed values are unimportant
function randomorbit()
    return (rand(Uniform(0.0, 100.0)) for i ∈ 1:8)
end

# ----------------------------------------------------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------------------------------------------------

## Test relationships between inverse constants
@testset "Constants" begin
    @test DirectOrbits.mas2rad == 1/DirectOrbits.rad2mas
    @test DirectOrbits.as2rad == 1/DirectOrbits.rad2as
    @test DirectOrbits.au2pc == 1/DirectOrbits.pc2au
    @test DirectOrbits.m2au == 1/DirectOrbits.au2m
    @test DirectOrbits.day2year == 1/DirectOrbits.year2day
    @test DirectOrbits.sec2year == 1/DirectOrbits.year2sec
    @test DirectOrbits.sec2day == 1/DirectOrbits.day2sec
end

## Test KeplerianElements attributes match required values
@testset "KeplerianElements Attributes" begin
    a, e, i, ω, Ω, τ, M, plx = randomparamsrad()
    elem = KeplerianElements(a, e, i, ω, Ω, τ, M, plx)
    @test elem.a ≈ a
    @test elem.e ≈ e 
    @test elem.i ≈ i
    @test elem.ω ≈ ω
    @test elem.Ω ≈ Ω 
    @test elem.τ ≈ τ
    @test elem.M ≈ M
    @test elem.plx ≈ plx 
    @test elem.dist ≈ 1000/plx * DirectOrbits.pc2au
    @test elem.T ≈ √(a^3/M) * DirectOrbits.year2day
    @test elem.n ≈ 2π/√(a^3/M)
    @test elem.ν_fact ≈ √((1 + e)/(1 - e))
    @test elem.p ≈ a*(1 - e^2)
    @test elem.cosi ≈ cos(i)
    @test elem.sini ≈ sin(i)
    @test elem.cosω ≈ cos(ω)
    @test elem.sinω ≈ sin(ω)
    @test elem.cosΩ ≈ cos(Ω)
    @test elem.sinΩ ≈ sin(Ω)
    @test elem.ecosω ≈ e*cos(ω)
    @test elem.esinω ≈ e*sin(ω)
    @test elem.cosi_cosΩ ≈ cos(i)*cos(Ω)
    @test elem.cosi_sinΩ ≈ cos(i)*sin(Ω)
    @test elem.J ≈ ((2π*a)/(elem.T*DirectOrbits.day2year)) * (1 - e^2)^(-1//2)
    @test elem.K ≈ elem.J*DirectOrbits.au2m*DirectOrbits.sec2year*sin(i)
    @test elem.A ≈ ((4π^2 * a)/(elem.T*DirectOrbits.day2year)^2) * (1 - e^2)^(-2)
end

## Test standard, keyword, and named tuple KeplerianElements are equal
@testset "KeplerianElements Input Styles" begin
    a, e, i, ω, Ω, τ, M, plx = randomparamsrad()
    nt = (a=a, e=e, i=i, ω=ω, Ω=Ω, τ=τ, M=M, plx=plx)
    elem = KeplerianElements(a, e, i, ω, Ω, τ, M, plx)
    elemkw = KeplerianElements(a=a, e=e, i=i, ω=ω, Ω=Ω, τ=τ, M=M, plx=plx)
    elemnt = KeplerianElements(nt)
    @test elem == elemkw
    @test elemkw == elemnt
    @test elemnt == elem
    @test astuple(elem) == nt
end

## Test KeplerianElementsDeg attributes match required values
@testset "KeplerianElementsDeg Attributes" begin
    a, e, i, ω, Ω, τ, M, plx = randomparamsdeg()
    elem = KeplerianElementsDeg(a, e, i, ω, Ω, τ, M, plx)
    @test elem.a ≈ a
    @test elem.e ≈ e 
    @test elem.i ≈ deg2rad(i) 
    @test elem.ω ≈ deg2rad(ω)
    @test elem.Ω ≈ deg2rad(Ω)
    @test elem.τ ≈ τ
    @test elem.M ≈ M
    @test elem.plx ≈ plx 
    @test elem.dist ≈ 1000/plx * DirectOrbits.pc2au
    @test elem.T ≈ √(a^3/M) * DirectOrbits.year2day
    @test elem.n ≈ 2π/√(a^3/M)
    @test elem.ν_fact ≈ √((1 + e)/(1 - e))
    @test elem.p ≈ a*(1 - e^2)
    @test elem.cosi ≈ cos(deg2rad(i))
    @test elem.sini ≈ sin(deg2rad(i))
    @test elem.cosω ≈ cos(deg2rad(ω))
    @test elem.sinω ≈ sin(deg2rad(ω))
    @test elem.cosΩ ≈ cos(deg2rad(Ω))
    @test elem.sinΩ ≈ sin(deg2rad(Ω))
    @test elem.ecosω ≈ e*cos(deg2rad(ω))
    @test elem.esinω ≈ e*sin(deg2rad(ω))
    @test elem.cosi_cosΩ ≈ cos(deg2rad(i))*cos(deg2rad(Ω))
    @test elem.cosi_sinΩ ≈ cos(deg2rad(i))*sin(deg2rad(Ω))
    @test elem.J ≈ ((2π*a)/(elem.T*DirectOrbits.day2year)) * (1 - e^2)^(-1//2)
    @test elem.K ≈ elem.J*DirectOrbits.au2m*DirectOrbits.sec2year*sin(deg2rad(i))
    @test elem.A ≈ ((4π^2 * a)/(elem.T*DirectOrbits.day2year)^2) * (1 - e^2)^(-2)
end

## Test standard, keyword, and named tuple KeplerianElementsDeg are equal
@testset "KeplerianElementsDeg Input Styles" begin
    a, e, i, ω, Ω, τ, M, plx = randomparamsdeg()
    nt = (a=a, e=e, i=i, ω=ω, Ω=Ω, τ=τ, M=M, plx=plx)
    elem = KeplerianElementsDeg(a, e, i, ω, Ω, τ, M, plx)
    elemkw = KeplerianElementsDeg(a=a, e=e, i=i, ω=ω, Ω=Ω, τ=τ, M=M, plx=plx)
    elemnt = KeplerianElementsDeg(nt)
    @test elem == elemkw
    @test elemkw == elemnt
    @test elemnt == elem
end

# ## Test OrbitSolution attributes match required values
# @testset "OrbitSolution Attributes" begin
#     x, y, ẋ, ẏ, ż, ẍ, ÿ = randomorbit()
#     o = OrbitSolution(x, y, ẋ, ẏ, ż, ẍ, ÿ)
#     @test o.x == x
#     @test o.y == y
#     @test o.ẋ == ẋ
#     @test o.ẏ == ẏ
#     @test o.ż == ż
#     @test o.ẍ == ẍ
#     @test o.ÿ == ÿ
# end

# ## Test standard, keyword, and named tuple OrbitSolution are equal
# @testset "OrbitSolution Input Styles" begin
#     x, y, ẋ, ẏ, ż, ẍ, ÿ = randomorbit()
#     nt = (x=x, y=y, ẋ=ẋ, ẏ=ẏ, ż=ż, ẍ=ẍ, ÿ=ÿ)
#     o = OrbitSolution(x, y, ẋ, ẏ, ż, ẍ, ÿ)
#     okw = OrbitSolution(x=x, y=y, ẋ=ẋ, ẏ=ẏ, ż=ż, ẍ=ẍ, ÿ=ÿ)
#     ont = OrbitSolution(nt)
#     @test o == okw
#     @test okw == ont 
#     @test ont == o
# end

# ## Test operations on OrbitSolution values
# @testset "OrbitSolution Operations" begin
#     x1, y1, ẋ1, ẏ1, ż1, ẍ1, ÿ1 = randomorbit()
#     x2, y2, ẋ2, ẏ2, ż2, ẍ2, ÿ2 = randomorbit()
#     o1 = OrbitSolution(x1, y1, ẋ1, ẏ1, ż1, ẍ1, ÿ1)
#     o1eps = OrbitSolution(x1 + rtol, y1 + rtol,
#                           ẋ1 + rtol, ẏ1 + rtol, ż1 + rtol,
#                           ẍ1 + rtol, ÿ1 + rtol)
#     o2 = OrbitSolution(x2, y2, ẋ2, ẏ2, ż2, ẍ2, ÿ2)
#     @test o1 == o1
#     @test o1 != o1eps
#     @test o1 ≈ o1eps rtol=rtol 
#     @test o1 != o2
#     @test o1 + o2 == OrbitSolution(x1 + x2, y1 + y2,
#                                    ẋ1 + ẋ2, ẏ1 + ẏ2, ż1 + ż2,
#                                    ẍ1 + ẍ2, ÿ1 + ÿ2)
#     @test o1 - o2 == OrbitSolution(x1 - x2, y1 - y2,
#                                    ẋ1 - ẋ2, ẏ1 - ẏ2, ż1 - ż2,
#                                    ẍ1 - ẍ2, ÿ1 - ÿ2)
#     @test -o1 == OrbitSolution(-x1, -y1, -ẋ1, -ẏ1, -ż1, -ẍ1, -ÿ1)
# end

## Idealized face-on Earth with circular orbit at 1 pc 
@testset "Earth, i = 0, e = 0, d = 1 pc" begin
    idealearth = KeplerianElements(
        a = 1.0,
        e = 0.0,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        τ = 0.0,
        M = 1.0,
        plx = 1000.0
    )

    # Test basic orbit properties
    @test period(idealearth) == DirectOrbits.year2day
    @test distance(idealearth) == 1.0
    @test meanmotion(idealearth) == 2π
    @test periastron(idealearth) == 58849
    @test semiamplitude(idealearth) == 0.0

    # Orbit solutions at quarters of the orbit
    oq1 = DirectOrbits.orbitsolve_ν(idealearth, 0.0)
    oq2 = DirectOrbits.orbitsolve_ν(idealearth, π/2)
    oq3 = DirectOrbits.orbitsolve_ν(idealearth, π)
    oq4 = DirectOrbits.orbitsolve_ν(idealearth, 3π/2)

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
    idealearth = KeplerianElements(
        a = 1.0,
        e = 0.0,
        i = π/2,
        ω = 0.0,
        Ω = 0.0,
        τ = 0.0,
        M = 1.0,
        plx = 1000.0
    )

    # Test basic orbit properties
    @test period(idealearth) == DirectOrbits.year2day
    @test distance(idealearth) == 1.0
    @test meanmotion(idealearth) == 2π
    @test periastron(idealearth) == 58849
    @test semiamplitude(idealearth) ≈ 29785.89 rtol=1e-3

    # Orbit solutions at quarters of the orbit
    oq1 = DirectOrbits.orbitsolve_ν(idealearth, 0.0)
    oq2 = DirectOrbits.orbitsolve_ν(idealearth, π/2)
    oq3 = DirectOrbits.orbitsolve_ν(idealearth, π)
    oq4 = DirectOrbits.orbitsolve_ν(idealearth, 3π/2)

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
    eccentric_1AU_1Msun_1pc = KeplerianElements(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        τ = 0.0,
        M = 1.0, # M_sun
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
        e = 0.5,
        i = 0.0,
        ω = 0.0,
        Ω = deg2rad(90),
        τ = 90.0,
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
    ecc_rot_ω = KeplerianElements(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = deg2rad(90.0),
        Ω = 0.0,
        τ = 0.0,
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
    ecc_rot_Ωτ = KeplerianElements(
        a = 1.0, # AU
        e = 0.5,
        i = 0.0,
        ω = deg2rad(-90),
        Ω = deg2rad(90),
        τ = 90.0,
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
    ecc09 = KeplerianElements(
        a = 1.0, # AU
        e = 0.9,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        τ = 0.0,
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
    ecc09 = KeplerianElements(
        a = 1.0, # AU
        e = 1-1e-3,
        i = 0.0,
        ω = 0.0,
        Ω = 0.0,
        τ = 0.0,
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
    k1(MA) = e->DirectOrbits.kepler_solver(MA, e)
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

## Test analytic derivatives match numeric derivatives
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
            e,
            i = 0.0,
            ω = 0.0,
            Ω = 0.0,
            τ = 0.0,
            M = 1.0,
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

        @test accra(elems, 100.0) ≈ ForwardDiff.derivative(
            t->pmra(elems, t),
            100.0
        )*DirectOrbits.year2day

        @test accdec(elems, 100.0) ≈ ForwardDiff.derivative(
            t->pmdec(elems, t),
            100.0
        )*DirectOrbits.year2day    
    end
end

# ----------------------------------------------------------------------------------------------------------------------