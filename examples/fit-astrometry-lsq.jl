
##

# @testset "Fitting" begin

# times = range(0, 365.24, length=20)
times = range(0, 365.24*25, length=10) |> collect
times .+= 100randn(size(times))



    
using Parameters
using ComponentArrays
using Distributions
using Plots
theme(:dao)

using Optim 

static = (;
    μ = 1.,
    # ω = 0,
    # Ω = 0,
    plx = 45.,
)

# Generate the points we will try to fit
truth = ComponentArray{Float64}(
    a = 15,
    i = 0.5,
    e = 0.0,
    τ = 0,
    ω = deg2rad(35),
    Ω = deg2rad(123)
)
truth_elements = KeplerianElements(;convert(NamedTuple,truth)..., static...)
points = hcat(raoff.(truth_elements, times), decoff.(truth_elements, times))
astrom_err = 10
points .+= rand(Normal(0, astrom_err), size(points))

# points = [
#     -860.348  -145.065
#     -151.847   995.722
#      871.223   158.71
#      126.111  -989.545
# ]

##
initial = (;
    a = 25.0,
    i = 0.9,
    e = 0.3,
    τ = 0.0,
    ω = deg2rad(24),
    Ω = deg2rad(100)
)
initial_elements = KeplerianElements(merge(initial, static)...)


##
# @unpack mle = DirectOrbits.fit_mle(points, times, static, initial)
@time @unpack mle, trace = DirectOrbits.fit_lsq(points, times, static, initial, trace=true)
mle

##
plot()
plot!(initial_elements, legend=:topleft, label="Initial");
plot!(trace, label="Trace", color=2)
plot!(mle, label="Converged", color=3)
scatter!(eachcol(points)..., color=:black, label="Astrometry")
# scatter!(raoff.(mle,times), decoff.(mle,times), color=:yellow, marker=(:cross, 5),  label="Converged",)
scatter!([0], [0], marker=(:star, :black, 5,), label="")

##
savefig(raw"C:\Users\William\OneDrive\Documents\2021\Orbit-Fitting\images\converge-optim-2.svg")
