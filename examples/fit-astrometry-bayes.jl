using DirectOrbits
using Distributions


# We will keep these elements constant
static = (;
    μ = 1,
    ω = 0,
    Ω = 0,
    plx = 45,
)

# Generate astrometry points using this template
# and we will try to fit the results
truth = (;
    a = 15,
    i = 0.5,
    e = 0.2,
    τ = 100,
)
truth_elements = KeplerianElements(merge(truth, static))
times = range(0, 365.24*25, length=10)
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Add some error to the measurements
astrom_err = 40.
points .+= rand(Normal(0, astrom_err), size(points))

# And this is also the input uncertainty
uncertainty = fill(astrom_err, size(points))

# Define our priors using any Distrubtions
priors = (;
    a = TruncatedNormal(25, 4, 0., Inf),
    i = Normal(0.7, 0.3),
    e = TruncatedNormal(0.5, 0.2, 0.0, 1.0),
    τ = Normal(250, 200.0),
)

## Run 

chains = DirectOrbits.fit_bayes(priors, static, points, times, uncertainty, burnin=30_000, numsamples_perwalker=20_000)

# displaying chains will give summary statistics on all fitted elements
# You can also use StatsPlots for traceplots, histograms, basic corner plots, etc.

## Corner Plot
# This step unfortunately requires the python corner library 
# until a nice one is made in Julia. Installing this might
# be tricky depending on your setup.
using PyCall
corner = pyimport("corner")
import PyPlot # Necessary for plots to auto-display

# Reorganize the samples, subset every 10th, and plot.
prepared = hcat(
    chains[:a][:], # a
    chains[:e][:], # e
    rad2deg.(chains[:i][:]), # i
    chains[:τ][:], # τ
)
figure = corner.corner(
    prepared,
    labels=["a - au", "ecc", "inc - °", "τ - days"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=true, title_kwargs=Dict("fontsize"=>12),
);
display(figure)

## Sample from posterior and make a nice plot
using Plots

N = 100
sampled = sample(KeplerianElements, chains, static, N)

plot()
scatter!(points[:,1], points[:,2], yerr=uncertainty[:,1], xerr=uncertainty[:,2], label="Astrometry", marker=(:black, :circle,1))
plot!(sampled, label="Posterior", color=1, alpha=0.05)
scatter!([0],[0], marker=(:star, :yellow,6),label="")