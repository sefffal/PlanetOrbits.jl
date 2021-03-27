ENV["MPLBACKEND"]="svg"
using DirectOrbits
using Distributions
using DirectImages

using ImageFiltering

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
    f = 20.,
    a = 15,
    i = 0.5,
    e = 0.2,
    τ = 100,
)
truth_elements = KeplerianElements(merge(truth, static))
times = range(0, period(truth_elements), length=10, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Create images at each time with those points
images = map(eachrow(points)) do (ra,dec)
    x = -ra
    y = dec

    # img = centered(20randn(2000,2000))
    img = centered(zeros(2001,2001))
    img[round(Int,x), round(Int,y)] += 50000

    imgf = imfilter(img, Kernel.gaussian(20), NA())
    map(imgf) do px
        px + randn()
    end
end
contrasts = contrast_interp.(images)

# Define our priors using any Distrubtions
priors = (;
    f = TruncatedNormal(20, 5, 0., Inf),
    a = TruncatedNormal(15, 4, 0., Inf),
    i = Normal(0.6, 0.3),
    e = TruncatedNormal(0.2, 0.2, 0.0, 1.0),
    τ = Normal(250, 200.0),
)
##
ll = DirectOrbits.make_ln_like_images(keys(priors), static, images, contrasts, times)
ll(truth)
# ll((f=40, a=15, i=0.5, e=0.2, τ=100))
## Run 

chains = DirectOrbits.fit_images(priors, static, images, contrasts, times, burnin=50_000, numsamples_perwalker=30_000)

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
# We have to use awkward python syntax to save the corner plot
figure.savefig("temp-corner.png", dpi=200)

## Sample from posterior and make a nice plot
using Plots
theme(:dao)
N = 100
sampled = sample(KeplerianElements, chains, static, N)

# plot(dpi=200, legend=:topright)
imshow(images[2],skyconvention=true)
xlims!(-1000,1000)
ylims!(-1000,1000)
scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
plot!(sampled, label="Posterior", color=2, alpha=0.05)
scatter!([0],[0], marker=(:star, :black,6),label="")
##
savefig("temp-plot.png")