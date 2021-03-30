ENV["MPLBACKEND"]="svg"
using DirectOrbits
using Distributions
using DirectImages

using ImageFiltering

# We will keep these elements constant
static = (;
    μ = 1,
    plx = 45,
    # τ = 100,
)

# Generate astrometry points using this template
# and we will try to fit the results
truth = (;
    f = 20.,
    a = 15,
    τ = 100,
    ω = 0,
    Ω = 0,
    e = 0.2,
    i = 0.5,
)
truth_elements = KeplerianElements(merge(truth, static))
times = range(0, period(truth_elements)/4, length=10, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

# Create synthetic images at each time with those points
images_contrasts = map(eachrow(points)) do (ra,dec)
    x = -ra
    y = dec

    img = zeros(201,201)
    r = imgsep(img)
    img[r .< 10] .= NaN


    img = map(zip(img, r)) do (px,r)
        px + 4000randn()/0.5r
    end

    img = centered(img)
    img[round(Int,x/10), round(Int,y/10)] += 5000
    # img = imfilter(img, Kernel.gaussian(5), NA())
    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img
end
contrasts = map(eachrow(points)) do _
    x = -ra
    y = dec

    img = zeros(201,201)
    r = imgsep(img)
    img[r .< 30] .= NaN


    img = map(zip(img, r)) do (px,r)
        px + 2000randn()/r
    end

    img = centered(img)
    img = imfilter(img, Kernel.gaussian(5), NA())

    contrast_interp(img)
end

# Define our priors using any Distrubtions
priors = (;
    f = TruncatedNormal(20, 5, 0., Inf),
    a = TruncatedNormal(15, 4, 0., Inf),
    i = Normal(0.6, 0.3),
    e = TruncatedNormal(0.2, 0.2, 0.0, 0.99),
    τ = Normal(150, 200.0),
    ω = Normal(0.0, 0.3),
    Ω = Normal(0.0, 0.3),
)

##


# Rule of thumb from Nienke
N_walk = 3length(priors)

chains = DirectOrbits.fit_images(priors, static, images, contrasts, times, platescale=10., burnin=15_000, numsamples_perwalker=15_000, numwalkers=N_walk)

# displaying chains will give summary statistics on all fitted elements
# You can also use StatsPlots for traceplots, histograms, basic corner plots, etc.

## Corner Plot
# This step unfortunately requires the python corner library 
# until a nice one is made in Julia. Installing this might
# be tricky depending on your setup.
# using PyCall
# corner = pyimport("corner")
# import PyPlot # Necessary for plots to auto-display

# # Reorganize the samples, subset every 10th, and plot.
# prepared = hcat(
#     chains[:f][:],
#     chains[:a][:],
#     chains[:e][:],
#     rad2deg.(chains[:i][:]), 
#     chains[:τ][:],
# )
# figure = corner.corner(
#     prepared,
#     labels=["f", "a - au", "ecc", "inc - °", "τ - days"],
#     quantiles=[0.16, 0.5, 0.84],
#     show_titles=true, title_kwargs=Dict("fontsize"=>12),
# );
# display(figure)
# We have to use awkward python syntax to save the corner plot
# figure.savefig("temp-corner-2.svg", dpi=200)
using DelimitedFiles
prepared = hcat(
    chains[:f][:],
    chains[:a][:],
    rad2deg.(chains[:i][:]),
    chains[:e][:],
    chains[:τ][:],
    rad2deg.(chains[:ω][:]),
    rad2deg.(chains[:Ω][:]),
)

writedlm("chains-2.txt", prepared)

##
run(`python examples/make-corner-plot.py`)

## Sample from posterior and make a nice plot

using Plots
theme(:dao)
N = 300
sampled = sample(KeplerianElements, chains, static, N)

# plot(dpi=200, legend=:topright)
i = DirectImage(images[5])
i.PLATESCALE = 10.
imshow(i, skyconvention=true, τ=20)
xlims!(-1000,1000)
ylims!(-1000,1000)
# scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
plot!(sampled, color=2, alpha=0.1, label="")
scatter!([0],[0], marker=(:star, :black,6),label="")


##
ra = raoff.(sampled, times')
dec = decoff.(sampled, times') 
plot()
histogram2d!(
    ra,
    dec,
    bins=(-1000:25:1000, -1000:25:1000),
    # cgrad=cgrad(:magma, scale=:log10),
    color=:viridis,
    label="Sampled points",
    xflip=true,
    xlims=(-1000,1000),
    ylims=(-1000,1000),
)

#savefig("temp-plot-2.png")

## Additional chains diagnostics
# using StatsPlots
# using MCMCChains: meanplot, traceplot

# meanplot(chains)
# traceplot(chains)