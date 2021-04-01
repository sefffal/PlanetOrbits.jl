import Random
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
Random.seed!(1234)
images_contrasts = map(eachrow(points)) do (ra,dec)
    x = -ra
    y = dec

    img = zeros(201,201)
    r = imgsep(img)
    img[r .< 2] .= NaN

    img = map(zip(img, r)) do (px,r)
        px + 2000randn()/0.5r
    end

    img = centered(img)

    # img = imfilter(img, Kernel.gaussian(5), NA())
    img_for_contrast = imfilter(img, Kernel.gaussian(5), "replicate")
    contrast = contrast_interp(img_for_contrast)

    img[round(Int,x/10), round(Int,y/10)] += 5000
    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]

# Define our priors using any Distrubtions
priors = (;
    f = TruncatedNormal(20, 5, 0., Inf),
    a = TruncatedNormal(15, 4, 0., Inf),
    i = Normal(0.6, 0.3),
    e = TruncatedNormal(0.2, 0.2, 0.0, 0.99),
    # τ = Normal(150, 200.0),
    # τ = Uniform(0,1),
    τ = Uniform(0,period(truth_elements)*3),
    ω = Normal(0.0, 0.3),
    Ω = Normal(0.0, 0.3),
)

##


# Rule of thumb from Nienke
@time chains = DirectOrbits.fit_images_emcee(
# @time chains = DirectOrbits.fit_images_RAM(
    priors,
    static,
    images,
    contrasts,
    times,
    platescale=10.,
    burnin=5_000_000,
    numsamples_perwalker=60_000,
    numwalkers=3length(priors)
)

# displaying chains will give summary statistics on all fitted elements
# You can also use StatsPlots for traceplots, histograms, basic corner plots, etc.

## Corner plot (python)
using DelimitedFiles
prepared = hcat(
    chains[:f][:],
    chains[:a][:],
    chains[:τ][:],
    rad2deg.(chains[:i][:]),
    rad2deg.(chains[:Ω][:]),
    chains[:e][:],
    rad2deg.(chains[:ω][:]),
)

writedlm("chains-2.txt", prepared)

##
run(`python examples/make-corner-plot.py`)

## Sample from posterior and make a nice plot

using Plots
theme(:dao)
N = 1000
sampled = sample(KeplerianElements, chains, static, N)

# plot(dpi=200, legend=:topright)
subplots = map(images[1:9]) do image
    i = DirectImage(image)
    i.PLATESCALE = 10.
    p = imshow(i, skyconvention=true, τ=6, colorbar=:none, margin=-1Plots.mm)
    xlims!(p, -1000,1000)
    ylims!(p, -1000,1000)
    # scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
    # plot!(p, sampled, color=:white, alpha=0.03, label="",)
    plot!(p, sampled, color=:white, alpha=0.01, label="",)
    plot!(p, truth_elements, color=:black, label="",)
    scatter!(p, [0],[0], marker=(:star, :black,6),label="")
    xlabel!(p, "")
    ylabel!(p, "")
    p
end
plot(subplots..., layout=(3,3), size=(700,700) )
# savefig("image-orbit-plot.svg")



##
N = 1500
sampled = sample(KeplerianElements, chains, static, N)
ra = raoff.(sampled, times')
dec = decoff.(sampled, times') 
plot()
histogram2d!(
    ra,
    dec,
    bins=(-1000:25:1000, -1000:25:1000),
    color=:plasma,
    label="Sampled points",
    xflip=true,
    xlims=(-1000,1000),
    ylims=(-1000,1000),
    background=:black,
    foreground=:white
)

#savefig("temp-plot-2.png")

## Additional chains diagnostics
# using StatsPlots
# using MCMCChains: meanplot, traceplot

# meanplot(chains)
# traceplot(chains)
##
using StatsBase
h = fit(Histogram, 
    (chains[:f][:], chains[:a][:],),
    nbins=(100,100)
)

h.weights |>ds9show