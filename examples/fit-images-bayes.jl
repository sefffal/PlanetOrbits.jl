import Random
using DirectOrbits
using Distributions
using DirectImages
using Plots, PairPlots
theme(:dao)

using MCMCChains: Chains
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
    a = 12,
    τ = 0.25,
    ω = 0,
    Ω = 0,
    e = 0.2,
    i = 0.5,
)
truth_elements = KeplerianElements(merge(truth, static))

truth2 = (;
    f = 15.,
    a = 18,
    τ = 0.75,
    ω = 0.1,
    Ω = 0.0,
    e = 0.25,
    i = 0.55,
)
truth_elements2 = KeplerianElements(merge(truth2, static))

times = range(0, period(truth_elements)/4, length=10, )
points = hcat(
    raoff.(truth_elements, times),
    decoff.(truth_elements, times)
)

times2 = range(0, period(truth_elements2)/4, length=10, )
points2 = hcat(
    raoff.(truth_elements2, times2),
    decoff.(truth_elements2, times2)
)

# Create synthetic images at each time with those points
Random.seed!(1234)
images_contrasts = map(zip(eachrow(points),eachrow(points2))) do ((ra,dec),(ra2,dec2))
    x = -ra
    y = dec
    x2 = -ra2
    y2 = dec2

    img = zeros(201,201)
    r = imgsep(img)
    img[r .< 2] .= NaN

    img = map(zip(img, r)) do (px,r)
        # px + 2000randn()/0.5r
        # px + 900randn()/0.5r
        px + 200randn()/0.5r
    end

    img = centered(img)

    # img = imfilter(img, Kernel.gaussian(5), NA())
    img_for_contrast = imfilter(img, Kernel.gaussian(5), "replicate")
    contrast = contrast_interp(img_for_contrast)

    img[round(Int,x/10), round(Int,y/10)] += 5000
    # img[round(Int,x2/10), round(Int,y2/10)] += 5000

    img = imfilter(img, Kernel.gaussian(5), "replicate")

    img, contrast
end
images = [img for (img,contrast) in images_contrasts]
contrasts = [contrast for (img,contrast) in images_contrasts]

# Define our priors using any Distrubtions
priors = (;
    f = TruncatedNormal(28, 5, 0., Inf),
    a = TruncatedNormal(16, 8, 4., Inf),
    i = Normal(0.6, 0.3),
    e = TruncatedNormal(0.21, 0.2, 0.0, 0.9999),
    τ = Uniform(0,1),
    ω = Normal(0.0, 0.3),
    Ω = Normal(0.0, 0.3),
)


## Visualze starting position
# initial = KeplerianElements(merge(namedtuple(keys(priors), DirectOrbits.find_starting_point(ln_post, priors)), static))
initial = map(namedtuple.(Ref(keys(priors)), eachcol(DirectOrbits.find_starting_walkers(ln_post, priors, 15)))) do nt
    KeplerianElements(merge(nt, static))
end


N = 250
sampled = sample(KeplerianElements, chains, static, N)

i = DirectImage(images[1])
i.PLATESCALE = 10.
# p = imshow(i, skyconvention=true, τ=6, legend=:topleft)
p = imshow(i, skyconvention=true, clims=(-10,maximum(i)), legend=:topleft)
xlims!(p, -1000,1000)
ylims!(p, -1000,1000)
# scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
plot!(p, initial, label="initial",)
plot!(p, truth_elements, color=:black, label="truth",)
# plot!(p, truth_elements2, color=:black, label="",)
scatter!(p, [0],[0], marker=(:star, :black,6),label="")

##


# Rule of thumb from Nienke
# @time chains = DirectOrbits.fit_images_emcee(
# @time chains = DirectOrbits.fit_images_RAM(
#     priors,
#     static,
#     images,
#     contrasts,
#     times,
#     platescale=10.,
#     burnin=150_000,
#     numsamples_perwalker=25000,
#     # numwalkers=3length(priors)
#     numwalkers=1
# )

# @time chains = DirectOrbits.fit_images_emcee(
#     priors,
#     static,
#     images,
#     contrasts,
#     times,
#     platescale=10.,
#     burnin=1000,
#     numsamples_perwalker=20_000,
#     numwalkers=100
# )

@time chains = DirectOrbits.fit_images_kissmcmc(
    priors,
    static,
    images,
    contrasts,
    times,
    platescale=10.,
    burnin=10_000,
    numwalkers=500,
    numsamples_perwalker=20_000,
)

# @time full_chains = DirectOrbits.fit_images_NUTS(
#     priors,
#     static,
#     images,
#     contrasts,
#     times,
#     platescale=10.,
#     numwalkers=5,
#     numsamples_perwalker=100000,
#     # numsamples_perwalker=1_000,
# )
# chains = full_chains
# chains = Chains(Array(full_chains[50000:end,:,:];), keys(full_chains))

# displaying chains will give summary statistics on all fitted elements
# You can also use StatsPlots for traceplots, histograms, basic corner plots, etc.


##
# Sample from posterior and make a nice plot

N = 250
sampled = sample(KeplerianElements, chains, static, N)

i = DirectImage(images[1])
i.PLATESCALE = 10.
p = imshow(i, skyconvention=true, clims=(-10,maximum(i)), legend=:topleft)
xlims!(p, -1000,1000)
ylims!(p, -1000,1000)
# scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
plot!(p, sampled, alpha=0.1, color=:white, label="",)
plot!([], [], alpha=1, color=:white, label="samples",)
plot!(p, truth_elements, color=:black, label="truth",)
# plot!(p, truth_elements2, color=:black, label="",)
scatter!(p, [0],[0], marker=(:star, :black,6),label="")


##
N = 500
sampled = sample(KeplerianElements, chains, static, N)
ra = raoff.(sampled, times')
dec = decoff.(sampled, times') 
plot()
plot!(p, sampled, alpha=0.1, color=:black, label="",)
histogram2d!(
    ra,
    dec,
    bins=(-1000:25:1000, -1000:25:1000),
    color=:plasma,
    label="Sampled points",
    xflip=true,
    xlims=(-1000,1000),
    ylims=(-1000,1000),
    # background=:black,
    # foreground=:white
)


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

writedlm("chains-6.txt", prepared)

##
#run(`python examples/make-corner-plot.py`)

# p = PairPlots.corner(chains, [
#     raw"f",
#     raw"a",
#     raw"i",
#     raw"e",
#     raw"\tau",
#     raw"\omega",
#     raw"\Omega",
# ], plotscatter=false, plotcontours=false)



chns = (chains[:f][:],
chains[:a][:],
chains[:e][:],
chains[:i][:],
chains[:τ][:],
chains[:ω][:],
chains[:Ω][:])
ss = std.(chns)
masks =  mapreduce(hcat, zip(chns,ss)) do (c,s)
    -2s .< (c .- median(c)) .< 2s
end
mask = prod(masks, dims=2)[:]
count(mask)


p = PairPlots.corner(
    (
        f=chains[:f][:][mask],
        a=chains[:a][:][mask],
        e=chains[:e][:][mask],
        i=chains[:i][:][mask],
        tau=chains[:τ][:][mask],
        # omega=chains[:ω][:][mask],
        # Omega=chains[:Ω][:][mask],
    ),
    plotscatter=false
)

##
function bp(kw) 
    kw2 = delete(kw, :inset)
    N = 250
    sampled = sample(KeplerianElements, chains, static, N)
    i = DirectImage(images[1])
    i.PLATESCALE = 10.
    DirectImages.imshow!(i; lims=1000, framestyle=:none, skyconvention=true, τ=6, legend=:topleft, colorbar=:none, kw...)
    ## scatter!(points[:,1], points[:,2], color=:black, label="Astrometry")
    plot!(sampled; alpha=1, color=:white, label="", kw2...)
    plot!([], []; alpha=1, color=:white, label="samples", kw2...)
    plot!(truth_elements; color=:black, label="truth", kw2...)
    scatter!([0],[0]; marker=(:star, :black,6),label="", kw2...)
end
corner(chains, [
    raw"f",
    raw"a",
    raw"\tau",
    raw"i",
    raw"\Omega",
    raw"e",
    raw"\omega"
], bonusplot=bp)




## Sample from posterior and make a nice plot

using Plots
theme(:dao)
N = 200
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
    plot!(p, sampled, color=:white, label="",)
    plot!(p, truth_elements, color=:black, label="",)
    scatter!(p, [0],[0], marker=(:star, :black,6),label="")
    xlabel!(p, "")
    ylabel!(p, "")
    p
end
plot(subplots..., layout=(3,3), size=(700,700) )
# savefig("image-orbit-plot.svg")



#savefig("temp-plot-2.png")

## Additional chains diagnostics
# using StatsPlots
using MCMCChains: meanplot, traceplot

# meanplot(chains)
traceplot(chains)
savefig("traceplot.png")
##
using StatsBase
h = fit(Histogram, 
    (chains[:f][:], chains[:a][:],),
    nbins=(100,100)
)

h.weights |>ds9show