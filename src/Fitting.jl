using ComponentArrays
using Optim
using Distributions: logpdf
using AffineInvariantMCMC
using RobustAdaptiveMetropolisSampler
using MCMCChains
using NamedTupleTools
using DirectImages: lookup_coord
using Base.Threads: @threads
# Least squares astrometry fitting

function least_squares_distance(elements, obs, times)

    total = zero(eltype(obs))
    for (pos,t) in zip(eachrow(obs), times)
        x = raoff(elements,t)
        y = decoff(elements,t)
        total += sqrt((x-pos[1])^2 + (y-pos[2])^2)
    end

    return total
end


function fit_lsq(points, times, static, initial; trace=false)
    initial_ca = ComponentArray{Float64}(initial);
    initial_elements = KeplerianElements(merge(initial, static))


    objective = let static=static, points=points, times=times
        function (params)
            el = KeplerianElements(merge(convert(NamedTuple, params), static))
            return log(DirectOrbits.least_squares_distance(el, points, times))
        end
    end


    results = Optim.optimize(objective, initial_ca, NelderMead(), Optim.Options(store_trace=trace, extended_trace=trace))

    optimized_elements =  KeplerianElements(merge(convert(NamedTuple, Optim.minimizer(results)), static))

    # Gather elements from trace
    if trace
        elements_trace = map(results.trace) do trace_step
            if haskey(trace_step.metadata, "centroid")
                KeplerianElements(;convert(NamedTuple, trace_step.metadata["centroid"])..., static...)
            else
                KeplerianElements(;convert(NamedTuple, trace_step.metadata["x"])..., static...)
            end
        end

        return (;mle=optimized_elements, trace=elements_trace)
    else
        return (;mle=optimized_elements)
    end
end



## Bayesian astrometry fitting

# This is a straight forward implementation that unfortunately is not type stable.
# This is because we are looping over a heterogeneous tuple
# function make_ln_prior(priors...)
#     function ln_prior(params)
#         lp = zero(first(params))
#         for i in eachindex(params)
#             pd = priors[i]
#             param = params[i]
#             lp += logpdf(pd, param)
#         end
#         return lp 
#     end
#     return ln_prior
# end

# Here is the same exact code, but manually unrolled.
# This does limit the max number of priors we can support arbitrarily but is ~7x faster
function make_ln_prior(p1, p2=nothing, p3=nothing, p4=nothing, p5=nothing, p6=nothing, p7=nothing, p8=nothing, p9=nothing)
    function ln_prior(params)
        lp = zero(first(params))
        lp += logpdf(p1, params[1])
        if !isnothing(p2)
            lp += logpdf(p2, params[2])
        end
        if !isnothing(p3)
            lp += logpdf(p3, params[3])
        end
        if !isnothing(p4)
            lp += logpdf(p4, params[4])
        end
        if !isnothing(p5)
            lp += logpdf(p5, params[5])
        end
        if !isnothing(p6)
            lp += logpdf(p6, params[6])
        end
        if !isnothing(p7)
            lp += logpdf(p7, params[7])
        end
        if !isnothing(p8)
            lp += logpdf(p8, params[8])
        end
        if !isnothing(p9)
            lp += logpdf(p9, params[9])
        end
        return lp
    end
    return ln_prior
end


function make_ln_like_astrom(props, static, points, times, uncertainties)
    function ln_like(params)
        nt = namedtuple(props, params)
        elements = KeplerianElements(merge(nt, static))

        ll = zero(first(params))
        for (point, errs, t) in zip(eachrow(points), eachrow(uncertainties), times)
            x, y = kep2cart(elements, t)
            residx = point[1] - x
            residy = point[2] - y
            σ²x = errs[1]^2
            σ²y = errs[2]^2
            χ²x = -0.5residx^2 / σ²x - log(sqrt(2π*σ²x))
            χ²y = -0.5residy^2 / σ²y - log(sqrt(2π*σ²y))
            ll += χ²x + χ²y
        end

        # I think there is a normalization term missing that is causing the likelihood
        # to completely dominate the priors.

        return ll
    end
    return ln_like
end

function make_ln_post_astrom(priors, static, points, times, uncertainties)
    
    ln_prior = make_ln_prior(priors...)
    ln_like = make_ln_like_astrom(keys(priors), static, points, times, uncertainties)

    ln_post(params) = ln_prior(params) + ln_like(params)

    return ln_post
end


function fit_bayes(
    priors, static, points, times, uncertainties;
    numwalkers=10,
    burnin = 10_000,
    thinning = 1,
    numsamples_perwalker = 10_000
    )

    # Initial values for the walkers are drawn from the priors
    initial_walkers = reduce(hcat, [rand.([priors...,]) for _ in 1:numwalkers])

    ln_post = make_ln_post_astrom(priors, static, points, times, uncertainties)

    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, numwalkers, initial_walkers, burnin, 1);
    @info "Burn in done"
    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, size(chain,2), chain[:, :, end], numsamples_perwalker, thinning);
    @info "Chains done"

    column_names = string.(collect(keys(priors)))
    chain = Chains(permutedims(chain, (3, 1, 2)),  column_names);

    return chain
end


struct TupleProtoHolder{T<:Any}
end


function make_ln_like_images(props, static, images, contrasts, times, platescale)
    # props is tuple of symbols
    # static is named tuple of values

    nt_proto = NamedTuple{props, NTuple{length(props),Float64}}

    if !(size(images)  == size(times) == size(contrasts))
        error("All values must have the same length")
    end

    template_holder = Val{nt_proto}()
    
    function function_barrier(template_holder::Val{T}) where T
        function ln_like(params)
            nt = T(params)
            merged = merge(nt, static)
            elements = KeplerianElements(merged)
            f = merged.f

            ll = zero(f)
            @inbounds for i in eachindex(times)
                image = images[i]
                contrast = contrasts[i]
                t = times[i]
                ra, dec = kep2cart(elements, t)
                x = -ra
                y = dec

                if !isfinite(x) || !isfinite(y)
                    continue
                end

                # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
                ix = round(Int, x/platescale)
                iy = round(Int, y/platescale)
                if ix ∈ axes(image,1) && iy ∈ axes(image,2)
                    f̃ₓ = image[ix,iy]
                else
                    continue
                end
                
                # f̃ₓ = lookup_coord(image, (x,y), platescale)

                r = √(x^2 + y^2)
                σₓ = contrast(r/platescale)

                # When we get a position that falls outside of our available
                # data (e.g. under the coronagraph) we cannot say anything
                # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
                # of zero.
                if !isfinite(σₓ) || !isfinite(f̃ₓ)
                    continue
                end

                l = -1/(2σₓ^2) * (f^2 - 2f*f̃ₓ)

                if !isfinite(l)
                    @show x y r σₓ f f̃ₓ l
                    error("Infinite log-likelihood encountered")
                end

                # Ruffio et al 2017, eqn 31
                ll += l
            end

            return ll
        end
    end
    return function_barrier(template_holder)
end



function make_ln_post_images(priors, static, images, contrasts, times, platescale)
    ln_prior = make_ln_prior(priors...)
    ln_like = make_ln_like_images(keys(priors), static, images, contrasts, times, platescale)
    ln_post(params) = ln_prior(params) + ln_like(params)
    return ln_post
end

function fit_images_emcee(
    priors, static, images, contrasts, times;
    platescale,
    numwalkers=10,
    burnin = 10_000,
    thinning = 1,
    numsamples_perwalker = 10_000
    )
    column_names = string.(collect(keys(priors)))

    # Initial values for the walkers are drawn from the priors
    initial_walkers = reduce(hcat, [rand.([priors...,]) for _ in 1:numwalkers])

    ln_post = make_ln_post_images(priors, static, images, contrasts, times, platescale)

    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, numwalkers, initial_walkers, burnin, 1);
    @info "Burn in done"
    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, size(chain,2), chain[:, :, end], numsamples_perwalker, thinning);
    @info "Chains done"

    chain = Chains(permutedims(chain, (3, 1, 2)),  column_names);

    return chain
end


function fit_images_RAM(
    priors, static, images, contrasts, times;
    platescale,
    numwalkers=10,
    burnin = 10_000,
    thinning = 1,
    numsamples_perwalker = 10_000
    )

    # Initial values for the walkers are drawn from the priors
    # initial_walkers = reduce(hcat, [rand.([priors...,]) for _ in 1:numwalkers])

    column_names = string.(collect(keys(priors)))

    ln_post = make_ln_post_images(priors, static, images, contrasts, times, platescale)

    
    ## RAM Sampler
    # using RobustAdaptiveMetropolisSampler
    chains = map(1:numwalkers) do walker_i
        initial_walker = rand.([priors...,])

        chain, accrate, S = RAM_sample(
            ln_post,
            initial_walker,
            0.5,
            burnin,
            show_progress=true
        )
        @show size(chain)
        @info "Burn in done $walker_i"
        chain, accrate, S = RAM_sample(
            ln_post,
            chain[end,:],
            0.5,
            numsamples_perwalker,
            show_progress=true
        )
        @info "Chains done $walker_i"
        return chain
    end
    chain = Chains(cat(chains...,dims=3), column_names)

    # chain = Chains(reshape(chain, size(chain,1), size(chain,2), 1),  column_names);

    return chain
end




# Analysis functions
import StatsBase: sample

# The stats base sample function makes it easy to get values from Chains
# but converting these values into KeplerianElements along with any static
# paramters takes a bit of work.
# Here we overload the sample function with a method just for this.
function sample(::Type{KeplerianElements}, chains::Chains, static, N=1)
    sampled = sample(chains, ceil(Int, N/size(chains,2)))
    out = KeplerianElements{Float64}[]

    proto = namedtuple(keys(chains))

    sizehint!(out, size(sampled,1)*size(sampled,2))
    for i in 1:size(sampled,1), j in 1:size(sampled,3)
        nt = proto(Array(sampled[i,:,j]))
        el = KeplerianElements(merge(nt, static))
        push!(out, el)
    end
    return out[begin:min(N,end)]
end
