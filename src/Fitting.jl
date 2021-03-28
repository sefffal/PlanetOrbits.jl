using ComponentArrays
using Optim
using Distributions: logpdf
using AffineInvariantMCMC
using MCMCChains

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


function make_ln_prior(priors)
    function ln_prior(params)
        lp = zero(first(params))
        for i in eachindex(params)
            pd = priors[i]
            param = params[i]
            lp += logpdf(pd, param)
        end
        return lp 
    end
    return ln_prior
end
function make_ln_like_astrom(props, static, points, times, uncertainties)
    function ln_like(params)
        expando = (;(prop=>param for (prop,param) in zip(props, params))...)
        elements = KeplerianElements(merge(expando, static))

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
    
    ln_prior = make_ln_prior(priors)
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



function make_ln_like_images(props, static, images, contrasts, times, platescale)
    function lnlike(params)
        expando = (;(prop=>param for (prop,param) in zip(props, params))...)
        merged = merge(expando, static)
        elements = KeplerianElements(merged)
        f = merged.f

        ll = zero(f)
        for (image, contrast, t) in zip(images, contrasts, times)
            ra, dec = kep2cart(elements, t)
            x = -ra
            y = dec

            ix = round(Int, x/platescale)
            iy = round(Int, y/platescale)
            if ix ∈ axes(image,1) && iy ∈ axes(image,2)
                f_img = image[ix,iy]
            else
                f_img = zero(eltype(image))
            end

            r = √(x^2 + y^2)
            σ = contrast(r)

            # Ruffio et al 2017, eqn 31
            ll += -1/(2σ^2) * (f^2 - 2f*f_img)

            # ll += 1/2σ * (f^2 - 2f*f_img)
            # l = 1/(√(2π)*σ) * exp(-1/2 * (f - f_img)^2/σ^2)
            # ll += log(l)
        end

        return ll
    end

        
        #         # Construct an orbit with these parameters
        #         orbit = Orbit(a, i, e, τ, μ, ω, Ω, plx)
        #         # Then get the liklihood by multipliying together
        #         # the liklihoods at each epoch (summing the logs)
        #         ln_post = 0.0
        #         # Go through each image
        #         for I in eachindex(convolved)
        #             # Find the current position in arcseconds
        #             pos = SVector(0., 0., 0.)
        #             try
        #                 pos = xyz(orbit, mjds[I])
        #             catch e
        #                 if !(typeof(e) <: Real)
        #                     @warn "error looking up orbit" exception = e maxlog = 2
        #                 end
        #             end
        #             # Get the value in the convolved images
        #             pos = SVector(pos[1], -pos[2])
        #             # I believe that the x-coordinate should be flipped because
        #             # image indices are opposite to sky coordinates
        
        #             # pos = SVector(-pos[2], pos[1])
        #             # pos = SVector(pos[2], -pos[1])
        
        #             phot_img = lookup_coord(convolved[I], pos, platescale)
        
        #             # NOTE: experimenting with flipped coords!
        #             # pos# .* SVector(1, 1, 1)
        #             # pos = SVector(-pos[2], pos[1])
        #             # phot_img = lookup_coord(convolved[I], pos, platescale)
        #             # Get the contrast at that location
        #             sep = sqrt(pos[1]^2 + pos[2]^2)
        #             σ = contrasts[I](sep / platescale)
        #             # Fallback if we fall off the edge of the image
        #             if isnan(σ)
        #                 σ = 1e1
        #             end
        #             # The liklihood function from JB
        #             # ln_post += 1/2σ * (phot^2 - 2phot*phot_img)
        
        #             # Seems to be negative?
        #             ln_post_add  = -1/2σ * (phot^2 - 2phot*phot_img)
        #             # ln_post_add  = 1/2σ * (phot^2 - 2phot*phot_img)
        
        #             # Me trying something
        #             if isfinite(ln_post_add)
        #                 ln_post += ln_post_add
        #             end
        #             # ln_post += -1/2σ * (phot^2 - 2phot*phot_img)
        
        #             # ln_post += log(
        #             #         1 / √(2π)σ * exp(
        #             #             -1 / 2((phot - phot_img)^2 / σ^2)
        #             #         )
        #             # )
        #         end
        #         # Fallback to a finite but bad value if we fall off the edge of the image
        #         # Is there a mathematically better way to express that we don't have this
        #         # information?
        #         if !isfinite(ln_post)
        #             @warn "non-finite posterior" maxlog = 20
        #             return Inf
        #         end
end



function make_ln_post_images(priors, static, images, contrasts, times, platescale)
    
    ln_prior = make_ln_prior(priors)
    ln_like = make_ln_like_images(keys(priors), static, images, contrasts, times, platescale)

    ln_post(params) = ln_prior(params) + ln_like(params)

    return ln_post
end

function fit_images(
    priors, static, images, contrasts, times;
    platescale,
    numwalkers=10,
    burnin = 10_000,
    thinning = 1,
    numsamples_perwalker = 10_000
    )

    # Initial values for the walkers are drawn from the priors
    initial_walkers = reduce(hcat, [rand.([priors...,]) for _ in 1:numwalkers])

    ln_post = make_ln_post_images(priors, static, images, contrasts, times, platescale)

    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, numwalkers, initial_walkers, burnin, 1);
    @info "Burn in done"
    chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, size(chain,2), chain[:, :, end], numsamples_perwalker, thinning);
    @info "Chains done"

    column_names = string.(collect(keys(priors)))
    chain = Chains(permutedims(chain, (3, 1, 2)),  column_names);

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
    sizehint!(out, size(sampled,1)*size(sampled,2))
    for i in 1:size(sampled,1), j in 1:size(sampled,2)
        nt = (;(k=>v for (k,v) in zip(keys(sampled), Array(sampled[i,:,j])))...)
        el = KeplerianElements(merge(nt, static))
        push!(out, el)
    end
    return out[begin:min(N,end)]
end


##
# 1


# module Fitting

# using ..DirectOrbits
# using DirectImages

# using Distributions
# using StaticArrays

# function planet_ln_like(convolved, priors, μ, plx, mjds, contrasts, platescale)

#     # Samplers normally must supply their arguments as a vector
#     function lnlike(phot, a, i, e, τ, ω, Ω)
#         return lnlike((phot, a, i, e, τ, ω, Ω))
#     end
    
#     function lnlike(params)
#         (phot, a, i, e, τ, ω, Ω) = params
#         # The prior is given by the input distributions.
#         # Sum their log pdf at this location
#         ln_prior = zero(phot)
#         for i in eachindex(params)
#             pd = priors[i]
#             param = params[i]
#             ln_prior += logpdf(pd, param)
#         end

#         # return ln_prior

#         # Construct an orbit with these parameters
#         orbit = Orbit(a, i, e, τ, μ, ω, Ω, plx)
#         # Then get the liklihood by multipliying together
#         # the liklihoods at each epoch (summing the logs)
#         ln_post = 0.0
#         # Go through each image
#         for I in eachindex(convolved)
#             # Find the current position in arcseconds
#             pos = SVector(0., 0., 0.)
#             try
#                 pos = xyz(orbit, mjds[I])
#             catch e
#                 if !(typeof(e) <: Real)
#                     @warn "error looking up orbit" exception = e maxlog = 2
#                 end
#             end
#             # Get the value in the convolved images
#             pos = SVector(pos[1], -pos[2])
#             # I believe that the x-coordinate should be flipped because
#             # image indices are opposite to sky coordinates

#             # pos = SVector(-pos[2], pos[1])
#             # pos = SVector(pos[2], -pos[1])

#             phot_img = lookup_coord(convolved[I], pos, platescale)

#             # NOTE: experimenting with flipped coords!
#             # pos# .* SVector(1, 1, 1)
#             # pos = SVector(-pos[2], pos[1])
#             # phot_img = lookup_coord(convolved[I], pos, platescale)
#             # Get the contrast at that location
#             sep = sqrt(pos[1]^2 + pos[2]^2)
#             σ = contrasts[I](sep / platescale)
#             # Fallback if we fall off the edge of the image
#             if isnan(σ)
#                 σ = 1e1
#             end
#             # The liklihood function from JB
#             # ln_post += 1/2σ * (phot^2 - 2phot*phot_img)

#             # Seems to be negative?
#             ln_post_add  = -1/2σ * (phot^2 - 2phot*phot_img)
#             # ln_post_add  = 1/2σ * (phot^2 - 2phot*phot_img)

#             # Me trying something
#             if isfinite(ln_post_add)
#                 ln_post += ln_post_add
#             end
#             # ln_post += -1/2σ * (phot^2 - 2phot*phot_img)

#             # ln_post += log(
#             #         1 / √(2π)σ * exp(
#             #             -1 / 2((phot - phot_img)^2 / σ^2)
#             #         )
#             # )
#         end
#         # Fallback to a finite but bad value if we fall off the edge of the image
#         # Is there a mathematically better way to express that we don't have this
#         # information?
#         if !isfinite(ln_post)
#             @warn "non-finite posterior" maxlog = 20
#             return Inf
#         end
#         # Multiply the liklihood by the prior
#         return ln_post + ln_prior
#         # return ln_prior
#         # return ln_post

#     end
#     return lnlike
# end

# end