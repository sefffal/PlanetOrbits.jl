using ComponentArrays
using Optim
using Distributions: mode, logpdf

import KissMCMC

using MCMCChains
using NamedTupleTools
using DirectImages: lookup_coord
using Base.Threads: @threads
import Random

using ThreadPools
using ProgressMeter
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
    # initial_ca = ComponentArray{Float64}(initial);
    initial_elements = KeplerianElements(merge(initial, static))

    props = keys(initial)
    objective = let static=static, points=points, times=times
        function (params)
            nt = namedtuple(props, params)
            el = KeplerianElements(merge(nt, static))
            dist =  DirectOrbits.least_squares_distance(el, points, times)
            return dist
        end
    end

    p = TwiceDifferentiable(objective, Float64.(collect(initial)), autodiff=:finite)
    # results = Optim.optimize(p, Float64.(collect(initial)), LBFGS(), Optim.Options(store_trace=trace, extended_trace=trace, f_tol=1e-6))
    results = Optim.optimize(p, Float64.(collect(initial)), NelderMead(), Optim.Options(store_trace=trace, extended_trace=trace, iterations=5000))

    # results = Optim.optimize(objective, initial_ca, NelderMead(), Optim.Options(store_trace=trace, extended_trace=trace))
    optimized_elements =  KeplerianElements(merge(namedtuple(props, Optim.minimizer(results)), static))

    # Gather elements from trace
    if trace
        elements_trace = map(results.trace) do trace_step
            if haskey(trace_step.metadata, "centroid")
                KeplerianElements(merge(namedtuple(props,trace_step.metadata["centroid"]), static))
            else
                KeplerianElements(merge(namedtuple(props,trace_step.metadata["x"]), static))
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
    burnin,
    thinning = 1,
    numsamples_perwalker
    )

    # Initial values for the walkers are drawn from the priors
    initial_walkers = reduce(hcat, [rand.([priors...,]) for _ in 1:numwalkers])

    ln_post = make_ln_post_astrom(priors, static, points, times, uncertainties)

    error("Not implemented -- adapt to kissmcmc instead of AffineInvariantMCMC").

    # chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, numwalkers, initial_walkers, burnin, 1);
    # @info "Burn in done"
    # chain, llhoodvals = AffineInvariantMCMC.sample(ln_post, size(chain,2), chain[:, :, end], numsamples_perwalker, thinning);
    # @info "Chains done"

    column_names = string.(collect(keys(priors)))
    chain = Chains(permutedims(chain, (3, 1, 2)),  column_names);

    return chain
end


struct TupleProtoHolder{T<:Any}
end


function make_ln_like_images(props, static, images, contrasts, times, platescale)
    # props is tuple of symbols
    # static is named tuple of values

    if !(size(images)  == size(times) == size(contrasts))
        error("All values must have the same length")
    end

    # template_holder = Val{props}()
    num_type = Float64
    template_holder = Val{NamedTuple{props,NTuple{length(props),num_type}}}()

    function function_barrier(template::Val{T}) where T
    # function function_barrier(template)# where T
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
                    # f̃ₓ = image[ix,iy]
                    f̃ₓ = lookup_coord(image, (x,y), platescale)
                else
                    continue
                end
                

                r = √(x^2 + y^2)
                σₓ = contrast(r/platescale)

                # When we get a position that falls outside of our available
                # data (e.g. under the coronagraph) we cannot say anything
                # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
                # of zero.
                if !isfinite(σₓ) || !isfinite(f̃ₓ)
                    continue
                end

                # Ruffio et al 2017, eqn 31
                l = -1/(2σₓ^2) * (f^2 - 2f*f̃ₓ)

                # At this point, a NaN or Inf log-likelihood implies
                # an error in preparing the inputs or in this code.
                if !isfinite(l)
                    @show x y r σₓ f f̃ₓ l
                    error("Non-finite log-likelihood encountered")
                end

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



function fit_images_lsq(
    priors,
    static, images, contrasts, times;
    platescale,
    )

    # Initial values for the walkers are drawn from the priors
    initial = rand.([priors...,])
    initial_elements = KeplerianElements(merge(initial, static))

    props = keys(initial)
    objective = let static=static, points=points, times=times
        function (params)
            nt = namedtuple(props, params)
            el = KeplerianElements(merge(nt, static))
            dist =  DirectOrbits.least_squares_distance(el, points, times)
            return dist
        end
    end

    results = Optim.optimize(p, Float64.(collect(initial)), NelderMead(), Optim.Options(store_trace=trace, extended_trace=trace, iterations=5000))

    # results = Optim.optimize(objective, initial_ca, NelderMead(), Optim.Options(store_trace=trace, extended_trace=trace))
    optimized_elements =  KeplerianElements(merge(namedtuple(props, Optim.minimizer(results)), static))

    return optimized_elements
end

function fit_images_kissmcmc(
    priors, static, images, contrasts, times;
    platescale,
    burnin,
    numwalkers=10,
    thinning = 1,
    numsamples_perwalker,
    squash=true
    )
    column_names = string.(collect(keys(priors)))

    ln_post = make_ln_post_images(priors, static, images, contrasts, times, platescale)

    
    @info "Finding starting point"
    # TODO: kissmcmc has a better method for creating the ball and rejecting some starting points
    initial_walkers = collect.(eachcol(find_starting_walkers(ln_post, priors, numwalkers)))
    initial_walkers = SVector{length(priors),Float64}.(initial_walkers)

    @time thetase, _accept_ratioe = KissMCMC.emcee(ln_post, initial_walkers; nburnin=burnin*numwalkers, use_progress_meter=true, nthin=thinning, niter=numsamples_perwalker*numwalkers);

    if squash
        @time thetase, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
    end

    # We can reinterpret the vector of SVectors as a matrix directly without copying!
    # This can save massive amounts of memory and time on large changes
    @time reinterptted = cat(
        [transpose(reinterpret(reshape, eltype(first(θ)), θ)) for θ in thetase]...,
        dims=3
    )
    return Chains(reinterptted, column_names)
end


# Function to get a maximum likelihood position to start the sampler from
function find_starting_point(ln_post, priors)
    initial_guess = mode.(collect(priors))
    for i in eachindex(initial_guess)
        if initial_guess[i] ==0
            initial_guess[i] += 0.1rand()
        end
    end
    objective(params) = -ln_post(params)

    # p = TwiceDifferentiable(objective, Float64.(collect(initial_guess)), autodiff=:forward)
    # results = Optim.optimize(objective, initial_guess, LBFGS(), Optim.Options(iterations=5000, f_tol=1e-9))
    results = Optim.optimize(objective, initial_guess, NelderMead(), Optim.Options(show_trace=false,iterations=5000, f_tol=1e-9))
    return results.minimizer
end


# Start walkers in a gaussian ball around the MLE, while ensuring we don't
# step outside the ranges defined by the priors
function find_starting_walkers(ln_post, priors, numwalkers)
    initial_walkers = mapreduce(hcat, 1:numwalkers) do _
        initial_position = find_starting_point(ln_post, priors)
        @showprogress "Finding initial positions" map(eachindex(initial_position)) do i
            p = NaN
            while !(minimum(priors[i]) < p < maximum(priors[i]))
                p = initial_position[i] + 0.01randn()*initial_position[i]
            end
            p
        end
    end
    return initial_walkers
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


function chain2elements(chains::Chains, static,)
    out = KeplerianElements{Float64}[]

    proto = namedtuple(keys(chains))

    sizehint!(out, size(chains,1)*size(chains,2))
    for i in 1:size(chains,1), j in 1:size(chains,3)
        nt = proto(Array(chains[i,:,j]))
        el = KeplerianElements(merge(nt, static))
        push!(out, el)
    end
    return out
end
