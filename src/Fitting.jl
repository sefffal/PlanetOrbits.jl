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

# Here is the same exact code, but manually unrolled.
# This does limit the max number of priors we can support arbitrarily but is ~7x faster
# function make_ln_prior(p1, p2=nothing, p3=nothing, p4=nothing, p5=nothing, p6=nothing, p7=nothing, p8=nothing, p9=nothing)
#     function ln_prior(params)
#         lp = zero(first(params))
#         lp += logpdf(p1, params[1])
#         if !isnothing(p2)
#             lp += logpdf(p2, params[2])
#         end
#         if !isnothing(p3)
#             lp += logpdf(p3, params[3])
#         end
#         if !isnothing(p4)
#             lp += logpdf(p4, params[4])
#         end
#         if !isnothing(p5)
#             lp += logpdf(p5, params[5])
#         end
#         if !isnothing(p6)
#             lp += logpdf(p6, params[6])
#         end
#         if !isnothing(p7)
#             lp += logpdf(p7, params[7])
#         end
#         if !isnothing(p8)
#             lp += logpdf(p8, params[8])
#         end
#         if !isnothing(p9)
#             lp += logpdf(p9, params[9])
#         end
#         return lp
#     end
#     return ln_prior
# end
# hierach_logpdf(p::Any, val) = logpdf(p, val)
# hierach_logpdf(ps::Tuple, vals) = sum(hierach_logpdf.(ps, vals))
# hierach_logpdf(ps::AbstractArray, vals) = sum(hierach_logpdf.(ps, vals))
# function make_ln_prior(p1)
#     function ln_prior(params)
#         return hierach_logpdf(p1, params[1])
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4, p5)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         lp += hierach_logpdf(p5, params[5])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4, p5, p6)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         lp += hierach_logpdf(p5, params[5])
#         lp += hierach_logpdf(p6, params[6])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4, p5, p6, p7)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         lp += hierach_logpdf(p5, params[5])
#         lp += hierach_logpdf(p6, params[6])
#         lp += hierach_logpdf(p7, params[7])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4, p5, p6, p7, p8)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         lp += hierach_logpdf(p5, params[5])
#         lp += hierach_logpdf(p6, params[6])
#         lp += hierach_logpdf(p7, params[7])
#         lp += hierach_logpdf(p8, params[8])
#         return lp
#     end
#     return ln_prior
# end
# function make_ln_prior(p1, p2, p3, p4, p5, p6, p7, p8, p9)
#     function ln_prior(params)
#         lp = hierach_logpdf(p1, params[1])
#         lp += hierach_logpdf(p2, params[2])
#         lp += hierach_logpdf(p3, params[3])
#         lp += hierach_logpdf(p4, params[4])
#         lp += hierach_logpdf(p5, params[5])
#         lp += hierach_logpdf(p6, params[6])
#         lp += hierach_logpdf(p7, params[7])
#         lp += hierach_logpdf(p8, params[8])
#         lp += hierach_logpdf(p9, params[9])
#         return lp
#     end
#     return ln_prior
# end


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


# TODO: more than one layer of nesting?
# TODO: This gets slow with more parameters; it hit's some threshold
@inline function hierach_merge(merged, i)
    if typeof(merged[1]) <: Number
        T = typeof(merged[1])
    else
        T = typeof(merged[1][1])
    end
    return NamedTuple{keys(merged),NTuple{length(keys(merged)),Float64}}(
        (length(vs) > 1 ? vs[i+1] * first(vs) : vs for vs in values(merged))
    )
end

function make_ln_like_images(priors, images, contrasts, times, platescale)
    # props is tuple of symbols
    # static is named tuple of values


    if !(all(size(images)  == size(times) == size(contrasts) == size(planet.epochs) for planet in priors.planets))
        error("All values must have the same length")
    end

    # Return our
    return function ln_like(θ)

        # The ln likelihood:
        ll = 0.0
        for (priors_planet, θ_planet) in zip(priors.planets, θ.planets)
            
            for i in eachindex(priors_planet.epochs)


                # TODO: this merging needs to be worked out at compile time, or at least when building the function!
                # TODO: see ComponentArrays.label2index
                # We already know the layout, can even just look up by index.                # Merge the three levels together. This gives us the deepest nested value for any given variable.
                θ_planet_epoch = θ_planet.epochs[i]
                # @time elements = KeplerianElements(merge(NamedTuple(θ), NamedTuple(θ_planet), NamedTuple(θ_planet_epoch)))
                elements = KeplerianElements((;θ.μ, θ.plx, θ_planet.i, θ_planet.Ω, θ_planet.ω, θ_planet.e, θ_planet.τ, θ_planet.a))
                f = θ_planet_epoch.f

                image = images[i]
                contrast = contrasts[i]
                t = times[i]
                ra, dec = kep2cart(elements, t)
                x = -ra
                y = dec

                # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
                f̃ₓ = lookup_coord(image, x,y, platescale)
                

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

                ll += l
            end

            # Spread in flux between epochs
            # ll += 
            # *exp(-0.5*(K0/L0 - Model_ratio)^2/model_spread^2)
        end

        # At this point, a NaN or Inf log-likelihood implies
        # an error in preparing the inputs or in this code.
        if !isfinite(ll)
            @show x y r σₓ f f̃ₓ l
            error("Non-finite log-likelihood encountered")
        end
        return ll
    end
end



function make_ln_post_images(priors, images, contrasts, times, platescale)
    ln_prior = make_ln_prior(priors)
    ln_like = make_ln_like_images(priors, images, contrasts, times, platescale)
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
    priors, images, contrasts, times;
    platescale,
    burnin,
    numwalkers=10,
    thinning = 1,
    numsamples_perwalker,
    squash=true
    )
    # column_names = string.(collect(keys(priors)))
    column_names = ComponentArrays.labels(priors)

    ln_post = make_ln_post_images(priors, images, contrasts, times, platescale)

    
    @info "Finding starting point"
    # TODO: kissmcmc has a better method for creating the ball and rejecting some starting points
    initial_walkers = find_starting_walkers(ln_post, priors, numwalkers)

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(;NamedTuple(cv)...)
        for cv in initial_walkers
    ]
    # initial_walkers = SVector{length(priors),Float64}.(initial_walkers)

    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin=burnin*numwalkers,
        use_progress_meter=true,
        nthin=thinning,
        niter=numsamples_perwalker*numwalkers
    )

    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(thetase′), thetase′)
    else
        # We can reinterpret the vector of SVectors as a matrix directly without copying!
        # This can save massive amounts of memory and time on large changes
        reinterptted = cat(
            [transpose(reinterpret(reshape, eltype(first(θ)), θ)) for θ in thetase]...,
            dims=3
        )
    end

    return Chains(reinterptted, column_names)
end

function find_starting_point(ln_post, priors)
    initial_guess = rand.(priors)
    i = 0
    while !isfinite(ln_post(initial_guess))
        i+=1
        initial_guess = rand.(priors)
        if i > 1000
            error("Could not find a starting point in the posterior that is finite by drawing from the priors after 1000 attempts")
        end
    end
    # return initial_guess


    objective(params) = -ln_post(params)

    # p = TwiceDifferentiable(objective, Float64.(collect(initial_guess)), autodiff=:forward)
    # results = Optim.optimize(objective, initial_guess, LBFGS(), Optim.Options(iterations=5000, f_tol=1e-9))
    results = Optim.optimize(objective, initial_guess, NelderMead(), Optim.Options(show_trace=false,iterations=5000, f_tol=1e-9))
    return results.minimizer
end

# # Function to get a maximum likelihood position to start the sampler from
# function find_starting_point(ln_post, priors)
#     initial_guess = mode.(collect(priors))
#     for i in eachindex(initial_guess)
#         if initial_guess[i] ==0
#             initial_guess[i] += 0.1rand()
#         end
#     end
#     objective(params) = -ln_post(params)

#     # p = TwiceDifferentiable(objective, Float64.(collect(initial_guess)), autodiff=:forward)
#     # results = Optim.optimize(objective, initial_guess, LBFGS(), Optim.Options(iterations=5000, f_tol=1e-9))
#     results = Optim.optimize(objective, initial_guess, NelderMead(), Optim.Options(show_trace=false,iterations=5000, f_tol=1e-9))
#     return results.minimizer
# end


# Start walkers in a gaussian ball around the MLE, while ensuring we don't
# step outside the ranges defined by the priors
function find_starting_walkers(ln_post, priors, numwalkers)
    # initial_walkers = mapreduce(hcat, 1:numwalkers) do _
    initial_walkers = map(1:numwalkers) do _
        initial_position = find_starting_point(ln_post, priors)
        # This used to intiialize the walkers in a Gaussian ball around the MAP.
        # But now we just draw starting points randomly from the priors, this isn't needed.
        @showprogress "Finding initial positions" map(eachindex(initial_position)) do i
            p = NaN
            while !(minimum(priors[i]) < p < maximum(priors[i]))
                p = (1+randn()) * initial_position[i]
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
