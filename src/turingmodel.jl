# wrapper for the the Turing model for a single tratigraphic interval
function sample_stratinterval(data_priors::StratInterval, iters, ϵ, τ_hmc, prior, postpredict)
    # dismantle the StratInterval object
    data = data_priors.data
    θ1_prior = data_priors.θ1_prior
    θ2_prior = data_priors.θ2_prior
    λ_prior = data_priors.λ_prior    
    # then declare the Turing model as a hierarchical Bayesian model with
    # priors on theta1, theta2, and lambda
    @model function stratint(τ, ::Type{T}=Float64) where {T}
        # initialise tau if missing for posterior predictive
        if τ === missing
            # Initialize `x` if missing
            τ = Vector{T}(undef, 1)
        end
        # priors are Float64 if fixed, Distribution otherwise
        # not all the priors can be fixed!
        if sum(map(x -> x isa Float64, (θ1_prior, θ2_prior, λ_prior))) == 3
            error("Not all the parameters can be fixed, at least one needs to be sampled. All parameters were Float64")
        end
        if θ1_prior isa Float64
            θ1 = θ1_prior
        else
            θ1 ~ θ1_prior
        end
        if θ2_prior isa Float64
            θ2 = θ2_prior
        else
            θ2 ~ θ2_prior
        end
        if λ_prior isa Float64
            λ = λ_prior
        else
            λ ~ λ_prior
        end
        # reject values outside the support of the parameters
        if θ1 > minimum(τ) || θ2 < maximum(τ) || θ1 < 0 || θ2 < 0
            Turing.@addlogprob! -Inf
            return nothing
        end
        # define the likelihood of the observed data
        N = length(τ)
        for i in 1:N
            τ[i] ~ BetaAdaptive(θ1, θ2, λ)
            #println(τ[i]~ BetaAdaptive(θ1, θ2, λ))
        end
    end
    # sample from prior? enter block and return chain
    if prior
        chain = sample(stratint(convert(Vector{Real}, data)),
                       Prior(),
                       iters)
        return chain
    end
    # else, sample from posterior and return chain
    chain = sample(stratint(convert(Vector{Real}, data)),
                   HMC(ϵ, τ_hmc),
                   iters)
    # if posterior predictive, sample from chain and using the model
    if postpredict
        pars = DataFrames.DataFrame(chain)[:, [:θ1, :θ2, :λ]]        
        # generate the posterior predictive from the samples in `chain`
        # by mapping a rand call to BetaAdaptive with each collection of values 
        # prepare the results vector
        postpredict = Vector(undef, iters)        
        # map over iterations to sample one value from each combination of params
        postpredict = map(x -> rand(BetaAdaptive(pars.θ1[x], pars.θ2[x], pars.λ[x]), 1)[1], 1:iters)
        return (chain, postpredict)
    else
        return chain
    end
end

# this method is just doing sequential estimation by looping over the vector of StratIntervals 
function sample_stratinterval(stratintervals::Vector{StratInterval}, iters, ϵ, τ_hmc, prior, postpredict)
    output = Vector(undef, length(stratintervals))
    counter = 1
    for i in stratintervals
        output[counter] = sample_stratinterval(i, iters, ϵ, τ_hmc, prior, postpredict)
        counter += 1
    end
    return output
end
