# wrapper for the the Turing model for a single tratigraphic interval
"""
    sample_stratinterval(data_priors::StratInterval, iters, sampler, prior, postpredict)

This function samples the posterior distributions of the stratigraphic interval using MCMC.

Input are the data and priors as constructed using StratInterval(). iters is the number of iterations during sampling,
sampler is the kind of sampler to use, e.g., HMC(epsilon, tau) or NUTS(), whether to sample from the prior,
and finally the last argument specifies whether we want to calculate the posterior predictive distribution of τ.

# examples

```@example
using Distributions
using Turing
using StratIntervals
my_strinterval = StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], Normal(10, 2), Exponential(1), Normal(0, 1))
sample_stratinterval(my_strinterval, 10000, NUTS(), true, false) # sample from prior
sample_stratinterval(my_strinterval, 10000, NUTS(), false, false) # sample from posterior
sample_stratinterval(my_strinterval, 10000, NUTS(), false, true) # sample from posterior and calculate posterior predictive of tau
```
"""
function sample_stratinterval(data_priors::StratInterval, iters, sampler, prior, postpredict)
    # dismantle the StratInterval object
    # if data is a vector of distributions, create the data vector by sampling one rand from each distribution
    if prod(map(x -> x isa T where T <: ContinuousUnivariateDistribution, data_priors.data))
        data = map(x -> rand(x, 1)[1], data_priors.data)
        data = convert(Vector{Real}, data)
    else
        # else just get the data from StratInterval
        data = data_priors.data
    end
    # then collect the priors    
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
        if θ2 > minimum(τ) || θ1 < maximum(τ) || θ1 < 0 || θ2 < 0
            Turing.@addlogprob! -Inf
            return nothing
        end
        # define the likelihood of the observed data
        N = length(τ)
        for i in 1:N
            τ[i] ~ ThreeParBeta(θ1, θ2, λ)
            #println(τ[i]~ ThreeParBeta(θ1, θ2, λ))
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
                   sampler,
                   iters)
    # if posterior predictive, sample from chain and using the model
    if postpredict
        # block here for repeating the fixed param value
        # during posterior predictive sampling
        if θ1_prior isa Float64
            θ1_trace = vcat(fill.(θ1_prior, iters))
        else
            θ1_trace = DataFrames.DataFrame(chain)[:, :θ1]
        end
        if θ2_prior isa Float64
            θ2_trace = vcat(fill.(θ2_prior, iters))
        else
            θ2_trace = DataFrames.DataFrame(chain)[:, :θ2]
        end
        if λ_prior isa Float64
            λ_trace = vcat(fill.(λ_prior, iters))
        else
            λ_trace = DataFrames.DataFrame(chain)[:, :λ]
        end
        # generate the posterior predictive from the samples in `chain`
        # by mapping a rand call to ThreeParBeta with each collection of values 
        # prepare the results vector
        postpredict = Vector(undef, iters)        
        # map over iterations to sample one value from each combination of params
        postpredict = map(x -> rand(ThreeParBeta(θ1_trace[x],
                                                 θ2_trace[x],
                                                 λ_trace[x]), 1)[1], 1:iters)
        return (chain, postpredict)
    else
        return chain
    end
end

# this method is just doing sequential estimation by looping over the vector of StratIntervals 
"""
    sample_stratinterval(data_priors::Vector{StratInterval}, iters, sampler, prior, postpredict)

This function samples the posterior distributions of a vector of  stratigraphic intervals using MCMC.

Input data is a vector of type StratInterval. iters is the number of iterations during sampling,
sampler is the kind of sampler to use, e.g., HMC(epsilon, tau) or NUTS(), whether to sample from the prior,
and finally the last argument specifies whether we want to calculate the posterior predictive distribution of τ.

# examples

```@example
using Distributions
using Turing
ndata = 100
iters = 10000
# first interval
true_lambda_1 = 0
true_theta2_1 = 10
true_theta1_1 = 22
data_1 = rand(ThreeParBeta(true_theta1_1, true_theta2_1, true_lambda_1), ndata)
# second interval
true_lambda_2 = 1
true_theta2_2 = 15
true_theta1_2 = 30
data_2 = rand(ThreeParBeta(true_theta1_2, true_theta2_2, true_lambda_2), ndata)
interval_1 = StratInterval(data_1, Normal(22, 3), Exponential(10), Normal(0,3))
interval_2 = StratInterval(data_2, Normal(30, 3), Exponential(15), Normal(1,3))
# construct the vector of StratIntervals
vecinterval = [interval_1, interval_2]
# MCMC sampling
mystratint_mcmc_vec = sample_stratinterval(vecinterval, iters, NUTS(), false, false)
# MCMC sampling and posterior predictive
mystratint_postpredict_vec = sample_stratinterval(vecinterval, iters, NUTS(), false, true)
```
"""
function sample_stratinterval(stratintervals::Vector{StratInterval}, iters, sampler, prior, postpredict)
    output = Vector(undef, length(stratintervals))
    counter = 1
    for i in stratintervals
        output[counter] = sample_stratinterval(i, iters, sampler, prior, postpredict)
        counter += 1
    end
    return output
end
