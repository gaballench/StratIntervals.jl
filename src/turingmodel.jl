# wrapper for the the Turing model for a single tratigraphic interval
function sample_stratinterval(data_priors::StratInterval, iters, ϵ, τ_hmc, prior)
    # then declare the Turing model as a hierarchical Bayesian model with
    # priors on theta1, theta2, and lambda
    @model function stratint(data_priors::StratInterval)
        # convert to real
        #data = convert(Vector{Real}, data_priors.data)
        # priors
        θ1 ~ data_priors.prior1
        θ2 ~ data_priors.prior2
        λ ~ data_priors.prior3
        # reject values outside the support of the parameters
        if θ1 > minimum(data_priors.data) || θ2 < maximum(data_priors.data) || θ1 < 0 || θ2 < 0
            Turing.@addlogprob! -Inf
            return nothing
        end
        # define the likelihood of the observed data
        N = length(data_priors.data)
        for i in 1:N
            τ[i] ~ BetaAdaptive(θ1, θ2, λ)
            #println(τ[i]~ BetaAdaptive(θ1, θ2, λ))
        end
    end
    # sample from prior? enter block and return chain
    if prior
        chain = sample(stratint(data_priors.data, θ1, θ2, λ), Prior(), HMC(ϵ, τ_hmc), iters)
        return chain
    end
    # else, sample from posterior and return chain
    chain = sample(stratint(data_priors.data, θ1, θ2, λ), HMC(ϵ, τ_hmc), iters)
    return chain
end

# sampling from a collection needs to be carried out by defining the likelihood by a conflation of the collection of distributions, but said conflation is univariate and we need to solve the issue of having multiple data vectors one for each distribution. Maybe after we sample, we calculate the conflation as a deterministic node y := conflation(d, x)
