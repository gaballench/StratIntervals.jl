# wrapper for the the Turing model for a single tratigraphic interval
function sample_stratinterval(data_priors::StratInterval, iters, ϵ, τ_hmc, prior)
    # dismantle the StratInterval object
    data = data_priors.data
    θ1_prior = data_priors.θ1_prior
    θ2_prior = data_priors.θ2_prior
    λ_prior = data_priors.λ_prior    
    # then declare the Turing model as a hierarchical Bayesian model with
    # priors on theta1, theta2, and lambda
    @model function stratint(τ::Vector{Real})
        # priors
        θ1 ~ θ1_prior
        θ2 ~ θ2_prior
        λ ~ λ_prior
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
        chain = sample(stratint(data),
                       Prior(),
                       iters)
        return chain
    end
    # else, sample from posterior and return chain
    chain = sample(stratint(data),
                   HMC(ϵ, τ_hmc),
                   iters)
    return chain
end
