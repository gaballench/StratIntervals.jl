# wrapper for the the Turing model
function sample_stratinterval(data_priors::StratInterval, iters, epsilon, tau_hmc, prior)
    # then declare the Turing model as a hierarchical Bayesian model with
    # priors on theta1, theta2, and lambda
    @model function stratint(data_priors::StratInterval)
        # convert to real
        #data = convert(Vector{Real}, data_priors.data)
        # priors
        theta1 ~ data_priors.prior1
        theta2 ~ data_priors.prior2
        lambda ~ data_priors.prior3
        # reject values outside the support of the parameters
        if theta1 > minimum(data_priors.data) || theta2 < maximum(data_priors.data) || theta1 < 0 || theta2 < 0
            Turing.@addlogprob! -Inf
            return nothing
        end
        # define the likelihood of the observed data
        N = length(data_priors.data)
        for i in 1:N
            tau[i] ~ BetaAdaptive(theta1, theta2, lambda)
            #println(tau[i]~ BetaAdaptive(theta1, theta2, lambda))
        end
    end
    # sample from prior? enter block and return chain
    if prior
        chain = sample(stratint(data_priors.data, theta1, theta2, lambda), Prior(), HMC(epsilon, tau_hmc), iters)
        return chain
    end
    # else, sample from posterior and return chain
    chain = sample(stratint(data_priors.data, theta1, theta2, lambda), HMC(epsilon, tau_hmc), iters)
    
    return chain
end
