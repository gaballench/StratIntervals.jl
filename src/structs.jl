
# struct for the input data into the model
# data are being used by the loglik inside the Turing model
#*_prior are passed on to the parameter priors inside the Turing model 
struct StratInterval
    data::Vector{Real}
    θ1_prior::ContinuousUnivariateDistribution
    θ2_prior::ContinuousUnivariateDistribution
    λ_prior::ContinuousUnivariateDistribution
end
