
# struct for the input data into the model
# data are being used by the loglik inside the Turing model
#*_prior are passed on to the parameter priors inside the Turing model 
struct StratInterval
    data::Vector{Real}
    θ1_prior::Union{ContinuousUnivariateDistribution,Float64}
    θ2_prior::Union{ContinuousUnivariateDistribution,Float64}
    λ_prior::Union{ContinuousUnivariateDistribution,Float64}
end
