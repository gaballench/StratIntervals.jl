
# struct for the input data into the model
# data are being used by the loglik inside the Turing model
#*_prior are passed on to the parameter priors inside the Turing model 
struct StratInterval
    data::Vector{Float64}
    theta1_prior::ContinuousUnivariateDistribution
    theta2_prior::ContinuousUnivariateDistribution
    lambda_prior::ContinuousUnivariateDistribution
end
