
# struct for the input data into the model
# data are being used by the loglik inside the Turing model
#*_prior are passed on to the parameter priors inside the Turing model 
"""
    StratInterval(data, θ1_prior, θ2_prior, λ_prior)

This structure allows to input the data for estimation of stratigraphic intervals.
Four attributes are available: data, a Vector{Real} with the time occurrences,
θ1_prior, a ContinuousUnivariateDistribution, specifying the prior on θ1,
θ2_prior, a ContinuousUnivariateDistribution, specifying the prior on θ2, and
λ_prior, a ContinuousUnivariateDistribution, specifying the prior on λ.

When sampling all these parameters, the priors need to be distributions. However,
it is possible to fix parameters to concrete values by specifying a Float64 instead of a 
ContinuousUnivariateDistribution. Care must be taken to avoid fixing all of the parameters,
at least one needs to be a ContinuousUnivariateDistribution.

Fixing a prior implies fixing its value in the model, so this needs to be used with extreme caution.

Fixing a prior essentially means calculating the posterior _conditional_ to that observed value,
and therefore does not correspond to the whole posterior. This is useful for simulation but
is discouraged in empirical applications, no matter how certain we seem to be about the value for
a given parameter.

An important exception to the rule above is when we have stratigraphic data for an extant lineage,
in this case we can (and in fact should) fix θ1_prior to 0.0

# examples

```@example
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], Normal(10, 2), Exponential(1), Normal(0, 1)) # use priors for all parameters
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], Normal(10, 2), 1.5, Normal(0, 1)) # fix θ2_prior to 1.5
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], 9.8, 1.5, Normal(0, 1)) # fix θ2_prior to 1.5 and θ1_prior to 9.8
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], 9.8, 1.5, 0.0) # this triggers an error
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], Normal(10, 2), Exponential(1), 0.0) # fix λ_prior to 0.0
StratInterval([2.0, 3.1, 3.2, 4.6, 6.77], Normal(10, 2), 0.0, Normal(0, 1)) # fix θ2_prior to 0.0 because the lineage is extant
```
"""
struct StratInterval
    data::Vector{Union{ContinuousUnivariateDistribution,Real}}
    θ1_prior::Union{ContinuousUnivariateDistribution,Float64}
    θ2_prior::Union{ContinuousUnivariateDistribution,Float64}
    λ_prior::Union{ContinuousUnivariateDistribution,Float64}
end
