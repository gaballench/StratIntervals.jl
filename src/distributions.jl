#using Distributions
#using SpecialFunctions # for the beta function gamma(a)gamma(b)/gamma(a+b) 
#using Base
#using Random

function fourpar_dbeta(τ,α,β,τ1,τ2) 
    f = ((τ-τ1)^(α - 1) * (τ2-τ)^(β - 1)) / ((τ2-τ1)^(α+β-1) * SpecialFunctions.beta(α,β))
    return(f)
end

### functions with reparam
function threepar_dbeta(τ,τ1,τ2,λ)
    if (λ <= 0)
        #println("Calculating under λ <= 0")
        f = ((τ2 - τ)^(-λ)) / ((τ2 - τ1)^(1-λ) * SpecialFunctions.beta(1,1-λ))
        #f = ((τ2 - τ)^(-λ)) / ((τ2 - τ1)^(1-λ) * Bnegative)
    else
        #println("Calculating under λ > 0")
        f = ((τ-τ1)^(λ)) / ((τ2 - τ1)^(1+λ) * SpecialFunctions.beta(1+λ,1))
        #f = ((τ-τ1)^(λ)) / ((τ2 - τ1)^(1+λ) * Bpositive)
    end
    return(f)    
end

####################################
### Implementing the distribution in Julia
####################################

### define a subtype
# this probably needs to be defined as imutable with the parameters θ. and λ
struct BetaAdaptive <: ContinuousUnivariateDistribution
    θ1::Real
    θ2::Real
    λ::Real
end

### declaring methods in a quick way, check whether this is correct
### define the pdf
function Distributions.pdf(d::BetaAdaptive, τ::Real)
    #threepar_dbeta(d.τ, d.θ1, d.θ2, d.λ)
    if (d.λ <= 0)
        #println("Calculating under λ <= 0")
        f = ((d.θ2 - τ)^(-d.λ)) / ((d.θ2 - d.θ1)^(1-d.λ) * SpecialFunctions.beta(1,1-d.λ))
        #f = ((t2 - τ)^(-λ)) / ((t2 - t1)^(1-λ) * Bnegative)
    else
        #println("Calculating under λ > 0")
        f = ((τ-d.θ1)^(d.λ)) / ((d.θ2 - d.θ1)^(1+d.λ) * SpecialFunctions.beta(1+d.λ,1))
        #f = ((τ-t1)^(λ)) / ((t2 - t1)^(1+λ) * Bpositive)
    end
    return(f)
end

### define the logpdf
function Distributions.logpdf(d::BetaAdaptive, τ::Real)
    log(pdf(d, τ))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object BetaAdaptive(θ1, θ2, λ) can be passed when calling Distributions.cdf as in the standard beta: Distributions.cdf(Beta(α, β), x) 
function Distributions.cdf(d::BetaAdaptive, τ::Real)
    # standardise to use the code for the Beta distribution
    x = (τ - d.θ1) / (d.θ2 - d.θ1)
    # calculate the cdf on the standard beta with the reparametrisation α,β -> λ
    if d.λ <= 0
        cdval = Distributions.cdf(Beta(1,1-d.λ), x)
    else
        cdval = Distributions.cdf(Beta(1+d.λ,1), x)
    end
    # de-standardisation is not necessary as we are returning y, not the cdvalue
    return cdval
end

### define quantile
function Distributions.quantile(d::BetaAdaptive, q::Real)
    # standardise to use the code for the Beta distribution
    #x = (τ - θ1) / (θ2 - θ1)
    # calculate the cdf on the standard beta with the reparametrisation α,β -> λ
    if d.λ <= 0
        x = Distributions.quantile(Beta(1,1-d.λ), q)
    else
        x = Distributions.quantile(Beta(1+d.λ,1), q)
    end
    # de-standardise in order to return the quantile in the support of the BetaAdaptive
    return x * (d.θ2 - d.θ1) + d.θ1
end

# minimum of the distribution
function Distributions.minimum(d::BetaAdaptive)
    return d.θ1
end

# maximum of the distribution
function Distributions.maximum(d::BetaAdaptive)
    return d.θ2
end

# test whether a value is in the support of the distribution
function Distributions.insupport(d::BetaAdaptive, τ::Real)
    return d.θ1 <= τ <= d.θ2
end

### sampler for the BetaAdaptive
struct BetaAdaptiveSampler <: Sampleable{Univariate, Continuous}
    distribution::Beta
    θ1::Real
    θ2::Real
    λ::Real
end

# rand method for the BetaAdaptive sampler
function Base.rand(rng::AbstractRNG, d::BetaAdaptiveSampler)
    if d.λ <= 0
        α = 1
        β = 1-d.λ
    else
        α = 1+d.λ
        β = 1
    end
    sample = rand(rng, d.distribution(α,β))
    return sample * (d.θ2 - d.θ1) + d.θ1
end
