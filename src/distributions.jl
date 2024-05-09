#using Distributions
#using SpecialFunctions # for the beta function gamma(a)gamma(b)/gamma(a+b) 
#using Base
#using Random

function fourpar_dbeta(tau,a,b,t1,t2) 
    f = ((tau-t1)^(a - 1) * (t2-tau)^(b - 1)) / ((t2-t1)^(a+b-1) * SpecialFunctions.beta(a,b))
    return(f)
end

### functions with reparam
function threepar_dbeta(tau,t1,t2,lambda)
    if (lambda <= 0)
        #println("Calculating under lambda <= 0")
        f = ((t2 - tau)^(-lambda)) / ((t2 - t1)^(1-lambda) * SpecialFunctions.beta(1,1-lambda))
        #f = ((t2 - tau)^(-lambda)) / ((t2 - t1)^(1-lambda) * Bnegative)
    else
        #println("Calculating under lambda > 0")
        f = ((tau-t1)^(lambda)) / ((t2 - t1)^(1+lambda) * SpecialFunctions.beta(1+lambda,1))
        #f = ((tau-t1)^(lambda)) / ((t2 - t1)^(1+lambda) * Bpositive)
    end
    return(f)    
end

####################################
### Implementing the distribution in Julia
####################################

### define a subtype
# this probably needs to be defined as imutable with the parameters theta. and lambda
struct BetaAdaptive <: ContinuousUnivariateDistribution
    theta1::Real
    theta2::Real
    lambda::Real
end

### declaring methods in a quick way, check whether this is correct
### define the pdf
function Distributions.pdf(d::BetaAdaptive, tau::Real)
    #threepar_dbeta(d.tau, d.theta1, d.theta2, d.lambda)
    if (d.lambda <= 0)
        #println("Calculating under lambda <= 0")
        f = ((d.theta2 - tau)^(-d.lambda)) / ((d.theta2 - d.theta1)^(1-d.lambda) * SpecialFunctions.beta(1,1-d.lambda))
        #f = ((t2 - tau)^(-lambda)) / ((t2 - t1)^(1-lambda) * Bnegative)
    else
        #println("Calculating under lambda > 0")
        f = ((tau-d.theta1)^(d.lambda)) / ((d.theta2 - d.theta1)^(1+d.lambda) * SpecialFunctions.beta(1+d.lambda,1))
        #f = ((tau-t1)^(lambda)) / ((t2 - t1)^(1+lambda) * Bpositive)
    end
    return(f)
end

### define the logpdf
function Distributions.logpdf(d::BetaAdaptive, tau::Real)
    log(pdf(d, tau))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object BetaAdaptive(theta1, theta2, lambda) can be passed when calling Distributions.cdf as in the standard beta: Distributions.cdf(Beta(alpha, beta), x) 
function Distributions.cdf(d::BetaAdaptive, tau::Real)
    # standardise to use the code for the Beta distribution
    x = (tau - d.theta1) / (d.theta2 - d.theta1)
    # calculate the cdf on the standard beta with the reparametrisation alpha,beta -> lambda
    if d.lambda <= 0
        cdval = Distributions.cdf(Beta(1,1-d.lambda), x)
    else
        cdval = Distributions.cdf(Beta(1+d.lambda,1), x)
    end
    # de-standardisation is not necessary as we are returning y, not the cdvalue
    return cdval
end

### define quantile
function Distributions.quantile(d::BetaAdaptive, q::Real)
    # standardise to use the code for the Beta distribution
    #x = (tau - theta1) / (theta2 - theta1)
    # calculate the cdf on the standard beta with the reparametrisation alpha,beta -> lambda
    if d.lambda <= 0
        x = Distributions.quantile(Beta(1,1-d.lambda), q)
    else
        x = Distributions.quantile(Beta(1+d.lambda,1), q)
    end
    # de-standardise in order to return the quantile in the support of the BetaAdaptive
    return x * (d.theta2 - d.theta1) + d.theta1
end

# minimum of the distribution
function Distributions.minimum(d::BetaAdaptive)
    return d.theta1
end

# maximum of the distribution
function Distributions.maximum(d::BetaAdaptive)
    return d.theta2
end

# test whether a value is in the support of the distribution
function Distributions.insupport(d::BetaAdaptive, tau::Real)
    return d.theta1 <= tau <= d.theta2
end

### sampler for the BetaAdaptive
struct BetaAdaptiveSampler <: Sampleable{Univariate, Continuous}
    distribution::Beta
    theta1::Real
    theta2::Real
    lambda::Real
end

# rand method for the BetaAdaptive sampler
function Base.rand(rng::AbstractRNG, d::BetaAdaptiveSampler)
    if d.lambda <= 0
        alpha = 1
        beta = 1-d.lambda
    else
        alpha = 1+d.lambda
        beta = 1
    end
    sample = rand(rng, d.distribution(alpha,beta))
    return sample * (d.theta2 - d.theta1) + d.theta1
end
