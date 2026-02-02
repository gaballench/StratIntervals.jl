
####################################
### Implementing the StratInterval distributions
####################################

################
### ThreeParBeta
################

### define a subtype
"""
    ThreeParBeta

Structure for constructing the PDF of stratigraphic interval estimation.
It is a subtype of `ContinuousUnivariateDitribution`. This structure has three fields,
corresponding to parameters: bounds on x [θ1,θ2] _in years before present_ and the preservation parameter λ. 

# examples

```@example
using Distributions
rand(ThreeParBeta(5.0, 1.0, 0.0)) # sample a random number from the ThreeParBeta with params θ1=5.0, θ2=1.0, λ=0.0
```
"""
struct ThreeParBeta <: ContinuousUnivariateDistribution
    θ1::Real
    θ2::Real
    λ::Real
end

### define the pdf
"""
    Distributions.pdf(d::ThreeParBeta, τ::Real)

This method extends the pdf function for the type ThreeParBeta so that the PDF is calculated
given the parameters θ1, θ2 and λ.
The function returns the value of the PDF at a given value of τ.

# examples

```@example
using Distributions
pdf(ThreeParBeta(5.0, 1.0, 0.0), 2.5)
```
"""
function Distributions.pdf(d::ThreeParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    if (d.λ <= 0)
        f = ((d.θ1 - τ)^(-d.λ)) / ((d.θ1 - d.θ2)^(1-d.λ) * (1/(1-d.λ)))
        #in real scale: f = ((t2 - τ)^(-λ)) / ((t2 - t1)^(1-λ) * Bnegative)
    else
        f = ((τ-d.θ2)^(d.λ)) / ((d.θ1 - d.θ2)^(1+d.λ) * (1/(1+d.λ)))
        #in real scale: f = ((τ-t1)^(λ)) / ((t2 - t1)^(1+λ) * Bpositive)
    end
    return(f)
end

### define the logpdf
"""
    Distributions.logpdf(d::ThreeParBeta, τ::Real)

This method extends the logpdf function for the type ThreeParBeta. See `Distributions.logpdf`.
The function returns the value of the logpdf at a given value of τ.

# examples

```@example
using Distributions
logpdf(ThreeParBeta(5.0, 1.0, 0.0), 2.5)
```
"""
function Distributions.logpdf(d::ThreeParBeta, τ::Real)
    log(pdf(d, τ))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object ThreeParBeta(θ1, θ2, λ) can be passed when calling Distributions.cdf as in the standard beta: Distributions.cdf(Beta(α, β), x) 
"""
    Distributions.cdf(d::ThreeParBeta, τ::Real)

This method extends the cdf function for the type ThreeParBeta. See `Distributions.cdf`.
The function returns the value of the cdf at a given value of τ.

# examples

```@example
using Distributions
cdf(ThreeParBeta(5.0, 1.0, 0.0), 2.5)
```
"""
function Distributions.cdf(d::ThreeParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # standardise to use the code for the Beta distribution
    x = (τ - d.θ2) / (d.θ1 - d.θ2)
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
"""
    Distributions.quantile(d::ThreeParBeta, q::Real)

This method extends the quantile function for the type ThreeParBeta. See `Distributions.quantile`
and `Distributions.quantile(d::UnivariateDistribution, q::Real)`.

# examples

```@example
using Distributions
quantile(ThreeParBeta(5.0, 1.0, 0.0), 0.5) # calculate the median of the distribution
```
"""
function Distributions.quantile(d::ThreeParBeta, q::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # standardise to use the code for the Beta distribution
    #x = (τ - θ1) / (θ2 - θ1)
    # calculate the cdf on the standard beta with the reparametrisation α,β -> λ
    if d.λ <= 0
        x = Distributions.quantile(Beta(1,1-d.λ), q)
    else
        x = Distributions.quantile(Beta(1+d.λ,1), q)
    end
    # de-standardise in order to return the quantile in the support of the ThreeParBeta
    return x * (d.θ1 - d.θ2) + d.θ2
end

# minimum of the distribution
"""
    Distributions.minimum(d::ThreeParBeta)

This method extends the minimum function for the type ThreeParBeta. See `Distributions.minimum`.
This function returns the _numerical_ minimum bound in the support of the ThreeParBeta distribution.
It is the parameter θ2.

# examples

```@example
using Distributions
minimum(ThreeParBeta(5.0, 1.0, 0.0))
```
"""
function Distributions.minimum(d::ThreeParBeta)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # returning the numerical minimum, which is actually the maximum in years before present
    return d.θ2
end

# maximum of the distribution
"""
    Distributions.maximum(d::ThreeParBeta)

This method extends the maximum function for the type ThreeParBeta. See `Distributions.maximum`.
This function returns the _numerical_ maximum bound in the support of the ThreeParBeta distribution.
It is the parameter θ1.

# examples

```@example
using Distributions
maximum(ThreeParBeta(5.0, 1.0, 0.0))
```
"""
function Distributions.maximum(d::ThreeParBeta)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # returning the numerical maximum, which is actually the maximum in years before present
    return d.θ1
end

# test whether a value is in the support of the distribution
"""
    Distributions.insupport(d::ThreeParBeta)

This method extends the insupport function for the type ThreeParBeta. See `Distributions.insupport`.
This function returns a Boolean being `true` if τ is in the support [θ1,θ2] of the ThreeParBeta distribution and `false` otherwise.

# examples

```@example
using Distributions
insupport(ThreeParBeta(5.0, 1.0, 0.0), 0.5) # false
insupport(ThreeParBeta(5.0, 1.0, 0.0), 3.2) # true
```
"""
function Distributions.insupport(d::ThreeParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    return d.θ2 <= τ <= d.θ1
end

### sampler for the ThreeParBeta
"""
    ThreeParBetaSampler

This structure provides the parameters and function type for the ThreeParBeta, to be used with `Distributions.rand`.
It is a subtype Sampleable of type Univariate, Continuous.
"""
struct ThreeParBetaSampler <: Sampleable{Univariate, Continuous}
    distribution::Beta
    θ1::Real
    θ2::Real
    λ::Real
end

# rand method for the ThreeParBeta sampler
"""
    Base.rand(rng::AbstractRNG, d::ThreeParBetaSampler)

This method extends the random sampler for the ThreeParBeta. It will return a single random Float64 from the ThreeParBeta
when called with only the distribution specification, or a vector of Float64 when specifying how many numbers to sample.

# examples

```@example
using Distributions
rand(ThreeParBeta(5.0, 1.0, 0.0)) # return a single random number
rand(ThreeParBeta(5.0, 1.0, 0.0), 10) # return a vector of 10 random numbers
```
"""
function Base.rand(rng::AbstractRNG, d::ThreeParBetaSampler)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    if d.λ <= 0
        α = 1
        β = 1-d.λ
    else
        α = 1+d.λ
        β = 1
    end
    sample = rand(rng, d.distribution(α,β))
    return sample * (d.θ1 - d.θ2) + d.θ2
end

################
### FourParBeta
################

### define a subtype
"""
    FourParBeta

Structure for constructing the PDF of stratigraphic interval estimation.
It is a subtype of `ContinuousUnivariateDitribution`. This structure has four fields,
corresponding to parameters: bounds on x [θ1,θ2] _in years before present_ and the preservation parameters α, and β. 

# examples

```@example
using Distributions
rand(FourParBeta(5.0, 1.0, 1.0, 1.0)) # sample a random number from the FourParBeta with params θ1=5.0, θ2=1.0, α=1.0, β=1.0
```
"""
struct FourParBeta <: ContinuousUnivariateDistribution
    θ1::Real
    θ2::Real
    α::Real
    β::Real
end

### define the pdf
"""
    Distributions.pdf(d::FouParBeta, τ::Real)

This method extends the pdf function for the type FourParBeta so that the PDF is calculated
given the parameters θ1, θ2, α, and β.
The function returns the value of the PDF at a given value of τ.

# examples

```@example
using Distributions
pdf(FourParBeta(5.0, 1.0, 1.0, 1.0), 2.5)
```
"""
function Distributions.pdf(d::FourParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    if (d.α <= 0.0 | d.β <= 0.0)
        error("Both α and β need to be larger than 0")
    end
    f = ((τ-d.θ2)^(d.α - 1) * (d.θ1-τ)^(d.β - 1)) / ((d.θ1-d.θ2)^(d.α+d.β-1) * SpecialFunctions.beta(d.α,d.β))
    return(f)
end

### define the logpdf
"""
    Distributions.logpdf(d::FourParBeta, τ::Real)

This method extends the logpdf function for the type FourParBeta. See `Distributions.logpdf`.
The function returns the value of the logpdf at a given value of τ.

# examples

```@example
using Distributions
logpdf(FourParBeta(5.0, 1.0, 1.0, 1.0), 2.5)
```
"""
function Distributions.logpdf(d::FourParBeta, τ::Real)
    log(pdf(d, τ))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object FourParBeta(θ1, θ2, α, β) can be passed when calling Distributions.cdf as in the standard beta: Distributions.cdf(Beta(α, β), x) 
"""
    Distributions.cdf(d::FourParBeta, τ::Real)

This method extends the cdf function for the type FourParBeta. See `Distributions.cdf`.
The function returns the value of the cdf at a given value of τ.

# examples

```@example
using Distributions
cdf(FourParBeta(5.0, 1.0, 1.0, 1.0), 2.5)
```
"""
function Distributions.cdf(d::FourParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # standardise to use the code for the Beta distribution
    x = (τ - d.θ2) / (d.θ1 - d.θ2)
    # calculate the cdf on the standard beta
    cdval = Distributions.cdf(Beta(d.α,d.β), x)
    # de-standardisation is not necessary as we are returning y, not the cdvalue
    return cdval
end

### define quantile
"""
    Distributions.quantile(d::FourParBeta, q::Real)

This method extends the quantile function for the type FourParBeta. See `Distributions.quantile`
and `Distributions.quantile(d::UnivariateDistribution, q::Real)`.

# examples

```@example
using Distributions
quantile(FourParBeta(5.0, 1.0, 1.0, 1.0), 0.5) # calculate the median of the distribution
```
"""
function Distributions.quantile(d::FourParBeta, q::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # standardise to use the code for the Beta distribution
    #x = (τ - θ1) / (θ2 - θ1)
    # calculate the cdf on the standard beta
    x = Distributions.quantile(Beta(d.α,d.β), q)
    # de-standardise in order to return the quantile in the support of the FourParBeta
    return x * (d.θ1 - d.θ2) + d.θ2
end

# minimum of the distribution
"""
    Distributions.minimum(d::FourParBeta)

This method extends the minimum function for the type FourParBeta. See `Distributions.minimum`.
This function returns the _numerical_ minimum bound in the support of the FourParBeta distribution.
It is the parameter θ2.

# examples

```@example
using Distributions
minimum(FourParBeta(5.0, 1.0, 1.0, 1.0))
```
"""
function Distributions.minimum(d::FourParBeta)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # returning the numerical minimum, which is actually the maximum in years before present
    return d.θ2
end

# maximum of the distribution
"""
    Distributions.maximum(d::FourParBeta)

This method extends the maximum function for the type FourParBeta. See `Distributions.maximum`.
This function returns the _numerical_ maximum bound in the support of the FourParBeta distribution.
It is the parameter θ1.

# examples

```@example
using Distributions
maximum(FourParBeta(5.0, 1.0, 1.0, 1.0))
```
"""
function Distributions.maximum(d::FourParBeta)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    # returning the numerical maximum, which is actually the maximum in years before present
    return d.θ1
end

# test whether a value is in the support of the distribution
"""
    Distributions.insupport(d::FourParBeta)

This method extends the insupport function for the type FourParBeta. See `Distributions.insupport`.
This function returns a Boolean being `true` if τ is in the support [θ1,θ2] of the FourParBeta distribution and `false` otherwise.

# examples

```@example
using Distributions
insupport(FourParBeta(5.0, 1.0, 1.0, 1.0), 0.5) # false
insupport(FourParBeta(5.0, 1.0, 1.0, 1.0), 3.2) # true
```
"""
function Distributions.insupport(d::FourParBeta, τ::Real)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    return d.θ2 <= τ <= d.θ1
end

### sampler for the FourParBeta
"""
    FourParBetaSampler

This structure provides the parameters and function type for the FourParBeta, to be used with `Distributions.rand`.
It is a subtype Sampleable of type Univariate, Continuous.
"""
struct FourParBetaSampler <: Sampleable{Univariate, Continuous}
    distribution::Beta
    θ1::Real
    θ2::Real
    α::Real
    β::Real
end

# rand method for the FourParBeta sampler
"""
    Base.rand(rng::AbstractRNG, d::FourParBetaSampler)

This method extends the random sampler for the FourParBeta. It will return a single random Float64 from the FourParBeta
when called with only the distribution specification, or a vector of Float64 when specifying how many numbers to sample.

# examples

```@example
using Distributions
rand(FourParBeta(5.0, 1.0, 1.0, 1.0)) # return a single random number
rand(FourParBeta(5.0, 1.0, 1.0, 1.0), 10) # return a vector of 10 random numbers
```
"""
function Base.rand(rng::AbstractRNG, d::FourParBetaSampler)
    if (d.θ2 > d.θ1)
        error("θ2 is greater than θ1, time needs to be in years before present, the other way around")
    end
    sample = rand(rng, d.distribution(d.α,d.β))
    return sample * (d.θ1 - d.θ2) + d.θ2
end

##############################
### RefOff distributions
### a.k.a. reflected-offset
##############################

################
### RefOffExponential
################

### define a subtype
"""
    RefOffExponential

Structure for constructing the PDF of a reflected and offset exponential distribution.
It is a subtype of `ContinuousUnivariateDitribution`. This structure has three fields,
corresponding to parameters: theta, offset, and reflect. 

# examples

```@example
using Distributions
rand(RefOffExponential(1, 10.0, 1)) # sample a random number from the RefOffExponential with params θ=1.0, o=10.0, ρ=1
```
"""
struct RefOffExponential <: ContinuousUnivariateDistribution
    θ::Real
    o::Real
    ρ::Real
end

### define the pdf
"""
    Distributions.pdf(d::RefOffExponential, x::Real)

This method extends the pdf for the type RefOffExponential so that the PDF is calculated
given the parameters θ, o, and ρ.
The function returns the value of the PDF at a given value of x.

# examples

```@example
using Distributions
pdf(RefOffExponential(1, 10.0, 1), 2.5)
```
"""
function Distributions.pdf(d::RefOffExponential, x::Real)
    f = pdf(Exponential(d.θ), d.ρ*(x - d.o))
    return(f)
end

### define the logpdf
"""
    Distributions.logpdf(d::RefOffExponential, x::Real)

This method extends the logpdf function for the type RefOffExponential. See `Distributions.logpdf`.
The function returns the value of the logpdf at a given value of x.

# examples

```@example
using Distributions
logpdf(RefOffExponential(1, 10.0, 1), 2.5)
```
"""
function Distributions.logpdf(d::RefOffExponential, x::Real)
    log(pdf(d, x))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object RefOffExponential(θ, o, ρ) can be passed when calling Distributions.cdf as in the standard exponential: Distributions.cdf(Exponential(θ), x) 
"""
    Distributions.cdf(d::RefOffExponential, x::Real)

This method extends the cdf function for the type RefOffExponential. See `Distributions.cdf`.
The function returns the value of the cdf at a given value of x.

# examples

```@example
using Distributions
cdf(RefOffExponential(1, 10.0, 1), 2.5)
```
"""
function Distributions.cdf(d::RefOffExponential, x::Real)
    # standardise to use the code for the Exponential distribution
    # first undo the offset
    x = x + d.o
    # calculate the cdf on the standard beta
    # if reflected, this is actually de-reflects and does 1-cdf
    if (d.ρ == 1)
        cdval = Distributions.cdf(Exponential(d.θ), x)
    else
        x = d.ρ*x
        cdval = 1 - Distributions.cdf(Exponential(d.θ), x)
    end
    # de-standardisation is not necessary as we are returning y, not the cdvalue
    return cdval
end

### define quantile
"""
    Distributions.quantile(d::RefOffExponential, q::Real)

This method extends the quantile function for the type RefOffExponential. See `Distributions.quantile`
and `Distributions.quantile(d::UnivariateDistribution, q::Real)`.

# examples

```@example
using Distributions
quantile(RefOffExponential(1, 10.0, 1), 0.5) # calculate the median of the distribution
```
"""
function Distributions.quantile(d::RefOffExponential, q::Real)
    # calculate the quantile on the standard beta, if reflected, calculate the complement on the reflected standard
    if (d.ρ == 1)
        x = d.o + Distributions.quantile(Exponential(d.θ), q)
    else
        x = d.o + Distributions.quantile(Exponential(d.θ), 1-q)
    end
    # re-reflect and re-offset if necessary
    
    return d.ρ*(x - d.o)
end

# minimum of the distribution
"""
    Distributions.minimum(d::RefOffExponential)

This method extends the minimum function for the type RefOffExponential. See `Distributions.minimum`.
This function returns the _numerical_ minimum bound in the support of the RefOffExponential distribution.
It depends on whether it is reflected and offset or not.

# examples

```@example
using Distributions
minimum(RefOffExponential(1, 10.0, 1))
```
"""
function Distributions.minimum(d::RefOffExponential)
    # if reflected, -Inf, else, offset
    if (d.ρ == -1)
        x = -Inf
    else
        x = d.o
    end    
    return x 
end

# maximum of the distribution
"""
    Distributions.maximum(d::RefOffExponential)

This method extends the maximum function for the type RefOffExponential. See `Distributions.maximum`.
This function returns the _numerical_ maximum bound in the support of the RefOffExponential distribution.
It depends on whether it is reflected and offset or not.

# examples

```@example
using Distributions
minimum(RefOffExponential(1, 10.0, 1))
```
"""
function Distributions.maximum(d::RefOffExponential)
    # standard, Inf, else, offset
    if (d.ρ == 1)
        x = Inf
    else
        x = d.o
    end    
    return x 
end

# test whether a value is in the support of the distribution
"""
    Distributions.insupport(d::RefOffExponential, x::Real)

This method extends the insupport function for the type RefOffExponential. See `Distributions.insupport`.
This function returns a Boolean being `true` if x is in the support (-Inf,offset] when ρ = -1, or [offset,Inf) when
standard, possible offset, and `false` otherwise.

# examples

```@example
using Distributions
insupport(RefOffExponential(1, 10.0, 1), 0.5) # false
insupport(RefOffExponential(1, 10.0, 1), 12.0) # true
```
"""
function Distributions.insupport(d::RefOffExponential, x::Real)
    if (d.ρ == 1)
        return d.o <= x < Inf
    end
    if (d.ρ == -1)
        return -Inf < x <= d.o
    end
end

### sampler for the RefOffExponential
"""
    RefOffExponentialSampler

This structure provides the parameters and function type for the RefOffExponential, to be used with `Distributions.rand`.
It is a subtype Sampleable of type Univariate, Continuous.
"""
struct RefOffExponentialSampler <: Sampleable{Univariate, Continuous}
    distribution::Exponential
    θ::Real
    o::Real
    ρ::Real
end

# rand method for the RefOffExponential sampler
"""
    Base.rand(rng::AbstractRNG, d::RefOffExponentialSampler)

This method extends the random sampler for the RefOffExponential. It will return a single random Float64 from the RefOffExponential
when called with only the distribution specification, or a vector of Float64 when specifying how many numbers to sample.

# examples

```@example
using Distributions
rand(RefOffExponential(1, 10.0, 1)) # return a single random number
rand(RefOffExponential(1, 10.0, 1), 10) # return a vector of 10 random numbers
```
"""
function Base.rand(rng::AbstractRNG, d::RefOffExponentialSampler)
    #sample = rand(rng, d.distribution(d.θ))
    #return d.ρ * (sample - d.o) # original
    return -1
end
