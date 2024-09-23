#using Distributions
#using SpecialFunctions # for the beta function gamma(a)gamma(b)/gamma(a+b) 
#using Base
#using Random

"""
    fourpar_dbeta(τ,θ1,θ2,α,β)

The four-parameter Beta PDF. This is not meant to be used directly but just for testing.
This is a Beta where the support is not bound to the interval [0,1] but rather two arbitrary
parameters [θ1,θ2]. Shape parameters α,β are the usual ones in the Beta PDF.
The function returns the value of the PDF at a given value of τ.

# examples

```jldoctest
julia> fourpar_dbeta(2.56, 1.0, 5.0, 1, 1) alpha and beta = 1
```
"""
function fourpar_dbeta(τ,θ1,θ2,α,β) 
    f = ((θ-θ1)^(α - 1) * (θ2-θ)^(β - 1)) / ((θ2-θ1)^(α+β-1) * SpecialFunctions.beta(α,β))
    return(f)
end

### functions with reparam
"""
    threepar_dbeta(θ,θ1,θ2,λ)

The three-parameter Beta PDF. This is not meant to be used directly but just for testing.
This is a Beta where the support is not bound to the interval [0,1] but rather two arbitrary
parameters [θ1,θ2]. Shape parameters are reparametrised to λ instead of the usual shape parameters
α,β in the Beta PDF.
The function returns the value of the PDF at a given value of τ.

# examples

```jldoctest
julia> threepar_dbeta(2.56, 1.0, 5.0, 2.0) # lambda = 2
```
"""
function threepar_dbeta(τ,θ1,θ2,λ)
    if (λ <= 0)
        #println("Calculating under λ <= 0")
        f = ((θ2 - θ)^(-λ)) / ((θ2 - θ1)^(1-λ) * SpecialFunctions.beta(1,1-λ))
        #f = ((θ2 - θ)^(-λ)) / ((θ2 - θ1)^(1-λ) * Bnegative)
    else
        #println("Calculating under λ > 0")
        f = ((θ-θ1)^(λ)) / ((θ2 - θ1)^(1+λ) * SpecialFunctions.beta(1+λ,1))
        #f = ((θ-θ1)^(λ)) / ((θ2 - θ1)^(1+λ) * Bpositive)
    end
    return(f)    
end

####################################
### Implementing the distribution in Julia
####################################

### define a subtype
# this probably needs to be defined as imutable with the parameters θ. and λ
"""
    ThreeParBeta

Structure for constructing the PDF of stratigraphic interval estimation.
It is a subtype of `ContinuousUnivariateDitribution`. This structure has three fields,
corresponding to parameters: bounds on x [θ1,θ2] and the preservation parameter λ. 

# examples

```jldoctest
julia> using Distributions
julia> rand(ThreeParBeta(1.0, 5.0, 0.0)) # sample a random number from the ThreeParBeta with params θ1=1.0, θ2=5.0, λ=0.0
```
"""
struct ThreeParBeta <: ContinuousUnivariateDistribution
    θ1::Real
    θ2::Real
    λ::Real
end

### declaring methods in a quick way, check whether this is correct
### define the pdf
"""
    Distributions.pdf(d::ThreeParBeta, τ::Real)

This method extends the pdf function for the type ThreeParBeta so that the PDF is calculated
given the parameters θ1, θ2 and λ.
The function returns the value of the PDF at a given value of τ.

# examples

```jldoctest
julia> using Distributions
julia> pdf(ThreeParBeta(1.0, 5.0, 0.0), 2.5)
```
"""
function Distributions.pdf(d::ThreeParBeta, τ::Real)
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
"""
    Distributions.logpdf(d::ThreeParBeta, τ::Real)

This method extends the logpdf function for the type ThreeParBeta. See `Distributions.pdf`.
The function returns the value of the logpdf at a given value of τ.

# examples

```jldoctest
julia> using Distributions
julia> logpdf(ThreeParBeta(1.0, 5.0, 0.0), 2.5)
```
"""
function Distributions.logpdf(d::ThreeParBeta, τ::Real)
    log(pdf(d, τ))
end

### define the CDF
# this needs to follow the usage of _delegate_statsfuns so that the object ThreeParBeta(θ1, θ2, λ) can be passed when calling Distributions.cdf as in the standard beta: Distributions.cdf(Beta(α, β), x) 
"""
    Distributions.cdf(d::ThreeParBeta, τ::Real)

This method extends the cdf function for the type ThreeParBeta. See `Distributions.pdf`.
The function returns the value of the cdf at a given value of τ.

# examples

```jldoctest
julia> using Distributions
julia> cdf(ThreeParBeta(1.0, 5.0, 0.0), 2.5)
```
"""
function Distributions.cdf(d::ThreeParBeta, τ::Real)
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
"""
    Distributions.quantile(d::ThreeParBeta, q::Real)

This method extends the quantile function for the type ThreeParBeta. See `Distributions.pdf`
and `Distributions.quantile(d::UnivariateDistribution, q::Real)`.

# examples

```jldoctest
julia> using Distributions
julia> quantile(ThreeParBeta(1.0, 5.0, 0.0), 0.5) # calculate the median of the distribution
```
"""
function Distributions.quantile(d::ThreeParBeta, q::Real)
    # standardise to use the code for the Beta distribution
    #x = (τ - θ1) / (θ2 - θ1)
    # calculate the cdf on the standard beta with the reparametrisation α,β -> λ
    if d.λ <= 0
        x = Distributions.quantile(Beta(1,1-d.λ), q)
    else
        x = Distributions.quantile(Beta(1+d.λ,1), q)
    end
    # de-standardise in order to return the quantile in the support of the ThreeParBeta
    return x * (d.θ2 - d.θ1) + d.θ1
end

# minimum of the distribution
"""
    Distributions.minimum(d::ThreeParBeta)

This method extends the minimum function for the type ThreeParBeta. See `Distributions.pdf`.
This function returns the minimum bound in the support of the ThreeParBeta distribution.
It is the parameter θ1.

# examples

```jldoctest
julia> using Distributions
julia> minimum(ThreeParBeta(1.0, 5.0, 0.0))
```
"""
function Distributions.minimum(d::ThreeParBeta)
    return d.θ1
end

# maximum of the distribution
"""
    Distributions.maximum(d::ThreeParBeta)

This method extends the maximum function for the type ThreeParBeta. See `Distributions.pdf`.
This function returns the maximum bound in the support of the ThreeParBeta distribution.
It is the parameter θ2.

# examples

```jldoctest
julia> using Distributions
julia> maximum(ThreeParBeta(1.0, 5.0, 0.0))
```
"""
function Distributions.maximum(d::ThreeParBeta)
    return d.θ2
end

# test whether a value is in the support of the distribution
"""
    Distributions.insupport(d::ThreeParBeta)

This method extends the insupport function for the type ThreeParBeta. See `Distributions.pdf`.
This function returns a Boolean being `true` if τ is in the support [θ1,θ2] of the ThreeParBeta distribution and `false` otherwise.

# examples

```jldoctest
julia> using Distributions
julia> insupport(ThreeParBeta(1.0, 5.0, 0.0), 0.5) # false
julia> insupport(ThreeParBeta(1.0, 5.0, 0.0), 3.2) # true
```
"""
function Distributions.insupport(d::ThreeParBeta, τ::Real)
    return d.θ1 <= τ <= d.θ2
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

```jldoctest
julia> using Distributions
julia> rand(ThreeParBeta(1.0, 5.0, 0.0)) # return a single random number
julia> rand(ThreeParBeta(1.0, 5.0, 0.0), 10) # return a vector of 10 random numbers
"""
function Base.rand(rng::AbstractRNG, d::ThreeParBetaSampler)
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
