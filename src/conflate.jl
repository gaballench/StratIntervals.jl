# code for implementing conflation of an arbitrary vector of PDFs

# calculate quantiles for conflated distributions

"""
    quant(taus::Array{Float64, 2}, p::Float64, lower_tail::Bool, min::Float64, max::Float64)

A function for calculating the quantile for a given probability value on a conflated distribution.
Whether to calculate lower or upper tail is available.
The method uses constrained numerical optimisation, and therefore a minimum and maximum need to be specified.

taus a matrix of stacked empirical distributions, one distribution per column
p the desired probability
lower_tail a boolean defining whether we want to calculate the probability from -Inf to quant, if false it calculates the prob from quant to Inf
min minimum value for the box optimisation, normally 0.0
max maximum value for the box optimisation, this should be something sensible outside of the bulk of the conflated PDF

# examples

```@example
using Random
using Distributions
using StatsPlots
using StratIntervals
Random.seed!(1507)

#empirical1
empirical1_vec = rand(Normal(0, 1), 20)
#empirical2
empirical2_vec = rand(Uniform(1, 2), 20)
# stack them together in a matrix, column-wise
empirical_matrix = [empirical1_vec empirical2_vec]

xx = -0.0:0.01:800.0
yy = map(x -> conflate(empirical_matrix, x), xx);

plot(xx, yy, label="Conflation")
density!(empirical1_vec, label="Empirical 1")
density!(empirical2_vec, label="Empirical 2")

quant(empirical_matrix, 0.025, true, 0.0, 800.0)
quant(empirical_matrix, 0.975, true, 0.0, 800.0)
quant(empirical_matrix, 0.025, false, 0.0, 800.0)
quant(empirical_matrix, 0.5, false, 0.0, 800.0)
```
"""
function quant(taus::Array{Float64, 2}, p::Float64, lower_tail::Bool, min::Float64, max::Float64)
    # from conflate(Array{Float64,2})
    interpolations = map(i -> InterpKDE(kde(taus[:, i])), 1:size(taus,2))
    # numerat not needed beforehand
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    denominat, err = quadgk(x -> prod(map(dd -> pdf(dd, x), interpolations)), -Inf, Inf, rtol=1e-8)
    # now we need to wrap the product pdf calculation
    function wrapper_conflate(x)
        prod(map(dd -> pdf(dd, x), interpolations))/denominat
    end
    # then depending on lower or upper tail, choose the proper minimiser function
    f(y) = quadgk(x -> wrapper_conflate(x), -Inf, y)[1]
    g(y) = quadgk(x -> wrapper_conflate(x), y, Inf)[1]
    if lower_tail
        res = optimize(x -> (f(x)-p)^2, min, max)
    end
    if !lower_tail
        res = optimize(x -> (g(x)-p)^2, min, max)
    end
    return [Optim.minimizer(res), res]
end

"""
    quant(taus::Vector{StratInterval}, p::Float64, lower_tail::Bool, min::Float64, max::Float64)

A function for calculating the quantile for a given probability value on a conflated distribution.
Whether to calculate lower or upper tail is available.
The method uses constrained numerical optimisation, and therefore a minimum and maximum need to be specified.

taus a vector which comes from calling sample_stratinterval(vector{StratInterval})
p the desired probability
lower_tail a boolean defining whether we want to calculate the probability from -Inf to quant, if false it calculates the prob from quant to Inf
min minimum value for the box optimisation, normally 0.0
max maximum value for the box optimisation, this should be something sensible outside of the bulk of the conflated PDF

# examples

```@example
using Random
using Distributions
using StatsPlots
using StratIntervals
Random.seed!(1507)

setprogress!(false)

ndata = 100;
iters = 10000;

# first interval
true_lambda_1 = 0;
true_theta2_1 = 10;
true_theta1_1 = 22;
data_1 = rand(ThreeParBeta(true_theta1_1, true_theta2_1, true_lambda_1), ndata);

# second interval
true_lambda_2 = 1;
true_theta2_2 = 15;
true_theta1_2 = 30;
data_2 = rand(ThreeParBeta(true_theta1_2, true_theta2_2, true_lambda_2), ndata);

# specifying each StratInterval object
interval_1 = StratInterval(data_1, Normal(22, 3), Exponential(10), Normal(0,3));
interval_2 = StratInterval(data_2, Normal(30, 3), Exponential(15), Normal(1,3));

# construct the vector of StratIntervals
vecinterval = [interval_1, interval_2];

# MCMC sampling and posterior predictive
mystratint_postpredict_vec = sample_stratinterval(vecinterval, iters, NUTS(), false, true);

# now calculate the conflation of posterior predictives for the collection interval_1, interval_2

xx = 0:0.01:40;
yy = map(x -> tau_collection(mystratint_postpredict_vec, x), xx);

plot(xx, yy, label="Conflation")
density!(mystratint_postpredict_vec[1][2], label="Postpredictive 1")
density!(mystratint_postpredict_vec[2][2], label="Postpredictive 2")

quant(mystratint_postpredict_vec, 0.025, true, 0.0, 800.0)
quant(mystratint_postpredict_vec, 0.975, true, 0.0, 800.0)
quant(mystratint_postpredict_vec, 0.025, false, 0.0, 800.0)
quant(mystratint_postpredict_vec, 0.5, false, 0.0, 800.0)
```
"""
function quant(taus::Vector{StratInterval}, p::Float64, lower_tail::Bool, min::Float64, max::Float64)
    # part of this comes from tau_collection, and part from conflate
    # from tau_collection
    postpredict_matrix = reduce(hcat, getindex.(taus, 2))
    # from conflate(Array{Float64,2})
    interpolations = map(i -> InterpKDE(kde(postpredict_matrix[:, i])), 1:size(postpredict_matrix,2))
    # numerat not needed beforehand
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    denominat, err = quadgk(x -> prod(map(dd -> pdf(dd, x), interpolations)), -Inf, Inf, rtol=1e-8)
    # now we need to wrap the product pdf calculation
    function wrapper_conflate(x)
        prod(map(dd -> pdf(dd, x), interpolations))/denominat
    end
    # then depending on lower or upper tail, choose the proper minimiser function
    f(y) = quadgk(x -> wrapper_conflate(x), -Inf, y)[1]
    g(y) = quadgk(x -> wrapper_conflate(x), y, Inf)[1]
    if lower_tail
        res = optimize(x -> (f(x)-p)^2, min, max)
    end
    if !lower_tail
        res = optimize(x -> (g(x)-p)^2, min, max)
    end
    return [Optim.minimizer(res), res]
end

"""
    Distributions.product_distribution(dists::Vector{Normal{Float64}})

Method for making product of Normals to behave as the generic one. The specialised one returning a MvNormal generates issues during conflation when all distribs are normal.

The methods associated to Product are going to be deprecated from Distributions. Prepare to modify this method accordingly.

"""
function Distributions.product_distribution(dists::Vector{Normal{Float64}})
    return Product(dists)
end


"""
    Distributions.pdf(d::Product, x::Float64)

Method for specifying the product of an arbitrary vector of pdfs to be evaluated at the same
x value.

Note: This method may need to be changed when Distributions deprecate Product.
See product_distribution(dists::Vector{Normal{Float64}}).

# examples

```@example
using Distributions
using Random
Random.seed!(4);
distribs = [Normal(0, 1), Normal(1, 1), Normal(2, 1)]
pdf(product_distribution(distribs), 1.5) # evaluate the product of three normals at x=1.5
```
"""
function Distributions.pdf(d::Product, x::Float64)
    xvec = [x]
    for i in 2:length(d)
        push!(xvec, x)
    end
    return pdf(d, xvec)
end

# conflation needs to be univariate wrt x. It can be vectorised by taking pdf. 

"""
    conflate(d::Vector{ContinuousUnivariateDitribution}, x::Float64)

Conflation method for a vector of distributions.

It takes a vector of pdfs and conflate then evaluating the resulting pdf at x.

This is a convenience method for allowing either input as vector or product_distribution
"""
function conflate(d::Vector{T} where T <: ContinuousUnivariateDistribution, x::Float64)
    return conflate(product_distribution(d), x)
end

"""
    conflate(d::Product, x::Float64)

Conflation method for a product distribution, which are coming from calling pdf on a vector of pdfs.

It takes a vector of pdfs and conflate then evaluating the resulting pdf at x.

For conflation of empirical distributions see `conflate(d::Array{Float64, 2}, x::Float64)`@ref
"""
function conflate(d::Product, x::Float64)
    numerat = d
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    # vector of values for product_distribution might be a metaprogramming like paste(rep('x', times =n), sep='*')
    # now integrating with the QuadGK package because Integrals simply did not work out of the box
    # integral, err = quadgk(x -> pdf(Normal(0,1),x), -Inf, Inf, rtol=1e-8)
    denominat, err = quadgk(x -> pdf(numerat, x), -Inf, Inf, rtol=1e-8)
    return pdf(numerat, x)/denominat
end

"""
    conflate(d::Array{Float64, 2}, x::Float64)

Conflation method for a matrix of floats, which are coming from a posterior predictive analysis,
one column per stratigraphic interval.

It takes a matrix of empirical values to smooth, column-wise, and then apply the KDE and interpolation on each column.

Then it applies the conflation to these empirical pdfs and evaluate it at x.
"""
function conflate(d::Array{Float64, 2}, x::Float64)
    # we need to wrap everything related to InterpKDE inside the method Array{Flot64, 2} as the type hierarchy does not allow to easily create a method for it
    # take a matrix of empirical values to smooth, column-wise, and then apply the KDE and interpolation on each column
    # then we will have a vector of elements of type InterpKDE which can be called using the pdf method.
    interpolations = map(i -> InterpKDE(kde(d[:, i])), 1:size(d,2))
    # numerat not needed beforehand
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    denominat, err = quadgk(x -> prod(map(dd -> pdf(dd, x), interpolations)), -Inf, Inf, rtol=1e-8)
    return prod(map(dd -> pdf(dd, x), interpolations))/denominat
end

"""
    tau_collection(taus, x::Float64)

Calculation of the pdf for tau as a conflation of StratIntervals.

This uses the output of `sample_interval` with a `Vector{StratInterval}` and where `postpredict` is `true`.

The value of x is the evaluated using the conflation of the collection of stratigraphic intervals, that is,
the tau for the co-occurrence of `StratInterval`s.
"""
function tau_collection(taus, x::Float64)
    # is taus a vector?
    if !(taus isa Vector{Any})
        error("taus must be of type Vector{Any}, the result of calling `sample_interval` on a Vector{StratInterval} where argument `postpredict` is true")
    end
    if length(taus) < 2
        error("taus must be at least of legnth = 2")
    end
    if !(taus[1][1] isa Chains) | !(taus[1][2] isa Vector{Float64})
        error("taus must have tuple elements of type Chains and Vector{Float64}, the result of calling `sample_interval` on a Vector{StratInterval} where argument `postpredict` is true")
    end
    # take all the postpredict vectors and build a matrix
    postpredict_matrix = reduce(hcat, getindex.(taus, 2))
    # return the call to conflate evaluated at the value of x
    return conflate(postpredict_matrix, x)
end
