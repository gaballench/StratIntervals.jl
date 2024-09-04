# code for implementing conflation of an arbitrary vector of PDFs

"""
    Distributions.pdf(d::Product, x::Float64)

Method for specifying the product of an arbitrary vector of pdfs to be evaluated at the same
x value.

# examples

```jldoctest
julia> using Distributions
julia> using Random
julia> Random.seed!(4);
julia> distribs = [Normal(0, 1), Normal(1, 1), Normal(2, 1)]
julia> pdf(distribs, 1.5) # evaluate the product of three normals at x=1.5
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
