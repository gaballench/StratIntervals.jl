# using QuadGK

# code for implementing conflation of an arbitrary vector of PDFs

# but first extend pdf to take a product and a float and simply recycle it
# for passing it further to pdf
function Distributions.pdf(d::Product, x::Float64)
    xvec = [x]
    for i in 2:length(d)
        push!(xvec, x)
    end
    return pdf(d, xvec)
end

# conflation needs to be univariate wrt x. It can be vectorised by taking pdf. but this needs to be done carefully in order to check whether Turing will behave correctly in this case
function conflate(d::Product, x::Float64)
    numerat = d
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    # vector of values for product_distribution might be a metaprogramming like paste(rep('x', times =n), sep='*')
    # now integrating with the QuadGK package because Integrals simply did not work out of the box
    # integral, err = quadgk(x -> pdf(Normal(0,1),x), -Inf, Inf, rtol=1e-8)
    denominat, err = quadgk(x -> pdf(numerat, x), -Inf, Inf, rtol=1e-8)
    return pdf(numerat, x)/denominat
end

# conflation method for a matrix of floats, which are coming from a posterior predictive analysis, one column per stratigraphic interval
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
