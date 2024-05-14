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
    numerator = d
    #evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
    # vector of values for product_distribution might be a metaprogramming like paste(rep('x', times =n), sep='*')
    # now integrating with the QuadGK package because Integrals simply did not work out of the box
    # integral, err = quadgk(x -> pdf(Normal(0,1),x), -Inf, Inf, rtol=1e-8)
    denominator, err = quadgk(x -> pdf(numerator, x), -Inf, Inf, rtol=1e-8)
    return pdf(numerator, x)/denominator
end
