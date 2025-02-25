module StratIntervals

# pick the function names in use in the code and just
# include these with using Module: fun1 fun2 etc
using Distributions
using Random
using StatsPlots
using Turing
using QuadGK
using DataFrames
using Interpolations
using KernelDensity
using SpecialFunctions

# functions to export for the users to call
export
    # distributions.jl
    ThreeParBeta,
    FourParBeta,
    RefOffExponential,
    pdf,
    logpdf,
    cdf,
    quantile,
    minimum,
    maximum,
    insupport,
    ThreeParBetaSampler,
    FourParBetaSampler,
    RefOffExponentialSampler,
    rand,
    # structs.jl
    StratInterval,
    # turingmodel.jl
    sample_stratinterval,
    # conflation.jl
    quant,
    product_distribution,
    conflate,
    pdf,
    tau_collection
    
include("distributions.jl")
include("structs.jl")
include("turingmodel.jl")
include("conflate.jl")
end
