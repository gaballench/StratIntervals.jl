module StratIntervals

# pick the function names in use in the code and just
# include these with using Module: fun1 fun2 etc
using Distributions
using Random
using SpecialFunctions
using StatsPlots
using Turing
using QuadGK

# functions to export for the users to call
export
    # distributions.jl
    foupar_dbeta,
    threepar_dbeta,
    BetaAdaptive,
    pdf,
    logpdf,
    cdf,
    quantile,
    minimum,
    maximum,
    insupport,
    BetaAdaptiveSampler,
    rand,
    # structs.jl
    StratInterval,
    # turingmodel.jl
    sample_stratinterval,
    # conflation.jl
    conflate,
    pdf

    
include("distributions.jl")
include("structs.jl")
include("turingmodel.jl")
include("conflate.jl")
end
