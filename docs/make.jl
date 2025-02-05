using Documenter

# copied from the make.jl of PhyloNetworks
# and modified whenever necessary
#using Pkg
#Pkg.add(PackageSpec(name="StratIntervals", rev="main"))

using StratIntervals
using Distributions
using Turing

push!(LOAD_PATH,"../src/")
#DocMeta.setdocmeta!(StratIntervals, :DocTestSetup, :(using StratIntervals); recursive=true)
#using PhyloPlots # to trigger any precompilation warning outside jldoctests

makedocs(
    sitename = "StratIntervals.jl",
    authors = "Gustavo A. Ballen",
    modules = [StratIntervals], # to list methods from StratIntervals only, not from Base etc.
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10, size_threshold_warn = 500 * 2^10), # 600 KiB
    # exception, so warning-only for :missing_docs. List all others:
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block,
        :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes,
                                 :linkcheck, :meta_block, :parse_error, :setup_block),
    pages = [
        "Home" => "index.md",
        "Stratigraphic intervals" => "stratinterval.md",
        "Distributions" => "distributions.md",
        "Conflation of PDFs" => "conflation.md",
        "Sampling a stratigraphic interval and posterior predictive" => "turingmodel.md",
        "Public  API" => "api.md",
    ],

)

deploydocs(
    repo = "github.com/gaballench/StratIntervals.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true,
)
