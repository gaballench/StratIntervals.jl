using Documenter

# copied from the make.jl of PhyloNetworks
# and modified whenever necessary
#using Pkg
#Pkg.add(PackageSpec(name="StratIntervals", rev="main"))

using StratIntervals
DocMeta.setdocmeta!(StratIntervals, :DocTestSetup, :(using StratIntervals); recursive=true)
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
        #"Manual" => [
        #    "Installation" => "man/installation.md",
        #    "Network manipulation" => "man/netmanipulation.md",
        #    "Input Data for SNaQ" => "man/inputdata.md",
        #    "TICR pipeline" => "man/ticr_howtogetQuartetCFs.md",
        #    "Network estimation and display" => "man/snaq_plot.md",
        #    "Network comparison and manipulation" => "man/dist_reroot.md",
        #    "Candidate Networks" => "man/fixednetworkoptim.md",
        #    "Extract Expected CFs" => "man/expectedCFs.md",
        #    "Bootstrap" => "man/bootstrap.md",
        #    "Multiple Alleles" => "man/multiplealleles.md",
        #    "Continuous Trait Evolution" => "man/trait_tree.md",
        #    "Parsimony on networks" => "man/parsimony.md",
        #    "Discrete Trait Evolution" => "man/fitDiscrete.md",
        #    "Neighbour Joining" => "man/nj.md",
        #],
        #"Library" => [
        #    "Public" => "lib/public.md",
        #    "Internals" => "lib/internals.md",
    ],

)

deploydocs(
    repo = "github.com/gaballench/StratIntervals.jl.git",
    push_preview = true,
)
