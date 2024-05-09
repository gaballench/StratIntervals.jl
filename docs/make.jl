using StratIntervals
using Documenter

DocMeta.setdocmeta!(StratIntervals, :DocTestSetup, :(using StratIntervals); recursive=true)

makedocs(;
    modules=[StratIntervals],
    authors="Gustavo A. Ballen",
    sitename="StratIntervals.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
