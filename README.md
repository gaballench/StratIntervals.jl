# StratIntervals <img src="docs/src/logo.png" height="140" align="right"></img>

[![Build Status](https://github.com/gaballench/StratIntervals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gaballench/StratIntervals.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![version](https://juliahub.com/docs/General/StratIntervals/stable/version.svg)](https://juliahub.com/ui/Packages/General/StratIntervals)

## Description

A Julia package for Bayesian inference of stratigraphic intervals using time occurrence data.

It also provides code for calculating the conflation of PDFs, both analytical and empirical.

## Is this package for me?

This package is useful for you if:

1. You have a series of occurrences in time, and want to estimate the time at which these events started and finished happening. Examples include a series of fossil occurrences for a given species for which we want to estimate the origination and extinction times.

2. You have a collection of multiple intervals as described above, and want to calculate the probability distribution for the time of co-occurrence of them, that is, the time at which these were co-existing. For example, we have stratigraphic intervals for different species that were found at the same fossil assemblage, and we want to estimate the time of the assemblage.

3. Want to use the distributions implemented in this package for other purposes, for instance the Three-Parameter Beta distribution, or the Reflected-Offset Exponential distribution. 

## Documentation

Please visit [https://gaballench.github.io/StratIntervals.jl/](https://gaballench.github.io/StratIntervals.jl/) for documentation and additional details on what this package can do for you.

## Installation

You may install the stable version of the package from the Julia Package General Registry with:

```julia
using Pkg
Pkg.add("StratIntervals")
```

Also, you may install the development version from GitHub:

```julia
using Pkg
Pkg.add("https://github.com/gaballench/StratIntervals.jl")
```

## Authorship

This package was created and maintained by Gustavo A. Ballen
[@gaballench](https://github.com/gaballench).

## Citation

Ballen, G.A. (2025). A flexible Bayesian method for estimating stratigraphic intervals and their co-occurrence in time. BioRxiv [https://doi.org/10.1101/2025.02.13.638199](https://doi.org/10.1101/2025.02.13.638199).

```
@article {Ballen2025.02.13.638199,
	author = {Ballen, Gustavo A.},
	title = {{A flexible Bayesian method for estimating stratigraphic intervals and their co-occurrence in time}},
	elocation-id = {2025.02.13.638199},
	year = {2025},
	doi = {10.1101/2025.02.13.638199},
	URL = {https://www.biorxiv.org/content/early/2025/02/17/2025.02.13.638199},
	eprint = {https://www.biorxiv.org/content/early/2025/02/17/2025.02.13.638199.full.pdf},
	journal = {bioRxiv}
}

```

## Problems?

If you find a bug or unexpected behavior, please [file an
issue](https://github.com/gaballench/StratIntervals.jl/issues).
