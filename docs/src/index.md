```@meta
CurrentModule = StratIntervals
```

# StratIntervals.jl

Documentation for the [StratIntervals](https://github.com/gaballench/StratIntervals.jl) package.
Please refer to Ballen (2025) for details as well as to the documentations herein.

## Is this package for me?

This package is useful for you if:

1. You have a series of occurrences in time, and want to estimate the time at which these events started and finished happening. Examples include a series of fossil occurrences for a given species for which we want to estimate the origination and extinction times.

2. You have a collection of multiple intervals as described above, and want to calculate the probability distribution for the time of co-occurrence of them, that is, the time at which these were co-existing. For example, we have stratigraphic intervals for different species that were found at the same fossil assemblage, and we want to estimate the time of the assemblage.

3. Want to use the distributions implemented in this package for other purposes, for instance the Three-Parameter Beta distribution, or the Reflected-Offset Exponential distribution. 

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

# Literature cited

The following references have been cited throughout the documentation.

Ballen, G.A. (2025). A flexible Bayesian method for estimating stratigraphic intervals and their co-occurrence in time. BioRxiv [https://doi.org/10.1101/2025.02.13.638199](https://doi.org/10.1101/2025.02.13.638199).

Ballen, G.A. & Reinales, S. (2025). tbea: tools for pre- and post-processing in Bayesian evolutionary analyses. Submitted to Evolutionary Journal of the Linnean Society. Preprint at BioRxiv [https://doi.org/10.1101/2024.06.18.599561](https://doi.org/10.1101/2024.06.18.599561).

Genest, C. and Zidek, J. V. (1986). Combining probability distributions: A critique and an annotated bibliography.
Statistical Science, 1(1):114–135.

Hill, T. P. (2011). Conflations of probability distributions. Transactions of the American Mathematical Society,
363(6):3351–3372.

Hill, T. P. and Miller, J. (2011). How to combine independent data sets for the same quantity. Chaos: An
Interdisciplinary Journal of Nonlinear Science, 21(3):033102.
