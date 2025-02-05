```@meta
CurrentModule = StratIntervals
```

# StratIntervals.jl

Documentation for the [StratIntervals](https://github.com/gaballench/StratIntervals.jl) package.
Please refer to Ballen (XXXX) for details as well as to the documentations herein.

## Is this package for me?

This package is useful for you if:

1. You have a series of occurrences in time, and want to estimate the time at which these events started and finished happening. Examples include a series of fossil occurrences for a given species for which we want to estimate the origination and extinction times.

2. You have a collection of multiple intervals as described above, and want to calculate the probability distribution for the time of oc-occurrence of them, that is, the time at which these were co-existing. For example, we have stratigraphic intervals for different species that were found at the same fossil assemblage, and we want to estimate the time of the assemblage.

3. Want to use the distributions implemented in this package for other purposes, for instance the Three-Parameter Beta distribution, or the Reflected-Offset Exponential distribution. 

# Literature cited

The following references have been cited throughout the documentation.

Ballen, G.A. (XXXX). A flexible Bayesian method for estimating stratigraphic intervals and
their co-occurrence in time. XXXX.

Ballen, G.A. & Reinales, S. (2025). tbea: tools for pre- and post-processing in Bayesian evolutionary analyses. Submitted to Evolutionary Journal of the Linnean Society. Preprint at BioRxiv https://www.biorxiv.org/content/10.1101/2024.06.18.599561.

Genest, C. and Zidek, J. V. (1986). Combining probability pistributions: A critique and an annotated bibliography.
Statistical Science, 1(1):114–135.

Hill, T. P. (2011). Conflations of probability distributions. Transactions of the American Mathematical Society,
363(6):3351–3372.

Hill, T. P. and Miller, J. (2011). How to combine independent data sets for the same quantity. Chaos: An
Interdisciplinary Journal of Nonlinear Science, 21(3):033102.
