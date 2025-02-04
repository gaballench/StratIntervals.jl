# Conflation: Combining distributions

Combining different distributions is not straitghforward (Genest and Zidek, 1986). However, the conflation of probability density functions is a useful procedure which combines them provided that each of them is independent (Hill, 2008, 2011).  Assuming that each distribution describing the posterior predictive distribution of each interval is independent, and simplifying notation so that ``\tilde{\tau} = \tilde{\tau} | \boldsymbol{\tau}``, we can define the composite distribution of ``\tilde{\tau}`` for ``M`` intervals as the conflation of individual posterior predictive distributions:

```math
\begin{equation}
  \begin{aligned}
    Q(\tilde{\tau}) = \frac{ \prod_{i=1}^{M} p_i(\tilde{\tau})}{\int_{-\infty}^{\infty} \prod_{i=1}^{M} p_i(x) \mathrm{d}x}
  \end{aligned}
\end{equation}
```

Such distribution is useful e.g. when we want to build credible intervals for the time ``\tau`` of co-occurrence of stratigraphic intervals. The conceptual interpretation as follows: ``M`` different lineages as represented by their stratigraphic intervals should coexist at most for some time interval when they all were alive. As the conflation of densities is a density itself, it can be used for asking questions on the probability of co-existence of lineages during some arbitrary time interval given the distribution.

The application of conflation is however more general than the specific case where we combine different posterior predictives as the distribution of the co-occurrence of stratigraphic intervals.
Ballen and Reinales (2025) have used the method for combining distributions into secondary calibrations as well as for estimating biogeographic events in divergence time estimation. As long as we have a vector of distributions, the calculation of their conflation is straightforward:

```@repl
using Distributions
using StatsPlots
using StratIntervals

dists = product_distribution([Normal(10, 2), Normal(11, 4), Normal(8, 5)])
xx = 5:0.0001:15
yy = map(x -> conflate(dists, x), xx)
plot(xx, yy, label="Conflation")
plot!(Normal(10, 2), label="Normal(10, 2)")
plot!(Normal(11, 4), label="Normal(11, 4)")
plot!(Normal(8, 5), label="Normal(8, 5)")
savefig("conflation_example.svg");

nothing
```

![](conflation_example.svg)
