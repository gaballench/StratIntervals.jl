# Sampling stratigraphic intervals, posterior predictive distribution, and the time of collections of stratigraphic intervals

This section describes how to sample from the model of stratigraphic interval once we have defined the data and priors.
What comes next is to use MCMC sampling for approximating the posterior distribution of parameters.

## Going Bayesian for stratigraphic intervals

We can define the posterior distribution of parameters describing a stratigraphic interval as follows:

```math
\begin{equation}
  \begin{aligned}
    f(\theta_1,\theta_2,\lambda | \boldsymbol{\tau}) = \frac{f(\theta_1) f(\theta_2) f(\lambda) \prod_{i = 1}^{N} f(\tau_i | \theta_1,\theta_2,\lambda)}{\int f(\theta_1) f(\theta_2) f(\lambda) \prod_{i = 1}^{N} f(\tau_i | \theta_1,\theta_2,\lambda) \mathrm{d}\theta_1 \mathrm{d}\theta_2 \mathrm{d}\lambda}
  \end{aligned}
\end{equation}
```

where ``f(\theta_1)``, ``f(\theta_2)``, and ``f(\lambda)`` are the priors of parameters, and ``\prod_{i = 1}^{N} f(\tau_i | \theta_1,\theta_2,\lambda)`` is the likelihood function, where data are fixed. Here, the likelihood function is the `ThreeParBeta` distribution.

Now we can use the function `sample_stratinterval` which takes an object of type `StratInterval` as input as well as settings for the MCMC sampling.
    
It can also define whether we want to sample from the prior, and whether we want to calculate the posterior predictive distribution.

Consider the following example, our vector of occurrences in time is `[2.0, 3.1, 3.2, 4.6, 6.77]`, and we want to use a Normal prior for ``\theta_1``, an Exponential for ``\theta_2``, and again a Normal prior for ``\lambda``:

```@repl
using Distributions
using Turing
using StratIntervals
using StatsPlots
using Random

Random.seed!(119)

setprogress!(false)

data = [2.0, 3.1, 3.2, 4.6, 6.77]
θ1_prior = Normal(10, 2)
θ2_prior = Exponential(1)
λ_prior = Normal(0, 1)
my_strinterval = StratInterval(data, θ1_prior, θ2_prior, λ_prior)
post_samples = sample_stratinterval(my_strinterval, 10000, NUTS(), false, false) # sample from posterior

plot(post_samples)
savefig("post_samples.svg");

nothing
```

![](post_samples.svg)

We note a few things about the analysis: First, it took quite little time to run (but this depends on processor architechture, so it's not generalisable). Second, by inspecting the effective sample sizes (ESSs), both bulk and tail, we see they are quite high, (say above 600), and therefore we can suspect we are reaching convergence.

The MCMC sampling, which is done via [Turing.jl](https://turinglang.org/), returns the samples, provides a nice summary of the sampling, and also allows to actually plot the chains and marginal distributions. From the summary we can also read the posterior mean and quantiles.

[Turing.jl](https://turinglang.org/) is a very powerful probabilistic language for Julia and allows to do many things with posterior samples, so make sure to take a look at their documentation for further details.

## Sampling from the prior

We can sample from the priors, which means that during MCMC sampling the value of the likelihood function is always 1.0. This way we can build an unconditional distribution to make sure it represents what we set as prior. We can sample from the prior setting the argument `prior = true` in `sample_stratinterval`:

```@repl
using Distributions
using Turing
using StratIntervals
using StatsPlots
using Random

Random.seed!(14)

setprogress!(false)

data = [2.0, 3.1, 3.2, 4.6, 6.77]
θ1_prior = Normal(10, 2)
θ2_prior = Exponential(1)
λ_prior = Normal(0, 1)
my_strinterval = StratInterval(data, θ1_prior, θ2_prior, λ_prior)
prior_samples = sample_stratinterval(my_strinterval, 10000, NUTS(), true, false) # sample from prior

plot(prior_samples)
savefig("prior_samples.svg");

nothing
```

![](post_samples.svg)

## Calculating the posterior predictive

Now that we have the samples from the posterior distribution of each parameter, we can integrate over them to build a distribution for the time ``\tau``, that is, the probability of observing occurrences through time integrating over uncertainty in parameter estimation. Such distribution is the posterior predictive.

We can estimate the posterior distributions and the posterior predictive in one shot, by just setting the argument `postpredict = true` as follows:

```@repl
using Distributions
using Turing
using StratIntervals
using StatsPlots
using LaTeXStrings
using Random

Random.seed!(14)

setprogress!(false)

data = [2.0, 3.1, 3.2, 4.6, 6.77]
θ1_prior = Normal(10, 2)
θ2_prior = Exponential(1)
λ_prior = Normal(0, 1)
my_strinterval = StratInterval(data, θ1_prior, θ2_prior, λ_prior)
post_predictive = sample_stratinterval(my_strinterval, 10000, NUTS(), false, true) # sample from the posterior predictive

density(post_predictive[2], label=L"$\tilde{\tau}$")
savefig("postpredict.svg");

nothing
```

![](postpredict.svg)

## Helping reach convergence: Suggestions for effective MCMC

Convergence is a very challenging aspect of Bayesian inference. In theory we can only claim that our approximation of the posterior distributions converged to the target distribution in very specific cases and as sample size grows to infinity. This is of course of little (if any) practical value, as we never can sample to infinity.

I like to use the no-U-turn sampler (`NUTS`), which is a very effective variant of the Hamiltonian Monte Carlo. Both these are very useful as all of our parameters are continuous. However, we can use any other sampler available from `Turing.jl` such as the Hamiltonian Monte Carlo (`HMC`), Gibbs (`Gibbs`), random-walk Metropolis-Hastings (`MH`), etc. We can try to use different samplers which may show different sampling efficiency in different situations and models. More efficient samplers mean better convergence.

Another way of rising the ESS values is by letting the chains run by longer. Consider the following example, and compare it with our first run above:

```@repl
using Distributions
using Turing
using StratIntervals
using StatsPlots
using Random

Random.seed!(44)

setprogress!(false)

data = [2.0, 3.1, 3.2, 4.6, 6.77]
θ1_prior = Normal(10, 2)
θ2_prior = Exponential(1)
λ_prior = Normal(0, 1)
my_strinterval = StratInterval(data, θ1_prior, θ2_prior, λ_prior)
post_samples = sample_stratinterval(my_strinterval, 1000, NUTS(), false, false) # sample from posterior

plot(post_samples)
savefig("lowess_samples.svg");

nothing
```

![](lowess_samples.svg)

As we can see, the ESS values are quite low compared to our original ones, and the trace plots show less "random" behaviour of the sampled parameter values. These two are signs of poor convergence. Just letting the sampling to run by longer had the effect of improving convergence. 

## Calculating ``\tau`` for a collection of stratigraphic intervals

Now that we can calculate posterior distributions and posterior predictives, let us combine these tools for collections of stratigraphic intervals to build a distribution of the co-occurrence time.

In this setting we have two stratigraphic intervals, each of which is a vector of time occurrences. We will set priors on each group of parameters, and sample from the posteriors. Then we will calculate the posterior predictive of each interval and finally combine them using conflation.

```@repl
using Distributions
using Turing
using StratIntervals
using StatsPlots
using Random

Random.seed!(67)

setprogress!(false)

ndata = 100;
iters = 10000;

# first interval
true_lambda_1 = 0;
true_theta2_1 = 10;
true_theta1_1 = 22;
data_1 = rand(ThreeParBeta(true_theta1_1, true_theta2_1, true_lambda_1), ndata);

# second interval
true_lambda_2 = 1;
true_theta2_2 = 15;
true_theta1_2 = 30;
data_2 = rand(ThreeParBeta(true_theta1_2, true_theta2_2, true_lambda_2), ndata);

# specifying each StratInterval object
interval_1 = StratInterval(data_1, Normal(22, 3), Exponential(10), Normal(0,3));
interval_2 = StratInterval(data_2, Normal(30, 3), Exponential(15), Normal(1,3));

# construct the vector of StratIntervals
vecinterval = [interval_1, interval_2];

# MCMC sampling and posterior predictive
mystratint_postpredict_vec = sample_stratinterval(vecinterval, iters, NUTS(), false, true);

# now calculate the conflation of posterior predictives for the collection interval_1, interval_2

xx = 0:0.01:40;
yy = map(x -> tau_collection(mystratint_postpredict_vec, x), xx);

plot(xx, yy, label="Conflation")
density!(mystratint_postpredict_vec[1][2], label="Postpredictive 1")
density!(mystratint_postpredict_vec[2][2], label="Postpredictive 2")

savefig("tau_collection.svg");

nothing
```

![](tau_collection.svg)
