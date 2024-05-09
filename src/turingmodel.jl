#using Turing, StatsPlots, Random, Distributions

#include("BetaAdaptive.jl")

# we want to check Turing by simulating data from known params and then sampling their posterior distributions
# set vars and true param values
ndata = 100
true_lambda = 0
true_theta1 = 0.5
true_theta2 = 22

data = rand(BetaAdaptive(true_theta1, true_theta2, true_lambda), ndata)

# then declare the Turing model as a hierarchical Bayesian model with
# priors on theta1, theta2, and lambda
@model function stratint(tau::Vector{Real})
    # priors
    theta1 ~ Normal(0.5,1)
    theta2 ~ Normal(22,1)
    lambda ~ Normal(0,1)
    # reject values outside the support of the parameters
    if theta1 > minimum(data) || theta2 < maximum(data) || theta1 < 0 || theta2 < 0
        Turing.@addlogprob! -Inf
        return nothing
    end
    # define the likelihood of the observed data
    N = length(data)
    for i in 1:N
        tau[i] ~ BetaAdaptive(theta1, theta2, lambda)
        #println(tau[i]~ BetaAdaptive(theta1, theta2, lambda))
    end
end

# settings for the MCMC sampler
iters = 50000
epsilon = 0.01
tau_hmc = 5

chain_prior = sample(stratint(convert(Vector{Real}, data)), Prior(), iters)
chain = sample(stratint(convert(Vector{Real}, data)), HMC(epsilon, tau_hmc), iters)
#println(summary(chain_prior))
#println(summary(chain))
plot(chain)
histogram(data)

########## defining the loglik as a product of distributions is NOT straightforward, BUT
# we can define a vector of X dists vector_dists =  [Dist1(), Dist2, ..., DistN()] and then
# sum(logpdf.(vector_dists, value)) will be the loglik as the sum of the value evaluated on each dist X  
# this allows to ignore the awkward behavour of both product_distribution and Product with seem to produce
# objects for Turing models which are sampleable but don't produce the actual product of the dists.
# I denfinitely need to look into these two functions better

# this is not working
plot(product_distribution([BetaAdaptive(0, 1, 0), BetaAdaptive(0.5, 1.5, 0)]))
# this is not working
plot(product_distribution([Normal(0, 1), Normal(1, 1)]))

############# checking what happens when n is large so that the loglik might run into issues with under-overflow in floating point arithmetics
###########3# the code below allows to reproduce the issue with big samples with a density other than my BetaAdaptive

# Set the true probability of heads in a coin.
mu = 0.0
sigma = 1.0
# Iterate from having seen 0 observations to 100 observations.
Ns = 1000000

# Draw data from a Bernoulli distribution, i.e. draw heads or tails.
Random.seed!(12)
largedata = rand(Normal(mu, sigma), last(Ns))
smalldata = rand(Normal(mu, sigma), 100)

# Declare our Turing model.
@model function normals(y)
    # Our prior belief about the probability of heads in a coin.
    sigmasq ~ InverseGamma(2, 1)
    mu ~ Normal(0, sqrt(sigmasq))
    #y ~ Normal(mu, sqrt(sigmasq))

    # The number of observations.
    N = length(y)
    for n in 1:N
        # Heads or tails of a coin are drawn from a Bernoulli distribution.
        y[n] ~ Normal(mu, sigma)
    end
end

# Settings of the Hamiltonian Monte Carlo (HMC) sampler.
iterations = 1000
ϵ = 0.05
τ = 10

# Start sampling.
chain_normals_small = sample(normals(smalldata), HMC(ϵ, τ), iterations)
chain_normals_large = sample(normals(largedata), HMC(ϵ, τ), iterations)
#plot(chain_coin)

# plot the sample and lik surface
#ENV["GKS_ENCODING"] = "utf-8" # Allows the use of unicode characters in Plots.jl
#using Plots
#using StatsPlots
#using Turing
#using Random
using Bijectors

# Set a seed.
Random.seed!(0)

# Define a strange model.
@model function gdemo(x)
    s² ~ InverseGamma(2, 1)
    m ~ Normal(0, sqrt(s²))
    for i in eachindex(x)
        x[i] ~ Normal(m, sqrt(s²))
    end
    return s², m
end

# Define our data points.
x = [1.5, 2.0, 13.0, 2.1, 0.0]

# Set up the model call, sample from the prior.
model = gdemo(largedata)

# Evaluate surface at coordinates.
evaluate(m1, m2) = logjoint(model, (m=m2, s²=invlink.(Ref(InverseGamma(2, 1)), m1)))

function plot_sampler(chain; label="")
    # Extract values from chain.
    val = get(chain, [:s², :m, :lp])
    ss = link.(Ref(InverseGamma(2, 1)), val.s²)
    ms = val.m
    lps = val.lp

    # How many surface points to sample.
    granularity = 100

    # Range start/stop points.
    spread = 0.5
    σ_start = minimum(ss) - spread * std(ss)
    σ_stop = maximum(ss) + spread * std(ss)
    μ_start = minimum(ms) - spread * std(ms)
    μ_stop = maximum(ms) + spread * std(ms)
    σ_rng = collect(range(σ_start; stop=σ_stop, length=granularity))
    μ_rng = collect(range(μ_start; stop=μ_stop, length=granularity))

    # Make surface plot.
    p = surface(
        σ_rng,
        μ_rng,
        evaluate;
        camera=(30, 65),
        #   ticks=nothing,
        colorbar=false,
        color=:inferno,
        title=label,
    )

    line_range = 1:length(ms)

    scatter3d!(
        ss[line_range],
        ms[line_range],
        lps[line_range];
        mc=:viridis,
        marker_z=collect(line_range),
        msw=0,
        legend=false,
        colorbar=false,
        alpha=0.5,
        xlabel="σ",
        ylabel="μ",
        zlabel="Log probability",
        title=label,
    )

    return p
end;

c = sample(model, Gibbs(HMC(0.01, 5, :s²), PG(20, :m)), 1000)
plot_sampler(c)
