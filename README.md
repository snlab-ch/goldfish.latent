
<!-- README.md is generated from README.Rmd. Please edit that file -->

# goldfish.latent

<!-- badges: start -->
<!-- badges: end -->

`goldfish.latent` extends the `goldfish` package with models that
include latent variables. Therefore, how to define the event sequence
data objects for modeling and the effects available for modeling uses
`goldfish`’s data objects definitions and statistics effects. Inferences
of the proposed models are made using Hamiltonian Chain Monte Carlo as
implemented in [Stan](https://mc-stan.org/).

## Installation

You can install the development version of goldfish.latent like so:

``` r
remotes::install_github("snlab-ch/goldfish.latent", build_vignettes = TRUE)
```

## Example

The first model introduced in `goldfish.latent` is the DyNAM with random
effects. Model formulation allows to include monadic statistic or
actors’ covariates to explain the variability of the random effects.

``` r
library(goldfish.latent)

# using cmdstanr for getting posterior samples and 
library(cmdstanr)

# using goldfish social evolution data set
library(goldfish)
data("Social_Evolution")
callNetwork <- defineNetwork(nodes = actors, directed = TRUE) |>
 linkEvents(changeEvent = calls, nodes = actors)
callsDependent <- defineDependentEvents(
 events = calls, nodes = actors, defaultNetwork = callNetwork
)
data2stan <- CreateData(
 randomEffects = list(inertia ~ 1),
 fixedEffects = callsDependent ~ recip + trans
)

stanCode <- CreateModelCode(data2stan)

mod01 <- cmdstan_model(stanCode)
mod01Samples <- mod01$sample(
  data = data2stan[["dataStan"]],
  parallel_chains = 4, chains = 4,  iter_warmup = 500, iter_sampling = 500,
  show_messages = FALSE
)
```

Using `cmdstanr` functionalities makes easier fit to use summary and
plotting function from packages that works with MCMC posterior samples.
For the summary of the posterior samples `cmdstanr` uses the
[`posterior`](https://mc-stan.org/posterior/) package. Plots from the
posterior samples can be done using, for example,
[`bayesplot`](https://mc-stan.org/bayesplot/) package.

``` r
mod01Samples$summary("beta")
#> # A tibble: 3 × 10
#>   variable   mean median    sd   mad     q5   q95  rhat ess_bulk ess_tail
#>   <chr>     <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 beta[1]   4.94   4.93  0.400 0.377  4.29  5.61   1.00    1107.    1433.
#> 2 beta[2]   1.64   1.63  0.206 0.210  1.31  1.97   1.00    2962.    1577.
#> 3 beta[3]  -0.275 -0.277 0.231 0.237 -0.656 0.101  1.00    3231.    1561.

  # names fixed effects coincides with colnames(data2stan$dataStan$X)
```

The `ComputeLogLikelihood()` function allows computing the marginal or
conditional log-likelihood for the model using MCMC samples from the
posterior distribution. Parallel computation using `parallel` package is
possible. Using the [`loo`](https://mc-stan.org/loo/) functionalities is
possible to compute the Watanabe Information Criterion `waic()` and the
Leave-One-Out approximation using Pareto smoothed importance sampling
`loo()`. Those information criteria can be used to compare models using
`loo_compare()`.

``` r
logLikMod01 <- ComputeLogLikelihood(mod01Samples, data2stan, spec = 4)

# loo and waic computation using the margina
library(loo)

relEff <- relative_eff(exp(logLikMod01))
looM01 <- loo(logLikMod01, r_eff = relEff)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
looM01
#> 
#> Computed from 2000 by 34 log-likelihood matrix
#> 
#>          Estimate    SE
#> elpd_loo   -679.9 112.1
#> p_loo        17.7   8.0
#> looic      1359.8 224.1
#> ------
#> Monte Carlo SE of elpd_loo is NA.
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. n_eff
#> (-Inf, 0.5]   (good)     27    79.4%   711       
#>  (0.5, 0.7]   (ok)        2     5.9%   199       
#>    (0.7, 1]   (bad)       0     0.0%   <NA>      
#>    (1, Inf)   (very bad)  5    14.7%   14        
#> See help('pareto-k-diagnostic') for details.
```

## Code of Conduct

Please note that the `goldfish.latent` project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
