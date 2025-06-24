
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
set_cmdstan_path()
cmdstan_version()

# using goldfish social evolution data set
library(goldfish)
data("Social_Evolution")
callNetwork <- make_network(nodes = actors, directed = TRUE) |>
 link_events(change_event = calls, nodes = actors)
callsDependent <- make_dependent_events(
 events = calls, nodes = actors, default_network = callNetwork
)
socialEvolutionData <- make_data(callsDependent)
data2stan <- make_data_re(
 random_effects = list(inertia ~ 1),
 fixed_effects = callsDependent ~ recip + trans,
 data = socialEvolutionData
)

stanCode <- make_model_code(data2stan)

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
mod01Samples$summary("betaChoice")
#> # A tibble: 3 × 10
#>   variable        mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
#>   <chr>          <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#> 1 betaChoice[1]  4.30   4.33  0.406 0.385  3.63  4.92    1.00     842.    1126.
#> 2 betaChoice[2]  1.65   1.65  0.200 0.203  1.33  1.98    1.00    2771.    1496.
#> 3 betaChoice[3] -0.251 -0.241 0.217 0.218 -0.622 0.0853  1.00    3489.    1445.

  # names fixed effects coincides with colnames(data2stan$dataStan$X)
```

The `compute_log_likelihood()` function allows computing the marginal or
conditional log-likelihood for the model using MCMC samples from the
posterior distribution. Parallel computation using `parallel` package is
possible. Using the [`loo`](https://mc-stan.org/loo/) functionalities is
possible to compute the Watanabe Information Criterion `waic()` and the
Leave-One-Out approximation using Pareto smoothed importance sampling
`loo()`. Those information criteria can be used to compare models using
`loo_compare()`.

``` r
logLikMod01 <- compute_log_likelihood(mod01Samples, data2stan, spec = 4)

# loo and waic computation using the marginal
library(loo)
#> This is loo version 2.8.0
#> - Online documentation and vignettes at mc-stan.org/loo
#> - As of v2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session.

relEff <- relative_eff(exp(logLikMod01))
looM01 <- loo(logLikMod01, r_eff = relEff)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
looM01
#> 
#> Computed from 2000 by 34 log-likelihood matrix.
#> 
#>          Estimate    SE
#> elpd_loo   -681.6 111.8
#> p_loo        18.1   8.2
#> looic      1363.3 223.5
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 1.3]).
#> 
#> Pareto k diagnostic values:
#>                          Count Pct.    Min. ESS
#> (-Inf, 0.7]   (good)     28    82.4%   155     
#>    (0.7, 1]   (bad)       2     5.9%   <NA>    
#>    (1, Inf)   (very bad)  4    11.8%   <NA>    
#> See help('pareto-k-diagnostic') for details.
```

## Code of Conduct

Please note that the `goldfish.latent` project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
