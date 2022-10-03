
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
#> 1 beta[1]   4.92   4.90  0.420 0.409  4.24  5.61   1.00     935.    1001.
#> 2 beta[2]   1.64   1.64  0.199 0.197  1.31  1.97   1.00    3588.    1525.
#> 3 beta[3]  -0.275 -0.272 0.232 0.230 -0.667 0.103  1.00    3634.    1726.

# names fixed effects coincides with colnames(data2stan$dataStan$X)

# Plots from the posterior distribution can be generated using `bayesplot`
# or any other package that handle MCMC posteriors
```

## Code of Conduct

Please note that the `goldfish.latent` project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
