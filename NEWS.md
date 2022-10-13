# goldfish.latent (development version)

# goldfish.latent 0.0.1

* Random-effects for DyNAM-choice model functionality is added:
  * `CreateData()` generates the expected data for `Stan`.
  * `CreateModelCode()` copies the `Stan` code to a temporal folder.
    `cmdstanr` uses the path to the new copied code to compile the model and
    generate sampled from the posterior distribution.
  * `ComputeLogLikelihood()` creates an array with the log-likelihood of either
    the marginal or the conditional version. Output is ready to be used by `loo`
    package.

# goldfish.latent 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
