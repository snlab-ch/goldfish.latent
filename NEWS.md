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

# goldfish.latent 0.0.2

* `scale` parameter in `CreateData()` allows to standardize variables and
  keep mean and standard deviation information.
* Internal function `RescaleCoef()` uses the standardization information to
  transform parameter values in the original variable scale.
* HMM-DyNAM includes a new Stan variant that allows to change temporal resolution,
  for example, the model considers changes the Hidden Markov happening after each
  day or other time window predefined.

