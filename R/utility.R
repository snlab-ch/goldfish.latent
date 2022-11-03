#' Create a Stan code file
#'
#' `goldfish.latent` offers working version of the models to work with Stan.
#' Users should have installed `cmdstanr` package and `CmdStan` before.
#' Follow the instructions instructions from `cmdstanr` documentation.
#'
#' `cmdstanr` functionalities makes possible for users to write a local copy
#' of the Stan code and modify it for other purposes.
#' The default in [cmdstanr::write_stan_file()] is to write the model in a
#' temporal folder. Using the `dir` argument is possible to write the model
#' in a folder specifies by the user.
#'
#' @param dataStan a `list` output of a [CreateData()] call.
#' @param ... additional arguments to be passed to
#'   [cmdstanr::write_stan_file()]
#'
#' @return The path to a file with `stan` extension.
#' It contains the code with the specification of data structure, priors,
#' and log-likelihood of the given model.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' library(cmdstanr)
#' data("Social_Evolution")
#' callNetwork <- defineNetwork(nodes = actors, directed = TRUE) |>
#'   linkEvents(changeEvent = calls, nodes = actors)
#' callsDependent <- defineDependentEvents(
#'   events = calls, nodes = actors, defaultNetwork = callNetwork
#' )
#' data2stan <- CreateDataModel(
#'   randomEffects = list(inertia ~ 1),
#'   fixedEffects = callsDependent ~ recip + trans
#' )
#'
#' stanCode <- CreateModelCode(data2stan)
#' }
CreateModelCode <- function(dataStan, ...) {
  stopifnot(inherits(dataStan, "goldfish.latent.data"))

  model <- attr(dataStan, "model")
  subModel <- attr(dataStan, "subModel")

  if (model == "DyNAMRE" && subModel == "choice") {
    if (dataStan[["dataStan"]][["Q"]] == 1) {
      fileModel <- "MCM_RE1.stan"
    } else stop("Not yet implemented for more than one random effect")
  } else if (model == "DyNAMSR") {
    if (subModel == "both" && dataStan[["dataStan"]][["hasIntercept"]]) {
      # stop("Not yet implemented, use independent submodels")
      # fileModel <- "DyNAMSR_both.stan"
      fileModel <- "DyNAM_both.stan"
    } else if (subModel == "both") {
      stop("Not yet implemented, use independent submodels")
      fileModel <- "DyNAMSR_both_ord.stan"
    } else if (subModel == "rate" && dataStan[["dataStan"]][["hasIntercept"]]) {
      fileModel <- "DyNAM_rate.stan"
    } else if (subModel == "choice") {
      fileModel <- "DyNAMSR_choice.stan"
    } else stop("not implemented yet")
  }

  stanCode <- readLines(
    system.file("stan", fileModel, package = "goldfish.latent")
  )

  if (requireNamespace("cmdstanr", quietly = TRUE) &&
      cmdstanr::cmdstan_version() >= "2.29.2") {
    model <- cmdstanr::write_stan_file(code = stanCode)
  } else
    stop(
      dQuote("cmdstanr"), " package and a working version of",
      dQuote("CmdStan"), "are required.",
      "\nPlease follow Stan documentation for instructions on how to install.")

  return(model)
}
