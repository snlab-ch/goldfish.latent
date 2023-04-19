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
      fileModel <- "DyNAMSR_rate.stan"
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


#' sample preprocessed data
#'
#' For arguments `methodChoiceSet` and `methodEvents` the available methods are
#' `"systematic"` and `"srswor"`. They correspond to
#'   systematic sampling and simple random sampling without replacement,
#'   respectively. It is possible to use an external function that has arguments
#'   `N` and `fraction` and return a numerical vector of the samples
#'   to keep is possible.
#' @param data output from [CreateDataSR()] or [CreateData()]
#' @param fractionChoiceSet numerical value that indicates the proportion of
#'   alternatives to sample from the choice set or the compiting actors in the
#'   choice and rate model, respectively.
#' @param methodChoiceSet character value indicating the function name used to
#'   generate the sample.
#' @param fractionEvents numerical value that indicates the proportion of
#'   events to sample. When the `data` object contains data for both sub-models,
#'   a pair sample is selected, i.e., the sample events contains information of
#'   the rate and choice sub-models.
#' @param methodEvents character value indicating the function name used to
#'   generate the sample.
#'
#' @return an object with the same meta information as `data` with the sampled
#' events and alternative choices.
#' @export
#'
#' @examples
#' sampledData <- SampleData(data)
SampleData <- function(
    data,
    fractionChoiceSet = 0.1, methodChoiceSet = "srswor",
    fractionEvents = 1, methodEvents = "systematic"
) {
  stopifnot(
    inherits(data, "goldfish.latent.data"),
    is.numeric(fractionChoiceSet) &&
      length(fractionChoiceSet) == 1 &&
      fractionChoiceSet > 0 && fractionChoiceSet <= 1,
    is.numeric(fractionEvents) &&
      length(fractionEvents) == 1 &&
      fractionEvents > 0 && fractionEvents <= 1,
    is.character(methodChoiceSet) && length(methodChoiceSet) == 1,
    is.character(methodEvents) && length(methodEvents) == 1
  )

  methodChoiceSet <- match.arg(methodChoiceSet, c("srswor", "systematic"))
  methodEvents <- match.arg(methodEvents, c("systematic", "srswor"))

  dataStan <- data$dataStan
  CollapseSample <- function(x) Reduce(f = rbind, x = x)

  if (attr(data, "model") %in% c("DyNAMSR") &&
      attr(data, "subModel") %in% c("both", "rate")) {
    dataStan <- within(
      dataStan,
      {
        expandedDF <- as.data.frame(Xrate)
        expandedDF$event <- rep(seq_len(Trate), endRate - startRate + 1)
        expandedDF$selected <- seq_len(Nrate) %in% choseRate

        if (fractionEvents < 1) {
          sampleEvents <- do.call(
            methodEvents,
            list(N = Trate, fraction = fractionEvents)
          )

          if (attr(data, "subModel") == "both") {
            if (Trate != Tchoice) {
              eventsChoice <- seq_len(Tchoice)
              eventsRate <- rep(0, Trate)
              eventsRate[isDependent] <- eventsChoice
              rm(eventsChoice)
            } else eventsRate <- seq_len(Trate)

            sampleEventsChoice <- eventsRate[sampleEvents]
            sampleEventsChoice <- sampleEventsChoice[sampleEventsChoice > 0]
            rm(eventsRate)
          }

          Trate <- length(sampleEvents)
          expandedDF <- subset(expandedDF, event %in% sampleEvents)
          isDependent <- isDependent[sampleEvents]
          timespan <- timespan[sampleEvents]
          rm(sampleEvents)
        }

        if (fractionChoiceSet < 1) {
          expandedDF <- by(
            expandedDF,
            expandedDF[, c("selected", "event")],
            \(x) {
              sample <- do.call(
                methodChoiceSet,
                list(N = nrow(x), fraction = fractionChoiceSet)
              )
              x[sample, ]
            }
          ) |> CollapseSample()

        }

        Nrate <- nrow(expandedDF)
        idxEvents <- tapply(seq_len(Nrate), expandedDF$event, range) |>
          simplify2array()
        startRate <- idxEvents[1, ]
        endRate <- idxEvents[2, ]
        Xrate <- as.matrix(expandedDF[, colnames(Xrate)])
        choseRate <- rep(0, Trate)
        choseRate[isDependent] <- which(expandedDF$selected)

        rm(idxEvents, expandedDF)
      })
  }

  if (attr(data, "model") %in% c("DyNAMSR") &&
             attr(data, "subModel") %in% c("both", "choice")) {
    dataStan <- within(
      dataStan,
      {
        expandedDF <- cbind(
          as.data.frame(Xchoice),
          data.frame(
            event = rep(seq_len(Tchoice), endChoice - startChoice + 1),
            selected = seq_len(Nchoice) %in% choseChoice
          )
        )

        if (fractionEvents < 1) {

          if (attr(data, "subModel") != "both") {
            sampleEventsChoice <- do.call(
              methodEvents,
              list(N = Tchoice, fraction = fractionEvents)
            )
          }

          expandedDF <- subset(expandedDF, event %in% sampleEventsChoice)
          Tchoice <- length(sampleEventsChoice)

          rm(sampleEventsChoice)
        }

        if (fractionChoiceSet < 1) {
          expandedDF <- by(
            expandedDF,
            expandedDF[, c("selected", "event")],
            \(x) {
              sample <- do.call(
                methodChoiceSet,
                list(N = nrow(x), fraction = fractionChoiceSet)
              )
              x[sample, ]
            }
          ) |> CollapseSample()
        }

        Nchoice <- nrow(expandedDF)
        idxEvents <- tapply(seq_len(Nchoice), expandedDF$event, range) |>
          simplify2array()
        startChoice <- idxEvents[1, ]
        endChoice <- idxEvents[2, ]
        Xchoice <- as.matrix(expandedDF[, colnames(Xchoice)])
        choseChoice <- which(expandedDF$selected)

        rm(idxEvents, expandedDF)
      }
    )
  }

  data$dataStan <- dataStan
  return(data)
}


AddColumns <- function(arrayTo, arrayFrom) {
  missCols <- setdiff(
    colnames(arrayFrom), colnames(arrayTo)
  )
  arrayTo <- as.data.frame(arrayTo)
  if (length(missCols) > 0)
    for (colAdd in missCols) {
      arrayTo[, colAdd] <- ""
    }

  return(arrayTo)
}

systematic <- function(N, fraction) {
  stopifnot(
    is.numeric(N) & length(N) == 1 && N >= 1,
    is.numeric(fraction) && length(fraction) == 1 &&
      fraction > 0 && fraction <= 1
  )

  if (N == 1) return(1)
  sampleSize <- ceiling(N * fraction)

  if (sampleSize >= N) return(seq_len(N))

  jump <- floor(1 / fraction)
  r <- sample(jump, 1)
  # cUpper <- N - sampleSize * jump
  # n <- (r <= cUpper) + jump

  return(r + (jump * (seq_len(sampleSize) - 1)))
}

srswor <- function(N, fraction) {
  stopifnot(
    is.numeric(N) & length(N) == 1 && N >= 1,
    is.numeric(fraction) && length(fraction) == 1 &&
      fraction > 0 && fraction <= 1
  )

  if (N == 1) return(1)
  sampleSize <- ceiling(N * fraction)

  if (sampleSize >= N) return(seq_len(N))

  sort(sample.int(N, ceiling(N * fraction)))
}

RescaleCoefs <- function(
    beta, scaleStats, offset = 0, isRate = TRUE
) {
  beta2 <- beta ## inherit names etc.

  if (isRate) {
    mu <- scaleStats$rate$`scaled:center`
    sigma <- scaleStats$rate$`scaled:scale`

    beta2[-1] <- beta[-1] / sigma
    beta2[1]  <- beta[1] + offset - sum(beta2[-1] * mu)
  } else {
    sigma <- scaleStats$choice$`scaled:scale`
    beta2 <- beta / sigma
  }
  return(beta2)
}
