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
#' @param data_stan a `list` output of a [CreateData()] call.
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
#' callNetwork <- make_network(nodes = actors, directed = TRUE) |>
#'   link_events(change_event = calls, nodes = actors)
#' callsDependent <- define_dependent_events(
#'   events = calls, nodes = actors, default_network = callNetwork
#' )
#' socialEvolutionData <- make_data(callsDependent)
#' data2stan <- make_data_re(
#'   random_effects = list(inertia ~ 1),
#'   fixed_effects = callsDependent ~ recip + trans,
#'   data = data
#' )
#'
#' stanCode <- make_model_code(data2stan)
#' }
# make_model_code <- function(data_stan, ...) {
#   UseMethod("make_model_code", data_stan)
# }

make_model_code <- function(data_stan, ...) {
  stopifnot(inherits(data_stan, "goldfish.latent.data"))

  model <- attr(data_stan, "model")
  subModel <- attr(data_stan, "subModel")

  if (model == "DyNAMRE" && subModel == "choice") {
    if (data_stan[["dataStan"]][["Qchoice"]] == 1) {
      fileModel <- "DNRE1_choice.stan"
    } else stop("Not yet implemented for more than one random effect")
  } else if (model == "DNHMM") {
    if (subModel == "both" && data_stan[["dataStan"]][["hasIntercept"]]) {
      # stop("Not yet implemented, use independent submodels")
      # fileModel <- "DyNAMSR_both.stan"
      fileModel <- "DNHMM_both.stan"
    } else if (subModel == "both") {
      stop("Not yet implemented, use independent submodels")
      fileModel <- "DNHMM_both_ord.stan"
    } else if (subModel == "rate" &&
               data_stan[["dataStan"]][["hasIntercept"]]) {
      fileModel <- "DNHMM_rate.stan"
    } else if (subModel == "choice") {
      fileModel <- "DNHMM_choice.stan"
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


# @rdname make_model_code
# @export
# make_model_code.DNRE <- function(
#     data_stan,
#     sub_model = attr(data_stan, "subModel"),
#     prior = c("normal", "t-student"),
#     prior_re = c("default", "gamma", "invWishart", "LKJ"),
#     generate_quantities = FALSE,
#     ...
# ) {
#   prior <- match.arg(prior)
#   priorRE <- match.arg(priorRE)
#   subModel <- match.arg(subModel, c("both", "choice", "rate"))
# 
#   Q <- ifelse(
#     !is.null(dataStan[["stan"]][["Qrate"]]),
#     dataStan[["stan"]][["Qrate"]],
#     0L
#   ) +
#     ifelse(
#       !is.null(dataStan[["stan"]][["Qchoice"]]),
#       dataStan[["stan"]][["Qchoice"]],
#       0L
#     )
# 
#   typeQ <- ifelse(Q > 1, "Qm", "Q1")
# 
#   parmsFE <- do.call(
#     utils::getS3method("ReadParms", subModel),
#     list(x = prior)
#   )
#   parmsRE <- c(
#     SplitJoinChunk(subModel, "DNRE", paste0("parms", typeQ), isCommon = TRUE),
#     SplitJoinChunk(subModel, "DNRE", paste0("parms", typeQ))
#   )
#   stanCode <- c(
#     "data {",
#     ReadStanChunkType(subModel, "DN", "data"),
#     "  int A; // number of senders",
#     ReadStanChunkType(subModel, "DNRE", "data"),
#     "}\nparameters {",
#     parmsFE[["pm"]],
#     parmsRE[["pm"]],
#     "}",
#     parmsRE[["tp"]],
#     "model {\n  //priors",
#     parmsFE[["pr"]],
#     parmsRE[["pr"]],
#     "  // loglikelihood",
#     parmsRE[["ll"]],
#     "}",
#     if (generateQuantities) parmsRE[["gq"]]
#   )
# 
#   cmdstanr::write_stan_file(code = stanCode)
# }
# 
# ReadStanChunk <- function(prefix, type, suffix, isType = TRUE) {
#   fileChunk <- paste0(
#     prefix, "_", if (isType) paste0(type, "_"), suffix, ".stan"
#   )
#   readLines(system.file("stan", fileChunk, package = "goldfish.latent"))
# }
# 
# ReadStanChunkType <- function(
#     x, prefix, suffix, isCommon = FALSE
# ) {
#   switch(
#     (!isCommon) * match(x, c("rate", "choice", "both")) + 1,
#     ReadStanChunk(prefix, "", suffix, isType = FALSE),
#     ReadStanChunk(prefix, "rt", suffix),
#     ReadStanChunk(prefix, "ch", suffix),
#     c(
#       ReadStanChunkType("rate", prefix, suffix, isCommon),
#       ReadStanChunkType("choice",  prefix, suffix, isCommon)
#     ),
#     stop("not recognize submodel")
#   )
# }
# 
# ReadStanChunkLL <- function(x, prefix) {
#   suffix <- "ll"
#   switch(
#     match(x, c("rate", "choice", "both")),
#     ReadStanChunk(prefix, "rt", suffix),
#     ReadStanChunk(prefix, "ch", suffix),
#     ReadStanChunk(prefix, "bt", suffix),
#     stop("not recognize submodel")
#   )
# }
# SplitJoinChunk <- function(
#     x, prefix, suffix, isCommon = FALSE
# ) {
#   codeLines <- ReadStanChunkType(x, prefix, suffix, isCommon)
#   chunksTitles <- grep("^\\h*// (\\w+)$", codeLines)
#   codeOrg <- list()
#   positions <- c(chunksTitles, length(codeLines))
#   for (ch in seq_along(chunksTitles)) {
#     title <- gsub("^\\h*// (\\w+)$", "\\1", codeLines[positions[ch]])
#     codeOrg[[title]] <- c(
#       codeOrg[[title]],
#       codeLines[seq.int(positions[ch] + 1, positions[ch + 1] - 1)]
#     )
#   }
#   return(codeOrg)
# }
# 
# ReadParms <- function(x, ...) {
#   UseMethod("ReadParms", x)
# }

# @noRd
# @export
# ReadParms.rate <- function(x) {
#   list(
#     parms = "  vector[Prate] betaRate;",
#     prior = switch(
#       x,
#       "normal" = c(
#         "  target += normal_lpdf(betaRate[1] | 0, 10);",
#         "  target += std_normal_lpdf(betaRate[2: ]);"
#       ),
#       "t-student" = c(
#         "  target += student_t_lpdf(betaRate[1] | 3, 0, 10);",
#         "  target += student_t_lpdf(betaRate[2: ] | 3, 0, 1);"
#       )
#     )
#   )
# }

# @noRd
# @export
# ReadParms.choice <- function(x) {
#   list(
#     parms = "  vector[Pchoice] betaChoice;",
#     prior = switch(
#       x,
#       "normal" = "  target += std_normal_lpdf(betaChoice);",
#       "t-student" = "  target += student_t_lpdf(betaChoice | 3, 0, 1);"
#     )
#   )
# }

# @noRd
# @export
# ReadParms.both <- function(x) {
#   rate <- ReadParms.rate(x)
#   choice <- ReadParms.choice(x)
#   list(
#     parms = c(rate[["parms"]], choice[["parms"]]),
#     prior = c(rate[["prior"]], choice[["prior"]])
#   )
# }



#' sample preprocessed data
#'
#' For arguments `methodChoiceSet` and `methodEvents` the available methods are
#' `"systematic"` and `"srswor"`. They correspond to
#'   systematic sampling and simple random sampling without replacement,
#'   respectively. It is possible to use an external function that has arguments
#'   `N` and `fraction` and return a numerical vector of the samples
#'   to keep.
#' @param data output from [make_data_hmm()] or [make_data_re()]
#' @param fraction_choice_set numerical value that indicates the proportion of
#'   alternatives to sample from the choice set or the compiting actors in the
#'   choice and rate model, respectively.
#' @param method_choice_set character value indicating the function name used to
#'   generate the sample.
#' @param fraction_events numerical value that indicates the proportion of
#'   events to sample. When the `data` object contains data for both sub-models,
#'   a pair sample is selected, i.e., the sample events contains information of
#'   the rate and choice sub-models.
#' @param method_events character value indicating the function name used to
#'   generate the sample.
#'
#' @return an object with the same meta information as `data` with the sampled
#' events and alternative choices.
#' @export
#'
#' @examples
#' sampledData <- sample_data(data)
sample_data <- function(
    data,
    fraction_choice_set = 0.1, method_choice_set = c("srswor", "systematic"),
    fraction_events = 1, method_events = c("systematic", "srswor")
) {
  
  method_choice_set <- match.arg(method_choice_set, c("srswor", "systematic"))
  method_events <- match.arg(method_events, c("systematic", "srswor"))
  
  stopifnot(
    inherits(data, "goldfish.latent.data"),
    is.numeric(fraction_choice_set) &&
      length(fraction_choice_set) == 1 &&
      fraction_choice_set > 0 && fraction_choice_set <= 1,
    is.numeric(fraction_events) &&
      length(fraction_events) == 1 &&
      fraction_events > 0 && fraction_events <= 1,
    is.character(method_choice_set) && length(method_choice_set) == 1,
    is.character(method_events) && length(method_events) == 1
  )

  dataStan <- data$dataStan
  model <- attr(data, "model")
  subModel <- attr(data, "subModel")
  CollapseSample <- function(x) Reduce(f = rbind, x = x)

  if (model %in% c("DNHMM", "DNCHMM", "DyNAMRE") &&
      subModel %in% c("both", "rate")) {
    dataStan <- within(
      dataStan,
      {
        expandedDF <- as.data.frame(Xrate)
        expandedDF$event <- rep(seq_len(Trate), endRate - startRate + 1)
        expandedDF$selected <- seq_len(Nrate) %in% choseRate

        if (fraction_events < 1) {
          sampleEvents <- do.call(
            method_events,
            list(N = Trate, fraction = fraction_events)
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

        if (fraction_choice_set < 1) {
          expandedDF <- by(
            expandedDF,
            expandedDF[, c("selected", "event")],
            \(x) {
              sample <- do.call(
                method_choice_set,
                list(N = nrow(x), fraction = fraction_choice_set)
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

  if (model %in% c("DNHMM", "DNCHMM", "DyNAMRE") &&
      subModel %in% c("both", "choice")) {
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

        if (fraction_events < 1) {

          if (attr(data, "subModel") != "both") {
            sampleEventsChoice <- do.call(
              method_events,
              list(N = Tchoice, fraction = fraction_events)
            )
          }

          expandedDF <- subset(expandedDF, event %in% sampleEventsChoice)
          Tchoice <- length(sampleEventsChoice)

          rm(sampleEventsChoice)
        }

        if (fraction_choice_set < 1) {
          expandedDF <- by(
            expandedDF,
            expandedDF[, c("selected", "event")],
            \(x) {
              sample <- do.call(
                method_choice_set,
                list(N = nrow(x), fraction = fraction_choice_set)
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

ChangeKR <- function(chrVec, kR) {
  whichLine <- grepl("%%", chrVec)
  chrVec[whichLine] <- gsub("%%kR%%", as.character(kR), chrVec[whichLine])
  return(chrVec)
}
