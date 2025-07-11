# Copyright (C) 2025, Alvaro Uzaheta - SNlab-ETH Zurich
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT
# License for more details.
#
# You should have received a copy of the MIT License along with this
# program. If not, see <https://opensource.org/licenses/MIT>.

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
  
  if (model == "DN_RE" && subModel == "choice") {
    if (data_stan[["data_stan"]][["Qchoice"]] == 1) {
      fileModel <- "DNRE1_choice.stan"
    } else {
      stop("Not yet implemented for more than one random effect")
    }
  } else if (model == "DNHMM") {
    if (subModel == "both" && data_stan[["data_stan"]][["hasIntercept"]]) {
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
    } else {
      stop("not implemented yet")
    }
  }
  
  stanCode <- readLines(
    system.file("stan", fileModel, package = "goldfish.latent")
  )
  
  if (requireNamespace("cmdstanr", quietly = TRUE) &&
      cmdstanr::cmdstan_version() >= "2.29.2") {
    model <- cmdstanr::write_stan_file(code = stanCode)
  } else {
    stop(
      dQuote("cmdstanr"), " package and a working version of",
      dQuote("CmdStan"), "are required.",
      "\nPlease follow Stan documentation for instructions on how to install."
    )
  }
  
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



#' Sample preprocessed data
#'
#' For arguments `method_set` and `method_events` the available methods are
#' `"systematic"` and `"srswor"`. They correspond to
#'   systematic sampling and simple random sampling without replacement,
#'   respectively. It is possible to use an external function that has arguments
#'   `N` and `fraction` and return a numerical vector of the samples
#'   to keep.
#'
#' @param data a `goldfish.latent.data` object.
#' @param fraction_set a numerical value between 0 and 1.
#'   It is the proportion of cases to sample. If `sample_size_set` is not
#'   `NULL`, it is ignored.
#' @param sample_size_set an integerish value greater than 0.
#'   It defines the number of cases to sample. Ignores `fraction_set` if
#'   `sample_size_set` is not `NULL`.
#' @param method_set a character string indicating the sampling method for the
#'   choice set.
#' @param fraction_events a numerical value between 0 and 1.
#'   It is the proportion of events to sample. If `sample_size_events` is not
#'   `NULL`, it is ignored.
#' @param sample_size_events an integerish value greater than 0.
#'   It defines the number of events to sample. Ignores `fraction_events` if
#'   `sample_size_events` is not `NULL`.
#' @param method_events a character string indicating the sampling method for
#'   the events.
#'
#' @return a `goldfish.latent.data` object with the sampled data and the
#'   selection probabilities of the events/cases sampled.
#'
#' @export
#'
#' @examples
#' sampledData <- sample_data(data)
sample_data <- function(
    data,
    fraction_set = 1,
    sample_size_set = NULL,
    method_set = c("srswor", "systematic"),
    fraction_events = 1,
    sample_size_events = NULL,
    method_events = c("systematic", "srswor")) {
  method_set <- match.arg(method_set, c("srswor", "systematic"))
  method_events <- match.arg(method_events, c("systematic", "srswor"))

  stopifnot(
    inherits(data, "goldfish.latent.data"),
    rlang::is_scalar_double(fraction_set) &&
      fraction_set > 0 && fraction_set <= 1,
    is.null(sample_size_set) ||
      (rlang::is_scalar_integerish(sample_size_set) && sample_size_set > 0),
    is.character(method_set) && length(method_set) == 1,
    rlang::is_scalar_double(fraction_events) &&
      fraction_events > 0 && fraction_events <= 1,
    is.null(sample_size_events) ||
      (rlang::is_scalar_integerish(sample_size_events) &&
         sample_size_events > 0),
    is.character(method_events) && length(method_events) == 1
  )

  dataStan <- data$data_stan
  model <- attr(data, "model")
  sub_model <- attr(data, "sub_model")

  # Event sampling
  sample_events_rate <- NULL
  should_sample_events <- fraction_events < 1 || !is.null(sample_size_events)
  if (should_sample_events && sub_model %in% c("both", "rate")) {
    sample_events_rate <- do.call(
      method_events,
      list(
        N = dataStan$T_rate,
        fraction = fraction_events,
        sample_size = sample_size_events
      )
    )
  }

  sample_events_choice <- NULL
  if (should_sample_events && sub_model %in% c("both", "choice")) {
    if (sub_model == "both" && !is.null(sample_events_rate)) {
      if (dataStan$T_rate != dataStan$T_choice) {
        events_choice <- seq_len(dataStan$T_choice)
        events_rate <- rep(0, dataStan$T_rate)
        events_rate[dataStan$is_dependent] <- events_choice
      } else {
        events_rate <- seq_len(dataStan$T_rate)
      }
      events_choice <- events_rate[sample_events_rate$sample]
      events_choice <- events_choice[events_choice > 0]
      sample_events_choice <- list(
        sample = events_choice,
        sel_prob = sample_events_rate$sel_prob
      )
    } else {
      sample_events_choice <- do.call(
        method_events,
        list(
          N = dataStan$T_choice,
          fraction = fraction_events,
          sample_size = sample_size_events
        )
      )
    }
  }

  if (sub_model %in% c("both", "rate")) {
    data_rate <- sample_set(
      data_stan = dataStan, model = model, sub_model = "rate",
      method_set = method_set,
      fraction_set = fraction_set, sample_size_set = sample_size_set,
      sample_events = sample_events_rate
    )
  } else {
    data_rate <- NULL
  }

  if (sub_model %in% c("both", "choice")) {
    data_choice <- sample_set(
      data_stan = dataStan, model = model, sub_model = "choice",
      method_set = method_set,
      fraction_set = fraction_set, sample_size_set = sample_size_set,
      sample_events = sample_events_choice
    )
  } else {
    data_choice <- NULL
  }

  if (sub_model == "both" && model %in% c("DN_RE", "DN_HMM_RE", "DN_CHMM_RE")) {
    data_choice[["A"]] <- data_choice[["sender"]] <- NULL
  }
  data$data_stan <- c(data_rate, data_choice)

  if (model %in% c("DN_RE", "DN_HMM_RE", "DN_CHMM_RE")) {
    if (sub_model %in% c("both", "choice")) {
      n_actors <- length(unique(data_choice[["sender"]]))
      if (n_actors != data_choice[["A"]]) {
        cli_abort(c(
          "Number of actors in choice data does not match number of actors in rate data.",
          "x" = "Number of actors in choice data is {n_actors}.",
          "y" = "Number of actors in rate data is {data_choice$A}."
        ))
      }
    } else if (sub_model %in% c("rate")) {
      n_actors <- length(unique(data_rate["sender"]))
      if (n_actors != data_rate[["A"]]) {
        cli_abort(c(
          "Number of actors in rate data does not match number of actors in choice data.",
          "x" = "Number of actors in rate data is {n_actors}.",
          "y" = "Number of actors in choice data is {data_rate$A}."
        ))
      }
    }
  }
  attr(data, "sample") <- TRUE
  return(data)
}

add_columns <- function(array_to, array_from) {
  miss_cols <- setdiff(
    colnames(array_from), colnames(array_to)
  )
  array_to <- as.data.frame(array_to)
  if (length(miss_cols) > 0) {
    for (col_add in miss_cols) {
      array_to[, col_add] <- ""
    }
  }

  array_to
}

systematic <- function(N, fraction, sample_size = NULL) {
  stopifnot(
    rlang::is_scalar_integerish(N) && N >= 1,
    rlang::is_scalar_double(fraction) && fraction > 0 && fraction <= 1,
    is.null(sample_size) || sample_size > 0
  )

  if (N == 1) {
    return(list(sample = 1, sel_prob = 1))
  }
  if (!is.null(sample_size)) {
    is_integer <- rlang::is_scalar_integerish(sample_size, finite = TRUE)
    if (!is_integer) {
      cli_abort(c(
        "Sample size must be an integer if greater than 1.",
        "x" = "Sample size is {sample_size}."
      ))
    } else {
      sample_size <- ceiling(N * fraction)
    }
  }

  if (sample_size >= N) {
    return(list(sample = seq_len(N), sel_prob = 1))
  }

  jump <- floor(N / sample_size)
  r <- sample(jump, 1)
  
  list(
    sample = r + (jump * (seq_len(sample_size) - 1)),
    sel_prob = sample_size / N
  )
}

srswor <- function(N, fraction, sample_size = NULL) {
  stopifnot(
    rlang::is_scalar_integerish(N) && N >= 1,
    rlang::is_scalar_double(fraction) && fraction > 0 && fraction <= 1,
    is.null(sample_size) || sample_size > 0
  )

  if (N == 1) {
    return(list(sample = 1, sel_prob = 1))
  }
  if (!is.null(sample_size)) {
    is_integer <- rlang::is_scalar_integerish(sample_size, finite = TRUE)
    if (!is_integer) {
      cli_abort(c(
        "Sample size must be an integer if greater than 1.",
        "x" = "Sample size is {sample_size}."
      ))
    }
  } else {
    sample_size <- ceiling(N * fraction)
  }

  if (sample_size >= N) {
    return(list(sample = seq_len(N), sel_prob = 1))
  }

  list(
    sample = sort(sample.int(N, sample_size)),
    sel_prob = sample_size / N
  )
}

rescale_coefs <- function(
    beta, scale_stats, offset = 0, is_rate = TRUE) {
  beta2 <- beta ## inherit names etc.

  if (is_rate) {
    mu <- scale_stats$rate$`scaled:center`
    sigma <- scale_stats$rate$`scaled:scale`

    beta2[-1] <- beta[-1] / sigma
    beta2[1] <- beta[1] + offset - sum(beta2[-1] * mu)
  } else {
    sigma <- scale_stats$choice$`scaled:scale`
    beta2 <- beta / sigma
  }
  return(beta2)
}

change_k_regimes <- function(chr_vec, k_regimes) {
  which_line <- grepl("%%", chr_vec)
  chr_vec[which_line] <-
    gsub("%%kR%%", as.character(k_regimes), chr_vec[which_line])
  chr_vec
}
