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

#' Create data for Stan Hidden Markov Model model
#'
#' The first step is create the data in the structure expected by the `Stan`
#' code designed for the HMM-DyNAM.
#' Additional information of the effects used during preprocessing is added to
#' the return object.
#'
#' @param rate_effects a `formula` specification as in [goldfish::estimate()].
#' @param choice_effects a `formula` specification as in [goldfish::estimate()].
#' @param model Current version only support `"DyNAM"` model.
#' @param scale logical value. Whether to standardize the effect stats at
#' the end.
#' @param data a `data.goldfish` object. Output of [goldfish::make_data()].
#' @inheritParams make_data_re
#'
#' @return a list with the following components.
#' \describe{
#'   \item{dataStan}{a list with the information necessary to run a HMC using
#'   Stan.}
#'   \item{namesEffects}{a character vector with terms in the random and fixed
#'   effects formulas and their final name.}
#'   \item{effectDescription}{an array with detailed and comprehensible
#'   information of the terms used in the random and fixed effects formulas.}
#'   \item{scale}{a list with the statistics used for the standardization when
#'   `scale = TRUE`.}
#' }
#' @importFrom goldfish set_preprocessing_opt
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' data("Social_Evolution")
#' callNetwork <- make_network(nodes = actors, directed = TRUE) |>
#'   link_events(change_event = calls, nodes = actors)
#' callsDependent <- define_dependent_events(
#'   events = calls, nodes = actors, default_network = callNetwork
#' )
#' socialEvolutionData <- make_data(callsDependent)
#' data2stan <- make_data_hmm(
#'   rate_effects = callsDependent ~ indeg + outdeg,
#'   choice_effects = callsDependent ~ recip + trans,
#'   data = socialEvolutionData
#' )
#' }
make_data_hmm <- function(
    rate_effects,
    choice_effects,
    model = c("DyNAM"),
    k_states = 3,
    support_constraint = NULL,
    control_preprocessing = NULL,
    scale = TRUE,
    progress = getOption("progress"),
    data = NULL) {
  ### 0. check parameters----
  model <- match.arg(model)

  stopifnot(
    is.null(progress) || inherits(progress, "logical"),
    is.null(control_preprocessing) ||
      inherits(control_preprocessing, "preprocessing_opt.goldfish"),
    is.null(rate_effects) || inherits(rate_effects, "formula"),
    is.null(choice_effects) || inherits(choice_effects, "formula"),
    is.null(support_constraint) ||
      inherits(support_constraint, "formula"),
    is.numeric(k_states) && length(k_states) == 1 && k_states >= 2
  )

  # setting initial values of some arguments
  if (is.null(progress)) progress <- FALSE

  if (is.null(control_preprocessing)) {
    control_preprocessing <- goldfish::set_preprocessing_opt()
  }

  if (!is.null(rate_effects) && !is.null(choice_effects)) {
    sub_type <- "both"


    if (as.character(choice_effects[[2]]) != as.character(rate_effects[[2]])) {
      stop("dependent event network needs to be the same for both formulas")
    }
  } else if (!is.null(rate_effects)) {
    sub_type <- "rate"
  } else if (!is.null(choice_effects)) {
    sub_type <- "choice"
  } else {
    stop(
      dQuote("rateEffects"), " and ", dQuote("choiceEffects"),
      " arguments are NULL objects. Specify at least one model."
    )
  }
  # process data
  scale_stats <- if (scale) list() else NULL
  if (!is.null(rate_effects)) {
    data_processed_rate <- gather_model_data(
      formula = rate_effects,
      model = model,
      sub_model = "rate",
      control_preprocessing = control_preprocessing,
      progress = progress,
      data = data
    )

    # nEvents <- length(dataProcessedRate$selected)
    # expandedDF <- cbind(
    #   as.data.frame(dataProcessedRate$stat_all_events),
    #   data.frame(
    #     event = rep(seq.int(nEvents), dataProcessedRate$n_candidates),
    #     selected = sequence(dataProcessedRate$n_candidates) ==
    #       rep(dataProcessedRate$selected, dataProcessedRate$n_candidates)
    #   )
    # )
    # idxEvents <- tapply(seq.int(nrow(expandedDF)), expandedDF$event, range) |>
    #   simplify2array()


    data_stan_rate <- within(data_processed_rate, {
      Trate <- length(n_candidates)
      Nrate <- nrow(stat_all_events)
      Prate <- ncol(stat_all_events)
      startRate <- cumsum(c(1, head(n_candidates, -1)))
      endRate <- cumsum(n_candidates)
      namesEffects <- gsub("\\$", "Of", namesEffects)
      Xrate <- setNames(stat_all_events, namesEffects)
      if (hasIntercept) {
        choseRate <- selected[, 1] + (startRate - 1) * isDependent
        # mean(timespan)
        offsetInt <- log(Trate / (sum(timespan) * mean(n_candidates))) 
      } else {
        choseRate <- selected[, 1] + (startRate - 1)
      }
      rm(
        stat_all_events, n_candidates, sender,
        namesEffects, effectDescription, selected
      )
    })

    if (scale) {
      Xrate <- scale(
        if (dataStanRate$hasIntercept) {
          dataStanRate$Xrate[, -1]
        } else {
          dataStanRate$Xrate
        }
      )
      scale_stats[["rate"]] <- attributes(Xrate)

      if (data_stan_rate$hasIntercept) {
        data_stan_rate$Xrate[, -1] <- Xrate[, ]
      } else {
        data_stan_rate$Xrate <- Xrate[, ]
      }
    }
    # table(idxEvents[1, ] == dataStanRate$startRate)
    # table(idxEvents[2, ] == dataStanRate$endRate)
    # chose <- which(expandedDF$selected)
    # table(chose == dataStanRate$choseRate[dataStanRate$isDependent])
  } else {
    data_processed_rate <- data_stan_rate <- NULL
  }

  if (!is.null(choice_effects)) {
    feTerms <- terms(choice_effects)
    cstrTerms <- if (!is.null(support_constraint)) {
      terms(support_constraint)
    } else {
      NULL
    }

    if (length(attr(cstrTerms, "term.labels")) > 1) {
      stop(dQuote("support_constraint"), " argument only works for one effect.")
    }

    if (any(attr(feTerms, "order") != 1)) {
      stop(
        dQuote("choice_effects"),
        "formula argument doesn't support interactions yet"
      )
    }

    termsDyNAM <- c(
      attr(feTerms, "term.labels"),
      attr(cstrTerms, "term.labels")
    ) |> unique()

    formulaDyNAM <- reformulate(
      termsDyNAM,
      response = as.character(choice_effects[[2]])
    )

    data_processed_choice <- gather_model_data(
      formula = formulaDyNAM,
      model = model,
      sub_model = "choice",
      preprocess_args = preprocess_args,
      progress = progress,
      data = data
    )

    nEvents <- length(data_processed_choice$sender)

    namesEffects <- setNames(
      gsub("\\$", "Of", data_processed_choice$namesEffects),
      termsDyNAM
    )

    expandedDF <- cbind(
      setNames(
        as.data.frame(data_processed_choice$stat_all_events),
        namesEffects
      ),
      data.frame(
        event = rep(seq_len(nEvents), data_processed_choice$n_candidates),
        selected = sequence(data_processed_choice$n_candidates) ==
          rep(data_processed_choice$selected, data_processed_choice$n_candidates)
      )
    )

    # subset if constraint
    if (!is.null(support_constraint)) {
      cstrName <- namesEffects[attr(cstrTerms, "term.labels")]
      keep <- expandedDF[, cstrName] == 1
      expandedDF <- expandedDF[keep, !names(expandedDF) %in% cstrName]
      effectDescription <-
        data_processed_choice$effectDescription[!namesEffects %in% cstrName, ]
      namesEffects <- namesEffects[!namesEffects %in% cstrName]
    } else {
      effectDescription <- data_processed_choice$effectDescription
    }

    # create objects for Stan
    nTotal <- nrow(expandedDF)

    idxEvents <- tapply(seq_len(nTotal), expandedDF$event, range) |>
      simplify2array()

    Xmat <- as.matrix(expandedDF[, namesEffects])

    data_stan_choice <- list(
      Tchoice = nEvents,
      Nchoice = nTotal,
      Pchoice = ncol(Xmat),
      startChoice = idxEvents[1, ],
      endChoice = idxEvents[2, ],
      Xchoice = Xmat,
      choseChoice = which(expandedDF[, "selected"])
    )

    if (scale) {
      Xchoice <- scale(Xmat)
      scale_stats[["choice"]] <- attributes(Xchoice)

      data_stan_choice$Xchoice <- Xchoice[, ]
    }
  } else {
    namesEffects <- effectDescription <- NULL
  }

  data_stan <- c(
    if (!is.null(rate_effects)) data_stan_rate,
    if (!is.null(choice_effects)) data_stan_choice,
    list(
      kS = k_states,
      alpha = c(floor((k_states - 1) / (1 - 0.8)), 1)
    )
  )

  if (!is.null(effectDescription) &
    !is.null(data_processed_rate$effectDescription)) {
    if (!setequal(
      colnames(effectDescription),
      colnames(data_processed_rate$effectDescription)
    )) {
      data_processed_rate$effectDescription <- AddColumns(
        data_processed_rate$effectDescription,
        effectDescription
      )

      effectDescription <- AddColumns(
        effectDescription,
        data_processed_rate$effectDescription
      )
    }
  }

  return(structure(
    list(
      data_stan = data_stan,
      namesEffects = c(
        data_processed_rate$namesEffects,
        namesEffects
      ),
      effectDescription = rbind(
        data_processed_rate$effectDescription,
        effectDescription
      ),
      scale_stats = scale_stats
    ),
    class = "goldfish.latent.data",
    model = "DNHMM",
    subModel = sub_type
  ))
}

#' Post process MCMC samples from a HMM
#'
#' The post-processing allows to compute for each posterior draw:
#' the most probable path of the Hidden States using the Viterbi algorithm,
#' the smoothed probabilities of a Hidden State given the observed sequence
#' using the Forward-Backward algorithm, and a sample path using the smoothed
#' probabilities `smoothProbsSt = "marginal"` or taking into account the
#' conditional probability given the next state  `smoothProbsSt = "joint"` as
#' implemented in Stan `hmm_latent_rng()`.
#' @param data_to_stan a object of class `"goldfish.latent.data"` outcome of
#' a call to [make_data_hmm()].
#' @param cmdstan_samples a object of class `"CmdStanMCMC"` with the posterior
#' samples of a HMM from the `data_to_stan`.
#' @param k_states an integer with the number of states used to get samples
#' from the posterior distribution in `cmdstan_samples`.
#' @param type a character specifying which algorithm to use for the post
#' processing of the posterior distribution samples.
#' @param smooth_probs_states a character specifying the way that state sample path
#' are drawn either using the smoothed probabilities `"marginal"` or
#' conditioning of the next sample state `"joint"`.
#' Only used when `type = "both"` or `type = "smoothProbs"`.
#' @param cl an optional `cluster` object to parallelize the
#' post-processing
#' @param n_nodes an integer specifying the number of nodes to use in the
#' parallelization.
#' @return a list with the most probable path when Viterbi algorithm is used
#' and/or the smoothed probabilities and sample paths when the
#' Forward-Backward algorithm is used.
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' data("Social_Evolution")
#' callNetwork <- make_network(nodes = actors, directed = TRUE) |>
#'   link_events(change_event = calls, nodes = actors)
#' callsDependent <- define_dependent_events(
#'   events = calls, nodes = actors, default_network = call_network
#' )
#' socialEvolutionData <- make_data(callsDependent)
#' data_to_stan <- make_data_hmm(
#'   rate_effects = callsDependent ~ indeg + outdeg,
#'   choice_effects = NULL,
#'   data = socialEvolutionData
#' )
#'
#' postProcess <- hmm_post_processing(data_to_stan, cmdstan_samples)
#' }
hmm_post_processing <- function(
    data_to_stan, cmdstan_samples,
    k_states = data_to_stan$data_stan$kS,
    type = c("both", "viterbi", "smoothProbs"),
    smooth_probs_states = c("joint", "marginal", "none"),
    cl = NULL, n_nodes = parallel::detectCores() - 1L) {
  type <- match.arg(type)
  smooth_probs_states <- match.arg(smooth_probs_states)

  model <- attr(data_to_stan, "model")
  stopifnot(
    inherits(cmdstan_samples, "CmdStanMCMC"),
    inherits(data_to_stan, "goldfish.latent.data"),
    is.numeric(k_states) && length(k_states) == 1 && k_states >= 2,
    model %in% c("DNHMM", "DNCHMM", "DNCHMMRE")
  )

  sub_model <- attr(data_to_stan, "subModel")

  # extract draws and reformat for posterior computations
  typeOutput <- ifelse(is.null(cl), "draws_matrix", "draws_array")
  typeOutput <- "draws_matrix"  # change parallelization

  drawsObject <- hmm_draws_to_label_switching(
    data_to_stan, cmdstan_samples,
    k_states = k_states, type = typeOutput
  )

  data_to_stan <- data_to_stan$data_stan

  if (is.null(cl)) {
    output <- do.call(
      switch(
        model,
         "DNHMM" = "hmm_pp_per_draws",
         "DNCHMM" = "chmm_pp_per_draws",
         "DNCHMMRE" = "chmm_re_pp_per_draws"
      ),
      args = list(
        chain_iter = NULL,
        draws_object = drawsObject,
        data_stan = data_to_stan,
        model = model,
        sub_model = sub_model,
        k_states = k_states,
        type = type,
        smooth_probs_state = smooth_probs_state
      )
    )
  } else {
    on.exit(parallel::stopCluster(cl))
    ignore <- parallel::clusterEvalQ(cl, {
      library(matrixStats)
      library(expm)
      NULL
    })

    nDraws <- nrow(drawsObject$draws)
    kFactor <- nDraws %% nnodes %% 500

    output <- parallel::clusterApplyLB(
      cl = cl,
      seq_len(nnodes * kFactor),
      fun = switch(
        model,
         "DNHMM" = hmm_pp_per_draws,
         "DNCHMM" = chmm_pp_per_draws,
         "DNCHMMRE" = chmm_re_pp_per_draws
      ),
      draws_object = drawsObject,
      data_stan = data_to_stan,
      model = model,
      sub_model = sub_model,
      k_states = k_states,
      type = type,
      smooth_probs_st = smooth_probs_st,
      nnodes = nnodes,
      k_factor = k_factor
    ) |>
      bind_pp_hmm(type = type, smooth_probs_state = smooth_probs_state)
  }

  return(output)
}

hmm_pp_per_draws <- function(
    chain_iter, draws_object, data_stan, model, sub_model,
    k_states, type, smooth_probs_states, n_nodes, k_factor) {
  # init output
  output <- list()

  draws <- draws_object$draws
  idx_theta <- draws_object$idxTheta
  idx_emission <- draws_object$idxEmission
  idx_beta_choice <- draws_object$idxBetaChoice
  idx_beta_rate <- draws_object$idxBetaRate

  if (!is.null(chain_iter)) {
    split_draws <- parallel::splitIndices(nrow(draws), n_nodes * k_factor)
    draws <- draws[split_draws[[chain_iter]], ]
  }


  n_draws <- nrow(draws)

  # define sizes
  n_events <- ifelse(sub_model %in% c("both", "rate"), "Trate", "Tchoice")
  n_events <- data_stan[[n_events]]

  is_res <- !is.null(data_stan$Nres)
  TT <- ifelse(is_res, data_stan$Nres, n_events)


  if (sub_model %in% c("choice")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (k_s in seq_len(k_states)) {
      xb <- tcrossprod(data_stan$Xchoice, draws[, idx_beta_choice[k_s, ]])

      for (event in seq_len(n_events)) {
        ll_event_state[event, , k_s] <-
          xb[data_stan$choseChoice[event], ] -
          colLogSumExps(
            xb,
            rows = seq(
              data_stan$startChoice[event],
              data_stan$endChoice[event]
            )
          )
      }
    }
  }

  if (sub_model %in% c("rate")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (k_s in seq_len(k_states)) {
      xb <- tcrossprod(data_stan$Xrate, draws[, idx_beta_rate[k_s, ]]) +
        data_stan$offsetInt

      for (event in seq_len(n_events)) {
        ll_event_state[event, , k_s] <-
          ifelse(data_stan$isDependent[event],
            xb[data_stan$choseRate[event], ], 0
          ) -
          data_stan$timespan[event] * exp(colLogSumExps(
            xb,
            rows = seq(
              data_stan$startRate[event],
              data_stan$endRate[event]
            )
          ))
      }
    }
  }

  if (sub_model %in% c("both")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (k_s in seq_len(k_states)) {
      xb_r <- tcrossprod(data_stan$Xrate, draws[, idx_beta_rate[k_s, ]]) +
        data_stan$offsetInt
      xb_c <- tcrossprod(data_stan$Xchoice, draws[, idx_beta_choice[k_s, ]])

      event_choice <- 1L
      for (event in seq_len(n_events)) {
        loglik <- -data_stan$timespan[event] * exp(colLogSumExps(
          xb_r,
          rows = seq(
            data_stan$startRate[event],
            data_stan$endRate[event]
          )
        ))
        if (data_stan$isDependent[event]) {
          loglik <- loglik + xb_r[data_stan$choseRate[event], ] +
            xb_c[data_stan$choseChoice[event_choice], ] -
            colLogSumExps(
              xb_c,
              rows = seq(
                data_stan$startChoice[event_choice],
                data_stan$endChoice[event_choice]
              )
            )
          event_choice <- event_choice + 1L
        }

        ll_event_state[event, , k_s] <- loglik
      }
    }
  }

  if (is_res) {
    ll_event_state <- array(apply(
      ll_event_state,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(data_stan$resA[y], data_stan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, n_draws, k_states))
  }

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, k_states, n_draws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, k_states, n_draws))

    # forward past computing most likely state from previous state
    # first observation: p(y_1| z_1) * p(z_1) (emission prob)
    delta[1, , ] <- t(log(draws[, idx_emission]) + ll_event_state[1, , ])

    for (tt in seq(2, TT)) {
      prev_delta <- t(delta[tt - 1, , ])
      for (k in seq_len(k_states)) {
        t_1_xj <- prev_delta + log(draws[, idx_theta[, k]]) +
          ll_event_state[tt, , k]

        delta[tt, k, ] <- apply(t_1_xj, 1, max)
        bpointer[tt, k, ] <- apply(t_1_xj, 1, which.max)
      }
    }

    # backward past
    z <- array(0L,
      dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT))
    )

    z[, TT] <- apply(delta[TT, , ], 2, which.max)

    for (tt in seq(TT - 1, 1)) {
      z[, tt] <- sapply(
        seq_len(n_draws),
        \(x) bpointer[tt + 1, z[x, tt + 1], x]
      )
    }

    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(n_draws, TT, k_states),
      dimnames = list(
        draws = seq_len(n_draws), time = seq_len(TT),
        states = seq_len(k_states)
      )
    )

    if (smooth_probs_states != "none") {
      z_sample <- array(0L, dim = c(n_draws, TT))
    }

    # forward -- filtering: alpha [tt, k] with running normalization
    #  11.2.2, Finite Mixture and Markov Switching models,
    #  Frühwirth, S., 2006

    # # filter at first observation:
    # # p(S_1 = k|y_0) = p(y_1| S_1) * p(S_1); (emission prob)
    p[, 1, ] <- log(draws[, idx_emission]) + ll_event_state[1, , ]
    # In Frühwirth: normalization is over llEventState, here follow Stan
    p[, 1, ] <- sweep(p[, 1, ], 1, apply(p[, 1, ], 1, max))
    # not need to convert to prob
    # p[, 1, ] <- sweep(p[, 1, ], 1, rowLogSumExps(p[, 1, ]))

    for (tt in seq(2, TT)) {
      # # one-step ahead prediction of S_t:
      # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
      for (k in seq_len(k_states)) {
        p[, tt, k] <- rowLogSumExps(
          log(draws[, idx_theta[, k]]) + p[, tt - 1, ]
        ) +
          # # filter for S_t: p(S_t = k| y_t) =
          # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
          # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) ; unnormalize enough
          ll_event_state[tt, , k]
      }
    }

    # Normalize last value, already smooth distribution
    p[, TT, ] <- exp(sweep(p[, TT, ], 1, rowLogSumExps(p[, TT, ])))

    if (smooth_probs_states != "none") { # sample last Hidden State (HS)
      z_sample[, TT] <- apply(
        p[, TT, ], 1,
        \(x) sample.int(k_states, 1, prob = x)
      )
    }

    # backward: smoother suggested in Hamilton expresses these as marginal
    # probabilities from the joint distribution of S_t and S_T | y
    # Implementation follows Stan hmm_hidden_state_prob()
    #

    # initial ending state ass as given (uniform)
    log_beta <- array(0, dim = c(n_draws, k_states))

    for (tt in seq(TT - 1, 1)) {
      # # Baum-Welch alg
      # #
      omega_beta <- log_beta + ll_event_state[tt + 1, , ] # element-wise product

      # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
      if (smooth_probs_states == "joint") {
        prob_last_hs <- p[, tt, ] +
          t(vapply(
            seq_len(n_draws),
            \(x){
              last_hs <- z_sample[x, tt + 1]
              log(draws[x, idx_theta[, last_hs]]) + omega_beta[x, last_hs]
            },
            numeric(k_states)
          ))
        prob_last_hs <- exp(sweep(prob_last_hs, 1, rowLogSumExps(prob_last_hs)))
        z_sample[, tt] <- apply(
          prob_last_hs, 1,
          \(x) sample.int(k_states, 1, prob = x)
        )
      }
      for (k in seq_len(k_states)) {
        log_beta[, k] <- rowLogSumExps(
          log(draws[, idx_theta[k, ]]) + omega_beta
        )
      }

      # running normalization
      log_beta <- sweep(log_beta, 1, apply(log_beta, 1, max))

      #
      gamma_t <- log_beta + p[, tt, ]
      p[, tt, ] <- exp(sweep(gamma_t, 1, rowLogSumExps(gamma_t)))

      # sample tt HS from marginal distribution
      if (smooth_probs_states == "marginal") {
        z_sample[, tt] <- apply(
          p[, tt, ], 1,
          \(x) sample.int(k_states, 1, prob = x)
        )
      }
    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smooth_probs_states != "none") {
      output[["smoothProbs"]][["zSample"]] <- z_sample
    }
  }

  return(output)
}

#' @importFrom expm expm
chmm_pp_per_draws <- function(
    chain_iter, draws_object, data_stan, model, sub_model,
    k_states, type, smooth_probs_states, n_nodes, k_factor) {
  # init output
  output <- list()

  draws <- draws_object$draws
  idx_theta <- draws_object$idx_theta
  idx_emission <- draws_object$idx_emission
  idx_beta_choice <- draws_object$idx_beta_choice
  idx_beta_rate <- draws_object$idx_beta_rate

  if (!is.null(chain_iter)) {
    split_draws <- parallel::splitIndices(nrow(draws), n_nodes * k_factor)
    draws <- draws[split_draws[[chain_iter]], ]
  }


  n_draws <- nrow(draws)

  # define sizes
  n_events <- ifelse(sub_model %in% c("both", "rate"), "Trate", "Tchoice")
  n_events <- data_stan[[n_events]]

  is_res <- !is.null(data_stan$Nres)
  TT <- ifelse(is_res, data_stan$Nres, n_events)

  log_crude_rate <- ifelse(
    sub_model %in% c("both", "rate"),
    log(data_stan$Trate / (data_stan$Nrate * mean(data_stan$timespan))), 0
  )

  if (sub_model %in% c("choice")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      xb <- tcrossprod(data_stan$Xchoice, draws[, idx_beta_choice[kS, ]])

      for (event in seq_len(n_events)) {
        ll_event_state[event, , kS] <-
          xb[data_stan$choseChoice[event], ] -
          colLogSumExps(
            xb,
            rows = seq(
              data_stan$startChoice[event],
              data_stan$endChoice[event]
            )
          )
      }
    }
  }

  if (sub_model %in% c("rate")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      xb <- tcrossprod(data_stan$Xrate, draws[, idx_beta_rate[kS, ]]) +
        log_crude_rate

      for (event in seq_len(n_events)) {
        ll_event_state[event, , kS] <-
          ifelse(data_stan$isDependent[event],
            xb[data_stan$choseRate[event], ], 0
          ) -
          data_stan$timespan[event] * exp(colLogSumExps(
            xb,
            rows = seq(
              data_stan$startRate[event],
              data_stan$endRate[event]
            )
          ))
      }
    }
  }

  if (sub_model %in% c("both")) {
    ll_event_state <- array(0, dim = c(n_events, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      xbR <- tcrossprod(data_stan$Xrate, draws[, idx_beta_rate[kS, ]]) +
        log_crude_rate
      xbC <- tcrossprod(data_stan$Xchoice, draws[, idx_beta_choice[kS, ]])

      event_choice <- 1L
      for (event in seq_len(n_events)) {
        loglik <- -data_stan$timespan[event] * exp(colLogSumExps(
          xbR,
          rows = seq(
            data_stan$startRate[event],
            data_stan$endRate[event]
          )
        ))
        if (data_stan$isDependent[event]) {
          loglik <- loglik + xbR[data_stan$choseRate[event], ] +
            xbC[data_stan$choseChoice[event_choice], ] -
            colLogSumExps(
              xbC,
              rows = seq(
                data_stan$startChoice[event_choice],
                data_stan$endChoice[event_choice]
              )
            )
          event_choice <- event_choice + 1L
        }

        ll_event_state[event, , kS] <- loglik
      }
    }
  }

  if (is_res) {
    ll_event_state <- array(apply(
      ll_event_state,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(data_stan$resA[y], data_stan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, n_draws, k_states))
  }

  # compute transition probabilities
  trans_probs <- vapply(
    data_stan$timespan,
    \(y)
    apply(
      draws[, idx_theta] * y,
      1,
      \(x){
        log(expm::expm(
          matrix(
            x,
            nrow = k_states, ncol = k_states
          )
        ))
      }
    ),
    matrix(0, nrow = k_states^2, ncol = n_draws)
  )

  idxTP <- idx_theta - min(idx_theta) + 1

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, k_states, n_draws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, k_states, n_draws))

    # forward past computing most likely state from previous state
    # first observation: p(y_1, z_1) = p(y_1| z_1) * p(z_1) (emission prob)
    delta[1, , ] <- t(log(draws[, idx_emission]) + ll_event_state[1, , ])

    for (tt in seq(2, TT)) {
      prevDelta <- t(delta[tt - 1, , ])
      for (k in seq_len(k_states)) {
        T_1_xj <- prevDelta + t(trans_probs[idxTP[, k], , tt]) +
          ll_event_state[tt, , k]

        delta[tt, k, ] <- apply(T_1_xj, 1, max)
        bpointer[tt, k, ] <- apply(T_1_xj, 1, which.max)
      }
    }

    # backward past
    z <- array(0L,
      dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT))
    )

    z[, TT] <- apply(delta[TT, , ], 2, which.max)

    for (tt in seq(TT - 1, 1)) {
      z[, tt] <- sapply(
        seq_len(n_draws),
        \(x) bpointer[tt + 1, z[x, tt + 1], x]
      )
    }

    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(n_draws, TT, k_states),
      dimnames = list(
        draws = seq_len(n_draws), time = seq_len(TT),
        states = seq_len(k_states)
      )
    )

    # log likelihood as sum of scaling factors, Kadhem 2015
    # those are the marginal probabilities of y_t
    logLikRec <- array(0, dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT)))

    # log likelihood conditional on the sample state, Kadhem 2015 and 2021
    logLikCond <- array(0, dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT)))

    # sampling of the last state, required for conditional logLik
    if (smooth_probs_states != "none") {
      zSample <- array(0L, dim = c(n_draws, TT))
    }

    # forward -- filtering: alpha [tt, k] with running normalization
    #  11.2.2, Finite Mixture and Markov Switching models,
    #  Frühwirth, S., 2006

    # # filter at first observation:
    # # p(S_1 = k, y_1) = p(y_1| S_1) * p(S_1); (emission prob)
    p[, 1, ] <- log(draws[, idx_emission]) + ll_event_state[1, , ]
    # Compute normalization constants Kadhem 2015 (It's just the marginal prob)
    # p(y_1) = \sum_k p(y_1, S_1 = k), marginal prob up to t, Frühwirth 2006
    logLikRec[, 1] <- rowLogSumExps(p[, 1, ])
    # Filtering for s_t: p(S_t = k| y_t) = p(y_t, S_t = k) / p(y_t)
    p[, 1, ] <- sweep(p[, 1, ], 1, logLikRec[, 1])

    # objects intermediate values
    ptt <- array(0, dim = c(n_draws, k_states))
    for (tt in seq(2, TT)) {
      # # one-step ahead prediction of S_t:
      # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
      for (k in seq_len(k_states)) {
        ptt[, k] <- rowLogSumExps(
          t(transProbs[idxTP[, k], ,tt]) + p[, tt - 1, ]) +
          # # filter for S_t: p(S_t = k| y_t) =
          # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
          # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1})
          ll_event_state[tt, , k]
      }
      # marginal prob (y) and filter probs P(S_t = k | y_t)
      marginalProb <- rowLogSumExps(ptt)
      logLikRec[, tt] <- marginalProb

      p[, tt, ] <- sweep(ptt, 1, marginalProb)
    }

    # Normalize last value, already smooth distribution
    p[, TT, ] <- exp(p[, TT, ])

    if (smooth_probs_states != "none") { # sample last Hidden State (HS)
      zSample[, TT] <- apply(
        p[, TT, ], 1,
        \(x) sample.int(k_states, 1, prob = exp(x))
      )

      # Compute conditional log likelihood
      logLikCond[, TT] <- vapply(
        seq_len(n_draws),
        \(x) ll_event_state[TT, x, zSample[x, TT]],
        numeric(1)
      )
    }

    # backward-smoothing the States: suggested in Hamilton expresses
    # these as marginal probabilities from the joint distribution
    # of S_t and S_T | y
    # Implementation follows Stan hmm_hidden_state_prob()
    #

    # # Baum-Welch alg
    # initial ending state is \beta_i(T) = 1
    logBeta <- array(0, dim = c(n_draws, k_states))

    for (tt in seq(TT - 1, 1)) {
      # \beta_j(t) * p(y_{t+1} | S_{t+1} = j)
      omegaBeta <- logBeta + ll_event_state[tt + 1, , ] # element-wise product

      # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
      # FFBS, a.k.a, Multi-move sampling
      if (smooth_probs_states == "joint") {
        probLastHS <- p[, tt, ] +
          t(vapply(
            seq_len(n_draws),
            \(x){
              lastHS <- zSample[x, tt + 1]
              # p(S_{t+1} = j|S_t = i) * p(y_{t+1}|S_{t+1} = j) * \beta_j(t)
              t(transProbs[idxTP[, lastHS], x, tt]) + omegaBeta[x, lastHS]
            },
            numeric(k_states)
          ))
        probLastHS <- exp(sweep(probLastHS, 1, rowLogSumExps(probLastHS)))
        zSample[, tt] <- apply(
          probLastHS, 1,
          \(x) sample.int(k_states, 1, prob = x)
        )

        # Compute conditional log likelihood
        logLikCond[, tt] <- vapply(
          seq_len(n_draws),
          \(x) ll_event_state[tt, x, zSample[x, tt]],
          numeric(1)
        )
      }
      for (k in seq_len(k_states)) {
        # p(S_{t+1} = j | S_t = i) * p(y_{t+1} | S_{t+1} = j) * \beta_j(t)
        logBeta[, k] <- rowLogSumExps(
          t(transProbs[idxTP[k, ], , tt]) + omegaBeta
        )
      }

      # running normalization
      logBeta <- sweep(logBeta, 1, apply(logBeta, 1, max))

      # Update probs using bayes rule: P(X_t|y) = P(y | X_t) * P(X_t | y) / P(y)
      # \beta_i(t) * \alpha_i(t)
      gammat <- logBeta + p[, tt, ]
      p[, tt, ] <- exp(sweep(gammat, 1, rowLogSumExps(gammat)))

      # sample tt HS from marginal distribution
      if (smooth_probs_states == "marginal") {
        zSample[, tt] <- apply(
          p[, tt, ], 1,
          \(x) sample.int(k_states, 1, prob = x)
        )
      }
    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smooth_probs_states != "none") {
      output[["smoothProbs"]][["zSample"]] <- zSample
    }

    output[["logLikRec"]] <- logLikRec
    output[["logLikCond"]] <- logLikCond
  }

  return(output)
}

chmm_re_pp_per_draws <- function(
    chain_iter, draws_object, data_stan, model, sub_model,
    k_states, type, smooth_probs_st, n_nodes, k_factor) {
  # init output
  output <- list()

  draws <- draws_object$draws
  idxTheta <- draws_object$idxTheta
  idxEmission <- draws_object$idxEmission
  idxBetaChoice <- draws_object$idxBetaChoice
  idxBetaRate <- draws_object$idxBetaRate

  idxGammaChoice <- draws_object$idxGammaChoice
  idxGammaRate <- draws_object$idxGammaRate
  idxLOmegaChoice <- draws_object$idxLOmegaChoice
  idxLOmegaRate <- draws_object$idxLOmegaRate
  idxLSigmaChoice <- draws_object$idxLSigmaChoice


  if (!is.null(chain_iter)) {
    split_draws <- parallel::splitIndices(nrow(draws), n_nodes * k_factor)
    draws <- draws[split_draws[[chain_iter]], ]
  }


  n_draws <- nrow(draws)

  # define sizes
  nEvents <- ifelse(sub_model %in% c("both", "rate"), "Trate", "Tchoice")
  nEvents <- data_stan[[nEvents]]

  isRes <- !is.null(data_stan$Nres)
  TT <- ifelse(isRes, data_stan$Nres, nEvents)

  logCrudeRate <- ifelse(
    sub_model %in% c("both", "rate"),
    log(data_stan$Trate / (data_stan$Nrate * mean(data_stan$timespan))), 0
  )

  # indices group
  startGroup <- c(data_stan$startGroup, TT + 1)
  group <- data_stan$groupChoice

  if (sub_model %in% c("choice")) {
    Xchoice <- data_stan$Xchoice
    Zchoice <- data_stan$Xchoice[, data_stan$V1choice]

    ll_event_state <- array(0, dim = c(nEvents, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(
          data_stan$startChoice[event], data_stan$endChoice[event]
        )
        eventChose <- data_stan$choseChoice[event] -
          data_stan$startChoice[event] + 1
        xb <- tcrossprod(Xchoice[eventIdx, ], draws[, idxBetaChoice[kS, ]]) +
          tcrossprod(
            Zchoice[eventIdx, ],
            draws[, idxGammaChoice[kS, group[eventIdx[1]], ]]
          )

        ll_event_state[event, , kS] <- xb[eventChose, ] - colLogSumExps(xb)
      }
    }
  }

  if (sub_model %in% c("rate")) {
    Xrate <- data_stan$Xrate
    Zrate <- data_stan$Xrate[, data_stan$V1rate]

    ll_event_state <- array(0, dim = c(nEvents, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(data_stan$startRate[event], data_stan$endRate[event])
        eventChose <- data_stan$choseRate[event] -
          data_stan$startRate[event] + 1
        xb <- tcrossprod(Xrate[eventIdx, ], draws[, idxBetaRate[kS, ]]) +
          tcrossprod(Zrate[eventIdx, ],
            draws[, idxGammaRate[kS, group[eventIdx[1]], ]]) +
            logCrudeRate
        ll_event_state[event, , kS] <-
          ifelse(data_stan$isDependent[event], xb[eventChose, ], 0) -
          data_stan$timespan[event] * exp(colLogSumExps(xb))
      }
    }
  }

  if (sub_model %in% c("both")) {
    Xchoice <- data_stan$Xchoice
    Zchoice <- data_stan$Xchoice[, data_stan$V1choice]
    Xrate <- data_stan$Xrate
    Zrate <- data_stan$Xrate[, data_stan$V1rate]

    ll_event_state <- array(0, dim = c(nEvents, n_draws, k_states))

    for (kS in seq_len(k_states)) {
      xbR <- tcrossprod(Xrate, draws[, idxBetaRate[kS, ]]) +
        logCrudeRate
      xbC <- tcrossprod(Xchoice, draws[, idxBetaChoice[kS, ]])

      eventChoice <- 1L
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(data_stan$startRate[event], data_stan$endRate[event])
        eventChose <- data_stan$choseRate[event] -
          data_stan$startRate[event] + 1
        xbR <- tcrossprod(Xrate[eventIdx, ], draws[, idxBetaRate[kS, ]]) +
          tcrossprod(Zrate[eventIdx, ],
            draws[, idxGammaRate[kS, group[eventIdx[1]], ]]) +
            logCrudeRate
        loglik <- -data_stan$timespan[event] * exp(colLogSumExps(xbR))
        if (data_stan$isDependent[event]) {
          loglik <- loglik + xbR[eventChose, ]

          eventIdx <- seq(
            data_stan$startChoice[eventChoice],
            data_stan$endChoice[eventChoice]
          )
          eventChose <- data_stan$choseChoice[eventChoice] -
            data_stan$startChoice[eventChoice] + 1

          xbC <- tcrossprod(Xchoice[eventIdx, ], draws[, idxBetaChoice[kS, ]]) +
            tcrossprod(Zchoice[eventIdx, ],
              draws[, idxGammaChoice[kS, group[eventIdx[1]], ]])

          loglik <- loglik + xbC[eventChose, ] - colLogSumExps(xbC)
          eventChoice <- eventChoice + 1L
        }

        ll_event_state[event, , kS] <- loglik
      }
    }
  }

  if (isRes) {
    ll_event_state <- array(apply(
      ll_event_state,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(data_stan$resA[y], data_stan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, n_draws, k_states))
  }

  # compute transition probabilities
  transProbs <- vapply(
    data_stan$timesSender,
    \(y)
    apply(
      draws[, idxTheta] * y,
      1,
      \(x){
        log(expm::expm(
          matrix(
            x,
            nrow = k_states, ncol = k_states
          )
        ))
      }
    ),
    matrix(0, nrow = k_states^2, ncol = n_draws)
  )

  idxTP <- idxTheta - min(idxTheta) + 1

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, k_states, n_draws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, k_states, n_draws))

    # forward past
    z <- array(0L,
      dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT))
    )

    for (gg in seq_len(data_stan$G)) {
      n_events_g <- startGroup[gg + 1] - startGroup[gg]
      eventIdx <- seq(startGroup[gg], startGroup[gg + 1] - 1)

      # forward past computing most likely state from previous state
      # first observation: p(y_1, z_1) = p(y_1| z_1) * p(z_1) (emission prob)
      delta[eventIdx[1], , ] <- t(
        log(draws[, idxEmission]) + ll_event_state[eventIdx[1], , ]
      )

      for (tt in tail(eventIdx, -1)) {
        prevDelta <- t(delta[tt - 1, , ])
        for (k in seq_len(k_states)) {
          T_1_xj <- prevDelta + t(transProbs[idxTP[, k], , tt]) +
            ll_event_state[tt, , k]

          delta[tt, k, ] <- apply(T_1_xj, 1, max)
          bpointer[tt, k, ] <- apply(T_1_xj, 1, which.max)
        }
      }

      # backward past
      z[, tail(eventIdx, 1)] <-
        apply(delta[tail(eventIdx, 1), , ], 2, which.max)

      for (tt in rev(head(eventIdx, -1))) {
        z[, tt] <- sapply(
          seq_len(n_draws),
          \(x) bpointer[tt + 1, z[x, tt + 1], x]
        )
      }
    }
    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(n_draws, TT, k_states),
      dimnames = list(
        draws = seq_len(n_draws), time = seq_len(TT),
        states = seq_len(k_states)
      )
    )

    # log likelihood as sum of scaling factors, Kadhem 2015
    # those are the marginal probabilities of y_t
    logLikRec <- array(0, dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT)))

    # log likelihood conditional on the sample state, Kadhem 2015 and 2021
    logLikCond <- array(0, dim = c(n_draws, TT),
      dimnames = list(draws = seq_len(n_draws), time = seq_len(TT)))

    # sampling of the last state, required for conditional logLik
    if (smooth_probs_states != "none") {
      zSample <- array(0L, dim = c(n_draws, TT))
    }

    for (gg in seq_len(data_stan$G)) {
      n_events_g <- startGroup[gg + 1] - startGroup[gg]
      eventIdx <- seq(startGroup[gg], startGroup[gg + 1] - 1)

      # forward -- filtering: alpha [tt, k] with running normalization
      #  11.2.2, Finite Mixture and Markov Switching models,
      #  Frühwirth, S., 2006

      # # filter at first observation:
      # # p(S_1 = k, y_1) = p(y_1| S_1) * p(S_1); (emission prob)
      p[, eventIdx[1], ] <- log(draws[, idxEmission]) +
        ll_event_state[eventIdx[1], , ]
      # Compute normalization constants Kadhem 2015 (It's just the marginal prob)
      # p(y_1) = \sum_k p(y_1, S_1 = k), marginal prob up to t, Frühwirth 2006
      logLikRec[, eventIdx[1]] <- rowLogSumExps(p[, eventIdx[1], ])
      # Filtering for s_t: p(S_t = k| y_t) = p(y_t, S_t = k) / p(y_t)
      p[, eventIdx[1], ] <- sweep(
        p[, eventIdx[1], ], 1, logLikRec[, eventIdx[1]]
      )

      # objects intermediate values
      ptt <- array(0, dim = c(n_draws, k_states))
      for (tt in tail(eventIdx, -1)) {
        # # one-step ahead prediction of S_t:
        # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
        for (k in seq_len(k_states)) {
          ptt[, k] <- rowLogSumExps(
            t(transProbs[idxTP[, k], ,tt]) + p[, tt - 1, ]) +
            # # filter for S_t: p(S_t = k| y_t) =
            # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
            # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1})
            ll_event_state[tt, , k]
        }
        # marginal prob (y) and filter probs P(S_t = k | y_t)
        marginalProb <- rowLogSumExps(ptt)
        logLikRec[, tt] <- marginalProb

        p[, tt, ] <- sweep(ptt, 1, marginalProb)
      }

      # Normalize last value, already smooth distribution
      p[, tail(eventIdx, 1), ] <- exp(p[, tail(eventIdx, 1), ])

      if (smooth_probs_states != "none") { # sample last Hidden State (HS)
        zSample[, tail(eventIdx, 1)] <- apply(
          p[, tail(eventIdx, 1), ], 1,
          \(x) sample.int(k_states, 1, prob = exp(x))
        )

        # Compute conditional log likelihood
        tail_eventIdx <- tail(eventIdx, 1)
        logLikCond[, tail_eventIdx] <- vapply(
          seq_len(n_draws),
          \(x) ll_event_state[tail_eventIdx, x, zSample[x, tail_eventIdx]],
          numeric(1)
        )
      }

      # backward-smoothing the States: suggested in Hamilton expresses
      # these as marginal probabilities from the joint distribution
      # of S_t and S_T | y
      # Implementation follows Stan hmm_hidden_state_prob()
      #

      # # Baum-Welch alg
      # initial ending state is \beta_i(T) = 1
      logBeta <- array(0, dim = c(n_draws, k_states))

      for (tt in rev(head(eventIdx, -1))) {
        # \beta_j(t) * p(y_{t+1} | S_{t+1} = j)
        # element-wise product
        omegaBeta <- logBeta + ll_event_state[tt + 1, , ]

        # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
        # FFBS, a.k.a, Multi-move sampling
        if (smooth_probs_states == "joint") {
          probLastHS <- p[, tt, ] +
            t(vapply(
              seq_len(n_draws),
              \(x){
                lastHS <- zSample[x, tt + 1]
                # p(S_{t+1} = j|S_t = i) * p(y_{t+1}|S_{t+1} = j) * \beta_j(t)
                t(transProbs[idxTP[, lastHS], x, tt]) + omegaBeta[x, lastHS]
              },
              numeric(k_states)
            ))
          probLastHS <- exp(sweep(probLastHS, 1, rowLogSumExps(probLastHS)))
          zSample[, tt] <- apply(
            probLastHS, 1,
            \(x) sample.int(k_states, 1, prob = x)
          )

          # Compute conditional log likelihood
          logLikCond[, tt] <- vapply(
            seq_len(n_draws),
            \(x) ll_event_state[tt, x, zSample[x, tt]],
            numeric(1)
          )
        }
        for (k in seq_len(k_states)) {
          # p(S_{t+1} = j | S_t = i) * p(y_{t+1} | S_{t+1} = j) * \beta_j(t)
          logBeta[, k] <- rowLogSumExps(
            t(transProbs[idxTP[k, ], , tt]) + omegaBeta
          )
        }

        # running normalization
        logBeta <- sweep(logBeta, 1, apply(logBeta, 1, max))

        # Update probs using bayes rule: P(X_t|y) = P(y | X_t) * P(X_t | y) / P(y)
        # \beta_i(t) * \alpha_i(t)
        gammat <- logBeta + p[, tt, ]
        p[, tt, ] <- exp(sweep(gammat, 1, rowLogSumExps(gammat)))

        # sample tt HS from marginal distribution
        if (smooth_probs_states == "marginal") {
          zSample[, tt] <- apply(
            p[, tt, ], 1,
            \(x) sample.int(k_states, 1, prob = x)
          )
        }
      }

    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smooth_probs_states != "none") {
      output[["smoothProbs"]][["zSample"]] <- zSample
    }

    output[["logLikRec"]] <- logLikRec
    output[["logLikCond"]] <- logLikCond
  }

  return(output)
}


#' Get Draws from the Posterior Distribution of a HMM-DyNAM
#'
#' @inheritParams hmm_post_processing
#' @param type a character specifying weather the draws are shaped for the
#' `label.switching` package or a simple matrix format.
#'
#' @return a list with draws and additional information
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' data("Social_Evolution")
#' callNetwork <- make_network(nodes = actors, directed = TRUE) |>
#'   link_events(changeEvent = calls, nodes = actors)
#' callsDependent <- define_dependent_events(
#'   events = calls, nodes = actors, defaultNetwork = callNetwork
#' )
#' socialEvolutionData <- make_data(callsDependent)
#' data_to_stan <- make_data_hmm(
#'   rate_effects = callsDependent ~ indeg + outdeg,
#'   choice_effects = NULL,
#'   data = socialEvolutionData
#' )
#'
#' drawsMCMC <- hmm_draws_to_label_switching(data_to_stan, cmdstanSamples)
#' }
hmm_draws_to_label_switching <- function(
    data_to_stan,
    cmdstan_samples,
    k_states = data_to_stan$dataStan$kS,
    type = c("2label.switching", "draws_matrix", "draws_array", "draws_df"),
    rescale = FALSE) {
  type <- match.arg(type)
  stopifnot(
    inherits(cmdstan_samples, "CmdStanMCMC"),
    inherits(data_to_stan, "goldfish.latent.data"),
    is.numeric(k_states) && length(k_states) == 1 && k_states >= 2
  )

  if (rescale && type %in% c("2label.switching", "draws_array")) {
    stop("Not rescaling available for ", dQuote("type"), " = ",
      type,
      call. = FALSE
    )
  }

  model <- attr(data_to_stan, "model")
  subModel <- attr(data_to_stan, "subModel")

  # extract draws and reformat for posterior computations
  theta <- ifelse(model == "DNHMM", "theta", "ta")

  parmsKeep <- c(
    "lp__", theta, "pi1",
    if (subModel %in% c("both", "choice")) "betaChoice",
    if (subModel %in% c("rate")) "beta",
    if (subModel %in% c("both")) "betaRate",
    if (model == "DNCHMMRE")
      c(
        if (subModel %in% c("both", "choice"))
          c("gammaChoice", "L_OmegaChoice", "L_sigmaChoice"),
        if (subModel %in% c("rate"))
          c("gamma", "L_Omega", "L_sigma"),
        if (subModel %in% c("both"))
          c("gammaRate", "L_OmegaRate", "L_sigmaRate")
      )
  )
  draws <- cmdstan_samples$draws(
    parmsKeep,
    format = if (type == "draws_df") "draws_df" else "matrix"
  )

  namesTheta <- grep(paste0("^", theta), colnames(draws))
  idxTheta <- matrix(
    namesTheta,
    nrow = k_states, ncol = length(namesTheta) / k_states
  )
  idxEmission <- grep("^pi1", colnames(draws))

  output <- list(
    draws = draws,
    idxEmission = idxEmission,
    idxTheta = idxTheta,
    chains = cmdstan_samples$num_chains(),
    iterSampling = cmdstan_samples$metadata()[["iter_sampling"]]
  )

  if (subModel %in% c("both", "choice")) {
    idxBetaChoice <- matrix(
      grep("^betaChoice[[]", colnames(draws)),
      nrow = k_states, ncol = data_to_stan$dataStan$Pchoice
    )
    output[["idxBetaChoice"]] <- idxBetaChoice
    kPchoice <- data_to_stan$dataStan$Pchoice

    if (model == "DNCHMMRE") {
      Qchoice <- data_to_stan$dataStan$Qchoice
      idxGammaChoice <- array(
        grep("^gammaChoice[[]", colnames(draws)),
        dim = c(k_states, data_to_stan$dataStan$G, Qchoice)
      )
      output[["idxGammaChoice"]] <- idxGammaChoice

      idxLOmegaChoice <- array(
        grep("^L_OmegaChoice[[]", colnames(draws)),
        dim = c(k_states, Qchoice, Qchoice)
      )
      output[["idxLOmegaChoice"]] <- idxLOmegaChoice

      idxLSigmaChoice <- matrix(
        grep("^L_sigmaChoice[[]", colnames(draws)),
        nrow = k_states, ncol = Qchoice
      )
      output[["idxLSigmaChoice"]] <- idxLSigmaChoice
    }
  } else {
    kPchoice <- 0
  }

  if (subModel %in% c("both", "rate")) {
    kPrate <- data_to_stan$dataStan$Prate
    idxBetaRate <- matrix(
      if (subModel == "rate") {
        grep("^beta[[]", colnames(draws))
      } else {
        grep("^betaRate[[]", colnames(draws))
      },
      nrow = k_states, ncol = kPrate
    )
    output[["idxBetaRate"]] <- idxBetaRate

    if (model == "DNCHMMRE") {
      Qrate <- data_to_stan$dataStan$Qrate
      idxGammaRate <- array(
        grep("^gammaRate[[]", colnames(draws)),
        dim = c(k_states, data_to_stan$dataStan$G, Qrate)
      )
      output[["idxGammaRate"]] <- idxGammaRate

      idxLOmegaRate <- array(
        grep("^L_OmegaRate[[]", colnames(draws)),
        dim = c(k_states, Qrate, Qrate)
      )
      output[["idxLOmegaRate"]] <- idxLOmegaRate

      idxLSigmaRate <- matrix(
        grep("^L_sigmaRate[[]", colnames(draws)),
        nrow = k_states, ncol = Qrate
      )
      output[["idxLSigmaRate"]] <- idxLSigmaRate
    }
  } else {
    kPrate <- 0
  }

  nDraws <- nrow(draws)

  # # rescale draws from betas if matrix of data frame and scale stats exists
  if (type %in% c("draws_matrix", "draws_df") &&
    !(rescale && !is.null(data_to_stan[["scaleStats"]]))) {
    return(output)
  } else if (type %in% c("draws_matrix", "draws_df")) {
    output[["draws"]] <- as.data.frame(output[["draws"]])
    if (subModel %in% c("both", "choice")) {
      for (state in seq_len(k_states)) {
        output[["draws"]][, idxBetaChoice[state, ]] <- t(apply(
          output[["draws"]][, idxBetaChoice[state, ]],
          1,
          RescaleCoefs,
          scaleStats = data_to_stan[["scaleStats"]],
          isRate = FALSE
        ))

        if (model == "DNCHMMRE") {
          Qchoice <- data_to_stan$dataStan$Qchoice
          V1choice <- data_to_stan$dataStan$V1choice
          stdDev <- data_to_stan$scaleStats$choice$`scaled:scale`[V1choice]
          for (gg in seq_len(data_to_stan$dataStan$G))
            output[["draws"]][, idxGammaChoice[state, gg, ]] <- sweep(
              output[["draws"]][, idxGammaChoice[state, gg, ]],
              2,
              stdDev,
              FUN = "/"
            )
# (\(x) {
#               omega <- diag(x[seq(Q^2 + 1, length(x))]) %*% matrix(x[seq(Q^2)], nrow = Q, ncol = Q)
#               diag(1 / stdDev) %*% omega %*% t(omega) %*% diag(1 / stdDev)
#               })(ver[1, ])

#           output[["draws"]][, idxLOmegaChoice[state, , ]] <- t(apply(
#             output[["draws"]][, c(idxLOmegaChoice[state, , ], idxLSigmaChoice[state, ])],
#             1,
#             \(x, stdDev = stdDev) {
#               omega <- diag(x[seq(Q^2 + 1, length(x))]) %*% matrix(x[seq(Q^2)], nrow = Q, ncol = Q)
#               diag(1 / stdDev) %*% omega %*% t(omega) %*% diag(1 / stdDev)
#               }
#           ))
#           output[["draws"]][, idxLSigmaChoice[state, , ]] <- t(apply(
#             output[["draws"]][, idxLSigmaChoice[state, , ]],
#             1,
#             RescaleCoefs,
#             scaleStats = data2Stan[["scaleStats"]],
#             isRate = FALSE
#           ))
        }
      }
    }
    if (subModel %in% c("both", "rate")) {
      offSetInt <- log(
        data_to_stan$dataStan$Trate /
          (data_to_stan$dataStan$Nrate *
             mean(data_to_stan$dataStan$timespan))
      ) - log(data_to_stan$dataStan$rescalingFactor)
      for (state in seq_len(k_states)) {
        output[["draws"]][, idxBetaRate[state, ]] <- t(apply(
          output[["draws"]][, idxBetaRate[state, ]],
          1,
          RescaleCoefs,
          scaleStats = data_to_stan[["scaleStats"]],
          offset = offSetInt,
          isRate = TRUE
        ))

        if (model == "DNCHMMRE") {
          Qrate <- data_to_stan$dataStan$Qrate
          V1rate <- data_to_stan$dataStan$V1rate
          stdDev <- data_to_stan$scaleStats$rate$`scaled:scale`[V1rate]
          for (gg in seq_len(data_to_stan$dataStan$G))
            output[["draws"]][, idxGammaRate[state, gg, ]] <- sweep(
              output[["draws"]][, idxGammaRate[state, gg, ]],
              2,
              stdDev,
              FUN = "/"
            )
        }
      }
      output[["draws"]][, idxTheta] <- output[["draws"]][, idxTheta] /
        data_to_stan$dataStan$rescalingFactor


    }
    return(output)
  }

  if (type == "draws_array") {
    output[["draws"]] <- cmdstan_samples$draws(
      parmsKeep,
      format = "draws_array"
    )
    return(output)
  }

  nParms <- 1 + kPchoice + kPrate # + kStates # if theta included

  drawsReshape <- array(
    0,
    dim = c(nDraws, k_states, nParms),
    dimnames = list(
      draws = NULL,
      state = NULL,
      parm = c(
        "pi1",
        if (kPchoice) colnames(draws)[idxBetaChoice[1, ]],
        if (kPrate) colnames(draws)[idxBetaRate[1, ]]
      )
    )
  )

  for (state in seq_len(k_states)) {
    colKeep <- c(
      idxEmission[state],
      # idxTheta[state, ], # it has the label switching problem
      if (kPchoice) idxBetaChoice[state, ],
      if (kPrate) idxBetaRate[state, ]
    )
    drawsReshape[, state, ] <- draws[, colKeep]
  }

  output[["draws"]] <- drawsReshape
  output[["sjwinit"]] <- which.max(draws[, "lp__"])
  return(output)
}

bind_pp_hmm <- function(
    list_outputs,
    type = c("both", "viterbi", "smoothProbs"),
    smooth_probs_states = c("joint", "marginal", "none")#,
    #model = c("DNHMM", "DNCHMM")
  ) {
  type <- match.arg(type)
  smooth_probs_states <- match.arg(smooth_probs_states)
  #model <- match.arg(model)

  output <- list()
  outputNames <- names(list_outputs[[1]])

  if (type %in% c("both", "viterbi")) {
    output[["viterbi"]] <- lapply(list_outputs, "[[", "viterbi") |>
      (\(x) Reduce(rbind, x = x))()
  }

  if (type %in% c("both", "smoothProbs")) {
    smoothProbs <- lapply(list_outputs, "[[", "smoothProbs") |>
      lapply("[[", "prob")
    dimensions <- sapply(smoothProbs, dim)
    nrow <- sum(dimensions[1, ])
    other <- apply(dimensions[-1, ], 1, max)

    prob <- array(0, dim = c(nrow, other))
    start <- 1L
    end <- dimensions[1, 1]

    for (ii in seq_len(length(smoothProbs))) {
      prob[seq.int(start, end), , ] <- smoothProbs[[ii]]
      start <- start + dimensions[1, ii]
      if (ii < length(smoothProbs)) end <- end + dimensions[1, ii + 1]
    }

    output[["smoothProbs"]] <- list(prob = prob)

    if (smooth_probs_states != "none") {
      output[["smoothProbs"]][["zSample"]] <- lapply(
        list_outputs, "[[", "smoothProbs"
      ) |>
        lapply("[[", "zSample") |>
        (\(x) Reduce(rbind, x = x))()
    }

    namesReduce <- outputNames[!outputNames %in% c("smoothProbs", "viterbi")]

    if (length(namesReduce) > 0) {
      outputReduce <- lapply(
        namesReduce,
        \(x) lapply(
          listOutputs,
          "[[", x
        ) |>
          (\(x) {
            if (is.array(x[[1]]) & inherits(x[[1]], "array")) {
              Reduce(rbind, x = x)
            } else if (inherits(x[[1]], "numeric")) {
              Reduce(c, x = x)
            }
          })()
      )
      names(outputReduce) <- namesReduce
      output <- c(output, outputReduce)
    }
  }

  return(output)
}

LabelSwitchingData <- function(
    data2Stan, cmdstanSamples, labelSwitching,
    kStates = data2Stan$dataStan$kS) {
  draws <- HMMDraws2LS(
    data2Stan, cmdstanSamples,
    kStates = kStates, type = "draws_df", rescale = TRUE
  )

  # model <- attr(data2Stan, "model")
  subModel <- attr(data2Stan, "subModel")
  methodsPer <- colnames(labelSwitching$similarity)
  methodsPer <- methodsPer[!methods %in% "groundTruth"]

  df <- lapply(
    methodsPer,
    \(x) {
      permute_mcmc_array(
        draws = draws, permutation = labelSwitching[["permutations"]][[x]],
        k_states = kStates, sub_model = subModel
      ) |>
        posterior::as_draws() |>
        bayesplot::mcmc_trace_data() |>
        within({
          permutation <- x
        })
    }
  )
  df <- c(
    df,
    list(draws[["draws"]] |>
      posterior::as_draws() |>
      bayesplot::mcmc_trace_data() |>
      within({
        permutation <- "original"
      }))
  ) |>
    (\(x) Reduce(rbind, x))()

  return(df)
}

plotLabelSwitching <- function(
    data2Stan, cmdstanSamples, labelSwitching,
    kStates = data2Stan$dataStan$kS) {
  df <- LabelSwitchingData(
    data2Stan = data2Stan, cmdstanSamples = cmdstanSamples,
    labelSwitching = labelSwitching, kStates = kStates
  )

  # nChain <- length(unique(df[, "chain"]))
  # colorsChains <- bayesplot::color_scheme_get() |> unlist()
  # if (length(colorsChains) < nChain) {
  #   seqColors <- length(colorsChains) - seq_len(nChain) + 1
  # } else
  #   seqColors <- rep_len(seq_along(colorsChains), nChain)
  # colorsChains <- colorsChains[seqColors]

  by(
    df,
    df$parameter,
    \(x) {
      plotTo <- ggplot2::ggplot(
        x,
        ggplot2::aes(
          x = .data$iteration, y = .data$value, color = .data$chain
        )
      ) +
        ggplot2::geom_line(linewidth = 1 / 3) +
        ggplot2::scale_color_brewer("Chain") +
        ggplot2::facet_wrap(vars(.data$permutation), scales = "free") +
        ggplot2::scale_x_continuous(breaks = pretty) +
        ggplot2::labs(x = "", tag = unique(x$parameter))
    }
  )
  # ) |> lapply(print)
}

HSData <- function(
    postProcessing,
    type = c("probRibbon", "stateSample"),
    probsQuant = c(0.1, 0.25, 0.75, 0.9),
    pointEst = c("median", "mean")) {
  type <- match.arg(type)
  pointEst <- match.arg(pointEst)

  stopifnot(length(type) == 1)

  if (type == "probRibbon") {
    dataPlot <- postProcessing[["smoothProbs"]][["prob"]]
    if (is.null(dataPlot)) {
      stop("no data for plotting probability ribbons", call. = FALSE)
    }
    dimData <- dim(dataPlot)
    if (is.null(dimnames(dataPlot))) {
      dimnames(dataPlot) <- list(
        draw = seq_len(dimData[1]),
        time = seq_len(dimData[2]),
        state = seq_len(dimData[3])
      )
    }

    dataPlot <- as.data.frame.table(dataPlot, responseName = "prob")

    dataSummary <- by(
      dataPlot,
      dataPlot[, c("time", "state")],
      \(x) {
        data.frame(
          unique(x[, c("time", "state")]),
          quantile(
            x[, "prob"],
            probs = probsQuant
          ) |> setNames(c("ll", "l", "h", "hh")) |> t(),
          m = switch(pointEst,
            mean = mean(x[, "prob"]),
            median = median(x[, "prob"])
          )
        )
      }
    ) |>
      (\(x) Reduce(rbind.data.frame, x))()

    return(dataSummary)
  } else if (type == "stateSample") {
    dataPlot <- postProcessing[["smoothProbs"]][["zSample"]]
    if (is.null(dataPlot)) {
      stop("no data for plotting probability states sample", call. = FALSE)
    }
    dimData <- dim(dataPlot)
    if (is.null(dimnames(dataPlot))) {
      dimnames(dataPlot) <- list(
        draw = seq_len(dimData[1]),
        time = seq_len(dimData[2])
      )
    }

    dataPlot <- as.data.frame.table(dataPlot, responseName = "zSample")

    dataZsample <- with(
      dataPlot,
      table(time, zSample)
    ) |>
      as.data.frame.table() |>
      subset(Freq > 0)

    return(dataZsample)
  }
}


plotHS <- function(
    postProcessing,
    type = c("probRibbon", "stateSample"),
    pointEst = c("median", "mean"),
    prob = 0.5,
    probOuter = 0.8) {
  type <- match.arg(type)
  pointEst <- match.arg(pointEst)

  probsQuant <- c(
    0.5 - probOuter / 2, 0.5 - prob / 2,
    0.5 + prob / 2, 0.5 + probOuter / 2
  )

  if (type == "probRibbon") {
    dataPlot <- HSData(
      postProcessing = postProcessing, type = type,
      probsQuant = probsQuant,
      pointEst = pointEst
    )
    posDodge <- ggplot2::position_dodge(
      width = 1 / length(unique(dataPlot$state))
    )
    plotRibbon <- ggplot2::ggplot(
      dataPlot,
      ggplot2::aes(
        x = .data$time, y = .data$m,
        group = .data$state, colour = .data$state
      )
    ) +
      ggplot2::geom_point(position = posDodge) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$ll, ymax = .data$hh),
        position = posDodge
      ) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$l, ymax = .data$h),
        position = posDodge, linewidth = 1.5
      ) +
      ggplot2::labs(x = "HMM step", y = "Prob. being in state")

    return(plotRibbon)
  } else if (type == "stateSample") {
    plotSample <- ggplot2::ggplot(
      HSData(postProcessing = postProcessing, type = type,
      pointEst = pointEst),
      ggplot2::aes(x = .data$time, y = .data$zSample)
    ) +
      ggplot2::geom_point(ggplot2::aes(size = .data$Freq, alpha = .data$Freq)) +
      ggplot2::scale_size(transform = "log10") +
      ggplot2::scale_alpha(transform = "log10") +
      ggplot2::labs(x = "State", y = "HMM step")

    return(plotSample)
  }
}

#' Transform MCMC Array
#'
#' @param data2Stan a list output from `make_data_hmm`
#' @param cmdstanSamples a list output from `CmdStanMCMC`
#' @param labelSwitching a list output from `label.switching`
#' @param postProcessing a list output from `hmm_post_processing`
#' @param kStates the number of states in the HMM
#' @param method a character specifying the method to use for the permutation
#' @param rescale a logical specifying if the draws should be rescaled
#'
#' @return a list with the permuted post processed statistics
#' @export
transformMCMCArray <- function(
    data2Stan, cmdstanSamples, labelSwitching, postProcessing,
    kStates = data2Stan$dataStan$kS, method, rescale = TRUE) {
  subModel <- attr(data2Stan, "subModel")
  draws <- HMMDraws2LS(
    data2Stan, cmdstanSamples,
    kStates = kStates, type = "draws_df", rescale = rescale
  )

  methodsPer <- colnames(labelSwitching$similarity)

  stopifnot(method %in% methodsPer)

  permutation <- labelSwitching[["permutations"]][[method]]

  drawsdf <- permute_mcmc_array(
    draws = draws, permutation = permutation,
    k_states = kStates, sub_model = subModel
  )

  output <- list(draws = drawsdf)

  if (!is.null(postProcessing)) {
    postProcessigPer <- postProcessing

    isViterbi <- !is.null(postProcessing[["viterbi"]])
    isSmoothProbs <- !is.null(postProcessing[["smoothProbs"]])
    if (isViterbi) output[["viterbi"]] <- postProcessing[["viterbi"]]
    if (isSmoothProbs) {
      output[["smoothProbs"]] <- postProcessing[["smoothProbs"]]
      isSample <- !is.null(postProcessing[["smoothProbs"]][["zSample"]])
    }

    iterPermute <- apply(permutation, 1, \(x) any(x != seq_len(kStates))) |>
      which()
    for (iter in iterPermute) {
      permuteProbs <- permutation[iter, ]
      permuteIter <- order(permuteProbs)
      if (isViterbi) {
        output[["viterbi"]][iter, ] <-
          permuteIter[postProcessing[["viterbi"]][iter, ]]
      }
      if (isSmoothProbs) {
        output[["smoothProbs"]][["prob"]][iter, , ] <-
          postProcessing[["smoothProbs"]][["prob"]][iter, , ][, permuteProbs]
        if (isSample) {
          output[["smoothProbs"]][["zSample"]][iter, ] <-
            permuteIter[postProcessing[["smoothProbs"]][["zSample"]][iter, ]]
        }
      }
    }
  }
  return(output)
}

#' Permute MCMC Array
#'
#' @param draws a list output of `HMMDraws2LS`
#' @param permutation a numeric vector output from the `label.switching`
#'  package with the permutation to apply to each draw
#' @param kStates the number of states in the HMM
#' @param subModel a character specifying the submodel of the HMM
#' @param hmmType a character specifying the type of HMM between
#'  discrete and continuous
#'
#' @return a data frame with the permuted draws
#' @export
permute_mcmc_array <- function(
    draws, permutation, k_states,
    model = c("DNHMM", "DNCHMM", "DNCHMMRE"), sub_model = c("choice", "rate", "both"),
    hmm_type = c("discrete", "continuous")) {

  model <- match.arg(model)
  sub_model <- match.arg(sub_model)
  hmm_type <- match.arg(hmm_type)
  output <- as.data.frame(draws$draws)
  draws$draws <- as.data.frame(draws$draws)

  iterPermute <- apply(permutation, 1, \(x) any(x != seq_len(k_states))) |>
    which()

  if (length(iterPermute) == 0) {
    return(output)
  }

  origin <- c(
    draws$idxEmission,
    draws$idxTheta |> as.vector(),
    if (subModel %in% c("both", "choice")) {
      c(draws$idxBetaChoice |> as.vector(),
      if (model == "DNCHMMRE") {
        c(draws$idxGammaChoice |> as.vector(),
        draws$idxLOmegaChoice |> as.vector(),
        draws$idxLSigmaChoice |> as.vector())
      }
      )
    },
    if (sub_model %in% c("both", "rate")) {
      c(draws$idxBetaRate |> as.vector(),
      if (model == "DNCHMMRE") {
        c(draws$idxGammaRate |> as.vector(),
        draws$idxLOmegaRate |> as.vector(),
        draws$idxLSigmaRate |> as.vector())
      }
      )
    }
  )
  for (iter in iterPermute) {
    if (hmm_type == "discrete") {
      colPermutation <-
        draws$idxTheta[permutation[iter, ], permutation[iter, ]] |> as.vector()
    } else {
      colPermutation <- draws$idxTheta
      permIter <- permutation[iter, ]

      for (row in seq_along(permIter)) {
        for (col in seq_along(permIter)) {
          if (row == col) next
          permCol <- permIter[col]
          jjPos <- ifelse(permIter[row] < permCol, permCol - 1, permCol)
          jjNewPos <- ifelse(row < col, col - 1, col)
          colPermutation[row, jjNewPos] <- draws$idxTheta[permIter[row], jjPos]
        }
      }
      colPermutation <- colPermutation |> as.vector()
    }

    if (sub_model %in% c("both", "choice")) {
      choicePermutation <- draws$idxBetaChoice[permutation[iter, ], ] |>
        as.vector()

      permIter <- permutation[iter, ]
      if (model == "DNCHMMRE") {
        choicePermutation <- c(
          choicePermutation,
          draws$idxGammaChoice[permIter, , ] |> as.vector(),
          draws$idxLOmegaChoice[permIter, , ] |> as.vector(),
          draws$idxLSigmaChoice[permIter, ] |> as.vector()
        )
      }
    } else choicePermutation <- NULL
    if (sub_model %in% c("both", "rate")) {
      ratePermutation <- draws$idxBetaRate[permIter, ] |> as.vector()
      if (model == "DNCHMMRE") {
        ratePermutation <- c(
          ratePermutation,
          draws$idxGammaRate[permIter, , ] |> as.vector(),
          draws$idxLOmegaRate[permIter, , ] |> as.vector(),
          draws$idxLSigmaRate[permIter, ] |> as.vector()
        )
      }
    } else ratePermutation <- NULL

    permuted <- c(
      draws$idxEmission[permutation[iter, ]],
      colPermutation,
      choicePermutation,
      ratePermutation
    )

    output[iter, origin] <- draws$draws[iter, permuted]
  }
  return(output)
}
