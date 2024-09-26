#' Create data for Stan Hidden Markov Model model
#'
#' The first step is create the data in the structure expected by the `Stan`
#' code designed for the HMM-DyNAM.
#' Additional information of the effects used during preprocessing is added to
#' the return object.
#'
#' @param rateEffects a `formula` specification as in [goldfish::estimate()].
#' @param choiceEffects a `formula` specification as in [goldfish::estimate()].
#' @param model Current version only support `"DyNAM"` model, enhancements
#' on the code would allow to use `"REM"` model too.
#' @param scale logical value. Whether to standardize the effect stats at
#' the end.
#' @param envir an `environment` where the formula objects are available.
#' @inheritParams CreateData
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
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' data("Social_Evolution")
#' callNetwork <- defineNetwork(nodes = actors, directed = TRUE) |>
#'   linkEvents(changeEvent = calls, nodes = actors)
#' callsDependent <- defineDependentEvents(
#'   events = calls, nodes = actors, defaultNetwork = callNetwork
#' )
#' data2stan <- CreateDataHMM(
#'   rateEffects = callsDependent ~ indeg + outdeg,
#'   choiceEffects = callsDependent ~ recip + trans
#' )
#' }
CreateDataHMM <- function(
    rateEffects,
    choiceEffects,
    model = c("DyNAM", "REM"),
    kStates = 3,
    supportConstraint = NULL,
    preprocessArgs = NULL,
    scale = FALSE,
    progress = getOption("progress"),
    envir = new.env()) {
  ### 0. check parameters----
  model <- match.arg(model)

  stopifnot(
    is.null(progress) || inherits(progress, "logical"),
    is.null(preprocessArgs) ||
      inherits(preprocessArgs, "list"),
    is.null(rateEffects) || inherits(rateEffects, "formula"),
    is.null(choiceEffects) || inherits(choiceEffects, "formula"),
    is.null(supportConstraint) ||
      inherits(supportConstraint, "formula"),
    is.numeric(kStates) && length(kStates) == 1 && kStates >= 2
  )

  # setting initial values of some arguments
  if (is.null(progress)) progress <- FALSE

  if (!is.null(rateEffects) && !is.null(choiceEffects)) {
    subType <- "both"


    if (as.character(choiceEffects[[2]]) != as.character(rateEffects[[2]])) {
      stop("dependent event network needs to be the same for both formulas")
    }
  } else if (!is.null(rateEffects)) {
    subType <- "rate"
  } else if (!is.null(choiceEffects)) {
    subType <- "choice"
  } else {
    stop(
      dQuote("rateEffects"), " and ", dQuote("choiceEffects"),
      " arguments are NULL objects. Specify at least one model."
    )
  }
  # process data
  scaleStats <- if (scale) list() else NULL
  if (!is.null(rateEffects)) {
    dataProcessedRate <- GatherPreprocessing(
      formula = rateEffects,
      model = model,
      subModel = "rate",
      preprocessArgs = preprocessArgs,
      progress = progress,
      envir = new.env(parent = envir)
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


    dataStanRate <- within(dataProcessedRate, {
      Trate <- length(n_candidates)
      Nrate <- nrow(stat_all_events)
      Prate <- ncol(stat_all_events)
      startRate <- cumsum(c(1, head(n_candidates, -1)))
      endRate <- cumsum(n_candidates)
      namesEffects <- gsub("\\$", "Of", namesEffects)
      Xrate <- setNames(stat_all_events, namesEffects)
      if (hasIntercept) {
        choseRate <- selected[, 1] + (startRate - 1) * isDependent
        offsetInt <- log(Trate / (sum(timespan) * mean(n_candidates))) # mean(timespan)
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
      scaleStats[["rate"]] <- attributes(Xrate)

      if (dataStanRate$hasIntercept) {
        dataStanRate$Xrate[, -1] <- Xrate[, ]
      } else {
        dataStanRate$Xrate <- Xrate[, ]
      }
    }
    # table(idxEvents[1, ] == dataStanRate$startRate)
    # table(idxEvents[2, ] == dataStanRate$endRate)
    # chose <- which(expandedDF$selected)
    # table(chose == dataStanRate$choseRate[dataStanRate$isDependent])
  } else {
    dataProcessedRate <- dataStanRate <- NULL
  }

  if (!is.null(choiceEffects)) {
    feTerms <- terms(choiceEffects)
    cstrTerms <- if (!is.null(supportConstraint)) {
      terms(supportConstraint)
    } else {
      NULL
    }

    if (length(attr(cstrTerms, "term.labels")) > 1) {
      stop(dQuote("supportConstraint"), " argument only works for one effect.")
    }

    if (any(attr(feTerms, "order") != 1)) {
      stop(
        dQuote("choiceEffects"),
        "formula argument doesn't support interactions yet"
      )
    }

    termsDyNAM <- c(
      attr(feTerms, "term.labels"),
      attr(cstrTerms, "term.labels")
    ) |> unique()

    formulaDyNAM <- reformulate(
      termsDyNAM,
      response = as.character(choiceEffects[[2]])
    )

    dataProcessedChoice <- GatherPreprocessing(
      formula = formulaDyNAM,
      model = model,
      subModel = "choice",
      preprocessArgs = preprocessArgs,
      progress = progress,
      envir = new.env(parent = envir)
    )

    nEvents <- length(dataProcessedChoice$sender)

    namesEffects <- setNames(
      gsub("\\$", "Of", dataProcessedChoice$namesEffects),
      termsDyNAM
    )

    expandedDF <- cbind(
      setNames(
        as.data.frame(dataProcessedChoice$stat_all_events),
        namesEffects
      ),
      data.frame(
        event = rep(seq_len(nEvents), dataProcessedChoice$n_candidates),
        selected = sequence(dataProcessedChoice$n_candidates) ==
          rep(dataProcessedChoice$selected, dataProcessedChoice$n_candidates)
      )
    )

    # subset if constraint
    if (!is.null(supportConstraint)) {
      cstrName <- namesEffects[attr(cstrTerms, "term.labels")]
      keep <- expandedDF[, cstrName] == 1
      expandedDF <- expandedDF[keep, !names(expandedDF) %in% cstrName]
      effectDescription <-
        dataProcessedChoice$effectDescription[!namesEffects %in% cstrName, ]
      namesEffects <- namesEffects[!namesEffects %in% cstrName]
    } else {
      effectDescription <- dataProcessedChoice$effectDescription
    }

    # create objects for Stan
    nTotal <- nrow(expandedDF)

    idxEvents <- tapply(seq_len(nTotal), expandedDF$event, range) |>
      simplify2array()

    Xmat <- as.matrix(expandedDF[, namesEffects])

    dataStanChoice <- list(
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
      scaleStats[["choice"]] <- attributes(Xchoice)

      dataStanChoice$Xchoice <- Xchoice[, ]
    }
  } else {
    namesEffects <- effectDescription <- NULL
  }

  dataStan <- c(
    if (!is.null(rateEffects)) dataStanRate,
    if (!is.null(choiceEffects)) dataStanChoice,
    list(
      kS = kStates,
      alpha = c(floor((kStates - 1) / (1 - 0.8)), 1)
    )
  )

  if (!is.null(effectDescription) &
    !is.null(dataProcessedRate$effectDescription)) {
    if (!setequal(
      colnames(effectDescription),
      colnames(dataProcessedRate$effectDescription)
    )) {
      dataProcessedRate$effectDescription <- AddColumns(
        dataProcessedRate$effectDescription,
        effectDescription
      )

      effectDescription <- AddColumns(
        effectDescription,
        dataProcessedRate$effectDescription
      )
    }
  }

  return(structure(
    list(
      dataStan = dataStan,
      namesEffects = c(
        dataProcessedRate$namesEffects,
        namesEffects
      ),
      effectDescription = rbind(
        dataProcessedRate$effectDescription,
        effectDescription
      ),
      scaleStats = scaleStats
    ),
    class = "goldfish.latent.data",
    model = "DNHMM",
    subModel = subType
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
#' @param data2Stan a object of class `"goldfish.latent.data"` outcome of a call
#' to [CreateDataHMM()].
#' @param cmdstanSamples a object of class `"CmdStanMCMC"` with the posterior
#' samples of a HMM from the `data2Stan`.
#' @param kStates an integer with the number of states used to get samples
#' from the posterior distribution in `cmdstanSamples`.
#' @param type a character specifying which algorithm to use for the post
#' processing of the posterior distribution samples.
#' @param smoothProbsSt a character specifying the way that state sample path
#' are drawn either using the smoothed probabilities `"marginal"` or
#' conditioning of the next sample state `"joint"`.
#' Only used when `type = "both"` or `type = "smoothProbs"`.
#'
#' @return a list with the most probable path when Viterbi algorithm is used
#' and/or the smoothed probabilities and sample paths when the
#' Forward-Backward algorithm is used.
#' @export
#'
#' @examples
#' \donttest{
#' library(goldfish)
#' data("Social_Evolution")
#' callNetwork <- defineNetwork(nodes = actors, directed = TRUE) |>
#'   linkEvents(changeEvent = calls, nodes = actors)
#' callsDependent <- defineDependentEvents(
#'   events = calls, nodes = actors, defaultNetwork = callNetwork
#' )
#' data2stan <- CreateDataHMM(
#'   rateEffects = callsDependent ~ indeg + outdeg,
#'   choiceEffects = NULL
#' )
#'
#' postProcess <- HMMPostProcessing(data2stan, cmdstanSamples)
#' }
HMMPostProcessing <- function(
    dataStan, cmdstanSamples,
    kStates = dataStan$dataStan$kS,
    type = c("both", "viterbi", "smoothProbs"),
    smoothProbsSt = c("joint", "marginal", "none"),
    cl = NULL) {
  type <- match.arg(type)
  smoothProbsSt <- match.arg(smoothProbsSt)

  model <- attr(dataStan, "model")
  stopifnot(
    inherits(cmdstanSamples, "CmdStanMCMC"),
    inherits(dataStan, "goldfish.latent.data"),
    is.numeric(kStates) && length(kStates) == 1 && kStates >= 2,
    model %in% c("DNHMM", "DNCHMM", "DNCHMMRE")
  )

  subModel <- attr(dataStan, "subModel")

  # extract draws and reformat for posterior computations
  typeOutput <- ifelse(is.null(cl), "draws_matrix", "draws_array")

  drawsObject <- HMMDraws2LS(
    dataStan, cmdstanSamples,
    kStates = kStates, type = typeOutput
  )

  dataStan <- dataStan$dataStan

  if (is.null(cl)) {
    output <- do.call(
      switch(
        model,
         "DNHMM" = "HMMPPperDraws",
         "DNCHMM" = "CHMMPPperDraws",
         "DNCHMMRE" = "CHMMREPPperDraws"
      ),
      args = list(
        chainIter = NULL,
        drawsObject = drawsObject,
        dataStan = dataStan,
        model = model,
        subModel = subModel,
        kStates = kStates,
        type = type,
        smoothProbsSt = smoothProbsSt
      )
    )
  } else {
    on.exit(parallel::stopCluster(cl))
    ignore <- parallel::clusterEvalQ(cl, {
      library(matrixStats)
      library(expm)
      NULL
    })
    output <- parallel::clusterApplyLB(
      cl = cl,
      seq_len(cmdstanSamples$num_chains()),
      fun = switch(
        model,
         "DNHMM" = HMMPPperDraws,
         "DNCHMM" = CHMMPPperDraws,
         "DNCHMMRE" = CHMMREPPperDraws
      ),
      drawsObject = drawsObject,
      dataStan = dataStan,
      model = model,
      subModel = subModel,
      kStates = kStates,
      type = type,
      smoothProbsSt = smoothProbsSt
    ) |>
      bindPPHMM(type = type, smoothProbsSt = smoothProbsSt)
  }

  return(output)
}

HMMPPperDraws <- function(
    chainIter, drawsObject, dataStan, model, subModel,
    kStates, type, smoothProbsSt) {
  # init output
  output <- list()

  draws <- drawsObject$draws
  idxTheta <- drawsObject$idxTheta
  idxEmission <- drawsObject$idxEmission
  idxBetaChoice <- drawsObject$idxBetaChoice
  idxBetaRate <- drawsObject$idxBetaRate

  if (!is.null(chainIter)) {
    draws <- draws[, chainIter, ]
  }


  nDraws <- nrow(draws)

  # define sizes
  nEvents <- ifelse(subModel %in% c("both", "rate"), "Trate", "Tchoice")
  nEvents <- dataStan[[nEvents]]

  isRes <- !is.null(dataStan$Nres)
  TT <- ifelse(isRes, dataStan$Nres, nEvents)


  if (subModel %in% c("choice")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xb <- tcrossprod(dataStan$Xchoice, draws[, idxBetaChoice[kS, ]])

      for (event in seq_len(nEvents)) {
        llEventState[event, , kS] <-
          xb[dataStan$choseChoice[event], ] -
          colLogSumExps(
            xb,
            rows = seq(
              dataStan$startChoice[event],
              dataStan$endChoice[event]
            )
          )
      }
    }
  }

  if (subModel %in% c("rate")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xb <- tcrossprod(dataStan$Xrate, draws[, idxBetaRate[kS, ]]) +
        dataStan$offsetInt

      for (event in seq_len(nEvents)) {
        llEventState[event, , kS] <-
          ifelse(dataStan$isDependent[event],
            xb[dataStan$choseRate[event], ], 0
          ) -
          dataStan$timespan[event] * exp(colLogSumExps(
            xb,
            rows = seq(
              dataStan$startRate[event],
              dataStan$endRate[event]
            )
          ))
      }
    }
  }

  if (subModel %in% c("both")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xbR <- tcrossprod(dataStan$Xrate, draws[, idxBetaRate[kS, ]]) +
        dataStan$offsetInt
      xbC <- tcrossprod(dataStan$Xchoice, draws[, idxBetaChoice[kS, ]])

      eventChoice <- 1L
      for (event in seq_len(nEvents)) {
        loglik <- -dataStan$timespan[event] * exp(colLogSumExps(
          xbR,
          rows = seq(
            dataStan$startRate[event],
            dataStan$endRate[event]
          )
        ))
        if (dataStan$isDependent[event]) {
          loglik <- loglik + xbR[dataStan$choseRate[event], ] +
            xbC[dataStan$choseChoice[eventChoice], ] -
            colLogSumExps(
              xbC,
              rows = seq(
                dataStan$startChoice[eventChoice],
                dataStan$endChoice[eventChoice]
              )
            )
          eventChoice <- eventChoice + 1L
        }

        llEventState[event, , kS] <- loglik
      }
    }
  }

  if (isRes) {
    llEventState <- array(apply(
      llEventState,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(dataStan$resA[y], dataStan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, nDraws, kStates))
  }

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, kStates, nDraws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, kStates, nDraws))

    # forward past computing most likely state from previous state
    # first observation: p(y_1| z_1) * p(z_1) (emission prob)
    delta[1, , ] <- t(log(draws[, idxEmission]) + llEventState[1, , ])

    for (tt in seq(2, TT)) {
      prevDelta <- t(delta[tt - 1, , ])
      for (k in seq_len(kStates)) {
        T_1_xj <- prevDelta + log(draws[, idxTheta[, k]]) +
          llEventState[tt, , k]

        delta[tt, k, ] <- apply(T_1_xj, 1, max)
        bpointer[tt, k, ] <- apply(T_1_xj, 1, which.max)
      }
    }

    # backward past
    z <- array(0L,
      dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT))
    )

    z[, TT] <- apply(delta[TT, , ], 2, which.max)

    for (tt in seq(TT - 1, 1)) {
      z[, tt] <- sapply(
        seq_len(nDraws),
        \(x) bpointer[tt + 1, z[x, tt + 1], x]
      )
    }

    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(nDraws, TT, kStates),
      dimnames = list(
        draws = seq_len(nDraws), time = seq_len(TT),
        states = seq_len(kStates)
      )
    )

    if (smoothProbsSt != "none") zSample <- array(0L, dim = c(nDraws, TT))

    # forward -- filtering: alpha [tt, k] with running normalization
    #  11.2.2, Finite Mixture and Markov Switching models,
    #  Frühwirth, S., 2006

    # # filter at first observation:
    # # p(S_1 = k|y_0) = p(y_1| S_1) * p(S_1); (emission prob)
    p[, 1, ] <- log(draws[, idxEmission]) + llEventState[1, , ]
    # In Frühwirth: normalization is over llEventState, here follow Stan
    p[, 1, ] <- sweep(p[, 1, ], 1, apply(p[, 1, ], 1, max))
    # not need to convert to prob
    # p[, 1, ] <- sweep(p[, 1, ], 1, rowLogSumExps(p[, 1, ]))

    for (tt in seq(2, TT)) {
      # # one-step ahead prediction of S_t:
      # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
      for (k in seq_len(kStates)) {
        p[, tt, k] <- rowLogSumExps(
          log(draws[, idxTheta[, k]]) + p[, tt - 1, ]
        ) +
          # # filter for S_t: p(S_t = k| y_t) =
          # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
          # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) ; unnormalize enough
          llEventState[tt, , k]
      }
    }

    # Normalize last value, already smooth distribution
    p[, TT, ] <- exp(sweep(p[, TT, ], 1, rowLogSumExps(p[, TT, ])))

    if (smoothProbsSt != "none") { # sample last Hidden State (HS)
      zSample[, TT] <- apply(
        p[, TT, ], 1,
        \(x) sample.int(kStates, 1, prob = x)
      )
    }

    # backward: smoother suggested in Hamilton expresses these as marginal
    # probabilities from the joint distribution of S_t and S_T | y
    # Implementation follows Stan hmm_hidden_state_prob()
    #

    # initial ending state ass as given (uniform)
    logBeta <- array(0, dim = c(nDraws, kStates))

    for (tt in seq(TT - 1, 1)) {
      # # Baum-Welch alg
      # #
      omegaBeta <- logBeta + llEventState[tt + 1, , ] # element-wise product

      # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
      if (smoothProbsSt == "joint") {
        probLastHS <- p[, tt, ] +
          t(vapply(
            seq_len(nDraws),
            \(x){
              lastHS <- zSample[x, tt + 1]
              log(draws[x, idxTheta[, lastHS]]) + omegaBeta[x, lastHS]
            },
            numeric(2)
          ))
        probLastHS <- exp(sweep(probLastHS, 1, rowLogSumExps(probLastHS)))
        zSample[, tt] <- apply(
          probLastHS, 1,
          \(x) sample.int(kStates, 1, prob = x)
        )
      }
      for (k in seq_len(kStates)) {
        logBeta[, k] <- rowLogSumExps(
          log(draws[, idxTheta[k, ]]) + omegaBeta
        )
      }

      # running normalization
      logBeta <- sweep(logBeta, 1, apply(logBeta, 1, max))

      #
      gammat <- logBeta + p[, tt, ]
      p[, tt, ] <- exp(sweep(gammat, 1, rowLogSumExps(gammat)))

      # sample tt HS from marginal distribution
      if (smoothProbsSt == "marginal") {
        zSample[, tt] <- apply(
          p[, tt, ], 1,
          \(x) sample.int(kStates, 1, prob = x)
        )
      }
    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smoothProbsSt != "none") {
      output[["smoothProbs"]][["zSample"]] <- zSample
    }
  }

  return(output)
}

#' @importFrom expm expm
CHMMPPperDraws <- function(
    chainIter, drawsObject, dataStan, model, subModel,
    kStates, type, smoothProbsSt) {
  # init output
  output <- list()

  draws <- drawsObject$draws
  idxTheta <- drawsObject$idxTheta
  idxEmission <- drawsObject$idxEmission
  idxBetaChoice <- drawsObject$idxBetaChoice
  idxBetaRate <- drawsObject$idxBetaRate

  if (!is.null(chainIter)) {
    draws <- draws[, chainIter, , drop = TRUE]
  }


  nDraws <- nrow(draws)

  # define sizes
  nEvents <- ifelse(subModel %in% c("both", "rate"), "Trate", "Tchoice")
  nEvents <- dataStan[[nEvents]]

  isRes <- !is.null(dataStan$Nres)
  TT <- ifelse(isRes, dataStan$Nres, nEvents)

  logCrudeRate <- ifelse(
    subModel %in% c("both", "rate"),
    log(dataStan$Trate / (dataStan$Nrate * mean(dataStan$timespan))), 0
  )

  if (subModel %in% c("choice")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xb <- tcrossprod(dataStan$Xchoice, draws[, idxBetaChoice[kS, ]])

      for (event in seq_len(nEvents)) {
        llEventState[event, , kS] <-
          xb[dataStan$choseChoice[event], ] -
          colLogSumExps(
            xb,
            rows = seq(
              dataStan$startChoice[event],
              dataStan$endChoice[event]
            )
          )
      }
    }
  }

  if (subModel %in% c("rate")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xb <- tcrossprod(dataStan$Xrate, draws[, idxBetaRate[kS, ]]) +
        logCrudeRate

      for (event in seq_len(nEvents)) {
        llEventState[event, , kS] <-
          ifelse(dataStan$isDependent[event],
            xb[dataStan$choseRate[event], ], 0
          ) -
          dataStan$timespan[event] * exp(colLogSumExps(
            xb,
            rows = seq(
              dataStan$startRate[event],
              dataStan$endRate[event]
            )
          ))
      }
    }
  }

  if (subModel %in% c("both")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xbR <- tcrossprod(dataStan$Xrate, draws[, idxBetaRate[kS, ]]) +
        logCrudeRate
      xbC <- tcrossprod(dataStan$Xchoice, draws[, idxBetaChoice[kS, ]])

      eventChoice <- 1L
      for (event in seq_len(nEvents)) {
        loglik <- -dataStan$timespan[event] * exp(colLogSumExps(
          xbR,
          rows = seq(
            dataStan$startRate[event],
            dataStan$endRate[event]
          )
        ))
        if (dataStan$isDependent[event]) {
          loglik <- loglik + xbR[dataStan$choseRate[event], ] +
            xbC[dataStan$choseChoice[eventChoice], ] -
            colLogSumExps(
              xbC,
              rows = seq(
                dataStan$startChoice[eventChoice],
                dataStan$endChoice[eventChoice]
              )
            )
          eventChoice <- eventChoice + 1L
        }

        llEventState[event, , kS] <- loglik
      }
    }
  }

  if (isRes) {
    llEventState <- array(apply(
      llEventState,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(dataStan$resA[y], dataStan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, nDraws, kStates))
  }

  # compute transition probabilities
  transProbs <- vapply(
    dataStan$timespan,
    \(y)
    apply(
      draws[, idxTheta] * y,
      1,
      \(x){
        log(expm::expm(
          matrix(
            x,
            nrow = kStates, ncol = kStates
          )
        ))
      }
    ),
    matrix(0, nrow = kStates^2, ncol = nDraws)
  )

  idxTP <- idxTheta - min(idxTheta) + 1

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, kStates, nDraws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, kStates, nDraws))

    # forward past computing most likely state from previous state
    # first observation: p(y_1, z_1) = p(y_1| z_1) * p(z_1) (emission prob)
    delta[1, , ] <- t(log(draws[, idxEmission]) + llEventState[1, , ])

    for (tt in seq(2, TT)) {
      prevDelta <- t(delta[tt - 1, , ])
      for (k in seq_len(kStates)) {
        T_1_xj <- prevDelta + t(transProbs[idxTP[, k], , tt]) +
          llEventState[tt, , k]

        delta[tt, k, ] <- apply(T_1_xj, 1, max)
        bpointer[tt, k, ] <- apply(T_1_xj, 1, which.max)
      }
    }

    # backward past
    z <- array(0L,
      dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT))
    )

    z[, TT] <- apply(delta[TT, , ], 2, which.max)

    for (tt in seq(TT - 1, 1)) {
      z[, tt] <- sapply(
        seq_len(nDraws),
        \(x) bpointer[tt + 1, z[x, tt + 1], x]
      )
    }

    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(nDraws, TT, kStates),
      dimnames = list(
        draws = seq_len(nDraws), time = seq_len(TT),
        states = seq_len(kStates)
      )
    )

    # log likelihood as sum of scaling factors, Kadhem 2015
    # those are the marginal probabilities of y_t
    logLikRec <- array(0, dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT)))

    # log likelihood conditional on the sample state, Kadhem 2015 and 2021
    logLikCond <- array(0, dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT)))

    # sampling of the last state, required for conditional logLik
    if (smoothProbsSt != "none") zSample <- array(0L, dim = c(nDraws, TT))

    # forward -- filtering: alpha [tt, k] with running normalization
    #  11.2.2, Finite Mixture and Markov Switching models,
    #  Frühwirth, S., 2006

    # # filter at first observation:
    # # p(S_1 = k, y_1) = p(y_1| S_1) * p(S_1); (emission prob)
    p[, 1, ] <- log(draws[, idxEmission]) + llEventState[1, , ]
    # Compute normalization constants Kadhem 2015 (It's just the marginal prob)
    # p(y_1) = \sum_k p(y_1, S_1 = k), marginal prob up to t, Frühwirth 2006
    logLikRec[, 1] <- rowLogSumExps(p[, 1, ])
    # Filtering for s_t: p(S_t = k| y_t) = p(y_t, S_t = k) / p(y_t)
    p[, 1, ] <- sweep(p[, 1, ], 1, logLikRec[, 1])

    # objects intermediate values
    ptt <- array(0, dim = c(nDraws, kStates))
    for (tt in seq(2, TT)) {
      # # one-step ahead prediction of S_t:
      # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
      for (k in seq_len(kStates)) {
        ptt[, k] <- rowLogSumExps(
          t(transProbs[idxTP[, k], ,tt]) + p[, tt - 1, ]) +
          # # filter for S_t: p(S_t = k| y_t) =
          # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
          # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1})
          llEventState[tt, , k]
      }
      # marginal prob (y) and filter probs P(S_t = k | y_t)
      marginalProb <- rowLogSumExps(ptt)
      logLikRec[, tt] <- marginalProb

      p[, tt, ] <- sweep(ptt, 1, marginalProb)
    }

    # Normalize last value, already smooth distribution
    p[, TT, ] <- exp(p[, TT, ])

    if (smoothProbsSt != "none") { # sample last Hidden State (HS)
      zSample[, TT] <- apply(
        p[, TT, ], 1,
        \(x) sample.int(kStates, 1, prob = exp(x))
      )

      # Compute conditional log likelihood
      logLikCond[, TT] <- vapply(
        seq_len(nDraws),
        \(x) llEventState[TT, x, zSample[x, TT]],
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
    logBeta <- array(0, dim = c(nDraws, kStates))

    for (tt in seq(TT - 1, 1)) {
      # \beta_j(t) * p(y_{t+1} | S_{t+1} = j)
      omegaBeta <- logBeta + llEventState[tt + 1, , ] # element-wise product

      # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
      # FFBS, a.k.a, Multi-move sampling
      if (smoothProbsSt == "joint") {
        probLastHS <- p[, tt, ] +
          t(vapply(
            seq_len(nDraws),
            \(x){
              lastHS <- zSample[x, tt + 1]
              # p(S_{t+1} = j|S_t = i) * p(y_{t+1}|S_{t+1} = j) * \beta_j(t)
              t(transProbs[idxTP[, lastHS], x, tt]) + omegaBeta[x, lastHS]
            },
            numeric(kStates)
          ))
        probLastHS <- exp(sweep(probLastHS, 1, rowLogSumExps(probLastHS)))
        zSample[, tt] <- apply(
          probLastHS, 1,
          \(x) sample.int(kStates, 1, prob = x)
        )

        # Compute conditional log likelihood
        logLikCond[, tt] <- vapply(
          seq_len(nDraws),
          \(x) llEventState[tt, x, zSample[x, tt]],
          numeric(1)
        )
      }
      for (k in seq_len(kStates)) {
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
      if (smoothProbsSt == "marginal") {
        zSample[, tt] <- apply(
          p[, tt, ], 1,
          \(x) sample.int(kStates, 1, prob = x)
        )
      }
    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smoothProbsSt != "none") {
      output[["smoothProbs"]][["zSample"]] <- zSample
    }

    output[["logLikRec"]] <- logLikRec
    output[["logLikCond"]] <- logLikCond
  }

  return(output)
}

CHMMREPPperDraws <- function(
    chainIter, drawsObject, dataStan, model, subModel,
    kStates, type, smoothProbsSt) {
  # init output
  output <- list()

  draws <- drawsObject$draws
  idxTheta <- drawsObject$idxTheta
  idxEmission <- drawsObject$idxEmission
  idxBetaChoice <- drawsObject$idxBetaChoice
  idxBetaRate <- drawsObject$idxBetaRate

  idxGammaChoice <- drawsObject$idxGammaChoice
  idxGammaRate <- drawsObject$idxGammaRate
  idxLOmegaChoice <- drawsObject$idxLOmegaChoice
  idxLOmegaRate <- drawsObject$idxLOmegaRate
  idxLSigmaChoice <- drawsObject$idxLSigmaChoice


  if (!is.null(chainIter)) {
    draws <- draws[, chainIter, , drop = TRUE]
  }


  nDraws <- nrow(draws)

  # define sizes
  nEvents <- ifelse(subModel %in% c("both", "rate"), "Trate", "Tchoice")
  nEvents <- dataStan[[nEvents]]

  isRes <- !is.null(dataStan$Nres)
  TT <- ifelse(isRes, dataStan$Nres, nEvents)

  logCrudeRate <- ifelse(
    subModel %in% c("both", "rate"),
    log(dataStan$Trate / (dataStan$Nrate * mean(dataStan$timespan))), 0
  )

  # indices group
  startGroup <- c(dataStan$startGroup, TT + 1)
  group <- dataStan$groupChoice

  if (subModel %in% c("choice")) {
    Xchoice <- dataStan$Xchoice
    Zchoice <- dataStan$Xchoice[, dataStan$V1choice]

    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(dataStan$startChoice[event], dataStan$endChoice[event])
        eventChose <- dataStan$choseChoice[event] - dataStan$startChoice[event] + 1
        xb <- tcrossprod(Xchoice[eventIdx, ], draws[, idxBetaChoice[kS, ]]) +
          tcrossprod(Zchoice[eventIdx, ], draws[, idxGammaChoice[kS, group[eventIdx[1]], ]])  

        llEventState[event, , kS] <- xb[eventChose, ] - colLogSumExps(xb)
      }
    }
  }

  if (subModel %in% c("rate")) {
    Xrate <- dataStan$Xrate
    Zrate <- dataStan$Xrate[, dataStan$V1rate]

    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(dataStan$startRate[event], dataStan$endRate[event])
        eventChose <- dataStan$choseRate[event] - dataStan$startRate[event] + 1
        xb <- tcrossprod(Xrate[eventIdx, ], draws[, idxBetaRate[kS, ]]) +
          tcrossprod(Zrate[eventIdx, ],
            draws[, idxGammaRate[kS, group[eventIdx[1]], ]]) +
            logCrudeRate
        llEventState[event, , kS] <-
          ifelse(dataStan$isDependent[event], xb[eventChose, ], 0) -
          dataStan$timespan[event] * exp(colLogSumExps(xb))
      }
    }
  }

  if (subModel %in% c("both")) {
    Xchoice <- dataStan$Xchoice
    Zchoice <- dataStan$Xchoice[, dataStan$V1choice]
    Xrate <- dataStan$Xrate
    Zrate <- dataStan$Xrate[, dataStan$V1rate]

    llEventState <- array(0, dim = c(nEvents, nDraws, kStates))

    for (kS in seq_len(kStates)) {
      xbR <- tcrossprod(dataStan$Xrate, draws[, idxBetaRate[kS, ]]) +
        logCrudeRate
      xbC <- tcrossprod(dataStan$Xchoice, draws[, idxBetaChoice[kS, ]])

      eventChoice <- 1L
      for (event in seq_len(nEvents)) {
        eventIdx <- seq(dataStan$startRate[event], dataStan$endRate[event])
        eventChose <- dataStan$choseRate[event] - dataStan$startRate[event] + 1
        xbR <- tcrossprod(Xrate[eventIdx, ], draws[, idxBetaRate[kS, ]]) +
          tcrossprod(Zrate[eventIdx, ],
            draws[, idxGammaRate[kS, group[eventIdx[1]], ]]) +
            logCrudeRate
        loglik <- -dataStan$timespan[event] * exp(colLogSumExps(xbR))
        if (dataStan$isDependent[event]) {
          loglik <- loglik + xbR[eventChose, ] 
          
          eventIdx <- seq(dataStan$startChoice[eventChoice], dataStan$endChoice[eventChoice])
          eventChose <- dataStan$choseChoice[eventChoice] - dataStan$startChoice[eventChoice] + 1
          
          xbC <- tcrossprod(Xchoice[eventIdx, ], draws[, idxBetaChoice[kS, ]]) +
            tcrossprod(Zchoice[eventIdx, ],
              draws[, idxGammaChoice[kS, group[eventIdx[1]], ]])

          loglik <- loglik + xbC[eventChose, ] - colLogSumExps(xbC)
          eventChoice <- eventChoice + 1L
        }

        llEventState[event, , kS] <- loglik
      }
    }
  }

  if (isRes) {
    llEventState <- array(apply(
      llEventState,
      3,
      \(x) lapply(
        seq_len(TT),
        \(y) colSums2(
          x,
          rows = seq(dataStan$resA[y], dataStan$resA[y + 1] - 1)
        )
      ) |> (\(x) Reduce(rbind, x = x))()
    ), dim = c(TT, nDraws, kStates))
  }

  # compute transition probabilities
  transProbs <- vapply(
    dataStan$timesSender,
    \(y)
    apply(
      draws[, idxTheta] * y,
      1,
      \(x){
        log(expm::expm(
          matrix(
            x,
            nrow = kStates, ncol = kStates
          )
        ))
      }
    ),
    matrix(0, nrow = kStates^2, ncol = nDraws)
  )

  idxTP <- idxTheta - min(idxTheta) + 1

  if (type %in% c("both", "viterbi")) {
    # back-pointer to the most likely previous state on the most probable path
    bpointer <- array(0, dim = c(TT, kStates, nDraws))
    # max prob for the sequence up to t that ends with an emission from state k
    delta <- array(0, dim = c(TT, kStates, nDraws))

    # forward past
    z <- array(0L,
      dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT))
    )

    for (gg in seq_len(dataStan$G)) {
      nEventsG <- startGroup[gg + 1] - startGroup[gg]
      eventIdx <- seq(startGroup[gg], startGroup[gg + 1] - 1)

      # forward past computing most likely state from previous state
      # first observation: p(y_1, z_1) = p(y_1| z_1) * p(z_1) (emission prob)
      delta[eventIdx[1], , ] <- t(log(draws[, idxEmission]) + llEventState[eventIdx[1], , ])

      for (tt in tail(eventIdx, -1)) {
        prevDelta <- t(delta[tt - 1, , ])
        for (k in seq_len(kStates)) {
          T_1_xj <- prevDelta + t(transProbs[idxTP[, k], , tt]) +
            llEventState[tt, , k]

          delta[tt, k, ] <- apply(T_1_xj, 1, max)
          bpointer[tt, k, ] <- apply(T_1_xj, 1, which.max)
        }
      }

      # backward past
      z[, tail(eventIdx, 1)] <- apply(delta[tail(eventIdx, 1), , ], 2, which.max)

      for (tt in rev(head(eventIdx, -1))) {
        z[, tt] <- sapply(
          seq_len(nDraws),
          \(x) bpointer[tt + 1, z[x, tt + 1], x]
        )
      }
    }
    output[["viterbi"]] <- z
  }

  if (type %in% c("both", "smoothProbs")) {
    p <- array(
      0,
      dim = c(nDraws, TT, kStates),
      dimnames = list(
        draws = seq_len(nDraws), time = seq_len(TT),
        states = seq_len(kStates)
      )
    )

    # log likelihood as sum of scaling factors, Kadhem 2015
    # those are the marginal probabilities of y_t
    logLikRec <- array(0, dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT)))

    # log likelihood conditional on the sample state, Kadhem 2015 and 2021
    logLikCond <- array(0, dim = c(nDraws, TT),
      dimnames = list(draws = seq_len(nDraws), time = seq_len(TT)))

    # sampling of the last state, required for conditional logLik
    if (smoothProbsSt != "none") zSample <- array(0L, dim = c(nDraws, TT))

    for (gg in seq_len(dataStan$G)) {
      nEventsG <- startGroup[gg + 1] - startGroup[gg]
      eventIdx <- seq(startGroup[gg], startGroup[gg + 1] - 1)

      # forward -- filtering: alpha [tt, k] with running normalization
      #  11.2.2, Finite Mixture and Markov Switching models,
      #  Frühwirth, S., 2006

      # # filter at first observation:
      # # p(S_1 = k, y_1) = p(y_1| S_1) * p(S_1); (emission prob)
      p[, eventIdx[1], ] <- log(draws[, idxEmission]) + llEventState[eventIdx[1], , ]
      # Compute normalization constants Kadhem 2015 (It's just the marginal prob)
      # p(y_1) = \sum_k p(y_1, S_1 = k), marginal prob up to t, Frühwirth 2006
      logLikRec[, eventIdx[1]] <- rowLogSumExps(p[, eventIdx[1], ])
      # Filtering for s_t: p(S_t = k| y_t) = p(y_t, S_t = k) / p(y_t)
      p[, eventIdx[1], ] <- sweep(p[, eventIdx[1], ], 1, logLikRec[, eventIdx[1]])

      # objects intermediate values
      ptt <- array(0, dim = c(nDraws, kStates))
      for (tt in tail(eventIdx, -1)) {
        # # one-step ahead prediction of S_t:
        # # p(S_t = k | y_{t-1}) = \sum_l \theta_{lk} p(S_{t-1} = l| y_{t-1})
        for (k in seq_len(kStates)) {
          ptt[, k] <- rowLogSumExps(
            t(transProbs[idxTP[, k], ,tt]) + p[, tt - 1, ]) +
            # # filter for S_t: p(S_t = k| y_t) =
            # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
            # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1})
            llEventState[tt, , k]
        }
        # marginal prob (y) and filter probs P(S_t = k | y_t)
        marginalProb <- rowLogSumExps(ptt)
        logLikRec[, tt] <- marginalProb

        p[, tt, ] <- sweep(ptt, 1, marginalProb)
      }

      # Normalize last value, already smooth distribution
      p[, tail(eventIdx, 1), ] <- exp(p[, tail(eventIdx, 1), ])

      if (smoothProbsSt != "none") { # sample last Hidden State (HS)
        zSample[, tail(eventIdx, 1)] <- apply(
          p[, tail(eventIdx, 1), ], 1,
          \(x) sample.int(kStates, 1, prob = exp(x))
        )

        # Compute conditional log likelihood
        logLikCond[, tail(eventIdx, 1)] <- vapply(
          seq_len(nDraws),
          \(x) llEventState[tail(eventIdx, 1), x, zSample[x, tail(eventIdx, 1)]],
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
      logBeta <- array(0, dim = c(nDraws, kStates))

      for (tt in rev(head(eventIdx, -1))) {
        # \beta_j(t) * p(y_{t+1} | S_{t+1} = j)
        omegaBeta <- logBeta + llEventState[tt + 1, , ] # element-wise product

        # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
        # FFBS, a.k.a, Multi-move sampling
        if (smoothProbsSt == "joint") {
          probLastHS <- p[, tt, ] +
            t(vapply(
              seq_len(nDraws),
              \(x){
                lastHS <- zSample[x, tt + 1]
                # p(S_{t+1} = j|S_t = i) * p(y_{t+1}|S_{t+1} = j) * \beta_j(t)
                t(transProbs[idxTP[, lastHS], x, tt]) + omegaBeta[x, lastHS]
              },
              numeric(kStates)
            ))
          probLastHS <- exp(sweep(probLastHS, 1, rowLogSumExps(probLastHS)))
          zSample[, tt] <- apply(
            probLastHS, 1,
            \(x) sample.int(kStates, 1, prob = x)
          )

          # Compute conditional log likelihood
          logLikCond[, tt] <- vapply(
            seq_len(nDraws),
            \(x) llEventState[tt, x, zSample[x, tt]],
            numeric(1)
          )
        }
        for (k in seq_len(kStates)) {
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
        if (smoothProbsSt == "marginal") {
          zSample[, tt] <- apply(
            p[, tt, ], 1,
            \(x) sample.int(kStates, 1, prob = x)
          )
        }
      }

    }

    output[["smoothProbs"]] <- list(
      prob = p
    )

    if (smoothProbsSt != "none") {
      output[["smoothProbs"]][["zSample"]] <- zSample
    }

    output[["logLikRec"]] <- logLikRec
    output[["logLikCond"]] <- logLikCond
  }

  return(output)
}


#' Get Draws from the Posterior Distribution of a HMM-DyNAM
#'
#' @inheritParams HMMPostProcessing
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
#' callNetwork <- defineNetwork(nodes = actors, directed = TRUE) |>
#'   linkEvents(changeEvent = calls, nodes = actors)
#' callsDependent <- defineDependentEvents(
#'   events = calls, nodes = actors, defaultNetwork = callNetwork
#' )
#' data2stan <- CreateDataHMM(
#'   rateEffects = callsDependent ~ indeg + outdeg,
#'   choiceEffects = NULL
#' )
#'
#' drawsMCMC <- HMMDraws2LS(data2stan, cmdstanSamples)
#' }
HMMDraws2LS <- function(
    data2Stan,
    cmdstanSamples,
    kStates = data2Stan$dataStan$kS,
    type = c("2label.switching", "draws_matrix", "draws_array", "draws_df"),
    rescale = FALSE) {
  type <- match.arg(type)
  stopifnot(
    inherits(cmdstanSamples, "CmdStanMCMC"),
    inherits(data2Stan, "goldfish.latent.data"),
    is.numeric(kStates) && length(kStates) == 1 && kStates >= 2
  )

  if (rescale && type %in% c("2label.switching", "draws_array")) {
    stop("Not rescaling available for ", dQuote("type"), " = ",
      type,
      call. = FALSE
    )
  }

  model <- attr(data2Stan, "model")
  subModel <- attr(data2Stan, "subModel")

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
  draws <- cmdstanSamples$draws(
    parmsKeep,
    format = if (type == "draws_df") "draws_df" else "matrix"
  )

  namesTheta <- grep(paste0("^", theta), colnames(draws))
  idxTheta <- matrix(
    namesTheta,
    nrow = kStates, ncol = length(namesTheta) / kStates
  )
  idxEmission <- grep("^pi1", colnames(draws))

  output <- list(
    draws = draws,
    idxEmission = idxEmission,
    idxTheta = idxTheta,
    chains = cmdstanSamples$num_chains(),
    iterSampling = cmdstanSamples$metadata()[["iter_sampling"]]
  )

  if (subModel %in% c("both", "choice")) {
    idxBetaChoice <- matrix(
      grep("^betaChoice[[]", colnames(draws)),
      nrow = kStates, ncol = data2Stan$dataStan$Pchoice
    )
    output[["idxBetaChoice"]] <- idxBetaChoice
    kPchoice <- data2Stan$dataStan$Pchoice

    if (model == "DNCHMMRE") {
      idxGammaChoice <- array(
        grep("^gammaChoice[[]", colnames(draws)),
        dim = c(kStates, data2Stan$dataStan$G, data2Stan$dataStan$Qchoice)
      )
      output[["idxGammaChoice"]] <- idxGammaChoice

      idxLOmegaChoice <- array(
        grep("^L_OmegaChoice[[]", colnames(draws)),
        dim = c(kStates, data2Stan$dataStan$Qchoice, data2Stan$dataStan$Qchoice)
      )
      output[["idxLOmegaChoice"]] <- idxLOmegaChoice

      idxLSigmaChoice <- matrix(
        grep("^L_sigmaChoice[[]", colnames(draws)),
        nrow = kStates, ncol = data2Stan$dataStan$Qchoice
      )
      output[["idxLSigmaChoice"]] <- idxLSigmaChoice
    }
  } else {
    kPchoice <- 0
  }

  if (subModel %in% c("both", "rate")) {
    idxBetaRate <- matrix(
      if (subModel == "rate") {
        grep("^beta[[]", colnames(draws))
      } else {
        grep("^betaRate[[]", colnames(draws))
      },
      nrow = kStates, ncol = data2Stan$dataStan$Prate
    )
    output[["idxBetaRate"]] <- idxBetaRate
    kPrate <- data2Stan$dataStan$Prate
    
    if (model == "DNCHMMRE") {
      idxGammaRate <- array(
        grep("^gammaRate[[]", colnames(draws)),
        dim = c(kStates, data2Stan$dataStan$G, data2Stan$dataStan$Qrate)
      )
      output[["idxGammaRate"]] <- idxGammaRate

      idxLOmegaRate <- array(
        grep("^L_OmegaRate[[]", colnames(draws)),
        dim = c(kStates, data2Stan$dataStan$Qrate, data2Stan$dataStan$Qrate)
      )
      output[["idxLOmegaRate"]] <- idxLOmegaRate

      idxLSigmaRate <- matrix(
        grep("^L_sigmaRate[[]", colnames(draws)),
        row = kStates, col = data2Stan$dataStan$Qrate
      )
      output[["idxLSigmaRate"]] <- idxLSigmaRate
    }
  } else {
    kPrate <- 0
  }

  nDraws <- nrow(draws)

  # # rescale draws from betas if matrix of data frame and scale stats exists
  if (type %in% c("draws_matrix", "draws_df") &&
    !(rescale && !is.null(data2Stan[["scaleStats"]]))) {
    return(output)
  } else if (type %in% c("draws_matrix", "draws_df")) {
    output[["draws"]] <- as.data.frame(output[["draws"]])
    if (subModel %in% c("both", "choice")) {
      for (state in seq_len(kStates)) {
        output[["draws"]][, idxBetaChoice[state, ]] <- t(apply(
          output[["draws"]][, idxBetaChoice[state, ]],
          1,
          RescaleCoefs,
          scaleStats = data2Stan[["scaleStats"]],
          isRate = FALSE
        ))

        if (model == "DNCHMMRE") {
          Q <- data2Stan$dataStan$Qchoice
          stdDev <- data2Stan$scaleStats$choice$`scaled:scale`[data2Stan$dataStan$V1choice]
          for (gg in seq_len(data2Stan$dataStan$G))
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
        data2Stan$dataStan$Trate / 
        (data2Stan$dataStan$Nrate *
          mean(data2Stan$dataStan$timespan))
      ) - log(data2Stan$dataStan$rescalingFactor)
      for (state in seq_len(kStates)) {
        output[["draws"]][, idxBetaRate[state, ]] <- t(apply(
          output[["draws"]][, idxBetaRate[state, ]],
          1,
          RescaleCoefs,
          scaleStats = data2Stan[["scaleStats"]],
          offset = offSetInt,
          isRate = TRUE
        ))

        if (model == "DNCHMMRE") {
          Q <- data2Stan$dataStan$Qrate
          stdDev <- data2Stan$scaleStats$rate$`scaled:scale`[data2Stan$dataStan$V1rate]
          for (gg in seq_len(data2Stan$dataStan$G))
            output[["draws"]][, idxGammaRate[state, gg, ]] <- sweep(
              output[["draws"]][, idxGammaRate[state, gg, ]],
              2,
              stdDev,
              FUN = "/"
            )
        }
      }
      output[["draws"]][, idxTheta] <- output[["draws"]][, idxTheta] /
        data2Stan$dataStan$rescalingFactor
      
      
    }
    return(output)
  }

  if (type == "draws_array") {
    output[["draws"]] <- cmdstanSamples$draws(
      parmsKeep,
      format = "draws_array"
    )
    return(output)
  }

  nParms <- 1 + kPchoice + kPrate # + kStates # if theta included

  drawsReshape <- array(
    0,
    dim = c(nDraws, kStates, nParms),
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

  for (state in seq_len(kStates)) {
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

bindPPHMM <- function(
    listOutputs,
    type = c("both", "viterbi", "smoothProbs"),
    smoothProbsSt = c("joint", "marginal", "none")#,
    #model = c("DNHMM", "DNCHMM")
  ) {
  type <- match.arg(type)
  smoothProbsSt <- match.arg(smoothProbsSt)
  #model <- match.arg(model)

  output <- list()
  outputNames <- names(listOutputs[[1]])

  if (type %in% c("both", "viterbi")) {
    output[["viterbi"]] <- lapply(listOutputs, "[[", "viterbi") |>
      (\(x) Reduce(rbind, x = x))()
  }

  if (type %in% c("both", "smoothProbs")) {
    smoothProbs <- lapply(listOutputs, "[[", "smoothProbs") |>
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

    if (smoothProbsSt != "none") {
      output[["smoothProbs"]][["zSample"]] <- lapply(
        listOutputs, "[[", "smoothProbs"
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
      permuteMCMCArray(
        draws = draws, permutation = labelSwitching[["permutations"]][[x]],
        kStates = kStates, subModel = subModel
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
      (\(x) {
        Reduce(rbind.data.frame, x)
      })()

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
#' @param data2Stan a list output from `CreateDataHMM`
#' @param cmdstanSamples a list output from `CmdStanMCMC`
#' @param labelSwitching a list output from `label.switching`
#' @param postProcessing a list output from `HMMPostProcessing`
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

  drawsdf <- permuteMCMCArray(
    draws = draws, permutation = permutation,
    kStates = kStates, subModel = subModel
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
permuteMCMCArray <- function(
    draws, permutation, kStates,
    model = c("DNHMM", "DNCHMM", "DNCHMMRE"), subModel = c("choice", "rate", "both"),
    hmmType = c("discrete", "continuous")) {
  
  model <- match.arg(model)
  subModel <- match.arg(subModel)
  hmmType <- match.arg(hmmType)
  output <- as.data.frame(draws$draws)
  draws$draws <- as.data.frame(draws$draws)

  iterPermute <- apply(permutation, 1, \(x) any(x != seq_len(kStates))) |>
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
    if (subModel %in% c("both", "rate")) {
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
    if (hmmType == "discrete") {
      colPermutation <-
        draws$idxTheta[permutation[iter, ], permutation[iter, ]] |> as.vector()
    } else {
      colPermutation <- draws$idxTheta
      permIter <- permutation[iter, ]

      for (row in seq_along(permIter)) {
        for (col in seq_along(permIter)) {
          if (row == col) next
          jjPos <- ifelse(permIter[row] < permIter[col], permIter[col] - 1, permIter[col])
          jjNewPos <- ifelse(row < col, col - 1, col)
          colPermutation[row, jjNewPos] <- draws$idxTheta[permIter[row], jjPos]
        }
      }
      colPermutation <- colPermutation |> as.vector()
    }

    if (subModel %in% c("both", "choice")) {
      choicePermutation <- draws$idxBetaChoice[permutation[iter, ], ] |> as.vector()
      if (model == "DNCHMMRE") {
        choicePermutation <- c(
          choicePermutation,
          draws$idxGammaChoice[permutation[iter, ], , ] |> as.vector(),
          draws$idxLOmegaChoice[permutation[iter, ], , ] |> as.vector(),
          draws$idxLSigmaChoice[permutation[iter, ], ] |> as.vector()
        )
      }
    } else choicePermutation <- NULL
    if (subModel %in% c("both", "rate")) {
      ratePermutation <- draws$idxBetaRate[permutation[iter, ], ] |> as.vector()
      if (model == "DNCHMMRE") {
        ratePermutation <- c(
          ratePermutation,
          draws$idxGammaRate[permutation[iter, ], , ] |> as.vector(),
          draws$idxLOmegaRate[permutation[iter, ], , ] |> as.vector(),
          draws$idxLSigmaRate[permutation[iter, ], ] |> as.vector()
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
