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
    kRegimes = 3,
    supportConstraint = NULL,
    preprocessArgs = NULL,
    scale = FALSE,
    progress = getOption("progress"),
    envir = new.env()
) {
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
    is.numeric(kRegimes) && length(kRegimes) == 1 && kRegimes >= 2
  )

  # setting initial values of some arguments
  if (is.null(progress)) progress <- FALSE

  if (!is.null(rateEffects) && !is.null(choiceEffects)) {
    subType <- "both"


    if (as.character(choiceEffects[[2]]) != as.character(rateEffects[[2]]))
      stop("dependent event network needs to be the same for both formulas")
  } else if (!is.null(rateEffects)) {
    subType <- "rate"
  } else if (!is.null(choiceEffects)) {
    subType <- "choice"
  } else stop(
    dQuote("rateEffects"), " and ", dQuote("choiceEffects"),
    " arguments are NULL objects. Specify at least one model."
  )
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


    dataStanRate = within(dataProcessedRate, {
      Trate <- length(n_candidates)
      Nrate <- nrow(stat_all_events)
      Prate <- ncol(stat_all_events)
      startRate <- cumsum(c(1, head(n_candidates, -1)))
      endRate <- cumsum(n_candidates)
      namesEffects <- gsub("\\$", "Of", namesEffects)
      Xrate <- setNames(stat_all_events, namesEffects)
      if (hasIntercept) {
        choseRate <- selected[, 1] + (startRate - 1) * isDependent
        offsetInt <- log(Trate / (sum(timespan) * mean(n_candidates)))
      } else choseRate <- selected[, 1] + (startRate - 1)
      rm(stat_all_events, n_candidates, sender,
         namesEffects, effectDescription, selected)
    })

    if (scale) {
      Xrate <- scale(
        if (dataStanRate$hasIntercept)
          dataStanRate$Xrate[, -1] else dataStanRate$Xrate
      )
      scaleStats[["rate"]] <- attributes(Xrate)

      if (dataStanRate$hasIntercept)
        dataStanRate$Xrate[, -1] <- Xrate[, ] else
          dataStanRate$Xrate <- Xrate[, ]
    }
    # table(idxEvents[1, ] == dataStanRate$startRate)
    # table(idxEvents[2, ] == dataStanRate$endRate)
    # chose <- which(expandedDF$selected)
    # table(chose == dataStanRate$choseRate[dataStanRate$isDependent])
  } else dataProcessedRate <- dataStanRate <- NULL

  if (!is.null(choiceEffects)) {
    feTerms <- terms(choiceEffects)
    cstrTerms <- if (!is.null(supportConstraint))
      terms(supportConstraint) else NULL

    if (length(attr(cstrTerms, "term.labels")) > 1)
      stop(dQuote("supportConstraint"), " argument only works for one effect.")

    if (any(attr(feTerms, "order") != 1))
      stop(dQuote("choiceEffects"),
           "formula argument doesn't support interactions yet")

    termsDyNAM <- c(attr(feTerms, "term.labels"),
                    attr(cstrTerms, "term.labels")) |> unique()

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
    } else effectDescription <- dataProcessedChoice$effectDescription

    # create objects for Stan
    nTotal <- nrow(expandedDF)

    idxEvents <- tapply(seq_len(nTotal), expandedDF$event, range) |>
      simplify2array()

    Xmat <- as.matrix(expandedDF[, namesEffects])

    dataStanChoice = list(
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
      kR = kRegimes,
      alpha = c(floor((kRegimes - 1) / (1 - 0.8)), 1)
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

  return(structure(list(
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
#' @param kRegimes an integer with the number of states used to get samples
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
    data2Stan, cmdstanSamples,
    kRegimes = data2Stan$dataStan$kR,
    type = c("both", "viterbi", "smoothProbs"),
    smoothProbsSt = c("joint", "marginal", "none")
    ) {
  type <- match.arg(type)
  smoothProbsSt <- match.arg(smoothProbsSt)

  stopifnot(
    inherits(cmdstanSamples, "CmdStanMCMC"),
    inherits(data2Stan, "goldfish.latent.data"),
    is.numeric(kRegimes) && length(kRegimes) == 1 && kRegimes >= 2
  )
  # init output
  output <- list()

  model <- attr(data2Stan, "model")
  subModel <- attr(data2Stan, "subModel")

  # extract draws and reformat for posterior computations
  drawsObject <- HMMDraws2LS(
    data2Stan, cmdstanSamples, kRegimes = kRegimes, type = "matrix"
  )
  draws <- drawsObject$draws
  idxTheta <- drawsObject$idxTheta
  idxEmission <- drawsObject$idxEmission
  idxBetaChoice <- drawsObject$idxBetaChoice
  idxBetaRate <- drawsObject$idxBetaRate

  nDraws <- nrow(draws)

  data2Stan <- data2Stan$dataStan
  # define sizes
  nEvents <- ifelse(subModel %in% c("both", "rate"), "Trate", "Tchoice")
  nEvents <- data2Stan[[nEvents]]

  isRes <- !is.null(data2Stan$Nres)
  TT <- ifelse(isRes, data2Stan$Nres, nEvents)


  if (subModel %in% c("choice", "both")) {
    llEventState <- array(0, dim = c(nEvents, nDraws, kRegimes))

    for (kR in seq_len(kRegimes)) {
      xb <- tcrossprod(data2Stan$Xchoice, draws[, idxBetaChoice[kR, ]])

      for (event in seq_len(nEvents))
        llEventState[event, , kR] <-
          xb[data2Stan$choseChoice[event], ] -
          colLogSumExps(
            xb,
            rows = seq(data2Stan$startChoice[event],
                       data2Stan$endChoice[event])
          )
    }

    if (isRes)
      llEventState <- array(apply(
        llEventState,
        3,
        \(x) lapply(
          seq_len(TT),
          \(y) colSums2(
            x,
            rows = seq(data2Stan$resA[y], data2Stan$resA[y + 1] - 1)
          )
        ) |> Reduce(rbind, x = _)
      ), dim = c(TT, nDraws, kRegimes))


    if (type %in% c("both", "viterbi")) {
      # back-pointer to the most likely previous state on the most probable path
      bpointer <- array(0, dim = c(TT, kRegimes, nDraws))
      # max prob for the sequence up to t that ends with an emission from state k
      delta <- array(0, dim = c(TT, kRegimes, nDraws))

      # forward past computing most likely state from previous state
      # first observation: p(y_1| z_1) * p(z_1) (emission prob)
      delta[1, , ] <- log(draws[, idxEmission]) + llEventState[1, , ]

      for (tt in seq(2, TT)) {
        logp <- array(sapply(
          seq_len(kRegimes),
          \(x) {
            T_1_xj = t(delta[tt - 1, , ]) + draws[, idxTheta[, x]] +
              llEventState[tt, , ]
            T_2_xj = apply(T_1_xj, 1, which.max)
            return(rbind(T_1_xj = apply(T_1_xj, 1, max), T_2_xj))
          }
        ), dim = c(2, kRegimes, nDraws))

        delta[tt, , ] <- logp[1, , ]
        bpointer[tt, , ] <- logp[2, , ]
      }

      # backward past
      z <- array(0L, dim = c(nDraws, TT))

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
      p <- array(0, dim = c(nDraws, TT, kRegimes))

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
        for (k in seq_len(kRegimes))
          p[, tt, k] <- rowLogSumExps(
            log(draws[, idxTheta[, k]]) + p[, tt - 1, ]
          ) +
          # # filter for S_t: p(S_t = k| y_t) =
          # #  p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) /
          # #  \sum_k p(y_t|S_t=k,y_{t-1}) p(S_t=k|y_{t-1}) ; unnormalize enough
            llEventState[tt, , k]
      }

      # Normalize last value, already smooth distribution
      p[, TT, ] <- exp(sweep(p[, TT, ], 1, rowLogSumExps(p[, TT, ])))

      if (smoothProbsSt != "none") # sample last Hidden State (HS)
        zSample[, TT] <- apply(
          p[, TT, ], 1,
          \(x) sample.int(kRegimes, 1, prob = x)
        )

      # backward: smoother suggested in Hamilton expresses these as marginal
      # probabilities from the joint distribution of S_t and S_T | y
      # Implementation follows Stan hmm_hidden_state_prob()
      #

      # initial ending state ass as given (uniform)
      logBeta <- array(0, dim = c(nDraws, kRegimes))

      for (tt in seq(TT - 1, 1)) {
        # # Baum-Welch alg
        # #
        omegaBeta <- logBeta + llEventState[tt + 1, , ] # element-wise product

        # intermezzo: sample the tt HS conditional on (tt+1)st HS as in Stan
        if (smoothProbsSt == "joint"){
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
            \(x) sample.int(kRegimes, 1, prob = x)
          )
        }
        for (k in seq_len(kRegimes))
          logBeta[, k] <- rowLogSumExps(
            log(draws[, idxTheta[k, ]]) + omegaBeta
          )

        # running normalization
        logBeta <- sweep(logBeta, 1, apply(logBeta, 1, max))

        #
        gammat <- logBeta + p[, tt, ]
        p[, tt, ] <- exp(sweep(gammat, 1, rowLogSumExps(gammat)))

        # sample tt HS from marginal distribution
        if (smoothProbsSt == "marginal")
          zSample[, tt] <- apply(
            p[, tt, ], 1,
            \(x) sample.int(kRegimes, 1, prob = x)
          )
      }

      output[["smoothProbs"]] <- list(
        prob = p
      )

      if (smoothProbsSt != "none")
        output[["smoothProbs"]][["zSample"]] <- zSample
    }

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
    kRegimes = data2Stan$dataStan$kR,
    type = c("2label.switching", "matrix")
) {
  type = match.arg(type)
  stopifnot(
    inherits(cmdstanSamples, "CmdStanMCMC"),
    inherits(data2Stan, "goldfish.latent.data"),
    is.numeric(kRegimes) && length(kRegimes) == 1 && kRegimes >= 2
  )

  model <- attr(data2Stan, "model")
  subModel <- attr(data2Stan, "subModel")

  # extract draws and reformat for posterior computations
  parmsKeep <- c(
    "lp__", "theta", "pi1",
    if (subModel %in% c("both", "choice")) "betaChoice",
    if (subModel %in% c("both", "rate")) "beta"
  )
  draws <- cmdstanSamples$draws(parmsKeep, format = "matrix")

  idxTheta <- matrix(
    grep("^theta", colnames(draws)),
    nrow = kRegimes, ncol = kRegimes
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
      nrow = kRegimes, ncol = data2Stan$dataStan$Pchoice
    )
    output[["idxBetaChoice"]] <- idxBetaChoice
    kPchoice <- data2Stan$dataStan$Pchoice
  } else kPchoice <- 0

  if (subModel %in% c("both", "rate")) {
    idxBetaRate <- matrix(
      grep("^beta[[]", colnames(draws)),
      nrow = kRegimes, ncol = data2Stan$Pchoice
    )
    output[["idxBetaRate"]] <- idxBetaRate
    kPrate <- data2Stan$dataStan$Prate
  } else kPrate <- 0

  nDraws <- nrow(draws)


  if (type == "matrix") return(output)

  nParms <- 1 + kPchoice + kPrate # + kRegimes # if theta included

  drawsReshape <- array(
    0,
    dim = c(nDraws, kRegimes, nParms),
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

  for (state in seq_len(kRegimes)) {
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
