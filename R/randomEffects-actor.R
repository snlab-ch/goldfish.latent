#' Create data for Stan
#'
#' The first step is create the data in the structure expected by the `Stan`
#' code designed for the DyNAM with actor random effects.
#' Additional information of the effects used during preprocessing is added to
#' the return object.
#'
#' The model formulation is composed of two parts: the fixed effect part
#' `fixedEffects` and the random effect part `randomEffects`.
#' The fixed effect formulation works as the formulation of models in `goldfish`.
#' All the effects used here in the right hand side are considered to be fixed,
#' hence, they won't have random effects.
#' The random effects formulation considers the possibility to have more than
#' one random effect, and that every random effect might be explain by actors'
#' monadic statistics or covariates.
#' A formula like `effect(network) ~ egoAlterInt(list(attrEgo, attrAlter))`
#' indicates that the `effect(network)` is added to the model having
#' random effects and those could be explain by the effects included on the
#' right hand side.
#'
#' @param randomEffects a `list`, each component is a `formula` and represents a
#' random effect to include in the model. Each `formula` has on the left hand
#' side the effect specification that plays the role of random effect, and
#' on the right hand side effects that would explain the variability of that
#' random effect.
#' @param fixedEffects a `formula` specification as in [goldfish::estimate()].
#' The effects include in the right hand side play the role of fixed effects in
#' the model.
#' @param model Current version only support `"DyNAM"` model, enhancements
#' on the code would allow to use `"REM"` model too.
#' @param subModel Current version only support `"choice"` sub-model.
#' @param supportConstraint a `formula` with only an effect that gives the
#' information of the alters available to received an event at each moment
#' of time. The left hand side of the formula is left empty.
#' It usually will correspond to a call of the `tie(network)` where `network`
#' is a binary adjacency matrix where the value 1 indicated that actor $j$ is
#' available to receive an event from actor $i$.
#' In the case that some actors left or join the process at any point of time is
#' better to use the `present` variable in the node data frame linking to it the
#' time varying changes of composition of the actors set.
#' @param preprocessArgs a list with additional preprocess arguments used by
#' [goldfish::GatherPreprocessing()] to compute the changes statistics of the
#' event sequence.
#' @param progress logical argument passed to [goldfish::GatherPreprocessing()]
#' to show the progress of the preprocessing of the events sequence.
#'
#' @return a list with the following components.
#' \describe{
#'   \item{dataStan}{a list with the information necessary to run a HMC using
#'   Stan.}
#'   \item{sendersIx}{a data frame with the label of the sender from the
#'   nodes data frame, and the index assign to the random coefficient.}
#'   \item{namesEffects}{a character vector with terms in the random and fixed
#'   effects formulas and their final name.}
#'   \item{effectDescription}{an array with detailed and comprehensible
#'   information of the terms used in the random and fixed effects formulas.}
#' }
#' @export
#' @importFrom stats terms setNames as.formula model.matrix reformulate
#' @importFrom goldfish GatherPreprocessing
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
#' data2stan <- CreateData(
#'   randomEffects = list(inertia ~ 1),
#'   fixedEffects = callsDependent ~ recip + trans
#' )
#' }
CreateData <- function(
  randomEffects,
  fixedEffects,
  model = c("DyNAM", "REM"),
  subModel = c("choice", "rate", "choice_coordination"),
  supportConstraint = NULL,
  preprocessArgs = NULL,
  progress = getOption("progress")
) {
  ### 0. check parameters----
  model <- match.arg(model)
  subModel <- match.arg(subModel)

  stopifnot(
    is.null(progress) || inherits(progress, "logical"),
    is.null(preprocessArgs) ||
      inherits(preprocessArgs, "list"),
    inherits(randomEffects, "list"),
    inherits(fixedEffects, "formula"),
    is.null(supportConstraint) ||
      inherits(supportConstraint, "formula")
  )

  if (!(model == "DyNAM" && subModel == "choice"))
    stop("Model and subModel combination not available yet, WIP")

  if (length(randomEffects) != 1)
    stop("The current implementation only consider one random effect.")

  # setting initial values of some arguments
  if (is.null(progress)) progress <- FALSE

  envir <- new.env()

  # formula treatment
  reTerms <- lapply(randomEffects, terms)
  feTerms <- terms(fixedEffects)
  cstrTerms <- if (!is.null(supportConstraint))
    terms(supportConstraint) else NULL

  if (length(attr(cstrTerms, "term.labels")) > 1)
    stop(dQuote("supportConstraint"), " argument only works for one effect.")

  if (any(attr(feTerms, "order") != 1))
    stop(dQuote("fixedEffects"),
         "formula argument doesn't support interactions yet")

  if (any(vapply(reTerms, \(x) any(attr(x, "order") != 1), logical(1))))
    stop(dQuote("fixedEffects"),
         "formula argument doesn't support interactions yet")

  # modify effects used to explain random effects to ego versions
  reTerms <- lapply(
    randomEffects, ModifyFormulaRE,
    fixedEffects = fixedEffects, envir = envir) |>
    lapply(terms)

  xDyNAM <- Reduce(
    c,
    lapply(reTerms, \(x) c(attr(x, "term.labels"), as.character(x[[2]])))
  )
  reDyNAM <- vapply(reTerms, \(x) as.character(x[[2]]), character(1))

  termsDyNAM <- c(xDyNAM, attr(feTerms, "term.labels"),
                  attr(cstrTerms, "term.labels")) |> unique()

  formulaDyNAM <- reformulate(
    termsDyNAM,
    response = as.character(fixedEffects[[2]])
  )

  # create a full matrix for filtering
    # get names from formula
  # parsedFormula <- goldfish:::parseFormula(formula = formulaDyNAM)
  # objectsEffectsLink <-
  #   goldfish:::getObjectsEffectsLink(rhsNames = parsedFormula$rhsNames)
  # effectDescription <- goldfish:::GetDetailPrint(
  #   objectsEffectsLink = objectsEffectsLink,
  #   parsedformula = parsedFormula
  # )
  # namesEffects <- goldfish:::CreateNames(
  #   effectDescription, sep = "_", joiner = "_"
  # )
    # process data
  dataProcessed <- GatherPreprocessing(
    formula = formulaDyNAM,
    model = if (model == "DyNAM" && subModel == "choice") "DyNAMRE" else model,
    subModel = subModel,
    preprocessArgs = preprocessArgs,
    progress = progress,
    envir = envir
  )

  nEvents <- length(dataProcessed$sender)

  namesEffects <- setNames(
    gsub("\\$", "Of", dataProcessed$namesEffects),
    termsDyNAM
  )

  expandedDF <- cbind(
    setNames(
      as.data.frame(dataProcessed$stat_all_events),
      namesEffects
    ),
    data.frame(
      event = rep(seq.int(nEvents), dataProcessed$n_candidates),
      selected = sequence(dataProcessed$n_candidates) ==
        rep(dataProcessed$selected, dataProcessed$n_candidates),
      sender = rep(dataProcessed$sender, dataProcessed$n_candidates)
    )
  )

  # subset if constraint
  if (!is.null(supportConstraint)) {
    cstrName <- namesEffects[attr(cstrTerms, "term.labels")]
    keep <- expandedDF[, cstrName] == 1
    expandedDF <- expandedDF[keep, !names(expandedDF) %in% cstrName]
    effectDescription <-
      dataProcessed$effectDescription[!namesEffects %in% cstrName, ]
    namesEffects <- namesEffects[!namesEffects %in% cstrName]
  } else effectDescription <- dataProcessed$effectDescription

  # create objects for Stan
  nTotal <- nrow(expandedDF)

  seqExDF <- seq.int(nTotal)

  idxEvents <- tapply(seq.int(nTotal), expandedDF$event, range) |>
    simplify2array()

  sendersIx <- data.frame(label = sort(unique(expandedDF$sender))) |>
    within(index <- seq.int(label))

  expandedDF[, "senderIx"] <-
    sendersIx[match(expandedDF[, "sender"], sendersIx[, "label"]), "index"]

  formulaDyNAMRE <- lapply(
    reTerms,
    \(x) paste(
      if (length(attr(x, "term.labels")) > 0)
        paste(
          namesEffects[match(as.character(x[[2]]), termsDyNAM)], "/",
          namesEffects[match(attr(x, "term.labels"), termsDyNAM)]
        ) else namesEffects[match(as.character(x[[2]]), termsDyNAM)],
      collapse = " + "
    )
  ) |>
    c(namesEffects[match(attr(feTerms, "term.labels"), termsDyNAM)]) |>
    paste(collapse = " + ")

  Xmat <- model.matrix(
    as.formula(paste("~ ", formulaDyNAMRE, " + 0")),
    data = expandedDF
  )

  return(structure(list(
    dataStan = list(
      T = nEvents,
      N = nTotal,
      P = ncol(Xmat),
      Q = length(randomEffects),
      A = nrow(sendersIx),
      start = idxEvents[1, ],
      end = idxEvents[2, ],
      sender = expandedDF[, "senderIx"],
      X = Xmat,
      Z = expandedDF[, namesEffects[match(reDyNAM, termsDyNAM)]],
      chose = which(expandedDF[, "selected"])
    ),
    sendersIx = sendersIx,
    namesEffects = namesEffects,
    effectDescription = effectDescription
  ),
  class = "goldfish.latent.data",
  model = "DyNAMRE",
  subModel = "choice"
  ))
}

ModifyFormulaRE <- function(reFormula, fixedEffects, envir = new.env()) {
  stopifnot(inherits(reFormula, "formula"))

  reForTerms <- terms(reFormula)
  effectsFormula <- attr(reForTerms, "term.labels")
  if (length(effectsFormula) == 0) {
    return(reFormula)
  }

  depName <- goldfish:::getDependentName(fixedEffects)
  defaultNetworkName <- attr(get(depName, envir = envir), "defaultNetwork")

  # modify calls
  effectForMod <- vapply(
    effectsFormula,
    \(x) {
      effectLang <- str2lang(x)
      if (length(effectLang) == 1) {
        return(deparse(as.call(
          list(effectLang, as.symbol(defaultNetworkName), type = "ego")
        )))
      } else if (deparse(effectLang[[1]]) != "ego") {
        return(deparse(as.call(
          c(list(effectLang[[1]]), as.list(effectLang[-1]), list(type = "ego"))
        )))
      } else return(x)
    },
    character(1)
  )

  # formula after modifications
  return(reformulate(effectForMod, response = reFormula[[2]]))
}

#' Compute log-likelihood using MCMC samples
#'
#' The function computes the log-likelihood for each MCMC sample from the
#' posterior distribution.
#'
#' Argument `type` allows to compute the conditional or marginal version of the
#' log-likelihood. The marginal version uses a Gauss-Hermite quadrature
#' approximation to integrate out the random effects
#' \insertCite{Merkle2019}{goldfish.latent}. The code is an adaptation from
#' the supplementary material of \insertCite{Merkle2019;textual}{goldfish.latent}.
#'
#' @param cmdstanrSamples a `draws_array` or a `CmdStanFit` object with MCMC
#'   samples from the posterior distribution. In the case of a `draws_array`
#'   object is expected to have three dimensions corresponding to iterations,
#'   chain and variables.
#' @param dataStan a `list` output of a [CreateData()] call.
#' @param type a `character` value. It indicates whether the log-likelihood
#'   computation should return the `"conditional"` or the `"marginal"` version.
#' @param nNodes an `integer`. The number of quadrature point to use in the
#'   Gauss-Hermite quadrature approximation use to integrate out the
#'   random-effects.
#' @param splitSize an `integer` or `NULL`. It is use when the `type` is
#'   `"conditional"`. When it is `NULL`, the `splitSize` is set to have
#'   roughly `4e4` rows sent to a processor when
#'   the matrix `X`, containing the change statistics, has more
#'   than `1e6`rows. If `X` has less than `1e6` rows, the `splitSize` is set
#'   in such way that every processor would have the same amount of rows to
#'   process.
#'   The default value is `NULL`.
#' @param spec A specification appropriate to the type of cluster, see
#'   [parallel::makeCluster()] for a detail description. In the simplest case,
#'   an integer defining the number of processors to use during the parallel
#'   computation.
#' @param ... Additional arguments and options to be passed to the
#'   [parallel::makeCluster()] call.
#'
#' @return An array with the log-likelihood for each event when
#'   `type = "conditional"` or the log-likelihood for each sender actor after
#'   the random effects are integrated out when `type = "marginal"`.
#'
#' @references
#' \insertRef{Merkle2019}{goldfish.latent}
#'
#' @export
#'
#' @importFrom Rdpack reprompt
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
#'
#' mod01 <- cmdstan_model(stanCode)
#' mod01Samples <- mod01$sample(
#'   data = data2stan[["dataStan"]],
#'   parallel_chains = 4, chains = 4,  iter_warmup = 500, iter_sampling = 500,
#'   show_messages = FALSE
#' )
#'
#' margLogLikMod01 <- ComputeLogLikelihood(mod01Samples, data2stan, spec = 4)
#' condLogLikMod01 <- ComputeLogLikelihood(mod01Samples, data2stan,
#'                                         type = "conditional", spec = 4)
#' }
ComputeLogLikelihood <- function(
  cmdstanrSamples,
  dataStan,
  type = c("marginal", "conditional"),
  nNodes = ifelse(type == "marginal", 11L, NULL),
  splitSize = NULL,
  spec = parallel::detectCores() - 1,
  ...
) {
  stopifnot(
    inherits(cmdstanrSamples, c("CmdStanFit", "draws")),
    inherits(dataStan, "goldfish.latent.data"),
    is.null(splitSize) || inherits(splitSize, "numeric") &&
      length(splitSize) == 1
  )

  type <- match.arg(type)


  if (dataStan[["dataStan"]][["Q"]] > 1)
    stop("Likelihood computation for a model with more than one random-effect",
         " is not yet available.")

  if (inherits(cmdstanrSamples, "CmdStanFit")) {
    draws <- cmdstanrSamples$draws("gamma_raw")
    drawsDimnames <- dimnames(draws)
    draws <- list(
      beta = cmdstanrSamples$draws("beta"),
      sigma = cmdstanrSamples$draws("sigma"),
      gamma_raw = draws
    ) |>
      lapply(\(x) apply(x, 3, rbind))
  } else if (length(dim(cmdstanrSamples)) == 3) {
    drawsDimnames <- dimnames(cmdstanrSamples)

    variableNames <- dimnames(cmdstanrSamples)[[3]]

    draws <- list(
      beta = cmdstanrSamples[, , grepl("^beta", variableNames)],
      sigma = cmdstanrSamples[, , grepl("^sigma$", variableNames)],
      gamma_raw = cmdstanrSamples[, , grepl("^gamma_raw", variableNames)]
    ) |>
      lapply(\(x) apply(x, 3, rbind))
  } else
    stop(
      dQuote("cmdstanrSamples"),
      " argument expects a three dimensional ", dQuote("draws"), " object."
    )

  if (is.null(splitSize) & type == "conditional") {
    if (!is.numeric(spec) || length(spec) != 1)
      stop(
        "Please provide an integer number for", dQuote("splitSize"),
        " parameter. It's not possible to assign a value with current value of",
        dQuote("spec")
      )

    splitSize <- if (spec == 1) NULL else
      ifelse(
        dataStan[["dataStan"]][["N"]] > 1e6,
        4e4 / dataStan[["dataStan"]][["A"]],
        dataStan[["dataStan"]][["T"]] / spec
      ) |> floor()

    eventsPerCore <- if (!is.null(splitSize)) {
      parallel::splitIndices(
        dataStan[["dataStan"]][["T"]],
        floor(dataStan[["dataStan"]][["T"]] / splitSize)
      )  |>
        lapply(range)
    } else NULL
  }

  # create cluster and initialize workers
  if (length(spec) == 1 & spec == 1) {
    cl <- NULL
  } else {
    cl <- parallel::makeCluster(spec = spec, ...)
    ignore <- parallel::clusterEvalQ(cl, {library(matrixStats);NULL})
    on.exit(parallel::stopCluster(cl))
  }
  #
  if (type == "conditional") {
    if (!is.null(eventsPerCore)) {
      logLik <- parallel::clusterApplyLB(
        cl = cl,
        seq.int(length(eventsPerCore)),
        fun = LogLikCondRE,
        draws = draws,
        dataList = dataStan[["dataStan"]],
        eventsPerCore = eventsPerCore
      )
      logLik <- Reduce(f = cbind, x = logLik)
    } else
      logLik <- LogLikCondRE(
        eventsIter = NULL,
        draws = draws,
        dataList = dataStan[["dataStan"]],
        eventsPerCore = NULL
      )

    drawsDimnames$variable <- sprintf("event[%d]", seq.int(ncol(logLik)))

  } else if (type == "marginal") {
    logLik <- mllDyNAMChoice(
      draws = draws,
      dataList = dataStan[["dataStan"]],
      nNodes = nNodes,
      cl = cl
    )

    drawsDimnames$variable <- sprintf("actor[%d]", seq.int(ncol(logLik)))
  }

  return(structure(
    logLik, class = c("draws_array", "draws", "array"),
    dim = sapply(drawsDimnames, length), dimnames = drawsDimnames
  ))
}

#' marginal likelihoods for the DyNAM choice
#'
#' Function to obtain marginal likelihoods with parallel processing.
#'
#' @param draws Data list. Draws from a Fitted Stan model converted to a list
#'   with three components: beta, gamma_raw and sigma in the case of a single
#'   random effect.
#' @param dataList Data list used in fitting the model
#' @param nNodes Number of adaptive quadrature nodes to use
#' @param cl  cluster defined for parallel computing
#'
#' @return A two dimensional array. Every row corresponds to posterior sample
#' iterations.Every column corresponds to a sending actor. Values are the
#' marginal log-likelihood after integrating out random effects.
#' @noRd
#' @importFrom matrixStats colLogSumExps rowLogSumExps
#'
#' @examples mllDyNAMChoice(draws, data2stan, 11)
mllDyNAMChoice <- function(draws, dataList, nNodes, cl = NULL) {

  # Get standard quadrature points
  quad <- statmod::gauss.quad.prob(nNodes, "normal", mu = 0, sigma = 1)
  # logarithm of adapted weights
  quad$logWA <- log(quad$weights) + log(2 * pi) / 2 + quad$nodes^2 / 2

  # draws <- extract(stan_fit, stan_fit@model_pars)
  # post_means <- better_posterior_means(draws)

  # Separate out draws for means and SD from MCMC samples
  gamma <- sweep(draws$gamma_raw, 1, draws$sigma, FUN = "*")
  mcmcStat <- list(
    means = matrixStats::colMeans2(gamma),
    sd = matrixStats::colSds(gamma)
  )

  nDraws <- nrow(gamma)

  # add helper data
  dataList$senderEvent <- dataList$sender[dataList$start]
  # Function to compute the approximate marginal log-lik for sender
  fMarginal <- function(sender, draws, mcmcStat, dataList, quad, nNodes) {
    events <- which(dataList$senderEvent == sender)

    #
    mcmcSd <- mcmcStat$sd[sender]
    adaptNodes <- mcmcStat$means[sender] + mcmcSd * quad$nodes
    #
    nEvents <- length(events)

    mll <- array(0, dim = c(nDraws, nNodes))
    # mll <- list()
    # contador <- 1
    for (event in events) {
      keep <- seq.int(dataList$start[event], dataList$end[event])
      xb <- tcrossprod(dataList$X[keep, ], draws$beta)
      choice <- which(dataList$chose[event] == keep)

      Z <- dataList$Z[keep]

      mll <- mll +
        sapply(
          seq.int(nNodes),
          function(i) {
            utility <- sweep(xb, 1, Z * adaptNodes[i], FUN = "+")
            # # the log of the prob is utility - logSumExp of the utility choice set
            utility[choice, ] - colLogSumExps(utility)

          }
        )
    }
    # # l_c + log(prob prior)
    mll <- mll + outer(draws$sigma[, 1], adaptNodes, function(x, y) dnorm(y, sd = x, log = TRUE))
    # dnorm(adaptNodes[i], sd = draws$sigma, log = TRUE)

    # log(\prod \sum_{qdr points} lik (qdr point)) = \sum logSumExp( log(log_lik (qdr point)))
    # l_c + log(prob prior) + log(adapted weight)
    rowLogSumExps(sweep(mll, 2, quad$logWA + log(mcmcSd), FUN = "+"))
  }

  # Parallel by sender
  if (!is.null(cl)) {
  parallel::parSapplyLB(
    cl,
    seq.int(dataList$A),
    fMarginal,
    draws = draws,
    mcmcStat = mcmcStat,
    dataList = dataList,
    quad = quad,
    nNodes = nNodes
  )
  } else
    sapply(
      seq.int(dataList$A),
      fMarginal,
      draws = draws,
      mcmcStat = mcmcStat,
      dataList = dataList,
      quad = quad,
      nNodes = nNodes
    )
}

# # function to compute the likelihood for mean draw for model w.o RE
LogLikCondWORE <- function(eventsIter, draws, dataList, eventsPerCore = NULL) {
  if (!is.null(eventsIter)) {
    events <- eventsPerCore[[eventsIter]]
    startI <- dataList$start[head(events, 1)]
    endI <- dataList$end[tail(events, 1)]
    X <- dataList$X[seq.int(startI, endI), ]

    seqEventsKeep <- seq.int(head(events, 1), tail(events, 1))
    start <- dataList$start[seqEventsKeep] - (startI - 1)
    end <- dataList$end[seqEventsKeep] - (startI - 1)
    chose <- dataList$chose[seqEventsKeep] - (startI - 1)

    nE <- length(seqEventsKeep)
  } else {
    X <- dataList$X
    start <- dataList$start
    end <- dataList$end
    chose <- dataList$chose
    nE <- dataList$T
  }

  xb <- tcrossprod(X, draws)

  ll <- array(0, dim = c(nrow(draws), nE))
  # ,
  #             dimnames = list(iteration = seq.int(nrow(draws)),
  #                             variable = sprintf('', seq())))

  for (event in seq.int(nE)) {
    ll[, event] <- xb[chose[event], ] -
      colLogSumExps(xb, rows = seq(start[event], end[event]))
  }

  return(ll)
}

LogLikCondRE <- function(eventsIter, draws, dataList, eventsPerCore = NULL) {
  if (!is.null(eventsIter)) {
    events <- eventsPerCore[[eventsIter]]
    startI <- dataList$start[head(events, 1)]
    endI <- dataList$end[tail(events, 1)]
    seqDataKeep <- seq.int(startI, endI)
    X <- dataList$X[seqDataKeep, ]
    Z <- dataList$Z[seqDataKeep]
    sender <- dataList$sender[seqDataKeep]

    seqEventsKeep <- seq.int(head(events, 1), tail(events, 1))
    start <- dataList$start[seqEventsKeep] - (startI - 1)
    end <- dataList$end[seqEventsKeep] - (startI - 1)
    chose <- dataList$chose[seqEventsKeep] - (startI - 1)

    nE <- length(seqEventsKeep)
  } else {
    X <- dataList$X
    Z <- dataList$Z
    sender <- dataList$sender
    start <- dataList$start
    end <- dataList$end
    chose <- dataList$chose
    nE <- dataList$T
  }

  gamma <- sweep(draws$gamma_raw, 1, draws$sigma, FUN = "*")

  xb <- tcrossprod(X, draws$beta) +  t(sweep(gamma[, sender], 2, Z, FUN = "*"))

  ll <- array(0, dim = c(nrow(gamma), nE))

  for (event in seq(nE)) {
    ll[, event] <- xb[chose[event], ] -
      colLogSumExps(xb, rows = seq(start[event], end[event]))
  }

  return(ll)
}

