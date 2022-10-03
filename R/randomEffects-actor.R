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
#' For the moment, it is only possible to add actors' covariates using the
#' effect `egoAlterInt()` and fixed the alter attribute to a variable with
#' constant value 1. Further enhancements would allow define it using
#' other monadic statistics like `indeg()`, `outdeg()`, `nodeTrans()` and using
#' the `ego()` effect to include actors' covariates.
#'
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
#' data2stan <- CreateDataModel(
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
  model = "DyNAM",
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
#' @param data2stan a `list` output of a [CreateData()] call.
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
#' #' \donttest{
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
CreateModelCode <- function(data2stan, ...) {
  stopifnot(inherits(data2stan, "goldfish.latent.data"))

  model <- attr(data2stan, "model")
  subModel <- attr(data2stan, "subModel")

  if (model == "DyNAM" && subModel == "choice") {
    if (data2stan[["dataStan"]][["Q"]] == 1) {
      stanCode <- readLines(
        system.file("stan", "MCM_RE1.stan", package = "goldfish.latent")
      )
    }
  } else stop("Not yet implemented for more than one random effect")

  if (requireNamespace("cmdstanr", quietly = TRUE) &&
      cmdstanr::cmdstan_version() >= "2.29.2") {
    model <- cmdstanr::write_stan_file(code = stanCode)
  } else
    stop(
      dQuote("cmdstanr"), " package and a working version of",
      dQuote("CmdStan"), "are required.",
      " Please follow Stan documentation for instructions on how to install.")

  return(model)
}
