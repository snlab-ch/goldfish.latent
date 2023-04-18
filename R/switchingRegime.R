#' Create data for Stan Switching Regime model
#'
#' The first step is create the data in the structure expected by the `Stan`
#' code designed for the DyNAM with Switching Regimes.
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
#' data2stan <- CreateDataSR(
#'   rateEffects = callsDependent ~ indeg + outdeg,
#'   choiceEffects = callsDependent ~ recip + trans
#' )
#' }
CreateDataSR <- function(
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
      inherits(supportConstraint, "formula")
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
      envir = envir
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
      envir = envir
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
  model = "DyNAMSR",
  subModel = subType
  ))
}

