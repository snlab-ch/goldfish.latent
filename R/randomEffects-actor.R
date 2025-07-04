#' Create data for Stan
#'
#' The first step is create the data in the structure expected by the `Stan`
#' code designed for the DyNAM with actor random effects.
#' Additional information of the effects used during preprocessing is added to
#' the return object.
#'
#' The model formulation is composed of two parts: the fixed effect part
#' `fixed_effects` and the random effect part `random_effects`.
#' The fixed effect formulation works as the formulation of models in
#' [goldfish::estimate_dynam()].
#' All the effects used here in the right hand side are considered to be fixed,
#' hence, they won't have random effects.
#' The random effects formulation considers the possibility to have more than
#' one random effect, and that every random effect might be explain by actors'
#' monadic statistics or covariates.
#' A formula like
#' `effect(network) ~ ego_alter_interaction(list(attr_ego, attr_alter))`
#' indicates that the `effect(network)` is added to the model having
#' random effects and those could be explain by the effects included on the
#' right hand side.
#'
#' @param random_effects a `list`, each component is a `formula` and
#' represents a random effect to include in the model. Each `formula` has on the
#' left hand side the effect specification that plays the role of random effect,
#' and on the right hand side effects that would explain the variability of
#' that random effect.
#' @param fixed_effects a `formula` specification as in
#' [goldfish::estimate_dynam()].
#' The effects include in the right hand side play the role of fixed effects in
#' the model.
#' @param model Current version only support `"DyNAM"` model, enhancements
#' on the code would allow to use `"REM"` model too.
#' @param sub_model Current version only support `"choice"` sub-model.
#' @param support_constraint a `formula` with only an effect that gives the
#' information of the restricted set to consider.
#' In the case of the `"choice"` sub-model, it corresponds to the choice set
#' available to received an event at each moment of time.
#' In the case of the `"rate"` sub-model, it corresponds to the competing set
#' available to send an event at each moment of time.
#' In the case of the `"REM"` model, it corresponds to the set of dyads
#' available to create an event at each moment of time.
#' The left hand side of the formula is left empty.
#' The effect should be coded in such a way that value 1 indicated
#' the availability and 0 the absence.
#' In the case that some actors left or join the process at any point of time is
#' better to use the `present` variable in the node data frame linking to it the
#' time varying changes of composition of the actors set.
#' @param control_preprocessing a `preprocessing_opt.goldfish` object
#'  output from a call to [goldfish::set_preprocessing_opt()] to compute
#' the update statistics of the event sequence.
#' @param progress logical argument passed to
#' [goldfish::gather_model_data()] to show the progress of the preprocessing
#' of the events sequence.
#'
#' @return an object of class `"goldfish.latent.data"` that contains
#' a list with the following components.
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
#' @importFrom goldfish gather_model_data
#' @importFrom cli cli_abort
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
#' data2stan <- make_data_re(
#'   random_effects = list(inertia ~ 1),
#'   fixed_effects = callsDependent ~ recip + trans,
#'   data = socialEvolutionData
#' )
#' }
make_data_re <- function(
  random_effects,
  fixed_effects,
  model = c("DyNAM", "REM"),
  sub_model = c("choice", "rate", "choice_coordination"),
  data = NULL,
  support_constraint = NULL,
  control_preprocessing = NULL,
  progress = getOption("progress")
) {
  ### 0. check parameters----
  model <- match.arg(model)
  sub_model <- match.arg(sub_model)

  stopifnot(
    is.null(progress) || inherits(progress, "logical"),
    is.null(control_preprocessing) ||
      inherits(control_preprocessing, "preprocessing_opt.goldfish"),
    inherits(random_effects, "list"),
    inherits(fixed_effects, "formula"),
    is.null(support_constraint) ||
      inherits(support_constraint, "formula")
  )

  # setting initial values of some arguments
  if (is.null(progress)) progress <- FALSE
  if (is.null(control_preprocessing)) {
    control_preprocessing <- goldfish::set_preprocessing_opt()
  }

  # dependent network
  dep_net_name <- attr(data[[deparse(fixed_effects[[2]])]], "default_network")

  # formula treatment
  set_re_as_ego <- model == "DyNAM" && sub_model != "rate"
  extended_formula <- modify_formula(
    formula = fixed_effects,
    support_constraint = support_constraint,
    random_effects = random_effects,
    re_as_ego = set_re_as_ego,
    dep_net_name = dep_net_name
  )

  # create a full matrix for filtering
  processed_data <- goldfish::gather_model_data(
    formula = extended_formula$dynam_formula,
    model = if (model == "DyNAM" && sub_model == "choice") "DyNAMRE" else model,
    sub_model = sub_model,
    control_preprocessing = control_preprocessing,
    progress = progress,
    data = data
  )

  names_effects <- setNames(
    gsub("\\$", "Of", processed_data$namesEffects),
    extended_formula$dynam_terms
  )

  cstr_data <- make_df_cstr(
    processed_data = processed_data,
    extended_formula = extended_formula,
    names_effects = names_effects
  )
  expanded_df <- cstr_data$expanded_df
  processed_data$effectDescription <- cstr_data$effect_description

  # create objects for Stan
  n_total <- nrow(expanded_df)
  seq_ex_df <- seq.int(n_total)
  idx_events <- tapply(seq.int(n_total), expanded_df$event, range) |>
    simplify2array()
  senders_ix <- data.frame(label = sort(unique(expanded_df$sender))) |>
    within(index <- seq.int(label))
  expanded_df[, "sender_ix"] <-
    senders_ix[match(expanded_df[, "sender"], senders_ix[, "label"]), "index"]

  fe_idx <- match(extended_formula$base_labels, extended_formula$dynam_terms)
  formula_dynam_re <- mapply(
    function(effect, explanatory, names_effects, terms_dynam) {
      if (length(explanatory) > 0) {
        paste(
          names_effects[match(effect, terms_dynam)], "/",
          names_effects[match(explanatory, terms_dynam)]
        )
      } else {
        names_effects[match(effect, terms_dynam)]
      }
    },
    extended_formula$random_labels$lhs,
    extended_formula$random_labels$rhs,
    MoreArgs = list(
      names_effects = names_effects,
      terms_dynam = extended_formula$dynam_terms
    )
  ) |>
    c(names_effects[fe_idx]) |>
    paste(collapse = " + ")

  X_mat <- model.matrix(
    as.formula(paste("~ ", formula_dynam_re, " + 0")),
    data = expanded_df
  )

  re_names <- names_effects[unlist(extended_formula$random_labels$lhs)]
  Z_mat <- expanded_df[, re_names, drop = FALSE] |>
    as.matrix()

  data_stan <- list(
    T = ncol(idx_events),
    N = n_total,
    P = ncol(X_mat),
    Q = ncol(Z_mat),
    A = nrow(senders_ix),
    start = idx_events[1, ],
    end = idx_events[2, ],
    sender = expanded_df[, "sender_ix"],
    X = X_mat,
    Z = Z_mat,
    chose = which(expanded_df[, "selected"])
    # event = expanded_df[, "event"],
    # selected = expanded_df[, "selected"]
  )

  names_stan <- names(data_stan)
  names_change <- !grepl("^A|sender$", names_stan)
  names_stan[names_change] <- glue("{names_stan[names_change]}_{sub_model}")
  names(data_stan) <- names_stan

  names_effects <- cstr_data$names_effects
  extended_formula[["dynam_re_terms"]] <- formula_dynam_re
   
  return(structure(
    list(
      data_stan = data_stan,
      senders_ix = senders_ix,
      names_effects = names_effects,
      effect_description = processed_data$effectDescription,
      extended_formula = extended_formula
    ),
    class = c("DN_RE", "goldfish.latent.data"),
    model = "DN_RE",
    sub_model = sub_model,
    sample = FALSE
  ))
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
#' the supplementary material of
#' \insertCite{Merkle2019;textual}{goldfish.latent}.
#'
#' @param cmdstan_samples a `draws_array` or a `CmdStanFit` object with MCMC
#'   samples from the posterior distribution. In the case of a `draws_array`
#'   object is expected to have three dimensions corresponding to iterations,
#'   chain and variables.
#' @param data_stan a `list` output of a [CreateData()] call.
#' @param type a `character` value. It indicates whether the log-likelihood
#'   computation should return the `"conditional"` or the `"marginal"` version.
#' @param n_nodes an `integer`. The number of quadrature point to use in the
#'   Gauss-Hermite quadrature approximation use to integrate out the
#'   random-effects.
#' @param split_size an `integer` or `NULL`. It is use when the `type` is
#'   `"conditional"`. When it is `NULL`, the `split_size` is set to have
#'   roughly `4e4` rows sent to a processor when
#'   the matrix `X`, containing the change statistics, has more
#'   than `1e6`rows. If `X` has less than `1e6` rows, the `split_size` is set
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
#' callNetwork <- make_network(nodes = actors, directed = TRUE) |>
#'   link_events(change_event = calls, nodes = actors)
#' callsDependent <- define_dependent_events(
#'   events = calls, nodes = actors, default_network = callNetwork
#' )
#' socialEvolutionData <- make_data(callsDependent)
#' data2stan <- make_data_re(
#'   random_effects = list(inertia ~ 1),
#'   fixed_effects = callsDependent ~ recip + trans,
#'   data = socialEvolutionData
#' )
#'
#' stanCode <- make_model_code(data2stan)
#'
#' mod01 <- cmdstan_model(stanCode)
#' mod01Samples <- mod01$sample(
#'   data = data2stan[["dataStan"]],
#'   parallel_chains = 4, chains = 4,  iter_warmup = 500, iter_sampling = 500,
#'   show_messages = FALSE
#' )
#'
#' margLogLikMod01 <- compute_log_likelihood(mod01Samples, data2stan, spec = 4)
#' condLogLikMod01 <- compute_log_likelihood(mod01Samples, data2stan,
#'                                         type = "conditional", spec = 4)
#' }
compute_log_likelihood <- function(
  cmdstan_samples,
  data_stan,
  type = c("marginal", "conditional"),
  n_nodes = ifelse(type == "marginal", 11L, NULL),
  split_size = NULL,
  spec = parallel::detectCores() - 1,
  ...
) {
  stopifnot(
    inherits(cmdstan_samples, c("CmdStanFit", "draws")),
    inherits(data_stan, "goldfish.latent.data"),
    is.null(split_size) || inherits(split_size, "numeric") &&
      length(split_size) == 1
  )

  type <- match.arg(type)


  if (data_stan[["data_stan"]][["Qchoice"]] > 1)
    stop("Likelihood computation for a model with more than one random-effect",
         " is not yet available.")

  if (inherits(cmdstan_samples, "CmdStanFit")) {
    draws <- cmdstan_samples$draws("gamma_raw")
    drawsDimnames <- dimnames(draws)
    draws <- list(
      beta = cmdstan_samples$draws("betaChoice"),
      sigma = cmdstan_samples$draws("sigma"),
      gamma_raw = draws
    ) |>
      lapply(\(x) apply(x, 3, rbind))
  } else if (length(dim(cmdstan_samples)) == 3) {
    drawsDimnames <- dimnames(cmdstan_samples)

    variableNames <- dimnames(cmdstan_samples)[[3]]

    draws <- list(
      beta = cmdstan_samples[, , grepl("^beta", variableNames)],
      sigma = cmdstan_samples[, , grepl("^sigma$", variableNames)],
      gamma_raw = cmdstan_samples[, , grepl("^gamma_raw", variableNames)]
    ) |>
      lapply(\(x) apply(x, 3, rbind))
  } else
    stop(
      dQuote("cmdstan_samples"),
      " argument expects a three dimensional ", dQuote("draws"), " object."
    )

  if (is.null(split_size) & type == "conditional") {
    if (!is.numeric(spec) || length(spec) != 1)
      stop(
        "Please provide an integer number for", dQuote("split_size"),
        " parameter. It's not possible to assign it value with",
        "the current value of", dQuote("spec")
      )

    split_size <- if (spec == 1) NULL else
      ifelse(
        data_stan[["data_stan"]][["Nchoice"]] > 1e6,
        4e4 / data_stan[["data_stan"]][["A"]],
        data_stan[["data_stan"]][["Tchoice"]] / spec
      ) |> floor()

    eventsPerCore <- if (!is.null(split_size)) {
      parallel::splitIndices(
        data_stan[["data_stan"]][["Tchoice"]],
        floor(data_stan[["data_stan"]][["Tchoice"]] / split_size)
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
        seq_len(length(eventsPerCore)),
        fun = LogLikCondRE,
        draws = draws,
        dataList = data_stan[["data_stan"]],
        eventsPerCore = eventsPerCore
      )
      logLik <- Reduce(f = cbind, x = logLik)
    } else
      logLik <- LogLikCondRE(
        eventsIter = NULL,
        draws = draws,
        dataList = data_stan[["data_stan"]],
        eventsPerCore = NULL
      )

    drawsDimnames$variable <- sprintf("event[%d]", seq.int(ncol(logLik)))

  } else if (type == "marginal") {
    logLik <- mllDyNAMChoice(
      draws = draws,
      dataList = data_stan[["data_stan"]],
      n_nodes = n_nodes,
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
  dataList$senderEvent <- dataList$senderChoice[dataList$startChoice]
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
      keep <- seq.int(dataList$startChoice[event], dataList$endChoice[event])
      xb <- tcrossprod(dataList$Xchoice[keep, ], draws$beta)
      choice <- which(dataList$choseChoice[event] == keep)

      Z <- dataList$Zchoice[keep]

      mll <- mll +
        sapply(
          seq.int(nNodes),
          function(i) {
            utility <- sweep(xb, 1, Z * adaptNodes[i], FUN = "+")
            # # the log of the prob is utility - logSumExp: utility choice set
            utility[choice, ] - colLogSumExps(utility)

          }
        )
    }
    # # l_c + log(prob prior)
    mll <- mll + outer(
      draws$sigma[, 1],
      adaptNodes,
      \(x, y) dnorm(y, sd = x, log = TRUE)
    )
    # dnorm(adaptNodes[i], sd = draws$sigma, log = TRUE)

    # log(\prod \sum_{qdr points} lik (qdr point)) =
    # \sum logSumExp( log(log_lik (qdr point)))
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
    startI <- dataList$startChoice[head(events, 1)]
    endI <- dataList$endChoice[tail(events, 1)]
    X <- dataList$Xchoice[seq.int(startI, endI), ]

    seqEventsKeep <- seq.int(head(events, 1), tail(events, 1))
    start <- dataList$startChoice[seqEventsKeep] - (startI - 1)
    end <- dataList$endChoice[seqEventsKeep] - (startI - 1)
    chose <- dataList$choseChoice[seqEventsKeep] - (startI - 1)

    nE <- length(seqEventsKeep)
  } else {
    X <- dataList$Xchoice
    start <- dataList$startChoice
    end <- dataList$endChoice
    chose <- dataList$choseChoice
    nE <- dataList$Tchoice
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
    startI <- dataList$startChoice[head(events, 1)]
    endI <- dataList$endChoice[tail(events, 1)]
    seqDataKeep <- seq.int(startI, endI)
    X <- dataList$Xchoice[seqDataKeep, ]
    Z <- dataList$Zchoice[seqDataKeep]
    sender <- dataList$senderChoice[seqDataKeep]

    seqEventsKeep <- seq.int(head(events, 1), tail(events, 1))
    start <- dataList$startChoice[seqEventsKeep] - (startI - 1)
    end <- dataList$endChoice[seqEventsKeep] - (startI - 1)
    chose <- dataList$choseChoice[seqEventsKeep] - (startI - 1)

    nE <- length(seqEventsKeep)
  } else {
    X <- dataList$Xchoice
    Z <- dataList$Zchoice
    sender <- dataList$senderChoice
    start <- dataList$startChoice
    end <- dataList$endChoice
    chose <- dataList$choseChoice
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

