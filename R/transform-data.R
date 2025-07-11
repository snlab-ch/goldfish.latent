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

#' @importFrom glue glue
make_cox_model_data <- function(data, sub_model = c("rate", "choice")) {
  sub_model <- match.arg(sub_model)
  data_stan <- data[["data_stan"]]
  model <- attr(data, "model")
  is_sample <- attr(data, "sample")

  expanded_df <- make_df_from_stan(data_stan, sub_model)

  if (model %in% c("DN_RE", "DN_HMM_RE", "DN_CHMM_RE")) {
    expanded_df$sender <- data_stan[["sender"]]
  }

  if (is_sample) {
    expanded_df$offset_sample <- data_stan[[glue("offset_set_{sub_model}")]]
  }

  if (sub_model == "choice") {
    expanded_df$uno <- 1
  }

  if (sub_model == "rate") {
    expanded_df$timespan <-
      data_stan[[glue("timespan_{sub_model}")]][expanded_df$event]
  }

  expanded_df
}

#' Create a data frame and apply support constraint
#'
#' @param processed_data a list with preprocess data according to the model
#'  formula. Output of [goldfish::gather_model_data()].
#' @param extended_formula a list with a formula used for preprocessing
#' that identifies the different components. Output of [modify_formula()].
#' @param names_effects a named character vector with the names of the effects
#' after preprocessing and the formula of the effects as names.
#'
#' @return a list with the following components:
#' \describe{
#'   \item{expanded_df}{a data frame with the data in long format.}
#'   \item{names_effects}{a character vector with the names of the effects
#'     after removing the support constraint.}
#'   \item{effect_description}{an array with detailed and comprehensible
#'   information of the terms used in the random and fixed effects formulas
#'   after removing the support constraint.}
#' }
#' @noRd
make_df_cstr <- function(
    processed_data,
    extended_formula,
    names_effects) {
  expanded_df <- with(processed_data, {
    n_events <- length(sender)

    expanded_df <- cbind(
      setNames(
        as.data.frame(stat_all_events),
        names_effects
      ),
      data.frame(
        event = rep(seq.int(n_events), n_candidates),
        selected = sequence(n_candidates) ==
          rep(selected, n_candidates),
        sender = rep(sender, n_candidates)
      )
    )

    # subset if constraint
    if (!is.null(extended_formula$cstr_label)) {
      cstr_name <- names_effects[extended_formula$cstr_label]
      keep <- expanded_df[, cstr_name] == 1
      expanded_df <- expanded_df[keep, !names(expanded_df) %in% cstr_name]
    }
    expanded_df
  })

  n_selected <- sum(expanded_df$selected)
  if (!is.null(processed_data$isDependent)) {
    n_origin <- length(processed_data$selected[processed_data$isDependent])
  } else {
    n_origin <- length(processed_data$selected)
  }
  if (n_selected != n_origin) {
    cli::cli_warn(c(
      cli_text(
        "There is a mismatch between the number of receivers ",
        "before and after constraint."),
      "i" = "receivers before: {n_origin}, receivers after: {n_selected}"
    ))
  }
  effect_description <- processed_data$effectDescription
  if (!is.null(extended_formula$cstr_label)) {
    cstr_name <- names_effects[extended_formula$cstr_label]
    effect_description <- effect_description[!names_effects %in% cstr_name, ]
    names_effects <- names_effects[!names_effects %in% cstr_name]
  }

  list(
    expanded_df = expanded_df,
    names_effects = names_effects,
    effect_description = effect_description
  )
}

#' Create an expanded data frame from stan data
#'
#' @param data_stan a list with data for Stan.
#' @param sub_model a character string indicating the sub_model,
#'   either "rate" or "choice".
#'
#' @return a data frame with the data in long format.
#' @noRd
make_df_from_stan <- function(data_stan, sub_model = c("rate", "choice")) {
  sub_model <- match.arg(sub_model)

  X_mat_name <- glue("X_{sub_model}")
  n_total_name <- glue("N_{sub_model}")
  n_events_name <- glue("T_{sub_model}")
  start_idx_name <- glue("start_{sub_model}")
  end_idx_name <- glue("end_{sub_model}")
  chose_vec_name <- glue("chose_{sub_model}")

  expanded_df <- as.data.frame(data_stan[[X_mat_name]])
  expanded_df$event <- rep(
    seq_len(data_stan[[n_events_name]]),
    data_stan[[end_idx_name]] - data_stan[[start_idx_name]] + 1
  )
  expanded_df$selected <- seq_len(data_stan[[n_total_name]]) %in%
    data_stan[[chose_vec_name]]

  expanded_df
}

# Helper function to process a submodel
sample_set <- function(
  data_stan,
  model, sub_model,
  method_set,
  fraction_set, sample_size_set = NULL,
  sample_events = NULL
) {
  data_sample <- list()
  x_names <- colnames(data_stan[[glue("X_{sub_model}")]])
  z_names <- colnames(data_stan[[glue("Z_{sub_model}")]])

  expanded_df <- make_df_from_stan(data_stan, sub_model)
  if (model %in% c("DN_RE", "DN_HMM_RE", "DN_CHMM_RE")) {
    expanded_df$sender <- data_stan[["sender"]]
  }

  is_dependent <- data_stan[["is_dependent"]]
  if (!is.null(sample_events) && sample_events[["sel_prob"]] < 1) {
    s_events <- sample_events$sample
    n_events_sample <- length(s_events)
    data_sample[[glue("T_{sub_model}")]] <- n_events_sample

    if (sub_model == "rate") {
      data_sample[["is_dependent"]] <- is_dependent[s_events]
      data_sample[["timespan"]] <- data_stan[["timespan"]][s_events]
    }
    expanded_df <- subset(expanded_df, event %in% s_events)
    data_sample[[glue("fraction_event_{sub_model}")]] <- sample_events$sel_prob
  } else {
    n_events_total <- data_stan[[glue("T_{sub_model}")]]
    data_sample[[glue("T_{sub_model}")]] <- n_events_total
    if (sub_model == "rate") {
      data_sample["is_dependent"] <- data_stan[["is_dependent"]]
      data_sample["timespan"] <- data_stan[["timespan"]]
    }
    data_sample[[glue("fraction_event_{sub_model}")]] <- 1
  }

  should_sample <- fraction_set < 1 || !is.null(sample_size_set)
  if (should_sample) {
    sample_df <- by(
      expanded_df,
      expanded_df[, "event"],
      function(x) {
        x[, "offset_sample"] <- 1
        selected_row <- x[x$selected, ]
        non_selected_rows <- x[!x$selected, ]
        if (nrow(non_selected_rows) > 0) {
          sample <- do.call(
            method_set,
            list(
              N = nrow(non_selected_rows),
              fraction = fraction_set,
              sample_size = sample_size_set
            ))
          non_selected_rows <- non_selected_rows[sample$sample, ]
          non_selected_rows[, "offset_sample"] <- 1 / sample$sel_prob
        }
        if (nrow(selected_row) == 0) selected_row <- NULL
        rbind(selected_row, non_selected_rows)
      }
    ) |>
    Reduce(f = rbind, x = _)
  }
  data_sample[[glue("offset_set_{sub_model}")]] <- sample_df$offset_sample

  data_sample[[glue("N_{sub_model}")]] <- nrow(sample_df)
  idx_events <- tapply(seq_len(nrow(sample_df)), sample_df$event, range) |>
    simplify2array()
  data_sample[[glue("start_{sub_model}")]] <- idx_events[1, ]
  data_sample[[glue("end_{sub_model}")]] <- idx_events[2, ]
  data_sample[[glue("X_{sub_model}")]] <- as.matrix(sample_df[, x_names])
  data_sample[[glue("P_{sub_model}")]] <- data_stan[[glue("P_{sub_model}")]]

  if (model %in% c("DN_RE", "DN_HMM_RE", "DN_CHMM_RE")) {
    data_sample[[glue("Z_{sub_model}")]] <- as.matrix(sample_df[, z_names])
    data_sample[[glue("Q_{sub_model}")]] <- data_stan[[glue("Q_{sub_model}")]]
    data_sample[["A"]] <- data_stan[["A"]]
    data_sample[["sender"]] <- sample_df[, "sender"]
  }

  chose_name <- glue("chose_{sub_model}")
  chose_idx <- which(sample_df$selected)
  if (sub_model == "rate") {
    data_sample[[chose_name]] <- rep(0, n_events_sample)
    data_sample[[chose_name]][data_sample[["is_dependent"]] == 1] <- chose_idx
  } else {
    data_sample[[chose_name]] <- chose_idx
  }

  data_sample
}
