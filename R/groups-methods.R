#' Gather data from multiple groups into a single object
#'
#' This function aggregates preprocessed data from a list of groups into a
#' single data structure suitable for a random effects by group model (using
#' Stan).
#' It combines matrices and vectors from each group, adjusting indices as
#' needed.
#'
#' @param groups_info A data frame with metadata about the groups. It must
#'   contain at least a column `ixGroup` with the group identifier.
#' @param data A list of preprocessed data objects. Each element of the list
#'   corresponds to a group and should be a `goldfish.latent.data` object
#'   containing a `data_stan` list.
#' @param keep_x An optional character vector of column names to keep from
#'   the `X_choice` matrix. If `NULL`, all columns are kept.
#' @param scale logical value. Whether to standardize the effect stats at
#'   the end.
#'
#' @return A `goldfish.latent.data` object containing the aggregated data in the
#'   `dataStan` list, along with metadata from the first group.
#' @export
gather_groups <- function(
    groups_info,
    data,
    keep_x = NULL,
    scale = FALSE
) {
  # Validate inputs
  if (!is.list(data)) {
    cli::cli_abort(c(
      "{.arg data} must be a list.",
      "x" = "You've supplied a {.cls {class(data)}}."
    ))
  }

  if (!is.data.frame(groups_info)) {
    cli::cli_abort(c(
      "{.arg groups_info} must be a data frame.",
      "x" = "You've supplied a {.cls {class(groups_info)}}."
    ))
  }

  if (!is.null(keep_x) && !is.character(keep_x)) {
    cli::cli_abort(c(
      "{.arg keep_x} must be a character vector or NULL.",
      "x" = "You've supplied a {.cls {class(keep_x)}}."
    ))
  }

  if (length(data) != nrow(groups_info)) {
    cli::cli_abort(c(
      "Length of {.arg data} must match number of rows in {.arg groups_info}.",
      "x" = "Length of {.arg data} is {.val {length(data)}},",
      "i" = "but {.arg groups_info} has {.val {nrow(groups_info)}} rows."
    ))
  }

  # set model variables and features
  model <- attr(data[[1]], "model")
  sub_model <- attr(data[[1]], "sub_model")
  has_sample <- attr(data[[1]], "sample")

  # Get dimensions and names from the first group to establish a template
  first_group_data <- data[[1]]$data_stan
  x_name <- glue("X_{sub_model}")
  z_name <- glue("Z_{sub_model}")
  n_name <- glue("N_{sub_model}")
  t_name <- glue("T_{sub_model}")
  start_name <- glue("start_{sub_model}")
  end_name <- glue("end_{sub_model}")
  chose_name <- glue("chose_{sub_model}")
  col_names_x <- colnames(first_group_data[[x_name]])

  # Calculate total dimensions for the aggregated data structures
  n_size <- sum(sapply(data, \(x) x$data_stan[[n_name]]))
  t_size <- sum(sapply(data, \(x) x$data_stan[[t_name]]))
  p_size <- length(col_names_x)

  # Optionally, filter covariates to keep
  if (!is.null(keep_x)) {
    if (!all(keep_x %in% col_names_x))
      cli::cli_abort(c(
        "Some variables in {.arg keep_x} are not in the design matrix.",
        "x" = "You've supplied {setdiff(keep_x, col_names_x)}."
      ))
    p_size <- length(keep_x)
    col_names_x <- keep_x
  }

  q_size <- max(sapply(data, \(x) x$data_stan[[glue("Q_{sub_model}")]]))
  has_z <- !is.null(first_group_data[[z_name]])

  # Initialize aggregated data structures
  x_matrix <- matrix(
    0, nrow = n_size, ncol = p_size,
    dimnames = list(NULL, col_names_x)
  )
  if (has_z) {
    z_matrix <- matrix(
      0, nrow = n_size, ncol = q_size,
      dimnames = list(NULL, colnames(first_group_data[[z_name]]))
    )
  }

  n_groups <- nrow(groups_info)

  start_agg <- integer(t_size)
  end_agg <- integer(t_size)
  chose_agg <- integer(t_size)
  sender_agg <- integer(n_size)
  start_group_agg <- integer(n_groups)

  if (sub_model == "rate") {
    timespan_agg <- integer(t_size)
    is_dependent_agg <- integer(t_size)
  }

  # Initialize offsets for rows in choice set (n) and events (t)
  n_offset <- 0
  t_offset <- 0

  # Loop through each group to populate the aggregated structures
  for (gr in seq_len(n_groups)) {
    group_data <- data[[gr]]$data_stan
    group_label <- groups_info$ixGroup[gr]

    # Check for consistency in design matrices across groups
    current_col_names <- colnames(group_data[[x_name]])
    if (!is.null(keep_x)) {
      missing_names <- setdiff(keep_x, current_col_names)
      if (length(missing_names) > 0) {
        cli::cli_abort(c(
          "Group {group_label} has missing variables defined in",
          "{.var keep_x}.",
          "x" = "Missing: {missing_names}"
        ))
      }
    } else {
      if (!identical(col_names_x, current_col_names)) {
        not_first <- setdiff(col_names_x, current_col_names)
        not_current <- setdiff(current_col_names, col_names_x)
        cli::cli_abort(c(
          "Variables do not match among groups:",
          "i" = "In first group but not in current: {not_first}",
          "i" = "In current but not in first: {not_current}"
        ))
      }
    }

    # Get a slice of the design matrix
    x_aux <- group_data[[x_name]][, col_names_x, drop = FALSE]

    # Determine row indices for the current group in the aggregated matrices
    n_indices <- (n_offset + 1):(n_offset + group_data[[n_name]])
    t_indices <- (t_offset + 1):(t_offset + group_data[[t_name]])

    # Assign group data to the correct slice of the aggregated matrices
    x_matrix[n_indices, ] <- x_aux
    if (has_z) z_matrix[n_indices, ] <- group_data[[z_name]]

    # Adjust and assign event data, shifting indices by the current offset
    start_agg[t_indices] <- group_data[[start_name]] + n_offset
    end_agg[t_indices] <- group_data[[end_name]] + n_offset
    chose_agg[t_indices] <- group_data[[chose_name]] + n_offset

    if (sub_model == "rate") {
      is_dependent <- group_data$is_dependent
      timespan_agg[t_indices] <- group_data$timespan
      is_dependent_agg[t_indices] <- is_dependent
      chose_agg[t_indices[!is_dependent]] <- 0
    }

    sender_agg[n_indices] <- groups_info$ixGroup[gr]
    start_group_agg[gr] <- t_indices[1]

    # Update offsets for the next iteration
    n_offset <- n_offset + group_data[[n_name]]
    t_offset <- t_offset + group_data[[t_name]]
  }

  # Assemble the final data object
  data_stan <- list(
    T = t_size,
    N = n_size,
    P = p_size,
    Q = q_size,
    start = start_agg,
    end = end_agg,
    X = x_matrix,
    Z = if (has_z) z_matrix else NULL,
    chose = chose_agg
  )
  names(data_stan) <- glue("{name}_{sub_model}", name = names(data_stan))
  data_stan[["A"]] <- n_groups
  data_stan[["sender"]] <- sender_agg
  data_stan[["start_group"]] <- start_group_agg

  if (sub_model == "rate") {
    data_stan[["timespan"]] <- timespan_agg
    data_stan[["is_dependent"]] <- is_dependent_agg
    mean_actors <- mean(end_agg - start_agg + 1)
    total_time <- sum(timespan_agg)
    data_stan[["offset_int_rate"]] <- log(t_size / (total_time * mean_actors))
  }

  data_gathered <- list(
    data_stan = data_stan,
    namesEffects = data[[1]]$namesEffects,
    effectsDescription = data[[1]]$effectsDescription,
    groupInfo = groups_info
  )

  # Preserve class and attributes from the original data object
  class(data_gathered) <- class(data[[1]])
  attr(data_gathered, "model") <- model
  attr(data_gathered, "subModel") <- sub_model
  attr(data_gathered, "sample") <- has_sample

  return(data_gathered)
}

#' Gather data from multiple groups to a JSON file
#' 
#' @param groups_info a data frame with the groups information.
#' @param data a list with the data from each group.
#' @param file_name the name of the file to save the data,
#'   with `.json` extension.
#' @param keep_x a character vector with the names of the variables to keep.
#' @param scale a logical indicating whether to scale the data.
#' @return a list with the data in long format.
#' @export
#' @importFrom stringr str_sub
#' @importFrom jsonlite toJSON
gather_groups_to_json <- function(
    groups_info,
    data,
    file_name,
    keep_x = NULL,
    scale = FALSE
) {
  if (!is.list(data)) {
    cli::cli_abort(c(
      "{.arg data} must be a list.",
      "x" = "You've supplied a {.cls {class(data)}}."
    ))
  }

  if (!is.data.frame(groups_info)) {
    cli::cli_abort(c(
      "{.arg groups_info} must be a data frame.",
      "x" = "You've supplied a {.cls {class(groups_info)}}."
    ))
  }

  if (!is.null(keep_x) && !is.character(keep_x)) {
    cli::cli_abort(c(
      "{.arg keep_x} must be a character vector or NULL.",
      "x" = "You've supplied a {.cls {class(keep_x)}}."
    ))
  }

  if (length(data) != nrow(groups_info)) {
    cli::cli_abort(c(
      "Length of {.arg data} must match number of rows in {.arg groups_info}.",
      "x" = "Length of {.arg data} is {.val {length(data)}},",
      "i" = "but {.arg groups_info} has {.val {nrow(groups_info)}} rows."
    ))
  }

  # set model variables and features
  model <- attr(data[[1]], "model")
  sub_model <- attr(data[[1]], "sub_model")
  has_sample <- attr(data[[1]], "sample")

  # Get dimensions and names from the first group to establish a template
  first_group_data <- data[[1]]$data_stan
  x_name <- glue("X_{sub_model}")
  z_name <- glue("Z_{sub_model}")
  n_name <- glue("N_{sub_model}")
  t_name <- glue("T_{sub_model}")
  start_name <- glue("start_{sub_model}")
  end_name <- glue("end_{sub_model}")
  chose_name <- glue("chose_{sub_model}")
  col_names_x <- colnames(first_group_data[[x_name]])

  # Calculate total dimensions for the aggregated data structures
  n_size <- sum(sapply(data, \(x) x$data_stan[[n_name]]))
  t_size <- sum(sapply(data, \(x) x$data_stan[[t_name]]))
  p_size <- length(col_names_x)

  # Optionally, filter covariates to keep
  if (!is.null(keep_x)) {
    if (!all(keep_x %in% col_names_x))
      cli::cli_abort(c(
        "Some variables in {.arg keep_x} are not in the design matrix.",
        "x" = "You've supplied {setdiff(keep_x, col_names_x)}."
      ))
    p_size <- length(keep_x)
    col_names_x <- keep_x
  }

  q_size <- max(sapply(data, \(x) x$data_stan[[glue("Q_{sub_model}")]]))
  has_z <- !is.null(first_group_data[[z_name]])

  # Initialize text containers
  n_groups <- nrow(groups_info)

  start_text <- character(n_groups + 1)
  end_text <- character(n_groups + 1)
  chose_text <- character(n_groups + 1)
  sender_text <- character(n_groups + 1)
  start_group_text <- character(n_groups + 1)
  x_text <- character(n_groups + 1)
  if (has_z) z_text <- character(n_groups + 1)

  if (sub_model == "rate") {
    timespan_text <- character(n_groups + 1)
    is_dependent_text <- character(n_groups + 1)
  }

  # Initialize offsets and total values
  n_offset <- 0
  t_offset <- 0
  total_actors <- 0
  total_time <- 0


  # Initialize containers
  start_text[1] <- glue('"start_{sub_model}": [')
  end_text[1] <- glue('"end_{sub_model}": [')
  chose_text[1] <- glue('"chose_{sub_model}": [')
  sender_text[1] <- glue('"sender_{sub_model}": [')
  start_group_text[1] <- glue('"start_group_{sub_model}": [')
  x_text[1] <- glue('"X_{sub_model}": [')
  if (has_z) z_text[1] <- glue('"Z_{sub_model}": [')

  if (sub_model == "rate") {
    timespan_text[1] <- '"timespan": ['
    is_dependent_text[1] <- '"is_dependent": ['
  }

  # Loop through each group to populate the aggregated structures
  for (gr in seq_len(n_groups)) {
    group_data <- data[[gr]]$data_stan
    group_label <- groups_info$ixGroup[gr]
    gr_idx <- gr + 1
    last_replace <- if (gr == n_groups) "\n]," else ","

    # Check for consistency in design matrices across groups
    current_col_names <- colnames(group_data[[x_name]])
    if (!is.null(keep_x)) {
      missing_names <- setdiff(keep_x, current_col_names)
      if (length(missing_names) > 0) {
        cli::cli_abort(c(
          "Group {group_label} has missing variables defined in",
          "{.var keep_x}.",
          "x" = "Missing: {missing_names}"
        ))
      }
    } else {
      if (!identical(col_names_x, current_col_names)) {
        not_first <- setdiff(col_names_x, current_col_names)
        not_current <- setdiff(current_col_names, col_names_x)
        cli::cli_abort(c(
          "Variables do not match among groups:",
          "i" = "In first group but not in current: {not_first}",
          "i" = "In current but not in first: {not_current}"
        ))
      }
    }

    # Get a slice of the design matrix
    x_aux <- group_data[[x_name]][, col_names_x, drop = FALSE]
    x_text[gr_idx] <- transform_json(
      x_aux,
      first_position = 3,
      last_position = -2, last_replace = last_replace
    )
    if (has_z) {
      z_text[gr_idx] <- transform_json(
        group_data[[z_name]],
        first_position = 3,
        last_position = -2,
        last_replace = last_replace
      )
    }

    # Adjust and assign event data, shifting indices by the current offset
    start_text[gr_idx] <- transform_json(
      group_data[[start_name]] + n_offset,
      last_replace = last_replace
    )
    end_text[gr_idx] <- transform_json(
      group_data[[end_name]] + n_offset,
      last_replace = last_replace
    )

    chose_tmp <- group_data[[chose_name]] + n_offset
    if (sub_model == "rate") {
      is_dependent <- group_data$is_dependent
      timespan_text[gr_idx] <- transform_json(
        group_data$timespan,
        last_replace = last_replace
      )
      is_dependent_text[gr_idx] <- transform_json(
        is_dependent,
        last_replace = last_replace
      )

      chose_tmp[!is_dependent] <- 0
      total_time <- total_time + sum(group_data$timespan)
      total_actors <- total_actors +
        sum(group_data[[end_name]] - group_data[[start_name]] + 1)
    }
    chose_text[gr_idx] <- transform_json(chose_tmp, last_replace = last_replace)

    sender_text[gr_idx] <- rep(groups_info$ixGroup[gr], group_data[[n_name]]) |>
      transform_json(last_replace = last_replace)
    start_group_text[gr_idx] <- glue("{t_offset + 1}{last_replace}")

    # Update offsets for the next iteration
    n_offset <- n_offset + group_data[[n_name]]
    t_offset <- t_offset + group_data[[t_name]]
  }

  # Assemble the final data object
  data_stan <- list(
    T = t_size,
    N = n_size,
    P = p_size,
    Q = q_size
  )
  names(data_stan) <- glue("{name}_{sub_model}", name = names(data_stan))
  data_stan[["A"]] <- n_groups
  
  if (sub_model == "rate") {
    data_stan[["offset_int_rate"]] <-
      log(t_size^2 / (total_time * total_actors))
  }

  data_text <- transform_json(
    data_stan,
    first_position = 3,
    last_replace = "}"
  )

  # Write the data to a json file
  conn <- file(file_name, open = "wb")
  writeLines("{", conn, useBytes = TRUE)
  writeLines(start_text, conn)
  writeLines(end_text, conn)
  writeLines(chose_text, conn)
  writeLines(sender_text, conn)
  writeLines(start_group_text, conn)
  writeLines(x_text, conn)
  if (has_z) writeLines(z_text, conn)
  if (sub_model == "rate") {
    writeLines(timespan_text, conn)
    writeLines(is_dependent_text, conn)
  }
  writeLines(data_text, conn)
  close(conn)


  data_gathered <- list(
    json_file = file_name,
    namesEffects = data[[1]]$namesEffects,
    effectsDescription = data[[1]]$effectsDescription,
    groupInfo = groups_info
  )

  # Preserve class and attributes from the original data object
  class(data_gathered) <- class(data[[1]])
  attr(data_gathered, "model") <- model
  attr(data_gathered, "subModel") <- sub_model
  attr(data_gathered, "sample") <- has_sample
  attr(data_gathered, "scale") <- scale
  attr(data_gathered, "json_file") <- TRUE

  return(data_gathered)
} 