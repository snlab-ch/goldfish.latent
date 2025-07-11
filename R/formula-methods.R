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

#' Modify formula
#'
#' This function modifies the formula to include the support constraint and
#' random effects.
#' @param formula a formula object
#' @param support_constraint a formula object
#' @param random_effects a list of formula objects
#' @param re_as_ego a logical value
#' @param dep_net_name the name of the dependent network
#' @return a list with the following elements:
#'   - dynam_formula: the modified formula to use in goldfish for preprocessing
#'   - dynam_terms: the terms of the modified formula
#'   - base_labels: the base labels of the modified formula
#'   - random_terms: the random terms of the modified formula
#'   - random_labels: the random labels of the modified formula
#'   - cstr_terms: the support constraint terms of the modified formula
#'   - cstr_label: the support constraint label of the modified formula
#' @noRd
modify_formula <- function(
    formula,
    support_constraint = NULL,
    random_effects = NULL,
    re_as_ego = FALSE,
    dep_net_name = "callsDependent") {
  formula_terms <- terms(formula)
  has_intercept <- has_explicit_intercept(formula)
  order_formula <- attr(formula_terms, "order")
  base_labels <- attr(formula_terms, "term.labels")
  if (any(order_formula != 1)) {
    cli_abort(c(
      "The {.var formula} argument doesn't support interactions yet",
      "x" = "There are {sum(order_formula != 1)} interactions in the formula."
    ))
  }

  if (!is.null(support_constraint)) {
    cstr_terms <- terms(support_constraint)
    cstr_label <- attr(cstr_terms, "term.labels")
    len_cstr <- length(cstr_label)
    if (len_cstr > 1) {
      cli_abort(c(
        "{.var support_constraint} must have only one effect.",
        "x" = "There are {len_cstr} terms in the formula."
      ))
    }
  } else {
    cstr_terms <- NULL
    cstr_label <- NULL
  }

  if (!is.null(random_effects)) {
    check_re <- vapply(random_effects, \(x) inherits(x, "formula"), logical(1))
    if (!all(check_re)) {
      class_re <- vapply(random_effects, class, character(1))
      cli_abort(c(
        "All elements of {.var random_effects} must be formula objects.",
        "x" = "You've have supplied a list with {.cls {class_re}} objects."
      ))
    }
    check_re <- vapply(
      random_effects,
      \(x) any(attr(terms(x), "order") != 1),
      logical(1)
    )
    if (any(check_re)) {
      cli_abort(c(
        "The {.var random_effects} argument doesn't support interactions yet",
        "x" = "There are {sum(check_re)} formulas with interactions."
      ))
    }
    if (re_as_ego) {
      random_effects <- lapply(
        random_effects,
        modify_formula_re,
        dep_net_name = dep_net_name
      )
    }

    random_terms <- lapply(random_effects, terms)
    random_labels <- list(
      lhs = lapply(random_terms, \(x) deparse(x[[2]])),
      rhs = lapply(random_terms, \(x) attr(x, "term.labels"))
    )
    extended_terms <- unlist(random_labels, use.names = FALSE)
  } else {
    random_terms <- NULL
    random_labels <- NULL
    extended_terms <- NULL
  }

  terms_dynam <- c(base_labels, cstr_label, extended_terms) |> unique()
  if (has_intercept) {
    terms_dynam <- c("1", terms_dynam)
  }

  formula_dynam <- stats::reformulate(
    terms_dynam,
    response = as.character(formula[[2]])
  )

  formula_extend <- list(
    dynam_formula = formula_dynam,
    dynam_terms = terms_dynam,
    base_labels = base_labels,
    random_terms = random_terms,
    random_labels = random_labels,
    cstr_terms = cstr_terms,
    cstr_label = cstr_label,
    has_intercept = has_intercept
  )

  return(formula_extend)
}

#' Modify the random effects formula
#'
#' This function modifies the random effects formula to include the dependent
#' network as the first term and set the type argument to `"ego"`.
#' @param re_formula a formula object
#' @param dep_net_name the name of the dependent network
#' @return a formula object
#' @noRd
modify_formula_re <- function(re_formula, dep_net_name) {
  stopifnot(inherits(re_formula, "formula"))

  re_terms <- terms(re_formula)
  effects_formula <- attr(re_terms, "term.labels")
  if (length(effects_formula) == 0) {
    return(re_formula)
  }

  # modify calls
  effect_for_mod <- vapply(
    effects_formula,
    \(x) {
      effect_lang <- str2lang(x)
      if (length(effect_lang) == 1) {
        deparse(as.call(
          list(effect_lang, as.symbol(dep_net_name), type = "ego")
        ))
      } else if (deparse(effect_lang[[1]]) != "ego") {
        deparse(as.call(
          c(
            list(effect_lang[[1]]),
            as.list(effect_lang[-1]),
            list(type = "ego")
          )
        ))
      } else {
        x
      }
    },
    character(1)
  )

  # formula after modifications
  reformulate(effect_for_mod, response = re_formula[[2]])
}

#' Check if a formula has an explicit intercept
#'
#' This function fails to identify that the intercept is repalced by the
#' inclusion of a factor covariate including all levels.
#' @param formula a formula object
#' @return a logical value
#' @noRd
has_explicit_intercept <- function(formula) {
  rhs <- deparse(formula[[3]])
  # Remove whitespace and collapse across lines if needed
  rhs <- gsub("\\s+", "", paste(rhs, collapse = ""))
  # Detect +1, starting 1 or bare 1 (only term)
  grepl("(^|\\+)1(\\W|$)", rhs)
}
