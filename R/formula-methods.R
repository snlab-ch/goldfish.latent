modify_formula <- function(
    formula,
    support_constraint = NULL,
    random_effects = NULL,
    re_as_ego = FALSE,
    dep_net_name = "callsDependent") {
  formula_terms <- terms(formula)
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
    extended_terms <- NULL
  }

  terms_dynam <- c(base_labels, cstr_label, extended_terms) |> unique()

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
    cstr_label = cstr_label
  )

  return(formula_extend)
}

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
