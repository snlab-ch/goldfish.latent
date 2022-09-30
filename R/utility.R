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
