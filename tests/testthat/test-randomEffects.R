# defining objects
actorsEx <- goldfish::defineNodes(actorsEx) |>
  goldfish::linkEvents(changeEvent = compChange, attribute = "present")

networkState <- goldfish::defineNetwork(
  matrix = networkState, nodes = actorsEx,
  directed = TRUE) |>
  goldfish::linkEvents(changeEvent = eventsIncrement, nodes = actorsEx)
depNetwork <- goldfish::defineDependentEvents(
  events = eventsIncrement,
  nodes = actorsEx,
  defaultNetwork = networkState)

# define goldfish objects
networkExog <- goldfish::defineNetwork(
  matrix = networkExog,
  nodes = actorsEx, directed = TRUE) |>
  goldfish::linkEvents(changeEvent = eventsExogenous, nodes = actorsEx)

test_that("warnings and stops", {
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip + trans,
      model = "dynam"
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip + trans,
      model = "REM",
      subModel = "choose"
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip + trans,
      model = "REM",
      subModel = "choice_coordination"
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1, recip ~ 1),
      fixedEffects = depNetwork ~ trans
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip * trans
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ ego(actorsEx$attr1) * outdeg),
      fixedEffects = depNetwork ~ recip + trans
    )
  )
  expect_error(
    CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip + trans,
      supportConstraint = ~ tie(networkExog) + recip(networkExog)
    )
  )
})
test_that("choice empty model RE", {
  res <- CreateData(
      randomEffects = list(inertia ~ 1),
      fixedEffects = depNetwork ~ recip + trans
  )
  expect_type(res, "list")
  expect_length(res, 4)
  expect_equal(
    res$sendersIx,
    data.frame(
      label = sprintf("Actor %d", 1:5),
      index = 1:5
    )
  )
  expect_equal(
    res$namesEffects,
    c("inertia" = "inertia_networkState", "recip" = "recip_networkState",
      "trans" = "trans_networkState")
  )
  expect_equal(res$dataStan$T, 12)
  expect_equal(res$dataStan$P, 3)
  expect_equal(res$dataStan$A, nrow(actorsEx))
})
test_that("choice RE with expl effects", {
  res <- CreateData(
    randomEffects = list(inertia ~ outdeg),
    fixedEffects = depNetwork ~ recip + trans
  )
  expect_type(res, "list")
  expect_length(res, 4)
  expect_equal(
    res$sendersIx,
    data.frame(
      label = sprintf("Actor %d", 1:5),
      index = 1:5
    )
  )
  expect_equal(
    res$namesEffects,
    c('outdeg(networkState, type = "ego")' = "outdeg_networkState_ego",
      "inertia" = "inertia_networkState", "recip" = "recip_networkState",
      "trans" = "trans_networkState")
  )
  expect_equal(res$dataStan$T, 12)
  expect_equal(res$dataStan$P, 4)
  expect_equal(res$dataStan$A, nrow(actorsEx))
})
