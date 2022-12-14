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
test_that("save code", {
  data <- structure(
    list(dataStan = list(Q = 2)),
    class = "goldfish.latent.data",
    model = "DyNAM", subModel = "rate"
  )
  expect_error(CreateModelCode(list()))
  expect_error(CreateModelCode(data))
  attr(data, "subModel") <- "choice"
  outCode <- CreateModelCode(data)
  expect_length(outCode, 1)
  expect_type(outCode, "character")
})
