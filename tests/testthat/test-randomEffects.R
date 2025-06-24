test_that("warnings and stops", {
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "dynam",
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "REM",
      sub_model = "choose",
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      model = "REM",
      sub_model = "choice_coordination",
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1, recip ~ 1),
      fixed_effects = depNetwork ~ trans,
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip * trans,
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ ego(actorsEx$attr1) * outdeg),
      fixed_effects = depNetwork ~ recip + trans,
      data = socialEvolutionData
    )
  )
  expect_error(
    make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      support_constraint = ~ tie(networkExog) + recip(networkExog),
      data = socialEvolutionData
    )
  )
})
test_that("choice empty model RE", {
  res <- make_data_re(
      random_effects = list(inertia ~ 1),
      fixed_effects = depNetwork ~ recip + trans,
      data = socialEvolutionData
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
  res <- make_data_re(
    random_effects = list(inertia ~ outdeg),
    fixed_effects = depNetwork ~ recip + trans,
    data = socialEvolutionData
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
  expect_error(make_model_code(list()))
  expect_error(make_model_code(data))
  attr(data, "subModel") <- "choice"
  outCode <- make_model_code(data)
  expect_length(outCode, 1)
  expect_type(outCode, "character")
})
