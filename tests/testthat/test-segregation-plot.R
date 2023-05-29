test_that("plotSegregation() catches errors", {

  x = nuclearPed()

  expect_silent(plotSegregation(x))
  expect_error(plotSegregation(x, proband = 1:2), "At most one proband is permitted")
  expect_error(plotSegregation(x, proband = 4), "Unknown proband:")

  expect_error(plotSegregation(x, affected = 4), "Unknown ID label:")
  expect_error(plotSegregation(x, unknown = 4), "Unknown ID label:")
  expect_error(plotSegregation(x, affected = 1, unknown = 1),
               "Individual specified as both affected and unknown")


  expect_error(plotSegregation(x, carriers = 4), "Unknown ID label:")
  expect_error(plotSegregation(x, noncarriers = 4), "Unknown ID label:")
  expect_error(plotSegregation(x, carriers = 1, noncarriers = 1),
               "Individual specified as both a carrier and a non-carrier")

  expect_error(plotSegregation(x, homozygous = 4), "Unknown ID label:")
  expect_error(plotSegregation(x, carriers = 1, homozygous = 1),
               "Individual specified as both a (heterozygous) carrier and homozygous", fixed = TRUE)
  expect_error(plotSegregation(x, noncarriers = 1, homozygous = 1),
               "Individual specified as both homozygous and a non-carrier")
})
