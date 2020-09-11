# Validate the examples of Thompson et al (2003)

PENETR = c(0.05, 0.7, 0.7)

test_that("Thompson's results for Pedigree 1 are correct", {
  x = relabel(swapSex(cousinPed(1), c(3,8)), "asPlot")

  quickFLB = function(carr, noncarr = NULL, aff, unkn = NULL) {
    bf = FLB(x, carriers = carr, noncarriers = noncarr,
             affected = aff, unknown = unkn,
             freq = 0.001,
             proband = "7",
             penetrances = PENETR,
             plot = F)
    round(bf, 2)
  }

  expect_equal(7.37, quickFLB(carr = c(4,5,7,8), aff = c(4,5,7,8)))

  expect_equal(7.30, quickFLB(carr = 7:8, aff = c(4,5,7,8)))
  expect_equal(5.28, quickFLB(carr = 7:8, aff = 7:8, unkn = 4:5))
  expect_equal(1.83, quickFLB(carr = 7:8, aff = 7:8))

  expect_equal(0.09, quickFLB(carr = 7, noncarr = 8, aff = c(4,5,7,8)))
  expect_equal(0.38, quickFLB(carr = 7, noncarr = 8, aff = 7:8, unkn = 4:5))
  expect_equal(0.88, quickFLB(carr = 7, noncarr = 8, aff = 7:8))
})


test_that("Thompson's results for Pedigree 2 are correct", {
  x = addChildren(nuclearPed(3, sex = c(1,2,2)), father = 3, nch = 2, sex = 2:1)
  # plot(x)

  bf = FLB(x,
           carriers = 7:8,
           noncarriers = 6,
           freq = 0.001,
           proband = 7,
           penetrances = PENETR,
           affected = c(4,5,6,7,8),
           plot = F)

  expect_equal(1.03, round(bf, 2))
})

