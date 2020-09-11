# Validate the examples of Bayrak-Toydemir et al (2008)

skip("Bayrak-Toydemir not working yet")

#--- SETUP ---

penetr = c(0.0015, 0.95, 0.95)
q = 0.001

# penetr = c(0, 1, 1)
# q = 1e-15

quickFLB = function(x, aff, unknown = NULL, carriers, noncarriers = NULL,
                    proband, plot = F) {
  res = FLB(x, affected = aff, freq = q, unknown = unknown,
            carriers = carriers, noncarriers = noncarriers,
            proband = proband, penetrances = penetr, plot = plot)
  # print(res)
  res[[1]][["FLB"]]
}

#------------

test_that("Bayrak's Family 4 is correct", {
  x = nuclearPed(3, sex = c(1,1,2))
  x = addDaughter(addSon(addChildren(x, fa = 3, nch = 2, sex = 2), 4), 5)
  bf = quickFLB(x, aff = c(3:5,7,8,10,12), carriers = c(4,7,8,10,12), proband = 4, pl=T, unk=1:2)
  # bf
  expect_equal(63.63, round(bf, 2))
})

test_that("Bayrak's Family 6 is correct", {
  x = nuclearPed(4, sex = c(1,1,1,2))
  x = addDaughter(addSon(x,4), 6)
  bf = quickFLB(x, aff = c(1,3,4,5,6,8,10), carriers = c(1,3,5,8,10), proband = 8, pl=T)
  # bf
  expect_equal(31.82, round(bf, 2))
})


test_that("Bayrak's Family 7 is correct", {
  x = swapSex(halfSibPed(nch1 = 1, nch2 = 2, sex1 = 2, sex2 = 2), 1)
  bf = quickFLB(x, aff = c(1,4,5,6), carriers = 4:6, proband = 4, pl=T)
  # bf
  expect_equal(7.98, round(bf, 2))
})

test_that("Bayrak's Family 8 is correct", {
  x = nuclearPed(3, sex = c(2,1,1))
  bf = quickFLB(x, aff = 2:5, carriers = 3:5, proband = 3, pl=T)
  # bf
  expect_equal(3.99, round(bf, 2))
})
