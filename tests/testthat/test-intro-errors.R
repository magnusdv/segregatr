test_that("FLB() catches errors", {

  expect_error(FLB(breakLoops(cousinPed(1, child = TRUE), verbose = FALSE)), "Pedigrees with pre-broken loops are not allowed as input to `FLB()`", fixed = TRUE)

  expect_error(FLB(nuclearPed(1), proband = ''), "A proband must be specified")
  expect_error(FLB(nuclearPed(1), proband = 0), "Unknown proband: 0")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 2), "The proband must be affected")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 2), "The proband must be a carrier")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, homozygous = 2), "The proband must be a carrier")

  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, unknown = 1, carriers = 1), "Individual specified as both affected and unknown: 1")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, noncarriers = 1), "Individual specified as both a carrier and a non-carrier: 1")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, homozygous = 1, noncarriers = 1), "Individual specified as both a carrier and a non-carrier: 1")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, homozygous = 1, Xchrom = TRUE), "Male individual specified as a homozygous carrier in an X-linked inheritance model: 1")

  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1), "An allele frequency must be specified")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = '0.1'), "The allele frequency must be a single number strictly between 0 and 1")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = c(0.1, 0.9)), "The allele frequency must be a single number strictly between 0 and 1")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 2), "The allele frequency must be a single number strictly between 0 and 1")

  expect_error(FLB(cousinPed(1, child = TRUE), proband = 9, aff = 9, carriers = 9, freq = 0.001, liability = rep(1,9)), "Liability classes are not yet implemented when `x` has loops")
  expect_error(FLB(cousinPed(1, child = TRUE), proband = 9, aff = 9, carriers = 9, freq = 0.001), "No penetrance values given")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001), "No penetrance values given")

  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = c(0, 0, 1), Xchrom = TRUE), "For X-linked models, `penetrances` must be a list with elements `male` and `female`")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = list(c(0, 1), c(0, 0, 1)), Xchrom = TRUE), "For X-linked models, `penetrances` must be a list with elements `male` and `female`")

  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = c(0, 0, 1), liability = c(1,1)), "Pedigree size (3) and assigned liability classes (2) must be equal", fixed = TRUE)
  expect_error(FLB(nuclearPed(3), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = c(0, 0, 1), liability = 1), "Pedigree size (5) and assigned liability classes (1) must be equal", fixed = TRUE)

  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = list(male = c(0, 1), female = c(0, 0, 1)), liability = c(1,1,2), Xchrom = TRUE), "Illegal liability class: male 2")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = list(male = c(0, 1), female = c(0, 0, 1)), liability = c(1,3,1), Xchrom = TRUE), "Illegal liability class: female 3")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = list(male = c(0, 1), female = c(0, 0, 1)), liability = c(2,2,2), Xchrom = TRUE), "Illegal liability class: male 2; female 2")
  expect_error(FLB(nuclearPed(1), proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = c(0, 0, 1), liability = c(1,1,2)), "Illegal liability class: 2")

})
