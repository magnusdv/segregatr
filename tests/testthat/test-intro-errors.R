
nuc = nuclearPed(1)
cous = cousinPed(1, child = TRUE)

test_that("FLB() catches errors", {

  # Loops
  expect_error(FLB(breakLoops(cous, verbose = FALSE), proband = 9, aff = 9, carriers = 9, freq = 0.001),
               "Pedigrees with pre-broken loops are not allowed", fixed = TRUE)
  #expect_error(FLB(cous, proband = 9, aff = 9, carriers = 9, freq = 0.001, liability = rep(1,9)),
  #             "Liability classes are not yet implemented when `x` has loops")

  # Proband
  expect_error(FLB(nuc, proband = ''), "A proband must be specified")
  expect_error(FLB(nuc, proband = 0), "Unknown proband: 0")
  expect_error(FLB(nuc, proband = 1, aff = 2, freq = 0.01), "The proband must be affected")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 2, freq = 0.01), "The proband must be a carrier")
  expect_error(FLB(nuc, proband = 1, aff = 1, homozygous = 2, freq = 0.01), "The proband must be a carrier")

  # Inconsistent input vectors
  expect_error(FLB(nuc, proband = 1, aff = 1, unknown = 1, carriers = 1),
               "Individual specified as both affected and unknown: 1")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, noncarriers = 1),
               "Individual specified as both a carrier and a non-carrier: 1")
  expect_error(FLB(nuc, proband = 1, aff = 1, homozygous = 1, noncarriers = 1),
               "Individual specified as both homozygous and a non-carrier: 1")
  expect_error(FLB(nuc, proband = 1, aff = 1, homozygous = 1, Xchrom = TRUE),
               "Male individual specified as a homozygous carrier in an X-linked inheritance model: 1")

  # Marker allele frequency
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1),
               "An allele frequency must be specified")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = '0.1'),
               "The allele frequency must be a single number strictly between 0 and 1")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = c(0.1, 0.9)),
               "The allele frequency must be a single number strictly between 0 and 1")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 2),
               "The allele frequency must be a single number strictly between 0 and 1")

  expect_error(FLB(cous, proband = 9, aff = 9, carriers = 9, freq = 0.001),
               "No penetrance values given")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001),
               "No penetrance values given")

  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = "AR", Xchrom = TRUE),
               "For X-linked models, `penetrances` must be a list with elements `male` and `female`")
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = list(m=c(0,1), f=c(0,1,1)), Xchrom = TRUE),
               "For X-linked models, `penetrances` must be a list with elements `male` and `female`")

  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penet = "AR", liability = c(1,1)),
               "Length of `liability` vector must be 1 or 3", fixed = TRUE)
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penet = "AR", liability = 1:2),
               "Length of `liability` vector must be 1 or 3", fixed = TRUE)
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penet = "AR", liability = c("1"=1, "2"=2, "4"=3)),
               "Unknown ID in names of liability vector: 4", fixed = TRUE)


  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = "XR", liability = c(1,1,2), Xchrom = TRUE),
               "Illegal liability class (males): 2", fixed = TRUE)
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = "XR", liability = c(1,3,1), Xchrom = TRUE),
               "Illegal liability class (females): 3", fixed = TRUE)
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = "XR", liability = c(2,2,2), Xchrom = TRUE),
               "Illegal liability class (males): 2", fixed = TRUE)
  expect_error(FLB(nuc, proband = 1, aff = 1, carriers = 1, freq = 0.001, penetrances = "AR", liability = c(1,1,2)),
               "Illegal liability class: 2", fixed = TRUE)
})
