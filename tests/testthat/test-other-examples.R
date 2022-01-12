test_that("Results for FLB() Example 1 are correct", {

  # pedigree
  x = nuclearPed(2)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = data.frame(f0 = c(0, 0.01, 0.05),
                 f1 = c(1, 0.99, 0.95),
                 f2 = c(1, 0.99, 0.95))
  l = rep(1, 4)

  # FLB call
  FLB_res = FLB(x, carriers = 3:4, aff = 3:4, unknown = 1:2,
                freq = b, penetrances = f, liability = l, proband = 3, details = TRUE)

  # calculations
  numer1 =
    1*a^2 * 1*a^2 * f[l[3],2]*0 * f[l[4],2]*0 + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],2]*1 * f[l[4],2]*1 + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],2]*1 * f[l[4],2]*1 + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],2]*0 * f[l[4],2]*0 + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5   # 1=a/b, 2=a/b
  likM =
    a^2 * a^2 * 0 * 0 + # 1=a/a, 2=a/a
    a^2 * b^2 * 1 * 1 + # 1=a/a, 2=b/b
    a^2 * 2*a*b * 0.5 * 0.5 + # 1=a/a, 2=a/b
    b^2 * a^2 * 1 * 1 + # 1=b/b, 2=a/a
    b^2 * b^2 * 0 * 0 + # 1=b/b, 2=b/b
    b^2 * 2*a*b * 0.5 * 0.5 + # 1=b/b, 2=a/b
    2*a*b * a^2 * 0.5 * 0.5 + # 1=a/b, 2=a/a
    2*a*b * b^2 * 0.5 * 0.5 + # 1=a/b, 2=b/b
    2*a*b * 2*a*b * 0.5 * 0.5   # 1=a/b, 2=a/b
  likDis =
    1*a^2 * 1*a^2 * sum(f[l[3],]*c(1,0,0)) * sum(f[l[4],]*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * sum(f[l[3],]*c(0,1,0)) * sum(f[l[4],]*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * sum(f[l[3],]*c(0.5,0.5,0)) * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * sum(f[l[3],]*c(0,1,0)) * sum(f[l[4],]*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * sum(f[l[3],]*c(0,0,1)) * sum(f[l[4],]*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * sum(f[l[3],]*c(0,0.5,0.5)) * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * sum(f[l[3],]*c(0.5,0.5,0)) * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * sum(f[l[3],]*c(0,0.5,0.5)) * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * sum(f[l[3],]*c(0.25,0.5,0.25)) * sum(f[l[4],]*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  denom2 =
    1*a^2 * 1*a^2 * f[l[3],2]*0 * sum(f[l[4],]*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],2]*1 * sum(f[l[4],]*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],2]*1 * sum(f[l[4],]*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],2]*0 * sum(f[l[4],]*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],2]*0.5 * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  likMproband = 2*a*b

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("Results for FLB() Example 1 using different liability classes are correct", {

  # pedigree
  x = nuclearPed(2)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = data.frame(f0 = c(0, 0.01, 0.05),
                 f1 = c(1, 0.99, 0.95),
                 f2 = c(1, 0.99, 0.95))
  l = c(1,3,2,2)

  # FLB call
  FLB_res = FLB(x, carriers = 3:4, aff = 3:4, unknown = 1:2,
                freq = b, penetrances = f, liability = l, proband = 3, details = TRUE)

  # calculations
  numer1 =
    1*a^2 * 1*a^2 * f[l[3],2]*0 * f[l[4],2]*0 + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],2]*1 * f[l[4],2]*1 + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],2]*1 * f[l[4],2]*1 + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],2]*0 * f[l[4],2]*0 + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],2]*0.5 * f[l[4],2]*0.5 + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],2]*0.5 * f[l[4],2]*0.5   # 1=a/b, 2=a/b
  likM =
    a^2 * a^2 * 0 * 0 + # 1=a/a, 2=a/a
    a^2 * b^2 * 1 * 1 + # 1=a/a, 2=b/b
    a^2 * 2*a*b * 0.5 * 0.5 + # 1=a/a, 2=a/b
    b^2 * a^2 * 1 * 1 + # 1=b/b, 2=a/a
    b^2 * b^2 * 0 * 0 + # 1=b/b, 2=b/b
    b^2 * 2*a*b * 0.5 * 0.5 + # 1=b/b, 2=a/b
    2*a*b * a^2 * 0.5 * 0.5 + # 1=a/b, 2=a/a
    2*a*b * b^2 * 0.5 * 0.5 + # 1=a/b, 2=b/b
    2*a*b * 2*a*b * 0.5 * 0.5   # 1=a/b, 2=a/b
  likDis =
    1*a^2 * 1*a^2 * sum(f[l[3],]*c(1,0,0)) * sum(f[l[4],]*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * sum(f[l[3],]*c(0,1,0)) * sum(f[l[4],]*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * sum(f[l[3],]*c(0.5,0.5,0)) * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * sum(f[l[3],]*c(0,1,0)) * sum(f[l[4],]*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * sum(f[l[3],]*c(0,0,1)) * sum(f[l[4],]*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * sum(f[l[3],]*c(0,0.5,0.5)) * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * sum(f[l[3],]*c(0.5,0.5,0)) * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * sum(f[l[3],]*c(0,0.5,0.5)) * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * sum(f[l[3],]*c(0.25,0.5,0.25)) * sum(f[l[4],]*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  denom2 =
    1*a^2 * 1*a^2 * f[l[3],2]*0 * sum(f[l[4],]*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],2]*1 * sum(f[l[4],]*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],2]*1 * sum(f[l[4],]*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],2]*0 * sum(f[l[4],]*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],2]*0.5 * sum(f[l[4],]*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],2]*0.5 * sum(f[l[4],]*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  likMproband = 2*a*b

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("Results for FLB() Example 2 are correct", {

  # pedigree
  x = nuclearPed(4)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = data.frame(f0 = c(0, 0.01, 0.05),
                 f1 = c(0, 0.01, 0.05),
                 f2 = c(1, 0.99, 0.95))
  l = rep(2, 6)

  # FLB call
  FLB_res = FLB(x, carriers = 4:5, homozygous = 3, noncarriers = 6, aff = 3, unknown = 1:2,
                freq = b, penetrances = f, liability = l, proband = 3, details = TRUE)

  # calculations
  numer1 =
    1*a^2 * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*0 * (1-f[l[5],2])*0 * (1-f[l[6],1])*1 + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],3]*0 * (1-f[l[4],2])*1 * (1-f[l[5],2])*1 * (1-f[l[6],1])*0 + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],3]*0 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.5 + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*1 * (1-f[l[5],2])*1 * (1-f[l[6],1])*0 + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],3]*1 * (1-f[l[4],2])*0 * (1-f[l[5],2])*0 * (1-f[l[6],1])*0 + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],3]*0.5 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0 + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.5 + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],3]*0.5 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0 + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],3]*0.25 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.25 # 1=a/b, 2=a/b
  likM =
    a^2 * a^2 * 0 * 0 * 0 * 1 + # 1=a/a, 2=a/a
    a^2 * b^2 * 0 * 1 * 1 * 0 + # 1=a/a, 2=b/b
    a^2 * 2*a*b * 0 * 0.5 * 0.5 * 0.5 + # 1=a/a, 2=a/b
    b^2 * a^2 * 0 * 1 * 1 * 0 + # 1=b/b, 2=a/a
    b^2 * b^2 * 1 * 0 * 0 * 0 + # 1=b/b, 2=b/b
    b^2 * 2*a*b * 0.5 * 0.5 * 0.5 * 0 + # 1=b/b, 2=a/b
    2*a*b * a^2 * 0 * 0.5 * 0.5 * 0.5 + # 1=a/b, 2=a/a
    2*a*b * b^2 * 0.5 * 0.5 * 0.5 * 0 + # 1=a/b, 2=b/b
    2*a*b * 2*a*b * 0.25 * 0.5 * 0.5 * 0.25 # 1=a/b, 2=a/b
  likDis =
    1*a^2 * 1*a^2 * sum(f[l[3],]*c(1,0,0)) * sum((1-f[l[4],])*c(1,0,0)) * sum((1-f[l[5],])*c(1,0,0)) * sum((1-f[l[6],])*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * sum(f[l[3],]*c(0,1,0)) * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * sum(f[l[3],]*c(0.5,0.5,0)) * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * sum(f[l[3],]*c(0,1,0)) * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * sum(f[l[3],]*c(0,0,1)) * sum((1-f[l[4],])*c(0,0,1)) * sum((1-f[l[5],])*c(0,0,1)) * sum((1-f[l[6],])*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * sum(f[l[3],]*c(0,0.5,0.5)) * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * sum(f[l[3],]*c(0.5,0.5,0)) * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * sum(f[l[3],]*c(0,0.5,0.5)) * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * sum(f[l[3],]*c(0.25,0.5,0.25)) * sum((1-f[l[4],])*c(0.25,0.5,0.25)) * sum((1-f[l[5],])*c(0.25,0.5,0.25)) * sum((1-f[l[6],])*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  denom2 =
    1*a^2 * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(1,0,0)) * sum((1-f[l[5],])*c(1,0,0)) * sum((1-f[l[6],])*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],3]*0 * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],3]*1 * sum((1-f[l[4],])*c(0,0,1)) * sum((1-f[l[5],])*c(0,0,1)) * sum((1-f[l[6],])*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],3]*0.5 * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0,5,0)) * sum((1-f[l[6],])*c(0.5,0,5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],3]*0.5 * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],3]*0.25 * sum((1-f[l[4],])*c(0.25,0.5,0.25)) * sum((1-f[l[5],])*c(0.25,0.5,0.25)) * sum((1-f[l[6],])*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
   likMproband = b^2

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("Results for FLB() Example 2 using different liability classes are correct", {

  # pedigree
  x = nuclearPed(4)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = data.frame(f0 = c(0, 0.01, 0.05),
                 f1 = c(0, 0.01, 0.05),
                 f2 = c(1, 0.99, 0.95))
  l = c(1,1,2,3,2,3)

  # FLB call
  FLB_res = FLB(x, carriers = 4:5, homozygous = 3, noncarriers = 6, aff = 3, unknown = 1:2,
                freq = b, penetrances = f, liability = l, proband = 3, details = TRUE)

  # calculations
  numer1 =
    1*a^2 * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*0 * (1-f[l[5],2])*0 * (1-f[l[6],1])*1 + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],3]*0 * (1-f[l[4],2])*1 * (1-f[l[5],2])*1 * (1-f[l[6],1])*0 + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],3]*0 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.5 + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*1 * (1-f[l[5],2])*1 * (1-f[l[6],1])*0 + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],3]*1 * (1-f[l[4],2])*0 * (1-f[l[5],2])*0 * (1-f[l[6],1])*0 + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],3]*0.5 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0 + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],3]*0 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.5 + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],3]*0.5 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0 + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],3]*0.25 * (1-f[l[4],2])*0.5 * (1-f[l[5],2])*0.5 * (1-f[l[6],1])*0.25 # 1=a/b, 2=a/b
  likM =
    a^2 * a^2 * 0 * 0 * 0 * 1 + # 1=a/a, 2=a/a
    a^2 * b^2 * 0 * 1 * 1 * 0 + # 1=a/a, 2=b/b
    a^2 * 2*a*b * 0 * 0.5 * 0.5 * 0.5 + # 1=a/a, 2=a/b
    b^2 * a^2 * 0 * 1 * 1 * 0 + # 1=b/b, 2=a/a
    b^2 * b^2 * 1 * 0 * 0 * 0 + # 1=b/b, 2=b/b
    b^2 * 2*a*b * 0.5 * 0.5 * 0.5 * 0 + # 1=b/b, 2=a/b
    2*a*b * a^2 * 0 * 0.5 * 0.5 * 0.5 + # 1=a/b, 2=a/a
    2*a*b * b^2 * 0.5 * 0.5 * 0.5 * 0 + # 1=a/b, 2=b/b
    2*a*b * 2*a*b * 0.25 * 0.5 * 0.5 * 0.25 # 1=a/b, 2=a/b
  likDis =
    1*a^2 * 1*a^2 * sum(f[l[3],]*c(1,0,0)) * sum((1-f[l[4],])*c(1,0,0)) * sum((1-f[l[5],])*c(1,0,0)) * sum((1-f[l[6],])*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * sum(f[l[3],]*c(0,1,0)) * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * sum(f[l[3],]*c(0.5,0.5,0)) * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * sum(f[l[3],]*c(0,1,0)) * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * sum(f[l[3],]*c(0,0,1)) * sum((1-f[l[4],])*c(0,0,1)) * sum((1-f[l[5],])*c(0,0,1)) * sum((1-f[l[6],])*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * sum(f[l[3],]*c(0,0.5,0.5)) * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * sum(f[l[3],]*c(0.5,0.5,0)) * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * sum(f[l[3],]*c(0,0.5,0.5)) * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * sum(f[l[3],]*c(0.25,0.5,0.25)) * sum((1-f[l[4],])*c(0.25,0.5,0.25)) * sum((1-f[l[5],])*c(0.25,0.5,0.25)) * sum((1-f[l[6],])*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  denom2 =
    1*a^2 * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(1,0,0)) * sum((1-f[l[5],])*c(1,0,0)) * sum((1-f[l[6],])*c(1,0,0)) + # 1=a/a, 2=a/a
    1*a^2 * 1*b^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=a/a, 2=b/b
    1*a^2 * 1*2*a*b * f[l[3],3]*0 * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0.5,0)) * sum((1-f[l[6],])*c(0.5,0.5,0)) + # 1=a/a, 2=a/b
    1*b^2 * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0,1,0)) * sum((1-f[l[5],])*c(0,1,0)) * sum((1-f[l[6],])*c(0,1,0)) + # 1=b/b, 2=a/a
    1*b^2 * 1*b^2 * f[l[3],3]*1 * sum((1-f[l[4],])*c(0,0,1)) * sum((1-f[l[5],])*c(0,0,1)) * sum((1-f[l[6],])*c(0,0,1)) + # 1=b/b, 2=b/b
    1*b^2 * 1*2*a*b * f[l[3],3]*0.5 * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=b/b, 2=a/b
    1*2*a*b * 1*a^2 * f[l[3],3]*0 * sum((1-f[l[4],])*c(0.5,0.5,0)) * sum((1-f[l[5],])*c(0.5,0,5,0)) * sum((1-f[l[6],])*c(0.5,0,5,0)) + # 1=a/b, 2=a/a
    1*2*a*b * 1*b^2 * f[l[3],3]*0.5 * sum((1-f[l[4],])*c(0,0.5,0.5)) * sum((1-f[l[5],])*c(0,0.5,0.5)) * sum((1-f[l[6],])*c(0,0.5,0.5)) + # 1=a/b, 2=b/b
    1*2*a*b * 1*2*a*b * f[l[3],3]*0.25 * sum((1-f[l[4],])*c(0.25,0.5,0.25)) * sum((1-f[l[5],])*c(0.25,0.5,0.25)) * sum((1-f[l[6],])*c(0.25,0.5,0.25)) # 1=a/b, 2=a/b
  likMproband = b^2

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("Results for FLB() Example 3 are correct", {

  # pedigree
  x = nuclearPed(3, sex = c(1, 1, 2))
  x = addChildren(x, mo = 5, sex = 1:2, verbose = FALSE)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = list(male   = data.frame(f0 = c(0, 0.01, 0.05),
                               f1 = c(1, 0.99, 0.95)),
           female = data.frame(f0 = c(0, 0.01, 0.05),
                               f1 = c(0, 0.01, 0.05),
                               f2 = c(1, 0.99, 0.95)))
  l = rep(1, 8)

  # FLB call
  FLB_res = FLB(x, carriers = c(3,7), nonc = 4, aff = c(3,7), unknown = 1:2,
      freq = b, penetrances = f, liability = l,
      proband = 7, Xchrom = TRUE, details = TRUE)

  # calculations (E-S algorithm (1,2,3,4-(5)-6,7,8))
  numer1_5_aa = (1-f$male[l[6],1])*a * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  numer1_5_bb = (1-f$male[l[6],1])*a * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  numer1_5_ab = (1-f$male[l[6],1])*a * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  numer1_5 = c(numer1_5_aa, numer1_5_ab, numer1_5_bb)
  numer1 =
    1*a * 1*a^2 * f$male[l[3],2]*0 * (1-f$male[l[4],1])*1 * sum((1-f$female[l[5],])*numer1_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * f$male[l[3],2]*1 * (1-f$male[l[4],1])*0 * sum((1-f$female[l[5],])*numer1_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * f$male[l[3],2]*0.5 * (1-f$male[l[4],1])*0.5 * sum((1-f$female[l[5],])*numer1_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * f$male[l[3],2]*0 * (1-f$male[l[4],1])*1 * sum((1-f$female[l[5],])*numer1_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * f$male[l[3],2]*1 * (1-f$male[l[4],1])*0 * sum((1-f$female[l[5],])*numer1_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * f$male[l[3],2]*0.5 * (1-f$male[l[4],1])*0.5 * sum((1-f$female[l[5],])*numer1_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likM_5_aa = a * 0 * 1 + b * 0 * 1 # 6=a + 6=b
  likM_5_bb = a * 1 * 1 + b * 1 * 1 # 6=a + 6=b
  likM_5_ab = a * 0.5 * 1 + b * 0.5 * 1 # 6=a + 6=b
  likM_5 = c(likM_5_aa, likM_5_ab, likM_5_bb)
  likM =
    1*a * 1*a^2 * 0 * 1 * sum(likM_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * 1 * 0 * sum(likM_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * 0.5 * 0.5 * sum(likM_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * 0 * 1 * sum(likM_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * 1 * 0 * sum(likM_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * 0.5 * 0.5 * sum(likM_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likDis_5_aa = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(1,0)) * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(1,0)) * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  likDis_5_bb = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(0,1)) * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(0,1)) * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  likDis_5_ab = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(0.5,0.5)) * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(0.5,0.5)) * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  likDis_5 = c(likDis_5_aa, likDis_5_ab, likDis_5_bb)
  likDis =
    1*a * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*likDis_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*likDis_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5)) * sum((1-f$female[l[5],])*likDis_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*likDis_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*likDis_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5))* sum((1-f$female[l[5],])*likDis_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  denom2_5_aa = (1-f$male[l[6],1])*a * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  denom2_5_bb = (1-f$male[l[6],1])*a * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  denom2_5_ab = (1-f$male[l[6],1])*a * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  denom2_5 = c(denom2_5_aa, denom2_5_ab, denom2_5_bb)
  denom2 =
    1*a * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*denom2_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*denom2_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5)) * sum((1-f$female[l[5],])*denom2_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*denom2_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*denom2_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5))* sum((1-f$female[l[5],])*denom2_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likMproband = b

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("Results for FLB() Example 3 using different liability classes are correct", {

  # pedigree
  x = nuclearPed(3, sex = c(1, 1, 2))
  x = addChildren(x, mo = 5, sex = 1:2, verbose = FALSE)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes
  f = list(male   = data.frame(f0 = c(0, 0.01, 0.05),
                               f1 = c(1, 0.99, 0.95)),
           female = data.frame(f0 = c(0, 0.01, 0.05),
                               f1 = c(0, 0.01, 0.05),
                               f2 = c(1, 0.99, 0.95)))
  l = c(1,1,2,2,2,2,3,3)

  # FLB call
  FLB_res = FLB(x, carriers = c(3,7), nonc = 4, aff = c(3,7), unknown = 1:2,
                freq = b, penetrances = f, liability = l,
                proband = 7, Xchrom = TRUE, details = TRUE)

  # calculations (E-S algorithm (1,2,3,4-(5)-6,7,8))
  numer1_5_aa = (1-f$male[l[6],1])*a * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  numer1_5_bb = (1-f$male[l[6],1])*a * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  numer1_5_ab = (1-f$male[l[6],1])*a * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  numer1_5 = c(numer1_5_aa, numer1_5_ab, numer1_5_bb)
  numer1 =
    1*a * 1*a^2 * f$male[l[3],2]*0 * (1-f$male[l[4],1])*1 * sum((1-f$female[l[5],])*numer1_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * f$male[l[3],2]*1 * (1-f$male[l[4],1])*0 * sum((1-f$female[l[5],])*numer1_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * f$male[l[3],2]*0.5 * (1-f$male[l[4],1])*0.5 * sum((1-f$female[l[5],])*numer1_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * f$male[l[3],2]*0 * (1-f$male[l[4],1])*1 * sum((1-f$female[l[5],])*numer1_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * f$male[l[3],2]*1 * (1-f$male[l[4],1])*0 * sum((1-f$female[l[5],])*numer1_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * f$male[l[3],2]*0.5 * (1-f$male[l[4],1])*0.5 * sum((1-f$female[l[5],])*numer1_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likM_5_aa = a * 0 * 1 + b * 0 * 1 # 6=a + 6=b
  likM_5_bb = a * 1 * 1 + b * 1 * 1 # 6=a + 6=b
  likM_5_ab = a * 0.5 * 1 + b * 0.5 * 1 # 6=a + 6=b
  likM_5 = c(likM_5_aa, likM_5_ab, likM_5_bb)
  likM =
    1*a * 1*a^2 * 0 * 1 * sum(likM_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * 1 * 0 * sum(likM_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * 0.5 * 0.5 * sum(likM_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * 0 * 1 * sum(likM_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * 1 * 0 * sum(likM_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * 0.5 * 0.5 * sum(likM_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likDis_5_aa = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(1,0)) * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(1,0)) * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  likDis_5_bb = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(0,1)) * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(0,1)) * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  likDis_5_ab = (1-f$male[l[6],1])*a * sum(f$male[l[7],]*c(0.5,0.5)) * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * sum(f$male[l[7],]*c(0.5,0.5)) * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  likDis_5 = c(likDis_5_aa, likDis_5_ab, likDis_5_bb)
  likDis =
    1*a * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*likDis_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*likDis_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5)) * sum((1-f$female[l[5],])*likDis_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*likDis_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*likDis_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5))* sum((1-f$female[l[5],])*likDis_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  denom2_5_aa = (1-f$male[l[6],1])*a * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(1,0,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0 * sum((1-f$female[l[8],])*c(0,1,0)) # 6=a + 6=b
  denom2_5_bb = (1-f$male[l[6],1])*a * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,1,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*1 * sum((1-f$female[l[8],])*c(0,0,1)) # 6=a + 6=b
  denom2_5_ab = (1-f$male[l[6],1])*a * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0.5,0.5,0)) + (1-f$male[l[6],2])*b * f$male[l[7],2]*0.5 * sum((1-f$female[l[8],])*c(0,0.5,0.5)) # 6=a + 6=b
  denom2_5 = c(denom2_5_aa, denom2_5_ab, denom2_5_bb)
  denom2 =
    1*a * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*denom2_5*c(1,0,0)) + # 1=a, 2=a/a
    1*a * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*denom2_5*c(0,1,0)) + # 1=a, 2=b/b
    1*a * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5)) * sum((1-f$female[l[5],])*denom2_5*c(0.5,0.5,0)) + # 1=a, 2=a/b
    1*b * 1*a^2 * sum(f$male[l[3],]*c(1,0)) * sum((1-f$male[l[4],])*c(1,0)) * sum((1-f$female[l[5],])*denom2_5*c(0,1,0)) + # 1=b, 2=a/a
    1*b * 1*b^2 * sum(f$male[l[3],]*c(0,1)) * sum((1-f$male[l[4],])*c(0,1)) * sum((1-f$female[l[5],])*denom2_5*c(0,0,1)) + # 1=b, 2=b/b
    1*b * 1*2*a*b * sum(f$male[l[3],]*c(0.5,0.5)) * sum((1-f$male[l[4],])*c(0.5,0.5))* sum((1-f$female[l[5],])*denom2_5*c(0,0.5,0.5)) # 1=b, 2=a/b
  likMproband = b

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

})

test_that("FLB() results for some recessive inheritance examples are correct", {

  dfreq = 0.00123

  # AR model: affected siblings
  expect_equal(
    FLB(nuclearPed(2), proband = 3, aff = 3:4, homozygous = 3:4, carriers = 1:2,
        freq = dfreq, penetrances = c(0,0,1)),
    4/(1-dfreq)^2
  )

  # XR model: affected brothers
  expect_equal(
    FLB(nuclearPed(2), proband = 3, aff = 3:4, carriers = 2:4, Xchrom = TRUE,
        freq = dfreq, penetrances = list(male = c(0,1), female = c(0,0,1))),
    2/(1-dfreq)
  )

})

test_that("Some quick checks for pedigrees with loops", {

  # pedigree1
  x = cousinPed(degree = 1, child = TRUE)

  # FLB call
  FLB_res = FLB(x, carriers = c(7,9), nonc = c(4,8), aff = c(3,7,9), unknown = 1,
                freq = 0.001, penetrances = c(0,1,1),
                proband = 9, details = TRUE)

  # calculations
  numer1 = denom2 = likDis =
    1*2*0.999*0.001 * 0.999^2 * 0.5 * 0.999^2 * 0.5 * 0.999^2 * 0.5 * 1 * 0.5 # 1=a/b, 2=a/a, 3=a/b, 5=a/a, 6=a/a
  likM = likelihood(setMarkers(x, marker(x, name = 'm', `4` = 'a/a', `7` = 'a/b', `8` = 'a/a', `9` = 'a/b', afreq = c('b' = 0.001, 'a' = 0.999))))
  likMproband = likelihood(setMarkers(x, marker(x, name = 'm', `9` = 'a/b', afreq = c('b' = 0.001, 'a' = 0.999))))

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))

  # pedigree2
  x = cousinPed(degree = 1, child = TRUE)

  # FLB call
  FLB_res = FLB(x, carriers = c(5,8,9), aff = c(5,9), freq = 0.001,
                penetrances = list(male = c(0, 1), female = c(0, 0, 1)),
                proband = 9, Xchrom = TRUE, details = TRUE)

  # calculations
  numer1 = denom2 = likDis =
    (1-0)*0.999 * (1-0)*2*0.999*0.001 * (1-0)*0.5 * (1-0)*0.999^2 * 1*0.5 * (1-0)*0.999^2 * (1-0)*1 * (1-0)*1 * 1*0.5 + # 1=a, 2=a/b, 3=a, 4=a/a, 6=a/a, 7=a
    (1-0)*0.999 * (1-0)*2*0.999*0.001 * (1-0)*0.5 * (1-0)*0.999^2 * 1*0.5 * (1-0)*2*0.999*0.001 * (1-0)*1 * (1-0)*0.5 * 1*0.5 + # 1=a, 2=a/b, 3=a, 4=a/a, 6=a/b, 7=a
    (1-0)*0.999 * (1-0)*2*0.999*0.001 * (1-0)*0.5 * (1-0)*2*0.999*0.001 * 1*0.5 * (1-0)*0.999^2 * (1-0)*0.5 * (1-0)*1 * 1*0.5 + # 1=a, 2=a/b, 3=a, 4=a/b, 6=a/a, 7=a
    (1-0)*0.999 * (1-0)*2*0.999*0.001 * (1-0)*0.5 * (1-0)*2*0.999*0.001 * 1*0.5 * (1-0)*2*0.999*0.001 * (1-0)*0.5 * (1-0)*0.5 * 1*0.5 # 1=a, 2=a/b, 3=a, 4=a/b, 6=a/b, 7=a
  likM = likelihood(setMarkers(x, marker(x, name = 'm', `5` = 'b', `8` = 'a/b', `9` = 'b', afreq = c('b' = 0.001, 'a' = 0.999), chrom = 'X')))
  likMproband = likelihood(setMarkers(x, marker(x, name = 'm', `9` = 'b', afreq = c('b' = 0.001, 'a' = 0.999), chrom = 'X')))

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likDis, unname(FLB_res[[3]]["likDis"]))
  expect_equal(denom2, unname(FLB_res[[2]]["denom2"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))
  expect_equal((numer1/likM)*(likMproband/denom2), unname(FLB_res[[1]]["FLB"]))


  # pedigree3
  x = cousinPed(degree = 1, child = TRUE)
  x = addChildren(x, fa = 1, mo = 2, sex = 2)

  # allele frequencies
  b = 0.0001
  a = 1 - b

  # penetrances and liability classes (not yet implemented)
  f = data.frame(f0 = c(0, 0.01, 0.05),
                 f1 = c(0, 0.01, 0.05),
                 f2 = c(1, 0.99, 0.95))
  l = rep(3, 10)

  # FLB call
  FLB_res = FLB(x, carriers = c(5,7,8), homozygous = c(9,10), noncarriers = c(4,6),
      aff = c(1,9,10), unknown = NULL, freq = b, penetrances = unlist(f[3,]),
      proband = 9, details = TRUE)

  # calculations
  numer1 =
    f[l[1],2]*2*a*b * (1-f[l[2],2])*2*a*b * (1-f[l[3],2])*0.50 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*0.50 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.25 + # 1=a/b, 2=a/b, 3=a/b
    f[l[1],2]*2*a*b * (1-f[l[2],2])*2*a*b * (1-f[l[3],3])*0.25 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*1.00 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.25 + # 1=a/b, 2=a/b, 3=b/b
    f[l[1],2]*2*a*b * (1-f[l[2],3])*b^2 * (1-f[l[3],2])*0.50 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*0.50 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.50 + # 1=a/b, 2=b/b, 3=a/b
    f[l[1],2]*2*a*b * (1-f[l[2],3])*b^2 * (1-f[l[3],3])*0.50 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*1.00 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.50 + # 1=a/b, 2=b/b, 3=b/b
    f[l[1],3]*b^2 * (1-f[l[2],2])*2*a*b * (1-f[l[3],2])*0.50 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*0.50 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.50 + # 1=b/b, 2=a/b, 3=a/b
    f[l[1],3]*b^2 * (1-f[l[2],2])*2*a*b * (1-f[l[3],3])*0.50 * (1-f[l[4],1])*a^2 * (1-f[l[5],2])*0.50 * (1-f[l[6],1])*a^2 * (1-f[l[7],2])*1.00 * (1-f[l[8],2])*0.50 * f[l[9],3]*0.25 * f[l[10],3]*0.50 # 1=b/b, 2=a/b, 3=b/b
  likM = likelihood(setMarkers(x, marker(x, name = 'm', `4` = 'a/a', `5` = 'a/b', `6` = 'a/a', `7` = 'a/b', `8` = 'a/b', `9` = 'b/b', `10` = 'b/b', afreq = c('b' = b, 'a' = a))))
  likMproband = likelihood(setMarkers(x, marker(x, name = 'm', `9` = 'b/b', afreq = c('b' = b, 'a' = a))))

  # test
  expect_equal(numer1, unname(FLB_res[[2]]["numer1"]))
  expect_equal(likM, unname(FLB_res[[3]]["likM"]))
  expect_equal(likMproband, unname(FLB_res[[3]]["likMproband"]))

})

