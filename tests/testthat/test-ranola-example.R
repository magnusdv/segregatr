# Validate the (benign) example of Ranola et al (2018)

PENET = read.table(header = T, sep = "\t", as.is = T, comment.char = "", text =
"f0	f1	comment
7.6e-08	1.1e-05	# 1: Male [0-20)
1.2e-06	0.00017	# 2: Male [20-30)
1.9e-05	0.0012	# 3: Male [30-40)
8.5e-05	0.003	# 4: Male [40-50)
0.00027	0.0062	# 5: Male [50-60)
0.00067	0.012	# 6: Male [60-70)
0.0012	0.018	# 7: Male [70-)
8.9e-07	0.001026	# 8: Female [0-20)
4.0997e-05	0.047524	# 9: Female [20-30)
0.00189916	0.18042	# 10: Female [30-40)
0.00878848	0.3736	# 11: Female [40-50)
0.0275136	0.5752	# 12: Female [50-60)
0.05646	0.6889	# 13: Female [60-70)
0.0793	0.785	# 14: Female [70-)")

PENET$f2 = PENET$f1
PENET = PENET[c("f0", "f1", "f2", "comment")]

test_that("Ranola's example pedigree evaluates correctly", {
  pedfile = system.file("extdata", "RanolaPedigree549.ped", package = "segregatr")
  x = readPed(pedfile)
  # plot(x)

  age = c(75, 60, 77, 70, 40, 70, 74, 57, 52, 53, 49, 49, 68, 43, 32, 66, 42, 38, 38, 80, 57, 57)
  liab = pmin(age %/% 10, 7) + 7*(getSex(x)-1)

  bf = FLB(x,
          carriers = c(3:5,11,14:15,17:18,22), #,2,6)
          noncarriers = c(7:10,12,13,16,19,20:21),
          affected = 5:6,
          proband = 5,
          penetrances = PENET[1:3],
          liability = liab,
          freq = 0.001,
          plot = F)

  expect_equal(0.43, round(bf, 2))
})
