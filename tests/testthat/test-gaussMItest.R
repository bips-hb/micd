test_that("gaussMItest equals gaussCItestMI", {
  
  n <- sample(200:500, 1)
  daten <- rnorm(5 * n)
  daten[sample(5 * n, trunc(5 * n / 10))] <- NA
  daten <- as.data.frame(matrix(daten, ncol = 5))
  
  imp <- mice(daten, method = "norm")
  suffMI <- getSuff(imp, test = "gaussMItest")
  
  gaussMItest(1, 2, c(4,5), suffStat = suffMI)
  gaussCItestMI(1, 2, c(4,5), data = imp)

  expect_equal(gaussMItest(1, 2, c(4,5), suffStat = suffMI), 
               gaussCItestMI(1, 2, c(4,5), data = imp))
})
