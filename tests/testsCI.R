# ----------------------------------------------------------------------------------------
# Compare test results for not imputed data
# ----------------------------------------------------------------------------------------
library(micd)

dat <- data.frame(x = rnorm(200), 
                  A = factor(rbinom(200, 1, 0.4)), 
                  y = rnorm(200), 
                  B = factor(rbinom(200, 2, c(0.1,0.4))),
                  z = rnorm(200), 
                  C = factor(rbinom(200, 2, c(0.2,0.3))),
                  Z = rnorm(200), 
                  D = factor(rbinom(200, 1, 0.5)))

# 2 variables ----------------------------------------------------------------------------
# 2 continuous ----------------
x <- 3
y <- 5
S <- NULL

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::gaussCItest(1, 2, NULL, suffStat = getSuff(dat[,c(x, y, S)], test = "gaussCItest"))

(eq1 <- all.equal(res1, res2, res3, tol = 1e-5))
test1 <- isTRUE(eq1)
rm(res1, res2, res3, x, y, S)


# 2 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- NULL

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::disCItest(1,2,NULL, 
                         suffStat = getSuff(dat[,c(x,y,S)], test = "disCItest", 
                                            adaptDF = TRUE, nlev = c(nlev_x, nlev_y)))

(eq2 <- all.equal(res1, res2, res3, tol = 1e-5))
test2 <- isTRUE(eq2)
rm(res1, res2, res3, x, y, S, nlev_x, nlev_y)


# 1 discrete, 1 continuous ----
x <- 3; 
y <- 6; nlev_y <- nlevels(dat[,y])
S <- NULL

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq3 <- all.equal(res1, res2, tol = 1e-5))
test3 <- isTRUE(eq3)
rm(res1, res2, x, y, S, nlev_y)




# 3 variables ----------------------------------------------------------------------------
# 3 continuous ----------------
x <- 3
y <- 5
S <- 1

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::gaussCItest(1, 2, NULL, suffStat = getSuff(dat[,c(x, y, S)], test = "gaussCItest"))

(eq4 <- all.equal(res1, res2, res3, tol = 1e-5))
test4 <- isTRUE(eq4)
rm(res1, res2, res3, x, y, S)


# 3 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- 6; nlev_S <- nlevels(dat[,S])

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::disCItest(1,2,3, 
                         suffStat =  getSuff(dat[,c(x,y,S)], test = "disCItest", 
                                     adaptDF = TRUE, nlev = c(nlev_x, nlev_y, nlev_S)))

(eq5 <- all.equal(res1, res2, res3, tol = 1e-5))
test5 <- isTRUE(eq5)
rm(res1, res2, res3, x, y, S, nlev_x, nlev_y, nlev_S)


# 1 discrete, 2 continuous ----
x <- 1
y <- 3
S <- 2; nlev_S <- nlevels(dat[,S])

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq6 <- all.equal(res1, res2, tol = 1e-5))
test6 <- isTRUE(eq6)
rm(res1, res2, x, y, S, nlev_S)






# more variables -------------------------------------------------------------------------
# 4 continuous ----------------
x <- 3
y <- 5
S <- c(1,7)

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::gaussCItest(1, 2, NULL, suffStat = getSuff(dat[,c(x, y, S)], test = "gaussCItest"))

(eq7 <- all.equal(res1, res2, res3, tol = 1e-5))
test7 <- isTRUE(eq7)
rm(res1, res2, res3, x, y, S)


# 4 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- c(6,8); nlev_S <- c(nlevels(dat[,S[1]]), nlevels(dat[,S[2]]))

res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
res3 <- pcalg::disCItest(1, 2, c(3,4), 
                         suffStat = getSuff(dat[,c(x,y,S)], test = "disCItest", 
                                    adaptDF = TRUE, nlev = c(nlev_x, nlev_y, nlev_S)))

(eq8 <- all.equal(res1, res2, res3, tol = 1e-5))
test8 <- isTRUE(eq8)
rm(res1, res2, res3, x, y, S, nlev_x, nlev_y, nlev_S)


# 2 discrete, 2 continuous ----
x <- 1
y <- 3
S <- c(7,8); nlev_S <- nlevels(dat[,8])

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq9 <- all.equal(res1, res2, tol = 1e-5))
test9 <- isTRUE(eq9)
rm(res1, res2, x, y, S, nlev_S)




# Result --------------------------------------------------------------------------------
if (!all(test1,test2,test3,test4,test5,test6,test7,test8,test9))
  stop("CI tests wrong!")
