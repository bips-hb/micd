# ----------------------------------------------------------------------------------------
# Compare test results
# ----------------------------------------------------------------------------------------
set.seed(346)
dat <- data.frame(x = rnorm(100), 
                  A = factor(rbinom(100, 1, 0.4)), 
                  y = rnorm(100), 
                  B = factor(rbinom(100, 2, c(0.1,0.4))),
                  z = rnorm(100), 
                  C = factor(rbinom(100, 2, c(0.2,0.3))),
                  Z = rnorm(100), 
                  D = factor(rbinom(100, 1, 0.5)))

dat.mis <- dat[-sample(100,10),]

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
rm(res1, res2, res3)


# 2 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- NULL

#res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
suffstat <- getSuff(dat[,c(x,y,S)], test = "disCItest", adaptDF = TRUE, nlev = c(nlev_x, nlev_y))
res3 <- pcalg::disCItest(1,2,NULL, suffStat = suffstat)

(eq2 <- all.equal(res1, res2, res3, tol = 1e-5))
test2 <- isTRUE(eq2)
rm(res1, res2, res3)


# 1 discrete, 1 continuous ----
x <- 3; 
y <- 6; nlev_y <- nlevels(dat[,y])
S <- NULL

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq3 <- all.equal(res1, res2, tol = 1e-5))
test3 <- isTRUE(eq3)
rm(res1, res2)










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
rm(res1, res2, res3)


# 3 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- 6; nlev_S <- nlevels(dat[,S])

#res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, S, suffStat = dat)
suffstat <- getSuff(dat[,c(x,y,S)], test = "disCItest", adaptDF = TRUE, 
                    nlev = c(nlev_x, nlev_y, nlev_S))
res3 <- pcalg::disCItest(1,2,NULL, suffStat = suffstat)

(eq5 <- all.equal(res1, res2, res3, tol = 1e-5))
test5 <- isTRUE(eq5)
rm(res1, res2, res3)


# 1 discrete, 2 continuous ----
x <- 1
y <- 3
S <- 2; nlev_S <- nlevels(dat[,S])

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq6 <- all.equal(res1, res2, tol = 1e-5))
test6 <- isTRUE(eq6)
rm(res1, res2)






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
rm(res1, res2, res3)


# 4 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- c(6,8); nlev_S <- c(nlevels(dat[,S[1]]), nlevels(dat[,S[2]]))

#res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
#res2 <- mixCItest(x, y, S, suffStat = dat)
suffstat <- getSuff(dat[,c(x,y,S)], test = "disCItest", adaptDF = TRUE, 
                    nlev = c(nlev_x, nlev_y, nlev_S))
res3 <- pcalg::disCItest(1,2,NULL, suffStat = suffstat)

(eq8 <- all.equal(res1, res2, res3, tol = 1e-5))
test8 <- isTRUE(eq8)
rm(res1, res2, res3)


# 2 discrete, 2 continuous ----
x <- 1
y <- 3
S <- c(7,8); nlev_S <- nlevels(dat[,8])

res1 <- flexCItest(x, y, NULL, suffStat = getSuff(dat, test="flexCItest"))
res2 <- mixCItest(x, y, NULL, suffStat = dat)

(eq9 <- all.equal(res1, res2, tol = 1e-5))
test9 <- isTRUE(eq9)
rm(res1, res2)




# Result --------------------------------------------------------------------------------
if (!all(test1,test2,test3,test4,test5,test6,test7,test8,test9))
  stop("CI tests wrong!")

















# ----------------------------------------------------------------
pcalg::gaussCItest(1, 2, 3, suffStat = getSuff(dat[,c(1,3,5)], test = "gaussCItest"))
mixCItest(1, 3, 5, suffStat = dat)
flexCItest(1, 3, 5, suffStat = getSuff(dat, test="flexCItest"))



suffstat <- getSuff(dat[,c(2,4,6)], test = "disCItest", adaptDF = TRUE, nlev = c(2,3,3))
pcalg::disCItest(1,2,3, suffStat = suffstat)
mixCItest(2, 4, 6, suffStat = dat)
flexCItest(2, 4, 6, suffStat = getSuff(dat, test="flexCItest")) ##!!!



imp  <- mice(dat.mis, method = "norm")
imp1 <- mice(dat.mis[,c(1,3,5)], method = "norm")
imp2 <- mice(dat.mis[,c(2,4,6)], method = "norm")


gaussMItest(1,2,3, suffStat = getSuff(imp1, test = 'gaussMItest'))  #!
mixMItest(1, 3,5, suffStat =  complete(imp, action = "all"))       #!
mixMItest(1, 3,5, suffStat =  getSuff(imp, test = "mixMItest"))    #!
flexMItest(1,3,5,  suffStat = getSuff(imp, test="flexMItest"))


disMItest(1,2,3, suffStat = getSuff(imp2, test = "disMItest"))
mixMItest(2,4,6, suffStat = complete(imp, action = "all"))
mixMItest(2,4,6, suffStat = getSuff(imp, test = "mixMItest"))
flexMItest(2,4,6,suffStat = getSuff(imp, test="flexMItest"))














gaussMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="gaussMItest"))
  mixMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="mixMItest"))
  mixMItest(1, 2, c(4,5), suffStat = mice::complete(imp2, action="all"))
flexMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="flexMItest"))
