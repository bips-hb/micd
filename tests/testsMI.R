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
imp  <- mice(dat.mis, method = "norm")

# 2 variables ----------------------------------------------------------------------------
# 2 continuous ----------------
x <- 3
y <- 5
S <- NULL

#warning: different
imp_c <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- gaussMItest(1,2, NULL, suffStat = getSuff(imp_c, test = 'gaussMItest'))
res2 <- mixMItest(x,y,S,  suffStat =  complete(imp, action = "all"))       
res3 <- mixMItest(x,y,S,  suffStat =  getSuff(imp, test = "mixMItest"))    
res4 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq1 <- all.equal(res1, res2, res3, res4, tol = 1e-5))
test1 <- isTRUE(eq1)
rm(res1, res2, res3, res4)


# 2 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- NULL

#res1 <- flexCItest(x, y, S, suffStat = getSuff(dat, test="flexCItest"))
imp_d <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- disMItest(1,2, NULL, suffStat = getSuff(imp_d, test = "disMItest"))
res2 <- mixMItest(x,y,S, suffStat = complete(imp, action = "all"))
res3 <- mixMItest(x,y,S, suffStat = getSuff(imp, test = "mixMItest"))
res4 <- flexMItest(x,y,S,suffStat = getSuff(imp, test="flexMItest"))

(eq2 <- all.equal(res1, res2, res3, res4, tol = 1e-5))
test2 <- isTRUE(eq2)
rm(res1, res2, res3, res4)


# 1 discrete, 1 continuous ----
x <- 3; 
y <- 6; nlev_y <- nlevels(dat[,y])
S <- NULL

res1 <- mixMItest(x,y,S, suffStat = getSuff(imp, test="mixMItest"))
res2 <- mixMItest(x,y,S, suffStat = mice::complete(imp, action="all"))
res3 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq3 <- all.equal(res1, res2, res3, tol = 1e-5))
test3 <- isTRUE(eq3)
rm(res1, res2)










# 3 variables ----------------------------------------------------------------------------
# 3 continuous ----------------
x <- 1
y <- 5
S <- 3

#warning: Unterchiedliche
imp_c <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- gaussMItest(1,2,3, suffStat = getSuff(imp_c, test = 'gaussMItest'))
res2 <- mixMItest(x,y,S,  suffStat =  complete(imp, action = "all"))       
res3 <- mixMItest(x,y,S,  suffStat =  getSuff(imp, test = "mixMItest"))    
res4 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq4 <- all.equal(res1, res2, res3, res4, tol = 1e-5))
test4 <- isTRUE(eq4)
rm(res1, res2, res3)

rm(res1, res2, res3, res4)


# 3 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- 6; nlev_S <- nlevels(dat[,S])

imp_d <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- disMItest(1,2,3, suffStat = getSuff(imp2, test = "disMItest"))
res2 <- mixMItest(x,y,S, suffStat = complete(imp, action = "all"))
res3 <- mixMItest(x,y,S, suffStat = getSuff(imp, test = "mixMItest"))
res4 <- flexMItest(x,y,S,suffStat = getSuff(imp, test="flexMItest"))

(eq5 <- all.equal(res1,res2, res3, res4, tol = 1e-5))
test5 <- isTRUE(eq5)
rm(res1, res2, res3)


# 1 discrete, 2 continuous ----
x <- 1
y <- 3
S <- 2; nlev_S <- nlevels(dat[,S])

res1 <- mixMItest(x,y,S, suffStat = getSuff(imp, test="mixMItest"))
res2 <- mixMItest(x,y,S, suffStat = mice::complete(imp, action="all"))
res3 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq6 <- all.equal(res1, res2, res3, tol = 1e-5))
test6 <- isTRUE(eq6)
rm(res1, res2, res3)






# more variables -------------------------------------------------------------------------
# 4 continuous ----------------
x <- 3
y <- 5
S <- c(1,7)

#warning: Unterchiedliche
imp_c <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- gaussMItest(1,2,3, suffStat = getSuff(imp_c, test = 'gaussMItest'))
res2 <- mixMItest(x,y,S,  suffStat =  complete(imp, action = "all"))       
res3 <- mixMItest(x,y,S,  suffStat =  getSuff(imp, test = "mixMItest"))    
res4 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq7 <- all.equal(res1, res2, res3, res4, tol = 1e-5))
test7 <- isTRUE(eq7)
rm(res1, res2, res3)


# 4 discrete ------------------
x <- 2; nlev_x <- nlevels(dat[,x])
y <- 4; nlev_y <- nlevels(dat[,y])
S <- c(6,8); nlev_S <- c(nlevels(dat[,S[1]]), nlevels(dat[,S[2]]))

imp_d <- mice(dat.mis[,c(x,y,S)], method = "norm")
res1 <- disMItest(1,2,c(3,4), suffStat = getSuff(imp_d, test = "disMItest"))
res2 <- mixMItest(x,y,S, suffStat = complete(imp, action = "all"))
res3 <- mixMItest(x,y,S, suffStat = getSuff(imp, test = "mixMItest"))
res4 <- flexMItest(x,y,S,suffStat = getSuff(imp, test="flexMItest"))

(eq8 <- all.equal(res1, res2, res3, res4, tol = 1e-5))
test8 <- isTRUE(eq8)
rm(res1, res2, res3, res4)


# 2 discrete, 2 continuous ----
x <- 1
y <- 3
S <- c(7,8); nlev_S <- nlevels(dat[,8])

res1 <- mixMItest(x,y,S, suffStat = getSuff(imp, test="mixMItest"))
res2 <- mixMItest(x,y,S, suffStat = mice::complete(imp, action="all"))
res3 <- flexMItest(x,y,S, suffStat = getSuff(imp, test="flexMItest"))

(eq9 <- all.equal(res1, res2, res3, tol = 1e-5))
test9 <- isTRUE(eq9)
rm(res1, res2)




# Result --------------------------------------------------------------------------------
if (!all(test1,test2,test3,test4,test5,test6,test7,test8,test9))
  stop("CI tests wrong!")
