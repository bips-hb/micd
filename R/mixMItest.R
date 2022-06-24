#' Likelihood Ratio Test for (Conditional) Independence between Mixed Variables
#' after Multiple Imputation
#'
#' A modified version of \code{\link{mixCItest}}, to be used within \code{pcalg::\link[pcalg]{skeleton}},
#' \code{pcalg::\link[pcalg]{pc}} or \code{pcalg::\link[pcalg]{fci}} when multiply imputed data sets are available.
#'
#' @param x,y,S    (integer) position of variable X, Y and set of variables S,
#' respectively, in \code{suffStat}. It is tested whether X and Y are
#' conditionally independent given the subset S of the remaining variables.
#' @param suffStat A list of \code{data.frame}s containing the multiply
#' imputed data sets. Usually obtained from a \code{mice::\link[mice:mids-class]{mids}}
#' object using \code{mice::\link[mice:complete.mids]{complete}} with argument \code{action="all"}.
#' Discrete variables must be coded as factors.
#' @param moreOutput If \code{TRUE}, the test statistic, its main components and
#'  the degrees of freedom are returned in addition to the p-value. Defaults to
#'  \code{FALSE}.
#'
#' @details See \code{\link{mixCItest}} for details on the assumptions of the
#' Conditional Gaussian likelihood ratio test. \code{CGtestMI} applies this test
#' to each \code{data.frame} in \code{suffStat}, then combines the results using
#' the rules in Meng & Rubin (1992).
#' @return A p-value. If \code{moreOutput=TRUE}, the test statistic, its main
#' components and the degrees of freedom are returned as well.
#'
#' @author Janine Witte
#'
#' @references Meng X.-L., Rubin D.B. (1992): Performing likelihood ratio tests
#' with multiply imputed data sets. \emph{Biometrika} 79(1):103-111.
#'
#' @examples
#'
#' ## load data (numeric and factor variables)
#' dat <- toenail2[1:400, ]
#'
#' ## delete some observations
#' set.seed(123)
#' dat[sample(400, 20), 2] <- NA
#' dat[sample(400, 30), 4] <- NA
#'
#' ## impute missing values using random forests
#' imp <- mice(dat, method = "rf")
#'
#' ## analyse data
#' # complete data:
#' mixCItest(2, 3, 5, suffStat = toenail2[1:400, ])
#' # multiple imputation:
#' suffMI <- complete(imp, action = "all")
#' mixMItest(2, 3, 5, suffStat =  suffMI)
#' # test-wise deletion:
#' mixCItwd(2, 3, 5, suffStat = dat)
#' # list-wise deletion:
#' sufflwd <- dat[complete.cases(dat), ]
#' mixCItest(2, 3, 5, suffStat = sufflwd)
#' 
#' ## use mixMItest within pcalg::pc
#' pc.fit <- pc(suffStat =  suffMI, indepTest = mixMItest, alpha = 0.01, p = 5)
#' pc.fit
#'
#' @export




mixMItest <- function(x, y, S = NULL, suffStat, moreOutput=FALSE) {
  # suffStat is a list of completed data sets

  # number of imputations / completed data sets
  M <- length(suffStat)

  # indices of discrete variables
  A_xyS <- which.is(suffStat[[1]][c(x,y,S)], "factor")
  A_yS <- which.is(suffStat[[1]][c(y,S)], "factor")
  A_xS <- which.is(suffStat[[1]][c(x,S)], "factor")
  # indeces of continuous variables
  X_xyS <- which.is(suffStat[[1]][c(x,y,S)], "numeric")
  X_yS <- which.is(suffStat[[1]][c(y,S)], "numeric")
  X_xS <- which.is(suffStat[[1]][c(x,S)], "numeric")
  # number of continuous variables
  k_xyS <- length(X_xyS)
  k_yS <- length(X_yS)
  k_xS <- length(X_xS)

  # degrees of freedom if we had complete data
  df_xyS <- dfCG(suffStat[[1]][c(x,y,S)], A_xyS, k_xyS)
  df_yS <- dfCG(suffStat[[1]][c(y,S)], A_yS, k_yS)
  df_xS <- dfCG(suffStat[[1]][c(x,S)], A_xS, k_xS)

  K <- df_xyS - df_yS - df_xS + 1


  if (!is.null(S)) {
    A_S <- which.is(suffStat[[1]][S], "factor")
    X_S <- which.is(suffStat[[1]][S], "numeric")
    k_S <- length(X_S)
    df_S <- dfCG(suffStat[[1]][S], A_S, k_S)
    K <- df_xyS - df_yS - df_xS + df_S
  }

  if (K <= 0) {K <- 1}

  # in each completed data set: calculate maximal log likelihood and output all relevant parameters

  ### xyS
  M1_xyS <- maxJoint(suffStat[[1]][c(x,y,S)], A_xyS, X_xyS, k_xyS, M)

  av_logL_xyS <- M1_xyS[[1]]
  av_Multinomial_xyS <- M1_xyS[[2]]
  av_Mu_xyS <- M1_xyS[[3]]
  av_Sigma_xyS <- M1_xyS[[4]]

  ### yS
  M1_yS <- maxJoint(suffStat[[1]][c(y,S)], A_yS, X_yS, k_yS, M)

  av_logL_yS <- M1_yS[[1]]
  av_Multinomial_yS <- M1_yS[[2]]
  av_Mu_yS <- M1_yS[[3]]
  av_Sigma_yS <- M1_yS[[4]]

  ### xS
  M1_xS <- maxJoint(suffStat[[1]][c(x,S)], A_xS, X_xS, k_xS, M)

  av_logL_xS <- M1_xS[[1]]
  av_Multinomial_xS <- M1_xS[[2]]
  av_Mu_xS <- M1_xS[[3]]
  av_Sigma_xS <- M1_xS[[4]]

  ### S
  if (!is.null(S)) {
    M1_S <- maxJoint(suffStat[[1]][c(S)], A_S, X_S, k_S, M)

    av_logL_S <- M1_S[[1]]
    av_Multinomial_S <- M1_S[[2]]
    av_Mu_S <- M1_S[[3]]
    av_Sigma_S <- M1_S[[4]]
  }
  for (i in 2:M) {
    Mi_xyS <- maxJoint(suffStat[[i]][c(x,y,S)], A_xyS, X_xyS, k_xyS, M)

    av_logL_xyS <- sum(av_logL_xyS, Mi_xyS[[1]])
    av_Multinomial_xyS <- listsum(av_Multinomial_xyS, Mi_xyS[[2]])
    av_Mu_xyS <- listsum(av_Mu_xyS, Mi_xyS[[3]])
    av_Sigma_xyS <- listsum(av_Sigma_xyS, Mi_xyS[[4]])

    Mi_yS <- maxJoint(suffStat[[i]][c(y,S)], A_yS, X_yS, k_yS, M)

    av_logL_yS <- sum(av_logL_yS, Mi_yS[[1]])
    av_Multinomial_yS <- listsum(av_Multinomial_yS, Mi_yS[[2]])
    av_Mu_yS <- listsum(av_Mu_yS, Mi_yS[[3]])
    av_Sigma_yS <- listsum(av_Sigma_yS, Mi_yS[[4]])

    Mi_xS <- maxJoint(suffStat[[i]][c(x,S)], A_xS, X_xS, k_xS, M)

    av_logL_xS <- sum(av_logL_xS, Mi_xS[[1]])
    av_Multinomial_xS <- listsum(av_Multinomial_xS, Mi_xS[[2]])
    av_Mu_xS <- listsum(av_Mu_xS, Mi_xS[[3]])
    av_Sigma_xS <- listsum(av_Sigma_xS, Mi_xS[[4]])

    if (!is.null(S)) {
      Mi_S <- maxJoint(suffStat[[i]][S], A_S, X_S, k_S, M)

      av_logL_S <- sum(av_logL_S, Mi_S[[1]])
      av_Multinomial_S <- listsum(av_Multinomial_S, Mi_S[[2]])
      av_Mu_S <- listsum(av_Mu_S, Mi_S[[3]])
      av_Sigma_S <- listsum(av_Sigma_S, Mi_S[[4]])
    }
  }

  # evaluate all likelihoods at the average parameter values
  eval_xyS <- sapply(suffStat, evalJoint, c(x,y,S), A_xyS, X_xyS, k_xyS,
                     av_Multinomial_xyS, av_Mu_xyS, av_Sigma_xyS)
  eval_yS  <- sapply(suffStat, evalJoint, c(y,S),   A_yS,  X_yS,  k_yS,
                     av_Multinomial_yS,  av_Mu_yS,  av_Sigma_yS)
  eval_xS  <- sapply(suffStat, evalJoint, c(x,S),   A_xS,  X_xS,  k_xS,
                     av_Multinomial_xS,  av_Mu_xS,  av_Sigma_xS)
  if (!is.null(S)) {
    eval_S   <- sapply(suffStat, evalJoint, S, A_S, X_S, k_S, av_Multinomial_S, av_Mu_S, av_Sigma_S)
  }



  ### mean log likelihood ratio statistic with individual parameters
  if (is.null(S)) {
    LR_bar <- 2* (av_logL_xyS - av_logL_yS - av_logL_xS)
  } else {
    LR_bar <- 2* (av_logL_xyS - av_logL_yS - av_logL_xS + av_logL_S)
  }

  if (is.na(LR_bar)) {return(NaN)}
  if (LR_bar==-Inf) {return(1)}
  if (LR_bar==Inf) {return(0)}

  ### mean log likelihood ratio statistic with averaged parameters
  if (is.null(S)) {
    LR_tilde <- 2* ( mean(eval_xyS) - mean(eval_yS) - mean(eval_xS) )
  } else {
    LR_tilde <- 2* ( mean(eval_xyS) - mean(eval_yS) - mean(eval_xS) + mean(eval_S) )
  }

  r3 <- ( M + 1 ) * ( LR_bar - LR_tilde )  / ( K * ( M - 1 ) )

  if (is.na(r3)) {return(0)}

  D3 <- LR_tilde / (K * (1 + r3))

  # degrees of freedom 2
  t <- K*(M-1)
  if (t > 4) {
    df <- 4 + (t-4) * (1 + (1-2/t) / r3)^2
  } else {
    df <- t * (1 + 1/K) * (1 + 1/r3)^2 /2
  }

  # p value
  pvalue <- pf(D3, K, df, lower.tail = FALSE)

  if (moreOutput) {
    return( c(D3=D3, LR_bar=LR_bar, LR_tilde=LR_tilde, K=K, df=df, pvalue=pvalue) )
  } else {
    return(pvalue)
  }
}




### sum list elements
plus <- function(a, b) {
  if (is.null(a)) {
    return(b)
  } else if (is.null(b)) {
    return(a)
  } else {
    return(a + b) 
  }
}

listsum <- function(l1, l2) {
  mapply(plus, l1, l2, SIMPLIFY=FALSE)
}



### degrees of freedom
dfCG <- function(dat, A, k) {
  fA <- df_f(A, dat)
  dof <-  fA * df_h(k) + fA
  return(dof)
}


# continous variables
df_h <- function(k){ k*(k+1) /2 + k }

# discrete variables
df_f <- function(A, dat) {
  if (length(A)==0) {return(1)}
  num_levels <- sapply(dat[A], function(i){nlevels(i)})
  prod(num_levels)
}



### maximum log Likelihood (Gaussian and multivariate component) in all cells
### 'Joint' because it is not conditional

maxJoint <- function(dat2, A, X, k, M) {
  N <- nrow(dat2)

  if (length(A)==0) {
    # no discrete variables
    maxMultinomial <- list(0)
    maxAll <- mvnorm.mle(Rfast::data.frame.to_matrix(dat2))
    maxMu <- list( maxAll$mu /M )
    logL <- maxAll$loglik /M
    maxSigma <- list( maxAll$sigma /M )
  } else {
    # at least one discrete variable
    x <- Rfast::data.frame.to_matrix(dat2[X])

    Sigma_all <- covm(x)
    # replace zero elements in Sigma_all by a small positive number
    if(length(Sigma_all[Sigma_all == 0]) != 0){
         Sigma_all[Sigma_all == 0] <- Sigma_all[Sigma_all == 0] + 1e-10
         }

    logL_cells <- by(dat2, dat2[A], maxCell, A, X, N, Sigma_all, k, M, simplify = FALSE)
    logL <- sum(unlist(sapply(logL_cells, '[[', 1)))
    maxMultinomial <- lapply(logL_cells, '[[', 2)
    maxMu <- lapply(logL_cells, '[[', 3)
    maxSigma <- lapply(logL_cells, '[[', 4)
  }

  return(list(logL, maxMultinomial, maxMu, maxSigma))
}


### maximum log likelihood (Gaussian and multivariate component) in one cell

maxCell <- function(dat_cell, A, X, N, Sigma_all, k, M) {
  a <- nrow(dat_cell)

  if (a==0) {return(list(0, 0, 0, 0))}

  # multinomial component of the log likelihood:
  maxP <- multinomialLikelihood(a, N) / M
  c1 <- a * maxP

  # multivariate Gaussian component of the log likelihood:
  c2 <- 0
  maxM <- 0
  maxS <- matrix(0)

  # make a numeric variable
  x <- Rfast::data.frame.to_matrix(dat_cell[X])

  if (k > 0) {
    if (a > (k + 5)) {
      Sigma <- covm(x)
    } else {
      Sigma <- Sigma_all
    }

    # replace zero elements in Sigma_all by a small positive number
    if(length(Sigma[Sigma == 0]) != 0){Sigma[Sigma == 0] <- Sigma[Sigma == 0] + 1e-10}

    dec <- tryCatch(chol(Sigma), error = function(e) e)
    if (inherits(dec, "error")) {
      c2 <- 0
      maxM <- apply(x, 2, mean) / M
      maxS <- Sigma / M

    } else {

    c2 <- sum(Rfast::dmvnorm(x = x,
                            mu = apply(x, 2, mean),  #colmeans(x)
                            sigma = matrix(Sigma, nrow = k),
                            logged = TRUE)) / M
    maxM <- apply(x, 2, mean) / M
    maxS <- Sigma / M
    }
  }

  lik <- c1 + c2

  return(list(lik, maxP, maxM, maxS))
}



evalJoint <- function(dat2, vars, A, X, k, maxMultinomial, maxMu, maxSigma) {

  dat2 <- dat2[ ,vars, drop = FALSE]

  if (length(A)==0) {
    logL <- 0
    dec <- tryCatch(chol(maxSigma[[1]]), error = function(e) e)
    if (!inherits(dec, "error")) {
      logL <- sum( Rfast::dmvnorm(Rfast::data.frame.to_matrix(dat2[X]),
                         mu = maxMu[[1]],
                         sigma = maxSigma[[1]], log = TRUE) )
    }
  } else {
    # at least one discrete variable
    cells <- split(dat2[X], dat2[A], drop=FALSE)
    Ncells <- length(cells)
    logL_cells <- mapply(evalCell, cells, rep(k, Ncells),
                         maxMultinomial,
                         maxMu,
                         maxSigma)
    logL <- sum(logL_cells)
  }

  return(logL)
}


evalCell <- function(dat_cell, k, maxMultinomial, maxMu, maxSigma) {
  a <- nrow(dat_cell)

  # if there are averaged parameters for this partition, but no data for this imputed data set
  # or if there are no continuous variables
  if (a==0) {return(0)}

  # if there are no averaged parameters for this partition b/c no imputed data set had data here
  # -> then this data set does not have data here either

  c1 <- a * maxMultinomial

  c2 <- 0
  if (k > 0) {
    dec <- tryCatch(chol(matrix(maxSigma, nrow=k)), error = function(e) e)
    if (!inherits(dec, "error")) {

    c2 <- sum( Rfast::dmvnorm(Rfast::data.frame.to_matrix(dat_cell),
                             mu = maxMu,
                             sigma = matrix(maxSigma, nrow=k),
                             logged=TRUE) )

    }

  }

  return(c1 + c2)
}


### maximum log likelihood multinomial component
multinomialLikelihood <- function(a, N) {
  log(a / N)
}


### Covm
covm <- function(dat) {
  n <- nrow(dat)
  covm <- Rfast::cova(Rfast::data.frame.to_matrix(dat)) *(n-1)/n
  if (anyNA(covm)) {
    covm <- matrix(1)
  }
  return(covm)
}
