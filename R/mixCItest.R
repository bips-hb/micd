#' Likelihood Ratio Test for (Conditional) Independence between Mixed Variables
#' 
#' A likelihood ratio test for (conditional) independence between mixed 
#' (continuous and unordered categorical) variables, to be used within 
#' \code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or \code{pcalg::\link[pcalg]{fci}}. 
#' It assumes that the variables in the test follow a Conditional Gaussian 
#' distribution, i.e. conditional on each combination of values of the discrete
#'  variables, the continuous variables are multivariate Gaussian. 
#'  Each multivariate Gaussian distribution is allowed to have its own mean 
#'  vector and covariance matrix.
#'
#' @param x,y,S    (integer) position of variable X, Y and set of variables S, 
#' respectively, in \code{suffStat}. It is tested whether X and Y are conditionally 
#' independent given the subset S of the remaining variables.
#' @param suffStat  \code{data.frame}. Discrete variables must be coded as factors.
#'
#' @details The implementation follows Andrews et al. (2018). The same test is also 
#' implemented in TETRAD and in the R-package rcausal, a wrapper for the 
#' TETRAD Java library. Small differences in the p-values returned by 
#' CGtest and the TETRAD/rcausal equivalent are due to differences in 
#' handling sparse or empty cells, and differences in covariance estimation.
#' 
#' Note that the test assumes a Conditional Gaussian distribution 
#' (see Lauritzen & Wermuth, 1989) for the variables \emph{in that specific test}, 
#' rather than for all variables in the graph / data set. Due to the 
#' non-collapsibility of the Conditional Gaussian distribution, this implies 
#' that the assumptions of two tests with overlapping sets of variables may 
#' contradict each other. Little is known about how this affects the quality of 
#' the estimated graph in practice.
#' 
#' @return A p-value.
#' 
#' @author Janine Witte
#' 
#' @references Andrews B., Ramsey J., Cooper G.F. (2018): Scoring Bayesian networks 
#' of mixed variables. \emph{International Journal of Data Science and Analytics} 6:3-18.
#' 
#' Lauritzen S.L., Wermuth N. (1989): Graphical models for associations between 
#' variables, some of which are qualitative and some quantitative. 
#' \emph{The Annals of Statistics} 17(1):31-57. 
#' 
#' Scheines R., Spirtes P., Glymour C., Meek C., Richardson T. (1998): 
#' The TETRAD project: Constraint based aids to causal model specification. 
#' \emph{Multivariate Behavioral Research} 33(1):65-117. 
#' http://www.phil.cmu.edu/tetrad/index.html
#'
#' @examples 
#' ## load data (numeric and factor variables)
#' dat <- toenail2[1:400, ]
#' ß
#' ## analyse data
#' mixCItest(5, 2, NULL, suffStat=dat)
#' mixCItest(2, 3, 4, suffStat=dat)
#'
#' @export


mixCItest <- function(x, y, S=NULL, suffStat) {
  # tests X _||_ Y | S
  
  logL_xyS <- likelihoodJoint(suffStat[c(x,y,S)])
  logL_yS <- likelihoodJoint(suffStat[c(y,S)])
  logL_xS <- likelihoodJoint(suffStat[c(x,S)])
  logL_S <- likelihoodJoint(suffStat[S])
  
  logLR <- 2* ( logL_xyS[1] - logL_yS[1] - logL_xS[1] + logL_S[1] )
  df <- logL_xyS[2] - logL_yS[2] - logL_xS[2] + logL_S[2]
  if (df<=0) {df <- 1}
  
  p <- stats::pchisq(logLR, df, lower.tail=FALSE)
  if (is.na(p)) {p <- 1}
  
  return(p)
}



# Maximum log Likelihood in all cells - 'Joint' because it is not conditional
likelihoodJoint <- function(dat2) {
  
  if (ncol(dat2)==0) {return(c(0,1))}
  
  N <- nrow(dat2)
  
  # indices of discrete variables
  A <- Rfast::which.is(dat2, "factor")
  # indices of continous variables
  X <- Rfast::which.is(dat2, "numeric")
  k <- length(X)
  
  if (length(A)==0) {
    # no discrete variables
    logL <- Rfast::mvnorm.mle(Rfast::data.frame.to_matrix(dat2))$loglik
  } else {
    # at least one discrete variable
    Sigma_all <- covm(dat2[X])
    logL_cells <- by(dat2, dat2[A], likelihoodCell, X, N, Sigma_all, k, simplify=FALSE)
    logL <- sum(unlist(logL_cells))
  }
  
  fA <- df_f(A, dat2)
  dof <-  fA * df_h(k) + fA
  
  return(c(logL, dof))
}




# Maximum log likelihood in one cell
likelihoodCell <- function(dat_cell, X, N, Sigma_all, k) {
  a <- nrow(dat_cell)
  
  # multinomial component of the log likelihood
  c1 <- a * multinomialLikelihood(a, N)
  
  # multivariate Gaussian component of the log likelihood
  c2 <- 0
  if (k > 0) {
    if (a > (k + 5)) {
      #if (a > (k)) {
      c2 <- Rfast::mvnorm.mle(Rfast::data.frame.to_matrix(dat_cell[X]))$loglik
    } else {
      dec <- tryCatch(chol(Sigma_all), error = function(e) e)
      if (!inherits(dec, "error")) {
        c2 <-sum(Rfast::dmvnorm(x = Rfast::data.frame.to_matrix(dat_cell[X]),
                        mu=apply(dat_cell[X], 2, mean),
                        sigma=Sigma_all,
                        logged=TRUE) )
      }
    }
  }
  
  return(c1 + c2)
}



# Maximum log likelihood multinomial component
multinomialLikelihood <- function(a, N) {
  log(a/N)
}


# Covm
covm <- function(dat) {
  n <- nrow(dat)
  covm <- Rfast::cova(Rfast::data.frame.to_matrix(dat)) *(n-1)/n
  if (anyNA(covm)) {
    covm <- matrix(1)
  }
  return(covm)
}

# Degrees of freedom for continous variables
df_h <- function(k){ k*(k+1) /2 }


# Degrees of freedom for discrete variables
df_f <- function(A, dat) {
  if (length(A)==0) {return(1)}
  num_levels <- sapply(dat[A], function(i){nlevels(i)})
  prod(num_levels)
}


