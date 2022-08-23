#' Test Conditional Independence of Gaussians via Fisher's Z using
#' multiple imputations
#'
#' A modified version of \code{pcalg::\link[pcalg:condIndFisherZ]{gaussCItest}}, 
#' to be used within
#' \code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or
#' \code{pcalg::\link[pcalg]{fci}} when multiply imputated data sets are available.
#'
#' @param x,y,S (integer) position of variable X, Y and set of variables S,
#'              respectively, in the adjacency matrix.
#'              It is tested, whether X and Y are conditionally independent
#'              given the subset S of the remaining nodes.
#' @param data  An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#' @param suffStat A list of length m+1, where m is the number of imputations; 
#'                 the first m elements are the covariance matrices of the m
#'                 imputed data sets, the m-th element is the sample size. Can
#'                 be obtained from a mids object by
#'                getSuff(mids, test="gaussMItest")
#'                
#' @details \code{gaussMItest} is faster, as it uses pre-calculated covariance matrices.
#'
#' @return  Returns the p-value of the test
#' 
#' @importFrom stats cov pt
#' @import mice
#'
#' @examples
#' ## load data (numeric variables)
#' dat <- as.matrix(windspeed)
#' 
#' ## delete some observations
#' set.seed(123)
#' dat[sample(1:length(dat), 260)] <- NA
#' 
#' ## Impute missing values under normal model
#' imp <- mice(dat, method = "norm")
#' 
#' ## analyse data
#' # complete data:
#' suffcomplete <- getSuff(windspeed, test = "gaussCItest")
#' gaussCItest(1, 2, c(4,5), suffStat = suffcomplete)
#' # multiple imputation:
#' suffMI <- getSuff(imp, test = "gaussMItest")
#' gaussMItest(1, 2, c(4,5), suffStat = suffMI)
#' gaussCItestMI(1, 2, c(4,5), data = imp)
#' # test-wise deletion:
#' gaussCItwd(1, 2, c(4,5), suffStat = dat)
#' # list-wise deletion:
#' dat2 <- dat[complete.cases(dat), ]
#' sufflwd <- getSuff(dat2, test = "gaussCItest")
#' gaussCItest(1, 2, c(4,5), suffStat = sufflwd)
#' 
#' ## use gaussMItest or gaussCItestMI within pcalg::pc
#' pc.fit <- pc(suffStat = suffMI, indepTest = gaussMItest, alpha = 0.01, p = 6)
#' pc.fit
#' pc.fit <- pc(suffStat = imp, indepTest = gaussCItestMI, alpha = 0.01, p = 6)
#' pc.fit
#' 
#' @export

gaussMItest <- function (x, y, S, suffStat) {
  # number of imputations
  M <- length(suffStat) - 1
  # sample size
  n <- suffStat[[M+1]]
  suffStat[[M+1]] <- NULL
  
  z <- sapply(suffStat, function(j) {
    zStatMI(x, y, S, C=j, n=n)
  })
  
  # 1. Average of M imputed data sets
  avgz <- mean(z)
  
  # 2. Average of completed-data variance
  W <- 1 / (n - length(S) - 3)
  
  # 3. Between variance
  B <- sum( ( z - avgz )^2 ) / (M-1)
  
  # 4. Total variance
  TV <- W + (1 + 1 / M) * B
  
  # 5. Test statistic
  ts <- avgz / sqrt(TV)
  
  # 6. Degrees of freedom
  df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  
  # 7. pvalue
  pvalue <- 2 * stats::pt(abs(ts), df = df, lower.tail = FALSE)
  
  return(pvalue)
}

#' @export
#' @rdname gaussMItest
gaussCItestMI <- function (x, y, S=NULL, data)
{
  M <- data$m
  z <- vector(mode = "list", length = M)
  n <- ifelse(mice::is.mids(data), nrow(data$data), nrow(data$data[[1]]))
  
  for (i in 1:M)
  {
    if(mice::is.mids(data)){
      data.i <- mice::complete(data, i)
    } else if(is.list(data$data)) {
      data.i <- data$data[[i]]
    } else stop("data is neither a list nor a mids object")
    
    suffStat <- list(C = stats::cov(scale(data.i)), n = nrow(data.i))
    z[[i]] <- zStatMI(x, y, S, C = suffStat$C, n = n)
  }
  
  
  # 1. Average of M imputed data sets
  avgz <- Reduce('+', z) / M
  
  # 2. Average of complete-data variance
  W <- 1 / (n - length(S) - 3)
  
  # 3. Between variance
  # Version 0.2.0: B.old <- mean(sapply(z, function(x) (x - avgz)^2))
  B <- sum(sapply(z, function(x) (x - avgz)^2)) / (M - 1)
  
  # 4. Total variance
  TV <- W + (1 + 1 / M) * B
  
  # 5. Test statistic
  ts <- avgz / sqrt(TV)
  
  # 6. Degrees of freedom
  df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  
  # 7. pvalue
  2 * stats::pt(abs(ts), df = df, lower.tail = FALSE)
  
}


zStatMI <- function (x, y, S, C, n)
{
  r <- pcalg::pcorOrder(x, y, S, C)
  res <- 0.5 * log_q1pm(r)
  if (is.na(res))
    0
  else res
}

log_q1pm <- function(r) log1p(2 * r / (1 - r))
