#' Test Conditional Independence of Gaussians via Fisher's Z using
#' multiple imputations
#'
#' This function is a modification of [functioname(pcalg::gaussCItest)]
#' to be used for multiple imputation.
#'
#' @param x,y,S (integer) position of variable X, Y and set of variables S,
#'              respectively, in the adjacency matrix.
#'              It is tested, whether X and Y are conditionally independent
#'              given the subset S of the remaining nodes.
#' @param data  An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#'
#' @return  Returns the p-value of the test
#'
#' @export
#' @examples
#' library(mice)
#' daten <- windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#'
#' ## Impute missing values
#' imp <- mice(daten)
#'
#' ## analyse data
#' gaussCItestMI(1,2,c(4,5), imp)
#' gaussCItest(1,2,c(4,5), suffStat = list(C = cor(windspeed), n = nrow(windspeed)))
#' gaussCItest(1,2,c(4,5), suffStat = list(C = cor(daten[complete.cases(daten),]),
#'             n = nrow(daten[complete.cases(daten),])))
#'
gaussCItestMI <- function (x, y, S, data)
{
  M <- data$m
  z <- vector(mode = "list", length = M)
  n <- ifelse(is.mids(data), nrow(data$data), nrow(data$data[[1]]))

  for (i in 1:M)
  {
   if(is.mids(data)){
      data.i <- mice::complete(data, i)
    } else if(is.list(data$data)) {
      data.i <- data$data[[i]]
    } else stop("data is neither a list nor a mids object")

    #suffStat <- list(C = cor(data.i), n = nrow(data.i))
    suffStat <- list(C = cov(scale(data.i)), n = nrow(data.i))
    z[[i]] <- zStatMI(x, y, S, C = suffStat$C, n = n)
  }

  # 1. Average of M imputed data sets
  avgz <- Reduce('+', z) / M

  # 2. Average of complete-data variance
  W <- 1 / (n - length(S) - 3)

  # 3. Between variance
  B <- mean(sapply(z, function(x) (x - avgz)^2))

  # 4. Total variance
  TV <- W + (1 + 1 / M) * B

  # 5. Test statistic
  ts <- avgz / sqrt(TV)

  # 6. Degrees of freedom
  df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2

  # 7. pvalue
  2 * pt(abs(ts), df = df, lower.tail = FALSE)

}



#' @export
zStatMI <- function (x, y, S, C, n)
{
    r <- pcalg::pcorOrder(x, y, S, C)
    res <- 0.5 * log.q1pm(r)
    if (is.na(res))
        0
    else res
}


#' @export
log.q1pm <- function(r) log1p(2*r/(1-r))
