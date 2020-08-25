#' Fisher's z-Test for (Conditional) Independence between Gaussian Variables with Missings
#' 
#' A wrapper for \code{\link[pcalg:condIndFisherZ]{gaussCItest}}, 
#' to be used within \code{\link[pcalg]{skeleton}}, \code{\link[pcalg]{pc}} or 
#' \code{\link[pcalg]{fci}} when the data contain missing values. Observations 
#' where at least one of the variables involved in the test is missing are 
#' deleted prior to performing the test (test-wise deletion).
#'
#' @param x,y,S (integer) position of variable X, Y and set of variables S, 
#' respectively, in each correlation matrix in \code{suffStat}. It is tested 
#' whether X and Y are conditionally independent given the subset S of the 
#' remaining variables.
#' @param suffStat \code{data.frame} containing the raw data.
#'
#' @return See \code{\link[pcalg:condIndFisherZ]{gaussCItest}} for details on 
#' Fisher's z-test. Test-wise deletion is valid if missingness does not jointly 
#' depend on X and Y.
#' 
#' @return A p-value.
#' 
#' @seealso \code{\link[pcalg]{gaussCItest}} for complete data, 
#' \code{\link[pcalg]{gaussMItest}} for multiply imputed data
#' 
#' @export
gaussCItwd <- function(x, y, S=NULL, suffStat) {
  
  miss <- apply(suffStat[ ,c(x,y,S)], 1, anyNA)
  C <- stats::cor(suffStat[!miss, c(x,y,S)])
  n <- sum(!miss)
  
  if (length(S)==0) {
    z <- pcalg::zStat(1, 2, NULL, C=C, n=n)
  } else {
    z <- pcalg::zStat(1, 2, 3:(2+length(S)), C=C, n=n)
  }
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}
