#' Likelihood Ratio Test for (Conditional) Independence between Mixed Variables with Missings
#' 
#' A version of \code{\link{mixCItest}}, to be used within \code{\link[pcalg]{skeleton}}, 
#' \code{\link[pcalg]{pc}} or \code{\link[pcalg]{fci}} when the data contain missing values. 
#' Observations where at least one of the variables involved in the test is missing 
#' are deleted prior to performing the test (test-wise deletion).
#'
#' @param x,y,S (integer) position of variable X, Y and set of variables S, 
#' respectively, in \code{suffStat}. It is tested whether X and Y are conditionally 
#' independent given the subset S of the remaining variables.
#' @param suffStat \code{data.frame}. Discrete variables must be coded as factors.
#'
#' @details See \code{\link{mixCItest}} for details on the assumptions of the 
#' Conditional Gaussian likelihood ratio test. Test-wise deletion is valid if 
#' missingness does not jointly depend on X and Y.
#' 
#' @return A p-value.
#' 
#' @seealso \code{\link{mixCItest}} for complete data, 
#' \code{\link{mixMItest}} for multiply imputed data
#' 
#' @export
mixCItwd <- function(x, y, S=NULL, suffStat) {
  
  miss <- apply(suffStat[ ,c(x,y,S)], 1, anyNA)
  mixCItest(x=x, y=y, S=S, suffStat=suffStat[!miss, ])
  
}
