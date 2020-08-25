#' G square Test for (Conditional) Independence between Discrete Variables with Missings
#' 
#' A wrapper for \code{\link[pcalg]{disCItest}}, to be used within 
#' \code{\link[pcalg]{skeleton}}, \code{\link[pcalg]{pc}} or 
#' \code{\link[pcalg]{fci}} when the data contain missing values. 
#' Observations where at least one of the variables involved in the test is 
#' missing are deleted prior to performing the test (test-wise deletion).
#'
#' @param x,y,S  (integer) position of variable X, Y and set of variables S, 
#'               respectively, in \code{suffStat}. It is tested whether X and Y 
#'               are conditionally independent given the subset S of the remaining variables. 
#' @param suffStat   a list with three elements, \code{"dm"}, \code{"nlev"}, 
#'                   \code{"adaptDF"}; each corresponding to the above arguments.
#'
#' @details
#' See \code{\link{disCItest}} for details the G square test. Test-wise deletion 
#' is valid if missingness does not jointly depend on X and Y.
#' 
#' @return A p-value.
#' 
#' @seealso \code{\link[pcalg]{disCItest}} for complete data, \code{\link[pcalg]{disMItest}}
#' for multiply imputed data
#' 
#' @export
disCItwd <- function(x, y, S=NULL, suffStat) {
  
  miss <- apply(suffStat$dm[ ,c(x,y,S)], 1, anyNA)
  suffStat$dm <- suffStat$dm[!miss,c(x,y,S)]
  if (length(S)>0) {S <- 3:ncol(suffStat$dm)}
  pcalg::disCItest(x=1, y=2, S=S, suffStat = suffStat)
}
