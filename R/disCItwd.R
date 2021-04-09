#' G square Test for (Conditional) Independence between Discrete Variables with Missings
#' 
#' A wrapper for \code{pcalg::\link[pcalg]{disCItest}}, to be used within 
#' \code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or 
#' \code{pcalg::\link[pcalg]{fci}} when the data contain missing values. 
#' Observations where at least one of the variables involved in the test is 
#' missing are deleted prior to performing the test (test-wise deletion).
#'
#' @param x,y,S  (integer) position of variable X, Y and set of variables S, 
#'               respectively, in \code{suffStat}. It is tested whether X and Y 
#'               are conditionally independent given the subset S of the remaining variables. 
#'               
#' @param dm  data matrix (rows: samples, columns: variables) with integer entries; the k levels for a given column must be coded by the integers 0,1,...,k-1. (see example)
#' 
#' @param nlev  optional vector with numbers of levels for each variable in \code{dm}.
#' 
#' @param adaptDF  logical specifying if the degrees of freedom should be lowered by one for each zero count. The value for the degrees of freedom cannot go below 1.
#'               
#' @param suffStat   a list with three elements, \code{"dm"}, \code{"nlev"}, 
#'                   \code{"adaptDF"}; each corresponding to the above arguments.
#'
#' @details
#' See \code{\link[pcalg]{disCItest}} for details on the G square test. Test-wise deletion 
#' is valid if missingness does not jointly depend on X and Y.
#' 
#' @return A p-value.
#' 
#' @seealso \code{pcalg::\link[pcalg]{disCItest}} for complete data, \code{\link{disMItest}}
#' for multiply imputed data
#' 
#' @examples
#' 
#' ## simulate variables (integers 0,...,k)
#' n <- 200
#' set.seed(123)
#' x <- sample(0:2, n, TRUE) # 3 levels
#' y <- sample(0:3, n, TRUE) # 4 levels
#' z <- sample(0:1, n, TRUE) # 2 levels
#' dat <- cbind(x,y,z)
#' 
#' ## delete some observations of y and z
#' dat[sample(1:n, 40), 2] <- NA
#' dat[sample(1:n, 40), 3] <- NA
#' 
#' ## analyse incomplete data
#' # test-wise deletion:
#' suffStat <- list(dm = dat, nlev = c(3,4,2), adaptDF = FALSE)
#' disCItwd(1, 3, NULL, suffStat = suffStat)
#' # list-wise deletion:
#' suffStat2 <- list(dm = dat[complete.cases(dat), ], nlev = c(3,4,2), adaptDF = FALSE)
#' disCItest(1, 3, NULL, suffStat = suffStat2)
#' 
#' @export


disCItwd <- function(x, y, S=NULL, suffStat) {
  miss <- apply(suffStat$dm[, c(x, y, S)], 1, anyNA)
  
  if (sum(!miss) < 2) { return(NA) }
  
  if (length(unique(suffStat$dm[ ,1]))==1) {return(NA)}
  if (length(unique(suffStat$dm[ ,2]))==1) {return(NA)}
  
  suffStat$dm <- suffStat$dm[!miss, c(x, y, S)]
  if (length(S) > 0) {
    S <- 3:ncol(suffStat$dm)
  }
  
  S1 <- apply(suffStat$dm[ ,S, drop=FALSE], 2, function(i){length(unique(i))==1})
  S <- S[!S1]
  if (length(S)==0) {S <- NULL}
  
  gSquareDis(x = 1, y = 2, S = S, dm=suffStat$dm, adaptDF=suffStat$adaptDF, n.min=-1)
}
