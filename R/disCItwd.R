#' G square Test for (Conditional) Independence between Discrete Variables with Missings
#'
#' A wrapper for \code{pcalg::\link[pcalg]{disCItest}}, to be used within
#' \code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or
#' \code{pcalg::\link[pcalg]{fci}} when the data contain missing values.
#' Observations where at least one of the variables involved in the test is
#' missing are deleted prior to performing the test (test-wise deletion).
#'
#' @param x,y,S  (Integer) position of variable X, Y and set of variables S,
#'               respectively, in \code{suffStat}. It is tested whether X and Y
#'               are conditionally independent given the subset S of the remaining variables.
#'
#' @param suffStat   A list with three elements, \code{"dm"}, \code{"nlev"},
#'                   \code{"adaptDF"}; each corresponding to the above arguments.
#'                   Can be obtained from a data.frame  of factor variables using
#'                   the \code{suffStat} function (see example section)
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
#' ## load data (200 observations)
#' data(gmD)
#' dat <- gmD$x[1:1000,]
#'
#' ## delete some observations of X2 and X3
#' set.seed(123)
#' dat[sample(1:1000, 50), 2] <- NA
#' dat[sample(1:1000, 50), 3] <- NA
#'
#' ## analyse incomplete data
#' # test-wise deletion ==========
#' sufftwd <- getSuff(dat, test = "disCItwd")
#' disCItwd(1, 3, NULL, suffStat = sufftwd)
#' 
#' # list-wise deletion ==========
#' dat2 <- dat[complete.cases(dat), ]
#' suffStat2 <- getSuff(dat2, test = "disCItest", adaptDF = FALSE)
#' disCItest(1, 3, NULL, suffStat = suffStat2)
#' 
#' ## use disCItwd within pcalg::pc ==========
#' pc.fit <- pc(suffStat = sufftwd, indepTest = disCItwd, alpha = 0.1, p = 5)
#' pc.fit
#' 
#' if(require("Rgraphviz", character.only = TRUE, quietly = TRUE)){
#' plot(pc.fit)
#' }
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

  S1 <- apply(suffStat$dm[ , S, drop = FALSE], 2, function(i){length(unique(i))==1})
  S <- S[!S1]
  if (length(S) == 0) {S <- NULL}

  pcalg::gSquareDis(x = 1, y = 2, S = S, dm=suffStat$dm, adaptDF=suffStat$adaptDF, n.min=-1)
}
