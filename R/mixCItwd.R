#' Likelihood Ratio Test for (Conditional) Independence between Mixed Variables with Missings
#'
#' A version of \code{\link{mixCItest}}, to be used within \code{pcalg::\link[pcalg]{skeleton}},
#' \code{pcalg::\link[pcalg]{pc}} or \code{pcalg::\link[pcalg]{fci}} when the data contain missing values.
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
#' @examples
#' ## load data (numeric and factor variables)
#' dat <- toenail2[1:400, ]
#'
#' ## delete some observations
#' set.seed(123)
#' dat[sample(400, 20), 2] <- NA
#' dat[sample(400, 30), 4] <- NA
#'
#' ## analyse data
#' # complete data:
#' mixCItest(2, 3, 5, suffStat=toenail2[1:400, ])
#' # test-wise deletion:
#' mixCItwd(2, 3, 5, suffStat = dat)
#' # list-wise deletion:
#' dat2 <- dat[complete.cases(dat), ]
#' mixCItest(2, 3, 5, suffStat = dat2)
#' 
#' ## use mixCItwd within pcalg::pc
#' pc.fit <- pc(suffStat = dat, indepTest = mixCItwd, alpha = 0.01, p = 5)
#' pc.fit
#'
#' @export



mixCItwd <- function(x, y, S = NULL, suffStat) {

  miss <- apply(suffStat[ ,c(x,y,S)], 1, anyNA)

  if (sum(!miss)==0) { return(NA) }

  mixCItest(x = x, y = y, S = S, suffStat = suffStat[!miss, ])
}
