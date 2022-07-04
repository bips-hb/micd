#' Wrapper for gaussMItest, disMItest and mixMItest
#'
#' A plug-in conditional independence test for \code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or
#' \code{pcalg::\link[pcalg]{fci}} when multiply imputed data sets are available. \code{flexMItest} detects whether
#' variables are continuous, discrete or mixed, and automatically switches between \code{\link{gaussMItest}} (continuous only),
#' \code{link{disMItest}} (discrete only) and \code{\link{mixMItest}} (mixed).
#'
#' @param x,y,S (integer) position of variable X, Y and set of variables S,
#' respectively, in the dataset. It is tested whether X and Y are conditionally
#' independent given the subset S of the remaining variables.
#' @param suffStat a list generated using \code{\link{getSuff}} with \code{test="flexMItest"}. See below for details.
#'
#' @details \code{suffStat} needs to be a list with four elements named \code{datlist}, \code{corlist},
#' \code{conpos} and \code{dispos}. \code{datlist} is the list of imputed datasets. \code{corlist}
#' is a list with M+1 elements, where M is the number of imputed datasets. For i=1,...,M, the
#' the i-th element of \code{corlist} is the correlation matrix of the continuous variables in the i-th imputed dataset;
#' the (M+1)-the element is the number of rows in each imputed dataset.
#' \code{conpos} is a vector containing the integer positions of the continuous variables in the original dataset.
#' \code{dispos} is a vector containing the integer positions of the discrete variables in the original dataset.
#'
#'
#' @return A p-value.
#'
#' @seealso \code{\link{gaussMItest}}, \code{\link{disMItest}} and \code{\link{mixMItest}}
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
#' ## impute missing values using random forests
#' imp <- mice(dat, method="rf")
#'
#' ## obtain correct input 'suffStat' for 'flexMItest'
#' suff <- getSuff(imp, test="flexMItest")
#'
#' ## analyse data
#' # continuous variables only
#' flexMItest(4,5,NULL, suffStat=suff)
#' implist <- complete(imp, action="all")
#' gaussSuff <- c(lapply(implist, function(i){cor(i[ ,c(4,5)])}), n=400)
#' gaussMItest(1,2,NULL, suffStat=gaussSuff)
#'
#' # discrete variables only
#' flexMItest(2,3,NULL, suffStat=suff)
#' disMItest(2,3,NULL, suffStat=complete(imp, action="all"))
#'
#' # mixed variables
#' flexMItest(2,3,4, suffStat=suff)
#' mixMItest(2,3,4, suffStat=complete(imp, action="all"))
#'
#' @export



flexMItestFlex <- function(x, y, S = NULL, suffStat) {

  if ( all(c(x,y,S) %in% suffStat$conpos) ) {
    x2 <- as.numeric(which(suffStat$conpos %in% x))
    y2 <- as.numeric(which(suffStat$conpos %in% y))

    if ( is.null(S) ) {S2 <- S} else {S2 <- which(suffStat$conpos %in% S)}
    pval <- gaussMItest(x=x2, y=y2, S=S2, suffStat=suffStat$corlist)
    pval <- ifelse(pval < 0.05, 0, .99)

  } else if ( all(c(x,y,S) %in% suffStat$dispos) ) {
    pval <- disMItest(x=x, y=y, S=S, suffStat=suffStat$datlist)
    pval <- ifelse(pval < 0.1, 0, .99)

  } else {
    pval <- mixMItest(x=x, y=y, S=S, suffStat=suffStat$datlist)
    pval <- ifelse(pval < 0.1, 0, .99)
  }

  return(pval)
}
