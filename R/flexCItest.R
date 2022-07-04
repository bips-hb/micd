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
#
#'
#' @export



flexCItest <- function(x, y, S = NULL, suffStat) {

  if (all(c(x,y,S) %in% suffStat$conpos)) {
    x2 <- as.numeric(which(suffStat$conpos %in% x))
    y2 <- as.numeric(which(suffStat$conpos %in% y))

    if ( is.null(S) ) {S2 <- S} else {S2 <- which(suffStat$conpos %in% S)}
    pval <- pcalg::gaussCItest(x=x2, y=y2, S=S2, suffStat=suffStat$corlist)

  } else if ( all(c(x,y,S) %in% suffStat$dispos) ) {
    x2 <- as.numeric(which(suffStat$dispos %in% x))
    y2 <- as.numeric(which(suffStat$dispos %in% y))

    if ( is.null(S) ) {S2 <- S} else {S2 <- which(suffStat$dispos %in% S)}
    pval <- pcalg::disCItest(x=x2, y=y2, S=S2, suffStat=suffStat$dmlist)

  } else {
    pval <- mixCItest(x=x, y=y, S=S, suffStat=suffStat$dat)
  }

  return(pval)
}

