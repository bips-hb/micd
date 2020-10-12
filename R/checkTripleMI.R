#' Check Consistency of Conditional Independence for a Triple of Nodes
#' used in FCI-MI for multiple imputed data sets of continuous data
#'
#' This function is a modification of \code{pcalg::\link[pcalg]{checkTriple}}
#' to be used for multiple imputation.
#'
#' @param a,b,c   (Integer) positions in adjacency matrix for nodes a, b, and c,
#'                 respectively.
#' @param nbrsA,nbrsC  (Integer) position in adjacency matrix for neighbors of
#'                      a and c, respectively.
#' @param sepsetA  Vector containing Sepset(a,c).
#' @param sepsetC  Vector containing Sepset(c,a).
#' @param data     An object of type mids, which stands for 'multiply imputed
#'                 data set', typically created by a call to function mice().
#' @param alpha    significance level of test.
#' @param version.unf (Integer) vector of length two. See
#'                    \code{pcalg::\link[pcalg]{checkTriple}} for more details.
#' @param maj.rule Logical indicating how the majority rule is applied. See
#'                    \code{pcalg::\link[pcalg]{checkTriple}} for more details.
#' @param verbose  Logical asking for detailed output of intermediate steps.
#'
#' @return See \code{pcalg::\link[pcalg]{checkTriple}} for details.
#'
#' @note This is a modified function of \code{pcalg::\link[pcalg]{checkTriple}}
#'       from the package 'pcalg' (Kalisch et al., 2012;
#'       http://www.jstatsoft.org/v47/i11/).

#' @export
checkTripleMI <- function (a, b, c, nbrsA, nbrsC, data, sepsetA, sepsetC,
                    alpha, version.unf = c(NA, NA), maj.rule = FALSE,
                    verbose = FALSE)
{
    nr.indep <- 0
    stopifnot(length(version.unf) == 2, version.unf %in% 1:2)
    tmp <- if (version.unf[2] == 2)
        (b %in% sepsetA || b %in% sepsetC)
    version <- 0
    if ((nn <- length(nbrsA)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        for (i in 1:nrow(allComb)) {
            S <- nbrsA[which(allComb[i, ] != 0)]
            pval <- gaussCItestMI(a, c, S, data)
            if (verbose)
                cat("a: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if ((nn <- length(nbrsC)) > 0) {
        allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
        for (i in 1:nrow(allComb)) {
            S <- nbrsC[which(allComb[i, ] != 0)]
            pval <- gaussCItestMI(a, c, S, data)
            if (verbose)
                cat("c: S =", S, " - pval =", pval, "\n")
            if (pval >= alpha) {
                nr.indep <- nr.indep + 1
                tmp <- c(tmp, b %in% S)
                version <- 1
            }
        }
    }
    if (version.unf[1] == 2 && nr.indep == 0) {
        version <- 2
    }
    if (is.null(tmp))
        tmp <- FALSE
    if (all(tmp)) {
        res <- 2
        if (b %nin% sepsetA)
            sepsetA <- c(sepsetA, b)
        if (b %nin% sepsetC)
            sepsetC <- c(sepsetC, b)
    }
    else {
        if (all(!tmp)) {
            res <- 1
            sepsetA <- setdiff(sepsetA, b)
            sepsetC <- setdiff(sepsetC, b)
        }
        else {
            if (!maj.rule) {
                res <- 3
            }
            else {
                if (sum(tmp)/length(tmp) < 0.5) {
                  res <- 1
                  sepsetA <- setdiff(sepsetA, b)
                  sepsetC <- setdiff(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) > 0.5) {
                  res <- 2
                  if (b %nin% sepsetA)
                    sepsetA <- c(sepsetA, b)
                  if (b %nin% sepsetC)
                    sepsetC <- c(sepsetC, b)
                }
                else if (sum(tmp)/length(tmp) == 0.5) {
                  res <- 3
                }
            }
        }
    }
    if (verbose && res == 3)
        cat("Triple ambiguous\n")
    lapply(list(decision = res, version = version, SepsetA = sepsetA,
        SepsetC = sepsetC), as.integer)
}


`%nin%` <- function (x, table) is.na(match(x, table))
