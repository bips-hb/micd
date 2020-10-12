#' Utility for conservative and majority rule in PC-MI and FCI-MI
#' applied for multiple imputed data sets of continuous data
#'
#' This function is a modification of \code{pcalg::\link[pcalg]{pc}}
#' to be used for multiple imputation.
#'
#' @param sk A skeleton object as returned from \code{\link{skeletonMI}}
#' @param data An object of type mids, which stands for 'multiply imputed
#'             data set', typically created by a call to function mice().
#' @param alpha Significance level for the conditional independence test.
#' @param version.unf (Integer) vector of length two. See
#'                    \code{pcalg::\link[pcalg]{checkTriple}} for more details.
#' @param maj.rule Logical indicating how the majority rule is applied. See
#'                 \code{pcalg::\link[pcalg]{checkTriple}} for more details.
#' @param verbose Logical asking for detailed output of intermediate steps.
#'
#' @return See \code{pcalg::\link[pcalg]{pc.cons.intern}} for more details.
#'
#' @note This is a modified function of \code{pcalg::\link[pcalg]{pc.cons.intern}}
#'       from the package 'pcalg' (Kalisch et al., 2012;
#'       http://www.jstatsoft.org/v47/i11/).
#'
#' @export
#'
pc.cons.internMI <- function (sk, data, alpha, version.unf = c(NA, NA),
                              maj.rule = FALSE, verbose = FALSE)
{
    g <- methods::as(sk@graph, "matrix")
    stopifnot(all(g == t(g)))
    p <- as.numeric(dim(g)[1])
    unfTripl <- vers <- rep(NA, min(p * p, 1e+05))
    counter <- 0
    if (sum(g) > 0) {
        ind <- which(g == 1, arr.ind = TRUE)
        tripleMatrix <- NULL
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            allC <- setdiff(which(g[b, ] == 1), a)
            newC <- allC[g[a, allC] == 0]
            tmpMatrix <- cbind(rep(a, length(newC)), rep(b, length(newC)),
                newC)
            tripleMatrix <- rbind(tripleMatrix, tmpMatrix)
            colnames(tripleMatrix) <- c("", "", "")
        }
        if ((m <- nrow(tripleMatrix)) > 0) {
            deleteDupl <- logical(m)
            for (i in seq_len(m)) if (tripleMatrix[i, 1] > tripleMatrix[i,
                3])
                deleteDupl[i] <- TRUE
            if (any(deleteDupl))
                tripleMatrix <- tripleMatrix[!deleteDupl, , drop = FALSE]
            for (i in seq_len(nrow(tripleMatrix))) {
                if (counter + 1L == length(unfTripl)) {
                  n.xtra <- min(p * p, 1e+05)
                  new.len <- counter + 1L + n.xtra
                  length(unfTripl) <- new.len
                  length(vers) <- new.len
                }
                a <- tripleMatrix[i, 1]
                b <- tripleMatrix[i, 2]
                c <- tripleMatrix[i, 3]
                nbrsA <- which(g[, a] != 0)
                nbrsC <- which(g[, c] != 0)
                if (verbose) {
                  cat("\nTriple:", a, b, c, "and sepset by skelet:",
                    unique(sk@sepset[[a]][[c]], sk@sepset[[c]][[a]]),
                    "\n")
                }

                r.abc <- checkTripleMI(a, b, c, nbrsA, nbrsC, data = data,
                                       sepsetA = sk@sepset[[a]][[c]],
                                       sepsetC = sk@sepset[[c]][[a]],
                                       alpha = alpha,
                                       version.unf = version.unf,
                                       maj.rule = maj.rule,
                                       verbose = verbose)
                if (r.abc$decision == 3) {
                  counter <- counter + 1
                  unfTripl[counter] <- pcalg::triple2numb(p, a, b, c)
                  vers[counter] <- r.abc$version
                }
                if ((version.unf[1] == 2) && (r.abc$version ==
                  2) && (r.abc$decision != 3)) {
                  counter <- counter + 1
                  unfTripl[counter] <- pcalg::triple2numb(p, a, b, c)
                  vers[counter] <- r.abc$version
                }
                sk@sepset[[a]][[c]] <- r.abc$SepsetA
                sk@sepset[[c]][[a]] <- r.abc$SepsetC
            }
        }
    }
    length(unfTripl) <- length(vers) <- counter
    list(unfTripl = unfTripl, vers = vers, sk = sk)
}
