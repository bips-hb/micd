#' Estimate Final Skeleton in the FCI-MI algorithm for multiple imputed
#' data sets of continuous data
#'
#' This function is a modification of [functioname(pcalg::pdsep)]
#' to be used for multiple imputation.
#'
#' @param skel Graph object returned by [functioname(micd::skeletonMI)]
#' @param data An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#' @param p Number of variables.
#' @param sepset See [functioname(pcalg::pdsep)] for more details.
#' @param alpha Significance level for the individual conditional independence tests.
#' @param pMax See [functioname(pcalg::pdsep)] for more details.
#' @param m.max Maximum size of the conditioning sets that are considered in the conditional independence tests.
#' @param pdsep.max See [functioname(pcalg::pdsep)] for more details.
#' @param NAdelete  See [functioname(pcalg::pdsep)] for more details.
#' @param unfVect See [functioname(pcalg::pdsep)] for more details.
#' @param biCC See [functioname(pcalg::pdsep)] for more details.
#' @param verbose Logical indicating that detailed output is to be provided.
#'
#' @return See [functioname(pcalg::pdsep)] for more details.
#' @export
#'
pdsepMI <- function (skel, data, p, sepset, alpha, pMax,
    m.max = Inf, pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
    biCC = FALSE, verbose = FALSE)
{
    G <- (methods::as(skel, "matrix") != 0)
    n.edgetests <- rep(0, 1000)
    ord <- 0L
    allPdsep.tmp <- vector("list", p)
    if (biCC)
        conn.comp <- lapply(RBGL::biConnComp(skel), as.numeric)
    if (any(G)) {
        amat <- G
        ind <- which(G, arr.ind = TRUE)
        storage.mode(amat) <- "integer"
        if (verbose)
            cat("\nCompute collider:\n")
        for (i in seq_len(nrow(ind))) {
            x <- ind[i, 1]
            y <- ind[i, 2]
            allZ <- setdiff(which(amat[y, ] != 0), x)
            for (z in allZ) {
                if (amat[x, z] == 0 && !(y %in% sepset[[x]][[z]] ||
                  y %in% sepset[[z]][[x]])) {
                  if (length(unfVect) == 0) {
                    amat[x, y] <- amat[z, y] <- 2
                    if (verbose)
                      cat("\n", x, "*->", y, "<-*", z, "\n")
                  }
                  else {
                    if (!any(unfVect == pcalg::triple2numb(p, x, y,
                      z), na.rm = TRUE) && !any(unfVect == pcalg::triple2numb(p,
                      z, y, x), na.rm = TRUE)) {
                      amat[x, y] <- amat[z, y] <- 2
                      if (verbose)
                        cat("\n", x, "*->", y, "<-*", z, "\n")
                    }
                  }
                }
            }
        }
        allPdsep <- lapply(1:p, pcalg::qreach, amat = amat)
        allPdsep.tmp <- vector("list", p)
        for (x in 1:p) {
            if (verbose)
                cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]],
                  "\n")
            if (any(an0 <- amat[x, ] != 0)) {
                tf1 <- setdiff(allPdsep[[x]], x)
                adj.x <- which(an0)
                for (y in adj.x) {
                  if (verbose)
                    cat(sprintf("\ny = %3d\n.........\n", y))
                  tf <- setdiff(tf1, y)
                  diff.set <- setdiff(tf, adj.x)
                  if (biCC) {
                    for (cci in conn.comp) {
                      if (x %in% cci && y %in% cci)
                        break
                    }
                    bi.conn.comp <- setdiff(cci, c(x, y))
                    tf <- intersect(tf, bi.conn.comp)
                    if (verbose) {
                      cat("There is an edge between", x, "and",
                        y, "\n")
                      cat("Possible D-Sep of", x, "intersected with the biconnected component of",
                        x, "and", y, "is:", tf, "\n")
                    }
                  }
                  allPdsep.tmp[[x]] <- c(tf, y)
                  if (length(tf) > pdsep.max) {
                    if (verbose)
                      cat("Size of Possible-D-SEP bigger than",
                        pdsep.max, ". Break the search for the edge between",
                        x, "and", y, "\n")
                  }
                  else if (length(diff.set) > 0) {
                    done <- FALSE
                    ord <- 0L
                    while (!done && ord < min(length(tf), m.max)) {
                      ord <- ord + 1L
                      if (verbose)
                        cat("ord = ", ord, "\n")
                      if (ord == 1) {
                        for (S in diff.set) {
                          pval <- gaussCItestMI(x, y, S, data)
                          n.edgetests[ord + 1] <- n.edgetests[ord +
                            1] + 1
                          if (is.na(pval))
                            pval <- as.numeric(NAdelete)
                          if (pval > pMax[x, y])
                            pMax[x, y] <- pval
                          if (pval >= alpha) {
                            amat[x, y] <- amat[y, x] <- 0
                            sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                            done <- TRUE
                            if (verbose)
                              cat("x=", x, " y=", y, " S=", S,
                                ": pval =", pval, "\n")
                            break
                          }
                        }
                      }
                      else {
                        tmp.combn <- utils::combn(tf, ord)
                        if (ord <= length(adj.x)) {
                          for (k in seq_len(ncol(tmp.combn))) {
                            S <- tmp.combn[, k]
                            if (!all(S %in% adj.x)) {
                              n.edgetests[ord + 1] <- n.edgetests[ord +
                                1] + 1
                              pval <- gaussCItestMI(x, y, S, data)
                              if (is.na(pval))
                                pval <- as.numeric(NAdelete)
                              if (pMax[x, y] < pval)
                                pMax[x, y] <- pval
                              if (pval >= alpha) {
                                amat[x, y] <- amat[y, x] <- 0
                                sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                                done <- TRUE
                                if (verbose)
                                  cat("x=", x, " y=", y, " S=",
                                    S, ": pval =", pval, "\n")
                                break
                              }
                            }
                          }
                        }
                        else {
                          for (k in seq_len(ncol(tmp.combn))) {
                            S <- tmp.combn[, k]
                            n.edgetests[ord + 1] <- n.edgetests[ord +
                              1] + 1
                            pval <- gaussCItestMI(x, y, S, data)
                            if (is.na(pval))
                              pval <- as.numeric(NAdelete)
                            if (pMax[x, y] < pval)
                              pMax[x, y] <- pval
                            if (pval >= alpha) {
                              amat[x, y] <- amat[y, x] <- 0
                              sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                              done <- TRUE
                              if (verbose)
                                cat("x=", x, " y=", y, " S=",
                                  S, ": pval =", pval, "\n")
                              break
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
        G[amat == 0] <- FALSE
        G[amat == 1] <- TRUE
        G[amat == 2] <- TRUE
    }
    list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp,
        max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
}
