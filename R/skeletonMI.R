#' Estimate (Initial) Skeleton of a DAG using the PC Algorithm for multiple
#' imputed data sets of continuous data
#'
#' This function is a modification of [pcalg::skeleton()]
#' to be used for multiple imputation.
#'
#' @param data  An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#' @param alpha  Significance level
#' @param labels (Optional) character vector of variable (or "node") names.
#'               Typically preferred to specifying p
#' @param p      (Optional) number of variables (or nodes). May be specified if
#'               labels are not, in which case labels is set to 1:p.
#' @param method  Character string specifying method; the default, "stable"
#'                provides an order-independent skeleton, see
#'                [pcalg::pc()] for details.
#' @param m.max  Maximal size of the conditioning sets that are considered in
#'               the conditional independence tests.
#' @param fixedGaps Logical symmetric matrix of dimension p*p. If entry \[i,j\]
#'                  is true, the edge i-j is removed before starting the
#'                  algorithm. Therefore, this edge is guaranteed to be absent
#'                  in the resulting graph.
#' @param fixedEdges A logical symmetric matrix of dimension p*p. If entry \[i,j\]
#'                   is true, the edge i-j is never considered for removal.
#'                   Therefore, this edge is guaranteed to be present in the
#'                   resulting graph.
#' @param NAdelete Logical needed for the case indepTest(*) returns NA.
#'                 If it is true, the corresponding edge is deleted, otherwise not.
#' @param verbose  If TRUE, detailed output is provided.
#'
#' @return See [pcalg::skeleton()] for more details.
#' 
#' @import pcalg
#'
#' @note This is a modified function of [pcalg::skeleton()]
#'       from the package 'pcalg' (Kalisch et al., 2012;
#'       http://www.jstatsoft.org/v47/i11/).
#'       
#' @author Original code by Markus Kalisch, Martin Maechler, Alain Hauser, and Diego Colombo.
#' Modifications by Ronja Foraita.         
#'
#' @export
#'
#' @examples
#' 
#' data(gmG)
#' n <- nrow(gmG8$x)
#' V <- colnames(gmG8$x) # labels aka node names
#' ## estimate Skeleton
#' data_mids <- mice(gmG8$x)
#' (skel.fit <- skeletonMI(data = data_mids, alpha = 0.01, labels = V, verbose = FALSE))

skeletonMI <- function (data, alpha, labels, p,
                        method = c("stable", "original"),
                        m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
                        NAdelete = TRUE, verbose = FALSE)
{
    cl <- match.call()
    if (!missing(p))
        stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
    if (missing(labels)) {
        if (missing(p))
            stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    }
    else {
        stopifnot(is.character(labels))
        if (missing(p))
            p <- length(labels)
        else if (p != length(labels))
            stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    }
    seq_p <- seq_len(p)
    method <- match.arg(method)
    if (is.null(fixedGaps)) {
        G <- matrix(TRUE, nrow = p, ncol = p)
    }
    else if (!identical(dim(fixedGaps), c(p, p)))
        stop("Dimensions of the dataset and fixedGaps do not agree.")
    else if (!identical(fixedGaps, t(fixedGaps)))
        stop("fixedGaps must be symmetric")
    else G <- !fixedGaps
    diag(G) <- FALSE
    if (any(is.null(fixedEdges))) {
        fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    }
    else if (!identical(dim(fixedEdges), c(p, p)))
        stop("Dimensions of the dataset and fixedEdges do not agree.")
    else if (!identical(fixedEdges, t(fixedEdges)))
        stop("fixedEdges must be symmetric")


    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list", p))
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)


    while (!done && any(G) && ord <= m.max)
    {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose){
          cat("Order=", ord, "; remaining edges:", remEdges, "\n", sep = "")
      }

      if (method == "stable") {
          G.l <- split(G, gl(p, p))
      }
      for (i in 1:remEdges) {
          if (verbose && (verbose >= 2 || i%%100 == 0)){
            cat("|i=", i, "|iMax=", remEdges, "\n")
          }
          x <- ind[i, 1]
          y <- ind[i, 2]

          if (G[y, x] && !fixedEdges[y, x]) {
            nbrsBool <- if (method == "stable") G.l[[x]]
                        else G[, x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)

            if (length_nbrs >= ord) {
              if (length_nbrs > ord){
                done <- FALSE
              }
              S <- seq_len(ord)

              repeat {
                n.edgetests[ord1] <- n.edgetests[ord1] + 1
# Ronja
                pval <- gaussCItestMI(x, y, nbrs[S], data)

                if (verbose){
                  cat("x=", x, " y=", y, " S=", nbrs[S],": pval =", pval, "\n")
                }
                if (is.na(pval)){
                  pval <- as.numeric(NAdelete)
                }
                if (pMax[x, y] < pval){
                  pMax[x, y] <- pval
                }

                if (pval >= alpha) {
                  G[x, y] <- G[y, x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                } else {
                  nextSet <- pcalg::getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast) break
                  S <- nextSet$nextSet
                }
              }
            }
          }
      }
      ord <- ord + 1L
  }
  for (i in 1:(p - 1)){
    for (j in 2:p){
      pMax[i, j] <- pMax[j, i] <- max(pMax[i,j], pMax[j, i])
    }
  }

 Gobject <- if (sum(G) == 0){methods::new("graphNEL", nodes = labels)}
            else {
              colnames(G) <- rownames(G) <- labels
              methods::as(G, "graphNEL")
            }
 methods::new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}





