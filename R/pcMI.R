#' Estimate the Equivalence Class of a DAG using the PC-MI algorithm for multiple
#' imputed data sets of continuous data
#'
#' This function is a modification of [functioname(pcalg::pc)]
#' to be used for multiple imputation.
#'
#' @param data  An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#' @param alpha Significance level (number in (0,1) for the conditional
#'              independence tests
#' @param labels (Optional) character vector of variable (or "node") names.
#'               Typically preferred to specifying p.
#' @param p  (Optional) number of variables (or nodes). May be specified if
#'            labels are not, in which case labels is set to 1:p.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i]
#'                  (or both) are TRUE, the edge i-j is removed before starting
#'                  the algorithm. Therefore, this edge is guaranteed to be
#'                  absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i]
#'                  (or both) are TRUE, the edge i-j is never considered for
#'                  removal. Therefore, this edge is guaranteed to be present
#'                  in the resulting graph
#' @param NAdelete If indepTest returns NA and this option is TRUE,
#'                 the corresponding edge is deleted. If this option is FALSE,
#'                 the edge is not deleted.
#' @param m.max Maximal size of the conditioning sets that are considered in the
#'              conditional independence tests.
#' @param u2pd  String specifying the method for dealing with conflicting
#'              information when trying to orient edges (see details below).
#' @param skel.method Character string specifying method; the default, "stable"
#'                provides an order-independent skeleton, see
#'                [functioname(pcalg::skeleton)] for details.
#' @param conservative Logical indicating if the conservative PC is used. See
#'                [functioname(pcalg::pc)] for details.
#' @param maj.rule Logical indicating that the triples shall be checked for
#'                 ambiguity using a majority rule idea, which is less strict
#'                 than the conservative PC algorithm. For more information, see
#'                 [functioname(pcalg::pc)].
#' @param solve.confl See [functioname(pcalg::pc)] for more details.
#' @param verbose If TRUE, detailed output is provided.
#'
#' @return An object of class "pcAlgo" (see pcAlgo) containing an estimate of
#'         the equivalence class of the underlying DAG.
#'
#' @return See [functioname(pcalg::pc)] for more details.
#'
#' @note This is a modified function of [functioname(pcalg::pc)]
#'       from the package 'pcalg' (Kalisch et al., 2012;
#'       http://www.jstatsoft.org/v47/i11/).
#'
#' @useDynLib micd
#' @importFrom Rcpp sourceCpp
#'
#' @export
#'
#' @examples
#' library(mice)
#' daten <- windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#'
#' ## Impute missing values
#' imp <- mice(daten)
#' pcMI(data = imp, label = colnames(imp$data), alpha = 0.01)
#'
pcMI <- function (data, alpha, labels, p, fixedGaps = NULL,
        fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, u2pd = c("relaxed",
        "rand", "retry"), skel.method = c("stable", "original"),
        conservative = FALSE, maj.rule = FALSE,
        solve.confl = FALSE, verbose = FALSE)
{
    cl <- match.call()
    if (!missing(p))
        stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
            1, p >= 2)
    if (missing(labels)) {
        if (missing(p))
            stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    }
    else {
        stopifnot(is.character(labels))
        if (missing(p)) {
            p <- length(labels)
        }
        else if (p != length(labels))
          stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
        else message("No need to specify 'p', when 'labels' is given")
    }
    u2pd <- match.arg(u2pd)
    skel.method <- match.arg(skel.method)
    if (u2pd != "relaxed") {
        if (conservative || maj.rule)
            stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
        if (solve.confl)
            stop("Versions of PC using lists for the orientation rules
                 (and possibly bi-directed edges)\n can only be run with
                  'u2pd = relaxed'")
    }
    if (conservative && maj.rule)
        stop("Choose either conservative PC or majority rule PC!")

    skel <- skeletonMI(data, alpha, labels = labels,
        method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
        NAdelete = NAdelete, m.max = m.max, verbose = verbose)
    skel@call <- cl

    if (!conservative && !maj.rule) {
        switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj,
            relaxed = udag2pdagRelaxed(skel, verbose = verbose,
                solve.confl = solve.confl))
    }
    else {
        pc. <- pc.cons.internMI(skel, data, alpha, version.unf = c(2, 1),
                              maj.rule = maj.rule, verbose = verbose)
        udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl,
            solve.confl = solve.confl)
    }
}




















