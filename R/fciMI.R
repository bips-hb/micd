#' Estimate a PAG by the FCI-MI Algorithm for multiple
#' imputed data sets of continuous data
#'
#' This function is a modification of [pcalg::fci()]
#' to be used for multiple imputation.
#'
#' @param data An object of type mids, which stands for 'multiply imputed
#'              data set', typically created by a call to function mice()
#' @param alpha Significance level (number in (0,1) for the conditional
#'              independence tests
#' @param labels (Optional) character vector of variable (or "node") names.
#'               Typically preferred to specifying p.
#' @param p (Optional) number of variables (or nodes). May be specified if
#'            labels are not, in which case labels is set to 1:p.
#' @param skel.method Character string specifying method; the default, "stable"
#'                provides an order-independent skeleton, see [pcalg::skeleton()] for details.
#' @param type Character string specifying the version of the FCI algorithm to be used. 
#'              See [pcalg::fci()] for details.
#' @param fixedGaps See [pcalg::fci()] for details.
#' @param fixedEdges See [pcalg::fci()] for details.
#' @param NAdelete See [pcalg::fci()] for details.
#' @param m.max Maximum size of the conditioning sets that are considered in
#'              the conditional independence tests.
#' @param pdsep.max See [pcalg::fci()] for details.
#' @param rules Logical vector of length 10 indicating which rules should be
#'              used when directing edges. The order of the rules is
#'              taken from Zhang (2008).
#' @param doPdsep See [pcalg::fci()] for details.
#' @param biCC See [pcalg::fci()] for details.
#' @param conservative See [pcalg::fci()] for details.
#' @param maj.rule See [pcalg::fci()] for details.
#' @param verbose If true, more detailed output is provided.
#'
#' @return See [pcalg::fci()] for details.
#' 
#' @author Original code by Diego Colombo, Markus Kalisch, and  Joris Mooij.
#' Modifications by Ronja Foraita.   
#' 
#' @importFrom methods as new
#' @export
#' @examples
#' 
#' daten <- windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#'
#' ## Impute missing values
#' imp <- mice(daten)
#' fc.res <- fciMI(data = imp, label = colnames(imp$data), alpha = 0.01)
#' 
#' if(require("Rgraphviz", character.only = TRUE, quietly = TRUE)){
#' plot(fc.res)
#' }
#'
fciMI <- function (data, alpha, labels, p, skel.method = c("stable",
    "original"), type = c("normal", "anytime", "adaptive"), fixedGaps = NULL, 
    fixedEdges = NULL, NAdelete = TRUE,
    m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10), doPdsep = TRUE,
    biCC = FALSE, conservative = FALSE, maj.rule = FALSE,
    verbose = FALSE)
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
    type <- match.arg(type)
    if (type == "anytime" && m.max == Inf)
        stop("To use the Anytime FCI you must specify a finite 'm.max'.")
    if (type == "adaptive" && m.max != Inf)
        stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
    if (conservative && maj.rule)
        stop("Choose either conservative FCI or majority rule FCI")
    cl <- match.call()
    if (verbose)
        cat("Compute Skeleton\n================\n")
     skel <- skeletonMI(data, alpha, labels = labels,
        method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
        NAdelete = NAdelete, m.max = m.max, verbose = verbose)
    skel@call <- cl
    G <- methods::as(skel@graph, "matrix")
    sepset <- skel@sepset
    pMax <- skel@pMax
    n.edgetestsSKEL <- skel@n.edgetests
    max.ordSKEL <- skel@max.ord
    allPdsep <- NA
    tripleList <- NULL
    if (doPdsep) {
        if (verbose)
            cat("\nCompute PDSEP\n=============\n")

        pc.ci <- micd:::pc.cons.internMI(skel, data, alpha = alpha,
            version.unf = c(1, 1), maj.rule = FALSE, verbose = verbose)

        pdsepRes <- micd:::pdsepMI(skel@graph, data, p = p,
                      sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                      m.max = if (type == "adaptive")
                max.ordSKEL
            else m.max, pdsep.max = pdsep.max, NAdelete = NAdelete,
            unfVect = pc.ci$unfTripl, biCC = biCC, verbose = verbose)
        G <- pdsepRes$G
        sepset <- pdsepRes$sepset
        pMax <- pdsepRes$pMax
        allPdsep <- pdsepRes$allPdsep
        n.edgetestsPD <- pdsepRes$n.edgetests
        max.ordPD <- pdsepRes$max.ord
        if (conservative || maj.rule) {
            if (verbose)
                cat("\nCheck v-structures conservatively\n=================================\n")
            tmp.pdsep <- methods::new("pcAlgo", graph = methods::as(G, "graphNEL"),
                call = cl, n = integer(0), max.ord = as.integer(max.ordSKEL),
                n.edgetests = n.edgetestsSKEL, sepset = sepset,
                pMax = pMax, zMin = matrix(NA, 1, 1))

            sk. <- micd:::pc.cons.internMI(tmp.pdsep, data, alpha, verbose = verbose,
                    version.unf = c(1, 1), maj.rule = maj.rule)
            tripleList <- sk.$unfTripl
            sepset <- sk.$sk@sepset
        }
    }
    else {
        n.edgetestsPD <- 0
        max.ordPD <- 0
        allPdsep <- vector("list", p)
        if (conservative || maj.rule) {
            if (verbose)
                cat("\nCheck v-structures conservatively\n=================================\n")

            nopdsep <- micd:::pc.cons.internMI(skel, data, alpha, verbose = verbose,
                         version.unf = c(2,1), maj.rule = maj.rule)
            tripleList <- nopdsep$unfTripl
            sepset <- nopdsep$sk@sepset
        }
    }
    if (verbose)
        cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
            "\nCompute collider:\n")
    res <- pcalg::udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
        verbose = verbose)
    colnames(res) <- rownames(res) <- labels
    methods::new("fciAlgo", amat = res, call = cl, n = integer(0), max.ord = as.integer(max.ordSKEL),
        max.ordPDSEP = as.integer(max.ordPD), n.edgetests = n.edgetestsSKEL,
        n.edgetestsPDSEP = n.edgetestsPD, sepset = sepset, pMax = pMax,
        allPdsep = allPdsep)
}
