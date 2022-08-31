#' Bootstrap Resampling for the PC-MI- and the FCI-MI-algorithm
#'
#' Generate R bootstrap replicates for the PC or FCI algorithm for data with
#' missing values.
#'
#'
#' @param data   Data.frame with missing values
#' @param select Variable of integers, indicating columns to select from a data frame;
#'               only continuous variables can be included in the model selection
#' @param method Character string specifying the algorithm for causal discovery
#'               from the package 'pcalg'.
#' @param method.mice Character string specifying imputation method; see [mice::mice()] for 
#'               more information.              
#' @param args   Arguments passed to `method`. NOTE: argument `labels` is set
#'               internally and should not be used!
#' @param R      A positive integer number of bootstrap replications.
#' @param m      Number of chains included in mice()`.
#' @param args.residuals (Optional) list containing vertices and confounders.
#'               May be specified when residuals for vertices should be calculated in each bootstrap
#'               data set. See [makeResiduals()] for more information
#' @param seed   A positive integer that is used as argument for set.seed().
#' @param quickpred   If true, mice uses quickpred to select predictors.
#' @param ...    Further arguments passed to the imputation function `mice()`.
#'
#' @return List of objects of class `pcalgo` (see [pcalg::pcAlgo-class])
#'         or of `fcmialgo` (see [pcalg::fciAlgo-class]).
#' 
#' @export
#' @examples
#' data(windspeed)
#' daten <- mice::ampute(windspeed)$amp
#'
#' \dontrun{
#' bgraph <- boot.graph(data = daten,
#'                      method = "pcMI",
#'                      args = "solve.confl = TRUE, alpha = 0.05",
#'                      R = 5)
#' }

boot.graph <- function(data, select = NULL, method = c("pcMI", "fciMI"),
                       method.mice = NULL,
                       args, R, m = 10, args.residuals = NULL,
                       seed = NA, quickpred = FALSE, ...)
{
  
    if(!is.data.frame(data)) stop("Data must be a data.frame object.")
    if (!is.na(seed)) set.seed(seed)
    if(R < 1) stop("R must be larger than 0.")
    call <- match.call()
  
    n <- nrow(data)
    samples <- vector(mode = "list", length = R)
    graphs <- vector(mode = "list", length = R)

    for (r in 1:R) samples[[r]] <- sample(1:n, n, replace = TRUE)
    for (g in 1:R) {
      cat(g, "of", R, "replications. \n")      
        if (!quickpred) {
            data.imp <- mice::mice(data[samples[[g]], ], method = method.mice, 
                                   m = m, printFlag = FALSE)
        }
        else {
            predictors <- mice::quickpred(data)
            data.imp <- mice::mice(data[samples[[g]], ], m = m, method = method.mice,
                                   pred = predictors, printFlag = FALSE)
        }
        if (!is.null(select)) {
            data.imp <- getsubsetcol(data.imp, var = select)
        }
        if (!is.null(args.residuals) == TRUE) {
            data.compl <- mice::complete(data.imp, "all", include = TRUE)
            data.res <- list(data = list(), 
                             original = data.compl[[1]],
                             m = m)
            for (ketten in 1:m) {
                data.res$data[[ketten]] <- makeResiduals(data.compl[[ketten +
                  1]], v = args.residuals$v, confounder = args.residuals$conf)
            }
            graphs[[g]] <- eval(parse(text = paste(method, "(data = data.res,",
                args, ", labels = colnames(data.res$data[[1]]))",
                sep = "")))
        }
        else {
            graphs[[g]] <- eval(parse(text = paste(method, "(data = data.imp,",
                args, ",\n                       labels = names(data.imp$imp))",
                sep = "")))
        }
    }
    list(graphs = graphs, call = call)
}


# extract subset out of the imputed data
getsubsetcol <- function(daten, var)
{
  newlist <- list()
  newlist$data <- daten$data[,var]
  newlist$imp <- daten$imp[var]
  newlist$m <- daten$m
  newlist$where <- daten$where[, var]
  newlist$nmis <- daten$nmis[var]

  oldClass(newlist) <- "mids"
  return(newlist)
}
