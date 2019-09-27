#' Bootstrap Resampling for the PC-MI- and the FCI-MI-algorithm
#'
#' Generate R bootstrap replicates for the PC or FCI algorithm for data with
#' missing values.
#'
#'
#' @param data   Data frame with missing values
#' @param select variable of integers, indicating columns to select from a data frame;
#'               only continous variables can be included in the model selection
#' @param method Character string specifying the algorithm for causal discovery
#'               from the package 'pcalg'.
#' @param args   Arguments passed to `method`. NOTE: argument `labels` is set
#'               internally and should not be used!
#' @param R      A positive integer number of bootstrap replications.
#' @param m      Number of chains included in mice()`.
#' @param args.residuals (optional) list containing vertices and confounders.
#'               May be specified when residuals for vertices should be calculated in each bootstrap
#'               data set. See [functioname(micd::makeResiduals)] for more information
#' @param seed   A positive integer that is used as argument for set.seed().
#' @param quickpred   If true, mice uses quickpred to select predictors.
#' @param ...    Further arguments passed to the imputation function `mice()`.
#'
#' @return List of objects of class `pcalgo` (see [functionname(pcalg::`pcAlgo-class`)])
#'         or of `fcmialgo` (see [functionname(pcalg::`fciAlgo-class`)]).
#' @export
#' @examples
#' library(mice)
#' daten <- windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#'
#' boot.graph <- boot.graph(
#'             data = daten,
#'             method = "pcMI",
#'             args = "solve.confl = TRUE, alpha = 0.05",
#'             R = 10)

boot.graph <- function(data, select = NULL, method = c("pcMI", "fciMI"),
                       args, R, m = 10, args.residuals = NULL,
                       seed = NA, quickpred = FALSE, ...)
{
  if (!is.na(seed)) set.seed(seed)

  n <- nrow(data)

  samples <- vector(mode = "list", length = R)
  graphs <- vector(mode = "list", length = R)

  for(r in 1:R) samples[[r]] <- sample(1:n, n, replace = TRUE)

  for(g in 1:R)
  {
    if(!quickpred){
       data.imp <- mice::mice(data[samples[[g]],], m = m)
    } else {
       predictors <- mice::quickpred(data)
         data.imp <- mice::mice(data[samples[[g]],], m = m, pred = predictors)
    }

    if(!is.null(select)){data.imp <- getsubsetcol(data.imp, var = select)}

    if(!is.null(args.residuals) == TRUE){
      data.compl <- mice::complete(data.imp, "all", include = TRUE)
      data.res <- list(data = list(), original = data.compl[[1]], m = m)

      for (ketten in 1:m){
        data.res$data[[ketten]] <- makeResiduals(data.compl[[ketten + 1]],
                                     v = args.residuals$v,
                                     confounder = args.residuals$conf)}

      graphs[[g]] <- eval(parse(text = paste(method, "(data = data.res,", args,
                        ", labels = colnames(data.res$data[[1]]))", sep = "")))

    } else {
      GE.imp <- getsubsetcol(data.imp, var = c(1:5))
      graphs[[g]] <- eval(parse(text = paste(method, "(data = GE.imp,", args, ",
                       labels = names(GE.imp$imp))", sep="")))
      #graphs[[g]] <- eval(parse(text = paste(method, "(data = data.imp,", args,
      #                  ", labels = names(data.imp$data))", sep = "")))
    }
  }
  list(graphs = graphs, call = call)
}


#' @export
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
