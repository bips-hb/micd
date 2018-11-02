#' Bootstrap Resampling for the PC-MI- and the FCI-MI-algorithm
#'
#' Generate R bootstrap replicates for the PC or FCI algorithm for data with
#' missing values.
#'
#'
#' @param data   Data frame with missing values
#' @param method Character string specifying the algorithm for causal discovery
#'               from the package 'pcalg'.
#' @param args   Arguments passed to `method`. NOTE: argument `labels` is set
#'               internally and should not be used!
#' @param R      A positive integer number of bootstrap replications.
#' @param seed   A positive integer that is used as argument for set.seed().
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

boot.graph <- function(data, method = c("pcMI", "fciMI"), args, R, m = 10,
                       seed = NA, quickpred = FALSE, ...)
{
  if (!is.na(seed)) set.seed(seed)

  n <- nrow(data)

  samples <- vector(mode = "list", length = R)
  graphs <- vector(mode = "list", length = R)

  for(r in 1:R) samples[[r]] <- sample(1:n, n, replace = TRUE)

  for(g in 1:R)
  {
    if(!qp){
       imp <- mice(data[samples[[g]],], m = m, ...)
    } else {
       predictors <- quickpred(data)
       imp <- mice(data[samples[[g]],], m = m, pred = predictors, ...)
    }
    GE.imp <- getsubsetcol(imp, var = c(1:5))
    graphs[[g]] <- eval(parse(text = paste(method, "(data = GE.imp,", args, ",
                       labels = names(GE.imp$imp))", sep="")))
  }
  list(graphs = graphs, call = call)
}


