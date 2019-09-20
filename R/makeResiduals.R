#' Generate residuals based on variables in imputed data sets
#'
#' @param x mids object
#' @param v variable of integers referring to the location in the data set
#' @param confounder variable of integers referring to the location in the data set
#' @param scaled whether the variables should be scaled; default is TRUE.
#'
#' @return
#' @export
#'
#' @examples
#' library(mice)
#' daten <- windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#' colnames(daten) <- names(windspeed)
#'
#' ## Impute missing values
#' imp <- mice(daten, m = 10)
#'
#' ## Build residuals
#' knoten <- 1:4
#' confounder <- 5:6
#' residuals <- list(data = list(), m = 10)
#' for (i in 1:10){
#'  data.i <- mice::complete(imp, i)
#'  residuals$data[[i]] <- makeResiduals(data.i,
#'                           v = knoten, confounder = confounder)}
#'
#' pc.res <- pcMI(data = residuals, p = ncol(residuals$data[[1]]), alpha = 0.05,
#'                maj.rule = TRUE, solve.confl = TRUE)
#' fci.res <- fciMI(data = imp, p = ncol(residuals$data[[1]]), alpha = 0.05)
#'
#' par(mfrow = c(1,2))
#' plot(pc.res)
#' plot(fci.res)
#'
makeResiduals <- function(x, v, confounder, scaled = TRUE)
{
  labels <- colnames(x)
  formeln <- paste0(labels[v],' ~ ', paste(labels[confounder], collapse = ' + '))

  daten <- matrix(ncol = length(v), nrow = nrow(x))

  for(node in 1:length(v))
  {
    daten[, node] <- lm(as.formula(formeln[node]), data = x)$residuals
  }
  colnames(daten) <- labels[v]
  if(scaled == TRUE){daten <- scale(daten)}
}


