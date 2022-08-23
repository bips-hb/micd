#' Generate residuals based on variables in imputed data sets
#'
#' @param data a data frame
#' @param v vector of integers referring to the location of the variable(s) in the data set
#' @param confounder vector of integers referring to the location of the variable(s) 
#' in the data set (confounders are not included in the network!)
#' @param method default method 'res' uses residuals, 'cc' uses complete cases
#'               and 'pd' uses pairwise deletion   
#'               
#' @return A data matrix of residuals.                                       
#'               
#' @importFrom stats na.omit as.formula lm
#'
#' @export
#'
#' @examples
#' data(windspeed)
#' daten <- mice::ampute(windspeed)$amp
#' 
#' # Impute missing values
#' imp <- mice(daten, m = 5)
#'
#' # Build residuals
#' knoten <- 1:4
#' confounder <- 5:6
#'
#' # Residuals based on dataset with missing values
#' res.pd <- makeResiduals(daten, v = knoten, confounder = confounder, method = "pd")
#'
#' # Residuals based in multiple imputed data
#' residuals <- list(data = list(), m = 5)
#' imp_c <- mice::complete(imp, "all")
#' for (i in 1:imp$m){
#'    residuals$data[[i]] <- makeResiduals(imp_c[[i]],
#'                           v = knoten, confounder = confounder)
#'  }
#'
#' pc.res <- pcMI(data = residuals, p = length(knoten), alpha = 0.05)
#' fci.res <- fciMI(data = imp, p = length(knoten), alpha = 0.05)
#'
#' if(require("Rgraphviz", character.only = TRUE, quietly = TRUE)){
#' par(mfrow = c(1,2))
#'   plot(pc.res)
#'   plot(fci.res)
#' par(mfrow = c(1,1))
#' }
#'
makeResiduals <- function (data, v, confounder, method = c("res","cc", "pd"))
{
  deletion.method <- match.arg(method)
  if(any(is.na(data[, c(v, confounder)])) & deletion.method == "res" )
     stop("Dataset contains missing values. Select method 'cc' or 'pd'.")
  
  labels <- colnames(data)
  rownames(data) <- 1:nrow(data)
  formeln <- paste0(labels[v], " ~ ", paste(labels[confounder], collapse = " + "))
  # browser()
  if (deletion.method == "cc"){data <- stats::na.omit(data)}
  daten <- matrix(ncol = length(v), nrow = nrow(data))
  
  for (node in 1:length(v)) {
    
    if (deletion.method == "pd") {
      tmp <- stats::lm(stats::as.formula(formeln[node]), data = data)$residuals
      daten[as.numeric(names(tmp)), node] <- tmp
    }
    else {
      daten[, node] <- stats::lm(stats::as.formula(formeln[node]), data = data)$residuals
    }
  }
  colnames(daten) <- labels[v]
  return(daten)
}
