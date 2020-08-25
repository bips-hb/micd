#' Evaluate causal graph discovery algorithm in multiple imputed datasets
#'
#' @param data An object of type mids, which stands for 'multiply imputed 
#'             data set', typically created by a call to function mice()
#' @param algo An algorithm for causal discovery from the package 'pcalg' 
#'             (see details).
#' @param args Additional arguments passed to the algo. Must be a string
#'             vector starting with comma, i.e. ", ..."
#' @param score Logical indicating whether a score-based or a constrained-based
#'              algorithm is applied.
#'              
#'              
#' @return A list object of S3 class \code{mira}.
#' @export
#'
#' @examples
#' \dontrun{
#' library(mice)
#' imp <- mice(windspeed)
#'  out.pc <- with_graph(data = imp, 
#'                       algo = "pc", 
#'                       args = ", indepTest = gaussCItest, solve.confl = TRUE,
#'                           labels = names(imp$imp), alpha = 0.05")
#' 
#'  out.fci <- with_graph(data = imp, 
#'                        algo = "fciPlus", 
#'                        args = ", indepTest = gaussCItest,
#'                               labels = names(imp$imp), alpha = 0.01")
#'                           
#'  out.ges <- with_graph(data = imp, 
#'                    algo = "ges", arg = NULL, 
#'                                      score = TRUE)
#' }                                      
with_graph <- function(data, algo = c("pc", "fci", "fciPlus", "ges"), 
                       args, score = FALSE)
{
  call <- match.call()
  analyses <- vector(mode = "list", length = data$m)

  for (i in seq_along(analyses)) {
    data.i <- mice::complete(data, i)
 
    if(score)
    {
      s <- methods::new("GaussL0penObsScore", data.i)
      analyses[[i]] <- eval(parse(text = paste(algo,"(s", args, ")", sep="")))

    } else {
      suff <- list(C = stats::cor(data.i), n = nrow(data.i))
      analyses[[i]] <- eval(parse(text = paste(algo,"(suff", args, ")", sep="")))
    }
  }
  object <- list(call = call, call1 = data$call, nmis = data$nmis, 
                 analyses = analyses)
  oldClass(object) <- c("mira", "matrix")
  return(object)
}
  