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
#' @return A list object of S3 class \code{mice::\link[mice]{mira-class}}.
#' @export
#'
#' @examples
#' dat <- as.matrix(windspeed)
#' 
#' ## delete some observations
#' set.seed(123)
#' dat[sample(1:length(dat), 260)] <- NA
#' 
#' ## Impute missing values under normal model
#' imp <- mice(dat, method = "norm")
#' 
#' out.pc <- with_graph(data = imp, 
#'                      algo = "pc", 
#'                      args = ", indepTest = gaussCItest, solve.confl = TRUE,
#'                      labels = names(imp$imp), alpha = 0.05")
#' 
#' out.fci <- with_graph(data = imp, 
#'                       algo = "fciPlus", 
#'                       args = ", indepTest = gaussCItest,
#'                       labels = names(imp$imp), alpha = 0.01")
#'                           
#'  out.ges <- with_graph(data = imp, algo = "ges", arg = NULL, score = TRUE)
#'  
#'  \dontrun{
#'  par(mfrow = c(1,3))
#'  plot(out.pc$res[[1]])
#'  plot(out.fci$res[[1]])
#'  plot(out.ges$res[[1]]$essgraph)
#'  par(mfrow = c(1,1))
#'  }
#'                                      
with_graph <- function(data, algo = c("pc", "fci", "fciPlus", "ges"), 
                       args, score = FALSE)
{
  call <- match.call()
  res <- vector(mode = "list", length = data$m)

  for (i in seq_along(res)) {
    data.i <- mice::complete(data, i)
 
    if(score)
    {
      s <- methods::new("GaussL0penObsScore", data.i)
      res[[i]] <- eval(parse(text = paste(algo,"(s", args, ")", sep="")))

    } else {
      suff <- list(C = stats::cor(data.i), n = nrow(data.i))
      res[[i]] <- eval(parse(text = paste(algo,"(suff", args, ")", sep="")))
    }
  }
  object <- list(call = call, call1 = data$call, nmis = data$nmis, 
                 res = res)
  oldClass(object) <- c("mira", "matrix")
  return(object)
}
  