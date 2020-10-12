#' Creates a \code{formulas} Argument
#' 
#' This helper function creates a valid \code{formulas} object. 
#' The \code{formulas} object is an argument to the \code{mice::\link[mice]{mice}} function. 
#' It is a list of formulas that specifies the target variables and the predictors 
#' by means of the standard ~ operator. In contrast to \code{mice::\link[mice]{make.formulas}}, 
#' which creates main effects formulas, \code{make.formulas.saturated} 
#' creates formulas including all possible interaction effects.
#'
#' @param data  A \code{data.frame} with the source data.
#' @param blocks  An optional specification for blocks of variables in the rows. 
#' The default assigns each variable in its own block.
#' @param predictorMatrix  A \code{predictorMatrix} specified by the user.
#'
#' @return A list of formulas.
#' 
#' @note A modification of \code{mice::\link[mice]{make.formulas}} by Stef van Buuren et al.
#' 
#' @seealso \code{mice::\link[mice]{make.formulas}}
#' 
#' @examples
#' ## main effecs model:
#' f1 <- make.formulas(nhanes)
#' f1
#' 
#' ## saturated model:
#' f2 <- make.formulas.saturated(nhanes)
#' f2
#' 
#' @export


make.formulas.saturated <- function (data, blocks = mice::make.blocks(data), 
                                     predictorMatrix = NULL) {
  #data <- mice:::check.dataform(data)
  data <- mice_check_dataform(data)
  formulas <- as.list(rep("~ 0", length(blocks)))
  names(formulas) <- names(blocks)
  for (h in names(blocks)) {
    y <- blocks[[h]]
    if (is.null(predictorMatrix)) {
      predictors <- colnames(data)
    }
    else {
      type <- predictorMatrix[h, ]
      predictors <- names(type)[type != 0]
    }
    x <- setdiff(predictors, y)
    formulas[[h]] <- paste(paste(y, collapse = "+"), "~ 0 +", 
                           paste(x, collapse = ":"))
  }
  formulas <- lapply(formulas, stats::as.formula)
  formulas
}


# hidden function from the mice package
mice_check_dataform <- function (data) 
{
  if (!(is.matrix(data) || is.data.frame(data))) 
    stop("Data should be a matrix or data frame", call. = FALSE)
  if (ncol(data) < 2) 
    stop("Data should contain at least two columns", 
         call. = FALSE)
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  df <- sapply(data, is.data.frame)
  if (any(mat)) 
    stop("Cannot handle columns with class matrix: ", 
         colnames(data)[mat])
  if (any(df)) 
    stop("Cannot handle columns with class data.frame: ", 
         colnames(data)[df])
  dup <- duplicated(colnames(data))
  if (any(dup)) 
    stop("Duplicate names found: ", paste(colnames(data)[dup], 
                                          collapse = ", "))
  data
}
