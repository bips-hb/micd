#' Obtain 'suffStat' for conditional independence testing
#'
#' A convenience function for transforming a multiply imputed data set into the 'suffStat' required
#' by [pcalg::gaussCItest()], [pcalg::disCItest()], [mixCItest()], [flexCItest()], [gaussMItest()], 
#' [disMItest()], [mixMItest()] or [flexMItest()].
#'
#' @param X  for __test="xxxCItest"__: a __data.frame__ or __matrix__;
#' for __test="xxxMItest"__: an object of class \code{mice::\link[mice:mids-class]{mids}},
#' or a list of data.frames containing the multiply imputed data sets.
#' @param test one of [gaussCItest()], [gaussMItest()], [disCItest()],
#' [disMItest()], [mixCItest()], [mixMItest()], [flexCItest()],
#' [flexMItest()].
#' @param adaptDF for discrete variables: logical specifying if the degrees of freedom 
#' should be lowered by one for each zero count. The value for the degrees of freedom 
#' cannot go below 1.
#' @param nlev (Optional) for discrete variables: vector with numbers of levels for each 
#' variable in the data.
#'
#' @return An R object that can be used as input to the specified conditional independence test:
#' 
#' @importFrom Rfast which.is
#' @importFrom mice complete
#' @importFrom stats cor
#'
#' @examples
#'
#' ## Example 1: continuous variables, no missing values =====================
#' ## load data
#' library(mice)
#' dat1 <- as.matrix(windspeed)
#'
#' ## analyse data
#' gaussCItest(1, 2, NULL, suffStat = getSuff(windspeed, test = "gaussCItest"))
#' mixCItest(1, 2, NULL, suffStat = windspeed)
#'
#' ## Example 2: continuous variables, multiple imputation ===================
#' ## load data
#' dat2 <- as.matrix(windspeed)
#'
#' ## delete some observations
#' set.seed(123)
#' dat2[sample(1:length(dat2), 260)] <- NA
#'
#' ## Impute missing values under normal model
#' imp2 <- mice(dat2, method = "norm")
#'
#' ## analyse imputed data
#' gaussMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="gaussMItest"))
#' mixMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="mixMItest"))
#' mixMItest(1, 2, c(4,5), suffStat = mice::complete(imp2, action="all"))
#' flexMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="flexMItest"))
#'
#' ## Example 3: discrete variables, multiple imputation =====================
#' ## simulate factor variables
#' n <- 200
#' set.seed(123)
#' x <- factor(sample(0:2, n, TRUE)) # factor, 3 levels
#' y <- factor(sample(0:3, n, TRUE)) # factor, 4 levels
#' z <- factor(sample(0:1, n, TRUE)) # factor, 2 levels
#' dat3 <- data.frame(x,y,z)
#'
#' ## delete some observations of z
#' dat3[sample(1:n, 40), 3] <- NA
#'
#' ## impute missing values under saturated model
#' form <- make.formulas.saturated(dat3)
#' imp3 <- mice::mice(dat3, method = "logreg", formulas = form)
#'
#' ## analyse imputed data
#' disMItest(1, 3, NULL, suffStat = getSuff(imp3, test="disMItest"))
#' disMItest(1, 3, NULL, suffStat = mice::complete(imp3, action = "all"))
#' mixMItest(1, 3, NULL, suffStat = getSuff(imp3, test="mixMItest"))
#' mixMItest(1, 3, NULL, suffStat = mice::complete(imp3, action = "all"))
#' flexMItest(1, 3, NULL, suffStat = getSuff(imp3, test="flexMItest"))
#'
#' ## Example 4: mixed variables, multiple imputation =========================
#' ## load data (numeric and factor variables)
#' dat4 <- toenail2[1:400, ]
#'
#' ## delete some observations
#' set.seed(123)
#' dat4[sample(400, 20), 2] <- NA
#' dat4[sample(400, 30), 4] <- NA
#'
#' ## impute missing values using random forests
#' imp4 <- mice(dat4, method="rf")
#' mixMItest(2, 3, 5, suffStat = getSuff(imp4, test="mixMItest"))
#' mixMItest(2, 3, 5, suffStat = mice::complete(imp4, action="all"))
#' flexMItest(2, 3, 5, suffStat = getSuff(imp4, test="flexMItest"))
#'
#' @export
getSuff <- function(X, test = c("gaussCItest", "gaussMItest", "disCItest", "disMItest",
                                "disCItwd", "mixCItest", "mixMItest", "flexMItest", "flexCItest"), 
                    adaptDF = NULL, nlev = NULL) {
  
  if (test=="gaussCItest") {
    C <- stats::cor(X)
    n <- nrow(X)
    list(C=C,n=n)
  } else if (test=="gaussMItest") {
    if(!(mice::is.mids(X) | is.list(X))) stop("data is neither a list nor a mids object")
    
    if (inherits(X, "mids")) {
      X <- mice::complete(X, action="all")
      if(length(which.is(X, "factor") > 0)) stop("data must be all numeric.")
    }
    C <- lapply(X, stats::cor)
    n <- nrow(X[[1]])
    c(C,n)
  } else if (test=="disCItest") {
    if (is.null(adaptDF)) { stop("'adaptDF' needs to be specified. See ?pcalg::disCItest") }
    if (!is.null(nlev)) { if (length(nlev) != ncol(X)) {stop("Something is wrong with nlev. Check ?pcalg::disCItest")} }
    for (i in 1:ncol(X)) {
      X[ ,i] <- as.numeric(X[ ,i]) - 1
    }

    if(is.null(adaptDF)) adaptDF = TRUE
    if (is.null(nlev)) {
      return(list(dm=as.matrix(X), adaptDF = adaptDF))
    } else {
      return(list(dm=as.matrix(X), nlev = nlev, adaptDF = adaptDF))
    }
  } else if (test=="disMItest") {
    if (inherits(X, "mids")) {
      X <- mice::complete(X, action="all")
    }
    X
  } else if (test=="disCItwd") {
    if (is.null(adaptDF)) adaptDF <- TRUE
    return(list(dm=as.matrix(X), adaptDF=adaptDF))    
  } else if (test=="mixCItest") {
    X
  }else if (test=="mixMItest") {
    if (inherits(X, "mids")) {
      X <- mice::complete(X, action="all")
    }
    X
  } else if (test=="flexCItest") {
    conpos <- Rfast::which.is(X, "numeric")
    dispos <- Rfast::which.is(X, "factor")
     
    corlist <- NULL
    datlist <- NULL
    
    if (length(conpos)==0) {
      # message("Note: The variables are all discrete (flexCI).\n")
    }
    if (length(dispos)==0) {
      # message("Note: The variables are all numeric (flexCI).\n")
    }
    
    if (length(conpos)>1 ) {
      C <- stats::cor(X[ ,conpos])
      n <- nrow(X)
      corlist <- list(C=C, n=n)
    }
    
    if (length(dispos)>1 ) {
      # if (is.null(adaptDF)) { stop("'adaptDF' needs to be specified. See ?pcalg::disCItest") }
      if (!is.null(nlev)) { if (length(nlev) != ncol(X)) {
        stop("Something is wrong with nlev. Check ?pcalg::disCItest")} }
      
      D <- X[ ,dispos]
      n_level <- vector(length = length(dispos))
      for (i in 1:ncol(D)) {
        n_level[i] <- nlevels(D[,i])
        # pcalg requires to start with 0
        D[ ,i] <- as.integer(D[, i]) - 1
        if(min(D[,i], na.rm = TRUE) != 0) D[,i] - min(D[,i], na.rm = TRUE)
      }
      if(is.null(adaptDF)) adaptDF = TRUE
      if (is.null(nlev)) {
        datlist <- list(dm = D, nlev = n_level, adaptDF = adaptDF)
      } else {
        datlist <- list(dm = D, nlev = nlev, adaptDF = adaptDF)
      }
    }
    
    return(list(dat = X, corlist = corlist, datlist = datlist, conpos = conpos, dispos = dispos))
    
  } else if (test=="flexMItest") {
    if (inherits(X, "mids")) {
      X <- mice::complete(X, action="all")
    }
    conpos <- Rfast::which.is(X[[1]], "numeric")
    dispos <- Rfast::which.is(X[[1]], "factor")
    
    if (length(conpos)==0) {
      # message("Note: The variables are all discrete (flexMI).\n")
    }
    if (length(dispos)==0) {
      # message("Note: The variables are all numeric (flexMI).\n")
    }
    
    corlist <- NULL
    
    if (length(conpos)>1 ) {
      C <- lapply(X, function(i){stats::cor(i[ ,conpos])})
      n <- nrow(X[[1]])
      corlist <- c(C, n)
    }
    
    return(list(datlist=X, corlist=corlist, conpos=conpos, dispos=dispos))
  }
}
