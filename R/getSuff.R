#' Obtain 'suffStat' for conditional independence testing
#'
#' A convenience function for transforming a (multiply imputed) data set into the 'suffStat' required
#' by \code{pcalg::\link[pcalg]{gaussCItest}}, \code{pcalg::\link[pcalg]{disCItest}},
#' \code{\link{mixCItest}}, \code{\link{gaussCItwd}}, \code{\link{disCItwd}},
#' \code{\link{mixCItwd}}, \code{\link{gaussMItest}}, \code{\link{disMItest}} or
#' \code{\link{mixMItest}}
#'
#' @param X  for \code{test="xxxCItest"} or \code{test="xxxCItwd"}: a \code{data.frame} or \code{matrix};
#' for \code{test="xxxMItest"}: an object of class \code{mice::\link[mice:mids-class]{mids}},
#' or a list of \code{data.frame}s containing the multiply imputed data sets.
#' @param test one of \code{"gaussCItest"}, \code{"gaussCItwd"}, \code{"gaussMItest"}, \code{"disCItest"},
#' \code{"disCItwd"}, \code{"disMItest"}, \code{"mixCItest"}, \code{"mixCItwd"}, \code{"mixMItest"}
#'
#' @return An R object that can be used as input to the specified conditional independence test
#'
#' @examples
#'
#' ## Example 1: continuous variables, no missing values
#' # load data
#' dat1 <- windspeed
#'
#' # analyse data
#' gaussCItest(1, 2, c(4,5), suffStat = getSuff(dat1, test = "gaussCItest"))
#' gaussCItest(1, 2, c(4,5), suffStat = list(C = cor(dat1), n = nrow(dat1)))
#'
#' ## Example 2: continuous variables, multiple imputation
#' # load data
#' dat2 <- as.matrix(windspeed)
#'
#' # delete some observations
#' set.seed(123)
#' dat2[sample(1:length(dat2), 260)] <- NA
#'
#' # impute missing values under normal model
#' imp2 <- mice(dat2, method = "norm")
#'
#' # analyse imputed data
#' gaussMItest(1, 2, c(4,5), suffStat = getSuff(imp2, test="gaussMItest"))
#' gaussMItest(1, 2, c(4,5), suffStat = 
#' c(lapply(complete(imp2, action = "all"), cor), n = nrow(imp2[[1]])))
#'
#' ## Example 3: discrete variables, multiple imputation
#' # load data
#' data(gmD)
#' dat3 <- gmD$x[1:200, ]
#' dat3[] <- lapply(dat3, as.factor)
#'
#' # delete some observations of X2 and X3
#' set.seed(123)
#' dat3[sample(1:200, 40), 2] <- NA
#' dat3[sample(1:200, 40), 3] <- NA
#'
#' # impute missing values under saturated model
#' form <- make.formulas.saturated(dat3)
#' imp3 <- mice(dat3, formulas = form)
#'
#' # analyse imputed data
#' disMItest(1, 3, NULL, suffStat = getSuff(imp3, test="disMItest"))
#' disMItest(1, 3, NULL, suffStat = complete(imp3, action = "all"))
#'
#' ## Example 4: mixed variables, multiple imputation
#' # load data (numeric and factor variables)
#' dat4 <- toenail2[1:400, ]
#'
#' # delete some observations
#' set.seed(123)
#' dat4[sample(400, 20), 2] <- NA
#' dat4[sample(400, 30), 4] <- NA
#'
#' # impute missing values using random forests
#' imp4 <- mice(dat4, method = "rf")
#' 
#' # analyse imputed data
#' mixMItest(2, 3, 5, suffStat = getSuff(imp4, test = "mixMItest"))
#' mixMItest(2, 3, 5, suffStat = complete(imp4, action = "all"))
#'
#' @export





# a function for conveniently obtaining the 'suffStat' required for the different conditional independence tests

getSuff <- function(X, test=c("gaussCItest", "gaussCItwd", "gaussMItest",
                              "disCItest", "disCItwd", "disMItest",
                              "mixCItest", "mixCItwd", "mixMItest"), adaptDF=NULL, nlev=NULL) {
  
  if ( test %in% c("gaussCItest", "gaussCItwd") ) {
    C <- cor(X)
    n <- nrow(X)
    list(C=C,n=n)
  } else if (test=="gaussMItest") {
    if (class(X)=="mids") {
      X <- mice::complete(X, action="all")
    }
    C <- lapply(X, cor)
    n <- nrow(X[[1]])
    c(C,n)
  } else if ( test %in% c("disCItest","disCItwd") ) {
    if (is.null(adaptDF)) { stop("'adaptDF' needs to be specified. See ?pcalg::disCItest") }
    if (!is.null(nlev)) { if (length(nlev)!=ncol(X)) {stop("Something is wrong with nlev. Check ?pcalg::disCItest")} }
    for (i in 1:ncol(X)) {
      X[ ,i] <- factor(X[ ,i])
      X[ ,i] <- as.numeric(X[ ,i]) - 1
    }
    if (is.null(nlev)) {
      return(list(dm=as.matrix(X), adaptDF=adaptDF))
    } else {
      return(list(dm=as.matrix(X), nlev=nlev, adaptDF=adaptDF))
    }
  } else if ( test=="disMItest" ) {
    if (class(X)=="mids") {
      X <- mice::complete(X, action="all")
    }
    X
  } else if ( test %in% c("mixCItest", "mixCItwd") ) {
    X
  }else if ( test=="mixMItest" ) {
    if (class(X)=="mids") {
      X <- mice::complete(X, action="all")
    }
    X
  }
}

