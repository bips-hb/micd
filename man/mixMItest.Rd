% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixMItest.R
\name{mixMItest}
\alias{mixMItest}
\title{Likelihood Ratio Test for (Conditional) Independence between Mixed Variables
after Multiple Imputation}
\usage{
mixMItest(x, y, S = NULL, suffStat, moreOutput = FALSE)
}
\arguments{
\item{x, y, S}{(integer) position of variable X, Y and set of variables S,
respectively, in \code{suffStat}. It is tested whether X and Y are
conditionally independent given the subset S of the remaining variables.}

\item{suffStat}{A list of \code{data.frame}s containing the multiply
imputed data sets. Usually obtained from a \code{mice::\link[mice:mids-class]{mids}}
object using \code{mice::\link[mice:complete.mids]{complete}} with argument \code{action="all"}.
Discrete variables must be coded as factors.}

\item{moreOutput}{(only for mixed of discrete variables) If \code{TRUE}, the test
statistic, its main components and
the degrees of freedom are returned in addition to the p-value. Defaults to
\code{FALSE}.}
}
\value{
A p-value. If \code{moreOutput=TRUE}, the test statistic, its main
components and the degrees of freedom are returned as well.
}
\description{
A modified version of \code{\link{mixCItest}}, to be used within \code{pcalg::\link[pcalg]{skeleton}},
\code{pcalg::\link[pcalg]{pc}} or \code{pcalg::\link[pcalg]{fci}} when multiply imputed data sets are available.
}
\details{
See \code{\link{mixCItest}} for details on the assumptions of the
Conditional Gaussian likelihood ratio test. \code{CGtestMI} applies this test
to each \code{data.frame} in \code{suffStat}, then combines the results using
the rules in Meng & Rubin (1992).
}
\examples{

## load data (numeric and factor variables)
data(toenail2)
dat <- toenail2[1:1000, ]

## delete some observations
set.seed(123)
dat[sample(1000, 20), 2] <- NA
dat[sample(1000, 30), 4] <- NA

## impute missing values using random forests (because of run time we just impute 2 chains)
imp <- mice(dat, method = "rf", m = 2, printFlag = FALSE)

## analyse data
# complete data:
mixCItest(2, 3, 5, suffStat = toenail2[1:1000, ])
# multiple imputation:
suffMI <- complete(imp, action = "all")
mixMItest(2, 3, 5, suffStat =  suffMI)
# test-wise deletion:
mixCItwd(2, 3, 5, suffStat = dat)
# list-wise deletion:
sufflwd <- dat[complete.cases(dat), ]
mixCItest(2, 3, 5, suffStat = sufflwd)

## use mixMItest within pcalg::pc
\donttest{
pc.fit <- pc(suffStat =  suffMI, indepTest = mixMItest, alpha = 0.01, p = 5)
pc.fit
}

}
\references{
Meng X.-L., Rubin D.B. (1992): Performing likelihood ratio tests
with multiply imputed data sets. \emph{Biometrika} 79(1):103-111.
}
\author{
Janine Witte
}
