% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussMItest.R
\name{gaussMItest}
\alias{gaussMItest}
\alias{gaussCItestMI}
\title{Test Conditional Independence of Gaussians via Fisher's Z Using
Multiple Imputations}
\usage{
gaussMItest(x, y, S, suffStat)

gaussCItestMI(x, y, S = NULL, data)
}
\arguments{
\item{x, y, S}{(Integer) position of variable X, Y and set of variables S,
respectively, in the adjacency matrix.
It is tested, whether X and Y are conditionally independent
given the subset S of the remaining nodes.}

\item{suffStat}{A list of length m+1, where m is the number of imputations;
the first m elements are the covariance matrices of the m
imputed data sets, the m-th element is the sample size. Can
be obtained from a mids object by
getSuff(mids, test="gaussMItest")}

\item{data}{An object of type mids, which stands for 'multiply imputed
data set', typically created by a call to function mice()}
}
\value{
A p-value.
}
\description{
A modified version of \code{pcalg::\link[pcalg:condIndFisherZ]{gaussCItest}},
to be used within
\code{pcalg::\link[pcalg]{skeleton}}, \code{pcalg::\link[pcalg]{pc}} or
\code{pcalg::\link[pcalg]{fci}} when multiply imputated data sets are available.
}
\details{
\code{gaussMItest} is faster, as it uses pre-calculated covariance matrices.
}
\examples{
## load data (numeric variables)
dat <- as.matrix(windspeed)

## delete some observations
set.seed(123)
dat[sample(1:length(dat), 260)] <- NA

## Impute missing values under normal model
imp <- mice(dat, method = "norm", printFlag = FALSE)

## analyse data
# complete data:
suffcomplete <- getSuff(windspeed, test = "gaussCItest")
gaussCItest(1, 2, c(4,5), suffStat = suffcomplete)
# multiple imputation:
suffMI <- getSuff(imp, test = "gaussMItest")
gaussMItest(1, 2, c(4,5), suffStat = suffMI)
gaussCItestMI(1, 2, c(4,5), data = imp)
# test-wise deletion:
gaussCItwd(1, 2, c(4,5), suffStat = dat)
# list-wise deletion:
dat2 <- dat[complete.cases(dat), ]
sufflwd <- getSuff(dat2, test = "gaussCItest")
gaussCItest(1, 2, c(4,5), suffStat = sufflwd)

## use gaussMItest or gaussCItestMI within pcalg::pc
(pc.fit <- pc(suffStat = suffMI, indepTest = gaussMItest, alpha = 0.01, p = 6))
(pc.fit <- pc(suffStat = imp, indepTest = gaussCItestMI, alpha = 0.01, p = 6))

}
