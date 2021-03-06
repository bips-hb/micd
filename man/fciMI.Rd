% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fciMI.R
\name{fciMI}
\alias{fciMI}
\title{Estimate a PAG by the FCI-MI Algorithm for multiple
imputed data sets of continuous data}
\usage{
fciMI(
  data,
  alpha,
  labels,
  p,
  skel.method = c("stable", "original"),
  type = c("normal", "anytime", "adaptive"),
  fixedGaps = NULL,
  fixedEdges = NULL,
  NAdelete = TRUE,
  m.max = Inf,
  pdsep.max = Inf,
  rules = rep(TRUE, 10),
  doPdsep = TRUE,
  biCC = FALSE,
  conservative = FALSE,
  maj.rule = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{data}{An object of type mids, which stands for 'multiply imputed
data set', typically created by a call to function mice()}

\item{alpha}{Significance level (number in (0,1) for the conditional
independence tests}

\item{labels}{(Optional) character vector of variable (or "node") names.
Typically preferred to specifying p.}

\item{p}{(Optional) number of variables (or nodes). May be specified if
labels are not, in which case labels is set to 1:p.}

\item{skel.method}{Character string specifying method; the default, "stable"
provides an order-independent skeleton, see
\code{pcalg::\link[pcalg]{skeleton}} for details.}

\item{type}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{fixedGaps}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{fixedEdges}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{NAdelete}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{m.max}{Maximum size of the conditioning sets that are considered in
the conditional independence tests.}

\item{pdsep.max}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{rules}{Logical vector of length 10 indicating which rules should be
used when directing edges. The order of the rules is
taken from Zhang (2008).}

\item{doPdsep}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{biCC}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{conservative}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{maj.rule}{See \code{pcalg::\link[pcalg]{fci}} for details.}

\item{verbose}{If true, more detailed output is provided.}
}
\value{
See \code{pcalg::\link[pcalg]{fci}} for details.
}
\description{
This function is a modification of \code{pcalg::\link[pcalg]{fci}}
to be used for multiple imputation.
}
\examples{

daten <- windspeed[,1]
for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
daten[sample(1:length(daten), 260)] <- NA
daten <- matrix(daten, ncol = 6)

## Impute missing values
imp <- mice(daten)
fciMI(data = imp, label = colnames(imp$data), alpha = 0.01)

}
