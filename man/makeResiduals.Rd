% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeResiduals.R
\name{makeResiduals}
\alias{makeResiduals}
\title{Generate residuals based on variables in imputed data sets}
\usage{
makeResiduals(data, v, confounder, method = c("res", "cc", "pd"))
}
\arguments{
\item{data}{A data.frame.}

\item{v}{Vector of integers referring to the location of the variable(s) in the data set}

\item{confounder}{Vector of integers referring to the location of the variable(s)
in the data set (confounders are not included in the network!)}

\item{method}{Default method 'res' uses residuals, 'cc' uses complete cases
and 'pd' uses pairwise deletion}
}
\value{
A data matrix of residuals.
}
\description{
Generate residuals based on variables in imputed data sets
}
\examples{
data(windspeed)
daten <- mice::ampute(windspeed)$amp

# Impute missing values
imp <- mice(daten, m = 5)

# Build residuals
knoten <- 1:4
confounder <- 5:6

# Residuals based on dataset with missing values
res.pd <- makeResiduals(daten, v = knoten, confounder = confounder, method = "pd")

# Residuals based in multiple imputed data
residuals <- list(data = list(), m = 5)
imp_c <- mice::complete(imp, "all")
for (i in 1:imp$m){
   residuals$data[[i]] <- makeResiduals(imp_c[[i]],
                          v = knoten, confounder = confounder)
 }

pc.res <- pcMI(data = residuals, p = length(knoten), alpha = 0.05)
fci.res <- fciMI(data = imp, p = length(knoten), alpha = 0.05)

if (requireNamespace("Rgraphviz", quietly = TRUE)){
oldpar <- par(mfrow = c(1,2)) 
  plot(pc.res)
  plot(fci.res)
par(oldpar)
}

}
