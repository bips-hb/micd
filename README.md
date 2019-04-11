# Multiple Imputation in Causal Graph Discovery (micd)

R functions based on the R package `pcalg` to handle multiple imputated data sets in causal graph discovery. 
At the moment, only PC and FCI algorithm for continuous data are supported. 

### Installation 
Install the packages `graph` and `RBGL` from Bioconductor. 
Then, simply type in R to install `micd`

```R
devtools::install_github("bips-hb/micd")
```

### References for pcalg

Markus Kalisch, Martin Maechler, Diego Colombo, Marloes H. Maathuis, Peter Buehlmann (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47(11), 1-26. URL http://www.jstatsoft.org/v47/i11/.

Alain Hauser, Peter Buehlmann (2012). Characterization and greedy learning of interventional Markov equivalence classes of directed acyclic graphs. Journal of Machine Learning Research, 13, 2409-2464. URL http://jmlr.org/papers/v13/hauser12a.html.
