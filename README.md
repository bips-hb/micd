# Multiple Imputation for Causal Graph Discovery (micd)

Add-on to the R package `pcalg` for handling missing data in contrataint-based causal graph discovery. Supports continuous, discrete and mixed data. Two options are available: 1) `gaussCItwd`, `disCItwd` and `mixedCItwd` perform test-wise deletion, where missing observations are deleted as necessary on a test-by-test basis; 2) `gaussMItest`, `disMItest` and `mixedMItest` perform conditional independence tests on multiply imputed data. Each of these functions can be used as a plug-in to  `pcalg::skeleton`, `pcalg::pc` or `pcalg::fci`.

Additionally, the package contains all functions required for replicating the analyses in Foraita et al. (2020).

### Installation 
Install the packages `graph` and `RBGL` from Bioconductor. 
Make sure that rtools40 is installed on your computer. 
Then, simply type in R to install `micd`

```R
devtools::install_github("bips-hb/micd")
```
(Note: `micd` was created using R 4.0.3)

### References

Ronja Foraita, Juliane Friemel, Kathrin Günther, Thomas Behrens, Jörn Bullerdiek, Rolf Nimzyk, Wolfgang Ahrens, Vanessa Didelez (2020). Causal discovery of gene regulation with incomplete data. Journal of the Royal Statistical Society: Series A (Statistics in Society), 183(4), 1747-1775. URL https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12565.

### References for pcalg

Markus Kalisch, Martin Mächler, Diego Colombo, Marloes H. Maathuis, Peter Bühlmann (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47(11), 1-26. URL www.jstatsoft.org/v47/i11/.

Alain Hauser, Peter Bühlmann (2012). Characterization and greedy learning of interventional Markov equivalence classes of directed acyclic graphs. Journal of Machine Learning Research, 13, 2409-2464. URL www.jmlr.org/papers/v13/hauser12a.html.

