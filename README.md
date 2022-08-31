# Multiple Imputation for Causal Graph Discovery (micd)

Add-on to the R package `pcalg` for handling missing data in contrataint-based causal graph discovery. Supports continuous, discrete and mixed data. Two options are available: 1) `gaussCItwd`, `disCItwd` and `mixedCItwd` perform test-wise deletion, where missing observations are deleted as necessary on a test-by-test basis; 2) `gaussMItest`, `disMItest` and `mixedMItest` perform conditional independence tests on multiply imputed data. Each of these functions can be used as a plug-in to  `pcalg::skeleton`, `pcalg::pc` or `pcalg::fci`.

Additionally, the package contains all functions required for replicating the analyses in Foraita et al. (2020).

### Installation 
Install the packages `graph` and `RBGL` from Bioconductor. 
Make sure that [rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed on your computer. 
Then, simply type in R to install `micd`

```R
devtools::install_github("bips-hb/micd")
```
(Note: The latest `micd` was created using R 4.2.1)

### References

Foraita R, Friemel J, Günther K, Behrens T, Bullerdiek J, Nimzyk R, Ahrens W, Didelez V (2020). Causal discovery of gene regulation with incomplete data. Journal of the Royal Statistical Society: Series A (Statistics in Society), 183(4), 1747-1775. URL https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12565.

Foraita R, Witte J, Börnhorst C, Gwozdz W, Pala V, Lissner L, Lauria F, Reisch L, Molnár D, De Henauw S, Moreno L, Veidebaum T, Tornaritis M, Pigeot I, Didelez V. A longitudinal causal graph analysis investigating modifiable risk factors and obesity in a European cohort of children and adolescents. 2022; medRxiv [2022.05.18.22275036](https://www.medrxiv.org/content/10.1101/2022.05.18.22275036v1).

Witte J, Foraita R, Didelez V (2022). Multiple imputation and test-wise deletion for causal discovery with incomplete cohort data. Statistics in Medicine, <https://doi.org/10.1002/sim.9535>.


### References for pcalg

Markus Kalisch, Martin Mächler, Diego Colombo, Marloes H. Maathuis, Peter Bühlmann (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47(11), 1-26. URL [www.jstatsoft.org/v47/i11/](https://www.jstatsoft.org/article/view/v047i11).

Alain Hauser, Peter Bühlmann (2012). Characterization and greedy learning of interventional Markov equivalence classes of directed acyclic graphs. Journal of Machine Learning Research, 13, 2409-2464. URL <https://www.jmlr.org/papers/v13/hauser12a.html>.

