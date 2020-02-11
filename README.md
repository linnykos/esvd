# Purpose

This repository contains all the data, functions, scripts to run simulations and analysis, and scripts to generate plots for the paper
"Exponential-family embedding with application to cell developmental trajectories 
  for single-cell RNA-seq data".

# Installation

This package can be installed through `devtools` in R.

```{r}
library("devtools")
devtools::install_github("linnylin92/esvd", subdir = "eSVD")
```
The package itself depends on several packages. These include `MASS', `foreach', `doMC', `princurve', `igraph', `clplite', `softImpute', `NMF', and `plot3D'
Warning: On Windows, to install the `doMC` package, use the following code in R.
```{r}
install.packages("doMC", repos="http://R-Forge.R-project.org")
```
See: http://stackoverflow.com/questions/16453625/package-domc-not-available-for-r-version-3-0-0-warning-in-install-packages

The above installation is only for the R package. To reproduce the entire simulation and analysis, you will need to pull/fork this entire repository.
You will need to install the Git Large File Storage system to do this (see below).

# Data 

The dataset used in this article is also included in the repository.
This is the Marques single-cell dataset collected by Marques et al. (2016). While the original dataset 
is publicly available on GEO (\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75330}),
we provide a locally preprocessed dataset, which was created to be amendable for our analysis in R.
This dataset is a 21 MB `.RData` file, and is synced onto GitHub using the Git Large File Storage system (https://git-lfs.github.com/). Please
install this system before proceeding.

# Reproducing the results

## Note

All the code below were run on a server with 15 cores. 

## Running the simulations

To reproduce the simulations (Section 5 and Appendix D of our paper), navigate to the `simulation` folder. 

* To reproduce Figure 2, run `illustration_example.R`.

* To reproduce Figure 3, run `factorization_suite.R` and `factorization_suite_postprocess.R`.

* To reproduce Figure 8, run `wasserstein_simulation.R` and `wasserstein_simulation_plot.R`.

* To reproduce Figure 9 through 13, run `factorization_suite_others.R` and `factorization_suite_others_postprocess.R`.


## Running the analysis

To reproduce the analysis (Section 1, Section 6, and Appendix G of our paper), navigate to the `main` folder. From this location, run the following lines in the command window. All the results and figures in these sections are reproduced by running `main.R`, which calls 7 different R scripts in succession. The figures are produced in the last two scripts, `step6_figures.R` and `step7_figures_embedding.R`.