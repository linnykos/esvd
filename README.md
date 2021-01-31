# Purpose

This repository contains all the data, functions, scripts to run simulations and analysis, and scripts to generate plots for the paper
"Exponential-family embedding with application to cell developmental trajectories 
  for single-cell RNA-seq data".

This code was developed and tested primarily on R 4.0.3 on a Macbook (macOS 11.1 Big Sur) equipped with an i5 processor.

# Installation

This package can be installed through `devtools` in R.

```{r}
library("devtools")
devtools::install_github("linnykos/esvd", subdir = "eSVD")
```
The package itself depends on several packages. These include `MASS`, `foreach`, `doMC`, `princurve`, `igraph`, `clplite`, `softImpute`, `RSpectra`, `plot3D`, `np`, `org.Mm.eg.db`, and `DBI`.
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
is publicly available on GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75330),
we provide a locally preprocessed dataset, which was created to be amendable for our analysis in R.
This dataset is a 21 MB `.RData` file, and is synced onto GitHub using the Git Large File Storage system (https://git-lfs.github.com/). Please
install this system before proceeding.

In the appendix, we investigate the Zeisel single-cell data collected by Zeisel et al. (2015). Similarly, 
while the original dataset 
is publicly available on GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361),
we provide a locally preprocessed dataset, which was created to be amendable for our analysis in R.
This dataset is a 13 MB `.RData` file, and is synced onto GitHub using the Git Large File Storage system 

# Reproducing the results

## Note

All the code below were run on a server with 15 cores. 

## Running the analysis

To reproduce the Marques (i.e., oligodendrocyte) analysis (Section 2, Section 7, and Appendix H of our paper, including Figures 1-3, 6-8, S.1, S.8-S.17), navigate to the `main` folder. From this location, run the R scripts command window. All the results and figures in these sections are reproduced by running `main.R`, which calls 16 different R scripts in succession, each producing `.RData` files that the next uses as input. The figures are produced in the last 8 script, `step8_figures_zz_data.R` through `step8_figures_zz_additional_analyses.R`.

Specifically: 

* `step8_figures_zz_data.R` produces Figures 3 and S.1.  

* `step8_figures_zz_training_testing.R` produces Figures 2, 7 and S.13. 

* `step8_figures_zz_2D_densities.R` produces Figures 1 and S.11. 

* `step8_figures_zz_2D_embedding.R` produces Figures S.8 and S.10.

* `step8_figures_zz_3D_embedding.R` produces Figures 6, 8, S.9 and S.12. 

* `step8_figures_zz_cascade.R` produces Figure S.14. 

* `step8_figures_zz_additional_analyses.R` produces Figures S.15-S.17.

To reproduce the Zeisel analysis (Appendix H.5), navigate to the `main_zeisel` folder. From this location, run the R scripts command window. All the results and figures in these sections are reproduced by running `main_zeisel.R`, which calls 5 different R scripts in succession, each producing `.RData` files that the next uses as input. Figure S.18 and Table S.1 are produced in the last script `step6_zeisel_metrics.R` and Figure S.19 is produced in the script, `step4_zeisel_analysis.R`.

## Running the simulations

To reproduce the simulations (Section 6 and Appendix D of our paper), navigate to the `simulation` folder. Below, we describe which scripts are associated with which files. Almost all these scripts depend on `factorization_generator.R` and `factorization_methods.R`.

* Figure 4: Run `illustration_example.R`.

* Figure 5: Run `factorization_suite_negbinom_esvd.R` and `factorization_suite_negbinom_rest.R`, followed by `factorization_suite_negbinom_postprocess.R`.

* Figure S.2: Run `consistency_simulation.R` followed by `consistency_simulation_plot.R`.

* Figure S.3: Run `factorization_suite_poisson_esvd.R` and `factorization_suite_poisson_rest.R`, followed by `factorization_suite_poisson_postprocess.R`.

* Figure S.4: Run `factorization_suite_curved_gaussian_esvd.R` and `factorization_suite_curved_gaussian_rest.R`, followed by `factorization_suite_curved_gaussian_postprocess.R`.

* Figure S.5: Run `factorization_suite_zinbwave_esvd.R` and `factorization_suite_zinbwave_rest.R`, followed by `factorization_suite_zinbwave_postprocess.R`.

* Figure S.6: Run `factorization_suite_tuning_zinbwave` followed by `factorization_suite_tuning_zinbwave_postprocess.R`.

* Figure S.7: Run `factorization_suite_pcmf_esvd.R` and `factorization_suite_pcmf_rest.R`, followed by `factorization_suite_pcmf_postprocess.R`.




