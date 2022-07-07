# Fully reproducible R code and data from Masselot et al. Constrained groupwise additive index models (2022) Biostatistics

Data and R code to reproduce the simulation study and real-world analyses of:

------

Pierre Masselot, Fateh Chebana, Céline Campagna, Éric Lavigne, Taha B M J Ouarda, Pierre Gosselin, **Constrained groupwise additive index models**, *Biostatistics*, 2022; https://doi.org/10.1093/biostatistics/kxac023

------

This work introduces Constrained groupwise additive index models (CGAIM), a model to construct indices from groups of variables, with indices non-linearly related to an outcome of interest. The paper describes the model, provides a fitting algorithm based on quadratic programming, and proposes inference and index selection procedures.

This repository contains the R code used to perform the four simulations studies and two applications includes in the paper. R implementation of the CGAIM is available in the R package `cgaim`.

## Details

The repository includes the following scripts.

0. These scripts contain functions implementing the methods to which the CGAIM is compared.
  1. *0.1_gMAVE.R*: Implements the gMAVE estimator that estimates groupwise indices based on local-linear regression.
  2. *0.2_FACTS.R*: Implements the FACTS model, a similar model to gMAVE with the addition of some constraints on the estimates.
  3. *0.3_PPR.R*: A simple wrapper for the `ppr` function in R.
  
1. The following scripts implement the four simulation experiments performed in the article.
  1. *1.1_Simulation_IndexEstimation.R*: Comparison of the RMSE of each algorithm in the estimation of index coefficients `alpha`.
  2. *1.2_Simulation_IndexSelection.R*: Evaluation of index selection by the proposed GCV criterion.
  3. *1.3_Simulation_Coverage.R*: Estimation of confidence intervals coverage from two inference procedures proposed in the paper.
  4. *1.4_Simulation_Exposome.R*: Evaluation of the CGAIM on a more complex setting mimicking exposome studies.
  
2. The following scripts implement the applications performed in the paper.
  1. *2.1_CaseStudy_Heat.R*: Fully reproducible code for the heat index construction in Montreal.
  2. *2.2_CaseStudy_Pollution.R*: Reproducible code for the air pollution index construction application found in Supplementary Materials.
  
3. Additional code creating illustrative plots found in Supplementary Materials.

In addition to these scripts, the repository contains a folder `Data` that includes data used in the simulation studies and applications. The folder `temp` includes files used to track the simulations performed in parallel and is not directly of interest.

## Notes

Some results might be slightly different than in the paper due to deprecated functions after the article was accepted. In particular, the function previously used to extract derivatives from `gam` objects has changed since the paper was accepted slightly impacting the results of unconstrained models.
