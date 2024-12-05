# TensorComBat
Tensor-ComBat is an image harmonization method which adjusts imaging outcomes while controlling for site effects and other covariates using Bayesian tensor response regression (BTRR). The model was motivated from a project which aimed to harmonize voxel-level cortical thickness measures to bolster downstream biological analyses and reproducibility using ADNI-1 neuroimaging data.

Relevant files with Tensor-ComBat functions are:
  **TC_MCMC.R** (MCMC Call with Rcpp),
  **TC_MCMC_c.cpp** (C++ functions to MCMC sampling),
  **TC_Harmonize.R** (Harmonization of 2D or 3D images, given group ids, and options for covariate control and longitudinal modeling) 
