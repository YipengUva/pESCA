---
title: "Demos and algorithms for the fitting and selection of the P-ESCA model"
author: "Yipeng Song, Biosystems Data Analysis Group, University of Amsterdam"
date: "February 16, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Make sure the current folder is in your Matlab path and run the help function from the Matlab console (e.g. help ESCA_group_concave) for more information on the input/output of the algorithms.

## Demos

The demos are used to show how to simulate multiple data sets with underlying global, local common and distinct structures according to the ESCA model, how to construct a P-ESCA model, and how to do the model selection. The explanation and the results of these demos can be found in the corresponding html files in the examples fold.

-**P-ESCA model on three quantitative data sets:**
   
- demo_PESCA_GGG.m: this demo shows how to do data simulation, model selection and fit the final P-ESCA model for multiple quantitative data sets.
   
- demo_PESCA_GGG_fullInfo.m: this demo is same as above, except that the simulated parameters are used to evaluate the model selection process. 

-**P-ESCA model on three binary data sets:**

- demo_PESCA_BBB.m: this demo shows how to do data simulation, model selection and fit the final P-ESCA model on multiple binary data sets.

- demo_PESCA_BBB_fullInfo.m: this demo is same as above, except that the simulated parameters are used to evaluate the model selection process.

-**P-ESCA model on mixed binary and quantitative data sets:**

- demo_PESCA_GBB.m: this demo shows how to do data simulation, model selection and fit the final P-ESCA model for mixed quantitative-binary-binary data sets.

- demo_PESCA_GGB.m: this demo shows how to do data simulation, model selection and fit the final P-ESCA model for mixed quantitative-quantitative-binary data sets.

## Examples:
- The html documents of the above demos published using Matlab are in the examples fold.

## Algorithms:
The algorithms are used to do data simulation, to fit a P-ESCA model, and to do model selection are in the algorithms fold.

-**dataSimulation:**

- This folder contains the functions used to do the data simulations.

-**functions:**

- This folder contains some functions used in the algorithm.

-**Dispersion parameter estimation:**

- alpha_estimation.m: the developed \alpha estimation procedure.

- SVD_CV_modelSelection.m: missing value based cross validation procedure for the rank selection of a SVD (PCA) model. Mainly used in alpha_estimation.m.

- SVD_missingValues.m: a SVD (PCA) algorithm, which is capable to tackle missing value problem.

-**Algorithms to fit a P-ESCA model:**

- ESCA_group_concave.m: the main algorithm developed to fit the P-ESCA model.
	
- ESCA_group_concave_composite.m: the algorithm used to fit a P-ESCA model when composite concave penalty is used in the P-ESCA model. This function is a part of ESCA_group_concave.m. 
		
- ESCA_group_concave_L1.m: the algorithm used to fit a P-ESCA model when concave L1-norm is used in the P-ESCA model. This function is a part of ESCA_group_concave.m.
		
- ESCA_group_concave_L2.m: the algorithm used to fit a P-ESCA model when concave L2-norm penalty is used in the P-ESCA model. This function is a part of ESCA_group_concave.m. This algorithm is used in the P-ESCA paper.

-**Algorithms to select the P-ESCA model:**

- ESCA_modelSelection_KCV.m: missing value based cross validation to select the P-ESCA model.

- ESCA_modelSelection_KCV_fullInfo.m: same as the above procedure. And the simulated parameters are used to evaluate the model selection procedure.

- ESCA_modelSelection_KCV_twoSteps_bg.m: the model selection procedure for mixed quantitative and binary data sets. 

- ESCA_modelSelection_KCV_twoSteps_Gaussian.m and ESCA_modelSelection_KCV_twoSteps_nonGaussian.m are used in ESCA_modelSelection_KCV_twoSteps_bg.m.

-**Others:**

- scaling_centering.m: a column centering and scaling procedure for the processing of the real data sets.


