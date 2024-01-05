# A sequential Monte Carlo algorithm for data assimilation problems in ecology

This repository hosts the code, supplementary materials and figures used to write the manuscript. To be able to run the scripts in this repository, the following R-packages must be installed in R: `nimbleSMC`, `nimble`, `sparta` and `nimMCMCSMCupdates`. The later can be installed from this [github repository](https://github.com/Peprah94/nimMCMCSMCupdates).

## Fitting models

The folders we describe below should contain the Rscipt and results from the analysis. However, some of the results files are relatively large, and could not be updated on GitHib. We have uploaded them to a [dropbox account](https://www.dropbox.com/scl/fo/27ha5bb7zbb5djabxjdsp/h?rlkey=jhdut9lrdxapmg6r73ksm0dc0&dl=0). 


### Linear Gaussian state space model (/linearGaussianSSM)

This folder contains R-script for fitting the linear Gaussian state space model described in simulation study one of the paper. The folder contains the files and folders: 

* `functionSimEstimation.R` - contains functions to simulate the data

* `simData.R` - simulate the data

* `Example1.R` - fit the models

* `simulatedDataEx1.RData` - simulated dataset

* `auxiliaryPF` - contains script to run the auxiliary particle filter

* `bootstrapPF` - contains script to run the auxiliary particle filter

 
