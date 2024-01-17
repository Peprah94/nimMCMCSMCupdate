# A sequential Monte Carlo algorithm for data assimilation problems in ecology



This repository hosts the code, supplementary materials and figures used to write the manuscript. To be able to run the scripts in this repository, the following R-packages must be installed in R: `nimbleSMC`, `nimble`, `sparta` and `nimMCMCSMCupdates`. The later can be installed from this [github repository](https://github.com/Peprah94/nimMCMCSMCupdates).

## Fitting models

The folders we describe below should contain the Rscipt and results from the analysis. However, some of the results files are relatively large, and could not be updated on GitHib. We have uploaded them to a [dropbox account](https://www.dropbox.com/scl/fo/27ha5bb7zbb5djabxjdsp/h?rlkey=jhdut9lrdxapmg6r73ksm0dc0&dl). 


### Linear Gaussian state space model (/linearGaussianSSM)

This folder contains R-script for fitting the linear Gaussian state space model described in simulation study one of the paper. The folder contains the files and folders: 

* `functionSimEstimation.R` - contains functions to simulate the data

* `simData.R` - simulate the data

* `Example1.R` - fit the models

* `simulatedDataEx1.RData` - simulated dataset

* `auxiliaryPF` - contains script to run the auxiliary particle filter

* `bootstrapPF` - contains script to run the auxiliary particle filter

### Dynamic Occupancy model (/dynamicOccupancyModel)

This folder contains R-script for fitting the dynamic occupancy model described in simulation study two of the paper. The folder contains the files and folders: 

* `auxiliaryEst.R` - Simulate data and fit the dynamic occupancy model to the data.

* `Archive.zip` - Results from the fitted model.

### Demographic SSM (/demographicSSM)

This folder contains R-script for fitting the demographic state-space model described in case study one in the paper. The folder contains the files and folders: 

* `populationDemographics.R` - Fit the demographic SSM to the data. The data is accessed from the R-package `AHMbook`.

* `Archive.zip` - Results from the fitted model.

### Occupancy model fitted using `sparta` (/spartaOccupancyModel)

This folder contains R-script for fitting the occupancy model described in case study two in the paper. The folder contains the files and folders: 

* `example3.R` - Fit the demographic SSM to the data. 

## Plotting results

The `plots.R` script is used to summarise the results obtained from the analysis. The figures and tables in the main paper are in the `Figures` folder.

## Supplementary Materials

This folder hosts three pdf files that provides additional information to the manuscript.


