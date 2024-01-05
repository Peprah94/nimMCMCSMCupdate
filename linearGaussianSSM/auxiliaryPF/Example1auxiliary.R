
## load the nimble library and set seed
load("linearGaussianSSM/simulatedDataEx1.RData")


# Load packages
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)
library(parallel)
library(doParallel)

pfTypeRun = "auxiliary" #tyoe of particle filter to run

# load function to fit the models with MCMC
source("linearGaussianSSM/functionSimEstimation.R")


thisCluster <- makeCluster(5)

# Fit an auxiliary PF
auxiliaryEstimates <- parallel::parLapply(cl = thisCluster,
                                          X = 16:30,
                                          fun = runFunction,
                                          simData = simData,
                                          iNodePrev = c(49, 45, 20, 10, 5),
                                          nIterations = 30000,
                                          nBurnin = 20000,
                                          nChains = 2,
                                          nThin= 1,
                                          nyears = 50,
                                          numParticles = 1000,
                                          pfTypeRun = "auxiliary")


#save results
save(auxiliaryEstimates, file = "Example1/auxiliaryPF/estimates2.RData")
stopCluster(thisCluster)

