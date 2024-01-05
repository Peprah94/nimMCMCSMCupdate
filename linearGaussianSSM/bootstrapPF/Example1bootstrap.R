## load the nimble library and set seed
load("linearGaussianSSM/simulatedDataEx1.RData")

# Load packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)
library(parallel)
library(doParallel)

pfTypeRun = "bootstrap" #type of particle fiter algorithm

# load data
source("linearGaussianSSM/functionSimEstimation.R")

# making cluster to run the functions parallely
thisCluster <- makeCluster(5)

# Fit the bootstrap PF
bootstrapEstimates <- parallel::parLapply(cl = thisCluster,
                                          X = 1:15,
                                          fun = runFunction,
                                          simData = simData,
                                          iNodePrev = c(49, 45, 20, 10, 5),
                                          nIterations = 30000,
                                          nBurnin = 20000,
                                          nChains = 2,
                                          nThin= 1,
                                          nyears = 50,
                                          numParticles = 1000,
                                          pfTypeRun = "bootstrap")



#save results
save(bootstrapEstimates, file = "Example1/bootstrapPF/estimates1.RData")
stopCluster(thisCluster)

