# This function is used to simulate the data for linear Gaussian SSM
source("linearGaussianSSM/functionSimEstimation.R")

#Simulate dataset
numSimulations <- 100

simData <- lapply(1:numSimulations, function(x){
  sim2(a = 0.5,
       #b = 0.5,
       c = 1,
      t = 50,
      x)
})

# save the data
save(simData, file = "linearGaussianSSM/simulatedDataEx1.RData")
