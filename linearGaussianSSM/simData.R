# This function is used to simulate the data for example 1
source("Example1/functionSimEstimation.R")

#Simulate dataset
numSimulations <- 30

simData <- lapply(1:numSimulations, function(x){
  sim2(a = 0.5,
       b = 0.5,
       c = 1,
      t = 50,
      x)
})

save(simData, file = "Example1/simulatedDataEx1.RData")
