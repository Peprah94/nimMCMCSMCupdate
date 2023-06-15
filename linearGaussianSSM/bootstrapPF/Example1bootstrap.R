## load the nimble library and set seed
load("Example1/simulatedDataEx1.RData")


library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)

pfTypeRun = "bootstrap"

# load data

source("Example1/functionSimEstimation.R")


library(parallel)
library(doParallel)
thisCluster <- makeCluster(5)
#doParallel::registerDoParallel(cl)
#setDefaultCluster(cl)
#clusterExport(cl, c("runFunction"))

# bootstrapEstimates <- pbapply::pblapply(simData[1:10], function(x){
#   parallelRun(simData = x,
#               iNodePrev = c(49, 45,  5),
#               nIterations = 10,
#               nBurnin = 2,
#               nChains = 2,
#               nThin= 1,
#               nyears = 50,
#               numParticles = 5,
#               pfTypeRun = "bootstrap")
# }, cl = cl)
bootstrapEstimates <- parallel::parLapply(cl = thisCluster,
                                          X = 1:15,
                                          fun = runFunction,
                                          simData = simData,
                                          iNodePrev = c(49, 45, 20, 10, 5),
                                          nIterations = 30000,
                                          nBurnin = 20000,
                                          nChains = 3,
                                          nThin= 1,
                                          nyears = 50,
                                          numParticles = 1000,
                                          pfTypeRun = "bootstrap")


#   foreach(iter = seq_along(simData[7:15]), .packages = c("nimble", "nimbleSMC", "nimMCMCSMCupdates")) %do% {
#   tryCatch({  runFunction(simData = simData[[iter]],
#                           iNodePrev = c(49, 45, 20, 10, 5),
#                           nIterations = 30000,
#                           nBurnin = 20000,
#                           nChains = 3,
#                           nThin= 1,
#                           nyears = 50,
#                           numParticles = 200,
#                           pfTypeRun = "bootstrap") }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }

#save results
save(bootstrapEstimates, file = "Example1/bootstrapPF/estimatesBFNew1.RData")
stopCluster(thisCluster)

