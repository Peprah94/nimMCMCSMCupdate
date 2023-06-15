# This is the general function that is used to fit the model with MCMC
# and simulate the required data

# Load the required packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)



# Function to fit replicated data

runFunction <- function(iter,
                        simData,
                        numParticles,
                        mcmcRun = TRUE,
                        pfTypeRun,
                        nIterations,
                        nBurnin,
                        nChains,
                        nThin,
                        nyears,
                        iNodePrev ){

  library('nimble')
  library(nimbleSMC)
  library(nimMCMCSMCupdates)

  # create empty lists to store the results
  example1ReducedModelTrue <- list()
  example1ReducedModelFalse <- list()
  example1UpdatedModelTrue <- list()
  example1UpdatedModelFalse <- list()
  baselineModel <- list()


  #load data
  simData <- simData[[iter]]

  # NIMBLE CODE
  stateSpaceCode <- nimbleCode({
    x[1] ~ dnorm(b/(1-a), sd = 1)
    #x[1] <- 0
    y[1] ~ dnorm(x[1], sd = 1)
    for(i in 2:t){
      x[i] ~ dnorm(a*x[i-1] + b, sd = 1)
      y[i] ~ dnorm(x[i] * c, 1)
    }
    a ~ dunif(0, 0.999)
    b ~ dnorm(0, 1)
    c ~ dnorm(1,1)
    #d ~ T(dnorm(0, 1), 0.001, 4)


    #estimate biases
    ahat <- a - aTrue
    bhat <- b - bTrue
    chat <- c - cTrue
  })

  # ## define data
  data <- list(
    y = simData$y
  )
  
  # define constants
  constants <- list(
    t = nyears,
    aTrue = 0.5,
    bTrue = 0.5,
    cTrue = 1,
    dTrue = 1
  )
  
  #define initial values
  inits <- list(
    a = 0.2,
    b = 0.8,
    # mu0= 0.8,
    c = 1
  )

  # ## build the model
  stateSpaceModel <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  for(iNodetag in seq_along(iNodePrev)){

    data <- list(
      y = simData$y[-c((iNodePrev[iNodetag]+1):50)]
    )
    constants <- list(
      t = iNodePrev[iNodetag],
      aTrue = 0.5,
      bTrue = 0.5,
      cTrue = 1,
      dTrue = 1
    )
    
    #Define model for reduced model
    newModelReduced <- nimbleModel(stateSpaceCode,
                                   data = data,
                                   constants = constants,
                                   inits = inits,
                                   check = FALSE)

  # Fit reduced model with MCMC
    example1ReducedModelTrue[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimWeights(model = newModelReduced,
                                                              latent = "x",
                                                              nParFiltRun = numParticles,
                                                              mcmc = TRUE,
                                                              block = FALSE,
                                                              pfType = pfTypeRun,
                                                              MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                       additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                       n.iter = nIterations,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = nBurnin,
                                                                                       n.thin = nThin)
    )

    # Fit reduced model with pMCMC
    example1ReducedModelFalse[[iNodetag]]  <- spartaNimWeights(model = newModelReduced,
                                                               latent = "x",
                                                               nParFiltRun = numParticles,
                                                               mcmc = FALSE,
                                                               pfType = pfTypeRun,
                                                               MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                        additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                        n.iter = nIterations,
                                                                                        n.chains = nChains,
                                                                                        n.burnin = nBurnin,
                                                                                        n.thin = nThin)
    )



    ################
    # Updated Model
    ################
    #message(paste("Running updated model for iNodePrev = ", iNodePrev[iNodetag], "and a = ", aVars[aVarstag]))

    # data for updated model
    data <- list(
      y = simData$y[c((iNodePrev[iNodetag]):50)]
    )
    
    # constants for updated model
    constants <- list(
      t = 50 - iNodePrev[iNodetag] + 1,
      aTrue = 0.5,
      bTrue = 0.5,
      cTrue = 1
    )

    # Define model for updated model
    newModelUpdated <- nimbleModel(stateSpaceCode,
                                   data = data,
                                   constants = constants,
                                   inits = inits,
                                   check = FALSE)

    # Fit model with results from fitted reduced MCMC model fit
    example1UpdatedModelTrue[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                              reducedModel = newModelReduced,
                                                              latent = "x", #latent variable
                                                              nParFiltRun = numParticles,
                                                              propCov = diag(3)* c(0.1, 1,1),
                                                              pfType = pfTypeRun,
                                                              MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                       additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                       n.iter = (nIterations - nBurnin)/nThin,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = 0,
                                                                                       n.thin = 1),  #saved loglikelihoods from reduced model
                                                              postReducedMCMC = example1ReducedModelTrue[[iNodetag]],# MCMC summary to use as initial values
                                                              pfControl = list(saveAll = TRUE,
                                                                               smoothing = TRUE,
                                                                               mcmc = TRUE,
                                                                               M = nyears - iNodePrev[iNodetag],
                                                                               iNodePrev = 1)
    )


    # Fit updated model with reduced model pMCMC
    example1UpdatedModelFalse[[iNodetag]] <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                              reducedModel = newModelReduced,
                                                              nParFiltRun = numParticles,
                                                              latent = "x", #latent variable
                                                              pfType = pfTypeRun,
                                                              propCov = diag(3)* c(0.1, 1,1),
                                                              MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                       additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                       n.iter = (nIterations - nBurnin)/nThin,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = 0,
                                                                                       n.thin = 1),  #saved loglikelihoods from reduced model
                                                              postReducedMCMC = example1ReducedModelFalse[[iNodetag]],# MCMC summary to use as initial values
                                                              pfControl = list(saveAll = TRUE,
                                                                               #lookahead = "mean",
                                                                               smoothing = TRUE,
                                                                               mcmc = FALSE,
                                                                               M = nyears - iNodePrev[iNodetag],
                                                                               iNodePrev = 1)
    )
  }


  ############
  # Baseline Model
  ##############
  # ## define data, constants, and initial values
  data <- list(
    y = simData$y
  )
  constants <- list(
    t = nyears,
    aTrue = 0.5,
    bTrue = 0.5,
    cTrue = 1
  )
  inits <- list(
    a = 0.2,
    b = 1,
    c = 1
  )
  #
  #
  # ## build the model
  stateSpaceModel <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)


  #Define a new model
  newModel <- stateSpaceModel$newModel(replicate = TRUE)

  # Function to run the baseline model
  
  # Fit baseline model with MCMC
  baselineModel[[1]] <- nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
                                                               latent = "x",
                                                               pfType = pfTypeRun,
                                                               nParFiltRun = numParticles,
                                                               MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                        additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                        n.iter = nIterations,
                                                                                        n.chains = nChains,
                                                                                        n.burnin = nBurnin,
                                                                                        n.thin = nThin))
  

  # Fit baseline model with SMC
  baselineModel[[2]]  <- nimMCMCSMCupdates::spartaNimWeights(model = newModel,
                                                         latent = "x",
                                                         nParFiltRun = numParticles,
                                                         mcmc = TRUE,
                                                         pfType = pfTypeRun,
                                                         MCMCconfiguration = list(target = c('a', 'b', 'c'),
                                                                                  additionalPars = c("x", "ahat", "bhat", "chat"),
                                                                                  n.iter = nIterations,
                                                                                  n.chains = nChains,
                                                                                  n.burnin = nBurnin,
                                                                                  n.thin = nThin)
  )

#return results
  ret <- c(baselineModel,
            example1ReducedModelTrue,
           example1ReducedModelFalse,
           example1UpdatedModelTrue,
           example1UpdatedModelFalse)

  return(ret)
}


#Function to simulate data
sim2 <- function(a, b, c, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, b/(1-a), sd = 1 )
  y[1] <- rnorm(1, x[1], 1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1] + b, sd = 1)
    y[k] <- rnorm(1, x[k]*c, 1)
  }
  return(list(x=x, y=y))
}
