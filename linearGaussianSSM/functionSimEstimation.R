# This is the general function that is used to fit the model with MCMC
# and simulate the required data

# Load the required packages
library('nimble')
library(nimbleSMC)
library(nimMCMCSMCupdates)



# Function to fit replicated data

runFunction <- function(iter, #iteration indez
                        simDataAll, #simulated data
                        numParticles, #number of particles
                        mcmcRun = TRUE, #whether to use MCMC (TRUE) or SMC (FALSE)
                        pfTypeRun, #type of particle filter algorithm
                        nIterations, #number of iterations
                        nBurnin, #number of burn-in samples
                        nChains, #number of chains
                        nThin, #thin in samples
                        nyears, #number of years (T)
                        iNodePrev #number of years used to fit the reduced models (t)
                        ){

  #load packages
  library('nimble')
  library(nimbleSMC)
  library(nimMCMCSMCupdates)

  # create empty lists to store the results
  example1ReducedModelTrue <- list()
  example1UpdatedModelTrueBF <- list()
  example1UpdatedModelTrueAux <- list()
  baselineModel <- list()


  #load data for the particular iteration index
  simData <- simDataAll[[iter]]

  # NIMBLE CODE
  stateSpaceCode <- nimbleCode({
    x[1] ~ dnorm(a, sd = 1)
    y[1] ~ dnorm(c*x[1], sd = 1)
    for(i in 2:t){
      x[i] ~ dnorm(a*x[i-1], sd = 1)
      y[i] ~ dnorm(x[i] * c, 1)
    }
    a ~ dunif(0, 0.999)
    c ~ dnorm(1,1)

    #estimate biases
    ahat <- a - aTrue
    chat <- c - cTrue
  })

  # ## define data
  data <- list(
    y = simData$y #observed data
  )

  # define constants
  constants <- list(
    t = nyears,
    aTrue = 0.5,
    cTrue = 1
  )

  #define initial values
  inits <- list(
    a = 0.2,
    c = 1
  )

  # ## build the model
  stateSpaceModel <- nimbleModel(stateSpaceCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  for(iNodetag in seq_along(iNodePrev)){

    ################
    # Reduced Model
    ################
    data <- list(
      y = simData$y[-c((iNodePrev[iNodetag]+1):50)] #select the subset of observed data
    )
    constants <- list(
      t = iNodePrev[iNodetag], #year used to fit the reduced model
      aTrue = 0.5,
      cTrue = 1
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
                                                              block = TRUE,
                                                              pfType = pfTypeRun,
                                                              MCMCconfiguration = list(target = c('a',  'c'),
                                                                                       additionalPars = c("x", "ahat", "chat"),
                                                                                       n.iter = nIterations,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = nBurnin,
                                                                                       n.thin = nThin)
    )




    ################
    # Updated Model
    ################

    # data for updated model
    data <- list(
      y = simData$y
    )

    # constants for updated model
    constants <- list(
      t = 50 ,
      aTrue = 0.5,
      cTrue = 1
    )

    inits <- list(
      a = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["a",1 ],
      c = example1ReducedModelTrue[[iNodetag]]$summary$all.chains["c",1 ]
    )

    # Define model for updated model
    newModelUpdated <- nimbleModel(stateSpaceCode,
                                   data = data,
                                   constants = constants,
                                   inits = inits,
                                   check = FALSE)

    # Fit model with results from fitted reduced MCMC model fit
    example1UpdatedModelTrueBF[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                              reducedModel = newModelReduced,
                                                              latent = "x", #latent variable
                                                              nParFiltRun = numParticles,
                                                              propCov = diag(2)* c(0.1,1),
                                                              pfType = "bootstrap",
                                                              extraVars = NULL,
                                                              MCMCconfiguration = list(target = c('a', 'c'),
                                                                                       additionalPars = c("x", "ahat", "chat"),
                                                                                       n.iter = (nIterations - nBurnin)/nThin,
                                                                                       n.chains = nChains,
                                                                                       n.burnin = 0,
                                                                                       n.thin = 1),  #saved loglikelihoods from reduced model
                                                              postReducedMCMC = example1ReducedModelTrue[[iNodetag]],# MCMC summary to use as initial values
                                                              pfControl = list(saveAll = TRUE,
                                                                               smoothing = TRUE,
                                                                               mcmc = TRUE,
                                                                               M = nyears - iNodePrev[iNodetag],
                                                                               iNodePrev = iNodePrev[iNodetag])
    )


    example1UpdatedModelTrueAux[[iNodetag]]  <- nimMCMCSMCupdates::spartaNimUpdates(model = newModelUpdated, #nimble model
                                                                                 reducedModel = newModelReduced,
                                                                                 latent = "x", #latent variable
                                                                                 nParFiltRun = numParticles,
                                                                                 propCov = diag(2)* c(0.1,1),
                                                                                 pfType = "auxiliary",
                                                                                 extraVars = NULL,
                                                                                 MCMCconfiguration = list(target = c('a', 'c'),
                                                                                                          additionalPars = c("x", "ahat", "chat"),
                                                                                                          n.iter = (nIterations - nBurnin)/nThin,
                                                                                                          n.chains = nChains,
                                                                                                          n.burnin = 0,
                                                                                                          n.thin = 1),  #saved loglikelihoods from reduced model
                                                                                 postReducedMCMC = example1ReducedModelTrue[[iNodetag]],# MCMC summary to use as initial values
                                                                                 pfControl = list(saveAll = TRUE,
                                                                                                  smoothing = TRUE,
                                                                                                  mcmc = TRUE,
                                                                                                  lookahead = "mean",
                                                                                                  M = nyears - iNodePrev[iNodetag],
                                                                                                  iNodePrev = iNodePrev[iNodetag])
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
                                                               pfType = "bootstrap",
                                                               nParFiltRun = numParticles,
                                                               pfControl = list(saveAll = TRUE,
                                                                                lookahead = "mean",
                                                                                smoothing = FALSE),
                                                               MCMCconfiguration = list(target = c('a', 'c'),
                                                                                        additionalPars = c("x", "ahat", "chat"),
                                                                                        n.iter = nIterations,
                                                                                        n.chains = nChains,
                                                                                        n.burnin = nBurnin,
                                                                                        n.thin = nThin))

  baselineModel[[2]] <- nimMCMCSMCupdates::baselineSpartaEstimation(model = newModel,
                                                                    latent = "x",
                                                                    pfType = "auxiliary",
                                                                    nParFiltRun = numParticles,
                                                                    pfControl = list(saveAll = TRUE,
                                                                                     lookahead = "mean",
                                                                                     smoothing = FALSE),
                                                                    MCMCconfiguration = list(target = c('a', 'c'),
                                                                                             additionalPars = c("x", "ahat", "chat"),
                                                                                             n.iter = nIterations,
                                                                                             n.chains = nChains,
                                                                                             n.burnin = nBurnin,
                                                                                             n.thin = nThin))


  # Fit baseline model with SMC
  baselineModel[[3]]  <- nimMCMCSMCupdates::spartaNimWeights(model = stateSpaceModel,
                                                         latent = "x",
                                                         nParFiltRun = numParticles,
                                                         mcmc = TRUE,
                                                         block = TRUE,
                                                         pfType = pfTypeRun,
                                                         MCMCconfiguration = list(target = c('a', 'c'),
                                                                                  additionalPars = c("x", "ahat", "chat"),
                                                                                  n.iter = nIterations,
                                                                                  n.chains = nChains,
                                                                                  n.burnin = nBurnin,
                                                                                  n.thin = nThin)
  )

#return results
  ret <- c(baselineModel,
            example1ReducedModelTrue,
           example1UpdatedModelTrueBF,
           example1UpdatedModelTrueAux
           )

  return(ret)
}


#Function to simulate data
sim2 <- function(a, c, t, seed){
  set.seed(seed)
  x <- y <- numeric(t)
  x[1] <- rnorm(1, a, sd = 1 )
  y[1] <- rnorm(1, x[1], sd=  1)

  for(k in 2:t){
    x[k] <- rnorm(1, a*x[k -1], sd = 1)
    y[k] <- rnorm(1, x[k]*c, sd= 1)
  }
  return(list(x=x, y=y))
}
