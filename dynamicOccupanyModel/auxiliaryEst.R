# Dynamic occupancy model
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)

# set the configurations for fitting the model
nyears = 30
nsites = 300
nvisits = 3
iNodePrevAll = c(29, 25) # the time steps used to fit the reduced models

# Set the configurations for MCMC
nIterations = 50000
nBurnin = 20000
nChains = 2
nThin = 15
numParticles = 10 # change this to 25, 50 and 100 for the particular number of particles you are interested in.

##################################################
# Function to simulate dynamic occupancy models  #
##################################################

dynOccupancyModels <- function(nyears, #number of years
                               nsites, #number of sites
                               nvisits, #number of visits
                               #intercept and covariate effect for initial occupancy probability
                               fixedPars = list(alphaPhi = - 2,
                                                betaPhi = 1.5),
                               # standard deviation hyperparameters
                               hyperParsSig = list(alphaP = 2,
                                                   betaP = 3,
                                                   alphaPsi = 2,
                                                   betaPsi = 3)){

  set.seed(1994)

  # Simulating covariates
  windSpeed <- array(runif(nyears*nsites*nvisits, -1, 1),
                     dim = c(nsites, nvisits, nyears))
  elevation <- runif(nsites, -1,1)
  springPrecipitation <- matrix(runif(nyears*nsites, -1, 1),
                                nrow = nsites,
                                ncol = nyears)
  sizeOfBeak <- matrix(runif(nyears*nsites, -1,1),
                       nrow = nsites,
                       ncol = nyears)

  # Simulating parameters
  alphaP <- rnorm(1, mean = 4, sd = hyperParsSig$alphaP)
  betaP <- rnorm(1, mean = -2, sd = hyperParsSig$betaP)
  alphaPsi <- rnorm(1, mean = 3, sd = hyperParsSig$alphaPsi)
  betaPsi <- rnorm(1, mean = -2.5, sd = hyperParsSig$betaPsi)

  # Detection Probability
  detectProb <- array(NA, dim = c(nsites, nvisits, nyears))
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        detectProb[site.tag, visit.tag, year.tag] <- plogis(alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag])
      }
    }
  }

  # Initial occupancy probability psi
  initOccuProb <- plogis(fixedPars$alphaPhi + fixedPars$betaPhi*elevation)

  # Persistence and colonisation probability
  persistenceProb  <- matrix(NA,  nrow = nsites, ncol = nyears)
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      persistenceProb[site.tag,  year.tag] <- plogis(alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag])
    }
  }

  colonisationProb <- 0.05

  # Latent state and observations
  y <- array(NA, dim = c(nsites, nvisits, nyears))
  z <- matrix(NA,nrow = nsites,ncol = nyears)

  # Initial Presence/Absence
  for(site.tag in 1:nsites){
    z[ site.tag, 1] <- rbinom(1, 1, initOccuProb[site.tag])
  }

  # True presence/absence
  for(year.tag in 2:nyears){
    for(site.tag in 1:nsites){
      z[ site.tag, year.tag] <- rbinom(1, 1, z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb)#[site.tag,  year.tag])
    }
  }

  #observations
  for(year.tag in 1:nyears){
    for(site.tag in 1:nsites){
      for(visit.tag in 1:nvisits){
        y[site.tag, visit.tag, year.tag] <- rbinom(1, 1,  z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag])
      }
    }
  }

  # Proportion of occupied sites
  psi.fs <- colMeans(z)

  # Return list
  retList <- list()
  retList$y = y
  retList$z = z
  retList$covariates = list(windSpeed = windSpeed,
                            elevation = elevation,
                            springPrecipitation = springPrecipitation,
                            sizeOfBeak = sizeOfBeak)
  retList$trueSigma = hyperParsSig
  retList$truePars = fixedPars
  retList$covariateEffects = list(alphaP = alphaP,
                                  betaP = betaP,
                                  alphaPsi = alphaPsi,
                                  betaPsi = betaPsi,
                                  alphaGamma = alphaGamma,
                                  betaGamma = betaGamma)

  retList$occSites = psi.fs

  return(retList)
}

# Simulating and saving data
simData <- dynOccupancyModels(nyears = 30,
                               nsites = 300,
                               nvisits = 3)

 save(simData, file = "simDataDynamicOccupancy.RData")

 # Uncomment this if you change the number of particles
#load("simDataDynamicOccupancy.RData")

# NIMBLE code
# Try to do this first before updating, so I can use it for the bootstrap updating
for(i in seq_along(iNodePrevAll)){
  #number of years to fit the reduced model
   iNodePrev <- iNodePrevAll[i]
  dynOccupancyCode <- nimbleCode({

    # Prior distributions of hyperparameters
    alphaPSig ~ T(dnorm(0, 0.1), 0.001, )
    betaPSig ~ T(dnorm(0, 0.1), 0.001, )
    alphaPsiSig ~ T(dnorm(0, 0.1), 0.001, )
    betaPsiSig ~ T(dnorm(0, 0.1), 0.001, )

    # Prior distribution for intercept and covariate effect
    alphaPhi ~ dnorm(0, sd = 10)
    betaPhi ~ dnorm(0, sd = 10)

    # Prior distributions
    alphaP ~ dnorm(mean = 4, sd = alphaPSig)
    betaP ~ dnorm(mean = -2, sd =  betaPSig)
    alphaPsi ~ dnorm(mean = 3, sd =  alphaPsiSig)
    betaPsi ~ dnorm(mean = -2.5, sd =  betaPsiSig)

    # Detection Probability
    for(year.tag in 1:nyears){
      for(site.tag in 1:nsites){
        for(visit.tag in 1:nvisits){
          logit(detectProb[site.tag, visit.tag, year.tag]) <- alphaP + betaP* windSpeed[site.tag, visit.tag, year.tag]
        }
      }
    }

    # Initial occupancy probability psi
    for(site.tag in 1:nsites){
      logit(initOccuProb[site.tag]) <- alphaPhi + betaPhi*elevation[site.tag]
    }

    # Persistence and colonisation probability
    for(year.tag in 1:nyears){
      for(site.tag in 1:nsites){
        logit(persistenceProb[site.tag,  year.tag]) <- alphaPsi + betaPsi* springPrecipitation[site.tag, year.tag]
        # logit(colonisationProb[site.tag,  year.tag]) <- alphaGamma[year.tag] + betaGamma[year.tag]* sizeOfBeak[site.tag, year.tag]
      }
    }

    colonisationProb ~ dunif(0.01, 1)


    # Initial Presence/Absence
    for(site.tag in 1:nsites){
      z[site.tag, 1] ~ dbin(prob = initOccuProb[site.tag], 1)
    }

    # True presence/absence
    for(year.tag in 2:nyears){
      for(site.tag in 1:nsites){
        z[site.tag, year.tag] ~ dbin(prob = (z[site.tag, (year.tag -1)] * persistenceProb[site.tag,  year.tag] + (1 - z[site.tag, (year.tag -1)])*colonisationProb), 1)
      }
    }

    #observations
    for(year.tag in 1:nyears){
      for(site.tag in 1:nsites){
        for(visit.tag in 1:nvisits){
          y[site.tag, visit.tag, year.tag] ~ dbin(prob = z[site.tag, year.tag] * detectProb[site.tag, visit.tag, year.tag], 1)
        }
      }
    }

    # Derived quantities
    for(year.tag in 1:nyears){
      psi.fs[year.tag] <- sum(z[1:nsites, year.tag])/nsites
    }

  })

## define data, constants, and initial values
  z <- apply(simData$y, c(1, 3), function(x){
    ifelse(sum(x) > 0, 1, NA)
  })

#####################
#   Reduced Model
########################

  data <- list(
    y = simData$y[,,-c((iNodePrev +1):nyears)],
    windSpeed = simData$covariates$windSpeed[,,-c((iNodePrev +1):nyears)],
    elevation = simData$covariates$elevation,
    springPrecipitation = simData$covariates$springPrecipitation[,-c((iNodePrev +1):nyears)],
    sizeOfBeak = simData$covariates$sizeOfBeak,
    z = z[,-c((iNodePrev +1):nyears)]
  )

  constants <- list(
    nyears = iNodePrev,
    nsites = nsites,
    nvisits = nvisits
  )
  inits <- list(
    alphaPSig = 1,
    betaPSig  = 1,
    alphaPsiSig = 1,
    betaPsiSig  = 1,
    alphaGammaSig = 1,
    betaGammaSig = 1,
    alphaPhi  = 0,
    betaPhi  = 0,
    alphaP =  rnorm(1, mean = 0, sd = 1),
    betaP =  rnorm(1, mean = 0, sd = 1),
    alphaPsi=  rnorm(1, mean = 0, sd = 1),
    betaPsi= rnorm(1, mean = 0, sd = 1),
    colonisationProb = 0.01,
    alphaPhi=  rnorm(1),
    betaPhi= rnorm(1)
  )

#   # NIMBLE model for reduced model
  newModelReduced <- nimbleModel(dynOccupancyCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

  # Fitting the model
  example2ReducedModelTrue <- spartaNimWeights(model = newModelReduced,
                                               latent = "z",
                                               nParFiltRun = numParticles,
                                               mcmc = TRUE,
                                               pfType = pfTypeRun,
                                               pfControl = list(saveAll = TRUE,
                                                                smoothing = TRUE,
                                                                timeIndex = 2),
                                               MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                   'alphaPsiSig', 'betaPsiSig',
                                                                                   'alphaPhi', 'betaPhi',
                                                                                   'alphaP', 'betaP',
                                                                                   'alphaPsi', 'betaPsi',
                                                                                   "colonisationProb"
                                               ),
                                               additionalPars = c( "psi.fs"
                                               ),
                                               n.iter = nIterations,
                                               n.chains = nChains,
                                               n.burnin = nBurnin,
                                               n.thin = nThin)
  )

  #save results
  save(example2ReducedModelTrue, file= paste0("example4ReducedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))

  ################
  # Updated Model
  ################
 load(paste0("example4ReducedBootstrapTrueM10Ind",iNodePrev,".RData"))

  data <- list(
    y = simData$y,
    windSpeed = simData$covariates$windSpeed,
    elevation = simData$covariates$elevation,
    springPrecipitation = simData$covariates$springPrecipitation,
    sizeOfBeak = simData$covariates$sizeOfBeak,
    z = z
  )

  constants <- list(
    nyears = nyears,
    nsites = nsites,
    nvisits = nvisits
  )


  inits <- list(
    alphaPSig = 1,
    betaPSig  = 1,
    alphaPsiSig = 1,
    betaPsiSig  = 1,
    alphaGammaSig = 1,
    betaGammaSig = 1,
    alphaPhi  = 0,
    betaPhi  = 0,
    alphaP =  rnorm(1, mean = 0, sd = 1),
    betaP =  rnorm(1, mean = 0, sd = 1),
    alphaPsi=  rnorm(1, mean = 0, sd = 1),
    betaPsi= rnorm(1, mean = 0, sd = 1),
    colonisationProb = 0.01
  )

  # NIMBLE model for updated model
  newModelUpdated <- nimbleModel(dynOccupancyCode,
                                 data = data,
                                 constants = constants,
                                 inits = inits,
                                 check = FALSE)

#   # Fitting the model
  example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                               reducedModel = newModelReduced,
                                               nParFiltRun = numParticles,
                                               pfType = "bootstrap",
                                               propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                               mcmcScale = 1,
                                               extraVars = NULL,
                                               MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                   'alphaPsiSig', 'betaPsiSig',
                                                                                   'alphaPhi', 'betaPhi',
                                                                                   'alphaP', 'betaP',
                                                                                   'alphaPsi', 'betaPsi',
                                                                                   "colonisationProb"),
                                                                        additionalPars = c("z", "psi.fs"),
                                                                        n.iter = (nIterations - nBurnin)/nThin,
                                                                        n.chains = nChains,
                                                                        n.burnin = 10,
                                                                        n.thin = 1
                                               ),  #saved loglikelihoods from reduced model
                                               postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                               pfControl = list(saveAll = TRUE,
                                                                timeIndex = 2,
                                                                smoothing = TRUE,
                                                                mcmc = TRUE,
                                                                M = nyears - iNodePrev,
                                                                iNodePrev = iNodePrev),
                                               nCores = 1
  )
#
   save(example2UpdatedModelTrue, file= paste0("example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))

#
  # Fitting the model
example2UpdatedModelTrue <- spartaNimUpdates(model = newModelUpdated, #nimble model
                                             reducedModel = newModelReduced,
                                             latent = "z", #latent variable
                                             nParFiltRun = numParticles,
                                             pfType = "auxiliary",
                                             propCov = c(1,1,1,1,1,1,0.01)*diag(7),
                                             mcmcScale = 1,
                                             extraVars = NULL,
                                             MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                                 'alphaPsiSig', 'betaPsiSig',
                                                                                 'alphaPhi', 'betaPhi',
                                                                                 'alphaP', 'betaP',
                                                                                 'alphaPsi', 'betaPsi',
                                                                                 "colonisationProb"),
                                                                      additionalPars = c("z", "psi.fs"),
                                                                      n.iter = (nIterations - nBurnin)/nThin,
                                                                      n.chains = nChains,
                                                                      n.burnin = 10,
                                                                      n.thin = 1),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              mcmc = TRUE,
                                                              M = nyears - iNodePrev,
                                                              iNodePrev = iNodePrev),
                                             nCores = 1
)

save(example2UpdatedModelTrue, file= paste0("example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
}


z <- apply(simData$y, c(1, 3), function(x){
  ifelse(sum(x) > 0, 1, NA)
})

data <- list(
  y = simData$y,
  windSpeed = simData$covariates$windSpeed,
  elevation = simData$covariates$elevation,
  springPrecipitation = simData$covariates$springPrecipitation,
  sizeOfBeak = simData$covariates$sizeOfBeak,
  z = z
)

constants <- list(
  nyears = nyears,
  nsites = nsites,
  nvisits = nvisits
)
inits <- list(
  alphaPSig = 1,
  betaPSig  = 1,
  alphaPsiSig = 1,
  betaPsiSig  = 1,
  alphaGammaSig = 1,
  betaGammaSig = 1,
  alphaPhi  = 0,
  betaPhi  = 0,
  alphaP =  rnorm(1, mean = 0, sd = 1),
  betaP =  rnorm(1, mean = 0, sd = 1),
  alphaPsi=  rnorm(1, mean = 0, sd = 1),
  betaPsi= rnorm(1, mean = 0, sd = 1),
  colonisationProb = 0.01
)

# NIMBLE model for reduced model
baselineModel <- nimbleModel(dynOccupancyCode,
                             data = data,
                             constants = constants,
                             inits = inits,
                             check = FALSE)

# Fitting the model
baselineModelEst <- spartaNimWeights(model = baselineModel,
                                     latent = "z",
                                     nParFiltRun = numParticles,
                                     mcmc = TRUE,
                                     pfType = pfTypeRun,
                                     pfControl = list(saveAll = TRUE,
                                                      smoothing = TRUE,
                                                      timeIndex = 2),
                                     MCMCconfiguration = list(target = c('alphaPSig', 'betaPSig',
                                                                         'alphaPsiSig', 'betaPsiSig',
                                                                         'alphaPhi', 'betaPhi',
                                                                         'alphaP', 'betaP',
                                                                         'alphaPsi', 'betaPsi',
                                                                         "colonisationProb"
                                     ),
                                     additionalPars = c( "psi.fs"
                                     ),
                                     n.iter = nIterations,
                                     n.chains = nChains,
                                     n.burnin = nBurnin,
                                     n.thin = nThin)
)
#save results
save(baselineModelEst, file= "baselineModel.RData")
