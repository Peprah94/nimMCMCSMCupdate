# This is the Rscript for case study one in the main article.
# The case study is on demographic SSM fitted to Created Tits data
#which is accessed from AHMbook package.

#load Packages needed to run the script
library(nimble)
library(nimbleSMC)
library(nimMCMCSMCupdates)
library(AHMbook)

# Type of PF to run. change if want to fit auxiliary PF
pfTypeRun = "bootstrap"

# Configurations for the MCMC
nIterations = 50000
nBurnin = 20000
nChains = 2
nThin = 6
numParticles = 100

# Define nimble code
popnGrowthModel <- nimbleCode({

  #Priors
  #Model for expected Initial abundance
  alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50) # mean for lambda

  # covariate coefficients
  for(v in 1:3){
    beta.lam[v] ~ dnorm(0,tau.lam)
    }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd

  # Model for immigration-free population growth rate
alpha.gam <- log(mean.gamma)
mean.gamma~dunif(0.9, 1.1)

for(v in 1:3){
  beta.gam[v] ~ dnorm(0,tau.gam)
}

tau.gam <- pow(sd.gam, -2)
sd.gam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd


#model for detection Probability
alpha.p <- logit(mean.p)
mean.p ~ dunif(0,1)
for(v in 1:3){
  beta.p[v] ~ dnorm(0, 0.1)
}


# Model for random immigration
for(t in 1:(T-1)){
  log(rho[t]) <- logrho[t]
  logrho[t] ~ dnorm(0, tau.rho)
}
tau.rho <- pow(sd.rho, -2)
sd.rho ~ T(dnorm(0, 0.5),0.001,) # Half-normal prior for sd

# Likelihood
#state Process
for(i in 1:M){
  #Initial xonditions
  N[i,1] ~ dpois(lambda[i])
  log(lambda[i]) <- loglam[i]
  loglam[i] <- alpha.lam + beta.lam[1]*elev[i] + beta.lam[2]*pow(elev[i], 2) + beta.lam[3]*forest[i]

  # Transition Model
  for(t in 2:T){
    N[i,t] ~ dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1])
    #N[i,t] ~ dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1])
    log(gamma[i, t-1]) <- loggam[i, t-1]
    loggam[i, t-1] <- alpha.gam + beta.gam[1]*elev[i] + beta.gam[2]*pow(elev[i], 2) + beta.gam[3]*forest[i]
  }

  # Observation process
  for(t in 1:T){
    C[i,t] ~ dbin(p[i,t], N[i,t])
    logit(p[i,t]) <- lp[i,t]
    lp[i,t] <- alpha.p + beta.p[1]*date[i,t] + beta.p[2]*pow(date[i,t], 2) + beta.p[3]*dur[i,t]
  }
}

# Derived quantities
for(t in 1:T){
  popindex[t] <- sum(N[1:M,t])
}

for(t in 1:(T-1)){
  gammaX[t] <- popindex[t+1]/popindex[t]
}

})

#load data
data("crestedTit")

#grab counts from 1999 to 2016 and select sites with non-zero counts
C <- as.matrix(crestedTit[ , 6:23])
nzero <- apply(C,1, function(x) {
  ret <-sum(x == 0, na.rm = TRUE)})

sel <- nzero <= 1  # select sites with non-zero counts
newM <- sum(sel)
year <- 1999:2016

# Get data for date and duration
nsite = nrow(C)
nyear = length(year)
datetmp <- as.matrix(crestedTit[, 24:77])
datefull <- array(datetmp, dim = c(nsite, 3, nyear))
durtmp <- as.matrix(crestedTit[, 78:131])
durfull <- array(durtmp, dim = c(nsite, 3, nyear))

#Get mean data and survey per year
date <- apply(datefull, c(1,3), mean, na.rm = TRUE)
dur <- apply(durfull, c(1,3), mean, na.rm = TRUE)
date[date == "NaN"] <- NA
dur[dur =="NaN"] <- NA


#scale the covariates
elev.sc <- standardize(crestedTit$elev)
forest.sc <- standardize(crestedTit$forest)
date.sc <- standardize(date)
date.sc[is.na(date.sc)] <- 0
dur.sc <- standardize(dur)
dur.sc[is.na(dur.sc)] <- 0

# Data for NIMBLE
data = list(C = C[sel, ],
            elev = elev.sc[sel],
            forest = forest.sc[sel],
            date = date.sc[sel,],
            dur = dur.sc[sel,])

# constants for NIMBLE
constants = list(M = newM,
                 T = ncol(C[sel, ]),
                 Nst = max(C[sel,], na.rm = TRUE)+ 10)

Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
Nst = Nst + 10

# set initial values
inits = list(mean.lambda = 30,
             beta.lam = c(0.1128, -0.0859, 0.4087),
             mean.gamma = 0.911,
             beta.gam = rnorm(3,0,1),
             mean.p = 0.3738,
             beta.p = c(-0.167, 0.068, 0.175),
             logrho = c(-0.2510079, -0.4308724, 0.7380799,
                        1.1775382,  0.7144803,  0.2239844, -0.2922550,
                        0.7901613, 0.3417639, -0.6317692, 0.6841965 , 1.0419031,
                        -0.8516481 , 0.1670149 ,0.4615348, -0.5136026,  0.3478865 ),
             sd.rho = 0.852,
             sd.lam = 0.416,
             sd.gam = 0.089,
             N = Nst)

# Define the nimble model
popnGrowthNimModel <- nimbleModel(popnGrowthModel,
                            data = data,
                            constants = constants,
                            inits = inits,
                            check = FALSE)

#####################
# Fit reduced model
##################

iNodePrev = 17 # Last year used to fit the reduced model

data = list(C = C[sel, 1:iNodePrev],
            elev = elev.sc[sel],
            forest = forest.sc[sel],
            date = date.sc[sel,1:iNodePrev],
            dur = dur.sc[sel,1:iNodePrev])

constants = list(M = newM,
                 T = ncol(C[sel, 1:iNodePrev]),
                 Nst = max(C[sel,], na.rm = TRUE)+ 10)

Nst <- C[sel,1:iNodePrev]
Nst[is.na(Nst)] <- 0
Nst = Nst+ 10

inits = list(mean.lambda = 10,
             beta.lam = rnorm(3, 0, 0.5),
             mean.gamma = 1,
             beta.gam = rnorm(3,0,1),
             mean.p = 0.5,
             beta.p = rnorm(3, 0,0.5),
             logrho = rnorm((constants$T -1), 2, 1),
             sd.rho = 1,
             sd.lam = 4,
             sd.gam = 1,
             N = Nst )

 reducedModel <- nimbleModel(popnGrowthModel,
                                data = data,
                                constants = constants,
                                inits = inits,
                                check = FALSE)

 #Fitting MCMC for reduced Model
#
# example2ReducedModelTrue <- spartaNimWeights(model = reducedModel,
#                                              latent = "N",
#                                              nParFiltRun = numParticles,
#                                              mcmc = TRUE,
#                                              pfType = pfTypeRun,
#                                              #pfType = pfTypeRun,
#                                              pfControl = list(saveAll = TRUE,
#                                                               smoothing = TRUE,
#                                                               timeIndex = 2),
#                                              MCMCconfiguration = list(target = c('mean.lambda',
#                                                                                  'beta.lam',
#                                                                                  'mean.gamma',
#                                                                                  'beta.gam',
#                                                                                  'mean.p',
#                                                                                  'beta.p',
#                                                                                  'logrho',
#                                                                                  'sd.rho',
#                                                                                  'sd.gam',
#                                                                                  'sd.lam'
#                                              ),
#                                              additionalPars = c("N", "gammaX", "popindex"
#                                              ),
#                                              n.iter = nIterations,
#                                              n.chains = nChains,
#                                              n.burnin = nBurnin,
#                                              n.thin = nThin)
# )
#
# save(example2ReducedModelTrue, file = "reducedModelResults.RData")
load("reducedModelResults.RData")

###################
# Updated Model
#################

popnGrowthModelUpdated <- nimbleCode({

  #Priors
  #Model for expected Initial abundance
  alpha.lam <- log(mean.lambda)
  mean.lambda ~ dunif(0, 50) # mean for lambda

  # covariate coefficients
  for(v in 1:3){
    beta.lam[v] ~ dnorm(0,tau.lam)
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd

  # Model for immigration-free population growth rate
  alpha.gam <- log(mean.gamma)
  mean.gamma~dunif(0.9, 1.1)

  for(v in 1:3){
    beta.gam[v] ~ dnorm(0,tau.gam)
  }

  tau.gam <- pow(sd.gam, -2)
  sd.gam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd


  #model for detection Probability
  alpha.p <- logit(mean.p)
  mean.p ~ dunif(0,1)
  for(v in 1:3){
    beta.p[v] ~ dnorm(0, 0.1)
  }


  # Model for random immigration
  for(t in 1:(T-1)){
    log(rho[t]) <- logrho[t]
    logrho[t] ~ dnorm(0, tau.rho)
  }
  tau.rho <- pow(sd.rho, -2)
  sd.rho ~ T(dnorm(0, 0.5),0.001,) # Half-normal prior for sd

  # Likelihood
  #state Process
  for(i in 1:M){
    #Initial xonditions
    N[i,1] ~ dpois(lambda[i])
    log(lambda[i]) <- loglam[i]
    loglam[i] <- alpha.lam + beta.lam[1]*elev[i] + beta.lam[2]*pow(elev[i], 2) + beta.lam[3]*forest[i]

    # Transition Model
    for(t in 2:T){
     # N[i,t] ~ T(dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1]), 1, Nst)
      N[i,t] ~ dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1])
      log(gamma[i, t-1]) <- loggam[i, t-1]
      loggam[i, t-1] <- alpha.gam + beta.gam[1]*elev[i] + beta.gam[2]*pow(elev[i], 2) + beta.gam[3]*forest[i]
    }

    # Observation process
    for(t in 1:T){
      C[i,t] ~ dbin(p[i,t], N[i,t])
      logit(p[i,t]) <- lp[i,t]
      lp[i,t] <- alpha.p + beta.p[1]*date[i,t] + beta.p[2]*pow(date[i,t], 2) + beta.p[3]*dur[i,t]
    }
  }

  # Derived quantities
  for(t in 1:T){
    popindex[t] <- sum(N[1:M,t])
  }

  for(t in 1:(T-1)){
    gammaX[t] <- popindex[t+1]/popindex[t]
  }

})

data = list(C = C[sel, ],
            Cst = C[sel, ],
            elev = elev.sc[sel],
            forest = forest.sc[sel],
            date = date.sc[sel,],
            dur = dur.sc[sel,])

constants = list(M = newM,
                 T = ncol(C[sel, ]),
                 Nst = max(C[sel,-c(1:16)], na.rm = TRUE))

Nst <- C[sel,]
Nst[is.na(Nst)] <- 0
Nst = Nst + 10

inits = list(mean.lambda = 30,
             beta.lam = c(0.1128, -0.0859, 0.4087),
             mean.gamma = 0.911,
             beta.gam = rnorm(3,0,1),
             mean.p = 0.3738,
             beta.p = c(-0.167, 0.068, 0.175),
             logrho = c(-0.2510079, -0.4308724, 0.7380799,
                        1.1775382,  0.7144803,  0.2239844, -0.2922550,
                        0.7901613, 0.3417639, -0.6317692, 0.6841965 , 1.0419031,
                        -0.8516481 , 0.1670149 ,0.4615348, -0.5136026,  0.3478865 ),
             sd.rho = 0.852,
             sd.lam = 0.416,
             sd.gam = 0.089,
             N = Nst)


popnGrowthNimModelUpd <- nimbleModel(popnGrowthModelUpdated,
                                  data = data,
                                  constants = constants,
                                  inits = inits,
                                  check = FALSE)

#  popnGrowthModelUpdated <- nimbleCode({
#
#    #Priors
#    #Model for expected Initial abundance
#    alpha.lam <- log(mean.lambda)
#    mean.lambda ~ dunif(0, 50) # mean for lambda
#
#    # covariate coefficients
#    for(v in 1:3){
#      beta.lam[v] ~ dnorm(0,tau.lam)
#    }
#    tau.lam <- pow(sd.lam, -2)
#    sd.lam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd
#
#    # Model for immigration-free population growth rate
#    alpha.gam <- log(mean.gamma)
#    mean.gamma~dunif(0.9, 1.1)
#
#    for(v in 1:3){
#      beta.gam[v] ~ dnorm(0,tau.gam)
#    }
#
#    tau.gam <- pow(sd.gam, -2)
#    sd.gam ~ T(dnorm(0, 1),0.001,) # Half-normal prior for sd
#
#
#    #model for detection Probability
#    alpha.p <- logit(mean.p)
#    mean.p ~ dunif(0,1)
#    for(v in 1:3){
#      beta.p[v] ~ dnorm(0, 0.1)
#    }
#
#
#    # Model for random immigration
#    for(t in 1:(T-1)){
#      log(rho[t]) <- logrho[t]
#      logrho[t] ~ dnorm(0, tau.rho)
#    }
#    tau.rho <- pow(sd.rho, -2)
#    sd.rho ~ T(dnorm(0, 0.5),0.001,) # Half-normal prior for sd
#
#    # Likelihood
#    #state Process
#    for(i in 1:M){
#      #Initial xonditions
#      N[i,1] ~ T(dpois(lambda[i]), 1, Nst)
#      log(lambda[i]) <- loglam[i]
#      loglam[i] <- alpha.lam + beta.lam[1]*elev[i] + beta.lam[2]*pow(elev[i], 2) + beta.lam[3]*forest[i]
#
#      # Transition Model
#      for(t in 2:T){
#        N[i,t] ~ T(dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1]),1, Nst )
#        #N[i,t] ~ dpois(N[i, t-1]*gamma[i, t-1] + rho[t-1])
#        log(gamma[i, t-1]) <- loggam[i, t-1]
#        loggam[i, t-1] <- alpha.gam + beta.gam[1]*elev[i] + beta.gam[2]*pow(elev[i], 2) + beta.gam[3]*forest[i]
#      }
#
#      # Observation process
#      for(t in 1:T){
#        C[i,t] ~ dbin(p[i,t], N[i,t])
#        logit(p[i,t]) <- lp[i,t]
#        lp[i,t] <- alpha.p + beta.p[1]*date[i,t] + beta.p[2]*pow(date[i,t], 2) + beta.p[3]*dur[i,t]
#      }
#    }
#
#    # Derived quantities
#    for(t in 1:T){
#      popindex[t] <- sum(N[1:M,t])
#    }
#
#    for(t in 1:(T-1)){
#      gammaX[t] <- popindex[t+1]/popindex[t]
#    }
#
#  })
#
# data = list(C = C,#[sel, -c(1:iNodePrev-1)],
#             elev = elev.sc,#[sel],
#             forest = forest.sc,#[sel],
#             date = date.sc,#[sel,-c(1:iNodePrev-1)],
#             dur = dur.sc#[sel,-c(1:iNodePrev-1)]
# )
#
# constants = list(M = newM,
#                  T = ncol(C),
#                  Nst = max(data$C, na.rm = TRUE)+ 10)
#
# Nst <- C#[sel,-c(1:iNodePrev-1)]
# Nst[is.na(Nst)] <- 0
#
#
# inits = list(mean.lambda = 10.977,
#              beta.lam = c(0.1128, -0.0859, 0.4087),
#              mean.gamma = 0.911,
#              beta.gam = rnorm(3,0,1),
#              mean.p = 0.3738,
#              beta.p = c(-0.167, 0.068, 0.175),
#              logrho = c(1.0419031,-0.8516481 , 0.1670149 ,0.4615348, -0.5136026),
#              sd.rho = 0.852,
#              sd.lam = 0.416,
#              sd.gam = 0.089,
#              N = Nst + 1)
#
# # Define model for updating process
# newModelUpdated <- nimbleModel(popnGrowthModel,
#                                data = data,
#                                constants = constants,
#                                inits = inits,
#                                check = FALSE)

# Fitting the updated model
#updatedModel <- popnGrowthNimModel$newModel(replicate = TRUE)
example2UpdatedModelTrue <- spartaNimUpdates(model = popnGrowthNimModelUpd , #nimble model
                                             reducedModel = reducedModel,
                                             latent = "N", #latent variable
                                             nParFiltRun = numParticles,
                                             #nParFiltRun = 1000,
                                             pfType = "bootstrap",
                                             extraVars = "logrho",
                                             mcmcScale = 1,
                                             MCMCconfiguration = list(target = c('mean.lambda',
                                                                                 'beta.lam',
                                                                                 'mean.gamma',
                                                                                 'beta.gam',
                                                                                 'mean.p',
                                                                                 'beta.p',
                                                                                 'logrho',
                                                                                 'sd.rho',
                                                                                 'sd.gam',
                                                                                 'sd.lam'
                                             ),
                                             additionalPars = c( "gammaX", "popindex"),
                                             n.iter = (nIterations - nBurnin)/nThin,
                                             n.chains = nChains,
                                             n.burnin = 10,#100,
                                             n.thin = 1),  #saved loglikelihoods from reduced model
                                             postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
                                             propCov1 =  c(0.02, 0.02, 0.02, 0.02) * diag(4),
                                             propCov =  c(rep(0.02, 7), 0.02,0.02, 0.02, 0.02) * diag(11),
                                             pfControl = list(saveAll = TRUE,
                                                              timeIndex = 2,
                                                              smoothing = TRUE,
                                                              #mcmc = TRUE,
                                                              M = 18 - iNodePrev,
                                                              iNodePrev = iNodePrev),
                                             nCores = 1
)


save(example2UpdatedModelTrue, file = "updatedModelResultsBPF.RData")
#
# example2UpdatedModelTrue <- spartaNimUpdates(model = popnGrowthNimModelUpd , #nimble model
#                                              reducedModel = reducedModel,
#                                              latent = "N", #latent variable
#                                              nParFiltRun = numParticles,
#                                              nCores =1,
#                                              #nParFiltRun = 1000,
#                                              pfType = "auxiliary",
#                                              extraVars = "logrho",
#                                              mcmcScale = 1,
#                                              MCMCconfiguration = list(target = c('mean.lambda',
#                                                                                  'beta.lam',
#                                                                                  'mean.gamma',
#                                                                                  'beta.gam',
#                                                                                  'mean.p',
#                                                                                  'beta.p',
#                                                                                  'logrho',
#                                                                                  'sd.rho',
#                                                                                  'sd.gam',
#                                                                                  'sd.lam'
#                                              ),
#                                              additionalPars = c( "gammaX", "popindex"),
#                                              n.iter = (nIterations - nBurnin)/nThin,
#                                              n.chains = nChains,
#                                              n.burnin = 10,#100,
#                                              n.thin = 1),  #saved loglikelihoods from reduced model
#                                              postReducedMCMC = example2ReducedModelTrue,# MCMC summary to use as initial values
#                                             #mcmcScale = 0.02,
#                                               propCov1 =  c(0.02, 0.02, 0.02, 0.02) * diag(4),
#                                              propCov =  c(rep(0.02, 7), 0.02,0.02, 0.02, 0.02) * diag(11),
#                                              pfControl = list(saveAll = TRUE,
#                                                               timeIndex = 2,
#                                                               smoothing = TRUE,
#                                                               #mcmc = TRUE,
#                                                               M = 18 - iNodePrev,
#                                                               iNodePrev = iNodePrev)
# )
#
#
# save(example2UpdatedModelTrue, file = "updatedModelResultsAPF.RData")

#############################
# Fit full model with MCMC
#############################
# baselineModelMCMC <- spartaNimWeights(model = popnGrowthNimModelUpd,
#                                              latent = "N",
#                                              nParFiltRun = numParticles,
#                                              mcmc = TRUE,
#                                              pfType = pfTypeRun,
#                                              #pfType = pfTypeRun,
#                                              pfControl = list(saveAll = TRUE,
#                                                               smoothing = TRUE,
#                                                               timeIndex = 2),
#                                              MCMCconfiguration = list(target = c('mean.lambda',
#                                                                                  'beta.lam',
#                                                                                  'mean.gamma',
#                                                                                  'beta.gam',
#                                                                                  'mean.p',
#                                                                                  'beta.p',
#                                                                                  'logrho',
#                                                                                  'sd.rho',
#                                                                                  'sd.gam',
#                                                                                  'sd.lam'
#                                              ),
#                                              additionalPars = c("N", "gammaX", "popindex"
#                                              ),
#                                              n.iter = nIterations,
#                                              n.chains = nChains,
#                                              n.burnin = nBurnin,
#                                              n.thin = nThin)
# )
#
# save(baselineModelMCMC, file = "baselineModelMCMCResults.RData")
#
# #############################
# # Fit full model with SMC
# #############################
#
# baselineModelSMC <- nimMCMCSMCupdates::baselineSpartaEstimation(model = popnGrowthNimModel,
#                                                                 latent = "N",
#                                                                 pfType = pfTypeRun,
#                                                                 nParFiltRun = numParticles,
#                                                                 pfControl = list(saveAll = TRUE,
#                                                                                  smoothing = FALSE,
#                                                                                  timeIndex = 2),
#                                                                 MCMCconfiguration = list(target = c('mean.lambda',
#                                                                                                     'beta.lam',
#                                                                                                     'mean.gamma',
#                                                                                                     'beta.gam',
#                                                                                                     'mean.p',
#                                                                                                     'beta.p',
#                                                                                                     'logrho',
#                                                                                                     'sd.rho',
#                                                                                                     'sd.gam',
#                                                                                                     'sd.lam'
#                                                                 ),
#                                                                 additionalPars = c("N", "gammaX", "popindex"
#                                                                 ),
#                                                                 n.iter = nIterations,
#                                                                 n.chains = nChains,
#                                                                 n.burnin = nBurnin,
#                                                                 n.thin = nThin
#                                                                 ))
#
# save(baselineModelSMC, file = "baselineModelSMCResults.RData")
#
