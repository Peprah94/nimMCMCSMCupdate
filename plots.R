# Results file sizes are large. We have attached the link to access our results
# https://www.dropbox.com/sh/9gdf426e2sw5yas/AADqlpM6jxvqRQMmccku3WxTa?dl=0


# Load data
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)
library(egg)


# set directory to the path the data was downloaded

###################
# Linear Gaussian SSM
####################

load("linearGaussianSSM/simulatedDataEx1.RData")

# Need to format the data from Example One
#write a function to do that
formatEx1 <- function(data){

  result <- lapply(data, function(x){
    BMC <- list(x[[1]], x[[2]], x[[3]])
    BSC <- list(x[[4]], x[[5]], x[[6]])
    putTogether <- list(BMC, BSC)
    others <- x[7:26]
    ret <- c(putTogether,
             others)
    return(ret)
  })

  return(result)
}

# Estimate root mean square errr
averageRMSE <- function(x,y){
  ret <- (mean((x-y)))
  return(ret)
}

## Load data for Aux PF
load("linearGaussianSSM/bootstrapPF/estimates1.RData")
pf1 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)
load("linearGaussianSSM/auxiliaryPF/estimates2.RData")
pf2 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)
load("linearGaussianSSM/bootstrapPF/estimates3.RData")
pf3 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)
load("linearGaussianSSM/auxiliaryPF/estimates4.RData")
pf4 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)
load("linearGaussianSSM/bootstrapPF/estimates5.RData")
pf5 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)
load("linearGaussianSSM/auxiliaryPF/estimates6.RData")
pf6 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)

## Put all results together

results1 <- c(pf1,
              pf2,
             pf3,
             pf4,
              pf5,
              pf6)

rm(pf1,
  pf2,
  pf3,
  pf4,
   pf5,
   pf6)


# Flatten list
results <- purrr::flatten(results1)


## Extract model Parameters

nodeNames = c("ahat", "chat")

nSamps <- 90 # Number of replicate samples

#estimate the bias of model parameter
biasModelPars <- lapply(results, function(x){
  red <- x[[2]]$all.chains
  output <-  red[rownames(red)[rownames(red) %in% nodeNames],"Median"]
  return(output)
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(t = rep(c(rep(50, 3),
                          rep(c(49, 45, 20, 10, 5), 3)), nSamps),
                model = rep(
                  c(rep(c("BBSC", "ABSC","BMC",
                          rep("RMC", 5),
                          rep("BUMC", 5),
                          rep("AUMC", 5)), nSamps))
                ))

## Extract results for auxiliary PF and plot
biasModelPars%>%
  ggplot(., aes(x = as.factor(t), y = ahat, fill = model))+
  geom_boxplot()


# Fill colors for plotting
 # Fill colors for plotting
  fill.colors <- c("BMC" = "#FF6600",
                   "ABSC" = "#FF3300",
                   "RMC" = "#00CCFF",
                   "BBSC" = "#E69F00",
                   "AUMC" = "#0033FF",
                   "BUMC" = "#3399FF",
                   "truth" = "black")

# extract results
biasModelPars1 <- biasModelPars%>%
  dplyr::group_by(t,model )%>%
  dplyr::summarise(across(c(ahat,chat), median))%>%
  reshape2::melt(., id.vars = c("t", "model"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))


#plot results
biasModelParsAhat <- biasModelPars1%>%
  filter(variable %in% c("ahat"))%>%
  filter( !t %in% c("45","50"))%>%
  ggplot(.,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][1], linetype = "dashed", col = "#FF3300", linewidth = 1)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][2], linetype = "dotdash", col = "#E69F00", linewidth = 1)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("ahat"),4][3], linetype = "dotted", col = "#FF6600", linewidth = 1)+
  facet_wrap(~variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(c( "ahat" = "a",
                                       #"bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = fill.colors)+
  xlab("")+
  ylab("Bias")

biasModelParsChat <- biasModelPars1%>%
  filter(variable %in% c("chat"))%>%
  filter( !t %in% c("45","50"))%>%
  ggplot(.,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][1], linetype = "dashed", col = "#FF3300", linewidth = 1)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][2], linetype = "dotdash", col = "#E69F00", linewidth = 1)+
  geom_hline(yintercept = biasModelPars1[biasModelPars1$t==50 & biasModelPars1$variable %in%c("chat"),4][3], linetype = "dotted", col = "#FF6600", linewidth = 1)+
  facet_wrap(~variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(c( "ahat" = "a",
                                       #"bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab(" ")

### Extract RMSE results and plot
# Estimating MCSE
mcseEstimates <- lapply(results, function(x){mcmcse::mcse.mat(as.matrix(rbind(x[[1]]$chain1,
                                                                              x[[1]]$chain2))[, c("a", "c")],
                                                              method = "bm",
                                                              #size = bacthSize,
                                                              g = NULL)%>%
    as.data.frame()%>%
    dplyr::select(se)%>%
    t()})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(t = rep(c(rep(50, 3),
                          rep(c(49, 45, 20, 10, 5), 3)), nSamps),
                model = rep(
                  c(rep(c("BBSC", "ABSC","BMC",
                          rep("RMC", 5),
                          rep("BUMC", 5),
                          rep("AUMC", 5)), nSamps))
                ))%>%
  reshape2::melt(., id.vars = c("t", "model"))%>%
  dplyr::group_by(model, t, variable)%>%
  #dplyr::group_by(t)%>%
  dplyr::summarise(value = median(value))

#select for auxPF
mcseEstimatesPars1 <- mcseEstimates%>%
  # filter(model %in% c("ABMC", "ABSC",
  #                     "ARMC",
  #                     "AUMC", "BUMC", "BBSC"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

mcseModelParsA <- mcseEstimatesPars1%>%
  filter( !t %in% c("45","50"))%>%
  filter(variable %in% c("a"))%>%
  ggplot(data = .,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
 geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[1], linetype = "dashed", col = "#FF3300",linewidth = 1)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[2], linetype = "dotdash", col = "#E69F00",linewidth = 1)+
  geom_hline(yintercept =c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("a"),4]$value)[3], linetype = "dotted", col = "#FF6600",linewidth = 1)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab("MCSE")

# Extract bootstrap PF
mcseModelParsB <- mcseEstimatesPars1%>%
  filter( !t %in% c("45","50"))%>%
  filter(variable %in% c("c"))%>%
  ggplot(data = .,
         mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
  geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[1], linetype = "dashed", col = "#FF3300",linewidth = 1)+
 geom_hline(yintercept = c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[2], linetype = "dotdash", col = "#E69F00",linewidth = 1)+
  geom_hline(yintercept =c(mcseEstimatesPars1[mcseEstimatesPars1$t==50 & mcseEstimatesPars1$variable %in%c("c"),4]$value)[3], linetype = "dotted", col = "#FF6600",linewidth = 1)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = fill.colors)+
  xlab(" ")+
  ylab(" ")

## Results for latent state distribution
# auxiliary latent states
allPF <- results1
auxLatentEst <- lapply(seq_along(allPF), function(i){

  x <- allPF[[i]]
  y <- simData[[i]]$x

  BBSC <- averageRMSE(x[[1]][[2]]$all.chains[!rownames(x[[1]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  ABSC <- averageRMSE(x[[2]][[2]]$all.chains[!rownames(x[[2]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  BMC <- averageRMSE(x[[3]][[2]]$all.chains[!rownames(x[[3]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)

  # BUMC
  BUMC_49 <- averageRMSE(c(x[[9]][[2]]$all.chains[!rownames(x[[9]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_45 <- averageRMSE(c(x[[10]][[2]]$all.chains[!rownames(x[[10]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_20 <- averageRMSE(c(x[[11]][[2]]$all.chains[!rownames(x[[11]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_10 <- averageRMSE(c(x[[12]][[2]]$all.chains[!rownames(x[[12]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  BUMC_5 <- averageRMSE(c(x[[13]][[2]]$all.chains[!rownames(x[[13]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  # AUMC
  AUMC_49 <- averageRMSE(c(x[[14]][[2]]$all.chains[!rownames(x[[14]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_45 <- averageRMSE(c(x[[15]][[2]]$all.chains[!rownames(x[[15]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_20 <- averageRMSE(c(x[[16]][[2]]$all.chains[!rownames(x[[16]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_10 <- averageRMSE(c(x[[17]][[2]]$all.chains[!rownames(x[[17]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  AUMC_5 <- averageRMSE(c(x[[18]][[2]]$all.chains[!rownames(x[[18]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"]), y)

  ret <- data.frame(BBSC, ABSC, BMC,
                    BUMC_49, BUMC_45,
                    BUMC_20, BUMC_10,
                    BUMC_5,AUMC_49, AUMC_45,
                    AUMC_20, AUMC_10,
                    AUMC_5
  )
  return(ret)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt()%>%
  dplyr::group_by(variable)%>%
  dplyr::summarise(mean = median(value),
                   sd = sd(value))%>%
  dplyr::mutate(t = c(rep(50, 3),
                      rep(c(49, 45, 20, 10, 5), 2)),
                model = c("BBSC", "ABSC","BMC",
                          rep("BUMC", 5),
                          rep("AUMC", 5))
  )%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

#Plot error plot
errorModelParsAuxPF <- auxLatentEst%>%
  filter( !t %in% c("45","50"))%>%
  ggplot(data = ., mapping = aes(x = as.factor(t), y = mean, col = model))+
  geom_point(position = position_dodge(width = 0.7), size = 2)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),size=1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = c(auxLatentEst[auxLatentEst$t==50,2]$mean)[1], linetype = "dashed", col = "#FF3300",linewidth = 1)+
 geom_hline(yintercept = c(auxLatentEst[auxLatentEst$t==50,2]$mean)[2], linetype = "dotdash", col = "#E69F00",linewidth = 1)+
  geom_hline(yintercept =c(auxLatentEst[auxLatentEst$t==50,2]$mean)[3], linetype = "dotted", col = "#FF6600",linewidth = 1)+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(legend.title = element_blank())+
  scale_color_manual(values = fill.colors)+
  xlab( " ")+
  ylab(expression(paste("Median ", "\u00B1", " SD")))


#Save results
biasPlots <- ggpubr::ggarrange(biasModelParsAhat,
                       biasModelParsChat,
                       mcseModelParsA,
                       mcseModelParsB,
                       errorModelParsAuxPF,
                       nrow = 3,
                       ncol = 2,
                       widths = c(1,1,2),
                       # labels = c("A)", "B)", "C"),
                       font.label = list(size = 25),
                       common.legend = TRUE)

auxPFpLots <- annotate_figure(biasPlots,
                              bottom = text_grob(expression(t)))


## Put plots together

ggsave(filename = "Figures/auxPFplotsAllBaseline.png",
       plot = auxPFpLots,
       width = 6,
       height = 6,
       units = "in"
       )


auxPFSubset <- allPF[[2]][c(1,2, 3, 9, 13, 14, 18)]

aPlots <- lapply(auxPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(legend.position = "bottom")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots[1:3],
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"))

aPlots1 <- annotate_figure(fig,
                           left = text_grob("Value", rot = 90),
                           bottom = text_grob("Iterations"))

ggsave(filename = "Figures/baselineModels.png",
       plot = aPlots1,
       width = 6,
       height = 7,
       units = "in")

fig <- ggpubr::ggarrange(plotlist = aPlots[4:7],
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"))
aPlots2 <- annotate_figure(fig,
                           left = text_grob("Value", rot = 90),
                           bottom = text_grob("Iterations"))

ggsave(filename = "Figures/updatedModels.png",
       plot = aPlots2,
       width = 8,
       height = 8,
       units = "in")

rm(list=ls())

######################
# Simulation study 2: Dynamic occupancy model
###################

library(dplyr)
#library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)
library(ggplot2)
theme_set(
  theme_light() + theme(legend.position = "top")
)

#fill.colors2 <- # Fill colors for plotting
fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")

load("dynamicOccupancyModel/auxiliaryPF/simDataDynamicOccupancy.RData")

#load baseline model results
load("dynamicOccupancyModel/auxiliaryPF/baselineModel.RData")
BMC <- baselineModelEst


#######
# Reduced Models... Same for all No of particles
#######
load("dynamicOccupancyModel/auxiliaryPF/example4ReducedBootstrapTrueM10Ind29.RData")
RMC <- example2ReducedModelTrue

load("dynamicOccupancyModel/auxiliaryPF/example4ReducedBootstrapTrueM10Ind25.RData")
RMC2 <- example2ReducedModelTrue
########
# M = 10
#####
numParticles <- 10
iNodePrev <- 29
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC <- example2UpdatedModelTrue

iNodePrev <- 25
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC2 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC2 <- example2UpdatedModelTrue

#########
# M = 25
######
# Note that the reduced model is the same irrespective of M
numParticles <- 25
iNodePrev <- 29
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC3 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC3 <- example2UpdatedModelTrue

iNodePrev <- 25
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC4 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC4 <- example2UpdatedModelTrue


#########
# M = 50
######
# Note that the reduced model is the same irrespective of M
numParticles <- 50
iNodePrev <- 29
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC5 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC5 <- example2UpdatedModelTrue

iNodePrev <- 25
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC6 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC6 <- example2UpdatedModelTrue


#########
# M = 100
######
# Note that the reduced model is the same irrespective of M
numParticles <- 100
iNodePrev <- 29
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC7 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC7 <- example2UpdatedModelTrue

iNodePrev <- 25
load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedBootstrapTrueM",numParticles,"Ind",iNodePrev,".RData"))
BUMC8 <- example2UpdatedModelTrue

load(paste0("dynamicOccupancyModel/auxiliaryPF/example4UpdatedAuxiliaryTrueM",numParticles,"Ind",iNodePrev,".RData"))
AUMC8 <- example2UpdatedModelTrue


#remove individual data results
rm(baselineModelEst,
   example2ReducedModelTrue,
   example2UpdatedModelTrue)

# T = 29
allModels <- list(BMC,
                  RMC,
                  BUMC,
                  AUMC,
                  BUMC3,
                  AUMC3,
                  BUMC5,
                  AUMC5,
                  BUMC7,
                  AUMC7)

# T = 25
allModels1000 <- list(BMC,
                      RMC2,
                      BUMC2,
                      AUMC2,
                      BUMC4,
                      AUMC4,
                      BUMC6,
                      AUMC6,
                      BUMC8,
                      AUMC8)


## Extract realised occupancy and plot
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Mean"]
})

ret1000 <- lapply(allModels1000, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Mean"]
})


allOccResults <- c(ret[-2],
                   ret1000[-2])%>%
  do.call("c", .)%>%
  cbind(.,
        rep(simData$occSites, 18),
        rep(1:30,18),
        rep(rep(c("BMC", "BUMC", "AUMC", "BUMC", "AUMC", "BUMC", "AUMC", "BUMC", "AUMC"), each = 30), 2),
        rep(c("t=29", "t=25"), each = 270),
        rep(rep(c("0", "10", "10", "25", "25", "50", "50", "100", "100"), each = 30), 2))%>%
  as.data.frame()

colnames(allOccResults) <- c("estimates",
                             "truth",
                             "year",
                             "model",
                             "type",
                             "nParticles")


#calculate correlation and annual growth rate
# New facet label names for supp variable
supp.labs <- c("Correlation", "Bias(Growth rate)")
names(supp.labs) <- c("corr", "perChange")

metricsToPlot <- allOccResults%>%
  group_by(model, type, nParticles)%>%
  mutate_at(., c("estimates","truth"), as.numeric)%>%
  mutate(corr = ifelse(type == "t=29",
                       cor(estimates[29:30], truth[29:30]),
                       cor(estimates[26:30], truth[26:30])),
         perChange = (((estimates[30]/estimates[1])^(1/20)-1)*100 - ((truth[30]/truth[1])^(1/30)-1)*100),
         biasChange = ((estimates[30] - truth[30])*100)
  )%>%
  summarise(corr = mean(corr),
            perChange = mean(perChange),
            biasChange = mean(biasChange))%>%
  ungroup()%>%
  reshape2::melt(., id.vars = c("model", "type", "nParticles"))

interceptVals <- metricsToPlot %>%
  dplyr::filter(variable %in% c("corr"))%>%
  dplyr::filter(nParticles %in% c("0"))

corrPlots <- metricsToPlot %>%
  dplyr::filter(variable %in% c("corr"))%>%
  dplyr::filter(!nParticles %in% c("0"))%>%
  ggplot(. , aes(x=as.factor(as.numeric(nParticles)), y = value, col = model))+
  geom_point( aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ type, labeller = labeller(variable = supp.labs),scales = "free_y")+
  scale_color_manual(values = fill.colors)+
  theme(
    strip.text.x = element_text( color = "black"
    )
  )+
  #theme_bw()+
  xlab("")+
  ylab("Correlation")


interceptVals <- metricsToPlot %>%
  dplyr::filter(variable %in% c("biasChange"))%>%
  dplyr::filter(nParticles %in% c("0"))

biasPlots <- metricsToPlot %>%
  dplyr::filter(variable %in% c("biasChange"))%>%
  dplyr::filter(!nParticles %in% c("0"))%>%
  ggplot(. , aes(x=as.factor(as.numeric(nParticles)), y = value, col = model))+
  geom_point( aes(shape = model), position = position_dodge(width = 0.9), size = 2)+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ type, labeller = labeller(variable = supp.labs),scales = "free_y")+
  scale_color_manual(values = fill.colors)+
  theme(
    strip.text.x = element_text( color = "black"
    )
  )+
  #theme_bw()+
  xlab("")+
  ylab("Bias")


fig <- ggpubr::ggarrange(corrPlots,
                         biasPlots,
                         ncol = 1,
                         common.legend = TRUE)

occPlots <- annotate_figure(fig,
                           #left = text_grob("Value", rot = 90),
                           bottom = text_grob("Number of particles (M)"))




ggsave(filename = "Figures/psifs.png",
       plot = occPlots,
       width = 6,
       height = 4,
       units = "in")

## Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})

ret1000 <- lapply(allModels1000, function(x){
  x[[3]]
})

## Model Parameters and plot them
pars <- c('alphaPSig', 'betaPSig',
          'alphaPsiSig', 'betaPsiSig',
          'alphaPhi', 'betaPhi',
          'alphaP', 'betaP',
          'alphaPsi', 'betaPsi',
          "colonisationProb")

pars <- sort(pars)
## Mean and SD of pars
## t = 48
# Mean of pars
retMean29 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, "Mean"]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()

retMean29 <- retMean29[sort(rownames(retMean29)),]%>%


  mutate(truth = c(5.44, -2, 7.40, 2, 2, 3.74, 1.5, -1.07,3,3, 0.05),
  metric = "mean",
  parameters = pars,
  t = 29)
colnames(retMean29)[1:10] <- c('BMC',
                         'RMC',
                         'BUMC_M10', 'AUMC_M10',
                         'BUMC_M25', 'AUMC_M25',
                         'BUMC_M50', 'AUMC_M50',
                         'BUMC_M100', 'AUMC_M100')
# sd of pars
retSD29 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, 3]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()

retSD29 <- retSD29[sort(rownames(retSD29)),]%>%
  mutate(truth = NA,
         metric = "sd",
         parameters = pars,
         t = 29)
colnames(retSD29)[1:10] <- c('BMC',
                             'RMC',
                             'BUMC_M10', 'AUMC_M10',
                             'BUMC_M25', 'AUMC_M25',
                             'BUMC_M50', 'AUMC_M50',
                             'BUMC_M100', 'AUMC_M100')


## t = 25
retMean25 <- lapply(allModels1000, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, "Mean"]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()

retMean25 <- retMean25[sort(rownames(retMean25)),]%>%
  mutate(truth = c(5.44, -2, 7.40, 2, 2, 3.74, 1.5, -1.07,3,3, 0.05),
         metric = "mean",
         parameters = pars,
         t = 25)
colnames(retMean25)[1:10] <- c('BMC',
                               'RMC',
                               'BUMC_M10', 'AUMC_M10',
                               'BUMC_M25', 'AUMC_M25',
                               'BUMC_M50', 'AUMC_M50',
                               'BUMC_M100', 'AUMC_M100')
# sd of pars
retSD25 <- lapply(allModels1000, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, 3]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()

retSD25 <- retSD25[sort(rownames(retSD25)),]%>%
  mutate(truth = NA,
         metric = "sd",
         parameters = pars,
         t = 25)
colnames(retSD25)[1:10] <- c('BMC',
                             'RMC',
                             'BUMC_M10', 'AUMC_M10',
                             'BUMC_M25', 'AUMC_M25',
                             'BUMC_M50', 'AUMC_M50',
                             'BUMC_M100', 'AUMC_M100')

# Put both estimates together
parEsts <- rbind(retMean29,
                 retSD29,
                 retMean25,
                 retSD25)
write.csv(parEsts, file = "Figures/occupancyParsEst.csv", row.names = FALSE)

## Effective sample size
ESSret <- lapply(allModels, function(x) {
  ret <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(ESSret)[1:10] <- c('BMC_0',
                           'RMC',
                           'BUMC_10', 'AUMC_10',
                           'BUMC_25', 'AUMC_25',
                           'BUMC_50', 'AUMC_50',
                           'BUMC_100', 'AUMC_100')

parsToPlot <- c( 'betaPsiSig', 'betaPsi', "colonisationProb")

supp.labs <- c("expression(sigma[alpha^p])",
               "expression(sigma[beta^p])",
               "expression(sigma[alpha^psi])",
               "expression(sigma[beta^psi])",
               "expression(alpha^phi)",
               "expression(beta^phi)",
               "expression(alpha^p)",
               "expression(beta^p)",
               "expression(alpha^Psi)",
               "expression(beta^Psi)",
               "expression(gamma)")
names(supp.labs) <- c('alphaPSig', 'betaPSig',
                      'alphaPsiSig', 'betaPsiSig',
                      'alphaPhi', 'betaPhi',
                      'alphaP', 'betaP',
                      'alphaPsi', 'betaPsi',
                      "colonisationProb")

ESSret29 <- ESSret%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  filter(!variable %in% c("RMC"))%>%
  filter(parameters %in%parsToPlot)%>%
  tidyr::separate_wider_delim(., variable, delim = "_", names = c("model", "M"))%>%
dplyr::mutate(Facets = rep(c("beta^Psi",
                             "sigma[beta^psi]",
                             "gamma"), 9))

interceptVals <- ESSret29%>%
  dplyr::filter(M %in% c("0"))

essFigPlots <- ESSret29%>%
  filter(!M %in% c("0"))%>%
  ggplot(., aes(x =as.factor(as.numeric(M)), y = value, col = model))+
  geom_point(aes(shape = model), size = 2,position = position_dodge(width = 0.9))+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ Facets, labeller = label_parsed)+
  scale_color_manual(values = fill.colors)+
  xlab("")+
  ylab("ESS")+
  theme(
    strip.text.x = element_text( color = "black"
    )
  )+
  labs(col = " ", shape = " ")


ESSret1000 <- lapply(allModels1000, function(x) {
  ret <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(ESSret1000)[1:10] <- c('BMC_0',
                                'RMC',
                                'BUMC_10', 'AUMC_10',
                                'BUMC_25', 'AUMC_25',
                                'BUMC_50', 'AUMC_50',
                                'BUMC_100', 'AUMC_100')

ESSret25 <- ESSret1000%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  filter(!variable %in% c("RMC"))%>%
  filter(parameters %in%parsToPlot)%>%
  tidyr::separate_wider_delim(., variable, delim = "_", names = c("model", "M"))%>%
  dplyr::mutate(Facets = rep(c("beta^Psi",
                               "sigma[beta^psi]",
                               "gamma"), 9))

interceptVals <- ESSret25%>%
  dplyr::filter(M %in% c("0"))

essFigPlots25 <- ESSret25%>%
  filter(!M %in% c("0"))%>%
  ggplot(., aes(x = as.factor(as.numeric(M)), y = value, col = model))+
  geom_point(aes(shape = model), size = 2,position = position_dodge(width = 0.9))+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ Facets, labeller = label_parsed)+
scale_color_manual(values = fill.colors)+
  theme(
  strip.text.x = element_text( color = "black"
  )
)+
  xlab("")+
  ylab("ESS")+
  labs(col = " ", shape = " ")


## Efficiency

EfficiencyRet <- lapply(allModels, function(x) {
  nDim <- length(x$timeRun)

  ess <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "secs")
  }else{
    as.numeric(x$timeRun$all.chains, units = "secs")
  }
print(timesRet)
  ret <- ess/timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(EfficiencyRet)[1:10] <- c('BMC_0',
                                   'RMC',
                                   'BUMC_10', 'AUMC_10',
                                   'BUMC_25', 'AUMC_25',
                                   'BUMC_50', 'AUMC_50',
                                   'BUMC_100', 'AUMC_100')



EfficiencyRet29 <- EfficiencyRet%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  filter(!variable %in% c("RMC"))%>%
  filter(parameters %in%parsToPlot)%>%
  tidyr::separate_wider_delim(., variable, delim = "_", names = c("model", "M"))%>%
  dplyr::mutate(Facets = rep(c("beta^Psi",
                               "sigma[beta^psi]",
                               "gamma"), 9))

interceptVals <- EfficiencyRet29%>%
  dplyr::filter(M %in% c("0"))

EfficiencyRet1 <- EfficiencyRet29%>%
  filter(!M %in% c("0"))%>%
  ggplot(., aes(x = as.factor(as.numeric(M)), y = value, col = model))+
  geom_point(aes(shape = model), size = 2,position = position_dodge(width = 0.9))+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ Facets, labeller = label_parsed)+
scale_color_manual(values = fill.colors)+  theme(
  strip.text.x = element_text( color = "black"
  )
)+
  xlab("")+
  ylab("Efficiency")+
  #labs(title = "A) M = 100")+
  labs(col = " ", shape = " ")


EfficiencyRet1000 <- lapply(allModels1000, function(x) {
  nDim <- length(x$timeRun)

  ess <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "secs")
  }else{
    as.numeric(x$timeRun$all.chains, units = "secs")
  }
print(timesRet)
  ret <- ess/timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(EfficiencyRet1000)[1:10] <- c('BMC_0',
                                       'RMC',
                                       'BUMC_10', 'AUMC_10',
                                       'BUMC_25', 'AUMC_25',
                                       'BUMC_50', 'AUMC_50',
                                       'BUMC_100', 'AUMC_100')


EfficiencyRet25 <- EfficiencyRet1000%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  filter(!variable %in% c("RMC"))%>%
  filter(parameters %in%parsToPlot)%>%
  tidyr::separate_wider_delim(., variable, delim = "_", names = c("model", "M"))%>%
  dplyr::mutate(Facets = rep(c("beta^Psi",
                               "sigma[beta^psi]",
                               "gamma"), 9))

interceptVals <- EfficiencyRet25%>%
  dplyr::filter(M %in% c("0"))

EfficiencyRet2 <- EfficiencyRet25%>%
  filter(!M %in% c("0"))%>%
  ggplot(., aes(x = as.factor(as.numeric(M)), y = value, col = model))+
  geom_point(aes(shape = model), size = 2,position = position_dodge(width = 0.9))+
  geom_hline(data = interceptVals, aes(yintercept = value),linetype = "dashed", col = "#FF3300")+
  facet_wrap( ~ Facets, labeller = label_parsed)+
  #theme_classic()+
  # scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]),
  #                             'betaPSig'= expression(sigma[beta^p]),
  #                             'alphaPsiSig'= expression(sigma[alpha^psi]),
  #                             'betaPsiSig'= expression(sigma[beta^psi]),
  #                             'alphaPhi'= expression(alpha^phi),
  #                             'betaPhi'= expression(beta^phi),
  #                             'alphaP'= expression(alpha^p),
  #                             'betaP'= expression(beta^p),
  #                             'alphaPsi'= expression(alpha^Psi),
  #                             'betaPsi'= expression(beta^Psi),
#                             "colonisationProb"= expression(gamma)))+
scale_color_manual(values = fill.colors)+  theme(
  strip.text.x = element_text( color = "black"
  )
)+
  xlab("")+
  ylab("Efficiency")+
  #labs(title = "A) M = 100")+
  labs(col = " ", shape = " ")


## Put all plots together
fig <- ggpubr::ggarrange(essFigPlots,
                         EfficiencyRet1,
                         ncol = 1, nrow = 2,
                         common.legend = TRUE,
                         legend = "top")

fig1 <- ggpubr::ggarrange(essFigPlots25,
                         EfficiencyRet2,
                         ncol = 1, nrow = 2,
                         common.legend = TRUE,
                         legend = "top")

effESSPlot <- annotate_figure(fig,
                              bottom = text_grob("Number of particles (M)"))

effESSPlot2 <- annotate_figure(fig1,
                              bottom = text_grob("Number of particles (M)"))


ggsave(filename = "Figures/essEffPlotEx2t29.png",
       plot = effESSPlot,
       width = 6,
       height = 5,
       units = "in")

ggsave(filename = "Figures/essEffPlotEx2t25.png",
       plot = effESSPlot2,
       width = 6,
       height = 5,
       units = "in")

# Convergence of parameters
# Auxiliary PF
aPlots <- lapply(allModels[1:4], function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[1:4])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExA.png",
       plot = aPlots,
       width = 6,
       height = 7,
       units = "in")


aPlots <- lapply(allModels[1:4], function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[5:8])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots[1:4],
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExB.png",
       plot = aPlots,
       width = 8,
       height = 8,
       units = "in")

aPlots <- lapply(allModels, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[9:11])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D"))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90),
                          bottom = text_grob("Iterations"))
ggsave(filename = "Figures/convergenceSecExC.png",
       plot = aPlots,
       width = 8,
       height = 8,
       units = "in")

# Extract the time in fitting the models.
timeRet29 <- lapply(allModels, function(x) {
  nDim <- length(x$timeRun)

  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "mins")
  }else{
    as.numeric(x$timeRun$all.chains, units = "mins")
  }

  ret <- timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)
colnames(timeRet29)[1:10] <- c('BMC_0',
                               'RMC',
                               'BUMC_10', 'AUMC_10',
                               'BUMC_25', 'AUMC_25',
                               'BUMC_50', 'AUMC_50',
                               'BUMC_100', 'AUMC_100')

timeRet25 <- lapply(allModels1000, function(x) {
  nDim <- length(x$timeRun)

  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "mins")
  }else{
    as.numeric(x$timeRun$all.chains, units = "mins")
  }

  ret <- timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)
colnames(timeRet25)[1:10] <- c('BMC_0',
                               'RMC',
                               'BUMC_10', 'AUMC_10',
                               'BUMC_25', 'AUMC_25',
                               'BUMC_50', 'AUMC_50',
                               'BUMC_100', 'AUMC_100')

allTimes <- rbind(timeRet29,
                  timeRet25)%>%
  as.data.frame()%>%
  dplyr::mutate(t = c("29", "25"))

write.csv(allTimes, file = "Figures/simOccupancyTimesRun.csv", row.names = FALSE)

rm(list=ls())

###########################################
# CASE STUDY ONE: Population Demographic model
##########################################

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

# Fill colors for plotting
fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")


load("demographicSSM/bootstrapPF/baselineModelMCMCResults.RData")
BMC <- baselineModelMCMC
rm(baselineModelMCMC)

load("demographicSSM/bootstrapPF/reducedModelResults.RData")
RMC <- example2ReducedModelTrue
rm(example2ReducedModelTrue)

load("demographicSSM/bootstrapPF/updatedModelResultsBPF.RData")
BUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)

load("demographicSSM/bootstrapPF/updatedModelResultsAPF.RData")
AUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)


# Put models together
allModels <- list(BMC,
                  RMC,
                  BUMC,
                  AUMC)

# Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})


# Extra growth rate and its standard deviation
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  c(x[[2]]$all.chains[extractNames, "Mean"], NA)
})


retSD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  c(x[[2]]$all.chains[extractNames, 3], "NA")
})


# Extract population index and its standard deviation
ret1 <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})

ret1SD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})

allOccResults <- c(ret1[-2],
                   ret[-2])%>%
  do.call("c", .)%>%
  cbind(.,
        rep(1999:2016,6),
        rep(rep(c("BMC", "BUMC", "AUMC"), each = 18), 2),
        rep(c("popnIndex", "gamma"), each = 54),
        "mean")%>%
  as.data.frame()

colnames(allOccResults) <- c("estimates",
                             "year",
                             "model",
                             "metric",
                             "pars")

allOccResultsSD <- c(ret1SD[-2],
                   retSD[-2])%>%
  do.call("c", .)%>%
  cbind(.,
        rep(1999:2016,6),
        rep(rep(c("BMC", "BUMC", "AUMC"), each = 18), 2),
        rep(c("popnIndex", "gamma"), each = 54),
        "sd")%>%
  as.data.frame()

colnames(allOccResultsSD) <- c("estimates",
                             "year",
                             "model",
                             "metric",
                             "pars")

resultsToPresent <- rbind(allOccResults,
                          allOccResultsSD)%>%
  filter(year %in% c("2015"))

write.csv(resultsToPresent , file = "Figures/demograhicInterAnnualEsts.csv", row.names = FALSE)

#calculate correlation and annual growth rate
popBase <- ret1[[1]][17:18]
gammaBase <- ret[[1]][16:17]
fig <- allOccResults%>%
  group_by(model, metric)%>%
  mutate_at(., c("estimates"), as.numeric)%>%
  mutate(corr = ifelse(metric == "popnIndex",
                       cor(estimates[17:18], popBase),
                       cor(estimates[16:17], gammaBase)),
         perChange = (((estimates[17]/estimates[1])^(1/50)-1))*100)%>%
  summarise(corr = mean(corr, na.rm = TRUE),
            perChange = mean(perChange, na.rm=TRUE))%>%
  ungroup()%>%
  reshape2::melt(., id.vars = c("model", "metric"))%>%
  ggplot(. , aes(x=metric, y = value, col = model))+
  geom_point( aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  facet_wrap( ~ variable)+
  scale_color_manual(values = fill.colors)+
  xlab("Metric")+
  ylab("Estimate")


ggsave(filename = "Figures/example2.png",
       plot = fig,
       width = 18,
       height = 10,
       dpi = 100)


pars <- c('mean.lambda',
          'beta.lam',
          'mean.gamma',
          'beta.gam',
          'mean.p',
          'beta.p',
          'logrho',
          'sd.rho',
          'sd.gam',
          'sd.lam'
)

pars <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[rownames(x[[2]]$all.chains) %in%pars]
  x[[2]]$all.chains[extractNames, 1]
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  mutate(model = c("BBMC",
                   "BRMC",
                   "BUMC",
                   "AUMC"))

write.csv(pars, file = "Figures/demograhicPars.csv", row.names = FALSE)

rm(list = ls())

#############
# Case Study 2: Sparta Model
###############

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

fill.colors <- c("BMC" = "#FF6600",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BBSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "BUMC" = "#3399FF",
                 "truth" = "black")
# load data
load("spartaOccupancyModel/reducedModelResults.RData")
RMC <- example2ReducedModelTrue

load("spartaOccupancyModel/baselineModelMCMC.RData")
BMC <- baselineModel

load("spartaOccupancyModel/updatedModelResultsBootstrap5.RData")
BUMC <- example2UpdatedModelTrue

load("spartaOccupancyModel/updatedModelResultsAuxiliary5.RData")
AUMC <- example2UpdatedModelTrue

rm(example2ReducedModelTrue, example2UpdatedModelTrue, baselineModel)

allModels <- list(RMC, BMC, BUMC, AUMC)


growthRate <- function(x){
  N <- length(x)
  ret <- ((x[N]/x[1])^(1/N) -1)*100
  return(ret)
}

## Extract psi.fs
ret <- lapply(allModels[-1], function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(parameter = "mean",
         yearID = seq(1,50))

colnames(ret[1:3]) <- c("BMC", "BUMC", "AUMC")





# growth <- lapply(ret, growthRate)
#
# x <- example2ReducedModelTrue
#   extractNames <- rownames(x$summary$all.chains)[grepl("eta", rownames(x$summary$all.chains))]
#   x$summary$all.chains[extractNames, 1]




retSD <- lapply(allModels[-1], function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})%>%
  do.call("cbind", .)%>%
as.data.frame()%>%
  mutate(parameter = "sd",
         yearID = seq(1,50))

colnames(retSD[1:3]) <- c("BMC", "BUMC", "AUMC")

allResults <- rbind(ret, retSD)
write.csv(allResults, file = "Figures/spartaModelUpdatedResults.csv", row.names = FALSE)

#extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})
