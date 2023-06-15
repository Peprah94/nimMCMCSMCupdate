# Results file sizes are large. We have attached the link to access our results
# https://www.dropbox.com/sh/9gdf426e2sw5yas/AADqlpM6jxvqRQMmccku3WxTa?dl=0


# Load data
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)


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
  ret <- sqrt(mean((x-y)^2))
  return(ret)
}

## Load data for Aux PF
load("linearGaussianSSM/auxiliaryPF/estimatesAFNew1.RData")
auxPF1 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)

load("linearGaussianSSM/auxiliaryPF/estimatesAFNew2.RData")
auxPF2 <- auxiliaryEstimates[1:15]
rm(auxiliaryEstimates)

# Put all dataset together
auxPF <- c(auxPF1, auxPF2)
rm(auxPF1, auxPF2)

## Load Bootstrap results
load("linearGaussianSSM/bootstrapPF/estimatesBFNew1.RData")
bootPF1 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)

load("linearGaussianSSM/bootstrapPF/estimatesBFNew2.RData")
bootPF2 <- bootstrapEstimates[1:15]
rm(bootstrapEstimates)

# Put all dataset together
bootPF <- c(bootPF1, bootPF2)
rm(bootPF1, bootPF2)

## Put all results together

results1 <- c(auxPF,
              bootPF)

# Flatten list
results <- purrr::flatten(results1)


## Extract model Parameters

nodeNames = c("ahat", "bhat", "chat")

#estimate the bias of model parameter
biasModelPars <- lapply(results, function(x){
  red <- x[[2]]$all.chains
  output <-  red[rownames(red)[rownames(red) %in% nodeNames],"Mean"]
  return(output)
})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(t = rep(c(rep(50, 2),
                          rep(c(49, 45, 20, 10, 5), 4)), 60),
                model = rep(
                  c(rep(c("ABSC", "ABMC",
                          rep("ARMC", 5),
                          rep("ARSC", 5),
                          rep("AUMC", 5),
                          rep("AUSC", 5)), 30),
                    rep(c("BBSC", "BBMC",
                          rep("BRMC", 5),
                          rep("BRSC", 5),
                          rep("BUMC", 5),
                          rep("BUSC", 5)), 30))
                ))%>%
  dplyr::group_by(model, t)%>%
  dplyr::summarise(ahat = mean(ahat),
                   bhat = mean(bhat),
                   chat = mean(chat))%>%
  reshape2::melt(., id.vars = c("t", "model"))

## Extract results for auxiliary PF and plot

# Fill colors for plotting
fill.colors <- c("BMC" = "#3399FF",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "ARSC" = "#E69F00",
                 "AUMC" = "#0033FF",
                 "AUSC" = "#FF6600",
                 "truth" = "black")

# extract results
biasModelPars1 <- biasModelPars%>%
  filter(model %in% c("ABMC", "ABSC",
                      "ARMC", "ARSC",
                      "AUMC", "AUSC"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

#plot results
biasModelParsAuxPF <- ggplot(data = biasModelPars1,
                             mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "red")+
  facet_wrap(~variable, ncol = 3, scales = "free_y", 
             labeller = as_labeller(c( "ahat" = "a", 
                                       "bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(expression(t[r]))+
  ylab("Bias")


## Extract bootstrap PF and plot
fill.colors1 <- c("BMC" = "#3399FF",
                  "BBSC" = "#FF3300",
                  "RMC" = "#00CCFF",
                  "BRSC" = "#E69F00",
                  "BUMC" = "#0033FF",
                  "BUSC" = "#FF6600",
                  "truth" = "black")

biasModelPars2 <- biasModelPars%>%
  filter(model %in% c("BBMC", "BBSC",
                      "BRMC", "BRSC",
                      "BUMC", "BUSC"))%>%
  dplyr::mutate(model = replace(model, model == "BBMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "BRMC", "RMC"))

biasModelParsBootPF <- ggplot(data = biasModelPars2,
                              mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "red")+
  facet_wrap(~variable, ncol = 3, scales = "free_y", 
             labeller = as_labeller(c( "ahat" = "a", 
                                       "bhat" = "b",
                                       "chat" = "c")))+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors1)+
  xlab(expression(t[r]))+
  ylab("Bias")



### Extract RMSE results and plot
# Estimating MCSE
mcseEstimates <- lapply(results, function(x){mcmcse::mcse.mat(as.matrix(rbind(x[[1]]$chain1,
                                                                              x[[1]]$chain2,
                                                                              x[[1]]$chain3))[, c("a", "b", "c")],
                                                              method = "bm",
                                                              #size = bacthSize,
                                                              g = NULL)%>%
    as.data.frame()%>%
    dplyr::select(se)%>%
    t()})%>%
  do.call("rbind", .)%>%
  as.data.frame()%>%
  dplyr::mutate(t = rep(c(rep(50, 2),
                          rep(c(49, 45, 20, 10, 5), 4)), 60),
                model = rep(
                  c(rep(c("ABSC", "ABMC",
                          rep("ARMC", 5),
                          rep("ARSC", 5),
                          rep("AUMC", 5),
                          rep("AUSC", 5)), 30),
                    rep(c("BBSC", "BBMC",
                          rep("BRMC", 5),
                          rep("BRSC", 5),
                          rep("BUMC", 5),
                          rep("BUSC", 5)), 30))
                ))%>%
  reshape2::melt(., id.vars = c("t", "model"))%>%
  dplyr::group_by(model, t, variable)%>%
  #dplyr::group_by(t)%>%
  dplyr::summarise(value = mean(value))

#select for auxPF
mcseEstimatesPars1 <- mcseEstimates%>%
  filter(model %in% c("ABMC", "ABSC",
                      "ARMC", "ARSC",
                      "AUMC", "AUSC"))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

mcseModelParsAuxPF <- ggplot(data = mcseEstimatesPars1,
                             mapping = aes(x = as.factor(t), y = value, col = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(expression(t[r]))+
  ylab("MCSE")

# Extract bootstrap PF
mcseEstimatesPars2 <- mcseEstimates%>%
  filter(model %in% c("BBMC", "BBSC",
                      "BRMC", "BRSC",
                      "BUMC", "BUSC"))%>%
  dplyr::mutate(model = replace(model, model == "BBMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "BRMC", "RMC"))

mcseModelParsBootPF <- ggplot(data = mcseEstimatesPars2,
                              mapping = aes(x = as.factor(t), y = value, color = model))+
  geom_point(aes(shape = model), position = position_dodge(width = 0.9), size = 8)+
  facet_wrap(~variable, ncol = 3, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors1)+
  xlab(expression(t[r]))+
  ylab("MCSE")


## Results for latent state distribution
# auxiliary latent states
auxLatentEst <- lapply(seq_along(auxPF), function(i){
  x <- auxPF[[i]]
  y <- simData[[i]]$x
  
  BBMC <- averageRMSE(x[[1]][[2]]$all.chains[!rownames(x[[1]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  BBSC <- averageRMSE(x[[2]][[2]]$all.chains[!rownames(x[[2]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  
  #tr= 49
  BRUMC_49 <- averageRMSE(c(x[[3]][[2]]$all.chains[!rownames(x[[3]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[13]][[2]]$all.chains[!rownames(x[[13]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_45 <- averageRMSE(c(x[[4]][[2]]$all.chains[!rownames(x[[4]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[14]][[2]]$all.chains[!rownames(x[[14]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_20 <- averageRMSE(c(x[[5]][[2]]$all.chains[!rownames(x[[5]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[15]][[2]]$all.chains[!rownames(x[[15]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_10 <- averageRMSE(c(x[[6]][[2]]$all.chains[!rownames(x[[6]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[16]][[2]]$all.chains[!rownames(x[[16]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_5 <- averageRMSE(c(x[[7]][[2]]$all.chains[!rownames(x[[7]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                           x[[17]][[2]]$all.chains[!rownames(x[[17]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  #SMC
  #tr= 49
  BRUSC_49 <- averageRMSE(c(x[[8]][[2]]$all.chains[!rownames(x[[8]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[18]][[2]]$all.chains[!rownames(x[[18]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_45 <- averageRMSE(c(x[[9]][[2]]$all.chains[!rownames(x[[9]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[19]][[2]]$all.chains[!rownames(x[[19]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_20 <- averageRMSE(c(x[[10]][[2]]$all.chains[!rownames(x[[10]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[20]][[2]]$all.chains[!rownames(x[[20]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_10 <- averageRMSE(c(x[[11]][[2]]$all.chains[!rownames(x[[11]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[21]][[2]]$all.chains[!rownames(x[[21]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_5 <- averageRMSE(c(x[[12]][[2]]$all.chains[!rownames(x[[12]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                           x[[22]][[2]]$all.chains[!rownames(x[[22]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  ret <- data.frame(BBMC, BBSC,
                    BRUMC_49, BRUMC_45,
                    BRUMC_20, BRUMC_10,
                    BRUMC_5,BRUSC_49, BRUSC_45,
                    BRUSC_20, BRUSC_10,
                    BRUSC_5
  )
  return(ret)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt()%>%
  dplyr::group_by(variable)%>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value))%>%
  dplyr::mutate( t = c(50, 50,
                       rep(c(49, 45, 20, 10, 5), 2)),
                 model = c("ABSC", "ABMC",
                           rep(c("AUMC", "AUSC"), each = 5)))%>%
  dplyr::mutate(model = replace(model, model == "ABMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "ARMC", "RMC"))

#Plot error plot
errorModelParsAuxPF <- ggplot(data = auxLatentEst, mapping = aes(x = as.factor(t), y = mean, col = model))+
  geom_point(position = position_dodge(width = 0.7), size = 8)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),size=1)+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors)+
  xlab(expression(t[r]))+
  ylab(expression(paste("Mean ", "\u00B1", " SE")))


## Bootstrap PF
bootLatentEst <-lapply(seq_along(bootPF), function(i){
  x <- bootPF[[i]]
  y <- simData[[i]]$y
  
  BBMC <- averageRMSE(x[[1]][[2]]$all.chains[!rownames(x[[1]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  BBSC <- averageRMSE(x[[2]][[2]]$all.chains[!rownames(x[[2]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"], y)
  
  #tr= 49
  BRUMC_49 <- averageRMSE(c(x[[3]][[2]]$all.chains[!rownames(x[[3]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[13]][[2]]$all.chains[!rownames(x[[13]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_45 <- averageRMSE(c(x[[4]][[2]]$all.chains[!rownames(x[[4]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[14]][[2]]$all.chains[!rownames(x[[14]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_20 <- averageRMSE(c(x[[5]][[2]]$all.chains[!rownames(x[[5]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[15]][[2]]$all.chains[!rownames(x[[15]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_10 <- averageRMSE(c(x[[6]][[2]]$all.chains[!rownames(x[[6]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[16]][[2]]$all.chains[!rownames(x[[16]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUMC_5 <- averageRMSE(c(x[[7]][[2]]$all.chains[!rownames(x[[7]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                           x[[17]][[2]]$all.chains[!rownames(x[[17]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  #SMC
  #tr= 49
  BRUSC_49 <- averageRMSE(c(x[[8]][[2]]$all.chains[!rownames(x[[8]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[18]][[2]]$all.chains[!rownames(x[[18]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_45 <- averageRMSE(c(x[[9]][[2]]$all.chains[!rownames(x[[9]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[19]][[2]]$all.chains[!rownames(x[[19]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_20 <- averageRMSE(c(x[[10]][[2]]$all.chains[!rownames(x[[10]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[20]][[2]]$all.chains[!rownames(x[[20]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_10 <- averageRMSE(c(x[[11]][[2]]$all.chains[!rownames(x[[11]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                            x[[21]][[2]]$all.chains[!rownames(x[[21]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  BRUSC_5 <- averageRMSE(c(x[[12]][[2]]$all.chains[!rownames(x[[12]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"],
                           x[[22]][[2]]$all.chains[!rownames(x[[22]][[2]]$all.chains) %in% c("a", "b", "c", "ahat", "bhat", "chat"), "Mean"][-1]), y)
  
  ret <- data.frame(BBMC, BBSC,
                    BRUMC_49, BRUMC_45,
                    BRUMC_20, BRUMC_10,
                    BRUMC_5,BRUSC_49, BRUSC_45,
                    BRUSC_20, BRUSC_10,
                    BRUSC_5
  )
  return(ret)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt()%>%
  dplyr::group_by(variable)%>%
  dplyr::summarise(mean = mean(value),
                   sd = sd(value))%>%
  dplyr::mutate( t = c(50, 50,
                       rep(c(49, 45, 20, 10, 5), 2)),
                 model = c("BBSC", "BBMC",
                           rep(c("BUMC", "BUSC"), each = 5)))%>%
  dplyr::mutate(model = replace(model, model == "BBMC", "BMC"))%>%
  dplyr::mutate(model = replace(model, model == "BRMC", "RMC"))

#Plot error plot
errorModelParsBootPF <- ggplot(data = bootLatentEst, mapping = aes(x = as.factor(t), y = mean, col = model))+
  geom_point(position = position_dodge(width = 0.7), size = 8)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5, position = position_dodge(width = 0.7),size=1)+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  scale_color_manual(values = fill.colors1)+
  xlab(expression(t[r]))+
  ylab(expression(paste("Mean ", "\u00B1", " SE")))


## Put plots together
auxPFpLots <- ggpubr::ggarrange(biasModelParsAuxPF,
                                mcseModelParsAuxPF ,
                                errorModelParsAuxPF,
                                ncol = 1,
                                common.legend = TRUE,
                                labels = c("A)", "B)", "C)"),
                                font.label = list(size = 25))
ggsave(filename = "Figures/auxPFplots.png",
       plot = auxPFpLots,
       width = 22,
       height = 20,
       dpi = 150)

bootPFpLots <- ggpubr::ggarrange(biasModelParsBootPF,
                                 mcseModelParsBootPF ,
                                 errorModelParsBootPF,
                                 ncol = 1,
                                 common.legend = TRUE,
                                 labels = c("A)", "B)", "C)"),
                                 font.label = list(size = 25))

ggsave(filename = "Figures/bootPFpLots.png",
       plot = bootPFpLots,
       width = 22,
       height = 20,
       dpi = 150)



### Convergence plots

## Auxiliary PF
#MCMC
auxPFSubset <- auxPF[[1]][c(2, 4, 6, 14, 16)]

aPlots <- lapply(auxPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "b", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20))
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"),
                         font.label = list(size = 25))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90, size = 25),
                          bottom = text_grob("Iterations", size = 25))

ggsave(filename = "Figures/auxPFMCMCconvergence.png",
       plot = aPlots,
       width = 22,
       height = 20,
       dpi = 150)

#SMC
auxPFSubset <- auxPF[[1]][c(1, 9, 11, 19, 21)]

aPlots <- lapply(auxPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "b", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20))
}
)

fig <- ggpubr::ggarrange(plotlist = aPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"),
                         font.label = list(size = 25))
aPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90, size = 25),
                          bottom = text_grob("Iterations", size = 25))

ggsave(filename = "Figures/auxPFSMCconvergence.png",
       plot = aPlots,
       width = 22,
       height = 20,
       dpi = 150)


## bootstrap PF

#MCMC
bootPFSubset <- bootPF[[1]][c(2, 4, 6, 14, 16)]


bPlots <- lapply(bootPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "b", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20))
}
)

fig <- ggpubr::ggarrange(plotlist = bPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"),
                         font.label = list(size = 25))
bPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90, size = 25),
                          bottom = text_grob("Iterations", size = 25))

ggsave(filename = "Figures/bootPFMCMCconvergence.png",
       plot = bPlots,
       width = 22,
       height = 20,
       dpi = 150)

#SMC
bootPFSubset <- bootPF[[1]][c(1, 9, 11, 19, 21)]

bPlots <- lapply(bootPFSubset, function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% c("a", "b", "c", "x[1]", "x[2]"))%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20))
}
)

fig <- ggpubr::ggarrange(plotlist = bPlots,
                         legend = "bottom",
                         common.legend = TRUE,
                         nrow = 1,
                         labels = c("A", "B", "C", "D", "E"),
                         font.label = list(size = 25))
bPlots <- annotate_figure(fig,
                          left = text_grob("Value", rot = 90, size = 25),
                          bottom = text_grob("Iterations", size = 25))

ggsave(filename = "Figures/bootPFSMCconvergence.png",
       plot = bPlots,
       width = 22,
       height = 20,
       dpi = 150)

rm(list=ls())




######################
# Simulation study 2: Dynamic occupancy model
###################

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

fill.colors2 <- c("ABSC" = "#3399FF",
                  "BBSC" = "#FF3300",
                  "RMC" = "#FF6600",
                  "BMC" = "blue",
                  "AUMC" = "#00CCFF",
                  "BUMC" = "#E69F00",
                  "truth" = "black")

load("~/Dropbox/Data for PHD/particleFilters/dynamicOccupanyModel/Bootstrap/simDataDynamicOccupancy.RData")

#load bootstrap results
load("dynamicOccupanyModel/Bootstrap/example4BaselineMCMC.RData")
BBMC <- baselineMCMC

load("dynamicOccupanyModel/Bootstrap/example4BaselineSMCboostrap.RData")
BBSC <- baselineModel

load("dynamicOccupanyModel/Bootstrap/example4ReducedBootstrapTrue.RData")
BRMC <- example2ReducedModelTrue

load("dynamicOccupanyModel/Bootstrap/example4UpdatedBootstrapTrue.RData")
BUMC <- example2UpdatedModelTrue

load("dynamicOccupanyModel/auxiliaryPF/example2UpdatedAuxiliaryTrue.RData")
AUMC <- example2UpdatedModelTrue

load("dynamicOccupanyModel/Bootstrap/example4UpdatedBootstrapTrueM1000.RData")
BUMC1000 <- example2UpdatedModelTrue

load("dynamicOccupanyModel/auxiliaryPF/example2UpdatedAuxiliaryTrueM1000.RData")
AUMC1000 <- example2UpdatedModelTrue

load("dynamicOccupanyModel/auxiliaryPF/example2BaselineSMCAuxiliary.RData")
ABSC <- baselineModel

#remove individual data results
rm(baselineModel, baselineMCMC, example2ReducedModelTrue, example2UpdatedModelTrue)

allModels <- list(BBSC, BBMC,
                  BRMC, BUMC,
                  ABSC, AUMC)

allModels1000 <- list(BBSC, BBMC,
                      BRMC, BUMC1000,
                      ABSC, AUMC1000)


## Extract realised occupancy and plot
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Mean"]
})

ret1000 <- lapply(allModels1000, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Mean"]
})

occSites <- data.frame(BBSC = ret[[1]],
                       BMC = ret[[2]],
                       BUMC = c(ret[[3]], ret[[4]][-1]),
                       ABSC = ret[[5]],
                       AUMC = c(ret[[3]], ret[[6]][-1]),
                       truth = simData$occSites,
                       year = 1:55)%>%
  reshape2::melt(id.vars = c("year"))%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = value, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 1),size=8)+
  geom_line(position = position_dodge(width = 1),linewidth = 2)+
  theme_classic()+
  #ylim(c(0.75, 1.6))+
  geom_vline(xintercept = 45, linetype = "dashed", col = "red")+
  scale_color_manual(name = "Model", values = fill.colors2)+
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50))+
  xlab("")+
  ylab("")+
  labs(title = "A) M = 100")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=30),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))

occSites1000 <- data.frame(BBSC = ret1000[[1]],
                           BMC = ret1000[[2]],
                           BUMC = c(ret1000[[3]], ret1000[[4]][-1]),
                           ABSC = ret1000[[5]],
                           AUMC = c(ret1000[[3]], ret1000[[6]][-1]),
                           truth = simData$occSites,
                           year = 1:55)%>%
  reshape2::melt(id.vars = c("year"))%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = value, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 1),size=8)+
  geom_line(position = position_dodge(width = 1),linewidth = 2)+
  theme_classic()+
  #ylim(c(0.75, 1.6))+
  geom_vline(xintercept = 45, linetype = "dashed", col = "red")+
  scale_color_manual(name = "Model", values = fill.colors2)+
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50))+
  xlab("")+
  ylab("")+
  labs(title = "B) M = 1000")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))

fig <- ggpubr::ggarrange(occSites, occSites1000,
                         nrow = 2, ncol = 1,
                         common.legend = TRUE,
                         legend = "top",
                         font.label = list(size = 25))

occPlots <- annotate_figure(fig,
                            left = text_grob("Average occupied sites", rot = 90, size = 30),
                            bottom = text_grob("Year", size = 30))


ggsave(filename = "Figures/psifs.png",
       plot = occPlots,
       width = 22,
       height = 16,
       dpi = 100)


## Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})

## Model Parameters and plot them
pars <- c('alphaPSig', 'betaPSig',
          'alphaPsiSig', 'betaPsiSig',
          'alphaPhi', 'betaPhi',
          'alphaP', 'betaP',
          'alphaPsi', 'betaPsi',
          "colonisationProb")

# Mean of pars
ret1 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, "Mean"]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(truth = c(1.6, 2,
                   2,
                   1.11, 4,
                   -0.493, 3,
                   1.5,
                   -1.69, 2,
                   0.05
  ))
colnames(ret1)[1:6] <- c('BBSC', 'BMC',
                         'RMC', 'BUMC',
                         'ABSC', 'AUMC')


# sd of pars
ret2 <- lapply(allModels, function(x){
  x[[2]]$all.chains[rownames(x[[2]]$all.chains) %in% pars, 3]
})%>%
  do.call("cbind", .)%>%
  as.data.frame()
colnames(ret2)[1:6] <- c('BBSC', 'BMC',
                         'RMC', 'BUMC',
                         'ABSC', 'AUMC')


## Effective sample size
ESSret <- lapply(allModels, function(x) {
  ret <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(ESSret)[1:6] <- c('BBSC', 'BMC',
                           'RMC', 'BUMC',
                           'ABSC', 'AUMC')

ESSret1 <- ESSret%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]), 
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]), 
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi), 
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p), 
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi), 
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors2)+
  xlab("")+
  ylab("ESS")+
  labs(title = "A) M = 100")+
  labs(col = "Model", shape = "Model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))



ESSret1000 <- lapply(allModels1000, function(x) {
  ret <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(ESSret1000)[1:6] <- c('BBSC', 'BMC',
                               'RMC', 'BUMC',
                               'ABSC', 'AUMC')

ESSret2 <- ESSret1000%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]), 
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]), 
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi), 
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p), 
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi), 
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors2)+
  xlab("")+
  ylab(" ")+
  labs(title = "B) M = 1000")+
  labs(col = "Model", shape = "Model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))


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
  
  ret <- ess/timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(EfficiencyRet)[1:6] <- c('BBSC', 'BMC',
                                  'RMC', 'BUMC',
                                  'ABSC', 'AUMC')

EfficiencyRet1 <- EfficiencyRet%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]), 
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]), 
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi), 
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p), 
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi), 
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors2)+
  xlab("")+
  ylab("Efficiency")+
  labs(col = "Model", shape = "Model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))



EfficiencyRet1000 <- lapply(allModels1000, function(x) {
  nDim <- length(x$timeRun)
  
  ess <- mcmcse::ess((as.matrix(x$samples))[,pars])%>%
    as.data.frame()
  timesRet <- if(nDim == 1){
    as.numeric(x$timeRun, units = "secs")
  }else{
    as.numeric(x$timeRun$all.chains, units = "secs")
  }
  
  ret <- ess/timesRet
})%>%
  do.call("cbind",.)%>%
  as.matrix(., rownames.force = FALSE)%>%
  as.data.frame(., row.names = NULL)

colnames(EfficiencyRet1000)[1:6] <- c('BBSC', 'BMC',
                                      'RMC', 'BUMC',
                                      'ABSC', 'AUMC')

EfficiencyRet2 <- EfficiencyRet1000%>%
  dplyr::mutate(parameters = c(pars))%>%
  reshape2::melt(., id.vars = c("parameters"))%>%
  ggplot(., aes(x = parameters, y = value, col = variable))+
  geom_point(aes(shape = variable), size = 8,position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete(labels = c('alphaPSig' = expression(sigma[alpha^p]), 
                              'betaPSig'= expression(sigma[beta^p]),
                              'alphaPsiSig'= expression(sigma[alpha^psi]), 
                              'betaPsiSig'= expression(sigma[beta^psi]),
                              'alphaPhi'= expression(alpha^phi), 
                              'betaPhi'= expression(beta^phi),
                              'alphaP'= expression(alpha^p), 
                              'betaP'= expression(beta^p),
                              'alphaPsi'= expression(alpha^Psi), 
                              'betaPsi'= expression(beta^Psi),
                              "colonisationProb"= expression(gamma)))+
  scale_color_manual(values = fill.colors2)+
  xlab("")+
  ylab("")+
  labs(col = "Model", shape = "model")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        plot.title = element_text(size = 35))


## Put all plots together
fig <- ggpubr::ggarrange(ESSret1, ESSret2,
                         EfficiencyRet1, EfficiencyRet2,
                         ncol = 2, nrow = 2,
                         common.legend = TRUE,
                         legend = "top")

effESSPlot <- annotate_figure(fig,
                              bottom = text_grob("Model Parameters", size = 25))


ggsave(filename = "Figures/essEffPlotEx2.png",
       plot = effESSPlot,
       width = 18,
       height = 10,
       dpi = 100)

# Convergence of parameters
# Auxiliary PF
#auxPFSubset <- auxPF[[1]][c(2, 4, 6, 14, 16)]
aPlots <- lapply(allModels[c(2,3,4,6)], function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[1:4])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=25),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 35))
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
       width = 22,
       height = 20,
       dpi = 150)


aPlots <- lapply(allModels[c(2,3,4,6)], function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[5:8])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=25),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 35))
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
ggsave(filename = "Figures/convergenceSecExB.png",
       plot = aPlots,
       width = 22,
       height = 20,
       dpi = 150)

aPlots <- lapply(allModels[c(2,3,4,6)], function(y){
  ggmcmc::ggs(y[[1]])%>%
    filter(Parameter %in% pars[9:11])%>%
    ggmcmc::ggs_traceplot()+
    theme_classic()+
    ylab("")+
    xlab("")+
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.title = element_text(size=25),
          legend.text = element_text(size=20),
          plot.title = element_text(size = 35))
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
       width = 22,
       height = 20,
       dpi = 150)


rm(list=ls())

###########################################
# Example 2
##########################################

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

fill.colors2 <- c("ABSC" = "#3399FF",
                   "BBSC" = "#FF3300",
                   "RMC" = "#FF6600",
                   "BMC" = "blue",
                   "AUMC" = "#00CCFF",
                   "BUMC" = "#E69F00",
                   "truth" = "black")

#bootstrap
load("demographicSSM/bootstrapPF/baselineModelSMCResults.RData")
BBSC <- baselineModelSMC
rm(baselineModelSMC)

load("demographicSSM/bootstrapPF/baselineModelMCMCResults.RData")
BBMC <- baselineModelMCMC
rm(baselineModelMCMC)

load("demographicSSM/bootstrapPF/reducedModelResults.RData")
BRMC <- example2ReducedModelTrue
rm(example2ReducedModelTrue)

load("demographicSSM/bootstrapPF/updatedModelResults.RData")
BUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)

load("demographicSSM/auxilaryPF/baselineModelSMCResults.RData")
ABSC <- baselineModelSMC
rm(baselineModelSMC)

load("demographicSSM/auxilaryPF/updatedModelResults.RData")
AUMC <- example2UpdatedModelTrue
rm(example2UpdatedModelTrue)


# Put models together
allModels <- list(BBMC,
                  BRMC, 
                  BUMC, 
                  AUMC)

# Extract time
ret <- lapply(allModels, function(x){
  x[[3]]
})


# Extra growth rate and its standard deviation
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})


retSD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("gammaX", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})


# Exteach population index and its standard deviation
ret1 <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})

ret1SD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("popindex", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})


# Plot population index
popnIndex <- data.frame(#BBSC = ret1[[1]],
  BMC = ret1[[1]],
  BUMC = c(ret1[[2]], ret1[[3]][-1]),
  #ABSC = ret1[[5]],
  AUMC = c(ret1[[2]], ret1[[4]][-1]),
  year = 1999:2016)%>%
  reshape2::melt(id.vars = c("year"),
                 value.name = "mean")

popnIndexSD <- data.frame(#BBSC = ret1[[1]],
  BMC = ret1SD[[1]],
  BUMC = c(ret1SD[[2]], ret1SD[[3]][-1]),
  #ABSC = ret1[[5]],
  AUMC = c(ret1SD[[2]], ret1SD[[4]][-1]),
  year = 1999:2016)%>%
  reshape2::melt(id.vars = c("year"))

popnIndex <- cbind(popnIndex, sd = popnIndexSD[,3])%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = mean, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 0.7),size=4)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  geom_vline(xintercept = 2011, linetype = "dashed", col = "red")+
  theme_classic()+
  xlab("Year")+
  ylab("Population Size")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))+
  scale_color_manual(values = fill.colors2)

growthRate <- data.frame(#BBSC = ret[[1]],
  BMC = ret[[1]],
  BUMC = c(ret[[2]], ret[[3]]),
  #ABSC = ret[[5]],
  AUMC = c(ret[[2]], ret[[4]]),
  year = 1999:2015)%>%
  reshape2::melt(id.vars = c("year"),
                 value.name = "mean")

growthRateSD <- data.frame(#BBSC = ret[[1]],
  BMC = retSD[[1]],
  BUMC = c(retSD[[2]], retSD[[3]]),
  #ABSC = ret[[5]],
  AUMC = c(retSD[[2]], retSD[[4]]),
  year = 1999:2015)%>%
  reshape2::melt(id.vars = c("year"))

growthRate <- cbind(growthRate, sd = growthRateSD[,3])%>%
  ggplot(data = ., mapping = aes(x = year,
                                 y = mean, col = variable, group = variable))+
  geom_point(position = position_dodge(width = 0.7),size=4)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  theme_classic()+
  ylim(c(0.7, 1.6))+
  geom_hline(yintercept = 1)+
  scale_color_manual(values = fill.colors2)+
  xlab("Year")+
  ylab("Growth rate")+
  geom_vline(xintercept = 2011, linetype = "dashed", col = "red")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

fig <- ggarrange(popnIndex,
                 growthRate,
                 ncol = 1,
                 nrow= 2,
                 common.legend = TRUE,
                 legend = "top",
                 labels = c("A)", "B)"),
                 font.label = list(size = 25))

ggsave(filename = "Figures/example2.png",
       plot = fig,
       width = 18,
       height = 10,
       dpi = 100)


lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("sd.gam", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 1]
})

rm(list = ls())


#############
# Example 3
###############

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(mcmcse)
library(coda)

fill.colors <- c("BMC" = "#3399FF",
                 "ABSC" = "#FF3300",
                 "RMC" = "#00CCFF",
                 "BRSC" = "#E69F00",
                 "BUSC" = "#0033FF",
                 "BUMC" = "#FF6600",
                 "truth" = "black")

# load data
load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/reducedModelResults.RData")
BRM <- example2ReducedModelTrue

load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/updatedModelResults.RData")
BURC <- example2UpdatedModelTrue

rm(example2ReducedModelTrue, example2UpdatedModelTrue)

allModels <- list(BRM, BURC)


## Extract psi.fs
ret <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, "Median"]
})%>%
  do.call("c", .)

ret <- ret[-50]


retSD <- lapply(allModels, function(x){
  extractNames <- rownames(x[[2]]$all.chains)[grepl("psi.fs", rownames(x[[2]]$all.chains))]
  x[[2]]$all.chains[extractNames, 3]
})%>%
  do.call("c", .)

retSD <- retSD[-50]

#load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/Example3/Lasius niger_z.rdata")
#psi.fs <- out$BUGSoutput$mean$psi.fs
load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/spartaOccupancyModel/baselineModel.RData")

names <- rownames(baselineModel[[2]]$all.chains)[grepl("psi.fs", rownames(baselineModel[[2]]$all.chains))]
psi.fs <- baselineModel[[2]]$all.chains[names, "Mean"]

psi.fsSD <- baselineModel[[2]]$all.chains[names, 3]

# Extract a

retA <- lapply(seq_along(allModels), function(x){
  if(x == 1){
    vals <- allModels[[x]]$summary$all.chains[1:49, 1]
  }else{
    vals <- allModels[[x]]$summary$all.chains[2:6, 1]
  }
})%>%
  do.call("c", .)

retAsd <- lapply(seq_along(allModels), function(x){
  if(x == 1){
    vals <- allModels[[x]]$summary$all.chains[1:49, 3]
  }else{
    vals <- allModels[[x]]$summary$all.chains[2:6, 3]
  }
})%>%
  do.call("c", .)

aOut <- baselineModel$summary$all.chains[1:54, 1]

aOutsd <- baselineModel$summary$all.chains[1:54, 1]
#aOut <- out$BUGSoutput$mean$a

year <- 1970:2023

extractedValuesA <- data.frame(a = c(retA, aOut),
                               psi.fs = c(ret, psi.fs),
                               #aSD = c(retAsd, aOutsd),
                               #psi.fsSD = c(retSD, psi.fsSD),
                               Model = rep(c("BUMC", "BMC"), each = 54),
                               year = rep(year, 2))%>%
  reshape2::melt(id.vars = c("Model", "year"),
                 value.name = "mean")


extractedValuesB <- data.frame(#a = c(retA, aOut),
  #psi.fs = c(ret, psi.fs),
  aSD = c(retAsd, aOutsd),
  psi.fsSD = c(retSD, psi.fsSD),
  Model = rep(c("BUMC", "BMC"), each = 54),
  year = rep(year, 2))%>%
  reshape2::melt(id.vars = c("Model", "year"),
                 value.name = "sd")

allData <- cbind(extractedValuesA, sd = extractedValuesB[,4 ])

extractedValues <- allData %>%
  ggplot(., mapping = aes(x = year, y = mean, col = Model))+
  geom_point(position = position_dodge(width = 0.7), size = 8)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.7),linewidth = 1)+
  geom_line(position = position_dodge(width = 0.7),linewidth = 1)+
  theme_classic()+
  ylab("Estimated value")+
  xlab("Year")+
  geom_vline(xintercept = 2018, linetype = "dashed")+
  facet_wrap( ~variable, ncol = 1, nrow = 2, scales = "free_y")+
  theme(legend.position = "top")+
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))+
  scale_color_manual(values = fill.colors)

ggsave(filename = "Figures/spartaPsifs.png",
       plot = extractedValues,
       width = 18,
       height = 10,
       dpi = 100)


#plot visits
load("/Users/kwakupa/Dropbox/Data for PHD/particleFilters/Example3/Lasius niger_z.rdata")
simData <- out$model$data()
hist(simData$Year)

str(simData)

ret <- lapply(allModels, function(x){
  x[[3]]
})
