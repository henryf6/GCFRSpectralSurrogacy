## ******************************************************************** ##
## master_surrogacy_hypothesis.R
##
## Author: Henry Frye
## Date Created: 2020-6-23
##
## Purpose:
## Analysis on surrogacy data
## ******************************************************************** ##

## ******************************************************************** ##
####Setup####
## ******************************************************************** ##
rm(list = ls())

#set working directories
if(Sys.info()['user']=='henryfrye') setwd('/Users/henryfrye/Dropbox/Intellectual_Endeavours/UConn/Research/ZA_Dimensions_Data/henry_phd')


# Load Packages
library(tidyverse)
library(ggthemes) 
library(ggfortify)
library(reshape2)
library(ggpubr)
library(DescTools)
library(RColorBrewer)
library(cowplot)

#read data
SpecSurr <- read.csv('clean_data/MasterSurrogacyData/SpectralSurrogacyAug2020.csv')

#Scale and center to address mismatch between variable scales.
SpecSurr$SpectralDistance.scale <- scale(SpecSurr$SpectralDistance, scale = TRUE, center = TRUE)
SpecSurr$SpRich.scale <- scale(SpecSurr$SpRich, scale = TRUE, center = TRUE)
SpecSurr$EffShan.scale <- scale(SpecSurr$EffShan, scale = TRUE, center = TRUE)
SpecSurr$PetcheyFD.scale <- scale(SpecSurr$PetcheyFD, scale = TRUE, center = TRUE)
SpecSurr$PD.scale <- scale(SpecSurr$PD, scale = TRUE, center = TRUE)

#Fynbo and Succulent Karoo Biome only models
FynKar <- SpecSurr %>% dplyr::filter(Biome  == "Fynbos" | Biome == "Succulent Karoo")

## ******************************************************************** ##
####Spectral Diversity Comparison via Cross Validation####
## ******************************************************************** ##

source('/Users/henryfrye/Dropbox/Intellectual_Endeavours/UConn/Research/Functions/cross_validation_functions.R')

predictors <- colnames(SpecSurr)[2:8] 
#cleaning from older data version, no longer needed, but good to check
master_1 <- SpecSurr[,c(2:11, 14)]  %>% na.omit()
infs <- which(is.infinite(rowSums(master_1)) == TRUE)
master_clean <- master_1[!is.infinite(rowSums(master_1)),]

#Effective Shannon's diversity Model prediction
shan_cv <- matrix(nrow =10 ,ncol = length(predictors))

for(i in 1:length(predictors)){
  shan_cv[,i] <- CVmaster(data = master_clean, stat.method ='lm', formula = paste0('EffShan ~ ', predictors[i],sep =""),
                          seed = 123, k = 10, response = "EffShan", error.method = 'rmse')
}

colnames(shan_cv) <- predictors
shan_gcfr_comp <- melt(shan_cv, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#Species richness
sprich_cv <- matrix(nrow =10 ,ncol = length(predictors))
for(i in 1:length(predictors)){ 
  sprich_cv[,i] <- CVmaster(data = master_clean, stat.method ='lm', formula = paste0( 'SpRich~ ', predictors[i],sep =""),
                            seed = 123, k = 10, response = "SpRich", error.method = 'rmse')
}

colnames(sprich_cv) <- predictors
sprich_gcfr_comp <- melt(sprich_cv, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#functional diversity
func_cv <- matrix(nrow =10 ,ncol = length(predictors))
for(i in 1:length(predictors)){
  func_cv[,i] <- CVmaster(data = master_clean, stat.method ='lm', formula = paste0('PetcheyFD ~ ', predictors[i],sep =""),
                          seed = 123, k = 10, response = "PetcheyFD", error.method = 'rmse')
}
colnames(func_cv) <- predictors
func_gcfr_comp <- melt(func_cv, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#phylogenetic diversity
phylo_cv <- matrix(nrow =10 ,ncol = length(predictors))
for(i in 1:length(predictors)){
  phylo_cv[,i] <- CVmaster(data = master_clean, stat.method ='lm', formula = paste0('PD ~ ', predictors[i],sep =""),
                           seed = 123, k = 10, response = "PD", error.method = 'rmse')
}
colnames(phylo_cv) <- predictors
phylo_gcfr_comp <- melt(phylo_cv, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

master_cv <- rbind(shan_gcfr_comp, sprich_gcfr_comp, func_gcfr_comp,phylo_gcfr_comp)
master_cv$diversity <- c( rep("exp(Shannons) Diversity", length(predictors)*10),
                          rep("Species Rich", length(predictors)*10),
                          rep("Petchey Functional", length(predictors)*10),
                          rep("Phylogenetic", length(predictors)*10))
master_cv <- master_cv %>% 
  group_by(diversity) %>%  # group_by(Method, diversity)
  mutate(NRMSE = (RMSE - min(RMSE))/ (max(RMSE)- min(RMSE)))

cv_comparison_all <- ggplot(master_cv, aes(x = Method, y = RMSE)) + geom_boxplot(notch = TRUE) + facet_wrap(diversity~., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Greater Cape Floristic Region")

cv_comparison_all

ggsave('Figures/GCFR_all_cross_valid_metric_comp.png',cv_comparison_all)

means_all <- master_cv %>% group_by(Method, diversity) %>% summarize(mean(RMSE))
sds_all <- master_cv %>% group_by(Method, diversity) %>% summarize(sd(RMSE))
NormMeansAll<- master_cv %>% group_by(Method, diversity) %>% summarize(mean(NRMSE))
NormSDAll <- master_cv %>% group_by(Method, diversity) %>% summarize(sd(NRMSE))

## ******************************************************************** ##
###Subregion Blocked Cross-validation####
## ******************************************************************** ##
source('/Users/henryfrye/Dropbox/Intellectual_Endeavours/UConn/Research/Functions/cross_validation_functions.R')

predictors <- colnames(SpecSurr)[2:8] 
#cleaning from older data version, no longer needed, but good to check
master_1 <- SpecSurr[,c(2:11, 14,19)]  %>% na.omit()
infs <- which(is.infinite(rowSums(master_1[,1:11])) == TRUE)
master_clean <- master_1[!is.infinite(rowSums(master_1[,1:11])),]

#Effective Shannon's diversity Model prediction
ShanCVBlock <- matrix(nrow =6 ,ncol = length(predictors))

for(i in 1:length(predictors)){
  ShanCVBlock[,i] <- CVFactor(data = master_clean, stat.method ='lm', formula = paste0('EffShan ~ ', predictors[i],sep =""), response = "EffShan", dependent = "Subregion", error.method = 'rmse')
}

colnames(ShanCVBlock) <- predictors
ShanCVBlockComp <- melt(ShanCVBlock, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#Species richness
SpRichCVBlock <- matrix(nrow =6 ,ncol = length(predictors))
for(i in 1:length(predictors)){ 
  SpRichCVBlock[,i] <- CVFactor(data = master_clean, stat.method ='lm', formula = paste0( 'SpRich~ ', predictors[i],sep =""), response = "SpRich",dependent = "Subregion", error.method = 'rmse')
}

colnames(SpRichCVBlock) <- predictors
SpRichCVBlockComp <- melt(SpRichCVBlock, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#functional diversity
FuncCVBlock <- matrix(nrow =6 ,ncol = length(predictors))
for(i in 1:length(predictors)){
  FuncCVBlock[,i] <- CVFactor(data = master_clean, stat.method ='lm', formula = paste0('PetcheyFD ~ ', predictors[i],sep =""), response = "PetcheyFD", dependent = "Subregion", error.method = 'rmse')
}
colnames(FuncCVBlock) <- predictors
FuncCVBlockComp <- melt(FuncCVBlock, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

#phylogenetic diversity
PhyloCVBlock <- matrix(nrow =6 ,ncol = length(predictors))
for(i in 1:length(predictors)){
  PhyloCVBlock[,i] <- CVFactor(data = master_clean, stat.method ='lm', formula = paste0('PD ~ ', predictors[i],sep =""), response = "PD", dependent = "Subregion", error.method = 'rmse')
}
colnames(PhyloCVBlock) <- predictors
PhyloCVBlockComp <- melt(PhyloCVBlock, measure.vars = predictors, value.names = c("RMSE", "Method")) %>% dplyr::select(-Var1) %>%
  dplyr::select("Method" = Var2,"RMSE" = value)

MasterCVBlock <- rbind(ShanCVBlockComp, SpRichCVBlockComp, FuncCVBlockComp,PhyloCVBlockComp)
MasterCVBlock$diversity <- c( rep("exp(Shannons) Diversity", length(predictors)*6),
                              rep("Species Rich", length(predictors)*6),
                              rep("Petchey Functional", length(predictors)*6),
                              rep("Phylogenetic", length(predictors)*6))
MasterCVBlock <- MasterCVBlock %>% 
  group_by(diversity) %>% 
  mutate(NRMSE = (RMSE - min(RMSE))/ (max(RMSE)- min(RMSE)))

CVBlockFigure <- ggplot(MasterCVBlock, aes(x = Method, y = RMSE)) + geom_boxplot(notch = TRUE) + facet_wrap(diversity~., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Greater Cape Floristic Region")

CVBlockFigure

## ******************************************************************** ##
####All GCFR surrogacy model####
## ******************************************************************** ##
#the y-axis here is unscaled to match the CV models (we want the RMSE to be on the original y unit scale for interpretation)
sprich.mod <- lm(SpRich.scale ~ SpectralDistance.scale, data= SpecSurr)
summary(sprich.mod)


shan.mod <- lm(EffShan.scale ~ SpectralDistance.scale, data= SpecSurr)
summary(shan.mod)
cook.test <- cooks.distance(shan.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(shan.mod) #diagnostics, some issue with non-homogenous variation and normality assumption

func.mod <- lm(PetcheyFD.scale ~ SpectralDistance.scale, data= SpecSurr)
summary(func.mod)
cook.test <- cooks.distance(func.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(func.mod) # some trend in residuals and issues with normality

phylo.mod <- lm(PD.scale ~ SpectralDistance.scale, data= SpecSurr)
summary(phylo.mod)
cook.test <- cooks.distance(phylo.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(phylo.mod) #similar to others

#Graphic for all surrogacy models
ShanSurrAll <- ggplot(SpecSurr, aes(x = SpectralDistance.scale, y = EffShan.scale)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Effective Shannon\'s Diversity') +
  xlab('Spectral Diversity') +
  theme_bw() + 
  theme(legend.position = 'none', aspect.ratio = 1,
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 8)) +
  annotate("text", x = 3.5, y = -1, label = paste0("R Squared = ", round( summary(shan.mod)$r.squared,3) ), size =2.5 )

SpRichAll <- ggplot(SpecSurr, aes(x = SpectralDistance.scale, y = SpRich.scale)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Species Richness') +
  xlab('Spectral Diversity') +
  theme_bw() + 
  theme(legend.position = 'none', aspect.ratio = 1,
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  annotate("text", x = 3.5, y = -1, label = paste0("R Squared = ", round( summary(sprich.mod)$r.squared,3) ), size =2.5 )


FuncSurrAll <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PetcheyFD.scale)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Functional Diversity') +
  
  theme_bw() + 
  theme(legend.position = 'top', aspect.ratio = 1,
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8))+
  annotate("text", x = 3, y = -1, label = paste0("R Squared = ", round( summary(func.mod)$r.squared,3) ), size =2.5 )

PhyloSurrAll <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PD.scale)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Phylogenetic Diversity') +
  theme_bw() + 
  theme(legend.position = 'top', aspect.ratio = 1,
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  annotate("text", x = 3, y = -1, label = paste0("R Squared = ", round( summary(phylo.mod)$r.squared,3) ), size =2.5) 

SurrogacyProofConcept<- ggarrange(ShanSurrAll, SpRichAll, FuncSurrAll ,
                                  PhyloSurrAll,
                                  labels = c('A','B','C','D') )
ggsave('Figures/Figure3.pdf',SurrogacyProofConcept, dpi =800, width = 110, height = 110, units= "mm")


## ******************************************************************** ##
####Biome Models####
## ******************************************************************** ##

shan.biome.mod <- lm(EffShan.scale ~ SpectralDistance.scale  * Biome,   data = FynKar)
summary(shan.biome.mod)
cook.test <- cooks.distance(shan.biome.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(shan.biome.mod) #pretty good outside of normality

func.biome.mod <- lm(PetcheyFD.scale ~ SpectralDistance.scale * Biome,   data = FynKar)
summary(func.biome.mod)
cook.test <- cooks.distance(func.biome.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(func.biome.mod)


phylo.biome.mod <- lm(PD.scale ~ SpectralDistance.scale * Biome,   data = FynKar)
summary(phylo.biome.mod)
cook.test <- cooks.distance(phylo.biome.mod) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE) #no outliers
autoplot(phylo.biome.mod)

## ******************************************************************** ##
####Subregion Models####
## ******************************************************************** ##
#independent contrasts
contrasts(SpecSurr$Subregion) <- "contr.sum"

shan.mod.region <- lm(EffShan.scale ~  SpectralDistance.scale * Subregion, data = SpecSurr )
summary(shan.mod.region)
cook.test <- cooks.distance(shan.mod.region) > 0.5 #suggested rule of thumb, take with a grain of salt
which(cook.test == TRUE)
autoplot(shan.mod.region) # model diagnostics look reasonable here; slight residual trend

func.mod.region <- lm(PetcheyFD.scale ~ SpectralDistance.scale* Subregion, data = SpecSurr )
summary(func.mod.region)
cook.test <- cooks.distance(func.mod.region) > 0.5 #suggested rule of thumb from stats dept
which(cook.test == TRUE)
autoplot(func.mod.region)

phylo.mod.region <- lm(PD.scale ~ SpectralDistance.scale* Subregion, data = SpecSurr)
summary(phylo.mod.region)
cook.test <- cooks.distance(phylo.mod.region) > 0.5 #suggested rule of thumb from stats dept
which(cook.test == TRUE)
autoplot(phylo.mod.region)

## ******************************************************************** ##
####Surrogacy Main Figure####
##Fig. 4 in manuscript
## ******************************************************************** ##
shan_surr <- ggplot(SpecSurr, aes(x = SpectralDistance.scale, y = EffShan.scale, color = Subregion)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Effective Shannon\'s Diversity') +
  xlab('Spectral Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'none', aspect.ratio = 1,
        axis.title.y = element_text(size = 8.5),
        axis.title.x = element_text(size = 10))

#Omitted for sake of clarity
# sp_rich_surr <-ggplot(SpecSurr, aes(SpectralDistance.scale, y = SpRich, color = Subregion)) +
#   geom_point(shape =21) +
#   geom_smooth(method ="glm", method.args = list(family = 
#                                                   "poisson")) +
#   xlab('Spectral Diversity') +
#   ylab('Species Richness') +
#   scale_color_brewer(palette = "Dark2") +
#   theme_classic() + theme(legend.position = 'none')

func_surr <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PetcheyFD.scale, color = Subregion)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Functional Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'top', aspect.ratio = 1/1,
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))

phylo_surr <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PD.scale, color = Subregion)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Phylogenetic Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'top', aspect.ratio = 1,
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))


#biome surrogacy graphics
table(SpecSurr$Subregion, SpecSurr$Biome)

shan_surr_bio <- ggplot(FynKar, aes(x = SpectralDistance.scale, y = EffShan.scale, color = Biome)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Effective Shannon\'s Diversity') +
  xlab('Spectral Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', aspect.ratio = 1,
        axis.title.y = element_text(size = 8.5),
        axis.title.x = element_text(size = 10))

#omitted for clarity sake
# sp_rich_surr_bio <-ggplot(FynKar, aes(SpectralDistance.scale, y = SpRich, color = Biome)) +
#   geom_point(shape =21) +
#   geom_smooth(method ="glm", method.args = list(family = 
#                                                   "poisson")) +
#   xlab('Dendrogram Length Spectral Diversity') +
#   ylab('Species Richness') +
#   theme_classic() + theme(legend.position = 'none')

func_surr_bio <- ggplot(FynKar, aes(SpectralDistance.scale, y = PetcheyFD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Functional Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw()+ 
  labs(color = "Biome") +
 theme(legend.position = 'top', aspect.ratio = 1,
       axis.title.y = element_text(size = 10),
       axis.title.x = element_text(size = 10))

phylo_surr_bio <- ggplot(FynKar, aes(SpectralDistance.scale, y = PD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Phylogenetic Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', aspect.ratio = 1,
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))

one<- ggarrange(shan_surr, func_surr,phylo_surr, common.legend = TRUE, nrow = 1, ncol =3, labels = c('D','E','F') )
two<-ggarrange(shan_surr_bio, func_surr_bio, phylo_surr_bio, nrow = 1, ncol =3, common.legend = TRUE, labels= c('A','B','C'))
master_surrogacy <- ggarrange(two, one, nrow = 2)

ggsave('Figures/Figure4.pdf',master_surrogacy,  height = 168, width = 168, units= "mm", dpi = 800 )

#presentation figures
shan_surr_bio2 <- ggplot(FynKar, aes(x = SpectralDistance.scale, y = EffShan.scale, color = Biome)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Effective Shannon\'s Diversity') +
  xlab('Spectral Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', text = element_text(size= 20) ,
        aspect.ratio = 1)

func_surr_bio2 <- ggplot(FynKar, aes(SpectralDistance.scale, y = PetcheyFD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Functional Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw()+ 
  labs(color = "Biome") +
  theme(legend.position = 'top', text = element_text(size= 20) ,
        aspect.ratio = 1)

phylo_surr_bio2 <- ggplot(FynKar, aes(SpectralDistance.scale, y = PD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Phylogenetic Diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', text = element_text(size= 20) ,
        aspect.ratio = 1)

##BiomePres<- ggarrange(shan_surr_bio2, func_surr_bio2, phylo_surr_bio2, nrow = 1, ncol =3, common.legend = TRUE, labels= NA)
ggsave('Figures/SurrogcyBiomePresentationShan.png',shan_surr_bio2,  height = 15, width = 15, units= "cm", dpi = 300 )
ggsave('Figures/SurrogcyBiomePresentationFunc.png',func_surr_bio2,  height = 15, width = 15, units= "cm", dpi = 300 )
ggsave('Figures/SurrogcyBiomePresentationPhylo.png',phylo_surr_bio2,  height = 15, width = 15, units= "cm", dpi = 300 )


shan_surr2 <- ggplot(SpecSurr, aes(x = SpectralDistance.scale, y = EffShan.scale, color = Subregion)) +
  geom_point(shape = 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  ylab('Effective Shannon\'s Diversity') +
  xlab('Spectral Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'none', text = element_text(size= 20) ,aspect.ratio = 1 )

func_surr2 <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PetcheyFD.scale, color = Subregion)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Functional Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'none', text = element_text(size= 20) ,aspect.ratio = 1)

phylo_surr2 <- ggplot(SpecSurr, aes(SpectralDistance.scale, y = PD.scale, color = Subregion)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Spectral Diversity') +
  ylab('Phylogenetic Diversity') +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(color = "Subregion") +
  theme(legend.position = 'top', text = element_text(size= 20) , aspect.ratio = 1)

ggsave('Figures/SurrogcySubregPresentationShan.png',shan_surr2,  height = 15, width = 15, units= "cm", dpi = 300 )
ggsave('Figures/SurrogcySubregPresentationFunc.png',func_surr2,  height = 15, width = 15, units= "cm", dpi = 300 )
ggsave('Figures/SurrogcySubregPresentationPhylo.png',phylo_surr2,  height = 15, width = 15, units= "cm", dpi = 300 )


## ******************************************************************** ##
####Spectral Redundancy Model####
## ******************************************************************** ##

#omitted; interesting but less applicable for broader audience
# spec.red.mod <- lm(SpectralDistance ~ SpRich.scale * Subregion,   data = SpecSurr)
# summary(spec.red.mod)
# 
# func.red.mod <- lm(PetcheyFD ~ SpRich.scale * Subregion,   data = SpecSurr)
# summary(func.red.mod)
# 
# phylo.red.mod <- lm(PD ~ SpRich.scale * Subregion,   data = SpecSurr)
# summary(phylo.red.mod)

func.red.mod.bio <- lm(PetcheyFD.scale ~ SpRich.scale * Biome,   data = FynKar)
summary(func.red.mod.bio)
autoplot(func.red.mod.bio)

phylo.red.mod.bio <- lm(PD.scale ~ SpRich.scale * Biome,   data = FynKar)
summary(phylo.red.mod.bio)
autoplot(phylo.red.mod.bio)

spec.red.mod.bio <- lm(SpectralDistance.scale ~ SpRich.scale * Biome,   data = FynKar)
summary(spec.red.mod.bio)
autoplot(spec.red.mod.bio)



## ******************************************************************** ##
####Spectral Redundancy Figure 5####
## ******************************************************************** ##

#again subregion level analysis here is outside of the main audience interest
# func_spe <- ggplot(SpecSurr, aes(x = SpRich.scale, y = PetcheyFD, color = Subregion)) +
#   geom_point(shape= 21) +
#   geom_smooth(method ="lm") +
#   ylab('Functional Diversity') +
#   xlab('Species Richness') +
#   scale_color_brewer(palette = "Dark2") +
#   theme_bw() + 
#   labs(color = "Subregion") +
#   theme(legend.position = 'top')
# 
# phylo_spe <- ggplot(SpecSurr, aes(x = SpRich.scale, y = Phylo, color = Subregion)) +
#   geom_point(shape= 21) +
#   geom_smooth(method ="lm") +
#   xlab('Species Richness') +
#   ylab('Phylogenetic diversity') +
#   scale_color_brewer(palette = "Dark2") +
#   theme_bw() + 
#   labs(color = "Subregion") +
#   theme(legend.position = 'top')
# 
# spectra_spe <- ggplot(SpecSurr, aes(x = SpRich.scale, y = SpectralDistance, color = Subregion)) +
#   geom_point(shape= 21) +
#   geom_smooth(method ="lm") +
#   xlab('Species Richness') +
#   ylab('Spectral diversity') +
#   scale_color_brewer(palette = "Dark2") +
#   theme_bw() + 
#   labs(color = "Subregion") +
#   theme(legend.position = 'top')

func_red_biom <- ggplot(FynKar, aes(x = SpRich.scale, y = PetcheyFD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  ylab('Functional Diversity') +
  xlab('Species Richness') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', axis.title.y = element_text(size = 10),
        axis.text.y =  element_text(size = 6, angle = 45), aspect.ratio = 1)

phylo_red_biom <- ggplot(FynKar, aes(x = SpRich.scale, y = PD.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  xlab('Species Richness') +
  ylab('Phylogenetic diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', axis.title.y = element_text(size = 10),
        axis.text.y =  element_text(size = 6, angle = 45), aspect.ratio = 1)

spectra_red_biom <- ggplot(FynKar, aes(x = SpRich.scale, y = SpectralDistance.scale, color = Biome)) +
  geom_point(shape= 21) +
  geom_smooth(method ="lm") +
  geom_abline() +
  xlab('Species Richness') +
  ylab('Spectral diversity') +
  scale_color_brewer(palette = "Set2") +
  theme_bw() + 
  labs(color = "Biome") +
  theme(legend.position = 'top', axis.title.y = element_text(size = 10),
        axis.text.y =  element_text(size = 6, angle = 45), aspect.ratio = 1)

#onered <- ggarrange(func_spe, phylo_spe,spectra_spe, common.legend = TRUE, nrow = 1, ncol =3, labels = c('A','B','C'))
twored <-ggarrange(spectra_red_biom + rremove('xlab'), func_red_biom + rremove('xlab'), phylo_red_biom, 
                   nrow = 3, ncol =1, common.legend = TRUE, labels= "AUTO")
#twored
#master_redundancy <- ggarrange(onered, twored, nrow = 2)

ggsave('Figures/Figure5.pdf',twored, height = 230, width = 79, units= "mm", dpi =800 )

## ******************************************************************** ##
####Evenness Comparison####
## ******************************************************************** ##
ggplot(FynKar, aes(x = Biome, y  =ShannonEven)) + geom_boxplot()
