## ******************************************************************** ##
## spetral_div_funcs.R
##
## Author: Henry Frye
## Date Created: 2019-1-18
##
## Purpose:
## Custom functions to 
## ******************************************************************** ##

library(geometry)
library(vegan)
#library(Momocs) #for hull volumes in ordination space

#Spectral Coefficient of Variation
spectral_diversityCV <- function(spe, spectra){
  #where spe is site by species abundance matrix
  #spectra is a species by spectra matrix
  spe_div <- matrix(nrow = dim(spe)[1], ncol = 1)
  
  for(i in 1:dim(spe)[1]){
    test <- t(spe[i,]) * spectra
    test2 <- as.matrix(spectra[-which(test == 0),])
    mean.spec <- t(as.matrix(apply(test2,2, mean)))
    sd.spec <- apply(test2,2, sd)
    spe_div[i] <- mean((sd.spec/mean.spec)*100)
    
  }
  return(spe_div)
}

#radian to degree
rad2deg <- function(rad) {(rad * 180) / (pi)}


#based on Kruse et al., 1993
spectral_angle <- function(spec, ref){
  rad2deg(acos( (sum(spec*ref)/
                   (sqrt(sum(spec*spec))*sqrt(sum(ref*ref))))
                
  ))
}

#Spectral angle diversity
#this takes the mean spectral angles for each plot, this wasn't clear in the gholizadeh paper
SAM_diversity <- function(spe, spectra){
  #where spe is site by species abundance matrix
  #spectra is a species by spectra matrix
  sam_div <- matrix(nrow = dim(spe)[1], ncol = 1)
  
  for(i in 1:dim(spe)[1]){
    site <-t(spe[i,])
    site_abun<- site[which(site > 0),]
    
    test <- ifelse( site > 0,1,0)* spectra
    plot_i <- test[which(test[,1] > 0),] * site_abun #This seems to be the most consistent way to abundance weight spectra
    mean_plot_spec <- apply(plot_i,2, mean) 
    sam_div[i]<- mean(apply(plot_i, 1, spectral_angle, ref = mean_plot_spec)) 
  }
  return(sam_div)
}

#Convex hull area inspired from Gholizadeh et al., 2018
CVAspectra <- function(spe, spectra, nwaves = 500){
  
  spe_cva <- matrix(nrow = nrow(spe), ncol = 1) 
  
  for(i in 1:nrow(spe)){
    test <- t(spe[i,]) * spectra
    test2 <- as.matrix(spectra[-which(test == 0),])
    mean.spec <- t(as.matrix(apply(test2,2, mean)))
    
    hullarea <- vector(length= nrow(test2))
    for(j in 1:nrow(test2)){
      ps <- matrix(test2[1,], mean.spec, ncol = 2, nrow = nwaves)
      hullarea[j] <- convhulln(ps, "FA")$area
    }
    
    
    spe_cva[i] <- sum(hullarea)/nrow(test2)
    
  }
  return(spe_cva)
}

#Average Spectral Information Divergence per plot
SIDspectra <- function(spe, spectra){
  
  sid <- matrix(nrow = nrow(spe), ncol = 1) 
  
  for(i in 1:nrow(spe)){
    
    test <- t(spe[i,]) * spectra
    test2 <- as.matrix(spectra[-which(test == 0),])
    mean.spec <- t(as.matrix(apply(test2,2, mean)))
    
    entropy_spp <- vector(length= nrow(test2)) #vector contianing entropy for species j in site i
    entropy_mean <- vector(length= nrow(test2)) # vector containing entropy for mean site
    
    for(j in 1:nrow(test2)){
      
      p = test2[j,] / sum(test2[j,])
      q = mean.spec/sum(mean.spec)
      
      entropy_spp[j] <- sum(p *log(p/q))
      entropy_spp[j] <- sum(q *log(q/p))
    }
    
    sid[i] <- mean(entropy_spp + entropy_mean)
  }
  
  return(sid)
  
}

#Convex hull volumes inspired byDahlin 2016
#So important note here is that you need the entire trait/spectra matrix to subset the PCA volumes
# otherwise you compared apples to oranges and get funky results 

CVVspectra <-function(spe, spectra){ #don't think I can abundance weight this...
  
  spe_cvv <- matrix(nrow = nrow(spe), ncol = 1)
  #div$sp_rich <- rowSums(spe_trimmed_spp > 0)
  #step1 is calculate a pca
  spectra.pca <- prcomp(spectra,scale = TRUE, center=TRUE)
  data.pca <- as.matrix(spectra.pca$x[,1:3])
  
  for(i in 1:nrow(spe)){
   # if(sum(spe[i,] > 0) > 3){ #have to use plots with at least 4 species to pull out 3 axes
    test <- t(spe[i,]) 
    test <- ifelse(test ==  0,0,1) #this turns it into pres/absence for now
    test2 <- as.matrix(test[-which(test == 0),]) #
    site_spp <- which(rownames(data.pca) %in% rownames(test2))
    hull_subset <- data.pca[site_spp,]
    spe_cvv[i] <- convhulln(hull_subset, "FA")$vol

  } 
  return(spe_cvv)
}



