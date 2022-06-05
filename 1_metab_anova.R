

############################################################
####### SCRIPT DESCRIPTION AND NOTES
############################################################

# Travis Hattery 2022
# Script uses data which has already been log10-transformed



############################################################
####### CONDUCT ANOVAS OF 2016&2017 3-ENVIRONMENT STUDY
############################################################


# Set working directory
setwd("C:/Users/travi/Desktop/Genetics/TJH GxE Manuscript")

input <- read.csv("2_widiv_3env_plotavg_log10_wthr50.csv")

traits <- colnames(input[11:118]) #9:116 - 108 total


#Create dataframe for ANOVA results
aresults <- data.frame(matrix(ncol = 541, nrow = 6))
labels <- read.csv("labels.csv")
terms <- c("growthenv", "genotype", "growthenv/rep", "genotype:growthenv",
           "growthenv/rep/block", "residual")


colnames(aresults) <- colnames(labels)
aresults[,1] <- terms

#Create dataframe for Heritability results
hresults <- data.frame(matrix(ncol=108, nrow=1))
colnames(hresults) <- traits

k=1

for (i in traits){
  
  #ANOVA full model on one metabolite trait at a time
  temp.aov <- aov(as.formula(paste(i,"~ growthenv+growthenv/rep+genotype+genotype:growthenv+growthenv/rep/block")), data=input)
  
  
  for (j in 2:6){
    
    #store ANOVA results to data frame
    aresults[,(j+(5*(k-1)))] <- summary(temp.aov)[[1]][, j-1]
  
    #store relevant MS values to placeholder values
    n_env <- 3
    n_rep <- 2
    MSG = summary(temp.aov)[[1]][,3][2] #Mean Square Genotype extracted from ANOVA
    MSGE = summary(temp.aov)[[1]][,3][4] #Mean Square Genotype*Environment extracted from ANOVA
    MSR = summary(temp.aov)[[1]][,3][6] #Mean Square Residual extracted from ANOVA
    
    H = MSG / (MSG + (MSGE/n_env) + (MSR/(n_env*n_rep)))
    
    hresults[i] <- H
    
  }
  
  k <- k+1
  }

write.csv(aresults,"widiv_aresults.csv", row.names = FALSE)
write.csv(hresults,"widiv_hresults.csv", row.names = FALSE)



############################################################
####### CONDUCT ANOVAS OF 2020 MULTI-ENVIRONMENT STUDY
############################################################

# Prepare Environment
library(car)
setwd("C:/Users/travi/Desktop/Genetics/WiDiv Full Study")
input <- read.csv("widiv_20_rowavg.csv")
traits <- colnames(input[11:117])

input$genotype <- as.factor(input$genotype)
input$location <- as.factor(input$location)
input$planting <- as.factor(input$planting)
input$environment <- as.factor(input$environment)
input$geno_loc_plant <- as.factor(input$geno_loc_plant)
input$row <- as.factor(input$row)
input$rep <- as.factor(input$rep)
input$block <- as.factor(input$block)
input$extract <- as.factor(input$extract)
input$batch <- as.factor(input$batch)

#Create dataframe for ANOVA results
#aresults <- data.frame(matrix(ncol = 536, nrow = 6))
aresults <- data.frame(matrix(ncol = 536, nrow = 8))
labels <- read.csv("anovalabels.csv")
terms <- c("genotype", "location", "planting", "genotype:location", "location/planting", "planting/rep", "genotype:location/planting","residual")


colnames(aresults) <- colnames(labels)
aresults[,1] <- terms

#Create dataframe for Heritability results
hresults <- data.frame(matrix(ncol=107, nrow=1))
colnames(hresults) <- traits

k <- as.integer(0)

for (i in traits){
  
  #ANOVA full model on one metabolite trait at a time
  #temp.aov <- aov(as.formula(paste(i,"~ genotype + location + planting + genotype:location + location/planting + planting/rep + genotype*location/planting")), data=input)
  
  temp1.aov <- glm(as.formula(paste(i,"~ genotype + location + planting + genotype:location + location/planting + planting/rep + genotype*location/planting")), data=input)
  temp.aov <- car::Anova(temp1.aov, test.statistic="F", type=2)
  
  #store ANOVA results to data frame
  dfcol = 2+(k*5)
  sscol = 3+(k*5)
  mscol = 4+(k*5)
  fcol = 5+(k*5)
  pcol = 6+(k*5)
  
  aresults[,dfcol]<-temp.aov[,2] #df
  aresults[,sscol]<-temp.aov[,1] #ss
  aresults[,fcol]<-temp.aov[,3] #F
  aresults[,pcol]<-temp.aov[,4] #p
  aresults[,mscol]<-aresults[,sscol]/aresults[,dfcol] #ms = ss/df
  
  
  #store relevant MS values to placeholder values
  n_env <- 12
  n_rep <- 2
  MSG = aresults[1,mscol] #Mean Square Genotype extracted from ANOVA
  MSGE = aresults[7,mscol] #Mean Square Genotype*Environment extracted from ANOVA
  MSR = aresults[8,mscol] #Mean Square Residual extracted from ANOVA
  H = MSG / (MSG + (MSGE/n_env) + (MSR/(n_env*n_rep)))
  hresults[k+1] <- H
  
  k <- k+1
  
}

write.csv(aresults,"anovaresults.csv", row.names = FALSE)
write.csv(hresults,"h2results.csv", row.names = FALSE)



