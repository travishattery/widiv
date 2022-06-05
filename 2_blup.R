## Original script - Travis Hattery
## 12 May 2022

## Conduct BLUP on metabolite data and save to .csv files 




############################################################
####### CLEAN ENVIRONMENT & PREPARE DATASETS ########
############################################################

# Clear global environment
rm(list=ls(all=TRUE))
gc()
library(lme4)

# Find dataset
# Dataset should have:
#  Column 1: genotype
#  Column 2: env
#  Column 3: rep
#  Column 4: block
#  Column 5 to n: Data
setwd("C:/Users/travi/Desktop/Genetics/TJH GxE Manuscript/2_Heatmaps/3env")
data <- read.csv("3env.csv", header=T, stringsAsFactors=F)


colnames(data)[1] <- "genotype"

# Setting variables as factors for genotype, env, rep, block
data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])
data[,4] <- as.factor(data[,4])

# Make sure identifying data is "Factor" and dataset is "num"
str(data)

# Save numerical data column names
traits <- colnames(data)[5:ncol(data)]

############################################################
####### BLUPS FOR EVERY TRAIT ########
############################################################

# just run the entire for-loop. It will save a bunch of stuff to the
# working directory, some analysis for each trait. One of those outputs
# is a file called "blups_trait.csv"

# In the middle of this for-loop, I put some space around a part where the
# actual BLUPping happens. It's a random-effects model, so I have a few options
# and only one is active. Check out the terms in your dataset, and decide
# which of the options is best, then use that one, and comment out the others.

# Each individual trait (column 5 to n in your dataset) will have its own
# BLUP values output. It shouldn't take too long to run this script on
# your own laptop - no supercomputing necessary.

for(trait in traits){
  pdf(paste0("plots_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(data[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across environments
  stripchart(data[,trait] ~data$env,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Environment",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(data[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality #SAMPLE SIZE MUST BE BELOW 5000
  normality <- shapiro.test(data[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA (switched : in Env:Rep and Rep:Block for nesting) #ALTERED FROM ORIGINAL
  model <- lm(get(trait) ~ genotype + rep + rep/block, data=data)
  #model <- lm(get(trait) ~ genotype + env + env/rep + rep/block + genotype:env, data=data)
  

  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  
  
  
  
  # Run a random effects model
  # THIS IS WHERE BLUPS ARE CALCULATED - CHOOSE ONLY ONE TO BE UNCOMMENTED OUT
  
  #model.1 <- lmer(get(trait) ~ (1|genotype) + (1|rep) + (1|rep/block), data = data, REML = TRUE)
  model.1 <- lmer(get(trait) ~ (1|genotype) + (1|env) + (1|env/rep) + (1|rep/block), data = data, REML = TRUE)
  #model.1 <- lmer(get(trait) ~ (1|genotype) + (1|env) + (1|env/rep) + (1|rep/block) + (1|genotype:env), data = data, REML = TRUE)
  
  
  
  
  
  
    # Decreasing stopping tolerances
  strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
  if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
  summary(model, correlation=FALSE)
  random_effects <- ranef(model)
  # Write out BLUPs for Genotypes
  write.table(random_effects$genotype, paste0("blups_", trait, ".csv"), col.names=F, row.names=T, sep=",")
  # Summary of random effects
  summary <- summary(model, correlation=FALSE)
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("resids_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  # Calculate heritability 
  model_variances <- as.data.frame(VarCorr(model))
  h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/5)+(model_variances$vcov[8]/10))
  out <- capture.output(h2)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  pdf(paste0("assumptions_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,3))
  
  # Model Fit with REML
  plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
  abline(h=0,col="red", lwd=1, lty=1)
  # histogram of residuals
  hist(residuals(model),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
  x=seq(-5e-15,9e-15,5e-15)
  curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
  # qq plot
  #qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
  
  dev.off()
  
}


############################################################
####### COMBINE BLUPS INTO MASTER FILES
############################################################

#Clear global environment
rm(list=ls(all=TRUE))
gc()

setwd("C:/Users/travi/Desktop/Genetics/TJH GxE Manuscript/2_Heatmaps/3env")

# because each trait has its own BLUP .csv output from the previous for-loop
# this script is designed to go find all of them and combine them into one .csv
# the first line here (data1) reads in 1 and pulls out some label info from it
# then it combined the rest by name. Make sure the working directory matches
# where all the blup .csv files are stored.
data1 <- read.csv("blups_c21hc.csv", header=F, stringsAsFactors=F)
rowlabels <- data1[,1]
blup_combined <- data.frame(matrix(ncol = 109, nrow = 468))
blup_combined[,1] <- rowlabels

collabels <- c("Taxa","c21hc","c22hc","c23hc","c24hc","c25hc","c26hc","c27hc","c28hc","c29hc",
               "c30hc","c31hc","c33hc","c35hc","c24_1_7hc","c24_1_9hc","c24_2_unkhc","c25_1_7hc",
               "c25_1_9hc","c27_1_7hc","c27_1_9hc","c27_2_unkhc","c29_1_7hc","c29_1_9hc",
               "c29_2_ahc","c29_2_bhc","c31_1_7hc","c31_1_9hc","c31_1_unkhc","c33_1_7hc","c33_1_9hc",
               "c33_1_unkhc","c35_1_7hc","c35_1_9hc","c16fatms","c17fatms","c18fatms","c20fatms","c22fatms",
               "c24fatms","c18_1_afa","c18_1_bfa","c18ohtms","c19ohtms","c25ohtms","c26ohtms","c27ohtms","c28ohtms",
               "sl","hc","hcsat","hcunsat","hceven","hcevensat","hcevenunsat","hcodd","hcoddsat",
               "hcoddunsat","fa","fasat","faunsat","faeven","faodd","oh","oheven","ohodd","7monoenes",
               "9monoenes","79alkenes","othermonoenes","hcmonoenes","dienes","hc_psl","hcsat_psl",
               "hcunsat_psl","hceven_psl","hcevensat_psl","hcevenunsat_psl","hcodd_psl","hcoddsat_psl",
               "hcoddunsat_psl","fa_psl","fasat_psl","faunsat_psl","faeven_psl","faodd_psl","oh_psl",
               "oheven_psl","ohodd_psl","7monoenes_psl","9monoenes_psl","hcsat_phc","hcunsat_phc",
               "hcevensat_phc","hcoddsat_phc","hcevenunsat_phc","hcoddunsat_phc","hcodd_phcsat",
               "hceven_phcsat","hcodd_phcunsat","hceven_phcunsat","hcmonoenes_phcunsat","dienes_phcunsat",
               "7monoenes_phcmonoenes","9monoenes_phcmonoenes","7monoenes_p79alkenes","9monoenes_p79alkenes",
               "21hc_p22fa","23hc_p24fa"
)

colnames(blup_combined) <- collabels


runlist <- c("c21hc","c22hc","c23hc","c24hc","c25hc","c26hc","c27hc","c28hc","c29hc",
             "c30hc","c31hc","c33hc","c35hc","c24_1_7hc","c24_1_9hc","c24_2_unkhc","c25_1_7hc",
             "c25_1_9hc","c27_1_7hc","c27_1_9hc","c27_2_unkhc","c29_1_7hc","c29_1_9hc",
             "c29_2_ahc","c29_2_bhc","c31_1_7hc","c31_1_9hc","c31_1_unkhc","c33_1_7hc","c33_1_9hc",
             "c33_1_unkhc","c35_1_7hc","c35_1_9hc","c16fatms","c17fatms","c18fatms","c20fatms","c22fatms",
             "c24fatms","c18_1_afa","c18_1_bfa","c18ohtms","c19ohtms","c25ohtms","c26ohtms","c27ohtms","c28ohtms",
             "sl","hc","hcsat","hcunsat","hceven","hcevensat","hcevenunsat","hcodd","hcoddsat",
             "hcoddunsat","fa","fasat","faunsat","faeven","faodd","oh","oheven","ohodd","X7monoenes",
             "X9monoenes","X79alkenes","othermonoenes","hcmonoenes","dienes","hc_psl","hcsat_psl",
             "hcunsat_psl","hceven_psl","hcevensat_psl","hcevenunsat_psl","hcodd_psl","hcoddsat_psl",
             "hcoddunsat_psl","fa_psl","fasat_psl","faunsat_psl","faeven_psl","faodd_psl","oh_psl",
             "oheven_psl","ohodd_psl","X7monoenes_psl","X9monoenes_psl","hcsat_phc","hcunsat_phc",
             "hcevensat_phc","hcoddsat_phc","hcevenunsat_phc","hcoddunsat_phc","hcodd_phcsat",
             "hceven_phcsat","hcodd_phcunsat","hceven_phcunsat","hcmonoenes_phcunsat","dienes_phcunsat",
             "X7monoenes_phcmonoenes","X9monoenes_phcmonoenes","X7monoenes_p79alkenes","X9monoenes_p79alkenes",
             "X21hc_p22fa","X23hc_p24fa"
)

j=2

for (i in runlist){

  filename1 <- ""
  filename <- ""
  filename1 <- paste("blups", i, sep = "_")
  filename <- paste(filename1, "csv", sep = ".")
  
  data <- read.csv(file=filename, header=F, stringsAsFactors=F)
  
  blup_combined[,j] <- data[,2]
  j <- j+1
  
}

setwd("C:/Users/travi/Desktop/Genetics/TJH GxE Manuscript/2_Heatmaps/3env")
write.csv(blup_combined,"allblups_3env.csv", row.names = FALSE)


