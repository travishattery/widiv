## Original script - Travis Hattery
## 8 May 2022

## Conduct Spearman correlations and save estimates and p-values to .csv files 



############################################################
####### CLEAN ENVIRONMENT & FIND DATASETS ########
############################################################

# Clear global environment
rm(list=ls(all=TRUE))
gc()

# Find dataset
# Dataset should have:
#  Column 1: genotype
#  Column 2 to n: Data
setwd("C:/<WORKING DIRECTORY>")
input <- read.csv("allblups_3env.csv", header=T, stringsAsFactors=F)


traits <- colnames(input[2:109])

############################################################
####### BUILD DATAFRAMES ########
############################################################

# Build enpty dataframes to store results to
# the ncol and nrow will change depending on the number of traits you have.
# An incorrectly sized dataframe will cause problems.
# too small, and it will hit a stop-error. Too large, and it will have missing data.
metab_est <- data.frame(matrix(ncol = 108, nrow = 108))
colnames(metab_est) <- traits
rownames(metab_est) <- traits
metab_p <- data.frame(matrix(ncol = 108, nrow = 108))
colnames(metab_p) <- traits
rownames(metab_p) <- traits


############################################################
####### SPEARMAN CORRELATION ########
############################################################

# Loop through all traits (i) by all traits (j) and store to grids
for(i in 1:108){
  
  for(j in 1:108){
    xname <- traits[i]
    yname <- traits[j]
    
    x <- input[,xname]
    y <- input[,yname]
    
    spearman <- cor.test(x,y,  method = "spearman")
    
    metab_est[i,j] = spearman$estimate
    metab_p[i,j] = spearman$p.value
    
    spearman
    
  }
}

# Store results
write.csv(metab_est,"C:/<WORKING DIRECTORY>/spearman_est_3env.csv", row.names = TRUE)
write.csv(metab_p,"C:/<WORKING DIRECTORY>/spearman_p_3env.csv", row.names = TRUE)
