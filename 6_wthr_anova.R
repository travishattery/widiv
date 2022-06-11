
#directory, libraries
setwd("C:/<WORKING DIRECTORY>")
library(car)


#input datasets and labelsets
input <- read.csv("2_widiv_3env_plotavg_log10_wthr50.csv")
wthranova1 <- read.csv("3_wthranova1.csv", header=TRUE, stringsAsFactors=FALSE)
collabels <- read.csv("wthrlabels.csv")
rowlabels <- c("genotype", "wthr_0_50","wthr_50_100","wthr_100_150","wthr_150_200",
           "wthr_200_250","wthr_250_300","wthr_300_350","wthr_350_400","wthr_400_450",
           "wthr_450_500","wthr_500_550","wthr_550_600","wthr_600_650","wthr_650_700",
           "wthr_700_750","wthr_750_800","wthr_800_850","wthr_850_900","wthr_900_950",
           "wthr_950_1000","wthr_1000_1050","wthr_1050_1100","wthr_1100_1150",
           "sig1","sig2","sig3","sig4","residual")


#Create dataframe to store ANOVA results
awthrresults <- data.frame(matrix(ncol = 541, nrow = 29))
colnames(awthrresults) <- colnames(collabels)
awthrresults[,1] <- rowlabels


sapply(input, typeof)
#rain=1:108 ... rad=109:216 ... tday=217:324 ... tnight=325:432 ... tsoil=433:540
#k values for loop storage: rain=1, rad=109, tday=217, tnight=325, tsoil=433
#1-10 DONE
for (i in 109:109){
  
  temp1.aov <- glm(as.formula(wthranova1[i,1]), data=input)
  temp2.aov <- car::Anova(temp1.aov, test.statistic="F")
  
  k=109
  dfcol = (i-k)*5+2
  sscol = (i-k)*5+3
  fcol = ((i-(k-1))*5)
  pcol = ((i-(k-1))*5)+1
  
  awthrresults[,dfcol]<-temp2.aov[,2] #df
  awthrresults[,sscol]<-temp2.aov[,1] #ss
  awthrresults[,fcol]<-temp2.aov[,3] #F
  awthrresults[,pcol]<-temp2.aov[,4] #p
  
}

write.csv(awthrresults,"widiv_aresults_wthr.csv", row.names = FALSE)
