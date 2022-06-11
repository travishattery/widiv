# Original TJH 210507
# Revised TJH 210602
# 1) LOESS fit data, extract values at each GDD point
# 2) Type 2 ANOVA on data for multiple metabolite traits
# 3) Run ANOVA using significant terms from previous ANOVA with Interactions









####################
####################
####################
####################
####################
##### Import weather data, create loess curve fit, extract new data at every GDD

##### Setup working directory and import original data
setwd("C:/<WORKING DIRECTORY>")
gdd_data_ISU16 <- read.csv("loess_ISU16.csv", header=T)
gdd_data_UMN16 <- read.csv("loess_UMN16.csv", header=T)
gdd_data_ISU17 <- read.csv("loess_ISU17.csv", header=T)

##### Create curves for 10% spans
loessMod10_ISU16 <- loess(tsoil ~ gdd_cmtv, data=gdd_data_ISU16, span=0.10) # 10% smoothing span
loessMod10_UMN16 <- loess(tsoil ~ gdd_cmtv, data=gdd_data_UMN16, span=0.10) # 10% smoothing span
loessMod10_ISU17 <- loess(tsoil ~ gdd_cmtv, data=gdd_data_ISU17, span=0.10) # 10% smoothing span

##### Find values of loess curves at original values of x
smoothed10_ISU16 <- predict(loessMod10_ISU16)
smoothed10_UMN16 <- predict(loessMod10_UMN16)
smoothed10_ISU17 <- predict(loessMod10_ISU17)

##### Plot smoothed loess curves on top of original data
plot(gdd_data_UMN16$tsoil, x=gdd_data_UMN16$gdd_cmtv, type="l", main="Loess Smoothing and Prediction", xlab="GDD", ylab="weather parameter")
lines(smoothed10_UMN16, x=gdd_data_UMN16$gdd_cmtv, col="green")

##### Find values of loess curves at each GDD value, write these values to .csv
gdd_trans_ISU16 <- predict(loessMod10_ISU16, data.frame(gdd_cmtv = seq(1, 1939, 1)), se = FALSE)
write.csv(gdd_trans_ISU16,"C:/<WORKING DIRECTORY>/gdd_trans_ISU16.csv", row.names = TRUE)

gdd_trans_UMN16 <- predict(loessMod10_UMN16, data.frame(gdd_cmtv = seq(1, 1640, 1)), se = FALSE)
write.csv(gdd_trans_UMN16,"C:/<WORKING DIRECTORY>/gdd_trans_UMN16.csv", row.names = TRUE)

gdd_trans_ISU17 <- predict(loessMod10_ISU17, data.frame(gdd_cmtv = seq(1, 2155, 1)), se = FALSE)
write.csv(gdd_trans_ISU17,"C:/<WORKING DIRECTORY>/gdd_trans_ISU17.csv", row.names = TRUE)











###### EVERYTHING BELOW HERE IS USELESS TO YOU PROBABLY #####



####################
####################
####################
####################
####################
##### PLOT DENDROGRAMS OF WEATHER PARAMETER GDD-WINDOWS
##### TJH 210806


#install.packages("heatmaply")
library("heatmaply")

setwd("C:/<WORKING DIRECTORY>")
metab_wthr <- read.csv("WiDiv_3env_log10_wthr50_trim.csv", header=T)

#rain, rad, tday, tnight, humid, tsoil

#IF 25-GDD-unit windows
weathercorr <- metab_wthr[ , c("Sample_.ID",
                               "tsoil_1100_1125",
                               "tsoil_1075_1100", "tsoil_1050_1075", 
                               "tsoil_1025_1050", "tsoil_1000_1025", 
                               "tsoil_975_1000", "tsoil_950_975", 
                               "tsoil_925_950", "tsoil_900_925", 
                               "tsoil_875_900", "tsoil_850_875", 
                               "tsoil_825_850", "tsoil_800_825", 
                               "tsoil_775_800", "tsoil_750_775", 
                               "tsoil_725_750", "tsoil_700_725", 
                               "tsoil_675_700", "tsoil_650_675", 
                               "tsoil_625_650", "tsoil_600_625", 
                               "tsoil_575_600", "tsoil_550_575", 
                               "tsoil_525_550", "tsoil_500_525", 
                               "tsoil_475_500", "tsoil_450_475",
                               "tsoil_425_450", "tsoil_400_425",
                               "tsoil_375_400", "tsoil_350_375",
                               "tsoil_325_350", "tsoil_300_325",
                               "tsoil_275_300", "tsoil_250_275",
                               "tsoil_225_250", "tsoil_200_225",
                               "tsoil_175_200", "tsoil_150_175",
                               "tsoil_125_150", "tsoil_100_125",
                               "tsoil_75_100", "tsoil_50_75",
                               "tsoil_25_50", "tsoil_0_25")]

#IF 50-GDD-unit windows
weathercorr <- metab_wthr[ , c("Sample_.ID", "tsoil_1100_1150",
                               "tsoil_1050_1100", "tsoil_1000_1050",
                               "tsoil_950_1000", "tsoil_900_950",
                               "tsoil_850_900", "tsoil_800_850",
                               "tsoil_750_800", "tsoil_700_750",
                               "tsoil_650_700", "tsoil_600_650",
                               "tsoil_550_600", "tsoil_500_550",
                               "tsoil_450_500", "tsoil_400_450",
                               "tsoil_350_400", "tsoil_300_350",
                               "tsoil_250_300", "tsoil_200_250",
                               "tsoil_150_200", "tsoil_100_150",
                               "tsoil_50_100", "tsoil_0_50")]

rownames(weathercorr) <- weathercorr[,1]
weathercorr[,1] <- NULL

my_cor <- cor(weathercorr)
heatmaply_cor(my_cor)




######Misc tests I'm no longer using
#Calculate Correlation Coefficients
res <- cor(weathercorr)
round(res, 2)

#Calculate Correlation Coefficients AND p-values, Save into new dataframes
res2 <- rcorr(as.matrix(weathercorr))
weathercorr_r <- res2$r
weathercorr_p <- res2$P

#Plot correlation matrix table
chart.Correlation(weathercorr, histogram=TRUE, pch=19)

#spearman test: -1 to +1
cor(metab_wthr$rain_100_125, metab_wthr$rain_150_175,  method = "spearman", use = "pairwise.complete.obs")





