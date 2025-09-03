
##REFLECTANCE PROFILES

library(rlang)
library(ggplot2)
library(scales)
library(pavo)
library(cowplot)
library(photobiology)
library(photobiologyWavebands)
library(ggspectra)
library(ggrepel)
setwd("E:/PhD_Organized/I. Science/E. Morphology/8. Reflectance/Data_reflectance_Heidi")
db <- read.csv("CoralRod_B.csv", sep =";", head = TRUE)
head(db)

#SHADE
names2 <- gsub("Shade[0-9]+.[0-9]+", "", names(db[,c(1:76)])[-1])
names2 <- c("wl", names2)
aggplot(db[,c(1:76)], names2, lcol = c("grey", "purple4"), shadecol =c("grey", "purple4"), alpha = .2, lty = c(2,1), lwd = 1.5, ylim = c(0,1.1))

#Optional: Plot using various error functions and options
#aggplot(sicalis, bysic, FUN.error = function(x) quantile(x, c(0.0275, 0.975)))

#CONTROL
names1 <- gsub(".Control[0-9]+.[0-9]+", "", names(db[,c(77:152)]))
names1 <- c("wl", names1)
aggplot(db[,c(1, 77:152)], names1, lcol = c("grey", "goldenrod"), shadecol =c("grey", "goldenrod"), alpha = .2, lty = c(2,1), lwd = 1.5, ylim = c(0,1.1))

#Statistical tests


#datatest <- 
#write.table(datatest, "data_reflectance.txt", sep="\t",
# row.names = FALSE)

datatest2 <- read.table("data_reflectance.txt",sep="\t", header = TRUE)

#re-formatting data
ref_df <- c(datatest2[,2],datatest2[,3],datatest2[,4],datatest2[,5])
ref_df_ok <- cbind(ref_df, c(rep("Tip control", nrow(datatest2)), rep("Body control", nrow(datatest2)), rep("Tip Shade", nrow(datatest2)), rep("Body shade", nrow(datatest2))))

ref_df_ok <- as.data.frame(ref_df_ok)

reflectance <- as.numeric(as.character(ref_df_ok$ref_df))

library(bestNormalize)
mytest <- lm(log(reflectance) ~ V2, data = ref_df_ok)
plot(mytest)
bestNormalize::bestNormalize(reflectance)
shapiro.test(mytest$residuals) #not normal
reflectance_ok <- orderNorm(reflectance)
shapiro.test(mytest2$residuals)  #not normal: kruskal-wallis


kruskal.test(log(reflectance)~V2, data= ref_df_ok)
wilcox.test(reflectance~V2, data= ref_df_ok) #faire pairwise

mytest2 <- lm(reflectance_ok$x.t ~ V2, data = ref_df_ok)
plot(mytest2)
dens <- ggplot(ref_df_ok, aes(x = 4+log(reflectance), fill = V2)) + 
  geom_density(alpha = .8)+
  theme_cowplot(12)+
  scale_fill_manual(values = c("goldenrod", "midnightblue", "white", "grey"))+
  ylab("Density (%)")+
  xlab("4 + Log(Refectance)")+
  labs(fill = "Sample type")+
  theme(legend.position = c(0.8,0.5))
dens



            