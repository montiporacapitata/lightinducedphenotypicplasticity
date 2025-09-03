###------------------------------------------------------------------------------
### Exploratory stats for 3D-based coral measures
###------------------------------------------------------------------------------

### Aim: to test assumptions of normality, homocedasticity and transform variables 

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## Load required libraries/packages

library(readr)
library(ggpubr)
library(stats)
library(car)
library(bestNormalize)
library(MVN)

## Set working directory: 
setwd("E:/PhD_Organized/I. Science/E. Morphology/2. Traits_extraction")

morphology <- as.data.frame(read_delim("data_raw_final.txt",
                                      "\t", escape_double = FALSE,
                                      na = "NA", trim_ws = TRUE))


#remove FD

ALLvars <- morphology %>% dplyr::select(-(`fractal_dimension`)) %>% na.omit()

##-------------------------------------------------------------------------------
## 2. Get log-ratios
##-------------------------------------------------------------------------------
## split df

df1_raw <- subset(ALLvars, time == "1yr")
df0_raw <- subset(ALLvars, time == "0yr")

df1_raw <- df1_raw[order(df1_raw$code),]
df0_raw <- df0_raw[order(df0_raw$code),]

#remove individuals that do not have both t1 and t0 scans due to high mortality or presence of NAs
df1 <- df1_raw[-c(37,39,41,44),] 
df0 <- df0_raw[-c(3,8,18,40,43),]

# get log ratios 
lr_rugosity <- log10(df1$rugosity/df0$rugosity)
lr_volume <- log10(df1$V/df0$V)
lr_SV_ratio <- orderNorm(log10(df1$SV_ratio/df0$SV_ratio))
lr_convexity <- log10(df1$convexity/df0$convexity)
lr_packing <- log10(df1$packing/df0$packing)
lr_toph <- log10(df1$top_heaviness/df0$top_heaviness)
lr_SA <- log10(df1$surface_area/df0$surface_area)

df_ratio <- as.data.frame(cbind(df1[,1:8],
                                lr_volume,
                                lr_SA,
                                lr_SV_ratio$x.t, 
                                lr_packing,
                                lr_rugosity,
                                lr_convexity,
                                lr_toph
                                ))

write.table(df_ratio, "habtoolsvarsTR_lr.txt",
            sep="\t",row.names = F)
##-------------------------------------------------------------------------------
## 3. Evaluate normality assumption of log-ratios
##-------------------------------------------------------------------------------

## A) Load extra packages to manipulate the data:
library(tidyr)
library(tidyselect)
library(tidyverse)
library(dplyr)
library(broom)

## B) Shapiro-Wilk test:
SWtest <- df_ratio %>% 
  gather(key = "variable_name", value = "value", 9:15) %>% 
  group_by(variable_name)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtest, "SWvarsTR_lr.txt",
            sep="\t",row.names = F)

##-------------------------------------------------------------------------------
## 4. Transform (TR) variables that violate the normality assumption
##-------------------------------------------------------------------------------
## For each variable with p-value <0.05: 
# A) Plot q-q plot and test normality using Shapiro-Wilk 
# B) Run bestNormalize for x5 and choose the most frequently suggested option
# C) Perform the suggested transformation
# D) Plot q-q plot again
# E) Plot q-q plot
# F) Perform Shapiro-Wilk test to check normality again

  
#top heaviness
ggqqplot(ALLvars$top_heaviness)
shapiro.test(ALLvars$top_heaviness)
bestNormalize(ALLvars$top_heaviness)
top_heaviness<-log10(ALLvars$top_heaviness)
ggqqplot(top_heaviness)
shapiro.test(top_heaviness)


#convexity
ggqqplot(ALLvars$convexity)
shapiro.test(ALLvars$convexity)
bestNormalize(ALLvars$convexity)
convexity<-log10(ALLvars$convexity)
ggqqplot(convexity)
shapiro.test(convexity)


#packing
ggqqplot(ALLvars$packing)
shapiro.test(ALLvars$packing)
bestNormalize(ALLvars$packing)
packing<-log10(ALLvars$packing)
ggqqplot(packing)
shapiro.test(packing)


#rugosity
ggqqplot(ALLvars$rugosity)
shapiro.test(ALLvars$rugosity)
bestNormalize(ALLvars$rugosity)
rugosity<-orderNorm(ALLvars$rugosity)
ggqqplot(rugosity$x.t)
shapiro.test(rugosity$x.t)


#surface area
ggqqplot(ALLvars$surface_area)
shapiro.test(ALLvars$surface_area)
bestNormalize(ALLvars$surface_area)
surface_area<-orderNorm(ALLvars$surface_area)
ggqqplot(surface_area$x.t)
shapiro.test(surface_area$x.t)

#V
ggqqplot(ALLvars$V)
shapiro.test(ALLvars$V)
bestNormalize(ALLvars$V)
Volume<-log10(ALLvars$V)
ggqqplot(Volume)
shapiro.test(Volume)

#SV_ratio
ggqqplot(ALLvars$SV_ratio)
shapiro.test(ALLvars$SV_ratio)
bestNormalize(ALLvars$SV_ratio)
SV_ratio<-center_scale(ALLvars$SV_ratio)
ggqqplot(SV_ratio$x.t)
shapiro.test(SV_ratio$x.t)



##-------------------------------------------------------------------------------
## 4. Create new data set with transformed variables (TR)
##-------------------------------------------------------------------------------

ALLvarsTR <- as.data.frame(cbind(ALLvars[,1:8],
                                Volume,
                                surface_area$x.t,
                                SV_ratio$x.t, 
                                packing,
                                rugosity$x.t,
                                convexity,
                                top_heaviness
))

#varsNOTnormal <- subset(SWtest, p.value < 0.05)
#varsNOTn_list <- as.factor(varsNOTnormal$variable_name) #List of TR variables
#varsNOTn_data <- dplyr::select(ALLvars, all_of(varsNOTn_list))
#ALLvarsTR <- dplyr::select(ALLvars,
 #                          -all_of(varsNOTn_list))

#ALLvarsTR <- cbind(ALLvarsTR, top_heaviness$x.t, convexity$x.t, 
#                   V$x.t,
#                   surface_area$x.t
#                  )

#ALLvarsTR <- ALLvarsTR[-c(16:19)]

# Save output
write.table(ALLvarsTR, "habtoolsvarsTR_log.txt",
            sep="\t", row.names = FALSE) # To export results 

##-------------------------------------------------------------------------------
## 6. Evaluate normality assumption in the transformed (TR) variables
##-------------------------------------------------------------------------------

SWtestTR <- ALLvarsTR %>% 
  gather(key = "variable_name", value = "value", 9:15) %>% 
  group_by(variable_name)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtestTR, "SWtestTR_log.txt",
            sep="\t", row.names = F)

SWtest_sppTR <- ALLvarsTR %>% 
  gather(key = "variable_name", value = "value", 9:15) %>% 
  group_by(variable_name, treatment_time)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtest_sppTR, "SWtest_sppTR.txt",
            sep="\t",row.names = F)

#BONUS: plot autocorrelation of non-transformed variables
#that checks for relationship with volume (ie size)
## A) Assess colors to treatments

Treatment <- factor(ALLvars$treatment_time,
                    levels = c("Control_0yr","Control_1yr","Shade_0yr","Shade_1yr"),
                    labels=c("C0","C1","S0","S1"))
cols <- c("grey", "goldenrod", "grey25", "midnightblue")
Genotype <- factor(ALLvars$genotype, levels = c("19","58","60","62","64","64Y","90","92","P1R2"))

## B) Check correlation and bi-plots

library(psych)

pairs.panels(ALLvars[,c(9,10,11,12,14,15)],
             gap = 0,
             hist.col="grey",
             bg = cols,
             pch = 21,
             ellipses=FALSE,
             density = TRUE,
             cex.cor=1,
             cex=1,
             stars = TRUE,
             smooth=FALSE)

