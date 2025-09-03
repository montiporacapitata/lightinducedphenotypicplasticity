###-------------------------------------------------------------------------------
### Inter-group differences 
###-------------------------------------------------------------------------------

### Aims: -Evaluate significant differences between light conditions at t0
###   -Make distribution boxplots of log-transformed variables for each group

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## A) Load required libraries/packages

library(readr)
library(stats)
library(rstatix)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)

## B) Set working directory: 
setwd("E:/PhD_Organized/I. Science/E. Morphology/2. Traits_extraction")

## C) Upload data set:
ALLvarsTR <- as.data.frame(read_delim("habtoolsvarsTR_log.txt", delim="\t"))

ALLvarsTR <- ALLvarsTR %>% dplyr::select(-(`packing`))

##-------------------------------------------------------------------------------
## 2. Assessing differences at t0
##-------------------------------------------------------------------------------

ALLvarsTR_0yr <- ALLvarsTR %>% filter(time == "0yr")

##Analysis of variance (ANOVA) between C0 and S0 (Supplementary Table S2)

# Store all formulae in a list
formulae <- lapply(colnames(ALLvarsTR_0yr)[9:ncol(ALLvarsTR_0yr)], function(x) as.formula(paste0("`",x,"`", " ~ treatment_time")))
ANOVA_resultsTR <- lapply(formulae, function(x) summary(aov(x, data = ALLvarsTR_0yr)))
names(ANOVA_resultsTR) <- format(formulae)

# Extract just p-values
p <- unlist(lapply(ANOVA_resultsTR, function(x) x[[1]]$"Pr(>F)"[1]))
dfpval <- cbind(variable=colnames(ALLvarsTR_0yr[,9:14]),p)

# Save output
write.table(dfpval, "ANOVA_results_log.txt", sep="\t",
            row.names = FALSE)


#Distribution boxplots (Fig. 2b)
ALLvarsTR <- ALLvarsTR %>%
  mutate(treatment_time = factor(treatment_time,
                          levels = c("Control_0yr","Control_1yr", "Shade_0yr","Shade_1yr"),
                          labels = c("C0","C1", "S0","S1"))) %>%
  dplyr::select(-Sample_ID,
                -`treatment`,
                -`code`,
                -genotype,
                -`rack`,
                -`time`,
                -`mortality`,
               ) 

ALLvarsTR_tibb <-  as_tibble(ALLvarsTR)
ALLvarsTR_tibb %>% sample_n(6) # Check

mydata.long <- ALLvarsTR_tibb %>%
  tidyr::pivot_longer(-treatment_time,
                      names_to = "variables",values_to = "value")

mydata.long %>% sample_n(6) # Check

mydata.long %>%
  group_by(variables, treatment_time) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value)
  ) %>%
  ungroup()

#stat.test <- mydata.long %>%
#  group_by(variables) %>%
#  t_test(value ~ treatment_time, p.adjust.method = "none")


# Remove unnecessary columns and display the outputs
#ttest <- stat.test %>%
#  dplyr::select(-.y., -statistic, -df) 

# Save output
#write.table(ttest, "t-test_Bonf-p-adj.txt",
#            sep="\t", row.names = F)

## D) OPTIONAL: Post hoc test with letters
#library(agricolae)

#List <- names(ALLvarsTR)[9:14] # select just the variables

#model1 <- lapply(List, function(x) {
#  lm(substitute(i~treatment_time, list(i = as.name(x))), data = ALLvarsTR)})

#lapply(model1, summary)

#letters = lapply(model1, function(m) HSD.test((m), "treatment_time", group = TRUE, console = TRUE))

#plot it
myplot <- ggboxplot(
  mydata.long, x = "treatment_time", y = "value",
  error.plot = "crossbar", notch=FALSE,
  fill = "treatment_time", legend = "none",
  ggtheme = theme_half_open(12)) +
  geom_jitter(stroke = 0.01, alpha = 0.2) +
  facet_wrap(~variables, scales="free_y") +
  ylab("Log10(value)")+ xlab("Group")+
  scale_fill_manual(name="treatment_time",
                    values = c("grey","goldenrod", "grey51", "midnightblue"),
                    labels = c("C0",
                               "C1",
                               "S0",
                               "S1")) 
myplot

##Add statistical test p-values
#stat.test <- stat.test %>% add_xy_position(x = "treatment_time")
#myplot + stat_pvalue_manual(stat.test, size = 2, label = "p.adj.signif")

