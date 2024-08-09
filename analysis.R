# get required libraries
library(foreign)
library(tidyverse)
library(INLA)
library(parallel)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)

# read functions
source("functions.R")

# fit models to matched woody vegetation data

# load data and add BA field
# format data fields
Data_Matched_Woody <- as_tibble(readRDS("input/joined_data_frames/woody_matched_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Matched_Woody <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Matched_Woody, family = "betabinomial", Ntrials = baseline)

Model_Matched_Woody <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + f(RID, model = "iid"), data = Data_Matched_Woody, family = "betabinomial", Ntrials = baseline)

# save the model object as an RDS file
saveRDS(Model_Matched_Woody, file = "output/model_matched_woody.rds")

# fit models to matched koala habitat data

# load data and add BA field
# format data fields
Data_Matched_Koala <- as_tibble(readRDS("input/joined_data_frames/khab_matched_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Matched_Koala <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Matched_Koala, family = "betabinomial", Ntrials = baseline)

# save the model object as an RDS file
saveRDS(Model_Matched_Koala, file = "output/model_matched_koala.rds")

# fit models to mixed woody vegetation data

# load data and add BA field
# format data fields
Data_Mixed_Woody <- as_tibble(readRDS("input/joined_data_frames/woody_mixed_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# randomly select 100,000 properties without replacement
Data_Mixed_Woody_Sub <- randomise(Data = Data_Mixed_Woody, N = 100000, Reps = 1)

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Mixed_Woody <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Mixed_Woody_Sub[[1]], family = "betabinomial", Ntrials = baseline)

# save the model object as an RDS file
saveRDS(Model_Mixed_Woody, file = "output/model_mixed_woody.rds")

# fit models to mixed koala habitat data

# load data and add BA field
# format data fields
Data_Mixed_Koala <- as_tibble(readRDS("input/joined_data_frames/khab_mixed_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# randomly select 100,000 properties without replacement
Data_Mixed_Koala_Sub <- randomise(Data = Data_Mixed_Koala, N = 100000, Reps = 1)

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Mixed_Koala <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Mixed_Koala_Sub[[1]], family = "betabinomial", Ntrials = baseline)

# save the model object as an RDS file
saveRDS(Model_Mixed_Koala, file = "output/model_mixed_koala.rds")

# do the predictions

# load models if needed
Model_Matched_Woody <- readRDS("output/model_matched_woody.rds")
Model_Matched_Koala <- readRDS("output/model_matched_koala.rds")
Model_Mixed_Woody <- readRDS("output/model_mixed_woody.rds")
Model_Mixed_Koala <- readRDS("output/model_mixed_koala.rds")

# get model summaries
Summary_Model_Matched_Woody <- summary(Model_Matched_Woody)$fixed
Summary_Model_Matched_Koala <- summary(Model_Matched_Koala)$fixed
Summary_Model_Mixed_Woody <- summary(Model_Mixed_Woody)$fixed
Summary_Model_Mixed_Koala <- summary(Model_Mixed_Koala)$fixed

# load property data 
zonal_woody_treat <- read_rds("input/analysis_data/woody_treat_zonal.rds")
zonal_woody_contr  <- read_rds("input/analysis_data/woody_contr_zonal.rds")
zonal_koala_treat <- read_rds("input/analysis_data/koala_treat_zonal.rds")
zonal_koala_contr <- read_rds("input/analysis_data/koala_contr_zonal.rds")

# get the treatment woody vegetation and koala habitat for mixed properties
Woody_Treat_Mixed <- zonal_woody_treat[which(zonal_woody_treat$baseline > 0 & zonal_woody_contr$baseline > 0), ] %>% filter(baseline > 0)  
Koala_Treat_Mixed <- zonal_koala_treat[which(zonal_koala_treat$baseline > 0 & zonal_koala_contr$baseline > 0), ] %>% filter(baseline > 0)

# get the treatment woody vegetation and koala habitat for matched properties 
Woody_Treat_Matched <- zonal_woody_treat[which((zonal_woody_treat$baseline > 0 & zonal_woody_contr$baseline == 0) | (zonal_woody_treat$baseline == 0 & zonal_woody_contr$baseline > 0)), ] %>% filter(baseline > 0)
Koala_Treat_Matched <- zonal_koala_treat[which((zonal_koala_treat$baseline > 0 & zonal_koala_contr$baseline == 0) | (zonal_koala_treat$baseline == 0 & zonal_koala_contr$baseline > 0)), ] %>% filter(baseline > 0)

# get the relevant coeficients from the models and calculate the impact of the treatment

# matched woody - no significant variables

# get current proportions cleared
p <- Woody_Treat_Matched %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021)

# get immediate effect coefficients
ImmediateLow <- 0
ImmediateMid <- 0
ImmediateHigh <- 0

# get trend effect coefficients
TrendLow <- c(0, 0, 0, 0, 0, 0, 0) * 0:6
TrendMid <- c(0, 0, 0, 0, 0, 0, 0) * 0:6
TrendHigh <- c(0, 0, 0, 0, 0, 0, 0) * 0:6

# matrix to hold new proportions cleared after removing the effect of the treatment
p_newLow <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newMid <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newHigh <- matrix(0, nrow = nrow(p), ncol = ncol(p))

# loop through propeties and time steps and calculate new proportions cleared after removing the effect of the treatment
for (i in 1:nrow(p)) {
  for (j in 1:ncol(p)) {
      if (p[i, j] == 1) {
        p_old <- 0.999999999999999 # dealing with proportions = 1
      } else if (p[i, j] == 0) {
        p_old <- 0.000000000000001 # dealing with proportions = 1  
      } else {
        p_old <- p[i, j]    
      }
      p_newLow[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j])))))
      p_newMid[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j])))))
      p_newHigh[i, j] <- as.numeric((((exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]))))))
  }
}

# sum values to get total effect and calculate areas
EffectLow_MatchedWoody <- (Woody_Treat_Matched[, 7:13] - p_newLow * matrix(rep(Woody_Treat_Matched$baseline, 7), nrow = nrow(Woody_Treat_Matched), ncol = 7)) * 25 * 25 / 10000
EffectMid_MatchedWoody <- (Woody_Treat_Matched[, 7:13] - p_newMid * matrix(rep(Woody_Treat_Matched$baseline, 7), nrow = nrow(Woody_Treat_Matched), ncol = 7)) * 25 * 25 / 10000
EffectHigh_MatchedWoody <- (Woody_Treat_Matched[, 7:13] - p_newHigh * matrix(rep(Woody_Treat_Matched$baseline, 7), nrow = nrow(Woody_Treat_Matched), ncol = 7)) * 25 * 25 / 10000

# matched koala - only immediate effect significant

# get current proportions cleared
p <- Koala_Treat_Matched %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021)

# get immediate effect coefficients
ImmediateLow <- Summary_Model_Matched_Koala[5, "0.025quant"]
ImmediateMid <- Summary_Model_Matched_Koala[5, "mean"]
ImmediateHigh <- Summary_Model_Matched_Koala[5, "0.975quant"]

# get trend effect coefficients
TrendLow <- c(0, 0, 0, 0, 0, 0, 0) * 0:6
TrendMid <- c(0, 0, 0, 0, 0, 0, 0) * 0:6
TrendHigh <- c(0, 0, 0, 0, 0, 0, 0) * 0:6

# matrix to hold new proportions cleared after removing the effect of the treatment
p_newLow <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newMid <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newHigh <- matrix(0, nrow = nrow(p), ncol = ncol(p))

# loop through propeties and time steps and calculate new proportions cleared after removing the effect of the treatment
for (i in 1:nrow(p)) {
  for (j in 1:ncol(p)) {
      if (p[i, j] == 1) {
        p_old <- 0.999999999999999 # dealing with proportions = 1
      } else if (p[i, j] == 0) {
        p_old <- 0.000000000000001 # dealing with proportions = 1  
      } else {
        p_old <- p[i, j]    
      }
      p_newLow[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j])))))
      p_newMid[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j])))))
      p_newHigh[i, j] <- as.numeric((((exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]))))))
  }
}

# sum values to get total effect and calculate areas
EffectLow_MatchedKoala <- (Koala_Treat_Matched[, 7:13] - p_newLow * matrix(rep(Koala_Treat_Matched$baseline, 7), nrow = nrow(Koala_Treat_Matched), ncol = 7)) * 25 * 25 / 10000
EffectMid_MatchedKoala <- (Koala_Treat_Matched[, 7:13] - p_newMid * matrix(rep(Koala_Treat_Matched$baseline, 7), nrow = nrow(Koala_Treat_Matched), ncol = 7)) * 25 * 25 / 10000
EffectHigh_MatchedKoala <- (Koala_Treat_Matched[, 7:13] - p_newHigh * matrix(rep(Koala_Treat_Matched$baseline, 7), nrow = nrow(Koala_Treat_Matched), ncol = 7)) * 25 * 25 / 10000

# mixed woody - immediate and trend effects significant

# get current proportions cleared
p <- Woody_Treat_Mixed %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021)

ImmediateLow <- Summary_Model_Mixed_Woody[5, "0.025quant"]
ImmediateMid <- Summary_Model_Mixed_Woody[5, "mean"]
ImmediateHigh <- Summary_Model_Mixed_Woody[5, "0.975quant"]

TrendLow <- Summary_Model_Mixed_Woody[8, "0.025quant"] * 0:6
TrendMid <- Summary_Model_Mixed_Woody[8, "mean"] * 0:6
TrendHigh <- Summary_Model_Mixed_Woody[8, "0.975quant"] * 0:6

# matrix to hold new proportions cleared after removing the effect of the treatment
p_newLow <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newMid <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newHigh <- matrix(0, nrow = nrow(p), ncol = ncol(p))

# loop through propeties and time steps and calculate new proportions cleared after removing the effect of the treatment
for (i in 1:nrow(p)) {
  for (j in 1:ncol(p)) {
      if (p[i, j] == 1) {
        p_old <- 0.999999999999999 # dealing with proportions = 1
      } else if (p[i, j] == 0) {
        p_old <- 0.000000000000001 # dealing with proportions = 1  
      } else {
        p_old <- p[i, j]    
      }
      p_newLow[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j])))))
      p_newMid[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j])))))
      p_newHigh[i, j] <- as.numeric((((exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]))))))
  }
}

# sum values to get total effect and calculate areas
EffectLow_MixedWoody <- (Woody_Treat_Mixed[, 7:13] - p_newLow * matrix(rep(Woody_Treat_Mixed$baseline, 7), nrow = nrow(Woody_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000
EffectMid_MixedWoody <- (Woody_Treat_Mixed[, 7:13] - p_newMid * matrix(rep(Woody_Treat_Mixed$baseline, 7), nrow = nrow(Woody_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000
EffectHigh_MixedWoody <- (Woody_Treat_Mixed[, 7:13] - p_newHigh * matrix(rep(Woody_Treat_Mixed$baseline, 7), nrow = nrow(Woody_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000

# mixed koala habitat - immediate and trend effects significant

# get current proportions cleared
p <- Koala_Treat_Mixed %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021)

ImmediateLow <- Summary_Model_Mixed_Koala[5, "0.025quant"]
ImmediateMid <- Summary_Model_Mixed_Koala[5, "mean"]
ImmediateHigh <- Summary_Model_Mixed_Koala[5, "0.975quant"]

TrendLow <- Summary_Model_Mixed_Koala[8, "0.025quant"] * 0:6
TrendMid <- Summary_Model_Mixed_Koala[8, "mean"] * 0:6
TrendHigh <- Summary_Model_Mixed_Koala[8, "0.975quant"] * 0:6

# matrix to hold new proportions cleared after removing the effect of the treatment
p_newLow <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newMid <- matrix(0, nrow = nrow(p), ncol = ncol(p))
p_newHigh <- matrix(0, nrow = nrow(p), ncol = ncol(p))

# loop through propeties and time steps and calculate new proportions cleared after removing the effect of the treatment
for (i in 1:nrow(p)) {
  for (j in 1:ncol(p)) {
      if (p[i, j] == 1) {
        p_old <- 0.999999999999999 # dealing with proportions = 1
      } else if (p[i, j] == 0) {
        p_old <- 0.000000000000001 # dealing with proportions = 1  
      } else {
        p_old <- p[i, j]    
      }
      p_newLow[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateLow - TrendLow[j])))))
      p_newMid[i, j] <- as.numeric(((exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateMid - TrendMid[j])))))
      p_newHigh[i, j] <- as.numeric((((exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]) / (1 + exp(log(p_old / (1 - p_old)) - ImmediateHigh - TrendHigh[j]))))))
  }
}

# sum values to get total effect and calculate areas
EffectLow_MixedKoala <- (Koala_Treat_Mixed[, 7:13] - p_newLow * matrix(rep(Koala_Treat_Mixed$baseline, 7), nrow = nrow(Koala_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000
EffectMid_MixedKoala <- (Koala_Treat_Mixed[, 7:13] - p_newMid * matrix(rep(Koala_Treat_Mixed$baseline, 7), nrow = nrow(Koala_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000
EffectHigh_MixedKoala <- (Koala_Treat_Mixed[, 7:13] - p_newHigh * matrix(rep(Koala_Treat_Mixed$baseline, 7), nrow = nrow(Koala_Treat_Mixed), ncol = 7)) * 25 * 25 / 10000

# get the total amounts cleared
Total_Koala_Low <- sum(EffectLow_MatchedKoala) + sum(EffectLow_MixedKoala)
Total_Koala_Mid <- sum(EffectMid_MatchedKoala) + sum(EffectMid_MixedKoala)
Total_Koala_High <- sum(EffectHigh_MatchedKoala) + sum(EffectHigh_MixedKoala)
Total_Woody_Low <- sum(EffectLow_MatchedWoody) + sum(EffectLow_MixedWoody)
Total_Woody_Mid <- sum(EffectMid_MatchedWoody) + sum(EffectMid_MixedWoody)
Total_Woody_High <- sum(EffectHigh_MatchedWoody) + sum(EffectHigh_MixedWoody)

# create some plots

# load models if needed
Model_Matched_Woody <- readRDS("output/model_matched_woody.rds")
Model_Matched_Koala <- readRDS("output/model_matched_koala.rds")
Model_Mixed_Woody <- readRDS("output/model_mixed_woody.rds")
Model_Mixed_Koala <- readRDS("output/model_mixed_koala.rds")

# get model summaries
Summary_Model_Matched_Woody <- summary(Model_Matched_Woody)$fixed
Summary_Model_Matched_Koala <- summary(Model_Matched_Koala)$fixed
Summary_Model_Mixed_Woody <- summary(Model_Mixed_Woody)$fixed
Summary_Model_Mixed_Koala <- summary(Model_Mixed_Koala)$fixed

# matched properties

# immediate effect

Coeff_EffectsA <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Matched_Woody["ba1:ci1", "mean"], Summary_Model_Matched_Koala["ba1:ci1", "mean"]), Lower = c(Summary_Model_Matched_Woody["ba1:ci1", "0.025quant"], Summary_Model_Matched_Koala["ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Matched_Woody["ba1:ci1", "0.975quant"], Summary_Model_Matched_Koala["ba1:ci1", "0.975quant"])))
Coeff_EffectsA <- Coeff_EffectsA %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotA <- ggplot(Coeff_EffectsA, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Immediate Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, by = 5))

# trend effect

Coeff_EffectsB <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Matched_Woody["time:ba1:ci1", "mean"], Summary_Model_Matched_Koala["time:ba1:ci1", "mean"]), Lower = c(Summary_Model_Matched_Woody["time:ba1:ci1", "0.025quant"], Summary_Model_Matched_Koala["time:ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Matched_Woody["time:ba1:ci1", "0.975quant"], Summary_Model_Matched_Koala["time:ba1:ci1", "0.975quant"])))
Coeff_EffectsB <- Coeff_EffectsB %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotB <- ggplot(Coeff_EffectsB, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Trend Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-10, 30), breaks = seq(-10, 30, by = 5))

# combined plot

Plot <- ggarrange(PlotA, PlotB, ncol = 1, nrow = 2, vjust = -1, hjust = 0, font.label = list(size = 18))
ggsave(Plot, file = "output/figures/effects_matched.jpg", width = 20, height = 30, units = "cm", dpi = 300)

# mixed properties

# immediate effect

Coeff_EffectsA <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["ba1:ci1", "mean"], Summary_Model_Mixed_Koala["ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.975quant"])))
Coeff_EffectsA <- Coeff_EffectsA %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotA <- ggplot(Coeff_EffectsA, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Immediate Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-0.05, 0.35), breaks = seq(-0.05, 0.35, by = 0.05))

# trend effect

Coeff_EffectsB <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "mean"], Summary_Model_Mixed_Koala["time:ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["time:ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["time:ba1:ci1", "0.975quant"])))
Coeff_EffectsB <- Coeff_EffectsB %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotB <- ggplot(Coeff_EffectsB, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Trend Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-0.15, 0), breaks = seq(-0.15, 0, by = 0.05))

# combined plot

Plot <- ggarrange(PlotA, PlotB, ncol = 1, nrow = 2, vjust = -1, hjust = 0, font.label = list(size = 18))
ggsave(Plot, file = "output/figures/effects_mixed.jpg", width = 20, height = 30, units = "cm", dpi = 300)

# plots of trends

# woody vegetation

EffectsWoody <- colSums(EffectMid_MixedWoody) + colSums(EffectMid_MatchedWoody)
EffectsWoody <- tibble(Clearing = EffectsWoody) %>% mutate(Year = 2015:2021)

Plot <- ggplot(EffectsWoody, aes(x = Year, y = Clearing)) + geom_line(color="red", linewidth = 1.5) + geom_point(size = 4) + theme_minimal() + labs(x = "Year", y = "Area Cleared (ha)") + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0) + scale_x_continuous(breaks = seq(2015, 2021, by = 1)) + scale_y_continuous(breaks = seq(-3500, 1000, by = 500))

ggsave(Plot, file = "output/figures/trend_woody.jpg", width = 30, height = 20, units = "cm", dpi = 300)

# koala habitat

EffectsKoala <- colSums(EffectMid_MixedKoala) + colSums(EffectMid_MatchedKoala)
EffectsKoala <- tibble(Clearing = EffectsKoala) %>% mutate(Year = 2015:2021)

Plot <- ggplot(EffectsKoala, aes(x = Year, y = Clearing)) + geom_line(color="red", linewidth = 1.5) + geom_point(size = 4) + theme_minimal() + labs(x = "Year", y = "Area Cleared (ha)") + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0) + scale_x_continuous(breaks = seq(2015, 2021, by = 1)) + scale_y_continuous(breaks = seq(0, 3000, by = 500))

ggsave(Plot, file = "output/figures/trend_koala.jpg", width = 30, height = 20, units = "cm", dpi = 300)
