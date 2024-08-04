# get required libraries
library(foreign)
library(tidyverse)
library(INLA)
library(parallel)

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

# save the model object as an RDS file
saveRDS(Model_Matched_Woody, file = "output/model_matched_woody.rds")

# get summary and statistical significance
Summary_Model_Matched_Woody <- summary(Model_Matched_Woody)$fixed

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

# get summary and statistical significance
Summary_Model_Matched_Koala <- summary(Model_Matched_Koala)$fixed

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

# get summary and statistical significance
Summary_Model_Mixed_Woody <- summary(Model_Mixed_Woody)$fixed

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

# get summary and statistical significance
Summary_Model_Mixed_Koala <- summary(Model_Mixed_Koala)$fixed

# create some plots

# immediate effect matched

Coeff_Effects <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Matched_Woody["ba1:ci1", "mean"], Summary_Model_Matched_Koala["ba1:ci1", "mean"]), Lower = c(Summary_Model_Matched_Woody["ba1:ci1", "0.025quant"], Summary_Model_Matched_Koala["ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Matched_Woody["ba1:ci1", "0.975quant"], Summary_Model_Matched_Koala["ba1:ci1", "0.975quant"])))
Coeff_Effects <- Coeff_Effects %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

Plot <- ggplot(Coeff_Effects, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 20),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm"))

ggsave(Plot, file = "output/figures/immediate_effects_matched.jpg", width = 20, height = 25, units = "cm", dpi = 300)

# trend effect matched

Coeff_Effects <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Matched_Woody["time:ba1:ci1", "mean"], Summary_Model_Matched_Koala["time:ba1:ci1", "mean"]), Lower = c(Summary_Model_Matched_Woody["time:ba1:ci1", "0.025quant"], Summary_Model_Matched_Koala["time:ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Matched_Woody["time:ba1:ci1", "0.975quant"], Summary_Model_Matched_Koala["ba1:ci1", "0.975quant"])))
Coeff_Effects <- Coeff_Effects %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

Plot <- ggplot(Coeff_Effects, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 20),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm"))

ggsave(Plot, file = "output/figures/trend_effects_matched.jpg", width = 20, height = 25, units = "cm", dpi = 300)

# immediate effect mixed

Coeff_Effects <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["ba1:ci1", "mean"], Summary_Model_Mixed_Koala["ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.975quant"])))
Coeff_Effects <- Coeff_Effects %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

Plot <- ggplot(Coeff_Effects, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 20),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm"))

ggsave(Plot, file = "output/figures/immediate_effects_mixed.jpg", width = 20, height = 25, units = "cm", dpi = 300)

# trend effect mixed

Coeff_Effects <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "mean"], Summary_Model_Mixed_Koala["time:ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["time:ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.975quant"])))
Coeff_Effects <- Coeff_Effects %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

Plot <- ggplot(Coeff_Effects, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 20),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm"))

ggsave(Plot, file = "output/figures/trend_effects_mixed.jpg", width = 20, height = 25, units = "cm", dpi = 300)
















# fit models to koala habitat for separated properties

# load data
Data2 <- as_tibble(read.dbf("input/khab_match_BACI.dbf")) %>% mutate(LOSS_count = ifelse(LOSS_count > FOREST_cou, FOREST_cou, LOSS_count), CI = ifelse(CI == 0, 1, 0))

#Data2 <- as_tibble(read.dbf("input/BACItime_khab_match_oldtime.dbf")) %>% mutate(LOSS_count = ifelse(LOSS_count > FOREST_cou, FOREST_cou, LOSS_count))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model2 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data2, family = "binomial", Ntrials = FOREST_cou)

# get summary and statistical significance
Summary_Model2 <- summary(Model2)$fixed

# save the model object as an RDS file
saveRDS(Model2, file = "output/Model2.rds")

# fit models to all woody vegetation data for mixed properties

# load data
Data3 <- as_tibble(read.dbf("input/mixed_BACI.dbf")) %>% mutate(LOSS_count = ifelse(LOSS_count > FOREST_cou, FOREST_cou, LOSS_count))

# first fit the model to the full data set

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model3 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data3, family = "binomial", Ntrials = FOREST_cou, verbose=TRUE)

Model3 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME, data = Data3, family = "binomial", Ntrials = FOREST_cou, verbose=TRUE)

# get summary and statistical significance
Summary_Model3 <- summary(Model3)$fixed

# save the model object as an RDS file
saveRDS(Model3, file = "output/Model3.rds")

# now sub-sample so we have the same sample size as for the separated properties

# here to ensure that we have the same sample size for the mixed properties as the separated properties we randomly select 3,665 poperties without replacement and do this multiple times
#BootData3 <- randomise(Data = Data3, N = 3665, Reps = 1000)

# make cluster
#cl <- makeCluster(5)

# export packages and nimble functions to cluster
#clusterEvalQ(cl,{
#  library(tidyverse)
#  library(INLA)})

# fit model
#Models3 <- parLapply(cl = cl, X = BootData3, fun = fit_inla)

#stop cluster
#stopCluster(cl)

# get fixed effect parameter estimates
#Effects3 <- lapply(Models3, function(x) as.matrix(summary(x)$fixed[, "mean"]))
# combine inot a single matrix
#Effects3 <- do.call(cbind, Effects3)

# create a shell to store results
#Summary_Models3 <- summary(Models3[[1]])$fixed
#Summary_Models3[1:nrow(Summary_Models3), 1:ncol(Summary_Models3)] <- 0

# get summary statistics
#Summary_Models3[, "mean"] <- t(apply(Effects3, 1, mean))
#Summary_Models3[, "sd"] <- t(apply(Effects3, 1, sd))
#Summary_Models3[, "0.025quant"] <- t(apply(Effects3, 1, quantile, probs = 0.025))
#Summary_Models3[, "0.5quant"] <- t(apply(Effects3, 1, quantile, probs = 0.5))
#Summary_Models3[, "0.975quant"] <- t(apply(Effects3, 1, quantile, probs = 0.975))
#Summary_Models3[, "mode"] <- NA # don't calculate the mode

# save the model object as an RDS file
#saveRDS(Models3, file = "output/Models3.rds")



# fit models to koala habitat for mixed properties

# load data
Data4 <- as_tibble(read.dbf("input/khab_mix_BACI.dbf")) %>% mutate(LOSS_count = ifelse(LOSS_count > FOREST_cou, FOREST_cou, LOSS_count))

# first fit the model to the full data set

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model4 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data4, family = "binomial", Ntrials = FOREST_cou)

# get summary and statistical significance
Summary_Model4 <- summary(Model4)$fixed

# save the model object as an RDS file
saveRDS(Model4, file = "output/Model4.rds")

# now sub-sample so we have the same sample size as for the separated properties

# here to ensure that we have the same sample size for the mixed properties as the separated properties we randomly select 1,735 poperties without replacement and do this multuiple times
#BootData4 <- randomise(Data = Data4, N = 1735, Reps = 1000)

# make cluster
#cl <- makeCluster(5)

# export packages and nimble functions to cluster
#clusterEvalQ(cl,{
#  library(tidyverse)
#  library(INLA)})

# fit model
#Models4 <- parLapply(cl = cl, X = BootData4, fun = fit_inla)

#stop cluster
#stopCluster(cl)

# get fixed effect parameter estimates
#Effects4 <- lapply(Models4, function(x) as.matrix(summary(x)$fixed[, "mean"]))
# combine inot a single matrix
#Effects4 <- do.call(cbind, Effects4)

# create a shell to store results
#Summary_Models4 <- summary(Models4[[1]])$fixed
#Summary_Models4[1:nrow(Summary_Models4), 1:ncol(Summary_Models4)] <- 0

# get summary statistics
#Summary_Models4[, "mean"] <- t(apply(Effects4, 1, mean))
#Summary_Models4[, "sd"] <- t(apply(Effects4, 1, sd))
#Summary_Models4[, "0.025quant"] <- t(apply(Effects4, 1, quantile, probs = 0.025))
#Summary_Models4[, "0.5quant"] <- t(apply(Effects4, 1, quantile, probs = 0.5))
#Summary_Models4[, "0.975quant"] <- t(apply(Effects4, 1, quantile, probs = 0.975))
#Summary_Models4[, "mode"] <- NA # don't calculate the mode

# save the model object as an RDS file
#saveRDS(Models4, file = "output/Models4.rds")

# CREATE FIGURES

SummaryMod1 <- summary(Model1)$fixed
SummaryMod2 <- summary(Model2)$fixed

Coeff_Effects <- as_tibble(cbind(Type = c("All Woody Vegetation", "Koala Habitat"), Est = c(SummaryMod1["BA:CI", "mean"], SummaryMod2["BA:CI", "mean"]), Lower = c(SummaryMod1["BA:CI", "0.025quant"], SummaryMod2["BA:CI", "0.025quant"]), Upper = c(SummaryMod1["BA:CI", "0.975quant"], SummaryMod2["BA:CI", "0.975quant"])))
Coeff_Effects <- Coeff_Effects %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% group_by(Type)

Plot1 <- ggplot(Coeff_Effects, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect of Self Assessment") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 20),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm"))

ggsave(Plot1, file = "./figures/immediate_effects_seperated.jpg", width = 20, height = 25, units = "cm", dpi = 300)



