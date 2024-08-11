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
Model_Matched_Woody <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Matched_Woody, family = "betabinomial", Ntrials = baseline, control.compute = list(config = TRUE))

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
Model_Matched_Koala <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Matched_Koala, family = "betabinomial", Ntrials = baseline, control.compute = list(config = TRUE))

# save the model object as an RDS file
saveRDS(Model_Matched_Koala, file = "output/model_matched_koala.rds")

# fit models to mixed woody vegetation data

# load data and add BA field
# format data fields
Data_Mixed_Woody <- as_tibble(readRDS("input/joined_data_frames/woody_mixed_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# randomly select 10,000 properties without replacement
Data_Mixed_Woody_Sub <- randomise(Data = Data_Mixed_Woody, N = 10000, Reps = 1)

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Mixed_Woody <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Mixed_Woody_Sub[[1]], family = "betabinomial", Ntrials = baseline, control.compute = list(config = TRUE))

# save the model object as an RDS file
saveRDS(Model_Mixed_Woody, file = "output/model_mixed_woody_new.rds")

# fit models to mixed koala habitat data

# load data and add BA field
# format data fields
Data_Mixed_Koala <- as_tibble(readRDS("input/joined_data_frames/khab_mixed_df.rds")) %>% mutate(ba = ifelse(time < 0, 0, 1)) %>% rename(ci = CI) %>% mutate(RID = as.integer(RID), baseline = as.integer(baseline), loss = as.integer(loss), time = as.integer(time), ba = as.factor(ba), ci = as.factor(ci))

# randomly select 10,000 properties without replacement
Data_Mixed_Koala_Sub <- randomise(Data = Data_Mixed_Koala, N = 10000, Reps = 1)

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for RID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model_Mixed_Koala <- inla(loss ~ time + ba + ci + ba:ci + ba:time + ci:time + ba:ci:time + f(RID, model = "iid"), data = Data_Mixed_Koala_Sub[[1]], family = "betabinomial", Ntrials = baseline, control.compute = list(config = TRUE))

# save the model object as an RDS file
saveRDS(Model_Mixed_Koala, file = "output/model_mixed_koala_new.rds")

# get the p-values

# load models if needed
Model_Matched_Woody <- readRDS("output/model_matched_woody.rds")
Model_Matched_Koala <- readRDS("output/model_matched_koala.rds")
Model_Mixed_Woody <- readRDS("output/model_mixed_woody.rds")
Model_Mixed_Koala <- readRDS("output/model_mixed_koala.rds")

print("Matched Woody")
inla.pmarginal(0, Model_Matched_Woody$marginals.fixed$`ba1:ci1`)
inla.pmarginal(0, Model_Matched_Woody$marginals.fixed$`time:ba1:ci1`)

print("Matched Koala")
inla.pmarginal(0, Model_Matched_Koala$marginals.fixed$`ba1:ci1`)
inla.pmarginal(0, Model_Matched_Koala$marginals.fixed$`time:ba1:ci1`)

print("Mixed Woody")
inla.pmarginal(0, Model_Mixed_Woody$marginals.fixed$`ba1:ci1`)
inla.pmarginal(0, Model_Mixed_Woody$marginals.fixed$`time:ba1:ci1`)

print("Mixed Koala")
inla.pmarginal(0, Model_Mixed_Koala$marginals.fixed$`ba1:ci1`)
inla.pmarginal(0, Model_Mixed_Koala$marginals.fixed$`time:ba1:ci1`)

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

# get posterior samples
# matched
Post_Matched_Woody <- inla.posterior.sample(10000, Model_Matched_Woody, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Matched_Woody_BACI <- lapply(Post_Matched_Woody, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Matched_Woody_BACITIME <- lapply(Post_Matched_Woody, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Matched_Koala <- inla.posterior.sample(10000, Model_Matched_Koala, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Matched_Koala_BACI <- lapply(Post_Matched_Koala, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Matched_Koala_BACITIME <- lapply(Post_Matched_Koala, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
# mixed
Post_Mixed_Woody <- inla.posterior.sample(10000, Model_Mixed_Woody, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Mixed_Woody_BACI <- lapply(Post_Mixed_Woody, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Mixed_Woody_BACITIME <- lapply(Post_Mixed_Woody, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Mixed_Koala <- inla.posterior.sample(10000, Model_Mixed_Koala, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Mixed_Koala_BACI <- lapply(Post_Mixed_Koala, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Mixed_Koala_BACITIME <- lapply(Post_Mixed_Koala, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()

# load property data
zonal_woody_treat <- read_rds("input/analysis_data/woody_treat_zonal.rds")
zonal_woody_contr  <- read_rds("input/analysis_data/woody_contr_zonal.rds")
zonal_koala_treat <- read_rds("input/analysis_data/koala_treat_zonal.rds")
zonal_koala_contr <- read_rds("input/analysis_data/koala_contr_zonal.rds")

# get the treatment woody vegetation and koala habitat for matched properties 
Woody_Treat_Matched <- zonal_woody_treat[which((zonal_woody_treat$baseline > 0 & zonal_woody_contr$baseline == 0) | (zonal_woody_treat$baseline == 0 & zonal_woody_contr$baseline > 0)), ] %>% filter(baseline > 0)
Koala_Treat_Matched <- zonal_koala_treat[which((zonal_koala_treat$baseline > 0 & zonal_koala_contr$baseline == 0) | (zonal_koala_treat$baseline == 0 & zonal_koala_contr$baseline > 0)), ] %>% filter(baseline > 0)

# get the treatment woody vegetation and koala habitat for mixed properties
Woody_Treat_Mixed <- zonal_woody_treat[which(zonal_woody_treat$baseline > 0 & zonal_woody_contr$baseline > 0), ] %>% filter(baseline > 0)  
Koala_Treat_Mixed <- zonal_koala_treat[which(zonal_koala_treat$baseline > 0 & zonal_koala_contr$baseline > 0), ] %>% filter(baseline > 0)

# get the relevant coeficient draws from the posteriors and calculate the impact of the treatment in hectares

# matched woody

# get actual proportions cleared
p <- Woody_Treat_Matched %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021) %>% mutate(p1415 = ifelse(p1415 == 1, 0.999999999999999, p1415)) %>% mutate(p1516 = ifelse(p1516 == 1, 0.999999999999999, p1516)) %>% mutate(p1617 = ifelse(p1617 == 1, 0.999999999999999, p1617)) %>% mutate(p1718 = ifelse(p1718 == 1, 0.999999999999999, p1718)) %>% mutate(p1819 = ifelse(p1819 == 1, 0.999999999999999, p1819)) %>% mutate(p1920 = ifelse(p1920 == 1, 0.999999999999999, p1920)) %>% mutate(p2021 = ifelse(p2021 == 1, 0.999999999999999, p2021)) %>% mutate(p1415 = ifelse(p1415 == 0, 0.000000000000001, p1415)) %>% mutate(p1516 = ifelse(p1516 == 0, 0.000000000000001, p1516)) %>% mutate(p1617 = ifelse(p1617 == 0, 0.000000000000001, p1617)) %>% mutate(p1718 = ifelse(p1718 == 0, 0.000000000000001, p1718)) %>% mutate(p1819 = ifelse(p1819 == 0, 0.000000000000001, p1819)) %>% mutate(p1920 = ifelse(p1920 == 0, 0.000000000000001, p1920)) %>% mutate(p2021 = ifelse(p2021 == 0, 0.000000000000001, p2021)) %>% as.matrix()

# get actual amounts cleared
c <- Woody_Treat_Matched[, 7:13] %>% as.matrix()

# get the baseline
b <- matrix(rep(Woody_Treat_Matched[, 2] %>% as.matrix(), 7), nrow = nrow(Woody_Treat_Matched), ncol = 7) %>% as.matrix()

# get time
t <- matrix(rep(c(0:6), nrow(p)), ncol = nrow(p), nrow = 7) %>% t()

# get the parameter draws
Params <- cbind(Post_Matched_Woody_BACI, Post_Matched_Woody_BACITIME) %>% as.matrix()

# get the predicted area effects 
Matched_Woody_Effect_Ha <- apply(Params, MARGIN = 1, FUN = function(x) {get_area_effect(p, c, b, t, x)})
Matched_Woody_Effect_Ha <- do.call("rbind", Matched_Woody_Effect_Ha)

# matched koala

# get actual proportions cleared
p <- Koala_Treat_Matched %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021) %>% mutate(p1415 = ifelse(p1415 == 1, 0.999999999999999, p1415)) %>% mutate(p1516 = ifelse(p1516 == 1, 0.999999999999999, p1516)) %>% mutate(p1617 = ifelse(p1617 == 1, 0.999999999999999, p1617)) %>% mutate(p1718 = ifelse(p1718 == 1, 0.999999999999999, p1718)) %>% mutate(p1819 = ifelse(p1819 == 1, 0.999999999999999, p1819)) %>% mutate(p1920 = ifelse(p1920 == 1, 0.999999999999999, p1920)) %>% mutate(p2021 = ifelse(p2021 == 1, 0.999999999999999, p2021)) %>% mutate(p1415 = ifelse(p1415 == 0, 0.000000000000001, p1415)) %>% mutate(p1516 = ifelse(p1516 == 0, 0.000000000000001, p1516)) %>% mutate(p1617 = ifelse(p1617 == 0, 0.000000000000001, p1617)) %>% mutate(p1718 = ifelse(p1718 == 0, 0.000000000000001, p1718)) %>% mutate(p1819 = ifelse(p1819 == 0, 0.000000000000001, p1819)) %>% mutate(p1920 = ifelse(p1920 == 0, 0.000000000000001, p1920)) %>% mutate(p2021 = ifelse(p2021 == 0, 0.000000000000001, p2021)) %>% as.matrix()

# get actual amounts cleared
c <- Koala_Treat_Matched[, 7:13] %>% as.matrix()

# get the baseline
b <- matrix(rep(Koala_Treat_Matched[, 2] %>% as.matrix(), 7), nrow = nrow(Koala_Treat_Matched), ncol = 7) %>% as.matrix()

# get time
t <- matrix(rep(c(0:6), nrow(p)), ncol = nrow(p), nrow = 7) %>% t()

# get the parameter draws
Params <- cbind(Post_Matched_Koala_BACI, Post_Matched_Koala_BACITIME) %>% as.matrix()

# get the predicted area effects 
Matched_Koala_Effect_Ha <- apply(Params, MARGIN = 1, FUN = function(x) {get_area_effect(p, c, b, t, x)})
Matched_Koala_Effect_Ha <- do.call("rbind", Matched_Koala_Effect_Ha)

# mixed woody

# get actual proportions cleared
p <- Woody_Treat_Mixed %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021) %>% mutate(p1415 = ifelse(p1415 == 1, 0.999999999999999, p1415)) %>% mutate(p1516 = ifelse(p1516 == 1, 0.999999999999999, p1516)) %>% mutate(p1617 = ifelse(p1617 == 1, 0.999999999999999, p1617)) %>% mutate(p1718 = ifelse(p1718 == 1, 0.999999999999999, p1718)) %>% mutate(p1819 = ifelse(p1819 == 1, 0.999999999999999, p1819)) %>% mutate(p1920 = ifelse(p1920 == 1, 0.999999999999999, p1920)) %>% mutate(p2021 = ifelse(p2021 == 1, 0.999999999999999, p2021)) %>% mutate(p1415 = ifelse(p1415 == 0, 0.000000000000001, p1415)) %>% mutate(p1516 = ifelse(p1516 == 0, 0.000000000000001, p1516)) %>% mutate(p1617 = ifelse(p1617 == 0, 0.000000000000001, p1617)) %>% mutate(p1718 = ifelse(p1718 == 0, 0.000000000000001, p1718)) %>% mutate(p1819 = ifelse(p1819 == 0, 0.000000000000001, p1819)) %>% mutate(p1920 = ifelse(p1920 == 0, 0.000000000000001, p1920)) %>% mutate(p2021 = ifelse(p2021 == 0, 0.000000000000001, p2021)) %>% as.matrix()

# get actual amounts cleared
c <- Woody_Treat_Mixed[, 7:13] %>% as.matrix()

# get the baseline
b <- matrix(rep(Woody_Treat_Mixed[, 2] %>% as.matrix(), 7), nrow = nrow(Woody_Treat_Mixed), ncol = 7) %>% as.matrix()

# get time
t <- matrix(rep(c(0:6), nrow(p)), ncol = nrow(p), nrow = 7) %>% t()

# get the parameter draws
Params <- cbind(Post_Mixed_Woody_BACI, Post_Mixed_Woody_BACITIME) %>% as.matrix()

# get the predicted area effects 
Mixed_Woody_Effect_Ha <- apply(Params, MARGIN = 1, FUN = function(x) {get_area_effect(p, c, b, t, x)})
Mixed_Woody_Effect_Ha <- do.call("rbind", Mixed_Woody_Effect_Ha)

# mixed koala

# get actual proportions cleared
p <- Koala_Treat_Mixed %>% mutate (p1415 = clear1415 / baseline, p1516 = clear1516 / baseline, p1617 = clear1617 / baseline, p1718 = clear1718 / baseline, p1819 = clear1819 / baseline, p1920 = clear1920 / baseline, p2021 = clear2021 / baseline) %>% select(p1415, p1516, p1617, p1718, p1819, p1920, p2021) %>% mutate(p1415 = ifelse(p1415 == 1, 0.999999999999999, p1415)) %>% mutate(p1516 = ifelse(p1516 == 1, 0.999999999999999, p1516)) %>% mutate(p1617 = ifelse(p1617 == 1, 0.999999999999999, p1617)) %>% mutate(p1718 = ifelse(p1718 == 1, 0.999999999999999, p1718)) %>% mutate(p1819 = ifelse(p1819 == 1, 0.999999999999999, p1819)) %>% mutate(p1920 = ifelse(p1920 == 1, 0.999999999999999, p1920)) %>% mutate(p2021 = ifelse(p2021 == 1, 0.999999999999999, p2021)) %>% mutate(p1415 = ifelse(p1415 == 0, 0.000000000000001, p1415)) %>% mutate(p1516 = ifelse(p1516 == 0, 0.000000000000001, p1516)) %>% mutate(p1617 = ifelse(p1617 == 0, 0.000000000000001, p1617)) %>% mutate(p1718 = ifelse(p1718 == 0, 0.000000000000001, p1718)) %>% mutate(p1819 = ifelse(p1819 == 0, 0.000000000000001, p1819)) %>% mutate(p1920 = ifelse(p1920 == 0, 0.000000000000001, p1920)) %>% mutate(p2021 = ifelse(p2021 == 0, 0.000000000000001, p2021)) %>% as.matrix()

# get actual amounts cleared
c <- Koala_Treat_Mixed[, 7:13] %>% as.matrix()

# get the baseline
b <- matrix(rep(Koala_Treat_Mixed[, 2] %>% as.matrix(), 7), nrow = nrow(Koala_Treat_Mixed), ncol = 7) %>% as.matrix()

# get time
t <- matrix(rep(c(0:6), nrow(p)), ncol = nrow(p), nrow = 7) %>% t()

# get the parameter draws
Params <- cbind(Post_Mixed_Koala_BACI, Post_Mixed_Koala_BACITIME) %>% as.matrix()

# get the predicted area effects 
Mixed_Koala_Effect_Ha <- apply(Params, MARGIN = 1, FUN = function(x) {get_area_effect(p, c, b, t, x)})
Mixed_Koala_Effect_Ha <- do.call("rbind", Mixed_Koala_Effect_Ha)

# get the total amounts cleared
Total_Woody <- Mixed_Woody_Effect_Ha + Matched_Woody_Effect_Ha
Total_Koala <- Mixed_Koala_Effect_Ha + Matched_Koala_Effect_Ha

# get median, mean, and upper and lower limits
Total_Woody_Mean <- apply(Total_Woody, MARGIN = 2, FUN = mean)
Total_Woody_Median <- apply(Total_Woody, MARGIN = 2, FUN = median)
Total_Woody_Lower <- apply(Total_Woody, MARGIN = 2, FUN = function(x) {quantile(x, probs = 0.025)})
Total_Woody_Upper <- apply(Total_Woody, MARGIN = 2, FUN = function(x) {quantile(x, probs = 0.975)})
Total_Koala_Mean <- apply(Total_Koala, MARGIN = 2, FUN = mean)
Total_Koala_Median <- apply(Total_Koala, MARGIN = 2, FUN = median)
Total_Koala_Lower <- apply(Total_Koala, MARGIN = 2, FUN = function(x) {quantile(x, probs = 0.025)})
Total_Koala_Upper <- apply(Total_Koala, MARGIN = 2, FUN = function(x) {quantile(x, probs = 0.975)})
Total_Woody_Mean_Sum <- sum(Total_Woody_Mean)
Total_Woody_Median_Sum <- sum(Total_Woody_Median)
Total_Woody_Lower_Sum <- sum(Total_Woody_Lower)
Total_Woody_Upper_Sum <- sum(Total_Woody_Upper)
Total_Koala_Mean_Sum <- sum(Total_Koala_Mean)
Total_Koala_Median_Sum <- sum(Total_Koala_Median)
Total_Koala_Lower_Sum <- sum(Total_Koala_Lower)
Total_Koala_Upper_Sum <- sum(Total_Koala_Upper)

# get the total proportion > 0
Total_Woody_SumYear <- apply(Total_Woody, MARGIN = 1, sum)
Total_Koala_SumYear <- apply(Total_Koala, MARGIN = 1, sum)

# get the total proportion > 0
Total_Woody_SumYearProp <- length(which(Total_Woody_SumYear > 0)) / length(Total_Woody_SumYear)
Total_Koala_SumYearProp <- length(which(Total_Koala_SumYear > 0)) / length(Total_Koala_SumYear)

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
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Trend Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-10, 20), breaks = seq(-10, 20, by = 5))

# combined plot

Plot <- ggarrange(PlotA, PlotB, ncol = 1, nrow = 2, vjust = -1, hjust = 0, font.label = list(size = 18))
ggsave(Plot, file = "output/figures/effects_matched.jpg", width = 20, height = 30, units = "cm", dpi = 300)

# mixed properties

# immediate effect

Coeff_EffectsA <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["ba1:ci1", "mean"], Summary_Model_Mixed_Koala["ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["ba1:ci1", "0.975quant"])))
Coeff_EffectsA <- Coeff_EffectsA %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotA <- ggplot(Coeff_EffectsA, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Immediate Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-0.2, 0.9), breaks = seq(-0.2, 0.9, by = 0.1))

# trend effect

Coeff_EffectsB <- as_tibble(cbind(Type = c("Woody Vegetation", "Koala Habitat"), Est = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "mean"], Summary_Model_Mixed_Koala["time:ba1:ci1", "mean"]), Lower = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.025quant"], Summary_Model_Mixed_Koala["time:ba1:ci1", "0.025quant"]), Upper = c(Summary_Model_Mixed_Woody["time:ba1:ci1", "0.975quant"], Summary_Model_Mixed_Koala["time:ba1:ci1", "0.975quant"])))
Coeff_EffectsB <- Coeff_EffectsB %>% mutate(Est = as.numeric(Est), Lower = as.numeric(Lower), Upper = as.numeric(Upper)) %>% mutate(Type = as.factor(Type)) %>% mutate(Type = relevel(Type, "Woody Vegetation")) %>% group_by(Type)

PlotB <- ggplot(Coeff_EffectsB, aes(x = Type, y = Est, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") + geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.7), width = 0.25) + labs(x = "Vegetation Type", y = "Effect Size") + theme_minimal() + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1)) + geom_hline(yintercept = 0, linetype = "solid", color = "black") + theme(plot.margin = margin(b = 0.5, unit = "cm")) + ggtitle("Trend Effect") + theme(plot.title = element_text(size = 22)) + scale_y_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, by = 0.1))

# combined plot

Plot <- ggarrange(PlotA, PlotB, ncol = 1, nrow = 2, vjust = -1, hjust = 0, font.label = list(size = 18))
ggsave(Plot, file = "output/figures/effects_mixed.jpg", width = 20, height = 30, units = "cm", dpi = 300)

# plots of trends

# woody vegetation

PlotData <- pivot_longer(Total_Woody, cols = 1:7, names_to = "Year", values_to = "Area") %>% mutate(Year = paste0("20", str_sub(Year, 8, 9))) %>% mutate(Year = as.factor(Year))
Plot <- ggplot(PlotData, aes(x = Year, y = Area)) + geom_violin(, fill="lightblue") + theme_minimal() + scale_y_continuous(limits = c(-20000, 20000), breaks = seq(-20000, 20000, by = 5000)) + stat_summary(fun = median, geom = "point", size=2, color="red") + geom_hline(yintercept = 0) + labs(x = "Year", y = "Impact on Area Cleared (ha)") + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1))

ggsave(Plot, file = "output/figures/trend_woody.jpg", width = 25, height = 20, units = "cm", dpi = 300)

# koala habitat

PlotData <- pivot_longer(Total_Koala, cols = 1:7, names_to = "Year", values_to = "Area") %>% mutate(Year = paste0("20", str_sub(Year, 8, 9))) %>% mutate(Year = as.factor(Year))
Plot <- ggplot(PlotData, aes(x = Year, y = Area)) + geom_violin(, fill="lightblue") + theme_minimal() + scale_y_continuous(limits = c(-10000, 10000), breaks = seq(-20000, 20000, by = 5000)) + stat_summary(fun = median, geom = "point", size=2, color="red") + geom_hline(yintercept = 0) + labs(x = "Year", y = "Impact on Area Cleared (ha)") + theme(legend.position = "none", axis.text = element_text(size = 18),  axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20, vjust = -1))

ggsave(Plot, file = "output/figures/trend_koala.jpg", width = 25, height = 20, units = "cm", dpi = 300)
