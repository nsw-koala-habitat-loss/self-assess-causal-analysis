# get required libraries
library(foreign)
library(tidyverse)
library(INLA)

# fit models to all woody vegetation data for separated properties

# load data
Data1 <- as_tibble(read.dbf("input/BACItime_v2.dbf"))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model1 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data1, family = "binomial", Ntrials = FOREST_cou)

# get summary and statistical significance
summary(Model1)

# save the model object as an RDS file
saveRDS(Model1, file = "output/Model1.rds")

# fit models to koala habitat for separated properties

# load data
Data2 <- as_tibble(read.dbf("input/BACItime_khab_match_oldtime.dbf"))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model2 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data2, family = "binomial", Ntrials = FOREST_cou)

# get summary and statistical significance
summary(Model2)

# save the model object as an RDS file
saveRDS(Model2, file = "output/Model2.rds")

# fit models to all woody vegetation data for mixed properties

# load data
Data3 <- as_tibble(read.dbf("input/MIXED_BACItime.dbf")) %>% mutate(CI = TREAT_inv)
Data3 <- Data3 %>% mutate(FOREST_cou = ifelse(FOREST_cou < LOSS_count, LOSS_count, FOREST_cou))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model3 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data3, family = "binomial", Ntrials = FOREST_cou, control.inla = list(control.vb = list(enable = FALSE)))

Model3 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME, data = Data3, family = "binomial", Ntrials = FOREST_cou, control.inla = list(control.vb = list(enable = FALSE)))

# get summary and statistical significance
summary(Model3)

# save the model object as an RDS file
saveRDS(Model3, file = "output/Model3.rds")

# fit models to koala habitat for mixed properties

# load data
Data4 <- as_tibble(read.dbf("input/BACItime_KHAB_MIXED_OLDTIME.dbf"))

# fit model
# binomial generalised linear model with counts of cells cleared as the dependent variable,
# a random-effect for property ID to control dependence in data points within properties
# and over-dispersion
# this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
Model4 <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = Data4, family = "binomial", Ntrials = FOREST_cou)

# get summary and statistical significance
summary(Model4)

# save the model object as an RDS file
saveRDS(Model4, file = "output/Model4.rds")
