randomise <- function(Data, N, Reps) {

  Out <- list()
  for (i in 1:Reps) {
    Out[[i]] <- Data %>% nest_by(RID) %>% ungroup() %>% slice_sample(n = N) %>% unnest(cols = c(data))
  }

  return(Out)
}

# fit an inla models - for use in fitting models in parallel
fit_inla <- function(X) {
  # fit model
  # binomial generalised linear model with counts of cells cleared as the dependent variable,
  # a random-effect for property ID to control dependence in data points within properties
  # and over-dispersion
  # this is a LOSS ~ TIME + BA + CI + (BA * CI) + (BA * TIME) + (CI * TIME) + (BA * CI * TIME) model
  Model <- inla(LOSS_count ~ TIME + BA + CI + BA:CI + BA:TIME + CI:TIME + BA:CI:TIME + f(RID, model = "iid"), data = X, family = "binomial", Ntrials = FOREST_cou, control.inla = list(control.vb = list(enable = FALSE)))

  return(Model)
}