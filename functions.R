get_area_effect <- function(p, c, b, t, Params) {
  # get area cleared
  Out <- ((c - ((exp(log(p / (1 - p)) - Params[1] - Params[2] * t) / (1 + exp(log(p / (1 - p)) - Params[1] - Params[2] * t)))) * b) * 25 * 25 / 10000) %>% as_tibble()

  OutSums <- Out %>% summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))

  return(OutSums)
}