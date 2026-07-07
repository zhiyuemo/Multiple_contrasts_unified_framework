# install.packages("remotes")
#remotes::install_github("zhiyuemo/MultiContrasts")

library(MultiContrasts)
library(ranger)
library(nnet)
library(tidyverse)

set.seed(2026)

available_sl_libraries()
moderna_bseronegative <- moderna_raw %>%
  # Include only Ptid with Perprotocol == 1
  filter(Perprotocol == 1) %>%
  # Include only Ptid with SubcohortInd == 1
  # filter(SubcohortInd == 1) %>% only use for antibody analysis
  # Include only Ptid with ph1.immuno == 1
  filter(ph1.immuno == 1) %>%
  # Filter for Bserostatus == 0
  filter(Bserostatus == 0) %>%
  # Select relevant columns and rename as needed
  select(
    X.1 = HighRiskInd,
    X.2 = Bstratum,
    X.3 = standardized_risk_score,
    X.4 = MinorityInd,
    A = Trt,
    Y = EventIndPrimaryD57
    # Y = Day57pseudoneutid50
  ) %>%
  # Remove rows with any NA values
  na.omit() 

# Check overlap between treatment groups
plot_covariate_diagnostics(
  data = moderna_bseronegative,
  treatment = "A",
  covariates = c("X.1", "X.2","X.3","X.4"),
  labels = c("Vaccine","Placebo")
)

# Estimate all 7 contrasts in one call
result <- aipw_estimate(
  data       = moderna_bseronegative,
  outcome    = "Y",
  treatment  = "A",
  n_folds = 5,
  covariates = c("X.1", "X.2","X.3","X.4"),
  sl_library_PS = c("SL.nnet","SL.glm"),
  sl_library_OR = c("SL.nnet","SL.glm")
)
result

# Plot all measures
plot_estimates(result)

# Plot specific measures and save
plot_estimates(result, measures = c("RD", "NNT"), save_dir = "figures/")
