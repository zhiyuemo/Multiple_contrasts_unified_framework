# Multiple contrasts unified framework: Replication materials

This repository contains replication materials for our manuscript, A Unified Doubly Robust Framework for Causal Contrasts, with Evidence on Randomized Trials of COVID-19 Vaccines.

Scripts can be opened and run in RStudio. Source code is available in the directory src/. The main scripts are:

- `1_clinicalgovAPI.R`
  - **Overview**: This script queries the ClinicalTrials.gov API for COVID-19 trials, filters to completed randomized Phase 3 RCTs, and exports their outcome measures and adverse events to Excel.
  - **NOTE**: You will need to sign up for a ClinicalTrials.gov API key to obtain the trials information. Please visit: https://clinicaltrials.gov/data-api/api for more information.
  - **Outputs**:
    - `COVID_rct_outcome_measures_multiple.xlsx`
    - `covid19_rct_outcome_measures_adverse_events.xlsx`

- `2_simulation.R`
  - **Overview**: This script generates the simulation results of a single experiment. We uploaded this to the OSG platform to run thousands times of simulations.
  - **NOTE**: If you would like to learn more about OSG, here is a link to the documentation for using OSG: https://portal.osg-htc.org/documentation/. To get started, please refer to the steps in https://portal.osg-htc.org/documentation/overview/account_setup/registration-and-login/
  - **Outputs**:
    - `simulation_results_1.csv`

- `3_multiContrast_moderna.R`
  - **Overview**: This script generates the results of Moderna COVID-19 data application.
  - **NOTE**: Moderna data is not publicly available.
  - **Input**:
    - dataframe `moderna_raw`
  - **Outputs**:
    - dataframe `result`

    
