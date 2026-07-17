## 1. Libraries ---------------------------------------------------------------
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(tidyverse)
library(writexl)

## 2. API base URL -------------------------------------------------------------
base_url <- "XXXXXX"

# COVID-19 vaccination RCTs
params2 <- list(
  "query.term" = "COVID-19 OR SARS-CoV-2 OR coronavirus disease 2019 OR COVID-19 vaccine OR SARS-CoV-2 vaccine" ,
  "query.cond" = "COVID-19",
  "pageSize" = 100)

raw_studies <- get_all_pages(base_url, params2)

safe_paste <- function(x, sep = "; ") {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  paste(x, collapse = sep)
}

trials_rct_completed <- raw_studies %>%
  filter(
    # Completed
    protocolSection.statusModule.overallStatus == "COMPLETED",
    
    # Interventional
    protocolSection.designModule.studyType == "INTERVENTIONAL",
    
    # Randomized
    protocolSection.designModule.designInfo.allocation == "RANDOMIZED"
  )

trials_rct_completed_p3 <- trials_rct_completed %>%
  filter(
    map_lgl(
      protocolSection.designModule.phases,
      ~ "PHASE3" %in% .x # Phase three
    )
  )

trials_rct_completed_p3 <- trials_rct_completed_p3 %>%
  filter(
    map_int(
      protocolSection.armsInterventionsModule.armGroups, # number of arms more than 2
      length
    ) >= 2
  )


## 7. Extract trial-level table -----------------------------------------------
trials <- trials_rct_completed_p3 %>%
  transmute(
    nct_id = protocolSection.identificationModule.nctId,
    title  = protocolSection.identificationModule.briefTitle,
    overall_status = protocolSection.statusModule.overallStatus,
    phase = map_chr(
      protocolSection.designModule.phases,
      safe_paste
    ),
    condition = map_chr(
      protocolSection.conditionsModule.conditions,
      safe_paste
    ),
    intervention = map_chr(
      protocolSection.armsInterventionsModule.interventions,
      ~ safe_paste(.x$name)
    ),
    start_date =
      protocolSection.statusModule.startDateStruct.date,
    completion_date =
      protocolSection.statusModule.completionDateStruct.date
  )

clean_trials <- trials_rct_completed_p3 %>% 
  select(
    protocolSection.identificationModule.nctId,
    protocolSection.identificationModule.briefTitle,
    protocolSection.identificationModule.officialTitle,
    protocolSection.sponsorCollaboratorsModule.leadSponsor.name,
    protocolSection.descriptionModule.briefSummary,
    protocolSection.descriptionModule.detailedDescription,
    resultsSection.outcomeMeasuresModule.outcomeMeasures,
    protocolSection.designModule.studyType,
    protocolSection.designModule.phases,
    resultsSection.adverseEventsModule.eventGroups
  ) %>% 
  rename(
    nct_id = protocolSection.identificationModule.nctId,
    brief_title = protocolSection.identificationModule.briefTitle,
    official_title = protocolSection.identificationModule.officialTitle,
    lead_sponsor = protocolSection.sponsorCollaboratorsModule.leadSponsor.name,
    brief_summary = protocolSection.descriptionModule.briefSummary,
    detailed_description = protocolSection.descriptionModule.detailedDescription,
    outcome_measures = resultsSection.outcomeMeasuresModule.outcomeMeasures,
    study_type = protocolSection.designModule.studyType,
    phases = protocolSection.designModule.phases,
    adverse_events = resultsSection.adverseEventsModule.eventGroups
  ) %>% 
  filter(!sapply(outcome_measures, is.null))

df_outcome_measure <- clean_trials %>%
  select(-adverse_events) %>% 
  unnest(outcome_measures) %>% 
  select(-paramType,-dispersionType) %>% 
  unnest(analyses)

aggre_measure_per_aim <- clean_trials %>%
  select(-adverse_events) %>% 
  unnest(outcome_measures) %>% 
  filter(!sapply(analyses, is.null))


df_adverse_event <- clean_trials %>%
  select(-outcome_measures)  %>%
  unnest(adverse_events)

df_outcome_measure_both <- df_outcome_measure %>%
  group_by(nct_id,brief_title,title,groupIds,nonInferiorityType,groupDescription) %>% 
  mutate(
    n_measures = n_distinct(paramType)) %>% 
  ungroup() %>% 
  group_by(nct_id) %>%
  mutate(
    n_aims = dense_rank(paste(title,groupIds,nonInferiorityType,groupDescription))) %>% 
  ungroup() %>% 
  arrange(nct_id,n_aims) %>% 
  mutate(groupIds = sapply(groupIds, function(x) paste(x, collapse = ", "))) %>%
  mutate(phases = sapply(phases, function(x) paste(x, collapse = ", ")))

# write to a excel file with two sheets
write_xlsx(list(
  outcome_measures = df_outcome_measure,
  adverse_events = df_adverse_event
), path = "covid19_rct_outcome_measures_adverse_events.xlsx")

write_xlsx(df_outcome_measure_both, path = "COVID_rct_outcome_measures_multiple.xlsx")




