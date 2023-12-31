# ============================================================================ #
#                           BREAKTHROUGH IMPUTED DATA                          #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #

ukb_covid <- readRDS(here(paste0(dir_data,'breakthrough.rds'))) %>%
  select('eid','event_dt','case_control_index','breakthrough_infection_type','severity_index' = 'origin')  %>%
  mutate('breakthrough_infection_index' = if_else(case_control_index == 'Control',0,1))

ukb_covar <- as_tibble(read_delim(paste0(dir_data,'ukb65397_covariates.csv')))

# Breakthrough infection, general covid ----------------------------------------
breakthrough <- ukb_covid %>%
  left_join(ukb_covar, by = 'eid') %>%
  filter(Caucasian == 1) %>% #  filter(is.na(Caucasian))
  filter(Sex == genetic_sex) %>%
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  filter(kinship_to_other_participants == 0) %>%
  mutate(Age = year(event_dt) - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID','IID') %>%
  select(FID,IID, bt_infection = breakthrough_infection_index, breakthrough_infection_type, severity_index,Sex, Age, Genetic_batch, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# breakthrough severity, general covid -----------------------------------------
breakthrough_severity <- breakthrough %>%
  filter(bt_infection == 1) %>%
  rename(severity = severity_index) %>%
  select(-bt_infection, -breakthrough_infection_type)

breakthrough <- breakthrough %>%
  select(-breakthrough_infection_type, -severity_index)

write.csv(breakthrough, paste0(dir_results,'/Cohorts/imputedData_breakthroughSusceptibility.csv'))
write.csv(breakthrough_severity, paste0(dir_results,'/Cohorts/imputedData_breakthroughSeverity.csv'))

write.phe(paste0(dir_results,'/Cohorts/imputedData_breakthroughSusceptibility.phe'), as.phe(breakthrough, pid = 'FID', id = 'IID'))
write.phe(paste0(dir_results,'/Cohorts/imputedData_breakthroughSeverity.phe'), as.phe(breakthrough_severity, pid = 'FID', id = 'IID'))


# Validation -------------------------------------------------------------------
breakthrough <- ukb_covid %>%
  left_join(ukb_covar, by = 'eid') %>%
  filter(!is.na(PC1), !is.na(PC2), !is.na(PC3), !is.na(PC4), !is.na(PC5),
         !is.na(PC6), !is.na(PC7), !is.na(PC8), !is.na(PC9), !is.na(PC10)) %>%
  filter(is.na(Caucasian)) %>% #  filter(is.na(Caucasian))
  filter(Sex == genetic_sex) %>%
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  filter(kinship_to_other_participants == 0) %>%
  mutate(Age = year(event_dt) - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID','IID') %>%
  select(FID,IID, bt_infection = breakthrough_infection_index, breakthrough_infection_type, severity_index,Sex, Age, Genetic_batch, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# breakthrough severity, covid -----------------------------------------
breakthrough_severity <- breakthrough %>%
  filter(bt_infection == 1) %>%
  rename(severity = severity_index) %>%
  select(-bt_infection, -breakthrough_infection_type)


breakthrough <- breakthrough %>%
  select(-breakthrough_infection_type, -severity_index)

write.phe(paste0(dir_results,'/Cohorts/imputedData_breakthroughSusceptibility_validation.phe'), as.phe(breakthrough, pid = 'FID', id = 'IID'))
write.phe(paste0(dir_results,'/Cohorts/imputedData_breakthroughSeverity_validation.phe'), as.phe(breakthrough_severity, pid = 'FID', id = 'IID'))

