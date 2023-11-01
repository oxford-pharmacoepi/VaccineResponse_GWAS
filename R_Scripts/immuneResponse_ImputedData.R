# ============================================================================ #
#                               IMMUNE RESPONSE                                #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
# Cohort -----------------------------------------------------------------------
ukb_covid <- as_tibble(read.delim(paste0(dir_data,"immune_response.tab"))) %>% 
  select("eid" = "f.eid",
         "Antibody_test_result" = "f.27981.0.0",
         "Date_antibody_test_performed" = "f.27982.0.0",
         "Received_first_COVID19_vaccination" = "f.27983.0.0",
         "Date_first_COVID19_vaccination" = "f.27984.0.0",
         "Received_second_COVID19_vaccination" = "f.27985.0.0",
         "Date_second_COVID19_vaccination" = "f.27986.0.0",
         "Covid_antibody_test_results" = "f.27990.0.0") %>%
  mutate(Antibody_test_result = if_else(Antibody_test_result == 2,1,Antibody_test_result)) %>%
  mutate(Date_first_COVID19_vaccination  = as.Date(Date_first_COVID19_vaccination,"%Y-%m-%d"),
         Year_first_COVID19_vaccination  = year(Date_first_COVID19_vaccination),
         Date_second_COVID19_vaccination = as.Date(Date_second_COVID19_vaccination,"%Y-%m-%d"),
         Year_second_COVID19_vaccination = year(Date_second_COVID19_vaccination),
         Date_antibody_test_performed    = as.Date(Date_antibody_test_performed,"%Y-%m-%d"),
         Year_antibody_test_performed    = year(Date_antibody_test_performed))

ukb_covar <- as_tibble(read_delim(paste0(dir_data,'ukb65397_covariates.csv')))

# General cohort
vaccine_cohort <- ukb_covid %>%
  filter(!is.na(Antibody_test_result)) %>%
  left_join(ukb_covar, by = "eid") %>%
  filter(Caucasian == 1) %>%
  filter(Sex == genetic_sex) %>%
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  filter(kinship_to_other_participants == 0) %>%
  filter(Received_first_COVID19_vaccination == 1) %>% 
  filter(Year_first_COVID19_vaccination >= 2021)


# ONE VACCINE COHORT -----------------------------------------------------------
one_vaccine_cohort <- vaccine_cohort %>% 
  filter(Received_second_COVID19_vaccination == 0) %>%
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_first_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 84) %>%
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) %>%
  filter(Year_antibody_test_performed >= 2021) %>%
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID') %>%
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
write.csv(one_vaccine_cohort,paste0(dir_results,'Cohorts/imputedData_oneDose.csv'))

# TWO VACCINE COHORT -----------------------------------------------------------
two_vaccine_cohort <- vaccine_cohort %>%
  filter(Received_second_COVID19_vaccination == 1) %>%
  filter(Year_second_COVID19_vaccination >= 2021) %>%
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_second_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 84) %>%
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) %>%
  filter(Year_antibody_test_performed >= 2021) %>%
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID') %>%
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
write.csv(two_vaccine_cohort,paste0(dir_results,'Cohorts/imputedData_twoDose.csv'))

# .PHE FILES
one_dose_cohort <- as.phe(one_vaccine_cohort, "FID", "IID") %>% 
  rename('FID' = 'pid',
         'IID' = 'id')
two_dose_cohort <- as.phe(two_vaccine_cohort, "FID", "IID") %>% 
  rename('FID' = 'pid',
         'IID' = 'id')

write.phe(paste0(dir_results,'/Cohorts/imputedData_oneDose.phe'), one_dose_cohort)
write.phe(paste0(dir_results,'/Cohorts/imputedData_twoDose.phe'), two_dose_cohort)


# Validation -------------------------------------------------------------------
# General cohort
vaccine_cohort <- ukb_covid %>%
  filter(Received_first_COVID19_vaccination == 1) %>% # GENERAL FILTERS
  filter(Year_first_COVID19_vaccination >= 2021) %>%
  left_join(ukb_covar,
            by = 'eid') %>% # QUALITY CONTROL
  filter(!is.na(PC1), !is.na(PC2), !is.na(PC3), !is.na(PC4), !is.na(PC5),
         !is.na(PC6), !is.na(PC7), !is.na(PC8), !is.na(PC9), !is.na(PC10)) %>%
  filter(is.na(Caucasian)) %>%
  filter(Sex == genetic_sex) %>%
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  filter(kinship_to_other_participants == 0)

# ONE VACCINE COHORT -----------------------------------------------------------
one_vaccine_cohort <- vaccine_cohort %>% 
  filter(Received_second_COVID19_vaccination == 0) %>%
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_first_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 84) %>%
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) %>%
  filter(Year_antibody_test_performed >= 2021) %>%
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID') %>%
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# TWO VACCINE COHORT -----------------------------------------------------------
two_vaccine_cohort <- vaccine_cohort %>%
  filter(Received_second_COVID19_vaccination == 1) %>%
  filter(Year_second_COVID19_vaccination >= 2021) %>%
  mutate(Days_since_vaccine = as.numeric(difftime(Date_antibody_test_performed,Date_second_COVID19_vaccination, units = "days"))) %>% 
  filter(Days_since_vaccine <= 84) %>%
  filter(Covid_antibody_test_results == 0 | is.na(Covid_antibody_test_results)) %>%
  filter(Year_antibody_test_performed >= 2021) %>%
  mutate(Age_antibody_test_performed = Year_antibody_test_performed - Year_of_birth) %>%
  mutate('FID' = eid, 'IID' = eid) %>% relocate('FID', 'IID') %>%
  select(FID, IID, 
         immuneResponse = Antibody_test_result, 
         Sex, 
         Age = Age_antibody_test_performed, 
         Genetic_batch, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# .PHE FILES
one_dose_cohort <- as.phe(one_vaccine_cohort, "FID", "IID") %>% 
  rename('FID' = 'pid',
         'IID' = 'id')
two_dose_cohort <- as.phe(two_vaccine_cohort, "FID", "IID") %>% 
  rename('FID' = 'pid',
         'IID' = 'id')

write.phe(paste0(dir_results,'/Cohorts/imputedData_oneDose_Validation.phe'), one_dose_cohort)
write.phe(paste0(dir_results,'/Cohorts/imputedData_twoDose_Validation.phe.phe'), two_dose_cohort)
