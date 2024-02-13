# ============================================================================ #
#                           BREAKTHROUGH IMPUTED DATA                          #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
# Cohort -----------------------------------------------------------------------
print("load databases")
ukb_covid <- readRDS(here(paste0(dir_data,'breakthrough.rds'))) %>%
  select('eid','event_dt','case_control_index','breakthrough_infection_type','severity_index' = 'origin')  %>%
  mutate('breakthrough_infection_index' = if_else(case_control_index == 'Control',0,1))

ukb_covar <- as_tibble(read_delim(paste0(dir_data,'ukb65397_covariates.csv')))

# Attrition --------------------------------------------------------------------
print("load attrition")
attr(ukb_covid, "cohort_attrition") <- tibble("N" = ukb_covid |> nrow(), 
                                              "Reason"   = "UK Biobank participants",
                                              "Excluded" = 0)
# Breakthrough infection, general covid ----------------------------------------
print("load breakthrough cohort")
breakthrough <- ukb_covid %>%
  recordAttrition("participants have received at least one COVID-19 vaccine") |>
  left_join(ukb_covar, by = 'eid') %>%
  filter(Caucasian == 1) %>%
  recordAttrition(reason = "participants with European ancestry") |>
  filter(Sex == genetic_sex) %>%
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  mutate(Age = year(event_dt) - Year_of_birth) %>%
  select(eid, bt_infection = breakthrough_infection_index, breakthrough_infection_type, severity_index, Sex, Age, Genetic_batch, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# Breakthrough susceptibility --------------------------------------------------
print("load breakthrough susceptibility")
breakthrough_susceptibility <- breakthrough |>
  select(eid, bt_infection, Sex, Age, Genetic_batch, starts_with("PC")) |>
  mutate("IID" = eid, "FID" = eid) |>
  relocate('FID', 'IID') 

attr(breakthrough_susceptibility,"cohort_attrition") <- attr(breakthrough_susceptibility,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_susceptibility |> filter(bt_infection == 1) |> tally() |> pull(),
                           breakthrough_susceptibility |> filter(bt_infection == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

# breakthrough severity --------------------------------------------------------
print("load breakthrough severity")
breakthrough_severity <- breakthrough %>%
  filter(bt_infection == 1) %>%
  rename(severity = severity_index) %>%
  recordAttrition(reason = "participants infected") |>
  select(-bt_infection, -breakthrough_infection_type) |>
  mutate("IID" = eid, "FID" = eid) |>
  relocate('FID', 'IID') 

attr(breakthrough_severity,"cohort_attrition") <- attr(breakthrough_severity,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_severity |> filter(severity == 1) |> tally() |> pull(),
                           breakthrough_severity |> filter(severity == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(breakthrough_susceptibility, paste0(dir_results,'/Cohorts/breakthrough_susceptibility.txt'), row.names = FALSE)
write.table(attr(breakthrough_susceptibility,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_susceptibility_attrition.txt'), row.names = FALSE)

write.table(breakthrough_severity, paste0(dir_results,'/Cohorts/breakthrough_severity.txt'), row.names = FALSE)
write.table(attr(breakthrough_severity,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_severity_attrition.txt'), row.names = FALSE)

# write .phe files
print("load .phe files")
breakthrough_susceptibility <- as.phe(breakthrough_susceptibility, "FID", "IID")
breakthrough_severity       <- as.phe(breakthrough_severity, "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/breakthrough_susceptibility.phe'), breakthrough_susceptibility)
write.phe(paste0(dir_results,'/Cohorts/breakthrough_severity.phe'), breakthrough_severity)


# Validation -------------------------------------------------------------------
# Breakthrough infection, general covid ----------------------------------------
print("load breakthrough validation cohort")
breakthrough <- ukb_covid %>%
  recordAttrition("participants have received at least one COVID-19 vaccine") |>
  left_join(ukb_covar, by = 'eid') %>%
  filter(is.na(Caucasian)) %>%
  recordAttrition(reason = "participants without European ancestry") |>
  filter(Sex == genetic_sex) %>%
  recordAttrition(reason = "participants with the same sex and genetic sex registered") |>
  filter(is.na(sex_chromosome_aneuploidy)) %>%
  recordAttrition(reason = "participants with no sex chromosome aneuploidy") |>
  mutate(Age = year(event_dt) - Year_of_birth) %>%
  select(eid, bt_infection = breakthrough_infection_index, breakthrough_infection_type, severity_index,Sex, Age, Genetic_batch, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

# Breakthrough susceptibility --------------------------------------------------
print("load breakthrough susceptibility")
breakthrough_susceptibility <- breakthrough |>
  select(eid, bt_infection, Sex, Age, Genetic_batch, starts_with("PC")) |>
  mutate("IID" = eid, "FID" = eid) |>
  relocate('FID', 'IID') 

attr(breakthrough_susceptibility,"cohort_attrition") <- attr(breakthrough_susceptibility,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_susceptibility |> filter(bt_infection == 1) |> tally() |> pull(),
                           breakthrough_susceptibility |> filter(bt_infection == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

# breakthrough severity --------------------------------------------------------
print("load breakthrough severity")
breakthrough_severity <- breakthrough %>%
  filter(bt_infection == 1) %>%
  rename(severity = severity_index) %>%
  recordAttrition(reason = "participants infected") |>
  select( -bt_infection, -breakthrough_infection_type) |>
  mutate("IID" = eid, "FID" = eid) |>
  relocate('FID', 'IID') 

attr(breakthrough_severity,"cohort_attrition") <- attr(breakthrough_severity,"cohort_attrition") |>
  union_all(tibble("N" = c(breakthrough_severity |> filter(severity == 1) |> tally() |> pull(),
                           breakthrough_severity |> filter(severity == 0) |> tally() |> pull()),
                   "Reason" = c("cases","controls"),
                   "Excluded" = 0))

write.table(breakthrough_susceptibility, paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation.txt'), row.names = FALSE)
write.table(attr(breakthrough_susceptibility,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation_attrition.txt'), row.names = FALSE)

write.table(breakthrough_severity, paste0(dir_results,'/Cohorts/breakthrough_severity_validation.txt'), row.names = FALSE)
write.table(attr(breakthrough_severity,"cohort_attrition"), paste0(dir_results,'/Cohorts/breakthrough_severity_validation_attrition.txt'), row.names = FALSE)


# write .phe files
print("load .phe files")
breakthrough_susceptibility <- as.phe(breakthrough_susceptibility |> select(-eid), "FID", "IID")
breakthrough_severity       <- as.phe(breakthrough_severity |> select(-eid), "FID", "IID")

write.phe(paste0(dir_results,'/Cohorts/breakthrough_susceptibility_validation.phe'), breakthrough_susceptibility)
write.phe(paste0(dir_results,'/Cohorts/breakthrough_severity_validation.phe'), breakthrough_severity)
