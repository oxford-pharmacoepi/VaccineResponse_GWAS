# ============================================================================ #
#                            CREATETABLE_VALIDATION                            #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #

CovidSeverity <- as_tibble(read.delim(paste0(dir_results,'/GWAS/breakthrough_ImputedData_covidSeverity_Validation.txt')))
Covid         <- as_tibble(read.delim(paste0(dir_results,'/GWAS/breakthrough_ImputedData_covidSusceptibility_Validation.txt')))
OneDose       <- as_tibble(read.delim(paste0(dir_results,'/GWAS/immuneResponse_ImputedData_one_dose_cohort_Validation.txt')))
TwoDose       <- as_tibble(read.delim(paste0(dir_results,'/GWAS/immuneResponse_ImputedData_two_dose_cohort_Validation.txt')))

OD <- OneDose %>% 
  filter(CHROM == 6) %>% filter(ID != "rs114903158" & ID != "rs3094106") %>%
  mutate(OR = round(exp(BETA), digits = 2)) %>%
  mutate('P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1),
         'Lower'   = exp(BETA-1.96*SE),
         'Upper'   = exp(BETA+1.96*SE)) %>%
  select(CHROM, ID, 'Validation_N' = 'N', 'Validation_OR' = 'OR', 'Validation_P Value' = 'P Value',
         'Validation_Lower' = 'Lower', 'Validation_Upper' = 'Upper') %>%
  mutate(Phenotype = 'One-dose antibody response') %>%
  left_join(
    read.delim(paste0(dir_results,'GWAS/immuneResponse_ImputedData_one_dose_cohort.txt')) %>%
      mutate('Main analysis_N' = N) %>%
      mutate('Main analysis_OR' = round(exp(BETA), digits = 2)) %>%
      mutate('Main analysis_P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
      mutate('Main analysis_Lower' = exp(BETA-1.96*SE)) %>%
      mutate('Main analysis_Upper' = exp(BETA+1.96*SE)) %>%
      select(ID, 'Main analysis_N', 'Main analysis_OR', 'Main analysis_P Value',
             'Main analysis_Lower', 'Main analysis_Upper'),
    by = 'ID'
  )

TD <- TwoDose %>% 
  filter(ID == "rs114903158" | ID == "rs3094106") %>%
  mutate(OR = round(exp(BETA), digits = 2)) %>%
  mutate('P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1),
         'Lower'   = exp(BETA-1.96*SE),
         'Upper'   = exp(BETA+1.96*SE)) %>%
  select(CHROM, ID, 'Validation_N' = 'N', 'Validation_OR' = 'OR', 'Validation_P Value' = 'P Value',
         'Validation_Lower' = 'Lower', 'Validation_Upper' = 'Upper') %>%
  mutate(Phenotype = 'Two-dose antibody response') %>%
  left_join(
    read.delim(paste0(dir_results,'GWAS/immuneResponse_ImputedData_two_dose_cohort.txt')) %>%
      mutate('Main analysis_N' = N) %>%
      mutate('Main analysis_OR' = round(exp(BETA), digits = 2)) %>%
      mutate('Main analysis_P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
      mutate('Main analysis_Lower' = exp(BETA-1.96*SE)) %>%
      mutate('Main analysis_Upper' = exp(BETA+1.96*SE)) %>%
      select(ID, 'Main analysis_N', 'Main analysis_OR', 'Main analysis_P Value',
             'Main analysis_Lower', 'Main analysis_Upper'),
    by = 'ID'
  )

BC <- Covid %>%
  filter(CHROM %in% c(3,19,10)) %>%
  mutate(OR = round(exp(BETA), digits = 2)) %>%
  mutate('P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1),
         'Lower'   = exp(BETA-1.96*SE),
         'Upper'   = exp(BETA+1.96*SE)) %>%
  select(CHROM, ID, 'Validation_N' = 'N', 'Validation_OR' = 'OR', 'Validation_P Value' = 'P Value',
         'Validation_Lower' = 'Lower', 'Validation_Upper' = 'Upper') %>%
  mutate(Phenotype = 'Breakthrough susceptibility') %>%
  left_join(
    read.delim(paste0(dir_results,'GWAS/breakthrough_ImputedData_covidSusceptibility.txt')) %>%
      mutate('Main analysis_N' = N) %>%
      mutate('Main analysis_OR' = round(exp(BETA), digits = 2)) %>%
      mutate('Main analysis_P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
      mutate('Main analysis_Lower' = exp(BETA-1.96*SE)) %>%
      mutate('Main analysis_Upper' = exp(BETA+1.96*SE)) %>%
      select(ID, 'Main analysis_N', 'Main analysis_OR', 'Main analysis_P Value',
             'Main analysis_Lower', 'Main analysis_Upper'),
    by = 'ID'
  )

BS <- CovidSeverity %>%
  filter(CHROM == 16) %>%
  mutate(OR = round(exp(BETA), digits = 2)) %>%
  mutate('P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1),
         'Lower'   = exp(BETA-1.96*SE),
         'Upper'   = exp(BETA+1.96*SE)) %>%
  select(CHROM, ID, 'Validation_N' = 'N', 'Validation_OR' = 'OR', 'Validation_P Value' = 'P Value',
         'Validation_Lower' = 'Lower', 'Validation_Upper' = 'Upper') %>%
  mutate(Phenotype = 'Breakthrough severity') %>%
  left_join(
    read.delim(paste0(dir_results,'GWAS/breakthrough_ImputedData_covidSeverity.txt')) %>%
      mutate('Main analysis_N' = N) %>%
      mutate('Main analysis_OR' = round(exp(BETA), digits = 2)) %>%
      mutate('Main analysis_P Value' = formatC(10^{-LOG10P}, format = "e", digit = 1)) %>%
      mutate('Main analysis_Lower' = exp(BETA-1.96*SE)) %>%
      mutate('Main analysis_Upper' = exp(BETA+1.96*SE)) %>%
      select(ID, 'Main analysis_N', 'Main analysis_OR', 'Main analysis_P Value',
             'Main analysis_Lower', 'Main analysis_Upper'),
    by = 'ID'
  )

gwas <- OD %>% full_join(TD) %>% full_join(BC) %>% full_join(BS)
write.table(gwas,paste0(dir_results,'/Validation/Validation.txt'))

gwas <- gwas %>%
  select(-CHROM) %>%
  relocate('Main analysis_N', .after = 'ID') %>%
  relocate('Main analysis_OR', .after = 'Main analysis_N') %>%
  relocate('Main analysis_P Value', .after = 'Main analysis_OR') %>%
  relocate('Phenotype', .before = 'ID')


gwas1 <- gwas %>%
  flextable() %>%
  span_header(sep = "_",) %>%
  align(align = "center",part = "all")
save_as_docx(gwas1, path=paste0(dir_results,"/Validation/Table_Validation.docx"), align = "center")





