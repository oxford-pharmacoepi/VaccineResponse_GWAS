rm(list=ls())
pacman::p_load('here','readr','qqman','dplyr','RColorBrewer','ggplot2','gridGraphics',
               'grid','ggplot2','gridExtra','tidyverse','egg','flextable','ftExtra','officer')

gwas_onedose <- read_delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/immuneResponse_ImputedData_one_dose_cohort.txt')) %>%
  filter(ID %in% c('rs9461694','rs2523496','6:32255712_CAGTT_C','rs9268465','6:32419074_CT_C',
                   'rs9268847','6:32440321_CTG_C','rs565122319','rs145945003','rs7763805',
                   'rs3094106','rs114903158',
                   'rs59776512','rs73062389','rs71322420','rs16861415','rs13097481','rs1977830','rs778809','rs11673136','rs681343',
                   'rs62038344')) %>%
  mutate(P = 10^(-LOG10P))  %>%
  select(CHR = CHROM, BP = GENPOS, SNP = ID, ALLELE1, 
         'Immune response_One dose_EAF' = A1FREQ, 
         'Immune response_One dose_Beta' = BETA,
         'Immune response_One dose_P Value' = P)

gwas_twodose <-  read_delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/immuneResponse_ImputedData_two_dose_cohort.txt')) %>%
  filter(ID %in% c('rs9461694','rs2523496','6:32255712_CAGTT_C','rs9268465','6:32419074_CT_C',
                   'rs9268847','6:32440321_CTG_C','rs565122319','rs145945003','rs7763805',
                   'rs3094106','rs114903158',
                   'rs59776512','rs73062389','rs71322420','rs16861415','rs13097481','rs1977830','rs778809','rs11673136','rs681343',
                   'rs62038344')) %>%
  mutate(P = 10^(-LOG10P))  %>%
  select(CHR = CHROM, BP = GENPOS, SNP = ID, ALLELE1, 
         'Immune response_Two dose_EAF' = A1FREQ, 
         'Immune response_Two dose_Beta' = BETA,
         'Immune response_Two dose_P Value' = P)

gwas_suscep <- read_delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/breakthrough_ImputedData_covidSusceptibility.txt')) %>%
  filter(ID %in% c('rs9461694','rs2523496','6:32255712_CAGTT_C','rs9268465','6:32419074_CT_C',
                   'rs9268847','6:32440321_CTG_C','rs565122319','rs145945003','rs7763805',
                   'rs3094106','rs114903158',
                   'rs59776512','rs73062389','rs71322420','rs16861415','rs13097481','rs1977830','rs778809','rs11673136','rs681343',
                   'rs62038344')) %>%
  mutate(P = 10^(-LOG10P))  %>%
  select(CHR = CHROM, BP = GENPOS, SNP = ID, ALLELE1, 
         'Breakthrough_Susceptibility_EAF' = A1FREQ, 
         'Breakthrough_Susceptibility_Beta' = BETA,
         'Breakthrough_Susceptibility_P Value' = P)

gwas_sever <-  read_delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/breakthrough_ImputedData_covidSeverity.txt')) %>%
  filter(ID %in% c('rs9461694','rs2523496','6:32255712_CAGTT_C','rs9268465','6:32419074_CT_C',
                   'rs9268847','6:32440321_CTG_C','rs565122319','rs145945003','rs7763805',
                   'rs3094106','rs114903158',
                   'rs59776512','rs73062389','rs71322420','rs16861415','rs13097481','rs1977830','rs778809','rs11673136','rs681343',
                   'rs62038344')) %>%
  mutate(P = 10^(-LOG10P))  %>%
  select(CHR = CHROM, BP = GENPOS, SNP = ID, ALLELE1, 
         'Breakthrough_Severity_EAF' = A1FREQ, 
         'Breakthrough_Severity_Beta' = BETA,
         'Breakthrough_Severity_P Value' = P)

mapping <- read_table(here('Mapping','snps_1.txt')) %>% select(rsID,nearestGene) %>%
  full_join(read_table(here('Mapping','snps_2.txt'))) %>% select(rsID,nearestGene) %>%
  full_join(read_table(here('Mapping','snps_3.txt')))  %>% select(rsID,nearestGene)%>%
  full_join(read_table(here('Mapping','snps_4.txt')))  %>% select(rsID,nearestGene)

gwas <- gwas_onedose %>%
  inner_join(gwas_twodose) %>% inner_join(gwas_suscep) %>% inner_join(gwas_sever) %>%
  left_join(mapping %>% rename(SNP = rsID)) %>%
  mutate('Immune response_One dose_Beta' = exp(`Immune response_One dose_Beta`),
         'Immune response_Two dose_Beta' = exp(`Immune response_Two dose_Beta`),
         'Breakthrough_Susceptibility_Beta' = exp(`Breakthrough_Susceptibility_Beta`),
         'Breakthrough_Severity_Beta' = exp(`Breakthrough_Severity_Beta`)) %>%
  arrange(CHR) %>%
  mutate(`Immune response_One dose_EAF` = round(`Immune response_One dose_EAF`, digits = 2),
         `Immune response_One dose_Beta` = round(`Immune response_One dose_Beta`, digits = 2),
         `Immune response_One dose_P Value` = formatC(`Immune response_One dose_P Value`, format = "e", digits = 1),
         `Immune response_Two dose_EAF` = round(`Immune response_Two dose_EAF`, digits = 2),
         `Immune response_Two dose_Beta` = round(`Immune response_Two dose_Beta`, digits = 2),
         `Immune response_Two dose_P Value` = formatC(`Immune response_Two dose_P Value`, format = "e", digits = 1),
         `Breakthrough_Susceptibility_EAF` = round(`Breakthrough_Susceptibility_EAF`, digits = 2),
         `Breakthrough_Susceptibility_Beta` = round(`Breakthrough_Susceptibility_Beta`, digits = 2),
         `Breakthrough_Susceptibility_P Value` = formatC(`Breakthrough_Susceptibility_P Value`, format = "e", digits = 1),
         `Breakthrough_Severity_EAF` = round(`Breakthrough_Severity_EAF`, digits = 2),
         `Breakthrough_Severity_Beta` = round(`Breakthrough_Severity_Beta`, digits = 2),
         `Breakthrough_Severity_P Value` = formatC(`Breakthrough_Severity_P Value`, format = "e", digits = 1)) %>%
  select(-BP) %>%
  rename("EA" = ALLELE1, 'Nearest gene' = nearestGene) %>%
  mutate("SNP" = if_else(SNP == "6:32419074_CT_C","rs2150392827",SNP)) %>%
  relocate(CHR) %>%
  relocate('Nearest gene', .after = SNP) %>%
  flextable() %>%
  span_header(sep = "_",) %>%
  align(i = 1, align = 'center', part = "header") %>%
  align(i = 2, align = 'center', part = "header") %>%
  align(i = 3, align = 'center', part = "header") %>%
  align(align = "center",part = "all") %>%
  vline(i = 1, j = c(4:8), part = "header") %>%
  vline(i = 2, j = seq(4,12,1), part = "header") %>%
  vline(i = 3, j = seq(4,15,3), part = "header") %>%
  vline(j = seq(4,15,3), part = "body") %>%
  bg(j = "Immune response_One dose_EAF", bg = "#EFEFEF", part = "all") %>%
  bg(j = "Immune response_One dose_P Value", bg = "#EFEFEF", part = "all") %>%
  bg(j = "Immune response_Two dose_Beta", bg = "#EFEFEF", part = "body") %>%
  bg(i = 3, j = c(9,11,13,15), bg = "#EFEFEF", part = "header") %>%
  bg(i = 2, j = c(11,13), bg = "#EFEFEF", part = "header") %>%
  bg(j = "Breakthrough_Susceptibility_EAF", bg = "#EFEFEF", part = "body") %>%
  bg(j = "Breakthrough_Susceptibility_P Value",bg = "#EFEFEF", part = "body") %>%
  bg(j = "Breakthrough_Severity_Beta",bg = "#EFEFEF", part = "body") %>%
  width(j = 1, width = 1, unit = "in") %>%
  bold(bold = TRUE, part="header") %>%
  width(j = "Immune response_One dose_EAF", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_Two dose_EAF", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_EAF", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_EAF", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_One dose_Beta", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_Two dose_Beta", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_Beta", width = 1.2, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_Beta", width = 1.2, unit = 'cm') %>%
  width(j = "Immune response_One dose_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Immune response_Two dose_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Breakthrough_Susceptibility_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "Breakthrough_Severity_P Value", width = 1.6, unit = 'cm') %>%
  width(j = "SNP", width = 2.7, unit = 'cm') %>%
  width(j = "CHR", width = 2.65, unit = 'cm') %>%
  width(j = "Nearest gene", width = 2.3, unit = 'cm') %>%
  width(j = "EA", width = 1.1, unit = 'cm') %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all")
save_as_docx(gwas, path=here("Table.docx"), align = "center")
