# ============================================================================ #
#                                 CODE TO RUN                                  #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
rm(list = ls())
pacman::p_load('dplyr','tibble','readr','here',
               'lubridate','pbatR','forcats','tableone',
               'qqman','RColorBrewer','ggplot2','gridGraphics','xlsx',
               'grid','gridExtra','tidyverse','egg','flextable','ftExtra','officer')

# Directory where the data is
dir_data <- 'C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/UKB_69741/GWAS_VaccineImmune/'
dir_results <- paste0(dir_data,'Results/')


# SCRIPTS TO MAKE THE COHORTS --------------------------------------------------
dir.create(paste0(dir_results,'Cohorts'))

# Immune response cohort:
source(here('R_Scripts','immuneResponse_ImputedData.R'))

# Breakthrough infection cohort:
source(here('R_Scripts','breakthrough_ImputedData.R'))


# SCRIPTS TO MAKE THE TABLES FROM THE PAPER ------------------------------------
snps <- c('rs9461694','rs2523496','6:32255712_CAGTT_C','rs9268465','6:32419074_CT_C',
          'rs9268847','6:32440321_CTG_C','rs565122319','rs145945003','rs7763805',
          'rs3094106','rs114903158',
          'rs59776512','rs73062389','rs71322420','rs16861415','rs13097481','rs1977830','rs778809','rs11673136','rs681343',
          'rs62038344')

dir.create(paste0(dir_results,'Tables'))
source(here('R_Scripts','CreateTable.R'))
source(here('R_Scripts','CreateTable_Validation.R'))


# hes  <- as_tibble(read.delim(paste0(dir_data,'hesin_diag.txt'), sep = "\t", quote = ""))
# hesD <- as_tibble(read.delim(paste0(dir_data,'hesin.txt'), sep = "\t", quote = ""))
# gp   <- as_tibble(read.delim(paste0(dir_data,'gp_clinical.txt'), sep = "\t", quote = ""))
# ukb <- as_tibble(read.delim(paste0(dir_data,'ukb669864_logistic.csv'), sep = "\t", quote = ""))
hes  <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/MR_Sclerostin/MR_Sclerostin_Data/UKB_old/hesin_diag.txt',sep = "\t", quote = ""))
hesD <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/MR_Sclerostin/MR_Sclerostin_Data/UKB_old/hesin.txt', sep = "\t", quote = ""))
gp   <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/MR_Sclerostin/MR_Sclerostin_Data/UKB_old/gp_clinical.txt',sep = "\t", quote = ""))
gp <- gp %>% filter(event_dt != "01/01/1900",
                    event_dt != "01/01/1901",
                    event_dt != "02/02/1902",
                    event_dt != "03/03/1903",
                    event_dt != "07/07/2037")
ukb  <- as_tibble(read.csv('C:/Users/martaa/Desktop/Projects/MR_Sclerostin/MR_Sclerostin_Data/UKB_old/ukb669864_logistic.csv'))


# SCRIPTS TO MAKE THE FIGURES FROM THE PAPER -----------------------------------
dir.create(paste0(dir_results,'Figures'))
source(here('R_Scripts','Figure.R'))
source(here('R_Scripts','CreateFigure_Validation.R'))




