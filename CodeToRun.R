# ============================================================================ #
#                                 CODE TO RUN                                  #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
rm(list = ls())
library("dplyr")
library("here")
library("lubridate")
library("pbatR")
library("tableone")
library("ggplot2")
library("stringr")
library("tidyverse")
library("flextable")
library("ftExtra")
library("coloc")
library("tidyverse")

source(here("R_Scripts/Functions.R"))

# Directory where the data is
dir_data    <- 'D:/Projects/VaccineResponse_GWAS/'
dir_ukb     <- 'D:/Projects/UKBiobank_65397/'
dir_results <- 'D:/Projects/VaccineResponse_GWAS/Results/'

# Load codes
source(here('R_Scripts','1-ImmuneResponse.R'))
source(here('R_Scripts','2-Breakthrough.R'))

# Run the GWAS in the RAP PLATFORM ---------------------------------------------
source(here("R_Scripts","3-ComputePVal.R"))

# Save the results within the "GWAS" folder under the names:
# imputedData_breakthroughSeverity.txt
# imputedData_breakthroughSeverity_Validation.txt


# Run FUMA ---------------------------------------------------------------------

# Colocalisation analysis ------------------------------------------------------
source(here("R_Scripts","4-Colocalization.R"))

# Validation -------------------------------------------------------------------

# SCRIPTS TO MAKE THE TABLES/FIGURES FROM THE PAPER ----------------------------
dir.create(paste0(dir_results,'Validation'))
source(here("R_Scripts","createTables.R"))
source(here("R_Scripts","CreateFigure_Validation.R"))
dir.create(paste0(dir_results,'Tables'))
source(here('R_Scripts','CreateTable.R'))
source(here('R_Scripts','CreateManhattanPlot.R'))


# SCRIPTS TO MAKE THE FIGURES FROM THE PAPER -----------------------------------
source(here('R_Scripts','Figure.R'))


t <- read.table("C:/Users/martaa/OneDrive - Nexus365/Marta/Projects/Breakthough_infection_GWAS/manuscript_versions/1-first_revision/Results/GWAS/oneDose_assoc.regenie.merged_new.txt", header = TRUE)
t1 <- read.table("D:/Projects/VaccineResponse_GWAS/Results_old/GWAS/imputedData_oneDose.txt", header = TRUE)


cohort_old <- read.table("D:/Projects/VaccineResponse_GWAS/Results_old/Cohorts/imputedData_oneDose.phe", sep = " ", header = TRUE) |> as_tibble()
cohort_new <- read.table("D:/Projects/VaccineResponse_GWAS/Results/Cohorts/one_dose_original.phe", sep = " ", header = TRUE) |> as_tibble()
