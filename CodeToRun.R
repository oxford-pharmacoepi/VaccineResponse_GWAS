# ============================================================================ #
#                                 CODE TO RUN                                  #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
rm(list = ls())
library("dplyr")
library("readr")
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

# Create directories
dir.create(paste0(dir_results,'Cohorts'))

# Load codes
source(here('R_Scripts','1-ImmuneResponse.R'))
source(here('R_Scripts','2-Breakthrough.R'))

# Run the GWAS in the RAP PLATFORM ---------------------------------------------
dir.create(paste0(dir_results,'GWAS'))
# Save the results within the "GWAS" folder under the names:
# imputedData_breakthroughSeverity.txt
# imputedData_breakthroughSeverity_Validation.txt

# Run FUMA ---------------------------------------------------------------------
# Run FUMA with the following specifications:
# Note: You may need to convert to .gz files the GWAS
# [INPUT FILES]
# - chrcol = CHROM
# - poscol = GENPOS
# - rsIDcol = ID
# - pcol = P
# - eacol = ALLELE1
# - neacol = ALLELE0
# - orcol = NA
# - becol = BETA
# - secol = SE
# - leadSNPsfile = NA
# - addleadSNPs = 1
# - regionsfile = NA
# [PARAMETERS]
# - GRCh38 = 0
# - N = NA
# - Ncol = N
# - exMHC = 0
# - MHCopt = NA
# - extMHC = NA
# - ensembl = v102
# - genetype = protein_coding
# - leadP = 5e-8
# - gwasP = 0.05
# - r2 = 0.6
# - r2_2 = 0.1
# - refpanel = 1KG/Phase3
# - pop = EUR
# - MAF = 0
# - refSNPs = 1
# - mergeDist = 250
# [MAGMA]
# - magma = 0
# [posMap]
# - posMap = 1
# - posMapWindowSize = 10
# - posMapAnnot = NA
# - posMapCADDth = 0
# - posMapRDBth = NA
# - posMapChr15 = NA
# - posMapChr15Max = NA
# - posMapChr15Meth = NA
# - posMapAnnoDs = NA
# - posMapAnnoMeth = NA
# [eqtlMap]
# - eqtlMap = 0
# [ciMap]
# - ciMap = 0

dir.create(paste0(dir_results,'FUMA'))
dir.create(paste0(dir_results,'FUMA/oneDose'))
dir.create(paste0(dir_results,'FUMA/twoDose'))
dir.create(paste0(dir_results,'FUMA/breakthroughSusceptibility'))
dir.create(paste0(dir_results,'FUMA/breakthroughSeverity'))
# Save the results within each folder. Files you may have to download include:
# - GenomicRisk.Loci.txt
# - leadSNPs.txt
# - IndSigSNPs.txt
# - annov.txt
# - annov.stats.txt
# - gwascatalog.txt
# - params.config
# - README.txt for more detailed information

# Colocalisation analysis: 
dir.create(paste0(dir_results,'Tables'))
dir.create(paste0(dir_results,'Figures'))
source(here("R_Scripts","3-Colocalization.R"))
source(here("R_Scripts","Colocalisation.R"))

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
