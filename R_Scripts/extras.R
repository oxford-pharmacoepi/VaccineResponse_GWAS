rm(list = ls())
pacman::p_load('dplyr','tibble','readr','here',
               'lubridate','pbatR','forcats','tableone',
               'qqman','RColorBrewer','ggplot2','gridGraphics','xlsx','stringr',
               'grid','gridExtra','tidyverse','egg','flextable','ftExtra','officer')

# Directory where the data is
dir_data <- 'C:/Users/martaa/Desktop/Projects/VaccineResponse_GWAS/'
dir_ukb  <- 'C:/Users/martaa/Desktop/Projects/UKBiobank_65397/'
dir_results <- paste0(dir_data,'Results/')


ch  <- c('oneDose', 
         'twoDose',
         'breakthroughSusceptibility',
         'breakthroughSeverity')
nam <- c('Immune response - One dose',
         'Immune response - Two dose',
         'Breakthrough susceptibility',
         'Breakthrough severity')

w <- 250e3


genRiskList <- list()
gwasList    <- list()
for(i in 1:4){
  genRiskList[[i]] <- read_delim(paste0(dir_results,'FUMA/',ch[i],'/GenomicRiskLoci.txt')) %>%
    select(rsID, chr, pos, p) %>%
    arrange(p)
  if(i == 1){
    genrisk <- genRiskList[[i]] %>% select(rsID,pos)
  }else{
    genrisk <- genrisk %>% union_all(genRiskList[[i]] %>% select(rsID,pos))
  }
  gwasList[[i]] <- read_delim(paste0(dir_results,'GWAS/imputedData_',ch[i],'.txt'))
}

i <- 1
gwas <- gwasList[[1]]

g1 <- gwas %>%
  filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)








getColocFormat <- function(gwas,genrisk,i,w){
  gwas1 <- gwas %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  d <- list(
    beta = gwas1$BETA,
    varbeta = gwas1$SE^2,
    type = 'cc',
    snp  = gwas1$ID,
    position = gwas1$GENPOS
  )
  
  return(d)
}
