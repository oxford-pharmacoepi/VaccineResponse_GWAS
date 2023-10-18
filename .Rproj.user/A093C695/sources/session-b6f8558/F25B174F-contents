# Colocalisation
if(!require("remotes"))
  install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
library(coloc)

genRiskList <- list()
gwasList    <- list()
for(i in 1:4){
  genRiskList[[i]] <- read_delim(paste0(dir_results,'FUMA/',ch[i],'/GenomicRiskLoci.txt')) %>%
    select(rsID, chr, pos, p) %>%
    arrange(p)
  gwasList[[i]] <- read_delim(paste0(dir_results,'GWAS/imputedData_',ch[i],'.txt'))
}

# Genomic risk loci - Immune response  One dose ================================
w <- 250e3
tableList_coloc_OD <- tibble(
  `Immune response - One dose genetic loci_SNP` = NA,
  `Immune response - One dose genetic loci_CHR` = NA,
  `Immune response - One dose genetic loci_POS` = NA,
  `Immune response - One dose genetic loci_Immune response - Two dose_N` = NA,
  `Immune response - One dose genetic loci_Immune response - Two dose_H4 Prob (%)` = NA,
  `Immune response - One dose genetic loci_Breakthrough susceptibility_N` = NA,
  `Immune response - One dose genetic loci_Breakthrough susceptibility_H4 Prob (%)` = NA,
  `Immune response - One dose genetic loci_Breakthrough severity_N` = NA,
  `Immune response - One dose genetic loci_Breakthrough severity_H4 Prob (%)` = NA
)

for (i in 1:nrow(genRiskList[[1]])){
  gwas1_1 <- gwasList[[1]] %>% 
    filter(GENPOS >= genrisk$pos[i]-w-1& GENPOS <= genrisk$pos[i]+w+1)
  
  gwas2_1 <- gwasList[[2]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas3_1 <- gwasList[[3]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas4_1 <- gwasList[[4]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  d1 <- list(
    beta = gwas1_1$BETA,
    varbeta = gwas1_1$SE^2,
    type = 'cc',
    snp  = gwas1_1$ID,
    position = gwas1_1$GENPOS
  )
  
  d2 <- list(
    beta = gwas2_1$BETA,
    varbeta = gwas2_1$SE^2,
    type = 'cc',
    snp  = gwas2_1$ID,
    position = gwas2_1$GENPOS
  )
  
  d3 <- list(
    beta = gwas3_1$BETA,
    varbeta = gwas3_1$SE^2,
    type = 'cc',
    snp  = gwas3_1$ID,
    position = gwas3_1$GENPOS
  )
  
  d4 <- list(
    beta = gwas4_1$BETA,
    varbeta = gwas4_1$SE^2,
    type = 'cc',
    snp  = gwas4_1$ID,
    position = gwas4_1$GENPOS
  )
  
  c12 <- coloc.abf(dataset1 = d1, dataset2 = d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c13 <- coloc.abf(dataset1 = d1, dataset2 = d3, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c14 <- coloc.abf(dataset1 = d1, dataset2 = d4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  tableList_coloc_OD <- tableList_coloc_OD %>% add_row(tibble(
    `Immune response - One dose genetic loci_SNP` = genRiskList[[1]]$rsID[[i]],
    `Immune response - One dose genetic loci_CHR` = genRiskList[[1]]$chr[[i]],
    `Immune response - One dose genetic loci_POS` = genRiskList[[1]]$pos[[i]],
    `Immune response - One dose genetic loci_Immune response - Two dose_N` = c12$summary[[1]],
    `Immune response - One dose genetic loci_Immune response - Two dose_H4 Prob (%)` = round(c12$summary[[6]]*100, digits = 2),
    `Immune response - One dose genetic loci_Breakthrough susceptibility_N` = c13$summary[[1]],
    `Immune response - One dose genetic loci_Breakthrough susceptibility_H4 Prob (%)` = round(c13$summary[[6]]*100, digits = 2),
    `Immune response - One dose genetic loci_Breakthrough severity_N` = c14$summary[[1]],
    `Immune response - One dose genetic loci_Breakthrough severity_H4 Prob (%)` = round(c14$summary[[6]]*100, digits = 2)
  ))
}

# Genomic risk loci - Immune response  Two dose ================================
w <- 250e3
tableList_coloc_TD <- tibble(
  `Immune response - Two dose genetic loci_SNP` = NA,
  `Immune response - Two dose genetic loci_CHR` = NA,
  `Immune response - Two dose genetic loci_POS` = NA,
  `Immune response - Two dose genetic loci_Immune response - One dose_N` = NA,
  `Immune response - Two dose genetic loci_Immune response - One dose_H4 Prob (%)` = NA,
  `Immune response - Two dose genetic loci_Breakthrough susceptibility_N` = NA,
  `Immune response - Two dose genetic loci_Breakthrough susceptibility_H4 Prob (%)` = NA,
  `Immune response - Two dose genetic loci_Breakthrough severity_N` = NA,
  `Immune response - Two dose genetic loci_Breakthrough severity_H4 Prob (%)` = NA
)

for (i in 1:nrow(genRiskList[[2]])){
  gwas1_1 <- gwasList[[1]] %>% 
    filter(GENPOS >= genrisk$pos[i]-w-1& GENPOS <= genrisk$pos[i]+w+1)
  
  gwas2_1 <- gwasList[[2]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas3_1 <- gwasList[[3]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas4_1 <- gwasList[[4]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  d1 <- list(
    beta = gwas1_1$BETA,
    varbeta = gwas1_1$SE^2,
    type = 'cc',
    snp  = gwas1_1$ID,
    position = gwas1_1$GENPOS
  )
  
  d2 <- list(
    beta = gwas2_1$BETA,
    varbeta = gwas2_1$SE^2,
    type = 'cc',
    snp  = gwas2_1$ID,
    position = gwas2_1$GENPOS
  )
  
  d3 <- list(
    beta = gwas3_1$BETA,
    varbeta = gwas3_1$SE^2,
    type = 'cc',
    snp  = gwas3_1$ID,
    position = gwas3_1$GENPOS
  )
  
  d4 <- list(
    beta = gwas4_1$BETA,
    varbeta = gwas4_1$SE^2,
    type = 'cc',
    snp  = gwas4_1$ID,
    position = gwas4_1$GENPOS
  )
  
  c12 <- coloc.abf(dataset1 = d2, dataset2 = d1, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c13 <- coloc.abf(dataset1 = d2, dataset2 = d3, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c14 <- coloc.abf(dataset1 = d2, dataset2 = d4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  tableList_coloc_TD <- tableList_coloc_TD %>% add_row(tibble(
    `Immune response - Two dose genetic loci_SNP` = genRiskList[[2]]$rsID[[i]],
    `Immune response - Two dose genetic loci_CHR` = genRiskList[[2]]$chr[[i]],
    `Immune response - Two dose genetic loci_POS` = genRiskList[[2]]$pos[[i]],
    `Immune response - Two dose genetic loci_Immune response - One dose_N` = c12$summary[[1]],
    `Immune response - Two dose genetic loci_Immune response - One dose_H4 Prob (%)` = round(c12$summary[[6]]*100, digits = 2),
    `Immune response - Two dose genetic loci_Breakthrough susceptibility_N` = c13$summary[[1]],
    `Immune response - Two dose genetic loci_Breakthrough susceptibility_H4 Prob (%)` = round(c13$summary[[6]]*100, digits = 2),
    `Immune response - Two dose genetic loci_Breakthrough severity_N` = c14$summary[[1]],
    `Immune response - Two dose genetic loci_Breakthrough severity_H4 Prob (%)` = round(c14$summary[[6]]*100, digits = 2)
  ))
}


# Genomic risk loci - BT Susceptibility ========================================
w <- 250e3
tableList_coloc_BTSU <- tibble(
  `Breakthrough susceptibility genetic loci_SNP` = NA,
  `Breakthrough susceptibility genetic loci_CHR` = NA,
  `Breakthrough susceptibility genetic loci_POS` = NA,
  `Breakthrough susceptibility genetic loci_Immune response - One dose_N` = NA,
  `Breakthrough susceptibility genetic loci_Immune response - One dose_H4 Prob (%)` = NA,
  `Breakthrough susceptibility genetic loci_Immune response - Two dose_N` = NA,
  `Breakthrough susceptibility genetic loci_Immune response - Two dose_H4 Prob (%)` = NA,
  `Breakthrough susceptibility genetic loci_Breakthrough severity_N` = NA,
  `Breakthrough susceptibility genetic loci_Breakthrough severity_H4 Prob (%)` = NA
)

for (i in 1:nrow(genRiskList[[3]])){
  gwas1_1 <- gwasList[[1]] %>% 
    filter(GENPOS >= genrisk$pos[i]-w-1& GENPOS <= genrisk$pos[i]+w+1)
  
  gwas2_1 <- gwasList[[2]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas3_1 <- gwasList[[3]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas4_1 <- gwasList[[4]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  d1 <- list(
    beta = gwas1_1$BETA,
    varbeta = gwas1_1$SE^2,
    type = 'cc',
    snp  = gwas1_1$ID,
    position = gwas1_1$GENPOS
  )
  
  d2 <- list(
    beta = gwas2_1$BETA,
    varbeta = gwas2_1$SE^2,
    type = 'cc',
    snp  = gwas2_1$ID,
    position = gwas2_1$GENPOS
  )
  
  d3 <- list(
    beta = gwas3_1$BETA,
    varbeta = gwas3_1$SE^2,
    type = 'cc',
    snp  = gwas3_1$ID,
    position = gwas3_1$GENPOS
  )
  
  d4 <- list(
    beta = gwas4_1$BETA,
    varbeta = gwas4_1$SE^2,
    type = 'cc',
    snp  = gwas4_1$ID,
    position = gwas4_1$GENPOS
  )
  
  c12 <- coloc.abf(dataset1 = d3, dataset2 = d1, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c13 <- coloc.abf(dataset1 = d3, dataset2 = d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c14 <- coloc.abf(dataset1 = d3, dataset2 = d4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  tableList_coloc_BTSU <- tableList_coloc_BTSU %>% add_row(tibble(
    `Breakthrough susceptibility genetic loci_SNP` = genRiskList[[3]]$rsID[[i]],
    `Breakthrough susceptibility genetic loci_CHR` = genRiskList[[3]]$chr[[i]],
    `Breakthrough susceptibility genetic loci_POS` = genRiskList[[3]]$pos[[i]],
    `Breakthrough susceptibility genetic loci_Immune response - One dose_N` = c12$summary[[1]],
    `Breakthrough susceptibility genetic loci_Immune response - One dose_H4 Prob (%)` = round(c12$summary[[6]]*100, digits = 2),
    `Breakthrough susceptibility genetic loci_Immune response - Two dose_N` = c13$summary[[1]],
    `Breakthrough susceptibility genetic loci_Immune response - Two dose_H4 Prob (%)` = round(c13$summary[[6]]*100, digits = 2),
    `Breakthrough susceptibility genetic loci_Breakthrough severity_N` = c14$summary[[1]],
    `Breakthrough susceptibility genetic loci_Breakthrough severity_H4 Prob (%)` = round(c14$summary[[6]]*100, digits = 2)
  ))
}

# Genomic risk loci - BT Severity ==============================================
w <- 250e3
tableList_coloc_BTSE <- tibble(
  `Breakthrough severity genetic loci_SNP` = NA,
  `Breakthrough severity genetic loci_CHR` = NA,
  `Breakthrough severity genetic loci_POS` = NA,
  `Breakthrough severity genetic loci_Immune response - One dose_N` = NA,
  `Breakthrough severity genetic loci_Immune response - One dose_H4 Prob (%)` = NA,
  `Breakthrough severity genetic loci_Immune response - Two dose_N` = NA,
  `Breakthrough severity genetic loci_Immune response - Two dose_H4 Prob (%)` = NA,
  `Breakthrough severity genetic loci_Breakthrough Susceptibility_N` = NA,
  `Breakthrough severity genetic loci_Breakthrough Susceptibility_H4 Prob (%)` = NA
)

for (i in 1:nrow(genRiskList[[4]])){
  gwas1_1 <- gwasList[[1]] %>% 
    filter(GENPOS >= genrisk$pos[i]-w-1& GENPOS <= genrisk$pos[i]+w+1)
  
  gwas2_1 <- gwasList[[2]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas3_1 <- gwasList[[3]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  gwas4_1 <- gwasList[[4]] %>%
    filter(GENPOS >= genrisk$pos[i]-w-1 & GENPOS <= genrisk$pos[i]+w+1)
  
  d1 <- list(
    beta = gwas1_1$BETA,
    varbeta = gwas1_1$SE^2,
    type = 'cc',
    snp  = gwas1_1$ID,
    position = gwas1_1$GENPOS
  )
  
  d2 <- list(
    beta = gwas2_1$BETA,
    varbeta = gwas2_1$SE^2,
    type = 'cc',
    snp  = gwas2_1$ID,
    position = gwas2_1$GENPOS
  )
  
  d3 <- list(
    beta = gwas3_1$BETA,
    varbeta = gwas3_1$SE^2,
    type = 'cc',
    snp  = gwas3_1$ID,
    position = gwas3_1$GENPOS
  )
  
  d4 <- list(
    beta = gwas4_1$BETA,
    varbeta = gwas4_1$SE^2,
    type = 'cc',
    snp  = gwas4_1$ID,
    position = gwas4_1$GENPOS
  )
  
  c12 <- coloc.abf(dataset1 = d4, dataset2 = d1, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c13 <- coloc.abf(dataset1 = d4, dataset2 = d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  c14 <- coloc.abf(dataset1 = d4, dataset2 = d3, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  tableList_coloc_BTSE <- tableList_coloc_BTSE %>% add_row(tibble(
    `Breakthrough severity genetic loci_SNP` = genRiskList[[4]]$rsID[[i]],
    `Breakthrough severity genetic loci_CHR` = genRiskList[[4]]$chr[[i]],
    `Breakthrough severity genetic loci_POS` = genRiskList[[4]]$pos[[i]],
    `Breakthrough severity genetic loci_Immune response - One dose_N` = c12$summary[[1]],
    `Breakthrough severity genetic loci_Immune response - One dose_H4 Prob (%)` = round(c12$summary[[6]]*100, digits = 2),
    `Breakthrough severity genetic loci_Immune response - Two dose_N` = c13$summary[[1]],
    `Breakthrough severity genetic loci_Immune response - Two dose_H4 Prob (%)` = round(c13$summary[[6]]*100, digits = 2),
    `Breakthrough severity genetic loci_Breakthrough Susceptibility_N` = c14$summary[[1]],
    `Breakthrough severity genetic loci_Breakthrough Susceptibility_H4 Prob (%)` = round(c14$summary[[6]]*100, digits = 2)
  ))

}

# Merge all the tables =========================================================
tableList_coloc <- list()
tableList_coloc[[1]] <- tableList_coloc_OD %>% 
  filter(!is.na(`Immune response - One dose genetic loci_SNP`)) %>%
  flextable() %>%
  span_header(sep = "_") %>%
  align(align = 'center', part = 'all') %>%
  bg(i = 1, bg = "#C8C8C8", part = "head") %>%
  bg(i = 2, j = c(4,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 3, j = c(4,6,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 1:2, j = c(4,6,8), bg = "#EFEFEF", part = "body") %>%
  bold(bold = TRUE, part="header") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all") %>%
  vline(j = seq(3,7,2))
  
tableList_coloc[[2]] <- tableList_coloc_TD %>% 
  filter(!is.na(`Immune response - Two dose genetic loci_SNP`)) %>%
  flextable() %>%
  span_header(sep = "_") %>%
  align(align = 'center', part = 'all') %>%
  bg(i = 1, bg = "#C8C8C8", part = "head") %>%
  bg(i = 2, j = c(4,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 3, j = c(4,6,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 1:2, j = c(4,6,8), bg = "#EFEFEF", part = "body") %>%
  bold(bold = TRUE, part="header") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all") %>%
  vline(j = seq(3,7,2))

tableList_coloc[[3]] <- tableList_coloc_BTSU %>% 
  filter(!is.na(`Breakthrough susceptibility genetic loci_SNP`)) %>%
  flextable() %>%
  span_header(sep = "_") %>%
  align(align = 'center', part = 'all') %>%
  bg(i = 1, bg = "#C8C8C8", part = "head") %>%
  bg(i = 2, j = c(4,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 3, j = c(4,6,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 1:2, j = c(4,6,8), bg = "#EFEFEF", part = "body") %>%
  bold(bold = TRUE, part="header") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all") %>%
  vline(j = seq(3,7,2))

tableList_coloc[[4]] <- tableList_coloc_BTSE %>% 
  filter(!is.na(`Breakthrough severity genetic loci_SNP`)) %>%
  flextable() %>%
  span_header(sep = "_") %>%
  align(align = 'center', part = 'all') %>%
  bg(i = 1, bg = "#C8C8C8", part = "head") %>%
  bg(i = 2, j = c(4,8), bg = "#EFEFEF", part = "head") %>%
  bg(i = 3, j = c(4,6,8), bg = "#EFEFEF", part = "head") %>%
  bg(j = c(4,6,8), bg = "#EFEFEF", part = "body") %>%
  bold(bold = TRUE, part="header") %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname='Calibri',part = "all") %>%
  vline(j = seq(3,7,2))
save_as_docx(values = tableList_coloc, path = paste0(dir_results,'Tables/Colocalisation.docx'))


