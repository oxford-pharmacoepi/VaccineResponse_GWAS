# ==============================================================================
#                              Validation analysis 
#                          Marta Alcalde-Herraiz, 2023
# ==============================================================================
ch  <- c('oneDose', 
         'twoDose',
         'breakthroughSusceptibility',
         'breakthroughSeverity')

nam <- c('Immune response - One dose',
         'Immune response - Two dose',
         'Breakthrough susceptibility',
         'Breakthrough severity')

leadSNPs <- tibble(
  "chr" = as.integer(),
  "rsID" = as.character()
)

for(i in 1:4){
  leadSNPs <- leadSNPs |>
    union_all(
      read.table(paste0(dir_results,"GWAS/FUMA_",ch[i],"/leadSNPs.txt"), 
                 header = TRUE, sep = "\t") |> as_tibble() |>
        select("chr", "rsID")
    )
} 

for(i in unique(leadSNPs$chr)){
  write_delim(leadSNPs |>
                filter(chr == i) |>
                select(rsID),
              paste0(dir_results,"/Validation/validation_", i,".txt"),
              col_names = FALSE,
              progress = TRUE,
              quote = "none")
}

