# ============================================================================ #
#                              Phenotyping HES                                 #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #

PhenotypingHes <- function(outc,hes,pop){
  phen <- read.xlsx(here("R_Scripts","PhenotypingR.xlsx"),
                    sheetName = outc)
  
  hes_codes <- phen %>% 
    select("HES") %>%
    filter(!is.na(HES))
  
  pop1 <- hes %>% 
    select(eid) %>% 
    filter(eid %in% pop$eid) %>%
    distinct()
  
  h <- hes %>%
    select("eid","diag_icd10") %>%
    filter(diag_icd10 %in% hes_codes$HES) %>% # Filter patients having a code
    select(eid) %>%
    distinct() %>% # No repeated patients
    mutate(state_hes = 1) %>%
    right_join( # add those patients without an outcome
      pop1,
      by = "eid"
    ) %>%
    mutate(state_hes = if_else(is.na(state_hes),0,state_hes)) %>%
    select(eid,state_hes)
  
  return(h)
}

