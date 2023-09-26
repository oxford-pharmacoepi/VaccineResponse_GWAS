# ============================================================================ #
#                                 Table One                                    #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #
t_onedose <- as_tibble(read.csv(paste0(dir_results,'Cohorts/immuneResponse_ImputedData_one_dose_cohort.csv'), row.names = 'X')) %>% rename('eid' = 'FID') %>% select(-IID, -PC1, -PC2, -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -Genetic_batch)
t_twodose <- as_tibble(read.csv(paste0(dir_results,'Cohorts/immuneResponse_ImputedData_two_dose_cohort.csv'), row.names = 'X')) %>% rename('eid' = 'FID') %>% select(-IID, -PC1, -PC2, -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -Genetic_batch)
t_bt      <- as_tibble(read.csv(paste0(dir_results,'Cohorts/breakthrough_ImputedData_covidSusceptibility.csv'), row.names = 'X')) %>% rename('eid' = 'FID') %>% select(-IID, -PC1, -PC2, -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -Genetic_batch)
t_bts     <- as_tibble(read.csv(paste0(dir_results,'Cohorts/breakthrough_ImputedData_covidSeverity.csv'), row.names = 'X')) %>% rename('eid' = 'FID') %>% select(-IID, -PC1, -PC2, -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -Genetic_batch)

tableList <- list(t_onedose, t_twodose, t_bt, t_bts)
outcomes  <- c('immuneResponse', 'immuneResponse', 'bt_infection', 'severity')
outcomes1 <- c('Immune response - One dose', 'Immune response - Two dose','Breakthrough susceptibility','Breakthrough severity')
tableOnes <- list()

for (i in 1:4){
  pop <- tableList[[i]] %>% select(eid)
  outc <- c('CAD','MI','IS','Hypertension','T2DM')
  for (j in outc){
    source(here('R_Scripts','PhenotypingHES.R'))
    hes_data <- PhenotypingHes(j,hes,pop)
    
    source(here("R_Scripts","PhenotypingGP.R"))
    gp_data <- PhenotypingGP(j,gp,pop)
    
    source(here("R_Scripts","PhenotypingUKB.R"))
    ukb_data <- as_tibble(PhenotypingUKB(j,pop,ukb))
    
    if (j == 'CAD'){
      ukb_data[1:nrow(ukb_data),'state_ukb'] <- NA
    }
    
    t <- pop %>% 
      left_join(hes_data, by = 'eid') %>% 
      left_join(gp_data, by = 'eid') %>%
      left_join(ukb_data %>% select('eid','state_ukb'), by = 'eid')
    
    t[is.na(t)] <- 2
    
    t <- t %>% mutate(state = if_else(state_hes == 1 | state_gp == 1 | state_ukb == 1,1,2))
    t[t$state == 2 & (t$state_hes == 0 | t$state_gp == 0), 'state'] <- 0
    t[t$state == 2, 'state'] <- NA
    
    tableList[[i]] <- tableList[[i]] %>% left_join(t %>% select(eid, state))
    names(tableList[[i]])[names(tableList[[i]]) == 'state'] <- j
  }
  
  
  tableList[[i]] <- tableList[[i]] %>% 
    mutate_at(vars('CAD','MI','IS','Hypertension','T2DM'), function(x) {as.factor(x)}) %>%
    rename('Coronary artery disease' = 'CAD',
           'Myocardial infarction' = 'MI',
           'Ischaemic stroke' = 'IS',
           'Type 2 diabetes mellitus' = 'T2DM',
           'Body mass index' = 'BMI')


  tableOnes[[i]] <- CreateTableOne(vars = colnames(tableList[[i]])[3:length(colnames(tableList[[i]]))], 
                       data = tableList[[i]],
                       strata = outcomes[i],
                       test = FALSE,
                       smd = TRUE,
                       includeNA = TRUE) %>%
    print(smd = TRUE, showAllLevels = TRUE)
  
  tableOnes[[i]] <- as_tibble(tableOnes[[i]]) %>%
    mutate(Variables = rownames(tableOnes[[i]])) %>%
    relocate(Variables) %>%
    rename("Controls" = "0",
           "Cases" = "1") %>% 
    mutate("level" = case_when(level == "0" ~ 'Control',
                               level == "1" ~ 'Cases',
                               is.na(level) ~ 'Missings',
                               .default = " ")) %>%
    rename(" " = "level") %>%
    flextable() %>%
    align(align = "center",part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    font(fontname='Calibri',part = "all") %>%
    set_caption(caption = outcomes1[i]) %>%
    bg(i = 1, bg = "#EFEFEF", part = "header") %>%
    hline(i = c(1,2,3,4,7,10,13,16), part = 'body') %>%
    width(j = c(3,4), 2.7, unit = 'cm') %>%
    width(j = 1, 3, unit = 'cm')
}


save_as_docx(values = tableOnes, path = paste0(dir_results,'Tables/TableOne.docx'), align = "center")




