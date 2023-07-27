rm(list=ls())
pacman::p_load('dplyr','tibble', 'flextable', 'ftExtra', 'here', 'ggplot2', 'forcats')

CovidSeverity <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/Validation/CovidSeverityValidation.txt'))
Covid         <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/Validation/CovidValidation.txt'))
OneDose       <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/Validation/OneDoseValidation.txt'))
TwoDose       <- as_tibble(read.delim('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/Validation/TwoDoseValidation.txt'))

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
    read.delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/immuneResponse_ImputedData_one_dose_cohort.txt')) %>%
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
    read.delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/immuneResponse_ImputedData_two_dose_cohort.txt')) %>%
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
    read.delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/breakthrough_ImputedData_covidSusceptibility.txt')) %>%
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
    read.delim(paste0('C:/Users/martaa/Desktop/Projects/GWAS_VaccineResponse/breakthrough_ImputedData_covidSeverity.txt')) %>%
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
write.table(gwas,here('validation.txt'))

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
save_as_docx(gwas1, path=here("Table_Validation.docx"), align = "center")

######################################33
gwas <- read.table(here('validation.txt'))

gwas <- gwas %>%
  mutate(ID = if_else(ID == '6:32419074_CT_C', 'rs2150392827', ID)) %>%
  mutate(ID = if_else(ID == '6:32440321_CTG_C', 'rs17202941', ID)) %>%
  mutate(Or = c(9,10,4,5,1,2,6,7,3,8,12,11,21,19,18,20,14,13,16,15,17,22)) %>%
  arrange(Or)

dodge <- 0.5
w <- 0.3
ww <- 0.6
# One dose plot
od <- gwas %>% filter(Phenotype == 'One-dose antibody response') %>%
  select(c('CHROM','ID','Phenotype','Or') | contains('Validation')) %>%
  mutate(Type = 'Validation') %>%
  rename('N' = 'Validation_N',
         'OR' = 'Validation_OR',
         'P.Value' = 'Validation_P.Value',
         'Lower' = 'Validation_Lower',
         'Upper' = 'Validation_Upper') %>%
  full_join(gwas %>% filter(Phenotype == 'One-dose antibody response') %>%
              select(c('CHROM','ID','Phenotype','Or') | contains('Main.analysis')) %>%
              mutate(Type = 'Main') %>%
              rename('N' = 'Main.analysis_N',
                     'OR' = 'Main.analysis_OR',
                     'P.Value' = 'Main.analysis_P.Value',
                     'Lower' = 'Main.analysis_Lower',
                     'Upper' = 'Main.analysis_Upper')) %>%
  arrange(desc(Or)) %>%
  mutate(Type = factor(Type, levels = c('Validation','Main')))

od_plot <- ggplot(od,aes(OR, fct_inorder(ID))) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper, color = Type),
    position = position_dodge(width = dodge), width = w, cex = ww
  ) +
  geom_point(aes(color = Type), size = 2, 
             position = position_dodge(width = dodge)) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("skyblue2","black")) +
  theme_gray() +
  facet_wrap(~Phenotype, strip.position = 'top') +
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_x_log10(breaks = c(0.5,0.75,1,1.25,1.5), limits = c(0.45, 1.55))
 

# Two dose plot
td <- gwas %>% filter(Phenotype == 'Two-dose antibody response') %>%
  select(c('CHROM','ID','Phenotype','Or') | contains('Validation')) %>%
  mutate(Type = 'Validation') %>%
  rename('N' = 'Validation_N',
         'OR' = 'Validation_OR',
         'P.Value' = 'Validation_P.Value',
         'Lower' = 'Validation_Lower',
         'Upper' = 'Validation_Upper') %>%
  full_join(gwas %>% filter(Phenotype == 'Two-dose antibody response') %>%
              select(c('CHROM','ID','Phenotype','Or') | contains('Main.analysis')) %>%
              mutate(Type = 'Main') %>%
              rename('N' = 'Main.analysis_N',
                     'OR' = 'Main.analysis_OR',
                     'P.Value' = 'Main.analysis_P.Value',
                     'Lower' = 'Main.analysis_Lower',
                     'Upper' = 'Main.analysis_Upper')) %>%
  arrange(Or) %>%
  mutate(Type = factor(Type, levels = c('Validation','Main')))
td_plot <- ggplot(td,
                  aes(OR, fct_rev(fct_inorder(ID)))) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper, color = Type),
    position = position_dodge(width = dodge), width = w, cex = ww
  ) +
  geom_point(aes(color = Type), size = 2, 
             position = position_dodge(width = dodge)) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("skyblue2","black")) +
  theme_gray() +
  facet_wrap(~Phenotype, strip.position = 'top') +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_x_log10(breaks = c(0.5,0.75,1,1.25,1.5), limits = c(0.45, 1.55))

# Breakthrough
bc <- gwas %>% filter(Phenotype == 'Breakthrough susceptibility') %>%
  select(c('CHROM','ID','Phenotype','Or') | contains('Validation')) %>%
  mutate(Type = 'Validation') %>%
  rename('N' = 'Validation_N',
         'OR' = 'Validation_OR',
         'P.Value' = 'Validation_P.Value',
         'Lower' = 'Validation_Lower',
         'Upper' = 'Validation_Upper') %>%
  full_join(gwas %>% filter(Phenotype == 'Breakthrough susceptibility') %>%
              select(c('CHROM','ID','Phenotype','Or') | contains('Main.analysis')) %>%
              mutate(Type = 'Main') %>%
              rename('N' = 'Main.analysis_N',
                     'OR' = 'Main.analysis_OR',
                     'P.Value' = 'Main.analysis_P.Value',
                     'Lower' = 'Main.analysis_Lower',
                     'Upper' = 'Main.analysis_Upper')) %>%
  arrange(Or) %>%
  mutate(Type = factor(Type, levels = c('Validation','Main')))
bc_plot <- ggplot(bc,
                  aes(OR, fct_rev(fct_inorder(ID)))) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper, color = Type),
    position = position_dodge(width = dodge), width = w, cex = ww
  ) +
  geom_point(aes(color = Type), size = 2, 
             position = position_dodge(width = dodge)) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("skyblue2","black")) +
  theme_gray() +
  facet_wrap(~Phenotype, strip.position = 'top') +
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") + 
  scale_x_log10(breaks = c(0.5,0.75,1,1.25,1.5), limits = c(0.45, 1.55))

# Breakthrough
bs <- gwas %>% filter(Phenotype == 'Breakthrough severity') %>%
  select(c('CHROM','ID','Phenotype','Or') | contains('Validation')) %>%
  mutate(Type = 'Validation') %>%
  rename('N' = 'Validation_N',
         'OR' = 'Validation_OR',
         'P.Value' = 'Validation_P.Value',
         'Lower' = 'Validation_Lower',
         'Upper' = 'Validation_Upper') %>%
  full_join(gwas %>% filter(Phenotype == 'Breakthrough severity') %>%
              select(c('CHROM','ID','Phenotype','Or') | contains('Main.analysis')) %>%
              mutate(Type = 'Main') %>%
              rename('N' = 'Main.analysis_N',
                     'OR' = 'Main.analysis_OR',
                     'P.Value' = 'Main.analysis_P.Value',
                     'Lower' = 'Main.analysis_Lower',
                     'Upper' = 'Main.analysis_Upper')) %>%
  arrange(Or) %>%
  mutate(Type = factor(Type, levels = c('Validation','Main'))) 
bs_plot <- ggplot(bs,
                  aes(OR, fct_rev(fct_inorder(ID)))) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper, color = Type),
    position = position_dodge(width = dodge), width = w, cex = ww
  ) +
  geom_point(aes(color = Type),  size = 2, 
             position = position_dodge(width = dodge)) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 1) +
  scale_color_manual(values = c("skyblue2","black")) +
  theme_gray() +
  facet_wrap(~Phenotype, strip.position = 'top') +
  xlab('Odds Ratio') +
  theme(strip.text.x = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none") + 
  scale_x_log10(breaks = c(0.5,0.75,1,1.25,1.5), limits = c(0.45, 1.55))

leg <- data.frame(
  or = c(0.1,0.1), lower = c(0,0), upper = c(0.2, 0.2),
  y = c(1,2), 
  type = c('Main', 'Validation')
)

legq <- ggplot(leg, aes(or, y)) + 
  geom_errorbar(aes(xmin = lower, xmax = upper, colour = type),
                width = w, cex = ww) +
  geom_point(aes(colour = type), shape = 15, size = 2) +
  theme_bw()+
  scale_color_manual(values = c("black","skyblue2"))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

a <- ggplot() + theme_void() + ylim(-2,45) + xlim(-1.2,1) +
  annotation_custom(ggplotGrob(legq), # Legend
                    ymin = 30.5, ymax = 45.5,
                    xmin = -.7, xmax = 1) +
  annotation_custom(ggplotGrob(od_plot), # One dose
                    ymin = 24.5, ymax = 44,
                    xmin = -.925,xmax = +1) +
  annotation_custom(ggplotGrob(td_plot), # Two doses
                    ymin = 19.25, ymax = 25.1,
                    xmin = -.767,xmax = +1) +
  annotation_custom(ggplotGrob(bc_plot),
                    ymin = 1.5, ymax = 19.85,
                    xmin = -.75, xmax = 1) +
  annotation_custom(ggplotGrob(bs_plot), # Breakthrough severity
                    ymin = -3, ymax = 2.1,
                    xmin = -.75, xmax = 1) 
a
ggsave('a.png', plot = a, width = 25, height = 36, dpi = 300, units = 'cm')



