# ============================================================================ #
#                           CREATEFIGURE_VALIDATION                            #
#                            Marta Alcalde Herraiz                             #
# ============================================================================ #

gwas <- read.table(paste0(dir_results,'Validation/Validation.txt'))

gwas <- gwas %>%
  mutate(ID = if_else(ID == '6:32419074_CT_C', 'rs2150392827', ID)) %>%
  mutate(ID = if_else(ID == '6:32440321_CTG_C', 'rs17202941', ID)) %>%
  mutate(Or = c(9,10,4,5,1,2,6,7,3,8,12,11,21,19,18,20,14,13,16,15,17,22)) %>%
  arrange(Or)

dodge <- 0.5
w <- 0.3
ww <- 0.6
shape_main <- 15
shape_val  <- 16
c_val <- c('#72ACB9','#D57100','#EFBA00')

name <- c('One-dose antibody response',
          'Two-dose antibody response',
          'Breakthrough susceptibility',
          'Breakthrough severity')
plots <- list()
vali  <- list(c(2,2,2,2,2,2,2,2,1,1,1,1,2,2,1,1,1,1,1,1),
              c(3,3,2,2),
              c(2,2,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1),
              c(2,2))

for(i in c(1:length(name))){
  tab <- gwas %>% filter(Phenotype == name[i]) %>%
    select(c('CHROM','ID','Phenotype','Or') | contains('Validation')) %>%
    mutate(Type = 'Validation') %>%
    rename('N' = 'Validation_N',
           'OR' = 'Validation_OR',
           'P.Value' = 'Validation_P.Value',
           'Lower' = 'Validation_Lower',
           'Upper' = 'Validation_Upper') %>%
    full_join(gwas %>% filter(Phenotype == name[i]) %>%
                select(c('CHROM','ID','Phenotype','Or') | contains('Main.analysis')) %>%
                mutate(Type = 'Main') %>%
                rename('N' = 'Main.analysis_N',
                       'OR' = 'Main.analysis_OR',
                       'P.Value' = 'Main.analysis_P.Value',
                       'Lower' = 'Main.analysis_Lower',
                       'Upper' = 'Main.analysis_Upper')) %>%
    arrange(desc(Or)) %>%
    mutate(Type = factor(Type, levels = c('Validation','Main'))) 
  num <- nrow(tab)
  tab <- tab %>%
    mutate(order = rep(seq(1, num/2), each = 2)) %>%
    mutate(Vali = factor(vali[[i]], levels = c(1,2,3)))
  
  
  plots[[i]] <- ggplot(tab, aes(x = order,y = OR)) +
    geom_errorbar(
      aes(ymin = Lower, ymax = Upper, color = Vali, shape = Type),
      position = position_dodge(width = dodge), width = w, cex = ww
    ) +
    geom_point(aes(y = OR, color = Vali, shape = Type, size = Type),
               position = position_dodge(width = dodge)) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
    scale_shape_manual(values = c(shape_main,shape_val)) +
    scale_color_manual(values = c("1" = '#72ACB9',"2"='#D57100',"3"='#EFBA00')) +
    scale_size_manual(values = c(2.5,2.5)) +
    theme_gray() +
    facet_wrap(~Phenotype, strip.position = 'top') +
    theme(strip.text.x = element_text(size = 11, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(),
          axis.text.y  = element_text(size = 10),
          axis.title.y = element_blank(),
          legend.position = "none") +
    scale_x_continuous(
      breaks = tab$order,
      labels = tab$ID
    ) +
    scale_y_log10(breaks = c(0.5,0.75,1,1.25,1.5), limits = c(0.45, 1.55)) +
    coord_flip()
  
}


leg1 <- data.frame(
  x = c(1,2,3),
  y = c(2,3,4),
  type =  factor(c('Fully validated','Partially validated','Not validated'), levels = c('Fully validated','Partially validated','Not validated'))
)
leg2 <- data.frame(
  or = c(0.1,0.1,0.1), lower = c(0,0,0), upper = c(0.2, 0.2,.2),
  y = factor(c(3,2,1)), 
  type = factor(c('Main','Validation','Main'), levels = c('Main','Validation')),
  vali = factor(c('Fully validated','Partially validated','Not validated'), levels = c('Fully validated','Partially validated','Not validated'))
) 


legq1 <- ggplot(leg1, aes(x,y, fill = type)) +
  geom_bar(stat="identity") +
  theme_bw()+
  scale_fill_manual(values = c('Fully validated' = '#72ACB9',
                               'Partially validated' = '#D57100',
                               'Not validated' = '#EFBA00'))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11))

legq2 <- ggplot(leg2, aes(or, y)) +
  geom_errorbar(aes(xmin = lower, xmax = upper, shape = type),
                width = w, cex = ww) +
  geom_point(aes(shape = type, size = type)) +
  theme_bw()+
  scale_shape_manual(values = c('Main' = 16,'Validation' = 15)) +
  scale_size_manual(values = c('Main' = 5, 'Validation' = 5)) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 11))

a <- ggplot() + theme_void() + ylim(-1,45) + xlim(-1.2,1) +
  annotation_custom(ggplotGrob(legq1), # Legend
                    ymin = 40, ymax = 45.5,
                    xmin = -.7, xmax = 1) +
  annotation_custom(ggplotGrob(legq2),
                    ymin = 40.5, ymax = 43.5,
                    xmin = -.7, xmax = 1) +
  annotation_custom(ggplotGrob(plots[[1]]), # One dose
                    ymin = 22.35, ymax = 42,
                    xmin = -.9335,xmax = +1) +
  annotation_custom(ggplotGrob(plots[[2]]), # Two doses
                    ymin = 18.35, ymax = 23.35,
                    xmin = -.769,xmax = +1) +
  annotation_custom(ggplotGrob(plots[[3]]),
                    ymin = 1, ymax = 19.35,
                    xmin = -.75, xmax = 1) +
  annotation_custom(ggplotGrob(plots[[4]]), # Breakthrough severity
                    ymin = -2, ymax = 2.1,
                    xmin = -.75, xmax = 1)

ggsave(paste0(dir_results,'Figures/Validation.png'), plot = a, width = 25, height = 36, dpi = 300, units = 'cm')