

dir <- c('imputedData_oneDose',
         'imputedData_twoDose',
         'imputedData_breakthroughSusceptibility',
         'imputedData_breakthroughSeverity')
tit <- c('A) Immune response - One dose', 'B) Immune response - Two dose',
         'C) Breakthrough susceptibility', 'D) Breakthrough severity')
ylim_MH <- c(30,18,56,24)
ylim_QQ <- c(14,14,14,14)
gws_y    <- c(-log10(1.5e-8),-log10(2.5e-8),-log10(3e-8),-log10(1e-8))
y_QQ_min <- c(15,9,30,12)
y_QQ_max <- c(33,20,63,27)
alph     <- 0.7
#colors <- c('#79A7A5','#A3C6BC','#CEE5D4')

for (ii in 1:length(dir)){
  gwas <- read_delim(paste0(dir_results,'GWAS/',paste0(dir[ii],'.txt'))) %>%
    mutate(P = 10^(-LOG10P))  %>%
    select(CHR = CHROM, BP = GENPOS, SNP = ID, LOG10P, P) %>%
    filter(!is.na(LOG10P))
  
  # gwas <- read_delim(paste0(dir[ii],'.txt'))
  
  source(here('R_Scripts',"getManhattanPlot.R"))
  plot1 <- getManhattanPlot(gwas,y_lim = ylim_MH[ii])
  
  source(here('R_Scripts',"getQQPlot.R"))
  plot2 <- getQQPlot(gwas, x_lim = 6.5, y_lim = ylim_QQ[ii])
  
  gwas_top <- gwas %>% filter(LOG10P > 1)
  gwas_low <- gwas %>% filter(LOG10P <= 1) %>% sample_frac(0.01)
  gwas1 <- gwas_top %>% full_join(gwas_low) 
  
  don <- gwas1 %>%
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwas1, ., by = c("CHR" ="CHR")) %>%
    # Add the cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum = BP+tot) 
  axisdf <- don %>%
    group_by(CHR) %>%
    summarize(center=(max(BPcum) + min(BPcum) ) / 2 )
  
  if (ii == 1){
    # plot1<- plot1 + 
    #   annotate('text',don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
    #            don %>% summarise(LOG10P=max(LOG10P)) %>% as.numeric()+1, 
    #            label = '6:32419074_CT_C', 
    #            size = 1.8, color = 'black', hjust=0.5, vjust = 0) +
    #   annotate('rect',
    #            xmin = don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() -1e8,
    #            xmax = don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() +1e8,
    #            ymin = 8.25, ymax = 9.75,
    #            alpha = alph, fill = "white")+
    #   annotate('text',don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),9,
    #            label = 'rs9461694', 
    #            size = 1.8, color = 'black', hjust=0.5, vjust = 0.5) +
    #   ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 2){
    # plot1 <- plot1 + 
    #   annotate('text',don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
    #            don %>% summarise(LOG10P=max(LOG10P)) %>% as.numeric()+.5, 
    #            label = 'rs114903158', 
    #            size = 1.8, color = 'black', hjust=0.5, vjust = 0) +
    #   annotate('rect',
    #            xmin = don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() -1e8,
    #            xmax = don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() +1e8,
    #            ymin = 7.6, ymax = 8.4,
    #            alpha = alph, fill = "white")+
    #   annotate('text',don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
    #            8,
    #            label = 'rs3094106',
    #            size = 1.8, color = 'black', hjust = 0.5, vjust = 0.5) +
    #   ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 3){
    plot1 <- plot1 +
      annotate('text',don %>% filter(LOG10P==max(LOG10P) & CHR == 3) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs73062389',
               size=1.8, color = 'black', hjust=0.5, vjust=0)
    don <- don %>% filter(LOG10P != max(LOG10P))
    plot1 <- plot1 +
      annotate('text',don %>% filter(LOG10P==max(LOG10P) & CHR == 3) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 3) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs16861415',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',don %>% filter(CHR == 10) %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 10) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs1977830',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',don %>% filter(CHR == 19) %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
               don %>% filter(CHR == 19) %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs778809',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      annotate('text',don %>% filter(CHR == 19) %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
               18.3,
               label = 'rs11673136',
               size=1.8, color = 'black', hjust=0.5, vjust=0)
    a <- don %>% filter(CHR==19) %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric()+1e6
    don1 <- don %>% filter(CHR == 19) %>% filter(BPcum > a)
    plot1 <- plot1 +
      annotate('rect',
               xmin = don1 %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() -1e8,
               xmax = don1 %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric() +1e8,
               ymin = don1  %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+0.75, ymax = don1  %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.75,
               alpha = alph, fill = "white")+
      annotate('text',don1 %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
               don1  %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+1.25,
               label = 'rs681343',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0,ylim_MH[ii],5)) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
  }else if(ii == 4){
    plot1 <- plot1 +
      annotate('text',don %>% filter(LOG10P==max(LOG10P)) %>% select(BPcum) %>% as.numeric(),
               don %>% summarise(LOG10=max(LOG10P)) %>% as.numeric()+.5,
               label = 'rs62038344',
               size=1.8, color = 'black', hjust=0.5, vjust=0) +
      ggtitle(tit[ii]) + theme(plot.title = element_text(size = 7, face = "bold"))
    
  }
  
  # Merge both plots -----------------------------------------------------------
  p <- plot1 + annotation_custom(ggplotGrob(plot2),
                                 ymin = y_QQ_min[ii], ymax = y_QQ_max[ii],
                                 xmin = axisdf$center[axisdf$CHR == 12],
                                 xmax = don %>% filter(BPcum == max(BPcum)) %>% summarise(BPcum) %>% as.numeric()) 
  
  ggsave(paste0(dir_results,'Figures/',dir[ii],'.png'), plot = p, width = 165, height = 54, dpi = 300, units = 'mm')
}




