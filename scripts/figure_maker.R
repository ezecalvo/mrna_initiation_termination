library(tidyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(BiocManager)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(data.table)
library(cowplot)
library(png)
library(grid)
library(gridExtra)
library(ggdendro)
library(dendextend)
library(ape)
library(VennDiagram)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(scales)
library(ragg)
library(grid)
library(gridGraphics)
library(MetBrewer)
library(ggdraw)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

g.legend <- function(a.plot){
  tmp <- ggplot_gtable(ggplot_build(a.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

#load in all processed and analyzed data for plotting data from the TSS_TES_short_read_analysis_code.R
##### Figure 1
## ------------- FIGURE 1A: schematic of first and last exons ------------- ##
img <- image_read("/projectnb/evolution/carroll8/APRIL2025_TSS_TES_final_analyses/gtex_analysis/figure1_output-1.png")
plot(img)
img_grob <- rasterGrob(as.raster(img), interpolate = TRUE)
figure1A <- img_grob


## -----FIGURE 1B: total number of annotated first and last exons per gene---- ## 
figure1B <- ggplot(figure1bdata, aes(x = Var1, y = Freq, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "exons per gene", y = "number of genes", fill = "exon type") +  # custom title
  theme_classic(base_size = 12) +
  scale_y_continuous(
    breaks = seq(0, max(figure1bdata$Freq), by = 2000)
  ) +
  scale_fill_manual(
    values = c("FE" = "black", "LE" = "grey"),
    labels = c("FE" = "first", "LE" = "last")  # <- legend labels
  ) +
  theme(
    legend.position = c(0.85, 0.75),  # adjust to taste
    legend.justification = c("right", "top"),
    legend.direction = "vertical",   # stacked
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )
figure1B


## ------ FIGURE 1C: distribution of pearson's r vals of # first~last exons across all gtex samples ------ ##

figure1C <- ggplot(data.frame(figure1cdata), aes(x = figure1cdata)) +
  geom_density(fill = "#72bcd5", alpha = 0.6, color = "#1e466e", linewidth = 1) +
  geom_vline(xintercept = r_mean, linetype = "dashed", color = "black") +
  labs(x = "pearson r", y = "density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlim(-1, 1) +
  theme_classic(base_size = 14)  


##  ------ FIGURE 1C INSET: one representative ovary sample, plotting all genes by number of FE/LEs expressed  ------##

fig1c_inset_data

# Define values for plotting
r <- fig1c_correlation_line 
x_start <- 2
y_start <- 2
x_end <- 6
y_end <- y_start + r * (x_end - x_start)

figure1C_inset <- ggplot(fig1c_inset_data, aes(Freq.x, Freq.y)) +
  geom_point(aes(color = n), size = 3) +
  scale_color_gradientn(
    name = "gene count",
    # colors = c("#aadce0", "#72bcd5", "#528fad", "#1e466e", "navy"),
    colors = c(
      "#bde2eb",  # pale blue, 
      "#4a9fc1",  # strong blue
      "#1e466e",  # navy
      "#081d40"   # deep navy
    ),
    trans = "log",
    breaks = c(1, 10, 100, 1000),
    labels = scales::number_format(accuracy = 1),
    guide = guide_colorbar(
      barwidth = 2,barheight = 0.3,title.position = "top",title.theme = element_text(size = 3),label.theme = element_text(size = 3)
    )
  ) +
  geom_segment(x = x_start, y = y_start, xend = x_end, yend = y_end, color = "black") +
  coord_cartesian(xlim = c(1, 5), ylim = c(1, 5)) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(breaks = 1:5) +
  theme_cowplot() +
  labs(
    x = "number of first exons", y = "number of last exons") +
  theme(
    legend.position = "top",legend.justification = "center", legend.box = "horizontal",axis.title = element_text(size = 10)  # shrink both x and y titles
  )

figure1C_inset


## ------ FIGURE 1DA: Correlation in AFEPSI~ALEPSI of a given ordinal position in genes w 2+ AFEs and ALEs ------ ##

figure1da_data <- unlist(figure1da_data)
AFEorder <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
ALEorder <- c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)

figure1da_data_df  <- data.frame( AFEorder = AFEorder, ALEorder = ALEorder, value = figure1da_data)

figure1DA <- ggplot(figure1da_data_df, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    # colours = c( "#1e466e", "#528fad", "#72bcd5","white", "#f7aa58","#ef8a47","#e76254"),
    breaks = c(-0.5, -.25, 0,.25,  0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    name = ""
  )  +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 5), size = 12),
    axis.title.y = element_text(margin = margin(r = 5, b = -20), size = 12),
    axis.text.x = element_text(margin = margin(t = -10), size = 12),
    axis.text.y = element_text(margin = margin(r = -10), size = 12)
  ) +
  labs(
    x = "AFE Order",
    y = "ALE Order"
  )


## ------ FIGURE 1DA: Correlation in RELATIVE USAGE between first and last exons of a given ordinal position in genes using exactly 3 AFEs and 3 ALES ------ ##

data <- figure1db_data 
data <- as.vector(as.matrix(data))
AFEorder <- c(1,1,1,2,2,2,3,3,3)
ALEorder <- c(1,2,3,1,2,3,1,2,3)
AFE_ALE_pearson_1DB <- data.frame(AFEorder, ALEorder, data)

figure1DB <- ggplot(AFE_ALE_pearson_1DB, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = data)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    # colours = c( "#1e466e", "#528fad", "#72bcd5","white", "#f7aa58","#ef8a47","#e76254"),
    breaks = c(-0.5, -.25, 0,.25,  0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    # breaks = c(-1, -.5, 0,.5, 1),
    # labels = c(-1, -.5, 0,.5, 1),
    # limits = c(-1, 1),
    name = ""
  )  +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 5), size = 12),
    axis.title.y = element_text(margin = margin(r = 5, b = -20), size = 12),
    axis.text.x = element_text(margin = margin(t = -10), size = 12),
    axis.text.y = element_text(margin = margin(r = -10), size = 12)
  ) +
  labs(
    x = "",
    y = ""
  )
figure1DB


## ------------------------------ Figure 1D legend  ------------------------ ##

Figure1DB_withlegend <- ggplot(figure1DB, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = data)) + 
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.5, -.25, 0,.25,  0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    name = "pearson r"
  ) +  
  theme(plot.background = element_rect(fill = "white", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "white"),  axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12),
        axis.title.y = element_text(margin = margin(r = 10, b = -20), size = 12),
        axis.text.x = element_text(margin = margin(t = -10), size = 12),
        axis.text.y = element_text(margin = margin(r = -10), size = 12),
        legend.title = element_text(hjust = 0.5)) +  labs(x = "AFE Order", y = "ALE Order")

figure1D_legend <- get_legend(Figure1DB_withlegend)


## ------------------------ FULL FIGURE 1 ------------------------------ ##

Figure1 <-  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  draw_plot(figure1A, x = 0.01, y = .8, width = .95, height = .15) +
  draw_plot(figure1B, x = 0.01, y = .36, width = .35, height = .38) +
  draw_plot(figure1C, x = 0.35, y = .35, width = .64, height = 0.38) +
  draw_plot(figure1C_inset, x = .48, y = .52, width = .25, height = 0.25) +
  draw_plot(figure1DA, x = 0.02, y = 0.01, width = 0.40, height = 0.32) +
  draw_plot(figure1DB, x = .46, y = 0.01, width = 0.4, height = 0.32) +
  draw_plot(figure1D_legend, x = .88, y = 0.04, width = 0.1, height = 0.3)



##### Figure 2
dotted_LRS <- myo10_example[myo10_example$dataset=='reads_dotted',]
LRS_loop <- myo10_example[myo10_example$dataset=='reads',]
polyA_loop <- myo10_example[myo10_example$dataset=='polyA',]
first_exons_loop <- myo10_example[myo10_example$dataset=='FirstExons',]
isoforms_loop <- myo10_example[myo10_example$dataset=='isoforms',]
dotted_isoforms_loop <- myo10_example[myo10_example$dataset=='isoforms_dotted',]

id_add_isoforms <- min(-1*LRS_loop$id)*0.01+min(-1*LRS_loop$id)#Add a coefficient to separate the reads and make the plot more readable
fig_2_a <- ggplot()+#These reads are from ENCODE sample ENCFF142LPL
  geom_segment(dotted_LRS,mapping=aes(y=(-1*id),x=start,xend=end,yend=(-1*id)),linetype='solid',size=0.000001,color='lightgray')+
  geom_segment(LRS_loop,mapping=aes(x=start,xend=end,y=(-1*id),yend=(-1*id)),size=0.5,color='#876FAC')+
  theme_minimal()+
  labs(x='genomic coordinate on chr5 (Mb)')+
  geom_segment(first_exons_loop,mapping = aes(x=start,xend=end,y=id+15,yend=id+15),size=5,color='#2F308C')+
  geom_segment(polyA_loop,mapping = aes(x=start,xend=end,y=id+6,yend=id+6),size=5,color='#BB4D51')+
  geom_segment(dotted_isoforms_loop,mapping = aes(y=(id+id_add_isoforms+207)/2,x=start,xend=end,yend=(id+id_add_isoforms+207)/2),linetype='solid',size=0.005,color='black')+
  geom_segment(isoforms_loop,mapping = aes(x=start,xend=end,y=(id+id_add_isoforms+207)/2,yend=(id+id_add_isoforms+207)/2),size=1,color='#EA993E')+
  scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  theme(legend.position = "none",panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_text(size = 13, face = 'bold'))+
  geom_segment(aes(x=16670534,xend=16671466,y=-195,yend=-195))+
  geom_text(aes(x=16670534,y=-195),label='FERM')+
  geom_segment(aes(x=16673774,xend=16675052,y=-195,yend=-195))+
  geom_text(aes(x=16673774,y=-195),label='MyTH4')+
  geom_segment(aes(x=16704673,xend=16818095,y=-195,yend=-195))+
  geom_text(aes(x=16704673,y=-195),label='Myosin head')

fig_2_b <- ggplot(myo10_example_scatter,aes(x=start,y=end))+#The coordinates are flipped (start is the real start)
  geom_point(alpha=0.4,size=5,color='#876FAC',fill='#876FAC',shape=16)+
  theme_minimal(base_size = 15)+
  geom_smooth(method='lm',fill='lightgrey',color='black',linetype='dashed')+
  labs(x='TSS position on chr 5 (Mb)',y='TES position on chr 5 (Mb)')+
  scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  scale_y_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.line = element_line(colour = 'black', size = 1))


fig_2_e_df <- encode_spearman_per_sample_df[encode_spearman_per_sample_df$biosample%in%c("neural crest cell","lower lobe of right lung","mucosa of descending colon","heart right ventricle","ovary"),]
fig_2_e_df$exon_type <- factor(fig_2_e_df$exon_type,levels = c('soloT','dualT'))

fig_2_e <- fig_2_e_df %>% ggplot()+
  ggridges::geom_density_ridges(aes(x=corr,y=factor(biosample,levels = c("neural crest cell","lower lobe of right lung","mucosa of descending colon","heart right ventricle","ovary"),labels = c("Neural","Lung","Colon","Heart","Ovary")),fill=exon_type,alpha=exon_type),size=0.1,scale=1.3)+
  geom_vline(xintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 15)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  labs(x="Spearman's rho",y='Density',fill='Terminal exon usage')+
  theme(legend.position = "none",strip.text.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13))+
  scale_alpha_manual(values = c(0.5,0.8))+
  scale_x_continuous(breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))


fig_2_f_df <- encode_auc[encode_auc$biosample%in%c("neural crest cell","lower lobe of right lung","mucosa of descending colon","heart right ventricle","ovary"),]

fig_2_f_df$variable <- factor(fig_2_f_df$variable,levels = c('soloT','dualT'))
fig_2_f_df <- fig_2_f_df[fig_2_f_df$sample!='ENCFF187BTK',]


fig_2_f <-fig_2_f_df %>% 
  ggplot(aes(x = value, y = biosample, fill = variable)) +
  stat_summary(geom = "bar", fun = mean, position = position_dodge(width = 0.9), alpha = 0.8) +
  theme_minimal(base_size = 15) +
  scale_fill_manual(values = c('#A2A2AA', met.brewer('Hiroshige')[10])) +
  xlim(-0.3, 0.3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(colour = 'black', size = 1),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.title.y = element_blank()) +
  labs(x = 'ΔAUC R>0') +
  stat_summary(geom = "errorbar",fun.data = mean_se,position = position_dodge(width = 0.9),color = '#75726f',size = 0.2,width = 0.5)

encode_auc$variable <- factor(encode_auc$variable,levels = c('soloT','dualT'),labels = c('single FE\nor polyA','multiple FE\n& polyA'))

fig_2_g <- ggplot(encode_auc)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- encode_auc[encode_auc$variable=='soloT',] %>% pull(value)
y <- encode_auc[encode_auc$variable=='dualT',] %>% pull(value)
t.test(x,y)

#pdf("~/Desktop/fig2_Rout.pdf",width=8.5,height=11)
pdf("Figure2.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig_2_a, 0, 0.66, 1, 0.33)+
  draw_plot(fig_2_b, 0, 0.4, 0.33, 0.25)+
  draw_plot(fig_2_e, 0, 0, 0.33, 0.37)+
  draw_plot(fig_2_f, 0.33, 0, 0.33, 0.37)+
  draw_plot(fig_2_g, 0.66, 0, 0.33, 0.37)
dev.off()

##### Figure 3

fig3_b <-ggplot(ENSG00000007923_reads_acc_tissues)+
  geom_segment(aes(x=start_FE,xend=end_LE,y=id,yend=id,color=factor(Biosample.term.name,labels =  c('astrocyte','iPSC','lung'))),size=0.1)+
  theme_cowplot()+
  scale_color_manual(values = met.brewer(name = 'Cross')[c(1,3,6)])+
  labs(x='genomic coordinate on chr1 (Mb)',color='tissue',title = 'DNAJC11')+
  scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_text(size = 13, face = 'bold'),axis.line = element_blank())+
  geom_segment(dotted_isoforms_ENSG00000007923,mapping = aes(y=id,x=start,xend=end,yend=id),linetype='solid',size=0.05,color='black')+
  geom_segment(isoforms_loop_ENSG00000007923,mapping = aes(x=start,xend=end,y=id,yend=id),size=3,color='#EA993E')


df_3d$variable <- factor(df_3d$variable,levels = c('soloT','dualT'),labels = c('Single FE\nor polyA','Multiple FE\n& polyA'))
fig3_c <- ggplot(df_3d,aes(y=variable,x=value,fill=variable))+
  stat_summary(geom = "bar", fun = mean, position = position_dodge(width = 0.9), alpha = 0.8)+
  stat_summary(geom = "errorbar",fun.data = mean_se,position = position_dodge(width = 0.9),color = '#75726f',size = 0.2,width = 0.5)+
  theme_cowplot()+
  labs(x='ΔAUC R>0')+
  scale_fill_manual(values=c('#A1A1A1',met.brewer('Hiroshige')[10]))+
  geom_vline(xintercept = 0,linetype='dashed')+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=11), 
        axis.title = element_text(size=13),axis.title.y = element_blank())+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_x_continuous(labels = c(-0.5,0,0.5),breaks = c(-0.5,0,0.5),limits = c(-0.5,0.5))


sample_phast$region <- factor(sample_phast$region, levels=c("upstart","downstart","upend","downend"))
sample_phast$position <- factor(sample_phast$position, levels=c("upstream","downstream"))
sample_phast$terminal <- factor(sample_phast$terminal, levels=c("TSS","TES"))
sample_phast$type <- factor(sample_phast$type, labels = c("no PITA","PITA"))

fig3_d <-ggplot(sample_phast, aes(y=mean, x=factor(position), fill=factor(type))) + 
  geom_boxplot(notch=T,alpha=0.8) + 
  labs(x="region", y="mean phastCons score", fill="type") +
  facet_wrap(~terminal) +
  theme_cowplot()+
  scale_fill_manual(values = c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position=c(0.41,0.85), axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.title = element_blank(), legend.background = element_rect(fill="white", color="lightgrey"),
        strip.background = element_rect(fill="gray95"),axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), strip.text=element_text(size=10),axis.line = element_line(colour = 'black', size = 1))


df_3a$classif_corr <- factor(df_3a$classif_corr,levels = c('MiPITA','Non-PITA','MaPITA'),labels = c('anti-PITA','no PITA','PITA'))

fig3_e <- df_3a %>% subset(domain_change=='Different domains') %>% ggplot()+
  geom_boxplot(aes(x=classif_corr,y=proportion_of_genes*100,fill=classif_corr),position='dodge2',alpha=0.8,notch = T)+
  theme_cowplot()+
  labs(x='',y='% genes with isoform-specific protein domains',fill='')+
  #scale_fill_manual(values = c(met.brewer('Hiroshige')[10],'#A2A2AA','#727278'))+
  scale_fill_manual(values =c(met.brewer('Hiroshige')[2],'#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))+
  scale_y_continuous(labels = c(0,25,50,75,100),breaks = c(0,25,50,75,100))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#####
#CRISPRa data

df_3f$Exon <- factor(df_3f$Exon,
                     levels = rev(c("FE1", "FE2", "LE1", "LE2"))  # reverse for top-to-bottom order
)

fig3_f <- ggplot(df_3f,aes(x = logFC, y = Exon)) +
  geom_bar(stat = "identity",aes(fill = bar_color, color = bar_outline),size = 0) +
  geom_errorbar(aes(xmin = logFC - SE, xmax = logFC + SE), width = 0.2) +
  geom_text(aes(x = logFC + ifelse(logFC >= 0, SE + 0.1, -SE - 0.1),label = stars),
            hjust = ifelse(all_genes$logFC >= 0, 0, 1),size = 5) +
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_hline(yintercept = 0.4, color = "black", size = 1) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=11), 
        axis.title = element_text(size=13),axis.title.y = element_blank())+
  facet_grid(.~factor(gene,levels = c('ZNF638','MAST1','SWI5'),labels = c(c('ZNF638\nactivation','MAST1\nactivation','SWI5\ninactivation'))))+
  coord_cartesian(xlim = c(-1,1))

pdf("Figure3.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig3_b, 0, 0.66, 0.33, 0.33)+ 
  draw_plot(fig3_c, 0.33, 0.66, 0.33, 0.33)+
  draw_plot(fig3_d, 0.66, 0.66, 0.33, 0.33)+
  draw_plot(fig3_e, 0, 0.33, 1, 0.33)+
  draw_plot(fig3_f, 0, 0, 0.95, 0.25)
dev.off()


##### Figure 4
df_4a$classif_for_plot <- factor(df_4a$classif_for_plot,levels = c('<-0.1','-0.1-0.1','0.1-0.25','0.25-0.5','0.5-0.75','0.75-1'))

fig4_a <- df_4a[df_4a$variable%in%c("Upstream-Downstream","Upstream-Upstream","Downstream-Downstream"),] %>% subset(value>0) %>% ggplot()+
  geom_boxplot(aes(x=factor(variable,levels=c("Upstream-Downstream","Upstream-Upstream","Downstream-Downstream"),labels=c('gene','PITA pre-mRNA','PITA pre-mRNA')),y=value,fill=classif_for_plot),outlier.color = 'black',outlier.size = 0.5)+
  theme_cowplot()+
  scale_y_log10(labels = label_log())+
  scale_fill_manual(values=met.brewer('Hiroshige')[c(2,3,5,6,8,10)])+
  labs(x='',y='Distance (nt)',fill="Spearman's rho")+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

fig_4ab_legend <- g.legend(fig4_a)
fig4_a <- fig4_a+theme(legend.position = 'none')


fig4_b <- df_4a[!df_4a$variable%in%c("Upstream-Downstream","Upstream-Upstream","Downstream-Downstream"),] %>% subset(value>0) %>% subset(value>0) %>% ggplot()+
  geom_boxplot(aes(x=factor(variable,levels=c("FE_distance","Downstream-Upstream","polyA_distance"),labels=c('TSS interval','internal pre-mRNA','TES interval')),y=value,fill=classif_for_plot),outlier.color = 'black',outlier.size = 0.5)+
  theme_cowplot()+
  scale_y_log10(labels = label_log())+
  scale_fill_manual(values=met.brewer('Hiroshige')[c(2,3,5,6,8,10)])+
  labs(x='',y='Distance (nt)',fill="Spearman's rho")+
  theme(legend.position = 'none',
        axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))


df_4b[df_4b$corr>=0.2,'classif_for_plot_new'] <- 'MaPITA'

fig4_d <- ggplot(mapping=aes(x=FE_distance,y=polyA_distance))+
  geom_point(data=subset(df_4b, classif_for_plot_new == "Non-MaPITA"), color='#A1A1A1', size=0.3,alpha=0.1) +
  geom_point(data=subset(df_4b, classif_for_plot_new == "MaPITA"), color=met.brewer('Hiroshige')[10], size=0.3, alpha=0.1) +
  scale_x_log10(labels = label_log())+
  scale_y_log10(labels = label_log())+
  labs(y='Max distance between polyA',x='Max distance between FE') +
  theme_cowplot()+
  theme(legend.position = 'none',
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        #axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))

#Export 4d to make it easier to open it in AI
png("Figure4d.png",width=8.5,height=11,units = 'in',res = 300)
ggdraw() + 
  draw_plot(fig4_d, 0, 0.28, 0.33, 0.28)
dev.off()


df_4c$classif_final <- factor(df_4c$classif_final,levels = c('Non-PITA both','MaPITA mouse','MaPITA human','MaPITA both'),labels = c('no PITA','PITA\nmouse','PITA\nhuman','PITA\nmouse&human'))
df_4c <- df_4c[df_4c$classif_final!='PITA\nmouse&human',]

fig4_c <- df_4c[df_4c$variable%in%c('Gene length','FE','polyA'),] %>% ggplot()+
  geom_boxplot(aes(x=classif_final,y=delta_normalized), notch=T,outlier.size = 0.5)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  facet_wrap(factor(variable,levels=c('Gene length','FE','polyA'),labels = c('gene','TSS interval','TES interval'))~.)+
  labs(y='deltalength (human-mouse')+
  coord_cartesian(ylim = c(-2.5,2.5))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))

x <- df_4c[df_4c$classif_final=='no PITA' & df_4c$variable=='Gene length',] %>% pull(delta_normalized)
y <- df_4c[df_4c$classif_final=='PITA\nmouse' & df_4c$variable=='Gene length',] %>% pull(delta_normalized)

ks.test(x,y)

#Bring 4c as a png
image_grob <- rasterGrob(readPNG("Figure4c.png"), interpolate = TRUE)

#4E
sorted_intervals <- df_4e$max_dist_between_sites[order(as.numeric(gsub("\\[|,.*", "", df_4e$max_dist_between_sites)))]
sorted_intervals <- sorted_intervals[!duplicated(sorted_intervals)]
df_4e$max_dist_between_sites <- factor(df_4e$max_dist_between_sites, levels = sorted_intervals)
df_4e$deltaAUC <- df_4e$area_high-df_4e$area_low
df_4e$terminal_site <- factor(df_4e$terminal_site,levels = c('FE','PAS'))
df_4e$SEM_max <- df_4e$deltaAUC+df_4e$SEM
df_4e$SEM_min <- df_4e$deltaAUC-df_4e$SEM
df_4e$upper_lim_bin <- sub(')','',unlist(lapply(strsplit(as.character(df_4e$max_dist_between_sites),split = ','),'[',2))) %>% as.numeric()
labels_for_x_axis <- df_4e %>% filter(upper_lim_bin<6e+04) %>% pull(upper_lim_bin) %>% sort()%>% unique()
labels_for_x_axis <- labels_for_x_axis/1000

fig4_e <- df_4e %>% filter(number_of_genes>200) %>% 
  filter(upper_lim_bin<6e+04) %>% 
  ggplot()+
  theme_minimal(base_size = 15)+
  geom_col(aes(x=max_dist_between_sites,y=deltaAUC,fill=terminal_site),position = 'dodge2')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x = element_line(colour = 'black', size = 1),
        axis.line.y = element_line(colour = 'black', size = 1))+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(x='distance between alternative sites',y='ΔAUC',fill='terminal site')+
  geom_errorbar(aes(x=max_dist_between_sites, ymin = percentile_5,ymax = percentile_95,color=terminal_site),position = 'dodge2',size=0.2)+
  scale_fill_manual(values = met.brewer(name = 'Cassatt2')[c(4,7)],guide='none')+
  scale_color_manual(values=c('gray75','gray75'),guide='none')+
  scale_x_discrete(labels=labels_for_x_axis)+
  background_grid(major="none",minor="y")


pdf("Figure4.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig4_a, 0, 0.66, 0.33, 0.33)+
  draw_plot(fig4_b, 0.33, 0.66, 0.66, 0.33)+
  draw_plot(fig_4ab_legend,0.5,0.55,0.3,0.1)+
  draw_plot(fig4_c, 0.33, 0.25, 0.66, 0.28)+
  draw_plot(image_grob, 0, 0.28, 0.33, 0.28)+
  draw_plot(fig4_e, 0, 0, 0.4, 0.25)
dev.off()

##### Figure5



fig5a_left <- ggplot(polII_mut_wt, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = results_cor)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.6, -.3, 0,.3, 0.6),
    labels = c(-0.6, -.3, 0,.3, 0.6),
    limits = c(-0.6, 0.6),
    # breaks = c(-1, -.5, 0,.5, 1),
    # labels = c(-1, -.5, 0,.5, 1),
    # limits = c(-1, 1),
    name = "")  +
  theme(
    legend.position = 'none',
    axis.line = element_blank(),
    #panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 10),
    axis.ticks = element_blank(),
  )+
  scale_y_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  scale_x_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  labs(x='AFE ordinal position',y='ALE ordinal position')


fig5a_mid <- ggplot(polII_mut_fast, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = results_cor)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.6, -.3, 0,.3, 0.6),
    labels = c(-0.6, -.3, 0,.3, 0.6),
    limits = c(-0.6, 0.6),
    # breaks = c(-1, -.5, 0,.5, 1),
    # labels = c(-1, -.5, 0,.5, 1),
    # limits = c(-1, 1),
    name = "")  +
  theme(
    legend.position = 'none',
    axis.line = element_blank(),
    #panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 10),
    axis.ticks = element_blank(),
  )+
  scale_y_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  scale_x_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  labs(x='AFE ordinal position',y='ALE ordinal position')

fig5a_right <- ggplot(polII_mut_slow, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = results_cor)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.6, -.3, 0,.3, 0.6),
    labels = c(-0.6, -.3, 0,.3, 0.6),
    limits = c(-0.6, 0.6),
    # breaks = c(-1, -.5, 0,.5, 1),
    # labels = c(-1, -.5, 0,.5, 1),
    # limits = c(-1, 1),
    name = ""
  )  +
  theme(
    #legend.position = 'none',
    axis.line = element_blank(),
    #panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 10),
    axis.ticks = element_blank(),
  )+
  scale_y_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  scale_x_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  labs(x='AFE ordinal position',y='ALE ordinal position',fill=expression(atop("Pearson R",paste("AFE",Psi,"~ALE",Psi))))

fig5a_legend <- g.legend(fig5a_right)
fig5a_right <- fig5a_right+theme(legend.position = 'none')

#
read_coordinates$afe <- factor(read_coordinates$afe, levels=sort(unique(read_coordinates$afe)),labels = c(1,2))
#FE_coord_loop_parsed$feature <- factor(FE_coord_loop_parsed$feature, levels=sort(unique(FE_coord_loop_parsed$feature)),labels = c(1,2))
isoforms_loop_TSS_example <- isoforms_loop_TSS_example[isoforms_loop_TSS_example$id%in%c(10,11,13,14),]
dotted_isoforms_TSS_example <- dotted_isoforms_TSS_example[dotted_isoforms_TSS_example$id%in%c(10,11,13,14),]

fig5C <- ggplot(read_coordinates)+
  geom_segment(data = subset(read_coordinates, type=="intron"), mapping=aes(x=start, xend=end, y=-1*id, yend=-1*id),color="grey75", linetype="dashed", size=0.09)+
  geom_segment(data = subset(read_coordinates, type=="exon"), aes(x= start, xend=end, y=-1*id, yend=-1*id, color=factor(afe)), size=0.4)+
  #geom_segment(data=exon_coord_loop, aes(x= start , xend= end, y= 1, yend= 1), color="gray", size= 2.1)+
  #geom_segment(data=FE_coord_loop_parsed, aes(x= start , xend= end, y= 1, yend= 1, color=factor(feature)), size= 2.1) +
  #coord_cartesian(c(min(read_coordinates$start)-1000, max(read_coordinates$end)+1000))+
  labs(x = "position on chr19(Mb), HNRNPL, - strand", y = " ", color = "AFE ordinal position")+
  scale_color_manual(values=paletteer_d("PrettyCols::Summer"))+
  theme_cowplot() +
  #scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=9),axis.title=element_text(size=9),axis.line = element_blank(),legend.text=element_text(size=5),legend.title = element_text(size=6))+
  geom_segment(dotted_isoforms_TSS_example,mapping=aes(x=start,xend=end,y=-3*id+50, yend=-3*id+50),linetype='dashed',color='gray50',size=0.2)+
  geom_segment(isoforms_loop_TSS_example,mapping=aes(x=start,xend=end,y=-3*id+50, yend=-3*id+50),size=0.8)+
  scale_x_reverse(labels = unit_format(unit = '',scale = 1e-6))
fig5C_legend <- g.legend(fig5E)
fig5C <- fig5E+theme(legend.position = 'none')


#
fig5C_inset <- ggplot(subset(df5E_inset, afe != 4), aes(x=factor(afe), y = mapped_length,fill=factor(afe))) + 
  geom_boxplot(notch=T,outlier.size = 0.1,size=0.1) +
  labs(x="AFE ordinal position",y="genomic\ndistance (kb)") +
  theme_cowplot()+
  scale_fill_manual(values=paletteer_d("PrettyCols::Summer"))+
  theme(panel.grid = element_blank(),panel.background = element_blank(),axis.text=element_text(size=6),axis.title=element_text(size=6),legend.text=element_text(size=6),legend.title = element_text(size=9),legend.position = 'none',axis.title.x = element_blank())+
  scale_x_discrete(labels=c("1","2")) +
  scale_y_continuous(breaks = c(100,3000),labels = c('0.1','3'),trans='log10')

#


fig5D <- ggplot(bothreps_direct10strand_subbed_gene, aes(x=factor(afe_reorder_edit), y = mean_mapped, fill=factor(pita))) + 
  #geom_boxplot(notch=T) + 
  geom_split_violin(color=NA, alpha=0.5) +
  geom_boxplot(width=0.25, 
               position = position_dodge2(padding=0.4),
               notch=T, outlier.shape =  NA, coef=0, size=0.1, color="white",alpha=0.8) +
  #scale_fill_manual(values=c("dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3","dodgerblue4"),guide=F) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  scale_fill_manual(values=c("#A2A2AA","#1E466E")) +
  labs(x="AFE ordinal position",y="mean genomic distance traveled (per gene)") +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        panel.grid = element_blank(),panel.background = element_blank(),axis.text=element_text(size=9),axis.title=element_text(size=9),legend.text=element_text(size=9)
  )
#
fig5E <- ggplot(combo_smoothed, aes(x=bin3, y=value, ymin=qinf, ymax=qsup))+
  geom_ribbon(aes(fill =factor(group)), alpha=0.1) +
  geom_line(aes(color=factor(group)), size=2,alpha=0.8) +
  # blocking off intermediate regions
  annotate("rect", xmin = 100, xmax = 100+offset, ymin = -Inf, ymax = Inf,fill= "white")  +
  annotate("rect", xmin = max(E_smoothed$bin2), xmax = combo_offset, ymin = -Inf, ymax = Inf,fill= "white")  +
  annotate("rect", xmin = 100+combo_offset-3, xmax = 100+combo_offset+offset+3, ymin = -Inf, ymax = Inf,fill= "white")  +
  ### adding exon boxes and lines at beginning
  annotate("segment", x = 50, xend = max(E_smoothed$bin2), y = -0.375, yend = -0.375,color= "lightgrey", linetype="solid")  +
  annotate("segment", x = combo_offset, xend = combo_offset+150+offset, y = -0.375, yend = -0.375,color= "lightgrey", linetype="solid")  +
  annotate("segment", x = max(E_smoothed$bin2)-1, xend = max(E_smoothed$bin2)+1, y = -0.475, yend = -0.275,color= "lightgrey", linetype="solid")  +
  annotate("segment", x = combo_offset-1, xend = combo_offset+1, y = -0.475, yend = -0.275,color= "lightgrey", linetype="solid")  +
  # TSS1
  annotate("segment", x = 50, xend = 50, y = -0.25, yend = 3,color= "lightgrey", linetype="dotted")  +
  annotate("segment", x = 50, xend = 50, y = -0.26, yend = -0.17,color= "grey25", linetype="solid")  +
  annotate("segment", x = 49.5, xend = 54, y = -0.17, yend = -0.17,color= "grey25", linetype="solid", arrow = arrow(type = "closed", length = unit(0.01, "npc")))  +
  annotate("rect", xmin = 50, xmax = 70, ymin = -0.5, ymax = -0.25,color= "grey50",fill="grey50")  +
  # TSS2
  annotate("segment", x = 150+offset, xend = 150+offset, y = -0.25, yend = 3,color= "lightgrey", linetype="dotted")  +
  annotate("segment", x = 150+offset, xend = 150+offset, y = -0.26, yend = -0.17,color= "grey25", linetype="solid")  +
  annotate("segment", x = 149.5+offset, xend = 154+offset, y = -0.17, yend = -0.17,color= "grey25", linetype="solid", arrow = arrow(type = "closed", length = unit(0.01, "npc")))  +
  annotate("rect", xmin = 150+offset, xmax = 170+offset, ymin = -0.5, ymax = -0.25,color= "grey50",fill="grey50")  +
  # PAS1
  annotate("segment", x = combo_offset+50, xend = combo_offset+50, y = -0.25, yend = 3,color= "lightgrey", linetype="dotted")  +
  annotate("point", x = combo_offset+50, y = -0.20, colour = "grey25", size = 1, shape=25,fill='grey25') +
  annotate("rect", xmin = combo_offset+50-20, xmax = combo_offset+50, ymin = -0.5, ymax = -0.25,color= "grey50",fill="grey50")  +
  # PAS2
  annotate("segment", x = combo_offset+150+offset, xend = combo_offset+150+offset, y = -0.25, yend = 3,color= "lightgrey", linetype="dotted")  +
  annotate("point", x = combo_offset+150+offset, y = -0.20, colour = "grey25", size = 1, shape=25,fill='grey25') +
  annotate("rect", xmin = combo_offset+150+offset-20, xmax = combo_offset+150+offset, ymin = -0.5, ymax = -0.25,color= "grey50",fill="grey50")  +
  # add text
  annotate("text", label=combo_labels_text, y = -0.375, color="white",
           x = c(50+10, x = 150+offset+10, combo_offset+50-10, combo_offset+150+offset-10), hjust=0.5,size=2.5)  +
  # attributes
  labs(y="elongation velocity",x="Distance to TSS",color="")+
  scale_fill_manual(values=c("#A2A2AA","#1E466E"),labels=c("no PITA","PITA"),guide=F) +
  scale_color_manual(values=c("#A2A2AA","#1E466E"),labels=c("no PITA","PITA")) +
  scale_x_continuous(breaks = combo_breaks, labels =combo_labels) +
  #facet_grid(~region) +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text=element_text(size=9),axis.title=element_text(size=10),
        legend.position = c(.87,.83))

#
#pairingschematic_grob <- rasterGrob('~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/LRSvelocity/final_sch.png', interpolate=T)
pairingschematic_grob <- rasterGrob(readPNG('~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/LRSvelocity/nascentelongation_schematic_final.png'), interpolate = TRUE)


pdf("Figure5.pdf",width=8.5,height=11)

ggdraw() + 
  draw_plot(fig5a_left, 0, 0.75, 0.3, 0.25)+
  draw_plot(fig5a_mid, 0.3, 0.75, 0.3, 0.25)+
  draw_plot(fig5a_right, 0.6, 0.75, 0.3, 0.25)+
  draw_plot(fig5a_legend, 0.93, 0.85, 0.05, 0.1)+
  draw_plot(pairingschematic_grob,0,0.585,0.66,0.165)+
  draw_plot(fig5C, 0, 0.415, 0.66, 0.165)+
  draw_plot(fig5C_legend, 0.05, 0.45, 0.1, 0.1)+
  draw_plot(fig5C_inset, 0.45, 0.45, 0.2, 0.1)+
  draw_plot(fig5D, 0.66, 0.415, 0.33, 0.33)+
  draw_plot(fig5E, 0, 0.1, 1, 0.33)+
  draw_plot_label(label=c("A","B","C","D","E"), x=c(0, 0, 0, 0.66, 0), y=c(1, 0.75, 0.585, 0.75, 0.415),family="Helvetica", color="black")
dev.off()


################
#Supp figures
#####extended figure 1A ####
# Save as a grob for cowplot
figureS1A <- grid.grabExpr({
  draw.pairwise.venn(
    area1 = 577, #these come from figures1adata
    area2 = 309,  #these come from figures1adata
    cross.area = 108, #these come from figures1adata
    category = c("genes expressing AFEs", "genes expressing ALEs"),
    fill = c("gray", "gray76"),
    cat.pos = c(0, 0),
    fontfamily = rep("sans", 3),         
    cat.fontface = rep("plain", 2),
    cat.fontfamily = rep("sans", 2)
  )
})

figureS1A


## ------------ FIGURE S1B : Mean Pearson's r by tissue type ------------  ##

figures1b_data

figureS1B <- ggplot(figures1b_data, aes(x = reorder(category, mean_correlation), y = mean_correlation)) +
  geom_col(fill = "#72bcd5") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  labs(
    title = "",
    x = "tissue category",
    y = "mean pearson r"
  )
)

## ------------ FIGURE S1C: Max AFEPSI and Max ALEPSI per gene  ------------  ##
figures1cdata
# Generate and capture heatmap as grob
psi_pheatmap <- pheatmap(
  figures1cdata,
  color = colorRampPalette(c("white", "#72bcd5", "#528fad", "#1e466e"))(100),
  breaks = c(seq(0, 0.94, length.out = 80), seq(0.95, 1, length.out = 20)),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  silent = TRUE  
)

png("/path/tosaving/data/as/figure/psi_pheatmap",
    width = 2400, height = 3200, res = 400)
grid::grid.draw(psi_pheatmap$gtable)
dev.off()

figureS1C <- ggdraw() +
  draw_image("/path/tosaving/data/as/figure/psi_pheatmap.png")


## ------------------------ FULL SUPPPLEMENTAL FIGURE 1 ------------------------------ ##

figureS1 <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  draw_plot(figureS1A, x = 0.00, y = 0.66, width = 0.33, height = 0.33) +
  draw_plot(figureS1B, x = 0.33, y = 0.66, width = 0.33, height = 0.33) +
  draw_plot(figureS1C,  x = 0.66, y = 0.66, width = 0.33, height = 0.33)

#Schematics fig 1
schematic_intra_mol_fig <- schematic_intra_mol %>% ggplot()+
  geom_segment(aes(x=start,xend=end,y=id,yend=id,color=isoform),size=0.2)+
  theme_cowplot()+
  scale_color_manual(values = met.brewer('Hiroshige',n=50)[c(30,3,10,40)])+
  theme(legend.position = "none",panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank(),axis.line = element_blank(),axis.title = element_blank())

schematic_inter_mol_fig <- schematic_inter_mol %>% ggplot()+
  geom_segment(aes(x=start,xend=end,y=id,yend=id,color=isoform),size=0.2)+
  theme_cowplot()+
  scale_color_manual(values = met.brewer('Hiroshige',n=50)[c(30,3,10,40)])+
  theme(legend.position = "none",panel.grid = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank(),axis.line = element_blank(),axis.title = element_blank())

pdf("inter_intra_cor.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(schematic_intra_mol_fig, 0, 0.33, 0.5, 0.33)+
  draw_plot(schematic_inter_mol_fig, 0.5, 0.33, 0.5, 0.33)
dev.off()

#### Figure S2

dotted_LRS <- myo10_example_supp_figure[myo10_example_supp_figure$dataset=='reads_dotted',]
LRS_loop <- myo10_example_supp_figure[myo10_example_supp_figure$dataset%in%c("Passed","No FE","No polyA","Same exon","No FE & polyA"),]
polyA_loop <- myo10_example_supp_figure[myo10_example_supp_figure$dataset=='polyA',]
first_exons_loop <- myo10_example_supp_figure[myo10_example_supp_figure$dataset=='FirstExons',]

dotted_LRS$id <- as.numeric(dotted_LRS$id)
exons_gene$id <- as.numeric(exons_gene$id)
polyA_loop$id <- as.numeric(polyA_loop$id)
first_exons_loop$id <- as.numeric(first_exons_loop$id)

LRS_loop$id <- as.numeric(LRS_loop$id)

LRS_loop$dataset <- factor(LRS_loop$dataset,levels = c("Passed","No FE","No polyA","Same exon","No FE & polyA"),labels = c("passed","no FE","no polyA","same exon","no FE & polyA"))

pdf("FigureS2.pdf",width=8.5,height=11)
print(ggplot(LRS_loop,mapping=aes(x=start,xend=end,y=(-1*id),yend=(-1*id)))+
        geom_segment(dotted_LRS,mapping=aes(y=(-1*id),x=start,xend=end,yend=(-1*id)),linetype='dashed',size=0.1)+
        geom_segment(size=3.5,aes(color=dataset))+
        scale_color_manual(values = c('#876FAC','#AFB86E','#74C2B9','#23C373','#BDB042'))+
        theme_minimal()+
        labs(x='genomic coordinate on chr5 (Mb)',title = 'MYO10')+
        geom_segment(first_exons_loop,mapping = aes(x=start,xend=end,y=id+9,yend=id+9),size=5,color='#2F308C')+
        geom_segment(polyA_loop,mapping = aes(x=start,xend=end,y=id+6,yend=id+6),size=5,color='#BB4D51')+
        scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
        theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_text(size = 13, face = 'bold')))
dev.off()

#### Figure S3
#2.1
figs2_1b <- ggplot(df_sup_2.1a)+
  geom_boxplot(aes(x=as.factor(sampling_factor),y=auc_sampling,fill=factor(variable,levels = c('soloT','dualT'))),alpha=0.8)+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(x='proportion of sampled library',y='ΔAUC R>0',fill='Exon type',title = '')+
  scale_x_discrete(breaks=seq(0,1,0.1))+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))


figs2_1c <- ggplot(df_sup_2.1b)+
  geom_boxplot(aes(x=factor(variable,levels = c('soloT','dualT')),y=value,fill=factor(variable,levels = c('soloT','dualT'))),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- df_sup_2.1b[df_sup_2.1b$variable=='soloT',] %>% pull(value)
y <- df_sup_2.1b[df_sup_2.1b$variable=='dualT',] %>% pull(value)
t.test(x,y)

#
df_sup_2.1_isoform_threshold$variable <- factor(df_sup_2.1_isoform_threshold$variable,levels=c('soloT','dualT'))

figs2_1d <- ggplot(df_sup_2.1_isoform_threshold)+
  geom_boxplot(aes(x=as.factor(isoform_exp_filter),y=value,fill=variable),alpha=0.8,outlier.size = 0.1,notch=T)+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y='ΔAUC',fill='Exon type',x='number of reads per isoform')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

figs2_1e <- ggplot(df_sup_2.1c)+
  geom_boxplot(aes(x=factor(bin_deciles,levels=c('[93,1.08e+03]','(1.08e+03,1.4e+03]','(1.4e+03,1.66e+03]','(1.66e+03,1.89e+03]','(1.89e+03,2.12e+03]','(2.12e+03,2.38e+03]','(2.38e+03,2.69e+03]','(2.69e+03,3.09e+03]','(3.09e+03,3.73e+03]','(3.73e+03,1.38e+04]'),labels = c('[93\n1.08e+03]','(1.08e+03\n1.4e+03]','(1.4e+03\n1.66e+03]','(1.66e+03\n1.89e+03]','(1.89e+03\n2.12e+03]','(2.12e+03\n2.38e+03]','(2.38e+03\n2.69e+03]','(2.69e+03\n3.09e+03]','(3.09e+03\n3.73e+03]','(3.73e+03\n1.38e+04]')),y=corr))+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y="Spearman's rho",x='read length')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))

figs2_1f <- ggplot(df_sup_2.1d,aes(x=factor(variable,levels = c("Unique FE or LE","Alternative FE & LE"),labels = c('soloT','dualT')),y=value,fill=factor(variable,levels = c("Unique FE or LE","Alternative FE & LE"),labels = c('soloT','dualT'))))+
  stat_summary(geom = "bar", fun = mean, position = position_dodge(width = 0.9), alpha = 0.8)+
  theme_cowplot()+
  scale_fill_manual(values=c('#808080',met.brewer('Hokusai3',n = 10)[10]))+
  stat_summary(geom = "errorbar",fun.data = mean_se,position = position_dodge(width = 0.9),color = '#75726f',size = 0.2,width = 0.5)+
  ylim(-0.2,0.2)+
  geom_vline(xintercept = 0,linetype='dashed')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))+
  labs(y='ΔAUC R>0',x='')


figs2_1g <- ggplot(df_sup_2.1f)+
  geom_density(aes(x=corr,fill=factor(exon_type,levels = c("soloT","dualT"),labels = c('soloT','dualT'))),alpha=0.8)+
  scale_fill_manual(values=c('#808080',met.brewer('Hokusai3',n = 10)[10]))+
  theme_cowplot()+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))+
  labs(x="Pearson's R")+
  xlim(-0.5,0.5)+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())

figs2_1h <- ggplot(df_sup_2.1e)+
  geom_boxplot(aes(x=factor(variable,levels = c('soloT','dualT')),y=value,fill=factor(variable,levels = c('soloT','dualT'))),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

figs2_1i <- ggplot(df_sup_2.1g)+
  geom_boxplot(aes(x=as.factor(threshold),y=value,fill=factor(variable,levels = c('soloT','dualT'))),alpha=0.8)+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y='ΔAUC',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

figs2_1j <- ggplot(df_sup_2.1h)+
  geom_boxplot(aes(x=same_different,y=value,fill=same_different))+
  theme_cowplot()+
  scale_fill_manual(values=c("#d4a0a4","#fbdcb1"))+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  labs(y='correlation between samples')

pdf("FigureS3.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs2_1b, 0.35, 0.66, 0.35, 0.33)+
  draw_plot(figs2_1c, 0.7, 0.66, 0.3, 0.33)+
  draw_plot(figs2_1d, 0, 0.33, 0.4, 0.33)+
  draw_plot(figs2_1f, 0.4, 0.33, 0.2, 0.33)+
  draw_plot(figs2_1e, 0.6, 0.33, 0.4, 0.33)+
  draw_plot(figs2_1g, 0, 0, 0.2, 0.33)+
  draw_plot(figs2_1h, 0.2, 0, 0.2, 0.33)+
  draw_plot(figs2_1i, 0.4, 0, 0.4, 0.33)+
  draw_plot(figs2_1j, 0.8, 0, 0.2, 0.33)
dev.off()

#### Figure S4
figs2_2a <- ggplot(df_sup2.2a)+
  geom_boxplot(aes(x=tech,y=value,fill=variable),alpha=0.8,outlier.size = 0.1)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  labs(x='',y='ΔAUC',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  ylim(-0.25,0.5)


x <- df_sup2.2a[df_sup2.2a$tech=='TIF-Seq2' & df_sup2.2a$variable=='soloT',] %>% pull(value)
y <- df_sup2.2a[df_sup2.2a$tech=='TIF-Seq2' & df_sup2.2a$variable=='dualT',] %>% pull(value)
t.test(x,y) 

x <- df_sup2.2a[df_sup2.2a$tech=='FLAM-seq' & df_sup2.2a$variable=='soloT',] %>% pull(value)
y <- df_sup2.2a[df_sup2.2a$tech=='FLAM-seq' & df_sup2.2a$variable=='dualT',] %>% pull(value)
t.test(x,y) 

figs2_2b <- ggplot(df_sup2.2b)+ #Our data through their pipeline but ours for exon classification
  geom_bar(aes(x=predominance_new,fill=factor(exon_type,levels = c('soloT','dualT'),labels = c('soloT','dualT'))),alpha=0.8)+
  theme_cowplot()+
  labs(fill='PiPITA exon\nclassification',x='LATER predominance')+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))


figs2_2c <- subset(df_sup2.2c,Var1=='PITA') %>% ggplot()+
  geom_boxplot(aes(x=Var2,y=mapita_norm))+
  theme_cowplot()+
  labs(x='',y='% of PITA genes')+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))

figs2_2d <- subset(df_sup2.2c,Var2=='predominant') %>% ggplot()+
  geom_boxplot(aes(x=Var1,y=dominant_norm))+
  theme_cowplot()+
  labs(x='',y='% of predominant genes')+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))


pdf("FigureS3.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs2_2a, 0, 0.66, 0.5, 0.33)+
  draw_plot(figs2_2b, 0.5, 0.66, 0.5, 0.33)+
  draw_plot(figs2_2c, 0, 0.33, 0.5, 0.33)+
  draw_plot(figs2_2d, 0.5, 0.33, 0.5, 0.33)
dev.off()



#### Figure S5
df_sup3.1a$classif <- factor(df_sup3.1a$classif,levels = c('MiPITA','Non-PITA','MaPITA'),labels = c('anti-PITA','no PITA','PITA'))
figs3_1a <- ggplot(df_sup3.1a)+
  geom_density(aes(x=n_rows,color=classif),size=2)+
  theme_cowplot()+
  xlim(0,15)+
  labs(x='Number of total protein domains per gene')+
  scale_color_manual(values =c(met.brewer('Hiroshige')[2],'#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',
        axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))

figs3_1_norm_exon_length$classif <- factor(figs3_1_norm_exon_length$classif,levels = c('MiPITA','Non-PITA','MaPITA'),labels = c('anti-PITA','no PITA','PITA'))

figs3_1_norm_exon_length <- ggplot(domains_length_norm)+
  geom_density(aes(x=domains_per_nt,color=classif),size=2)+
  theme_cowplot()+
  labs(x='Domains per exon nt')+
  scale_color_manual(values =c(met.brewer('Hiroshige')[2],'#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',
        axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))+
  xlim(0,0.01)

df_sup3.1b$classif_corr <- factor(df_sup3.1b$classif_corr,levels = c('MiPITA','Non-PITA','MaPITA'),labels = c('anti-PITA','no PITA','PITA'))

figs3_1b <- df_sup3.1b[df_sup3.1b$domain_change=="Different domains",] %>% ggplot(aes(x=classif_corr,y=proportion_of_genes,fill=classif_corr))+
  stat_summary(geom = "errorbar",fun.data = mean_se,position = position_dodge(width = 0.9),color = '#75726f',size = 0.2,width = 0.5)+
  stat_summary(geom = "bar", fun = mean, position = position_dodge(width = 0.9), alpha = 0.8)+
  theme_cowplot()+
  labs(y='% genes with isoform−specific\nprotein domains across tissues',fill='')+
  scale_fill_manual(values =c(met.brewer('Hiroshige')[2],'#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(colour = 'black', size = 1),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13))

fig_s3_legend <- g.legend(figs3_1b)
figs3_1b <- figs3_1b+theme(legend.position = 'none')


tissue_phast$region <- factor(tissue_phast$region, levels=c("upstart","downstart","upend","downend"))
tissue_phast$position <- factor(tissue_phast$position, levels=c("upstream","downstream"))
tissue_phast$terminal <- factor(tissue_phast$terminal, levels=c("TSS","TES"))
tissue_phast$type <- factor(tissue_phast$type, levels=c('nonPITA','PITA'),labels = c('no PITA','PITA'))

figs3_1c <- ggplot(tissue_phast, aes(y=mean, x=factor(position), fill=factor(type))) + 
  geom_boxplot(notch=T,alpha=0.8) + 
  labs(x="region", y="mean phastCons score", fill="type") +
  theme_cowplot()+
  facet_wrap(~terminal) +
  theme(legend.position=c(0.41,0.85), axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.title = element_blank(), legend.background = element_rect(fill="white", color="lightgrey"),
        strip.background = element_rect(fill="gray95"),axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), strip.text=element_text(size=10),axis.line = element_line(colour = 'black', size = 1))+
  scale_fill_manual(values = c('#A2A2AA',met.brewer('Hiroshige')[10]))


#When adding the exon length norm I remade everything manually on AI
pdf("FigureS5.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs3_1a, 0.3, 0.66, 0.35, 0.33)+
  draw_plot(figs3_1_norm_exon_length, 0.65, 0.66, 0.35, 0.33)+
  draw_plot(figs3_1c, 0, 0.66, 0.3, 0.33)
dev.off()

#### Figure S6

dfS6_qPCR <- rbind(TP53_dataframe[,-c(6,8)],SYT9_dataframe[,-c(6,8)],LEPR_dataframe)
dfS6_exon_coordinates <- exon_coordinates[exon_coordinates$gene%in%c('TP53','SYT9','LEPR'),]


tp53_qpcr <- dfS6_qPCR %>% filter(gene=='TP53') %>% 
  ggplot(aes(y = factor(Exon,levels = c('LE2','LE1','FE2','FE1')), x = Relative_Change)) +
  geom_bar(aes(fill = fill_color), stat = "identity", color = NA) +
  geom_errorbar(aes(xmin = Relative_Change - SEM, xmax = Relative_Change + SEM), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-1,1))+
  labs(y = "", x = "") +
  theme_classic() +
  theme(
    axis.line.x = element_line(size = 1),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    #axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=10),
    strip.text = element_blank(),
    legend.position = "none"
  )+
  #facet_grid(.~factor(gene,levels = c('ZNF638','MAST1','SWI5'),labels = c(c('ZNF638\nactivation','MAST1\nactivation','SWI5\ninactivation'))),scales = 'free')+
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_hline(yintercept = 0.4, color = "black", size = 1)

tp53_exon <- dfS6_exon_coordinates %>% filter(gene=='TP53') %>% 
  ggplot()+
  geom_hline(yintercept = 1,linetype='dashed',color='lightgray')+
  geom_segment(mapping=aes(x=start,xend=end,y=1,yend=1,color=factor(exon_type)),size=10)+
  
  theme_cowplot()+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )+
  scale_color_manual(values = c('#72bcd5','#A2A2AA','gray90','#72bcd5','#A2A2AA'))


#
syt9_qpcr <- dfS6_qPCR %>% filter(gene=='SYT9') %>% 
  ggplot(aes(y = factor(Exon,levels = c('LE2','LE1','FE2','FE1')), x = Relative_Change)) +
  geom_bar(aes(fill = fill_color), stat = "identity", color = NA) +
  geom_errorbar(aes(xmin = Relative_Change - SEM, xmax = Relative_Change + SEM), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-50,50))+
  labs(y = "", x = "") +
  theme_classic() +
  theme(
    axis.line.x = element_line(size = 1),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    #axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=10),
    strip.text = element_blank(),
    legend.position = "none"
  )+
  #facet_grid(.~factor(gene,levels = c('ZNF638','MAST1','SWI5'),labels = c(c('ZNF638\nactivation','MAST1\nactivation','SWI5\ninactivation'))),scales = 'free')+
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_hline(yintercept = 0.4, color = "black", size = 1)

syt9_exon <- dfS6_exon_coordinates %>% filter(gene=='SYT9') %>% 
  ggplot()+
  geom_hline(yintercept = 1,linetype='dashed',color='lightgray')+
  geom_segment(mapping=aes(x=start,xend=end,y=1,yend=1,color=factor(exon_type)),size=10)+
  
  theme_cowplot()+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )+
  scale_color_manual(values = c('#A2A2AA','#72bcd5','gray90','#A2A2AA','#72bcd5'))

#
LEPR_qpcr <- dfS6_qPCR %>% filter(gene=='LEPR') %>% 
  ggplot(aes(y = factor(Exon,levels = c('LE2','LE1','FE2','FE1')), x = Relative_Change)) +
  geom_bar(aes(fill = fill_color), stat = "identity", color = NA) +
  geom_errorbar(aes(xmin = Relative_Change - SEM, xmax = Relative_Change + SEM), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-42,42))+
  labs(y = "", x = "") +
  theme_classic() +
  theme(
    axis.line.x = element_line(size = 1),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    #axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=10),
    strip.text = element_blank(),
    legend.position = "none"
  )+
  #facet_grid(.~factor(gene,levels = c('ZNF638','MAST1','SWI5'),labels = c(c('ZNF638\nactivation','MAST1\nactivation','SWI5\ninactivation'))),scales = 'free')+
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_hline(yintercept = 0.4, color = "black", size = 1)

LEPR_exon <- dfS6_exon_coordinates %>% filter(gene=='LEPR') %>% 
  ggplot()+
  geom_hline(yintercept = 1,linetype='dashed',color='lightgray')+
  geom_segment(mapping=aes(x=start,xend=end,y=1,yend=1,color=factor(exon_type)),size=10)+
  
  theme_cowplot()+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )+
  scale_color_manual(values = c('#A2A2AA','#72bcd5','gray90','#A2A2AA','#72bcd5'))


znf_qpcr <- df3F_qPCR %>% filter(gene=='ZNF638') %>% 
  ggplot(aes(y = factor(Exon,levels = c('LE2','LE1','FE2','FE1')), x = Relative_Change)) +
  geom_bar(aes(fill = fill_color), stat = "identity", color = NA) +
  geom_errorbar(aes(xmin = Relative_Change - SEM, xmax = Relative_Change + SEM), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-.25,.25))+
  labs(y = "", x = "") +
  theme_classic() +
  theme(
    axis.line.x = element_line(size = 1),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    #axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=10),
    strip.text = element_blank(),
    legend.position = "none"
  )+
  #facet_grid(.~factor(gene,levels = c('ZNF638','MAST1','SWI5'),labels = c(c('ZNF638\nactivation','MAST1\nactivation','SWI5\ninactivation'))),scales = 'free')+
  geom_vline(xintercept = 0, color = "black", size = .5) +
  geom_hline(yintercept = 0.4, color = "black", size = 1)

znf_exon <- df3F_exon_coordinates %>% filter(gene=='ZNF638') %>% 
  ggplot()+
  geom_hline(yintercept = 1,linetype='dashed',color='lightgray')+
  geom_segment(mapping=aes(x=start,xend=end,y=1,yend=1,color=factor(exon_type)),size=10)+
  
  theme_cowplot()+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )+
  scale_color_manual(values = c('#72bcd5','#A2A2AA','gray90','#72bcd5','#A2A2AA'))

pdf("FigureS6.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(tp53_qpcr, 0, 0, 0.33, 0.25)+
  draw_plot(tp53_exon, 0, 0.15, 0.33, 0.25)+
  draw_plot(syt9_qpcr, 0.33, 0, 0.33, 0.25)+
  draw_plot(syt9_exon, 0.33, 0.15, 0.33, 0.25)+
  draw_plot(LEPR_qpcr, 0.66, 0, 0.33, 0.25)+
  draw_plot(LEPR_exon, 0.66, 0.15, 0.33, 0.25)+
  draw_plot(znf_qpcr, 0.5, 0.66, 0.33, 0.25)+
  draw_plot(znf_exon, 0.66, 0.43, 0.33, 0.25)
dev.off()

#### Figure S7
#Expanded cor <-0.1 for for human
distances_full_neg_values$classif_for_plot_full <- factor(distances_full_neg_values$classif_for_plot_full, levels = c("-1 - -0.75","-0.75 - -0.5","-0.5 - -0.25","-0.25 - -0.1","-0.1 - 0.1","0.1 - 0.25","0.25 - 0.5","0.5 - 0.75","0.75 - 1"))

S7A <- distances_full_neg_values %>% subset(value>0) %>% ggplot()+
  geom_boxplot(aes(x=factor(variable,levels=c("Upstream-Downstream","Upstream-Upstream",'mRNA',"Downstream-Downstream","FE_distance","Downstream-Upstream","polyA_distance"),labels=c('gene','PITA pre-mRNA','mRNA','PITA pre-mRNA','TSS interval','internal pre−mRNA','TES interval')),y=value,fill=classif_for_plot_full),outlier.color = 'black',outlier.size = 0.5)+
  theme_cowplot()+
  scale_y_log10(labels = label_log())+
  scale_fill_manual(values=met.brewer('Hiroshige')[c(1,2,3,4,5,6,7,8,10)])+
  labs(x='',y='Distance (nt)',fill="Spearman's rho")+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))

#####
no_pi_FE <- subset(df_sup4.1a, new_classif == "no PITA")
pi_FE <- subset(df_sup4.1a, new_classif == "PITA")

cor.test(no_pi_FE$polyA_distance,no_pi_FE$gene_length,method = 'pearson')
cor.test(pi_FE$polyA_distance,pi_FE$gene_length,method = 'pearson')

#Delta AUC mouse
S7B <- ggplot(df_sup4.1c)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  labs(y='ΔAUC R>0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- df_sup4.1c[df_sup4.1c$variable=='soloT',] %>% pull(value)
y <- df_sup4.1c[df_sup4.1c$variable=='dualT',] %>% pull(value)

t.test(x,y)

#Distances mouse
S7C <- df_sup4.1d %>% subset(value>0) %>% ggplot()+
  geom_boxplot(aes(x=variable,y=value,fill=classif_for_plot),outlier.color = 'black',outlier.size = 0.5)+
  theme_cowplot()+
  scale_y_log10(labels = label_log())+
  scale_fill_manual(values=met.brewer('Hiroshige')[c(2,3,5,6,8,10)])+
  labs(x='',y='Distance (nt)',fill="Spearman's rho")+
  theme(
    legend.position = 'none',
    axis.line = element_line(colour = 'black', size = 1),
    axis.text = element_text(size=9), 
    axis.title = element_text(size=10),
    legend.text=element_text(size=9), 
    legend.title=element_text(size=9),
    axis.title.x = element_blank(),
    strip.background=element_rect(fill='gray95'),
    strip.text=element_text(size=10))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

S7D <- df_4c[!df_4c$variable%in%c('Gene length','FE','polyA'),] %>% ggplot()+
  geom_boxplot(aes(x=classif_final,y=delta_normalized), notch=T,outlier.size = 0.5)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  facet_wrap(factor(variable,levels=c('Internal','MaPITA (down-down)','MaPITA (up-up)'),labels = c('internal pre-mRNA','PITA pre-mRNA','PITA pre-mRNA'))~.)+
  labs(y='human length - mouse length')+
  coord_cartesian(ylim = c(-2.5,2.5))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))

#Proportion of PITA types per delta length bin
df_s_4g$bin_for_plot <- factor(df_s_4g$bin_for_plot,levels = c("<-1.5","-1.5 - -0.5","-0.5 - 0.5","0.5 - 1.5",">1.5"))
S7E <- ggplot(df_s_4g[df_s_4g$classif_final!='no PITA' & df_s_4g$classif_final!='Non-PITA both' & !is.na(df_s_4g$bin_for_plot),], aes(x = bin_for_plot, alpha=classif_final)) +
  geom_bar(position = "fill") +
  labs(x = "delta distance (human-mouse)", y = "proportion of genes", fill = "") +
  scale_fill_discrete(name = "classif_final") +
  theme_cowplot()+
  facet_wrap(.~factor(variable,levels=c("FE","Gene length","polyA" ,"Internal","MaPITA (down-down)","MaPITA (up-up)"),labels = c('TSS interval','gene','TES interval','internal pre−mRNA','PITA pre−mRNA','PITA pre−mRNA')),nrow = 1)+
  scale_alpha_manual(values = c(0.5,1))+
  theme(#axis.title.x = element_blank(),
    #legend.position = 'none',
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.line = element_line(colour = 'black', size = 1),
    panel.background = element_rect(color='lightgray'),
    axis.text = element_text(size=9), 
    axis.title = element_text(size=10),
    legend.text=element_text(size=9), 
    legend.title=element_text(size=9),
    strip.background=element_rect(fill='gray95'),
    strip.text=element_text(size=10))

#### Figure S8
S8A <- ggplot(mapping=aes(x=FE_distance,y=gene_length))+
  geom_point(data=subset(df_sup4.1a, new_classif == "no PITA"), color='#A1A1A1', size=0.3,alpha=0.1) +
  geom_point(data=subset(df_sup4.1a, new_classif == "PITA"), color=met.brewer('Hiroshige')[10], size=0.3, alpha=0.1) +
  scale_x_log10(labels = label_log())+
  scale_y_log10(labels = label_log())+
  labs(x='max distance between FE',y='gene length') +
  theme_cowplot()+
  theme(
    axis.line = element_line(colour = 'black', size = 1),
    #panel.background = element_rect(color='lightgray'),
    axis.text = element_text(size=9), 
    axis.title = element_text(size=10),
    legend.text=element_text(size=9), 
    legend.title=element_text(size=9),
    #axis.title.x = element_blank(),
    strip.background=element_rect(fill='gray95'),
    strip.text=element_text(size=10))

S8B <- ggplot(mapping=aes(x=polyA_distance,y=gene_length))+
  geom_point(data=subset(df_sup4.1a, new_classif == "no PITA"), color='#A1A1A1', size=0.3,alpha=0.1) +
  geom_point(data=subset(df_sup4.1a, new_classif == "PITA"), color=met.brewer('Hiroshige')[10], size=0.3, alpha=0.1) +
  scale_x_log10(labels = label_log())+
  scale_y_log10(labels = label_log())+
  labs(x='max distance between polyA',y='gene length') +
  theme_cowplot()+
  theme(legend.position = 'none',
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        #axis.title.x = element_blank(),
        strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))
#Export as a png to avoid AI from crashing
png("FigureS8A_B.png",width=8.5,height=11,units = 'in',res = 300)
ggdraw() + 
  draw_plot(S8A, 0, 0.66, 0.5, 0.33)+
  draw_plot(S8B, 0.5, 0.66, 0.5, 0.33)
dev.off()

S8C <- ggplot(df_s4.1_f,aes(x = pita_type, y = rho,alpha=pita_type)) +
  geom_col() +
  geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci), width = 0.2)+
  theme_cowplot()+
  scale_alpha_manual(values = c(1,0.5,1))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        #strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))+
  labs(y="Spearman's rho deltaTSS-deltaTES (human-mouse)")


sorted_intervals <- df_s7d$max_dist_between_sites[order(as.numeric(gsub("\\[|,.*", "", df_s7d$max_dist_between_sites)))]
sorted_intervals <- sorted_intervals[!duplicated(sorted_intervals)]
df_s7d$max_dist_between_sites <- factor(df_s7d$max_dist_between_sites, levels = sorted_intervals)
df_s7d$deltaAUC <- df_s7d$area_high-df_s7d$area_low

labels_for_x_axis_fig_s7d <- sub(".*,(.*)\\).*", "\\1", df_s7d$max_dist_between_sites) %>% as.numeric() %>% unique() %>% sort()
labels_for_x_axis_fig_s7d <- as.character(labels_for_x_axis_fig_s7d/1000)

S8D <- df_s7d %>% filter(number_of_genes>200) %>% 
  ggplot()+
  geom_col(aes(x=max_dist_between_sites,y=deltaAUC,fill=terminal_3_classif),position='dodge2')+
  theme_minimal(base_size = 15)+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        #axis.line.x = element_line(colour = 'black', linewidth = 1)
  )+
  labs(x='distance between PAS (kb)',y='ΔAUC',fill='')+
  geom_errorbar(aes(x=max_dist_between_sites, ymin = percentile_5,ymax = percentile_95,color=terminal_3_classif),position = 'dodge2',size=0.2)+
  scale_fill_manual(values = met.brewer(name = 'Cassatt2',n = 20)[c(12,14)])+
  scale_color_manual(values=c('gray75','gray75'),guide='none')+
  scale_x_discrete(labels=labels_for_x_axis_fig_s7d)

image_grob <- rasterGrob(readPNG("FigureS8A_B.png"), interpolate = TRUE)


pdf("FigureS8.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(image_grob,0,0.66,1,0.33)+
  draw_plot(S8C, 0.3, 0.33, 0.7, 0.33)+
  draw_plot(S8D, 0.3, 0.33, 0.7, 0.33)
dev.off()

##### Figure S9
#Stackup plots were made by Sergey Venev (github.com/dekkerlab/pita_project)

S9B <- ggplot(all_bins[all_bins$cell_line=='esc',])+
  geom_line(aes(x=bin,y=value,color=pita_classif_psi),size=1.3)+
  theme_cowplot()+
  facet_grid(.~site)+
  scale_x_continuous(breaks = c(0,100,200),labels = c('-7.5kb','TSS','7.5kb'))+
  scale_color_manual(values = c("#d3d3d3","#A2A2AB","#4D87A8","#1E466E"),labels=c(
    expression(paste("no PITA, ",Psi," < 0.5",sep = "")),
    expression(paste("no PITA, ",Psi," > 0.5",sep = "")),
    expression(paste("PITA, ",Psi," < 0.5",sep = "")),
    expression(paste("PITA, ",Psi," > 0.5",sep = ""))
  ))+
  labs(x='distance from TSS',y='obs/exp')+
  theme(legend.position = c(.37,.2),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_blank(),
        strip.background=element_rect(fill='white'),
        strip.text=element_text(size=10))+
  background_grid(major="none",minor="y")

S9C_D <- micro_c_heatmap_insulations_esc %>% ggplot()+
   geom_tile(aes(x=variable,y=ypos,fill=log(value)))+
   scale_fill_gradientn(
     colours = rev(met.brewer('Hiroshige',n = 100)[10:90]),
     na.value = 'white',
     breaks=c(log(1/1.5),0,log(1.5)),
     labels = c(0.67,1,1.5),
     limits = c(log(1/1.5), log(1.5)))+
   theme_cowplot()+
   scale_x_continuous(breaks = c(max(micro_c_heatmap_insulations$variable)/3,2*max(micro_c_heatmap_insulations$variable)/3),labels = c('TSS','TES'))+
   labs(fill='obs/exp')+
   facet_grid(treatment~isoform)+
   theme(axis.title.x = element_blank(),
         strip.background = element_rect(fill="white"),
         axis.text = element_text(size=9),
         axis.title.y = element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y =element_blank(),
         strip.text=element_blank(),
         axis.line = element_blank())
 micro_c_figure_legend <- g.legend(micro_c_figure)
 micro_c_figure <- micro_c_figure+theme(legend.position = 'none')

##### Figure S10
#Stackup plots were made by Sergey Venev (github.com/dekkerlab/pita_project)
 S10B <- ggplot(all_bins[all_bins$cell_line=='hff',])+
   geom_line(aes(x=bin,y=value,color=pita_classif_psi),size=1.3)+
   theme_cowplot()+
   facet_grid(.~site)+
   scale_x_continuous(breaks = c(0,100,200),labels = c('-7.5kb','TSS','7.5kb'))+
   scale_color_manual(values = c("#d3d3d3","#A2A2AB","#4D87A8","#1E466E"),labels=c(
     expression(paste("no PITA, ",Psi," < 0.5",sep = "")),
     expression(paste("no PITA, ",Psi," > 0.5",sep = "")),
     expression(paste("PITA, ",Psi," < 0.5",sep = "")),
     expression(paste("PITA, ",Psi," > 0.5",sep = ""))
   ))+
   labs(x='distance from TSS',y='obs/exp')+
   theme(legend.position = c(.37,.2),
         axis.line = element_line(colour = 'black', size = 1),
         #panel.background = element_rect(color='lightgray'),
         axis.text = element_text(size=8),
         axis.title = element_text(size=9),
         legend.text=element_text(size=9),
         legend.title=element_blank(),
         strip.background=element_rect(fill='white'),
         strip.text=element_text(size=10))+
   background_grid(major="none",minor="y")
 
 S10C_D <- micro_c_heatmap_insulations_hff %>% ggplot()+
  geom_tile(aes(x=variable,y=ypos,fill=log(value)))+
  scale_fill_gradientn(
    colours = rev(met.brewer('Hiroshige',n = 100)[10:90]),
    na.value = 'white',
    breaks=c(log(1/1.5),0,log(1.5)),
    labels = c(0.67,1,1.5),
    limits = c(log(1/1.5), log(1.5)))+
  theme_cowplot()+
  scale_x_continuous(breaks = c(max(micro_c_heatmap_insulations$variable)/3,2*max(micro_c_heatmap_insulations$variable)/3),labels = c('TSS','TES'))+
  labs(fill='obs/exp')+
  facet_grid(treatment~isoform)+
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text = element_text(size=9),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        strip.text=element_blank(),
        axis.line = element_blank()
  )

#### Figure S11
fig_s11a_a <- ggplot(mouse_polII_mut_wt, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = results_cor)) + 
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.5, -.25, 0,.25, 0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
  ) +  
  theme(
    axis.line = element_blank(),
    #panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    #legend.text = element_text(size = 9),
    #legend.title = element_text(size = 9),
    strip.text = element_text(size = 10),
    axis.ticks = element_blank(),
  )+
  scale_y_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  scale_x_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  labs(x='AFE ordinal position',y='ALE ordinal position',fill=expression(atop("Pearson R",paste("AFE",Psi,"~ALE",Psi))))

fig_s11a_legend <- g.legend(fig_s10a_a)
fig_s11a_a <- fig_s10a_a+theme(legend.position = 'none')

fig_s11a_b <- ggplot(mouse_polII_mut_slow, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = results_cor)) + 
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.5, -.25, 0,.25, 0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    name = ""
  ) +  
  theme(
    legend.position = 'none',
    axis.line = element_blank(),
    #panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    strip.text = element_text(size = 10),
    axis.ticks = element_blank(),
  )+
  scale_y_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  scale_x_continuous(expand=c(0, 0),breaks =c(1,2,3))+
  labs(x='AFE ordinal position',y='ALE ordinal position')


#
sp_with_el_index$corr_classif <- factor(sp_with_el_index$corr_classif,levels = c("anti-PITA","no PITA","PITA"),labels =  c("no PITA","no PITA","PITA"))

fig_s11b <- ggplot(sp_with_el_index,aes(x = ave_elongation.velocity_ctrl_DMSO1h, fill=corr_classif,color = corr_classif)) +
  #stat_ecdf() +
  geom_line(stat = 'ecdf',size = 1.3)+
  theme_cowplot() +
  scale_x_log10() +
  coord_cartesian(xlim = c(6,700))+
  labs(y = 'cumulative proportion', x = 'elongation velocity') +
  scale_color_manual(values = c('#A2A2AA', met.brewer('Hiroshige')[10])) +
  theme(
    legend.position = 'none',
    axis.line = element_line(colour = 'black', size = 1),
    panel.background = element_rect(color = 'lightgray'),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    #strip.background = element_rect(fill = 'gray95'),
    strip.text = element_text(size = 10)
  )

fig_s11b_insert <- ggplot(sp_with_el_index,aes(x=corr_classif,y=ave_elongation.velocity_ctrl_DMSO1h,fill=corr_classif))+
  geom_boxplot(notch=T, show.legend = F,outlier.size = 0.2)+
  theme_cowplot()+
  scale_y_log10()+
  labs(y='elongation velocity')+
  scale_fill_manual(values =c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=7),
        axis.title = element_text(size=8),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        #strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))
no_pita <- sp_with_el_index[sp_with_el_index$corr_classif=='no PITA',]
pita <- sp_with_el_index[sp_with_el_index$corr_classif=='PITA',]

ks.test(no_pita$ave_elongation.velocity_ctrl_DMSO1h,pita$ave_elongation.velocity_ctrl_DMSO1h)

#


elong_index_per_length_bin$bin_log2 <- factor(elong_index_per_length_bin$bin_log2, levels = c("[1.14e+04,4.14e+04]","(4.14e+04,7.14e+04]","(7.14e+04,1.01e+05]","(1.01e+05,1.31e+05]","(1.31e+05,1.61e+05]","(1.61e+05,1.91e+05]","(1.91e+05,2.21e+05]","(2.21e+05,2.51e+05]","(2.51e+05,2.81e+05]","(2.81e+05,3.11e+05]","(3.11e+05,3.41e+05]","(3.41e+05,3.71e+05]","(3.71e+05,4.01e+05]","(4.01e+05,4.31e+05]","(4.31e+05,4.61e+05]","(5.51e+05,5.81e+05]"),labels = c("[1.14e+04\n4.14e+04]","(4.14e+04\n7.14e+04]","(7.14e+04\n1.01e+05]","(1.01e+05\n1.31e+05]","(1.31e+05\n1.61e+05]","(1.61e+05\n1.91e+05]","(1.91e+05\n2.21e+05]","(2.21e+05\n2.51e+05]","(2.51e+05\n2.81e+05]","(2.81e+05\n3.11e+05]","(3.11e+05\n3.41e+05]","(3.41e+05\n3.71e+05]","(3.71e+05\n4.01e+05]","(4.01e+05\n4.31e+05]","(4.31e+05\n4.61e+05]","(5.51e+05\n5.81e+05]"))

fig_s11c <- ggplot(elong_index_per_length_bin, aes(x = bin_log2, y = mean_value, fill = corr_classif)) +
  geom_bar(stat = "identity", position = position_dodge2(),alpha=0.8) +
  geom_errorbar(aes(ymin = mean_value - sem, ymax = mean_value + sem), 
                width = 0.25, position = position_dodge(width=0.9)) +
  theme_cowplot() +
  scale_y_log10() +
  labs(x = 'gene length bin', y = 'mean RNAPII velocity') +
  scale_fill_manual(values = c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',
        axis.line = element_line(colour = 'black', size = 1),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        #panel.background = element_rect(color='lightgray'),
        axis.text.y  = element_text(size=9),
        axis.text.x.top = element_text(size=7),
        axis.text.x.bottom = element_text(size=7),
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        #strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))


#
df_s11D$region[which(df_s11D$region == "non+anti-PITA")] <- "noPITA" 
# edited plot
nodata <- ggplot(subset(df_s11D, region=="noPITA")) + 
  stat_smooth(aes(x = bin, y = qinf), method = "loess", se = FALSE) +
  stat_smooth(aes(x = bin, y = qsup), method = "loess", se = FALSE)
pitadata <- ggplot(subset(df_s10D, region=="PITA")) + 
  stat_smooth(aes(x = bin, y = qinf), method = "loess", se = FALSE) +
  stat_smooth(aes(x = bin, y = qsup), method = "loess", se = FALSE)
# build plot object for rendering 
gg_no <- ggplot_build(nodata)
gg_pita <- ggplot_build(pitadata)
# extract data for the loess lines from the 'data' slot
ribbon_no <- data.frame(x = gg_no$data[[1]]$x,
                        ymin = gg_no$data[[1]]$y,
                        ymax = gg_no$data[[2]]$y)
ribbon_no$region = "noPITA"
ribbon_pita <- data.frame(x = gg_pita$data[[1]]$x,
                          ymin = gg_pita$data[[1]]$y,
                          ymax = gg_pita$data[[2]]$y) 
ribbon_pita$region = "PITA"
ribbon_data <- rbind(ribbon_no, ribbon_pita)

fig_s11D <- ggplot(df_s11D)+
  # add ribbon
  #geom_ribbon(aes_string(fill ="region"), alpha=0.3) +
  geom_ribbon(data = ribbon_data, aes(x=x, ymin=ymin, ymax=ymax, fill=factor(region)), alpha=0.1) +
  # add values
  geom_line(aes(x=bin, y=value, color=factor(region)),stat = "smooth",size=1.5,alpha=0.8) +
  # attributes
  #ylim(0,2.1)+
  scale_fill_manual(values=c("#A2A2AA","#1E466E"),guide=F) +
  #scale_color_manual(values=c("darkgrey","1E466E"),labels=c("no PITA","PITA")) +
  scale_color_manual(values =c('#A2A2AA',met.brewer('Hiroshige')[10]),labels=c("no PITA","PITA"))+
  labs(y="elongation velocity",x="Gene body")+
  scale_x_continuous(breaks = c(0,100), labels =c("TSS","TES")) +
  ### dna lines
  annotate("segment", x = -10, xend = 110, y = 0, yend = 0,color= "grey", linetype="solid")  +
  ### gene 
  annotate("segment", x = 1, xend = 1, y = 0.05, yend = 0.1,color= "grey25", linetype="solid")  +
  annotate("segment", x = 1, xend = 5, y = 0.1, yend = 0.1,color= "grey25", linetype="solid", arrow = arrow(type = "closed", length = unit(0.01, "npc")))  +
  annotate("text", label=expression(TSS[1]), x = 1, y = 0.14, color="grey25") +
  annotate("rect", xmin = 1, xmax = 100, ymin = -0.05, ymax = 0.05,color= "grey50",fill="grey50")  +
  annotate("point", x = 100, y = 0.075, colour = "grey25", size = 5, shape="\u25BC") +
  annotate("text", label=expression(PAS[n]), x = 100, y = 0.14, color="grey25") +
  theme_cowplot()+
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank())


fig_s11D_legend <- g.legend(fig_s11D)
fig_s11D <- fig_s11D+theme(legend.position = 'none')

#
df_s11E <- df_s11E[!duplicated(df_s11E),]
fig_s11E <- df_s11E %>%  filter(pita_classif!='anti PITA') %>%
  ggplot()+
  geom_smooth(aes(x=proportional_position,y=logodds,color=pita_classif),fill='gray90')+
  scale_x_continuous(breaks = c(0,0.5,1),labels = c('upstream','intermediate','downstream'))+
  theme_cowplot()+
  scale_color_manual(values =c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(axis.line.x=element_blank(),
        legend.position = 'none')+
  labs(y='predicted PAS strength (logit)',x='PAS position')

pdf("FigureS11.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig_s11a_a, 0, 0.66, 0.4, 0.33)+
  draw_plot(fig_s11a_b, 0.6, 0.66, 0.4, 0.33)+
  draw_plot(fig_s11a_legend,0.47,0.8,0.1,0.1)+
  draw_plot(fig_s11b,0,0.33,0.4,0.33)+
  draw_plot(fig_s11b_insert,0.19,0.37,0.2,0.17)+
  draw_plot(fig_s11c,0.42,0.33,0.58,0.33)+
  draw_plot(fig_s11D,0,0,0.5,0.33)+
  draw_plot(fig_s11E,0.5,0,0.5,0.33)+
  draw_plot_label(label=c("A","B","C","D","E"), x=c(0, 0, 0.42, 0,0.5), y=c(1,0.66,0.66,0.33,0.33),family="Helvetica", color="black")

dev.off()


#### Figure S12
## LRS replicate correlations
figS12A <- ggplot(combo_repcorr_gene, aes(x=rep1_mean, y=rep2_mean)) + 
  geom_point(size=1, alpha=0.15) + 
  geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
  scale_x_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  labs(x="mean genomic distance, rep 1",
       y="mean genomic distance, rep 2") +
  annotate("text", x = 1e2+1500, y=1e5, label="Pearson R = 0.76", size=2.5, color="red") +
  annotate("text", x = 1e2+1500, y=1e5-30000, label="3,373 AFEs across 1,462 genes", size=2.5, color="red") +
  theme_cowplot() +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8))

# Pearson corr = 0.76
cor(combo_repcorr_gene$rep1_mean, combo_repcorr_gene$rep2_mean, use = "complete")
# number of AFEs = 3373
nrow(unique(subset(combo_repcorr_gene, !is.na(rep1_mean) & !is.na(rep2_mean))))
# number of genes = 1462
length(unique(subset(combo_repcorr_gene, !is.na(rep1_mean) & !is.na(rep2_mean))$gene))
## violin plots across all reads, each AFE separate
figS12B <- ggplot(subset(bothreps_10m_stranded_subsampled, mapped_length >= 100 & afe_reorder <= 6), aes(x=factor(afe_reorder), y = mapped_length)) + 
  geom_violin() + 
  geom_boxplot(width=0.15, notch=T, outliers = F) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  labs(x="AFE ordinal position",y="genomic distance (nt), direct RNA-seq") +
  theme_cowplot() +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8))
## global LRS results using cDNA
figS12C <- ggplot(subset(cDNA_unstranded_10m_subbed, genomic_length >= 100 & afe_reorder <= 5), aes(x=factor(afe_reorder), y = genomic_length)) + 
  geom_violin() + geom_boxplot(width=0.15, notch=T, outliers = F) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  labs(x="AFE ordinal position",y="genomic distance (nt), cDNA-seq") +
  theme_cowplot() +
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8))


figS12D <- ggplot(CpG_plot_data, aes(x = position, y = fraction, fill = prediction)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group) +
  scale_fill_manual(
    values = c("No" = "black", "Yes" = "grey70"),
    labels = c("No" = "without CpG island", "Yes" = "with CpG island"),
    breaks = c("Yes", "No"),  # This sets legend order (gray, then black),
    name = NULL  # removes legend title
  ) +
  labs(
    x = "promoter ordinal position",
    y = "fraction of promoters"
  ) +
  theme_cowplot() +
  theme(axis.text = element_text(size=7),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 8)
  )


figS12E <- ggplot(promoter_prediction_plot_data, aes(x = position, y = fraction, fill = prediction)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group) +
  scale_fill_manual(
    values = c("No" = "black", "Yes" = "grey70"),
    breaks = c("Yes", "No"),  # <- legend order,
    labels = c("No" = "no predicted promoter", "Yes" = "predicted promoter"),
    name = NULL  # removes legend title
    
  ) +
  labs(
    x = "promoter ordinal position",
    y = "fraction of promoters"
  ) +
  theme_cowplot() +
  theme(axis.text = element_text(size=7),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 8)
  )


## combine
pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/FigureS12.pdf",width=8,height=11.5, useDingbats = F)
ggdraw() + # 
  draw_plot(figS12A, 0, 0.75, 0.33, 0.25) +
  draw_plot(figS12B, 0.33, 0.75, 0.33, 0.25) +
  draw_plot(figS12C, 0.66, 0.75, 0.33, 0.25) +
  draw_plot(figS12D, 0, 0.5, 0.33, 0.25) +
  draw_plot(figS12E, 0.33, 0.5, 0.33, 0.25) +
  
  cowplot::draw_plot_label(label=c("A","B","C",'D','E'), x=c(0, 0.33, 0.66,0,0.33), y=c(1, 1, 1,0.75,0.75), 
                           family="Helvetica", color="black")
dev.off()



























