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
library(magick)
library(ggdendro)
library(dendextend)
library(ape)
library(VennDiagram)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(ggdist)
library(gridExtra)
library(ComplexUpset)
library(ggthemes)
library(scales)
library(ragg)
library(grid)
library(gridGraphics)
library(MetBrewer)
library(ggdraw)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
####FIGURE1A####

pdf_path <- "/projectnb/encore/carroll8/carroll8/2023_HITindexoutputsforPITApaper/Final_Code/Figure1A_850by100.pdf"
pdf_image <- image_read(pdf_path)
plot(pdf_image)
png_path <- "/projectnb/encore/carroll8/carroll8/2023_HITindexoutputsforPITApaper/Final_Code/29NOV2023_Figure1A_850by100.png"  # Specify the path to save the PNG image
image_write(pdf_image, path = png_path)
Figure1A <- rasterGrob(as.raster(pdf_image))

### FIGURE 1B #####
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA")
load("FIGURE1B_Data.RData")
FIGURE1Bdata$type <- c("first", "first", "first", "first", "last", "last", "last", "last")
FIGURE1Bdata$group <- c(1,2,3,4,1,2,3,4)
sci_labels <- function(x) {
  parse(text = sprintf("10^{%d}", log10(x)))
}

mycols  <- c("gray22", "gray76")
breaks = 10**(1:5 * 1)


Figure1B_final <- ggplot(FIGURE1Bdata, aes(x = group, y = mean, fill = type)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), 
                position = position_dodge(width = 0.9), 
                width = 0.2) +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(
    trans = 'log10',
    breaks = c(1, 10, 100, 1000, 10000),
    labels = sci_labels
  ) +
  theme_cowplot() +
  theme(
    axis.line = element_line(colour = 'black', linewidth = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 11),
    axis.text.y = element_text(angle = 0, vjust = 0, hjust = .5, size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 20),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white", colour = "white"),
    legend.position = c(0.99, 0.93),  # Adjust the legend position (top-right)
    legend.justification = c(1, 1),  # Adjust the justification for top-right corner
    legend.title.align = 0.5
  ) +
  labs(
    x = "exons per gene",
    y = "mean exon count",
    fill = "exon type"
  ) +
  guides(fill = guide_legend(title = "exon type"))


Figure1B <- Figure1B_final + theme(plot.margin = unit(c(0,0,0,0), "cm"))

### FIGURE1C ###

setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
load("FIGURE1C_Data.RData")
eachgene <- FIGURE1Cdata[[1]]
numberofobservations <- FIGURE1Cdata[[2]]
pearsoncor <- FIGURE1Cdata[[3]]

Figure1C <- ggplot(eachgene, aes(Freq.x, Freq.y)) +
  geom_point(data = numberofobservations, aes(color = n), size = 4) +
  scale_color_gradientn(name = "genecount",  # Stack title
                        colors = c("#aadce0" , "#679AA5", "#1c5062","#1e466e"),
                        # colors = c("#aadce0", "#72bcd5", "#376795", "#1e466e"),
                        trans = "log",  # Apply log scale
                        breaks = c(1, 10, 100),  # Specify breaks for log scale
                        labels = scales::number_format(accuracy = 1)
  ) +
  geom_smooth(data = eachgene, aes(Freq.x, Freq.y), color = 'black', method = "lm", se = FALSE) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  theme_cowplot() +
  # ... (rest of your theme settings)
  labs(
    x = "number of first exons",
    y = "number of last exons",
    fill = NULL  # Remove the fill legend title
  ) +
  guides(fill = guide_legend(title = NULL)) +  # Remove the fill legend title
  theme(
    legend.title.align = -0.5,  # Center the legend title
    text = element_text(size = 12),  # Set default text size
    axis.text = element_text(size = 11),  # Set axis label text size
    axis.title = element_text(size = 12),  # Set axis title text size
    legend.text = element_text(size = 11),  # Set legend text size
    legend.title = element_text(size = 12),  # Set legend title text size
    legend.position = "top",  # Place legend above the plot
    legend.justification = "center",  # Center the legend within the plot area
    legend.box = "horizontal"  # Arrange legend items horizontally
  )

Figure1C

### FIGURE1D ###


setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
load("FIGURE1D_Data.RData")

FIGURE1D_matrix <- as.matrix(FIGURE1Ddata)
col_fun = colorRamp2(c(.4, .5, .6, .7), c("#aadce0", "#72bcd5", "#528fad", "#1e466e"))

Figure1D_heatmap <- Heatmap(FIGURE1D_matrix, name = "Pearson R", col = col_fun, show_row_names = TRUE, width = unit(3, "cm"),
                            row_dend_width = unit(1.5, "cm"), column_labels = "", row_names_gp = gpar(fontsize = 11),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = list(direction = "horizontal", legend_position = "right"))

Figure1D <- draw(Figure1D_heatmap, heatmap_legend_side = "top")
Figure1D <- grid.grabExpr(draw(Figure1D))

### FIGURE 1EA ### 
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
load("FIGURE1E_Data.RData")

AFEorder <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
ALEorder <- c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
AFE_ALE_pearson <- data.frame(AFEorder, ALEorder, FIGURE1Edata)

Figure1EA <- ggplot(AFE_ALE_pearson, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = FIGURE1Edata)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.5, -.25, 0,.25,  0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    name = ""
  ) +  theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10), size = 12),
    axis.title.y = element_text(margin = margin(r = 10, b = -20), size = 12),
    axis.text.x = element_text(margin = margin(t = -10), size = 12),
    axis.text.y = element_text(margin = margin(r = -10), size = 12)
  ) +
  labs(
    x = "AFE Order",
    y = "ALE Order"
  )

Figure1EA


### FIGURE1EB ###

AFEorder <- c(1,1,1,2,2,2,3,3,3)
ALEorder <- c(1,2,3,1,2,3,1,2,3)
AFE_ALE_pearson <- data.frame(AFEorder, ALEorder, FIGURE1EBdata)

Figure1EB <- ggplot(AFE_ALE_pearson, aes(x = AFEorder, y = ALEorder)) +
  geom_tile(aes(fill = FIGURE1EBdata)) +
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
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
    x = "AFE Order",
    y = "ALE Order"
  )

Figure1EB



Figure1EB_withlegend <- ggplot(AFE_ALE_pearson, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = FIGURE1EBdata)) + 
  scale_fill_gradientn(
    colours = c("#e76254", "#ef8a47", "#f7aa58", "white", "#72bcd5", "#528fad", "#1e466e"),
    breaks = c(-0.5, -.25, 0,.25,  0.5),
    labels = c(-0.5, -.25, 0,.25, 0.5),
    limits = c(-0.5, 0.5),
    name = ""
  ) +  
  theme(plot.background = element_rect(fill = "white", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "white"),  axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10), size = 12),
        axis.title.y = element_text(margin = margin(r = 10, b = -20), size = 12),
        axis.text.x = element_text(margin = margin(t = -10), size = 12),
        axis.text.y = element_text(margin = margin(r = -10), size = 12),
        legend.title = element_text(hjust = 0.5)) +  labs(x = "AFE Order", y = "ALE Order")

legend_1E <- get_legend(Figure1EB_withlegend)


figure1 <-  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  draw_plot(Figure1A, x = 0.01, y = .75, width = .95, height = .15) +
  draw_plot(Figure1B, x = 0.02, y = .4, width = .29, height = .31) +
  draw_plot(Figure1C, x = 0.35, y = .4, width = .29, height = 0.31) +
  draw_plot(Figure1D, x = .68, y = .425, width = .27, height = 0.30) +
  draw_plot(Figure1EA, x = 0.02, y = 0.01, width = 0.44, height = 0.35) +
  draw_plot(Figure1EB, x = .49, y = 0.015, width = 0.415, height = 0.35) +
  draw_plot(legend_1E, x = .92, y = 0.02, width = 0.08, height = 0.4)

setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
ggsave("07DEC2023_figure1_updates.pdf", plot = figure1, width = 8.5, height = 10)



#####extended figure 1A ####
#the numbers are directly inputted from the "means" output  
x <- draw.pairwise.venn(1447, 851, 259, 
                        category = c("AFEs", "ALEs"), 
                        fill = c("gray22", "gray76"),
                        cat.pos = c(0, 0),
                        fontfamily = rep("Helvetica", 3),
                        cat.fontface = rep("plain", 2), 
                        cat.fontfamily = rep("Helvetica", 2))

EXTFIG1A <- x

###extended figure 1B#####


# Rename columns with unique names in order to allow for proper heatmap generation
colnames(master_df) <- c("gene", paste0("", seq_along(master_df)[-1]))
master_dataframe <- master_df
rownames(master_dataframe) <- master_dataframe[,1]
master_dataframe_NAs <- master_dataframe
master_dataframe_NAs <- as.matrix(master_dataframe_NAs[,-1])
master_dataframe <- as.matrix(master_dataframe[,-1])
master_dataframe[is.na(master_dataframe)] <- 0 #in order to do hierarchical clustering, set NA to 0 (NA expression means no expression) 
min(master_dataframe_NAs, na.rm = TRUE) #determine minimum PSI value in order to define heatmap color range  
ggplot_data <- melt(master_dataframe)
#extract legend only 
col_fun_NAs <- colorRamp2(c(.1, .55,  1), c("#aadce0","#72bcd5","#1e466e"))  # Define the color gradient
extract_legend <- ggplot(ggplot_data, aes(x = Var2, y = value)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn( colours = c("#aadce0","#72bcd5","#1e466e"), breaks = c(.1, .55, 1), labels = c(.1, .55, 1), limits = c(.1, 1), guide = guide_colorbar(title = NULL))
maxPSI_legend <- get_legend(extract_legend)
#extract heatmap only 
col_fun <- colorRamp2(c(0, .09, .1, .55,  1), c("white" , "white", "#aadce0","#72bcd5", "#1e466e"))  # Define the color gradient
maxPSI_heatmap <- Heatmap(master_dataframe, show_row_names = FALSE, col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE, km = 3, heatmap_legend_param = list(title = "maximum PSI", at = c(0.1, .55, 1), labels = c("0.1", ".55", "1"), limits = c(0.1, 1)), show_heatmap_legend = FALSE, row_title = "genes")
maxPSI_heatmap <- draw(maxPSI_heatmap)
maxPSI_heatmap <- grid.grabExpr(draw(maxPSI_heatmap))

extended_figure1<-  ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  draw_plot(EXTFIG1A, x = 0.01, y = .65, width = .32, height = .25) +
  draw_plot(maxPSI_heatmap, x = 0.34, y = .55, width = .57, height = .41) +
  draw_plot(maxPSI_legend, x = 0.92, y = .69, width = .07, height = .12)

#Load all the data for long read seq
load("~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/Paper/figures/data_for_figures.RData")

cage_color <- met.brewer('Archambault',n = 40)[5]
polyA_color <- met.brewer('Archambault',n = 20)[10]

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

g.legend <- function(a.plot){
  tmp <- ggplot_gtable(ggplot_build(a.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##### Figure 2

dotted_LRS <- myo10_example[myo10_example$dataset=='reads_dotted',]
LRS_loop <- myo10_example[myo10_example$dataset=='reads',]
polyA_loop <- myo10_example[myo10_example$dataset=='polyA',]
first_exons_loop <- myo10_example[myo10_example$dataset=='FirstExons',]
isoforms_loop <- myo10_example[myo10_example$dataset=='isoforms',]
dotted_isoforms_loop <- myo10_example[myo10_example$dataset=='isoforms_dotted',]

id_add_isoforms <- min(-1*LRS_loop$id)*0.01+min(-1*LRS_loop$id)

#These reads are from sample ENCFF142LPL
fig_2_a <- ggplot()+
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
fig_2_e_df$exon_type <- factor(fig_2_e_df$exon_type,levels = c('Unique FE or polyA','Alternative FE & polyA'))

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

fig_2_f_df$variable <- factor(fig_2_f_df$variable,levels = c('Unique FE or polyA','Alternative FE & polyA'))
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


encode_auc$variable <- factor(encode_auc$variable,levels = c('Unique FE or polyA','Alternative FE & polyA'),labels = c('single FE\nor polyA','multiple FE\n& polyA'))

fig_2_g <- ggplot(encode_auc)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC >0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- encode_auc[encode_auc$variable=='Unique FE or polyA',] %>% pull(value)
y <- encode_auc[encode_auc$variable=='Alternative FE & polyA',] %>% pull(value)
t.test(x,y)

#pdf("~/Desktop/fig2_Rout.pdf",width=8.5,height=11)
pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure2.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig_2_a, 0, 0.66, 1, 0.33)+
  draw_plot(fig_2_b, 0, 0.4, 0.33, 0.25)+
  draw_plot(fig_2_e, 0, 0, 0.33, 0.37)+
  draw_plot(fig_2_f, 0.33, 0, 0.33, 0.37)+
  draw_plot(fig_2_g, 0.66, 0, 0.33, 0.37)
dev.off()

##### Figure 3

df_3a$classif_corr <- factor(df_3a$classif_corr,levels = c('MiPITA','Non-PITA','MaPITA'),labels = c('anti-PITA','no PITA','PITA'))

fig3_a <- df_3a %>% subset(domain_change=='Different domains') %>% ggplot()+
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


sample_phast$region <- factor(sample_phast$region, levels=c("upstart","downstart","upend","downend"))
sample_phast$position <- factor(sample_phast$position, levels=c("upstream","downstream"))
sample_phast$terminal <- factor(sample_phast$terminal, levels=c("TSS","TES"))
sample_phast$type <- factor(sample_phast$type, labels = c("no PITA","PITA"))

fig3_b <-ggplot(sample_phast, aes(y=mean, x=factor(position), fill=factor(type))) + 
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

df_3d$variable <- factor(df_3d$variable,levels = c('Single FE or polyA','Multiple FE & polyA'),labels = c('Single FE\nor polyA','Multiple FE\n& polyA'))
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


fig3_d <-ggplot(ENSG00000007923_reads_acc_tissues)+
  geom_segment(aes(x=start_FE,xend=end_LE,y=id,yend=id,color=factor(Biosample.term.name,labels =  c('astrocyte','iPSC','lung'))),size=0.1)+
  theme_cowplot()+
  scale_color_manual(values = met.brewer(name = 'Cross')[c(1,3,6)])+
  labs(x='genomic coordinate on chr1 (Mb)',color='tissue',title = 'DNAJC11')+
  scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
  theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_text(size = 13, face = 'bold'),axis.line = element_blank())+
  geom_segment(dotted_isoforms_ENSG00000007923,mapping = aes(y=id,x=start,xend=end,yend=id),linetype='solid',size=0.05,color='black')+
  geom_segment(isoforms_loop_ENSG00000007923,mapping = aes(x=start,xend=end,y=id,yend=id),size=3,color='#EA993E')

pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure3_inprogress.pdf",width=8.5,height=11)
#pdf("~/Desktop/Figure3_inprogress.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig3_a, 0, 0.66, 0.33, 0.33)+
  draw_plot(fig3_b, 0.33, 0.66, 0.33, 0.33)+
  draw_plot(fig3_c, 0.66, 0.66, 0.33, 0.33)+
  draw_plot(fig3_d, 0, 0.33, 1, 0.33)
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

fig4_c <- ggplot(mapping=aes(x=FE_distance,y=polyA_distance))+
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


no_pi_FE <- subset(df_4b, classif_for_plot_new == "Non-MaPITA")
pi_FE <- subset(df_4b, classif_for_plot_new == "MaPITA")

cor.test(no_pi_FE$FE_distance,no_pi_FE$polyA_distance,method = 'pearson')
cor.test(pi_FE$polyA_distance,pi_FE$gene_length,method = 'pearson')



#Export 4c to make it easier to open it in AI
png("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure4c.png",width=8.5,height=11,units = 'in',res = 300)
ggdraw() + 
  draw_plot(fig4_c, 0, 0.28, 0.33, 0.28)
dev.off()


df_4c$classif_final <- factor(df_4c$classif_final,levels = c('Non-PITA both','MaPITA mouse','MaPITA human','MaPITA both'),labels = c('no PITA','PITA\nmouse','PITA\nhuman','PITA\nmouse&human'))
df_4c <- df_4c[df_4c$classif_final!='PITA\nmouse&human',]

fig4_d <- df_4c[df_4c$variable%in%c('Gene length','FE','polyA'),] %>% ggplot()+
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
image_grob <- rasterGrob(readPNG("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure4c.png"), interpolate = TRUE)

pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure4_inprogress.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(fig4_a, 0, 0.66, 0.33, 0.33)+
  draw_plot(fig4_b, 0.33, 0.66, 0.66, 0.33)+
  draw_plot(fig_4ab_legend,0.5,0.55,0.3,0.1)+
  #draw_plot(image_grob, 0, 0.28, 0.33, 0.28)+
  draw_plot(fig4_d, 0.33, 0.25, 0.66, 0.28)
dev.off()
#####

micro_c_heatmap_insulations$treatment <- factor(micro_c_heatmap_insulations$treatment,levels = c('ctrl','pita'),labels = c('no PITA','PITA'))
micro_c_heatmap_insulations$isoform <- factor(micro_c_heatmap_insulations$isoform,levels = c('full','up','down'),labels = c('upstreamTSS-downstreamTES','upstreamTSS-upstreamTES','downstreamTSS-downstreamTES'))
micro_c_heatmap_insulations_esc <- micro_c_heatmap_insulations[micro_c_heatmap_insulations$cell_line=='esc',]
micro_c_heatmap_insulations_hff <- micro_c_heatmap_insulations[micro_c_heatmap_insulations$cell_line=='hff',]

micro_c_figure <- micro_c_heatmap_insulations_esc %>% ggplot()+
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

micro_c_figure_legend <- g.legend(micro_c_figure)
micro_c_figure <- micro_c_figure+theme(legend.position = 'none')


all_bins$site <- factor(all_bins$site,levels = c("upstreamTSS","downstreamTSS"))

f5_b_final <- ggplot(all_bins)+
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
        strip.text=element_text(size=10),
  )+
  background_grid(major="none",minor="y")


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

fig5_legend <- g.legend(fig5a_right)
fig5a_right <- fig5a_right+theme(legend.position = 'none')

sp_with_el_index$corr_classif <- factor(sp_with_el_index$corr_classif,levels = c("anti-PITA","no PITA","PITA"))

fig5b <- ggplot(sp_with_el_index,aes(x = ave_elongation.velocity_ctrl_DMSO1h, fill=corr_classif,color = corr_classif)) +
  #stat_ecdf() +
  geom_line(stat = 'ecdf',size = 1.3)+
  theme_cowplot() +
  scale_x_log10() +
  coord_cartesian(xlim = c(6,700))+
  labs(y = 'cumulative proportion', x = 'elongation velocity') +
  scale_color_manual(values = c(met.brewer('Hiroshige')[2], '#A2A2AA', met.brewer('Hiroshige')[10])) +
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

fig5b_insert <- ggplot(sp_with_el_index,aes(x=corr_classif,y=ave_elongation.velocity_ctrl_DMSO1h,fill=corr_classif))+
  geom_boxplot(notch=T, show.legend = F,outlier.size = 0.2)+
  theme_cowplot()+
  scale_y_log10()+
  labs(y='elongation velocity')+
  scale_fill_manual(values =c(met.brewer('Hiroshige')[2],'#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=7), 
        axis.title = element_text(size=8),
        legend.text=element_text(size=7), 
        legend.title=element_text(size=7),
        #strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))


anti_pita <- sp_with_el_index[sp_with_el_index$corr_classif=='anti-PITA',]
no_pita <- sp_with_el_index[sp_with_el_index$corr_classif=='no PITA',]
pita <- sp_with_el_index[sp_with_el_index$corr_classif=='PITA',]

ks.test(anti_pita$ave_elongation.velocity_ctrl_DMSO1h,no_pita$ave_elongation.velocity_ctrl_DMSO1h)
ks.test(anti_pita$ave_elongation.velocity_ctrl_DMSO1h,pita$ave_elongation.velocity_ctrl_DMSO1h)
ks.test(no_pita$ave_elongation.velocity_ctrl_DMSO1h,pita$ave_elongation.velocity_ctrl_DMSO1h)

#Fold change calculation
mean(pita$ave_elongation.velocity_ctrl_DMSO1h)/mean(no_pita$ave_elongation.velocity_ctrl_DMSO1h)


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figure5_inprogress_new.pdf",width=8.5,height=11)
#pdf("~/Desktop/Figure5_inprogress_new.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(micro_c_figure,0, 0.8, 1, 0.20)+
  draw_plot(micro_c_figure_legend,1, 0.66, 0.08, 0.23)+
  
  draw_plot(f5_b_final, 0, 0.33, 0.6, 0.33)+
  
  
  draw_plot(fig5a_left, 0, 0, 0.3, 0.33)+
  draw_plot(fig5a_mid, 0.3, 0, 0.3, 0.33)+
  draw_plot(fig5a_right, 0.6, 0, 0.3, 0.33)+
  draw_plot(fig5_legend, 0.91, 0.10, 0.08, 0.20)+
  draw_plot(fig5b, 0.585, 0.33, 0.41, 0.33)+
  draw_plot(fig5b_insert, 0.79, 0.38, 0.20, 0.15)
dev.off()




################
#Supp figures

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

pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/inter_intra_cor.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(schematic_intra_mol_fig, 0, 0.33, 0.5, 0.33)+
  draw_plot(schematic_inter_mol_fig, 0.5, 0.33, 0.5, 0.33)
dev.off()

#2.1
figs2_1a <- ggplot(df_sup_2.1a)+
  geom_boxplot(aes(x=as.factor(sampling_factor),y=auc_sampling,fill=factor(variable,levels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8)+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(x='proportion of sampled library',y='ΔAUC R>0',fill='Exon type',title = '')+
  scale_x_discrete(breaks=seq(0,1,0.1))+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))


figs2_1b <- ggplot(df_sup_2.1b)+
  geom_boxplot(aes(x=factor(variable,levels = c('single FE or polyA','multiple FE & polyA')),y=value,fill=factor(variable,levels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC >0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- df_sup_2.1b[df_sup_2.1b$variable=='single FE or polyA',] %>% pull(value)
y <- df_sup_2.1b[df_sup_2.1b$variable=='multiple FE & polyA',] %>% pull(value)
t.test(x,y)

figs2_1c <- ggplot(df_sup_2.1c)+
  geom_boxplot(aes(x=factor(bin_deciles,levels=c('[93,1.08e+03]','(1.08e+03,1.4e+03]','(1.4e+03,1.66e+03]','(1.66e+03,1.89e+03]','(1.89e+03,2.12e+03]','(2.12e+03,2.38e+03]','(2.38e+03,2.69e+03]','(2.69e+03,3.09e+03]','(3.09e+03,3.73e+03]','(3.73e+03,1.38e+04]'),labels = c('[93\n1.08e+03]','(1.08e+03\n1.4e+03]','(1.4e+03\n1.66e+03]','(1.66e+03\n1.89e+03]','(1.89e+03\n2.12e+03]','(2.12e+03\n2.38e+03]','(2.38e+03\n2.69e+03]','(2.69e+03\n3.09e+03]','(3.09e+03\n3.73e+03]','(3.73e+03\n1.38e+04]')),y=corr))+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y="Spearman's rho",x='read length')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))

figs2_1d <- ggplot(df_sup_2.1d,aes(x=factor(variable,levels = c("Unique FE or LE","Alternative FE & LE"),labels = c('single FE or polyA','multiple FE & polyA')),y=value,fill=factor(variable,levels = c("Unique FE or LE","Alternative FE & LE"),labels = c('single FE or polyA','multiple FE & polyA'))))+
  stat_summary(geom = "bar", fun = mean, position = position_dodge(width = 0.9), alpha = 0.8)+
  theme_cowplot()+
  scale_fill_manual(values=c('#808080',met.brewer('Hokusai3',n = 10)[10]))+
  stat_summary(geom = "errorbar",fun.data = mean_se,position = position_dodge(width = 0.9),color = '#75726f',size = 0.2,width = 0.5)+
  ylim(-0.2,0.2)+
  geom_vline(xintercept = 0,linetype='dashed')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  labs(y='ΔAUC R>0')


figs2_1e <- ggplot(df_sup_2.1f)+
  geom_density(aes(x=corr,fill=factor(exon_type,levels = c("Unique FE or polyA","Alternative FE & polyA"),labels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8)+
  scale_fill_manual(values=c('#808080',met.brewer('Hokusai3',n = 10)[10]))+
    theme_cowplot()+
    theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13))+
    labs(x="Pearson's R")+
    xlim(-0.5,0.5)+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())

figs2_1f <- ggplot(df_sup_2.1e)+
  geom_boxplot(aes(x=factor(variable,levels = c('single FE or polyA','multiple FE & polyA')),y=value,fill=factor(variable,levels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_minimal(base_size = 25,base_line_size = 0.8)+
  labs(y='ΔAUC >0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

figs2_1g <- ggplot(df_sup_2.1g)+
  geom_boxplot(aes(x=as.factor(threshold),y=value,fill=factor(variable,levels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8)+
  theme_cowplot()+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(y='ΔAUC >0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

figs2_1h <- ggplot(df_sup_2.1h)+
  geom_boxplot(aes(x=same_different,y=value,fill=same_different))+
  theme_cowplot()+
  #scale_fill_manual(values=c("#d4a0a4","#4b7b7e"))
  #scale_fill_manual(values=c("#d4a0a4","#b18e5e"))
  scale_fill_manual(values=c("#d4a0a4","#fbdcb1"))+
  #scale_fill_manual(values=c("#d4a0a4","#4b7b7e"))
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  labs(y='correlation between samples')

x <- df_sup_2.1h[df_sup_2.1h$same_different=='non-replicates',] %>% pull(value)
y <- df_sup_2.1h[df_sup_2.1h$same_different=='replicates',] %>% pull(value)
wilcox.test(x,y)


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/fig_sup2_1.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs2_1a, 0, 0.66, 0.66, 0.33)+
  draw_plot(figs2_1b, 0.66, 0.66, 0.33, 0.33)+
  draw_plot(figs2_1c, 0, 0.33, 0.6, 0.33)+
  draw_plot(figs2_1d, 0.6, 0.33, 0.2, 0.33)+
  draw_plot(figs2_1e, 0.8, 0.33, 0.2, 0.33)+
  draw_plot(figs2_1g, 0, 0, 0.6, 0.33)+
  draw_plot(figs2_1f, 0.6, 0, 0.2, 0.33)+
  draw_plot(figs2_1h, 0.8, 0, 0.2, 0.33)
dev.off()

#2.2
figs2_2a <- ggplot(df_sup2.2a)+
  geom_boxplot(aes(x=tech,y=value,fill=variable),alpha=0.8,outlier.size = 0.1)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  labs(x='',y='ΔAUC >0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  ylim(-0.25,0.5)


x <- df_sup2.2a[df_sup2.2a$tech=='TIF-Seq2' & df_sup2.2a$variable=='single FE or polyA',] %>% pull(value)
y <- df_sup2.2a[df_sup2.2a$tech=='TIF-Seq2' & df_sup2.2a$variable=='multiple FE & polyA',] %>% pull(value)
t.test(x,y) 

x <- df_sup2.2a[df_sup2.2a$tech=='FLAM-seq' & df_sup2.2a$variable=='single FE or polyA',] %>% pull(value)
y <- df_sup2.2a[df_sup2.2a$tech=='FLAM-seq' & df_sup2.2a$variable=='multiple FE & polyA',] %>% pull(value)
t.test(x,y) 

figs2_2b <- ggplot(df_sup2.2b)+ #Our data through their pipeline but ours for exon classification
  geom_bar(aes(x=predominance_new,fill=factor(exon_type,levels = c('Unique FE or polyA','Alternative FE & polyA'),labels = c('single FE or polyA','multiple FE & polyA'))),alpha=0.8)+
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


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/fig_sup2_2.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs2_2a, 0, 0.66, 0.5, 0.33)+
  draw_plot(figs2_2b, 0.5, 0.66, 0.5, 0.33)+
  draw_plot(figs2_2c, 0, 0.33, 0.5, 0.33)+
  draw_plot(figs2_2d, 0.5, 0.33, 0.5, 0.33)
dev.off()

#2.3

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

pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/fig_sup2_3.pdf",width=8.5,height=11)
print(ggplot(LRS_loop,mapping=aes(x=start,xend=end,y=(-1*id),yend=(-1*id)))+
        geom_segment(dotted_LRS,mapping=aes(y=(-1*id),x=start,xend=end,yend=(-1*id)),linetype='dashed',size=0.1)+
        geom_segment(size=3.5,aes(color=dataset))+
        scale_color_manual(values = c('#876FAC','#AFB86E','#74C2B9','#23C373','#BDB042'))+
        theme_minimal()+
        labs(x='genomic coordinate on chr5 (Mb)',title = 'MYO10')+
        geom_segment(first_exons_loop,mapping = aes(x=start,xend=end,y=id+9,yend=id+9),size=5,color='#2F308C')+
        geom_segment(polyA_loop,mapping = aes(x=start,xend=end,y=id+6,yend=id+6),size=5,color='#BB4D51')+
        #geom_segment(dotted_isoforms_loop,mapping = aes(y=(id+id_add_isoforms+207)/2,x=start,xend=end,yend=(id+id_add_isoforms+207)/2),linetype='dashed',size=0.05)+
        #geom_segment(isoforms_loop,mapping = aes(x=start,xend=end,y=(id+id_add_isoforms+207)/2,yend=(id+id_add_isoforms+207)/2),size=2.5,color='#EA993E')+
        scale_x_continuous(labels = unit_format(unit = '',scale = 1e-6))+
        theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.background = element_blank(),axis.ticks.x = element_blank(),axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_text(size = 13, face = 'bold')))
dev.off()

#3.1
load('~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/Paper/figures/ext5_data.RData')


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
pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/ext_data_5_R_out.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs3_1a, 0.3, 0.66, 0.35, 0.33)+
  draw_plot(figs3_1_norm_exon_length, 0.65, 0.66, 0.35, 0.33)+
  #draw_plot(figs3_1b, 0.5, 0.66, 0.5, 0.33)+
  #draw_plot(fig_s3_legend, 0.33, 0.66, 0.2, 0.15)+
  draw_plot(figs3_1c, 0, 0.66, 0.3, 0.33)
dev.off()

#4.1 supp
#Expanded cor <-0.1 for for human
distances_full_neg_values$classif_for_plot_full <- factor(distances_full_neg_values$classif_for_plot_full, levels = c("-1 - -0.75","-0.75 - -0.5","-0.5 - -0.25","-0.25 - -0.1","-0.1 - 0.1","0.1 - 0.25","0.25 - 0.5","0.5 - 0.75","0.75 - 1"))

figs_4_1_A <- distances_full_neg_values %>% subset(value>0) %>% ggplot()+
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

figs4_1b <- ggplot(mapping=aes(x=FE_distance,y=gene_length))+
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

figs4_1c <- ggplot(mapping=aes(x=polyA_distance,y=gene_length))+
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

#Delta AUC mouse
figs4_1d <- ggplot(df_sup4.1c)+
  geom_boxplot(aes(x=variable,y=value,fill=variable),alpha=0.8,outlier.size = 0.1,notch=T)+
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_cowplot()+
  labs(y='ΔAUC R>0',fill='Exon type')+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),axis.line = element_line(colour = 'black', size = 1),axis.text=element_text(size=11),axis.title=element_text(size=13),axis.title.x = element_blank())+
  ylim(-0.3,0.3)+
  scale_fill_manual(values=c('#A2A2AA',met.brewer('Hiroshige')[10]))

x <- df_sup4.1c[df_sup4.1c$variable=='single FE or polyA',] %>% pull(value)
y <- df_sup4.1c[df_sup4.1c$variable=='multiple FE & polyA',] %>% pull(value)

t.test(x,y)

#Distances mouse
figs4_1e <- df_sup4.1d %>% subset(value>0) %>% ggplot()+
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

#fig_4ab_legend <- g.legend(fig4_a)
#fig4_a <- fig4_a+theme(legend.position = 'none')

df_4c$classif_final <- factor(df_4c$classif_final,levels=c('Non-PITA both','MaPITA mouse','MaPITA human','MaPITA both'),labels=c('no PITA','PITA mouse','PITA human','PITA both'))

figs4_1f <- df_4c[!df_4c$variable%in%c('Gene length','FE','polyA'),] %>% ggplot()+
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
figs4_1g <- ggplot(df_s_4g[df_s_4g$classif_final!='no PITA' & df_s_4g$classif_final!='Non-PITA both' & !is.na(df_s_4g$bin_for_plot),], aes(x = bin_for_plot, alpha=classif_final)) +
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

#Calculate the correlation for the delta TSS and TES between human and mouse
figs4_1h <- ggplot(df_s4.1_f,aes(x = pita_type, y = rho,alpha=pita_type)) +
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


png("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/Figures4_1AB.png",width=8.5,height=11,units = 'in',res = 300)
ggdraw() + 
  draw_plot(figs4_1b, 0, 0.66, 0.5, 0.33)+
  draw_plot(figs4_1c, 0.5, 0.66, 0.5, 0.33)
dev.off()


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/fig_sup4_1_rout.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs_4_1_A, 0, 0.66, 1, 0.33)+
  draw_plot(figs4_1d, 0, 0, 0.33, 0.33)+
  draw_plot(figs4_1e, 0.33, 0, 0.66, 0.33)
dev.off()

pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/fig_sup4_2_rout.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs4_1g, 0, 0.66, 1, 0.33)+
  draw_plot(figs4_1f, 0, 0.33, 0.5, 0.33)+
  draw_plot(figs4_1h, 0.5, 0.33, 0.5, 0.33)
dev.off()


#####
#Supp 5

#Ext 8

ext8_A <- micro_c_heatmap_insulations_hff %>% ggplot()+
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

ext8_A_legend <- g.legend(ext8_A)
ext8_A <- ext8_A+theme(legend.position = 'none')

micro_c_site_specific_insulations$treatment <- factor(micro_c_site_specific_insulations$treatment,levels = c('ctrl','pita'),labels = c('no PITA','PITA'))
micro_c_site_specific_insulations$isoform <- factor(micro_c_site_specific_insulations$isoform,levels = c('tss','tes'),labels = c('TSS','TES'))
micro_c_site_specific_insulations_esc <- micro_c_site_specific_insulations[micro_c_site_specific_insulations$cell_line=='esc',]
micro_c_site_specific_insulations_hff <- micro_c_site_specific_insulations[micro_c_site_specific_insulations$cell_line=='hff',]

ext8_f <- micro_c_site_specific_insulations_hff %>% ggplot()+
  geom_tile(aes(x=variable,y=ypos,fill=log(value)))+
  scale_fill_gradientn(
    colours = rev(met.brewer('Hiroshige',n = 100)[10:90]),
    na.value = 'white',
    breaks=c(log(1/1.5),0,log(1.5)),
    labels = c(0.67,1,1.5),
    limits = c(log(1/1.5), log(1.5)))+
  theme_cowplot()+
  scale_x_continuous(breaks = c(max(micro_c_heatmap_insulations$variable)/3,2*max(micro_c_heatmap_insulations$variable)/3),labels = c('site1','site2'))+
  labs(fill='obs/exp')+
  facet_grid(treatment~isoform)+
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text = element_text(size=9),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y =element_blank(),
        strip.text=element_blank(),
        axis.line = element_blank()
  )



average_insulation_hff$site <- factor(average_insulation_hff$site,levels = c("upstreamTSS","downstreamTSS"))
ext8_g <- ggplot(average_insulation_hff)+
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
        strip.text=element_text(size=10),
  )+
  background_grid(major="none",minor="y")


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/ext_data_8_Rout.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(ext8_A,0, 0.8, 1, 0.20)+
  draw_plot(ext8_A_legend,1, 0.66, 0.08, 0.23)+
  draw_plot(ext8_f,0, 0, 0.66, 0.20)+
  draw_plot(ext8_g,0.66, 0, 0.33, 0.20)
dev.off()

#Ext 9

ext9_B <- micro_c_site_specific_insulations_esc %>% ggplot()+
  geom_tile(aes(x=variable,y=ypos,fill=log(value)))+
  scale_fill_gradientn(
    colours = rev(met.brewer('Hiroshige',n = 100)[10:90]),
    na.value = 'white',
    breaks=c(log(1/1.5),0,log(1.5)),
    labels = c(0.67,1,1.5),
    limits = c(log(1/1.5), log(1.5)))+
  theme_cowplot()+
  scale_x_continuous(breaks = c(max(micro_c_heatmap_insulations$variable)/3,2*max(micro_c_heatmap_insulations$variable)/3),labels = c('site1','site2'))+
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

ext9_B_legend <- g.legend(ext9_B)
ext9_B <- ext9_B+theme(legend.position = 'none')


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/ext_data_9_Rout.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(ext9_B,0, 0, 0.66, 0.20)+
  draw_plot(ext9_B_legend,1, 0.15, 0.08, 0.23)
dev.off()




#Ext 10
figs5_1a <- isoform_freqs %>% 
  filter(Var3=='Alternative FE & polyA',Var4=='no PITA') %>% 
  ggplot(aes(x=factor(Var1,levels = c("upstream","middle point","downstream")),y=factor(Var2,levels = c("upstream","middle point","downstream"))))+
  geom_tile(aes(fill=normalized_counts))+
  labs(x='position FE',y='position TES',fill='proportion of\ngenes')+
  scale_fill_gradientn(#colors = c("#E76254", "#EF8A47", "#F7AA58", "white", "#72BCD5", "#528FAD", "#1E466E"),
                      colors = c("white", "#72BCD5", "#528FAD", "#1E466E"),
                       breaks = c(0,.25,  0.5),
                       labels = c(0,.25, 0.5),
                       limits = c(0, 0.5))+
  theme_cowplot()+
  theme(
        axis.line = element_blank(),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        axis.ticks = element_blank(),
        strip.text=element_text(size=10))

figs5_1a_legend <- g.legend(figs5_1a)
figs5_1a <- figs5_1a+theme(legend.position = 'none')




figs5_1b <- isoform_freqs %>% 
  filter(Var3=='Alternative FE & polyA',Var4=='PITA') %>% 
  ggplot(aes(x=factor(Var1,levels = c("upstream","middle point","downstream")),y=factor(Var2,levels = c("upstream","middle point","downstream"))))+
  geom_tile(aes(fill=normalized_counts))+
  labs(x='position TES',y='position FE',fill='proportion of genes')+
  scale_fill_gradientn(#colors = c("#E76254", "#EF8A47", "#F7AA58", "white", "#72BCD5", "#528FAD", "#1E466E"),
                      colors = c("white", "#72BCD5", "#528FAD", "#1E466E"),
                       breaks = c(0,.25,  0.5),
                       labels = c(0,.25, 0.5),
                       limits = c(0, 0.5))+
  theme_cowplot()+
  theme(axis.line = element_blank(),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        strip.text=element_text(size=10),
        axis.ticks = element_blank(),
        legend.position = 'none')


elong_index_per_length_bin$bin_log2 <- factor(elong_index_per_length_bin$bin_log2, levels = c("[1.14e+04,4.14e+04]","(4.14e+04,7.14e+04]","(7.14e+04,1.01e+05]","(1.01e+05,1.31e+05]","(1.31e+05,1.61e+05]","(1.61e+05,1.91e+05]","(1.91e+05,2.21e+05]","(2.21e+05,2.51e+05]","(2.51e+05,2.81e+05]","(2.81e+05,3.11e+05]","(3.11e+05,3.41e+05]","(3.41e+05,3.71e+05]","(3.71e+05,4.01e+05]","(4.01e+05,4.31e+05]","(4.31e+05,4.61e+05]","(5.51e+05,5.81e+05]"))

figs5_1c <- ggplot(elong_index_per_length_bin, aes(x = bin_log2, y = mean_value, fill = corr_classif)) +
  geom_bar(stat = "identity", position = position_dodge2(),alpha=0.8) +
  geom_errorbar(aes(ymin = mean_value - sem, ymax = mean_value + sem), 
                width = 0.25, position = position_dodge(width=0.9)) +
  theme_cowplot() +
  scale_y_log10() +
  labs(x = 'gene length bin', y = 'mean RNAPII velocity') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_fill_manual(values = c('#A2A2AA',met.brewer('Hiroshige')[10]))+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = 'black', size = 1),
        #panel.background = element_rect(color='lightgray'),
        axis.text = element_text(size=9), 
        axis.title = element_text(size=10),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=9),
        #strip.background=element_rect(fill='gray95'),
        strip.text=element_text(size=10))


ext10_d <- ggplot(mouse_polII_mut_wt, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = results_cor)) + 
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

ext10_d_e_legend <- g.legend(ext10_d)
ext10_d <- ext10_d+theme(legend.position = 'none')

ext10_e <- ggplot(mouse_polII_mut_slow, aes(x = AFEorder, y = ALEorder)) + geom_tile(aes(fill = results_cor)) + 
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


pdf("~/Dropbox (UMass Medical School)/PaiLab/collabprojects/TSS_TES/figures/paper/ext_data_10_Rout.pdf",width=8.5,height=11)
ggdraw() + 
  draw_plot(figs5_1a, 0, 0.66, 0.45, 0.33)+
  draw_plot(figs5_1b, 0.55, 0.66, 0.45, 0.33)+
  draw_plot(figs5_1a_legend,0.46,0.81,0.1,0.1)+
  draw_plot(figs5_1c, 0, 0.4, 1, 0.25)+
  draw_plot(ext10_d, 0, 0.1, 0.35, 0.264)+
  draw_plot(ext10_e, 0.55, 0.1, 0.35, 0.264)+
  draw_plot(ext10_d_e_legend,0.4,0.2,0.1,0.1)
dev.off()
























