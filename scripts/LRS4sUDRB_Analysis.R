library(ggplot2)
library(cowplot)
library(R.utils)
library(data.table)
library(tidyr)
library(dplyr)
library(ggridges)
library(ggbeeswarm)
library(MetBrewer)

##### functions #####

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

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

filterReads <- function(data, afereads, genecounts){
  data_geneinfo <- data %>% 
    # count number of total reads per gene and number of afes per gene
    group_by(gene) %>% mutate(n_reads = n(), n_afe = n_distinct(afe)) %>% 
    # count number of reads per afe
    group_by(gene, afe) %>% mutate(n_afe_reads = n())
  # filter for at least 10 reads/gene, 2+ afes, minimum number of reads/afe
  data_geneinfo <- subset(data_geneinfo, n_reads >= genecounts & n_afe >= 2 & n_afe_reads >= afereads)
  # using filtered data, count number of afes that meet read filter
  data_geneinfo <- data_geneinfo %>% group_by(gene) %>% mutate(n_afe_filtered = n_distinct(afe))
  # filter for 2+ afes after conditioning on afe read filter
  data_geneinfo <- subset(data_geneinfo, n_afe_filtered >= 2)
  return(data_geneinfo)
}

subsampleReads <- function(data, afereads, samples){
  genes <- unique(data$gene)
  datahere <- c()
  for(i in genes){
    genedata <- subset(data, gene == i)
    genedata <- genedata[order(genedata$afe),]
    afes <- sort(unique(genedata$afe))
    afenum <- c(1:length(afes))
    for(n in 1:length(afes)){
      afedata <- subset(genedata, afe == afes[n])
      afedata$afe_reorder <- afenum[n]
      for(s in c(1:samples)){
        afedata_sub <- sample_n(afedata, afereads)
        afedata_sub$sample <- s
        datahere <- rbind(datahere, afedata_sub)
      }
    }
  }
  return(datahere)
}

geneSummarize <- function(data){
	datahere <- data %>% group_by(gene, afe_reorder, sample) %>% 
	summarize(mean_mapped = mean(mapped_length), mean_read = mean(read_length),
		median_mapped = median(mapped_length), median_read = median(read_length),
		max_mapped = max(mapped_length), max_read = max(read_length),
		n_gene_reads = unique(n_reads), n_afe_reads = unique(n_afe_reads))
	return(datahere)
}

filterReadsDNA <- function(data, afereads, genecounts){
  data_geneinfo <- data %>% 
    # count number of total reads per gene and number of afes per gene
    group_by(gene_id_FE) %>% mutate(n_reads = n(), n_afe = n_distinct(FE_index)) %>% 
    # count number of reads per afe
    group_by(gene_id_FE, FE_index) %>% mutate(n_afe_reads = n())
  # filter for at least 10 reads/gene, 2+ afes, minimum number of reads/afe
  data_geneinfo <- subset(data_geneinfo, n_reads >= genecounts & n_afe >= 2 & n_afe_reads >= afereads)
  # using filtered data, count number of afes that meet read filter
  data_geneinfo <- data_geneinfo %>% group_by(gene_id_FE) %>% mutate(n_afe_filtered = n_distinct(FE_index))
  # filter for 2+ afes after conditioning on afe read filter
  data_geneinfo <- subset(data_geneinfo, n_afe_filtered >= 2)
  return(data_geneinfo)
}

subsampleReadsDNA <- function(data, afereads, samples){
  genes <- unique(data$gene_id_FE)
  datahere <- c()
  for(i in genes){
    genedata <- subset(data, gene_id_FE == i)
    genedata <- genedata[order(genedata$FE_index),]
    afes <- sort(unique(genedata$FE_index))
    afenum <- c(1:length(afes))
    for(n in 1:length(afes)){
      afedata <- subset(genedata, FE_index == afes[n])
      afedata$afe_reorder <- afenum[n]
      for(s in c(1:samples)){
        afedata_sub <- sample_n(afedata, afereads)
        afedata_sub$sample <- s
        datahere <- rbind(datahere, afedata_sub)
      }
    }
  }
  return(datahere)
}

getReadsGrouped <- function(bed, g, minlen = 200){
  # Remove reads with no readname
  bed_for_genes <- bed[!is.na(bed$read_name), ]
  #subset genes we are interested in 
  bed <- subset(bed_for_genes, gene==g & mapped_length >= minlen)
  
  # get the frequency of each read
  readfreq <- as.data.frame(table(bed$read_name))
  # add the order number of exons
  bed$nread <- readfreq$Freq[match(bed$read_name, readfreq$Var1)]
  bed <- bed %>% data.frame() %>% group_by(read_name) %>% dplyr::mutate(count = row_number())
  ## loop through number of introns
  counter = 0
  #Group reads by FE index, and order them in ascending genomic order 
  forloop <- unique(bed$afe)
  forloop_inorder <- sort(forloop)
  LRS_loop <- c()
  
  for(i in forloop_inorder){
    print(i)
    bed_here <- subset(bed, afe == i)
    #Correct read order for strand
    strand <- bed_here$strand[1]
    
    if (strand=='+'){
      if  (any(readfreq$Freq > 1)) { 
        start_temp <- aggregate(bed_here$start, by = list(bed_here$read_name), min) %>% `colnames<-`(c("read_name", "start"))
        end3prime_temp <- aggregate(bed_here$start, by = list(bed_here$read_name), max) %>% `colnames<-`(c("read_name", "start"))
        end_temp <- aggregate(bed_here$end, by = list(bed_here$read_name), max)%>% `colnames<-`(c("read_name", "end"))
        merged <- merge(start_temp,end3prime_temp,by="read_name")
        colnames(merged)<- c('read_name', 'start', 'start_old')
        merged <- merge(merged,end_temp,by="read_name")
      } else {
        
        merged <- data.frame(read_name=bed_here$read_name, 
                             start= bed_here$start, 
                             start_old= bed_here$start, 
                             end= bed_here$end)
      }
      colnames(merged)<- c('read_name', 'start', 'start_old', 'end')
      merged <- merged[order(merged$start_old,merged$start, merged$end),] %>% transform(id=match(read_name, unique(read_name)))
      
    }else{
      
      if  (any(readfreq$Freq > 1)) { 
        start_temp <- aggregate(bed_here$end, by = list(bed_here$read_name), min) %>% `colnames<-`(c("read_name", "start"))
        end_short_temp <- aggregate(bed_here$start, by = list(bed_here$read_name), min) %>% `colnames<-`(c("read_name", "old_end"))
        end_temp <- aggregate(bed_here$start, by = list(bed_here$read_name), max)%>% `colnames<-`(c("read_name", "end"))
        merged <- merge(start_temp,end_short_temp,by="read_name")
        colnames(merged)<- c('read_name', 'start', 'old_end')
        merged <- merge(merged,end_temp,by="read_name")
      } else {
        merged <- data.frame(read_name=bed_here$read_name, 
                             start= bed_here$end,
                             old_end= bed_here$start, 
                             end= bed_here$start)
      }
      colnames(merged)<- c('read_name', 'start', 'old_end', 'end')
      merged <- merged %>% mutate(length= merged$end-merged$start)
      merged <- merged[order(merged$start,merged$length, merged$end,decreasing = T),]%>% transform(id=match(read_name, unique(read_name)))
    }
    
    LRS_loop_here <- merge(merged[,c('read_name','id')],bed_here, by="read_name")
    LRS_loop_here$id <- LRS_loop_here$id + counter
    LRS_loop <- rbind(LRS_loop, LRS_loop_here)
    counter = length(unique(LRS_loop$read_name))
    print(counter)
  }
  
  # get introns
  
  countsegments <- as.data.frame(table(LRS_loop$read_name))
  splicedreads <- subset(countsegments, Freq >= 2)$Var1
  LRS_loop_spliced <- subset(LRS_loop, read_name %in% splicedreads)
  
  if  (any(readfreq$Freq > 1)) {
    dotted_LRS <- as.data.frame(aggregate(LRS_loop_spliced$end, by = list(LRS_loop_spliced$id), min))
    dotted_LRS$end <- as.data.frame(aggregate(LRS_loop_spliced$start, by = list(LRS_loop_spliced$id), max))[,2]
    colnames(dotted_LRS) <- c('id','start','end')
    dotted_LRS_expand <- data.frame(read_name = NA,
                                    id = dotted_LRS$id,
                                    chr=NA,
                                    start = dotted_LRS$start,
                                    end = dotted_LRS$end, 
                                    gene=NA,
                                    strand = NA,
                                    afe = NA,
                                    mapped_length = NA,
                                    nread = NA, 
                                    count = NA,
                                    type = "intron")
    
    # combine dataframes
    LRS_loop$type = "exon"
    LRS_combo <- rbind(LRS_loop, dotted_LRS_expand)
  } else {
    LRS_loop$type = "exon"
    LRS_combo <- LRS_loop
  }
  
  
  
  return(LRS_combo)
}


##### direct RNA data #####

direct_10m_stranded_rep1 <- fread("dRNA_rep1_K562_10m_nonrRNA_sorted_20250214_stranded_parsedReads.info.gz")
direct_10m_stranded_rep2 <- fread("K562_rep2_10m_intersectiwthAFEandALE_stranded.bam_parsedReads.info.gz")

bothreps_10m_stranded <- rbind(direct_10m_stranded_rep1, direct_10m_stranded_rep2)

### filter reads for gene expression, # of AFEs, and # of reads/afe

direct_10m_stranded_geneinfo_rep1 <- filterReads(direct_10m_stranded_rep1, 2, 10)
direct_10m_stranded_geneinfo_rep2 <- filterReads(direct_10m_stranded_rep2, 2, 10)

bothreps_10m_stranded_geneinfo <- filterReads(bothreps_10m_stranded, 2, 10)

### subsample reads

direct_10m_stranded_subsampled_rep1 <- subsampleReads(direct_10m_stranded_geneinfo_rep1, 2, 1)
direct_10m_stranded_subsampled_rep2 <- subsampleReads(direct_10m_stranded_geneinfo_rep2, 2, 1)

bothreps_10m_stranded_subsampled <- subsampleReads(bothreps_10m_stranded_geneinfo, 2, 1)

### add PITA classification
pitaclassif <- read.table("pita_classif.tsv",sep="\t",header=T)

direct_10m_stranded_subsampled_rep1$pita <- pitaclassif$classification[match(direct_10m_stranded_subsampled_rep1$gene, pitaclassif$gene_id)]
direct_10m_stranded_subsampled_rep2$pita <- pitaclassif$classification[match(direct_10m_stranded_subsampled_rep2$gene, pitaclassif$gene_id)]

bothreps_10m_stranded_subsampled$pita <- pitaclassif$classification[match(bothreps_10m_stranded_subsampled$gene, pitaclassif$gene_id)]

### summarize by gene

direct_10m_stranded_subsampled_gene_rep1 <- geneSummarize(direct_10m_stranded_subsampled_rep1)
direct_10m_stranded_subsampled_gene_rep1$pita <- pitaclassif$classification[match(direct_10m_stranded_subsampled_gene_rep1$gene, pitaclassif$gene_id)]

direct_10m_stranded_subsampled_gene_rep2 <- geneSummarize(direct_10m_stranded_subsampled_rep2)
direct_10m_stranded_subsampled_gene_rep2$pita <- pitaclassif$classification[match(direct_10m_stranded_subsampled_gene_rep2$gene, pitaclassif$gene_id)]

bothreps_10m_stranded_subsampled_gene <- geneSummarize(bothreps_10m_stranded_subsampled)
bothreps_10m_stranded_subsampled_gene$pita <- pitaclassif$classification[match(bothreps_10m_stranded_subsampled_gene$gene, pitaclassif$gene_id)]

### modify summarized both reps

# make NAs into no PITA
bothreps_10m_stranded_subsampled_gene$pita[which(is.na(bothreps_10m_stranded_subsampled_gene$pita))] <- "no PITA"

# combine adjacent AFEs together
bothreps_10m_stranded_subsampled_gene$afe_reorder_edit <- "1"
bothreps_10m_stranded_subsampled_gene$afe_reorder_edit[which(bothreps_10m_stranded_subsampled_gene$afe_reorder >= 2)] <- "2-3"
bothreps_10m_stranded_subsampled_gene$afe_reorder_edit[which(bothreps_10m_stranded_subsampled_gene$afe_reorder >= 4)] <- "4-5"
bothreps_10m_stranded_subsampled_gene$afe_reorder_edit[which(bothreps_10m_stranded_subsampled_gene$afe_reorder >= 6)] <- "6+"

table(bothreps_10m_stranded_subsampled_gene$afe_reorder_edit)
table(bothreps_10m_stranded_subsampled_gene$pita, bothreps_10m_stranded_subsampled_gene$afe_reorder_edit)

write.table(bothreps_10m_stranded_subsampled_gene, file="bothreps_directRNA_10mstranded_subsampled_gene.txt",sep="\t", quote=F, row.names=F, col.names=T)

### correlate replicates

# 2, 10: 3373 AFEs, 1462 genes, mean R = 0.7573591
# 5, 10: 1344 AFEs, 629 genes, mean R = 0.856845
#10, 20: 432 AFEs, 206 genes, mean R = 0.8881038

direct_10m_stranded_geneinfo_highconf_rep1 <- filterReads(direct_10m_stranded_rep1, 2, 10)
direct_10m_stranded_geneinfo_highconf_rep2 <- filterReads(direct_10m_stranded_rep2, 2, 10)

# all afe reads, no subsampling
direct_10m_stranded_gene_highconf_rep1 <- direct_10m_stranded_geneinfo_highconf_rep1 %>% group_by(gene, afe) %>%
											summarize(mean_mapped = mean(mapped_length),
													  median_mapped = median(mapped_length),
													  n_afe_reads = unique(n_afe_reads))
direct_10m_stranded_gene_highconf_rep2 <- direct_10m_stranded_geneinfo_highconf_rep2 %>% group_by(gene, afe) %>%
											summarize(mean_mapped = mean(mapped_length),
													  median_mapped = median(mapped_length),
													  n_afe_reads = unique(n_afe_reads))

direct_10m_stranded_gene_highconf_rep1$repgenename <- paste(direct_10m_stranded_gene_highconf_rep1$gene, direct_10m_stranded_gene_highconf_rep1$afe, sep="-")
direct_10m_stranded_gene_highconf_rep2$repgenename <- paste(direct_10m_stranded_gene_highconf_rep2$gene, direct_10m_stranded_gene_highconf_rep2$afe, sep="-")

genenames <- unique(c(direct_10m_stranded_gene_highconf_rep1$repgenename, direct_10m_stranded_subsampled_gene_highconf_rep2$repgenename))
rep1_match = match(genenames, direct_10m_stranded_gene_highconf_rep1$repgenename)
rep2_match = match(genenames, direct_10m_stranded_gene_highconf_rep2$repgenename)

combo_repcorr_gene <- data.frame(gene = unlist(lapply(strsplit(genenames, split="-"), "[", 1)),
                                 afe_reorder = unlist(lapply(strsplit(genenames, split="-"), "[", 2)),
                                 rep1_nreads = direct_10m_stranded_gene_highconf_rep1$n_afe_reads[rep1_match],
                                 rep2_nreads = direct_10m_stranded_gene_highconf_rep2$n_afe_reads[rep2_match],
                                 rep1_mean = direct_10m_stranded_gene_highconf_rep1$mean_mapped[rep1_match],
                                 rep2_mean = direct_10m_stranded_gene_highconf_rep2$mean_mapped[rep2_match],
                                 rep1_median = direct_10m_stranded_gene_highconf_rep1$median_mapped[rep1_match],
                                 rep2_median = direct_10m_stranded_gene_highconf_rep2$median_mapped[rep2_match])

cor(combo_repcorr_gene$rep1_mean, combo_repcorr_gene$rep2_mean, use = "complete.obs")
cor(combo_repcorr_gene$rep1_median, combo_repcorr_gene$rep2_median, use = "complete.obs")

dim(subset(combo_repcorr_gene, !is.na(rep1_mean) & !is.na(rep2_mean))) #afes
length(unique(subset(combo_repcorr_gene, !is.na(rep1_mean) & !is.na(rep2_mean))$gene)) #genes

ggplot(combo_repcorr_gene, aes(x=rep1_mean, y=rep2_mean)) + 
  geom_point(size=3, alpha=0.25) + 
  scale_x_log10() + scale_y_log10() + 
  theme_cowplot()

ggplot(combo_repcorr_gene, aes(x=rep1_median, y=rep2_median)) + 
  geom_point(size=3, alpha=0.25) + 
  scale_x_log10() + scale_y_log10() + 
  theme_cowplot()

##### cDNA data #####

cDNA_unstranded <- fread("cDNA_bothtimepoints_allreads_NOTsummarizedbygene.tsv.gz")
cDNA_unstranded_10m <- subset(cDNA_unstranded, sample="10m")

# filter reads
cDNA_unstranded_10m_geneinfo <- filterReadsDNA(cDNA_unstranded_10m, 2, 10)

# subsample reads
cDNA_unstranded_10m_subbed <- subsampleReadsDNA(cDNA_unstranded_10m_geneinfo, 2, 1)

table(cDNA_unstranded_10m_subbed$afe_reorder)
table(cDNA_unstranded_10m_subbed$afe_reorder, cDNA_unstranded_10m_subbed$PITA)

ggplot(subset(cDNA_unstranded_10m_subbed, genomic_length >= 100 & afe_reorder <= 5), aes(x=factor(afe_reorder), y = genomic_length)) + 
  geom_violin() + geom_boxplot(width=0.15, notch=T, outliers = F) +
  scale_y_log10(breaks=c(1e2,1e3,1e4,1e5),labels=c("0.1kb","1kb","10kb","100kb")) +
  labs(x="AFE ordinal position",y="genomic distance traveled") +
  theme_cowplot()

##### PLOT EXAMPLE FROM DIRECT RNA DATA #####

### load split bed files of reads
reads_complete_split_rep1 <- fread("plot_raw_reads/dRNA_rep1_K562_10m_nonrRNA_sorted_20250214_stranded_parsedReads_split.bed.gz", sep = "\t")
reads_complete_split_rep2 <- fread("plot_raw_reads/dRNA_rep2_K562_10m_nonrRNA_sorted_20250323_stranded_parsedReads_split.bed.gz", sep = "\t")
colnames(reads_complete_split_rep1) <- colnames(reads_complete_split_rep2) <- c('chr','start','end','readname','score', 'strand')

### load read info
reads_withinfo_rep1 <- as.data.frame(read.table("plot_raw_reads/dRNA_rep1_K562_10m_nonrRNA_sorted_20250214_stranded_parsedReads.info.gz",header = TRUE, sep="\t",stringsAsFactors=FALSE))
reads_withinfo_rep2 <- as.data.frame(read.table("plot_raw_reads/dRNA_rep2_K562_10m_nonrRNA_sorted_20250323_stranded_parsedReads.info.gz",header = TRUE, sep="\t",stringsAsFactors=FALSE))

reads_complete_rep1 <- merge(reads_complete_split_rep1,reads_withinfo_rep1, by="readname", all=TRUE)
reads_complete_rep2 <- merge(reads_complete_split_rep2,reads_withinfo_rep2, by="readname", all=TRUE)

colnames(reads_complete_rep1) <- colnames(reads_complete_rep2) <- c("read_name",  "chr" , "start" , "end",  "score" , "strand", "gene", "afe", "mapped_length", "read_length" )
bothreads_complete <- rbind(reads_complete_rep1, reads_complete_rep2)

subset_reads_withindex<- bothreads_complete[,c("read_name","chr", "start","end","strand","gene","afe", "mapped_length")]

### load feature coordinates
exon_coordinates <- fread("plot_raw_reads/HG38_exons.bed.gz", sep ="\t")
colnames(exon_coordinates) <- c("chr", "start", "end", "gene", "feature", "strand")
FE_coordinates<- fread("../DRB-LRS_python/first_exonstmp.bed.gz", sep = "\t")
colnames(FE_coordinates) <- c("chr", "start", "end", "gene", "feature", "strand")

# plot examples
install.packages("paletteer")
library(paletteer)

plotLRSexample <- function(genehere, exon_coordinates, FE_coordinates, read_coordinates){
  # gene title
  gene_title_withstrand <- paste(genehere, read_coordinates$strand[1], sep=" | strand ")
  # exon, FE coordinates for gene
  exon_coord_loop<- filter(exon_coordinates, gene == genehere)
  FE_coord_loop<- filter(FE_coordinates, gene == genehere)
  # format reads for plotting
  read_coordinates$afe <- factor(read_coordinates$afe, levels=sort(unique(read_coordinates$afe)))
  FE_coord_loop_parsed <- subset(FE_coord_loop, feature %in% unique(read_coordinates$afe))
  # don't plot genes that have 12+ FEs used
  if(length(unique(read_coordinates$afe)) > 12){ next }
  # plot genes
  print(ggplot(read_coordinates)+
    # plot reads
    geom_segment(data = subset(read_coordinates, type=="intron"), mapping=aes(x=start, xend=end, y=-1*id, yend=-1*id),color="grey75", linetype="dashed", linewidth=1)+
    geom_segment(data = subset(read_coordinates, type=="exon"), aes(x= start, xend=end, y=-1*id, yend=-1*id, color=factor(afe)), linewidth=1.5)+
    geom_segment(data=exon_coord_loop, aes(x= start , xend= end, y= 1, yend= 1), color="gray", linewidth= 3)+
    geom_segment(data=FE_coord_loop_parsed, aes(x= start , xend= end, y= 1, yend= 1, color=factor(feature)), linewidth= 3) +
    #reduce xlimits, as most reads are only starting to be transcribed
    coord_cartesian(c(min(read_coordinates$start)-1000, max(read_coordinates$end)+1000))+
    labs(x = " ", y = " ", title=gene_title_withstrand, color = "Genomic order of FE")+
    scale_color_manual(values=paletteer_d("PrettyCols::Summer"))+
    theme_cowplot() +
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
}

read_coordinates <- getReadsGrouped(subset_reads_withindex, gene, minlen = 200)
plotLRSexample(genehere, exon_coordinates, FE_coordinates, read_coordinates)

bestgenes <- c("ENSG00000072864", "ENSG00000104824", "ENSG00000104957", "ENSG00000105323", "ENSG00000147677", "ENSG00000159131", "ENSG00000160959", "ENSG00000166925",
               "ENSG00000171552", "ENSG00000089009", "ENSG00000126561")

candidategenes <- c("ENSG00000048405", "ENSG00000063245", "ENSG00000064607", "ENSG00000084207", "ENSG00000092531", "ENSG00000099194", "ENSG00000100335", "ENSG00000100991",
                    "ENSG00000105699", "ENSG00000105819", "ENSG00000108582", "ENSG00000112584", "ENSG00000119760", "ENSG00000120253", "ENSG00000124535", "ENSG00000131236",
                    "ENSG00000137312", "ENSG00000147162", "ENSG00000156515", "ENSG00000160211", "ENSG00000171853", "ENSG00000183010", "ENSG00000204227", "ENSG00000065183", 
                    "ENSG00000072310", "ENSG00000103275", "ENSG00000109534")

# boxplot for gene
g = "ENSG00000104824"
genehere <- subset(bothreps_10m_stranded_geneinfo, gene == g & mapped_length >= 100)
titletext = paste0(g, "\nreads = ", nrow(genehere))
ggplot(subset(genehere, afe != 4), aes(x=factor(afe), y = mapped_length)) + 
  geom_boxplot(notch=T) +
  scale_x_discrete(labels=c("1","2")) +
  labs(x="AFE ordinal position",y="genomic distance (nt)") +
  scale_y_log10() + theme_cowplot()


