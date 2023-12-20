##############
#This scripts analyzes ENCODE long-read seq data to find TSS-TES correlations and the analysis of sequencing depth with subsampling.

#Input: bam2bed split files where each feature of each read is one individual row
#Output: Spearman correlations between TSS-TES, filtered reads matching FE and polyA databases or any other databases that want to be checked

#Requirements: besides the loaded packages, Bedtools intersect is called from the script

library(tidyverse)
library(stringr)
library(data.table)

setDTthreads(threads = 8)

#####
#Gets first and last exons
firs_last_fun <- function(reads){
  reads_pos <- reads[reads$strand=='+',]
  reads_neg <- reads[reads$strand=='-',]
  
  #Find the first feature of each read (first exon)
  reads_all_pos_first <- reads_pos[reads_pos[, .I[start == min(start)], by=read_name]$V1] %>% as.data.frame()
  reads_all_pos_last <- reads_pos[reads_pos[, .I[start == max(start)], by=read_name]$V1] %>% as.data.frame()
  
  #Find the last feature of each read (last exon)
  reads_all_neg_first <- reads_neg[reads_neg[, .I[end == max(end)], by=read_name]$V1] %>% as.data.frame()
  reads_all_neg_last <- reads_neg[reads_neg[, .I[end == min(end)], by=read_name]$V1] %>% as.data.frame()
  
  reads_all_first <- rbind(reads_all_pos_first,reads_all_neg_first)
  reads_all_last <- rbind(reads_all_pos_last,reads_all_neg_last)
  
  return(list(reads_all_first=reads_all_first,reads_all_last=reads_all_last))
}

#Runs intersect for read start and ends and the references
intersect_function <- function(reads_all_list){
  
  #Split starts and ends from firs_last_fun
  read_starts <- reads_all_list[["reads_all_first"]][,c('chr','start','end','read_name','gene_id','strand')]
  read_ends <- reads_all_list[["reads_all_last"]][,c('chr','start','end','read_name','gene_id','strand')]
  
  # Write the location information to a temp file
  fwrite(read_starts, "read_start_tmp.bed", sep = "\t", append = FALSE, col.names = FALSE,row.names = F,quote = F)
  fwrite(read_ends, "read_end_tmp.bed", sep = "\t", append = FALSE, col.names = FALSE,row.names = F,quote = F)
  
  # Run bedtools to extract the sequence -> reeplace first_exons_reference.bed & polyA_reference.bed with the desired databases
  system("bedtools intersect -s -wao -a read_start_tmp.bed -b first_exons_reference.bed > starts_intersect.bed")
  system("bedtools intersect -s -wao -a read_end_tmp.bed -b polyA_reference.bed > ends_intersect.bed")
  
  starts_intersect <- fread('starts_intersect.bed',sep = '\t') %>% as.data.frame()
  ends_intersect <- fread('ends_intersect.bed',sep = '\t') %>% as.data.frame()
  
  #Classify reads according to their FE and LE matching the databases
  reads_no_start <- subset(starts_intersect,V10=='.') #Don't match first
  reads_no_end <- subset(ends_intersect,V10=='.') #Don't match polyA
  reads_matching_start <- starts_intersect[starts_intersect$V5==starts_intersect$V10,] #Match first
  reads_matching_end <- ends_intersect[ends_intersect$V5==ends_intersect$V10,] #Match polyA
  reads_matching_end <- setDT(reads_matching_end)[reads_matching_end[, .I[V11 == max(V11)], by=V4]$V1] %>% as.data.frame() #Keep the most downstream polyA in case there's overlap between peaks (could happen for tandemUTRs)
  
  #Merge reads using their names if they passed FE and polyA filters
  reads_passed_filters <- merge(reads_matching_start[,c('V1','V2','V3','V4','V5','V6','V11')],reads_matching_end[,c('V2','V3','V4','V10','V11')],by='V4')
  reads_passed_filters <- reads_passed_filters[reads_passed_filters$V5==reads_passed_filters$V10,]
  colnames(reads_passed_filters) <- c('read_name','chr','start_FE','end_FE','gene_id_FE','strand','FE_index','start_LE','end_LE','gene_id_LE','polyA_index')
  
  
  return(list(reads_no_start=reads_no_start,reads_no_end=reads_no_end,reads_passed_filters=reads_passed_filters))
  
}

#####
#Get first exon reference (using HG38 coordinates)
first_exons_hitindex <- read.table('first_exons_hitindex_gtex_hg38.txt',header = T)

#Assign a FE ordinal position (id) per gene going from 1-n of FE (being 1 the most upstream FE)
first_exons_hitindex_pos <- first_exons_hitindex %>% subset(strand=='+') %>% dplyr::arrange(start) %>% group_by(gene_name) %>% dplyr::mutate(id = row_number(gene_name)) %>% as.data.frame()
first_exons_hitindex_neg <- first_exons_hitindex %>% subset(strand=='-') %>% dplyr::arrange(desc(end)) %>% group_by(gene_name) %>% dplyr::mutate(id = row_number(gene_name)) %>% as.data.frame()
first_exons_hitindex <- rbind(first_exons_hitindex_pos,first_exons_hitindex_neg)
first_exons_hitindex <- first_exons_hitindex[,c('chr','start','end','gene_name','id','strand')]
fwrite(first_exons_hitindex, "first_exons_reference.bed", sep = "\t", append = FALSE, col.names = FALSE,row.names = F,quote = F)

#Get polyA reference, for example using https://www.polyasite.unibas.ch/atlas
polya_db <- read.table('polyA_db.bed',header = T)
#Assign a polyA ordinal position (id) per gene going from 1-n of polyA (being 1 the most upstream polyA)
polya_db <- polya_db[,c('chr','start','end','gene_name','id','strand')]
polya_db_pos <- polya_db %>% subset(strand=='+') %>% dplyr::arrange(start) %>% group_by(gene_name) %>% dplyr::mutate(id = row_number(gene_name)) %>% as.data.frame()
polya_db_neg <- polya_db %>% subset(strand=='-') %>% dplyr::arrange(desc(end)) %>% group_by(gene_name) %>% dplyr::mutate(id = row_number(gene_name)) %>% as.data.frame()
polya_db <- rbind(polya_db_pos,polya_db_neg)
polya_db <- polya_db[,c('chr','start','end','gene_name','id','strand')]
fwrite(polya_db, "polyA_reference.bed", sep = "\t", append = FALSE, col.names = FALSE,row.names = F,quote = F)

#####
#Initiate variables

#All bam files to analyze, in our case one per sample
files_all <-  list.files(path = './all_samples/',pattern = '*bed.gz')

spearman_per_sample <- list() #Correlations per gene per sample
gene_list_no_FE <- list() #Reads not matching FE per sample
gene_list_no_polyA <- list() #Reads not matching polyA per sample
reads_across_datasets <- list() #Reads that passed all filters per sample
reads_across_datasets_3utr <- list() #Reads only covering one feature per sample (single exon genes or non spliced)
truncation_metrics_list <- list() #Metrics per sample
file_counter <- 1 #Just a counter

for (i in 1:length(files_all)){
  print(paste('Starting sample',i,'of',length(files_all)))
  
  reads_all <- fread(paste('./all_samples/',files_all[i],sep=''),fill = T,sep = '\t') %>% as.data.frame()
  reads_all <- reads_all[,c(1,2,3,4,6,10)]
  
  file_name <- gsub("\\..*","",files_all[i])
  colnames(reads_all) <- c('chr','start','end','read_name','strand','gene_id')
  reads_all$end <- as.numeric(reads_all$end)
  reads_all$start <- as.numeric(reads_all$start)
  
  #Find reads with only one feature (not covering a splice site, starting and ending on the same exon). This is very important an a normal truncation type on long reads, especially on ONT. It's recommended to run this for all bed files: for i in *bed.gz;do gunzip -cd $i | cut -f 4 | sort | uniq -c > counts_per_read/$i;done.

  counts_per_read <- fread(paste('./counts_per_read/',files_all[i],sep=''),col.names = c('Freq','read_name')) %>% as.data.frame()
  terminal_reads <- counts_per_read[counts_per_read$Freq==1,'read_name'] #Reads with only one feature are terminal/truncated
  non_terminal_reads <- counts_per_read[counts_per_read$Freq!=1,'read_name']
  reads_across_datasets_3utr[[file_name]] <- reads_all[reads_all$read_name%in%terminal_reads,]
  reads_all <- reads_all[reads_all$read_name%in%non_terminal_reads,]
  
  #Split read starts and ends to make bed file 
  reads_all_list <- firs_last_fun(reads = setDT(reads_all))
  
  #Run intersect and store outputs
  intersect_function_out <- intersect_function(reads_all_list)
  gene_list_no_FE[[file_name]] <- intersect_function_out[["reads_no_start"]] #Reads not matching FE
  gene_list_no_polyA[[file_name]] <- intersect_function_out[["reads_no_end"]] #Reads not matching polyA
  reads_passed_filters <- intersect_function_out[["reads_passed_filters"]] #Reads that passed all filters
  
  #Get filtering metrics
  counts_per_gene_raw <- table(reads_all$V10) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_no_filter"))
  terminal_truncation <- table(reads_across_datasets_3utr[[file_name]]$gene_id) %>% as.data.frame()   %>% `colnames<-`(c("gene_id", "freq_terminal_reads"))
  truncation_5 <- table(reads_no_start$V5) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_5end_filter"))
  truncation_3 <- table(reads_no_end$V5) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_3end_filter"))
  counts_per_gene_after_filters <- table(reads_passed_filters$gene_id_FE) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "final_freq"))
  truncation_metrics <- merge(merge(merge(merge(counts_per_gene_raw,terminal_truncation,by='gene_id',all.x=T),truncation_5,by='gene_id',all.x=T),truncation_3,by='gene_id',all.x=T),counts_per_gene_after_filters,,by='gene_id',all.x=T)
  truncation_metrics_list[[file_name]] <- truncation_metrics
  
  #Keep reads with TSS and TES pairs (using their ordinal position) that are used at least three times for each gene. This filter is highly important to remove spurious reads
  reads_passed_filters <- reads_passed_filters %>% dplyr::group_by(gene_id_FE,FE_index) %>% dplyr::mutate(TSS_counts=n())%>% filter(TSS_counts>2)
  reads_passed_filters <- reads_passed_filters %>% dplyr::group_by(gene_id_FE,polyA_index) %>% dplyr::mutate(polyA_counts=n()) %>% filter(polyA_counts>2) %>% as.data.frame()
  
  #Store filtered reads on a list
  reads_across_datasets[[file_name]] <- reads_passed_filters
  
  #Make vectors of genes using AFE & ALE or Unique F | L exons
  data_count_1 <- aggregate(data = reads_passed_filters,FE_index ~ gene_id_FE,function(x) length(unique(x)))
  data_count_1 <- merge(data_count_1,aggregate(data = reads_passed_filters,polyA_index ~ gene_id_FE,function(x) length(unique(x))),by='gene_id_FE')
  alt_genes <- subset(data_count_1,FE_index>1 & polyA_index>1) %>% pull(gene_id_FE)
  unique_genes <- subset(data_count_1,FE_index==1 | polyA_index==1) %>% pull(gene_id_FE)
  
  #Calculate Spearman Rs
  spearman_all <- by(reads_passed_filters, reads_passed_filters$gene_id_FE, FUN = function(X) cor(X$start_FE, X$end_LE, method = "spearman")) #From the start of the first exon to the end of the last
  spearman_all <- data.frame(gene_id = dimnames(spearman_all)[[1]],corr = as.vector(spearman_all))
  
  #Assign alternative/unique exons
  spearman_all$exon_type[spearman_all$gene_id %in% alt_genes] <- 'Alternative FE & polyA'
  spearman_all$exon_type[spearman_all$gene_id %in% unique_genes] <- 'Unique FE or polyA'
  
  #Get counts per gene used for Spearman
  reads_all_counts <- counts_per_gene_after_filters %>% dplyr::mutate(cpm=10^6*final_freq/sum(final_freq))
  spearman_all <- merge(spearman_all,reads_all_counts,by = 'gene_id') #Add counts (after filtering) information
  
  #Add useful information
  spearman_all <- merge(spearman_all,data_count_1,by.x='gene_id',by.y='gene_id_FE')
  
  colnames(spearman_all)[c(6,7)] <- c('number_of_used_FE','number_of_used_polyA')
  spearman_per_sample[[file_name]] <- data.frame(spearman_all,sample=file_name)
}













