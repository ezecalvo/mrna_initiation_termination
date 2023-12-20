library(tidyr)
library(dplyr)
library(reshape)

##### data list for gtex #####
setwd("/projectnb/encore/carroll8/PITA/gtex/original_exons/")
getwd()
AFEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/gtex/original_exons/", pattern = "conservative.AFEPSI", recursive = FALSE)
ALEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/gtex/original_exons/", pattern = "conservative.ALEPSI", recursive = FALSE)
AFEPSI_list <- lapply(AFEPSI, read.delim, sep = "\t")
ALEPSI_list <- lapply(ALEPSI, read.delim, sep = "\t")
desired_order <-  c(8,9,3,14,13,20,6,7,16,18,17,22,10,19,1,21,4,12,5,11,2,15) #order tissues in defined manner for downstream processing 
tissue <- c("brain", "breast", "colon", "heart", "kidney", "liver", "lung", "ovary", "prostate", "stomach", "testis")
# Reorder the list of dataframes
AFEPSI_list <- AFEPSI_list[desired_order]
ALEPSI_list <- ALEPSI_list[desired_order]

### FIGURE1B ###

#code for determining mean number of FE/LEs per gene across tissues 
FELEpergene_upto4 <- function(a,b){
#cut off of calling FE/LE exon is PSI =< .05 
a <- a[a$AFEPSI >= .05,]
b <- b[b$ALEPSI >= .05,]
#group genes by # of FE/LE, group genes with 4 or more 
atable <- as.data.frame(table(a$gene))
btable <- as.data.frame(table(b$gene))
atable$Freq <- ifelse(atable$Freq > 4, 4, atable$Freq) 
btable$Freq <- ifelse(btable$Freq > 4, 4, btable$Freq) 
#calculate frequency of genes per # of FE/LEs
FE <- atable %>% group_by(Freq) %>% count("Freq") %>% select(-1)
LE <- btable %>% group_by(Freq) %>% count("Freq") %>% select(-1)
values <- rbind(FE,LE)
}

#apply function to both lists of dataframes and calculate mean, std_dev, and SEM
samples_df <- as.data.frame(mapply(FELEpergene_upto4, AFEPSI_list, ALEPSI_list))
samples_df$mean <- samples_df %>% rowMeans() 
samples_df$std_dev <- apply(samples_df, 1, sd)
samples_df <- samples_df %>% mutate(SEM = (samples_df$std_dev/ sqrt(22))) # fill in number of observations
FIGURE1Bdata <- samples_df
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
save(FIGURE1Bdata, file = "FIGURE1B_Data.RData")

###### FIGURE1C #####

#here, we will create 2 dfs- first will identify # FE/LEs per gene, second will be # observations per combination
AFEvsALE <- function(a, b) {
  #remove exons with less than 0.05 PSI values and then subset for genes with alternative AFEs and ALEs
  a <- a[a$AFEPSI >= 0.05,]
  b <- b[b$ALEPSI >= 0.05,]
  a_table <- as.data.frame(table(a$gene))
  a_table <- a_table[a_table$Freq != 1,]
  b_table <- as.data.frame(table(b$gene))
  b_table <- b_table[b_table$Freq != 1,]
  AFEALEgenes <- merge(a_table, b_table, by = "Var1", all = FALSE)
  AFEALEgenes <- AFEALEgenes[AFEALEgenes$Freq.x <= 6, ]
  AFEALEgenes <- AFEALEgenes[AFEALEgenes$Freq.y <= 6, ]
  AFEALEgenes_count <- AFEALEgenes %>%group_by(Freq.x, Freq.y) %>% tally()
  pearson_cor <- cor(AFEALEgenes$Freq.x, AFEALEgenes$Freq.y, method = c("pearson"))
  rounded_pearson_cor <- round(pearson_cor, digits = 3) 
  return(list(AFEALEgenes, AFEALEgenes_count, rounded_pearson_cor))
}

AFE_vs_ALE <- as.data.frame(mapply(AFEvsALE, AFEPSI_list, ALEPSI_list)) 
FIGURE1Cdata <- AFEvsALE(AFEPSI_list[[16]], ALEPSI_list[[16]])
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
save(FIGURE1Cdata, file = "FIGURE1C_Data.RData")


##### FIGURE 1D ####
#we want to extract the pearson correlation determined above in order to summarize pearson correlations across tissues 
pearsons <- AFE_vs_ALE[3,]
pearsons <- as.data.frame(t(pearsons))
pearsons$group <- rep(1:11, each = 2)
names(pearsons)[1] <- "value"
pearsons$value <- as.numeric(pearsons$value)
mean <- aggregate(pearsons, by = list(group = pearsons$group), FUN = mean)
cor <- mean[, -c(1,3)]
tissue <- c("brain", "breast", "colon", "heart", "kidney", "liver", "lung", "ovary", "prostate", "stomach", "testis")
data <- data.frame(tissue,cor)
row.names(data) <- data$tissue
data <- data[-1]
data
FIGURE1Ddata <- data
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
save(FIGURE1Ddata, file = "FIGURE1D_Data.RData")

### Figure 1EA ###

result_list <- list()
PSIbyGO_upto5 <- function(a,b){
  a <- a[a$AFEPSI >= .05,]
  b <- b[b$ALEPSI >= .05,]
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  a_5 <- a_table[a_table$Freq >=2 & a_table$Freq <=5,]
  b_5 <- b_table[b_table$Freq >=2 & b_table$Freq <=5,]
  genes <- merge(a_5, b_5, by = "Var1", all = FALSE)
  a_5 <- a %>%filter(gene %in% genes$Var1)
  b_5 <- b %>%filter(gene %in% genes$Var1)
  a_5 <- separate(a_5,exon,into=c("chromosome", "start", "end"))
  b_5 <- separate(b_5,exon,into=c("chromosome", "start", "end"))
  a_5_plus <- subset(a_5, strand == "+")
  b_5_plus <- subset(b_5, strand == "+")
  a_5_minus <- subset(a_5, strand == "-")
  b_5_minus <- subset(b_5, strand == "-")
  #order brain AFEplus ascending 
  a_5_plus <- a_5_plus[order(a_5_plus$start),]
  #order brain AFEE plus ascending 
  b_5_plus <- b_5_plus[order(b_5_plus$start),]
  #order brain AFEminus ascending 
  a_5_minus <- a_5_minus[order(a_5_minus$start, decreasing = TRUE),]
  #order brain AFEE minus ascending 
  b_5_minus <- b_5_minus[order(b_5_minus$start, decreasing = TRUE),]
  #add in the genomic order of the brain exons 
  a_5_plus <-  transform(a_5_plus, Order = ave(1:nrow(a_5_plus), gene,FUN = seq_along))
  b_5_plus <-  transform(b_5_plus , Order = ave(1:nrow(b_5_plus ), gene,FUN = seq_along))
  a_5_minus <-  transform(a_5_minus, Order = ave(1:nrow(a_5_minus), gene,FUN = seq_along))
  b_5_minus <-  transform(b_5_minus, Order = ave(1:nrow(b_5_minus), gene,FUN = seq_along))
  ##
  AFEs <- rbind(a_5_plus, a_5_minus)
  ALEs <- rbind(b_5_plus, b_5_minus)
  #sort by genomic order
  go_AFE_list <- split(AFEs, f = AFEs$Order)
  go_ALE_list <- split(ALEs, f = ALEs$Order)
  #master data list with all pair-wise combinations 
  for (i in 1:5) {
    for (j in 1:5) {
      result_list[[paste0("x_", i, "v", j)]] <- merge(go_AFE_list[[i]], go_ALE_list[[j]], by = "gene", all = FALSE)
    }
    combined_df <- do.call(rbind, result_list)
    print(combined_df)
  }
  return(combined_df)
}

alldata <- data.frame()
data <- data.frame()
start_id <- 0

# Loop through each pair of dataframes in the list
for (i in 1:(length(AFEPSI_list))) {
  # Extract the current and next dataframes
  data <- PSIbyGO_upto5(AFEPSI_list[[i]], ALEPSI_list[[i]])
  start_id <- start_id + length(unique(data$gene))
  data$gene <- as.numeric(factor(data$gene, levels = unique(data$gene))) + start_id
  alldata <- rbind(alldata, data)
  #print(start_id)
  #print(data)
}

alldata

results_cor_list <- vector("list", length = 25)

for (i in 1:5) {
  for (j in 1:5) {
    subset_df <- alldata[alldata$Order.x == i & alldata$Order.y == j, ]
    results_cor_list[[(i - 1) * 5 + j]] <- cor(subset_df$AFEPSI, subset_df$ALEPSI)
  }
}


list <-results_cor_list
FIGURE1EAdata <- unlist(list)
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
save(FIGURE1EAdata, file = "FIGURE1EA_Data.RData")


#### FIGURE 1EB ####
result_list <- list()
PSIbyGO_exactly3 <- function(a,b){
  a <- a[a$AFEPSI >= .05,]
  b <- b[b$ALEPSI >= .05,]
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  a_3 <- a_table[a_table$Freq ==3,]
  b_3 <- b_table[b_table$Freq ==3,]
  genes <- merge(a_3, b_3, by = "Var1", all = FALSE)
  print(genes)
  a_3 <- a %>%filter(gene %in% genes$Var1)
  b_3 <- b %>%filter(gene %in% genes$Var1)
  a_3 <- separate(a_3,exon,into=c("chromosome", "start", "end"))
  b_3 <- separate(b_3,exon,into=c("chromosome", "start", "end"))
  a_3_plus <- subset(a_3, strand == "+")
  b_3_plus <- subset(b_3, strand == "+")
  a_3_minus <- subset(a_3, strand == "-")
  b_3_minus <- subset(b_3, strand == "-")
  #order brain AFEplus ascending 
  a_3_plus <- a_3_plus[order(a_3_plus$start),]
  #order brain AFEE plus ascending 
  b_3_plus <- b_3_plus[order(b_3_plus$start),]
  #order brain AFEminus ascending 
  a_3_minus <- a_3_minus[order(a_3_minus$start, decreasing = TRUE),]
  #order brain AFEE minus ascending 
  b_3_minus <- b_3_minus[order(b_3_minus$start, decreasing = TRUE),]
  #add in the genomic order of the brain exons 
  a_3_plus <-  transform(a_3_plus, Order = ave(1:nrow(a_3_plus), gene,FUN = seq_along))
  b_3_plus <-  transform(b_3_plus , Order = ave(1:nrow(b_3_plus ), gene,FUN = seq_along))
  a_3_minus <-  transform(a_3_minus, Order = ave(1:nrow(a_3_minus), gene,FUN = seq_along))
  b_3_minus <-  transform(b_3_minus, Order = ave(1:nrow(b_3_minus), gene,FUN = seq_along))
  ##
  AFEs <- rbind(a_3_plus, a_3_minus)
  ALEs <- rbind(b_3_plus, b_3_minus)
  #
  #sort by genomic order
  go_AFE_list <- split(AFEs, f = AFEs$Order)
  go_ALE_list <- split(ALEs, f = ALEs$Order)
  #master data list with all pair-wise combinations 
  for (i in 1:3) {
    for (j in 1:3) {
      result_list[[paste0("x_", i, "v", j)]] <- merge(go_AFE_list[[i]], go_ALE_list[[j]], by = "gene", all = FALSE)
    }
    combined_df <- do.call(rbind, result_list)
  }
  return(combined_df)
}

# 
# # Loop through each pair of dataframes in the list
# for (i in 1:(length(AFEPSI_list))) {
#   # Extract the current and next dataframes
#   data <- PSIbyGO_exactly3(AFEPSI_list[[i]], ALEPSI_list[[i]])
#   # start_id <- start_id + length(unique(data$gene))
#   # data$gene <- as.numeric(factor(data$gene, levels = unique(data$gene))) + start_id
#   alldata <- rbind(alldata, data)
#   #print(start_id)
#   #print(data)
# }


alldata <- data.frame()
data <- data.frame()
start_id <- 0

for (i in 1:length(AFEPSI_list)) {
  # Use tryCatch to handle potential errors in the loop
  tryCatch({
    # Extract the current and next dataframes
    data <- PSIbyGO_exactly3(AFEPSI_list[[i]], ALEPSI_list[[i]])
    # Append the data to the result dataframe
    alldata <- rbind(alldata, data)
  }, error = function(e) {
    # Print an error message if an error occurs
    message(paste("Error occurred in iteration", i, ": ", e$message))
  })
}

results_cor_list <- vector("list", length = 9)
for (i in 1:3) {
  for (j in 1:3) {
    subset_df <- alldata[alldata$Order.x == i & alldata$Order.y == j, ]
    results_cor_list[[(i - 1) * 3 + j]] <- cor(subset_df$AFEPSI, subset_df$ALEPSI)
  }
}

list <-results_cor_list
FIGURE1EBdata <- unlist(list)
setwd("/projectnb/encore/carroll8/CODE_FOR_PITA/FINAL/DATA/")
save(FIGURE1EBdata, file = "FIGURE1EB_Data.RData")

PSIbyGO_exactly3_genes <- function(a,b){
  a <- a[a$AFEPSI >= .05,]
  b <- b[b$ALEPSI >= .05,]
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  a_3 <- a_table[a_table$Freq ==3,]
  b_3 <- b_table[b_table$Freq ==3,]
  genes <- merge(a_3, b_3, by = "Var1", all = FALSE)
  totalgenes <- nrow(genes)
  print(totalgenes)
}

totalgenes <- data.frame()
# Loop through each pair of dataframes in the list
for (i in 1:(length(AFEPSI_list))) {
  # Extract the current and next dataframes
  data <- PSIbyGO_exactly3_genes(AFEPSI_list[[i]], ALEPSI_list[[i]])
  totalgenes <- rbind(totalgenes, data)
  #print(start_id)
  #print(data)
}
sum(totalgenes)

### ext fig 1a ###

geneswithAFEs <- function(a, b) {
  a <- a[a$AFEPSI >= 0.05, ]
  atable <- as.data.frame(table(a$gene))
  alternative_first <- atable[atable$Freq >= 2, ]
  b <- b[b$ALEPSI >= 0.05, ]
  btable <- as.data.frame(table(b$gene))
  alternative_last <- btable[btable$Freq >= 2, ]
  both <- merge(alternative_last, alternative_first, by = "Var1", all = FALSE)
  altboth <- both[both$Freq.x >= 2 & both$Freq.y >= 2,]
  allgenes <- merge(atable, btable, by = "Var1", all = TRUE)
  w <- nrow(allgenes)
  x <- nrow(altboth)
  y <- nrow(alternative_first)
  z <- nrow(alternative_last)
  # Create a dataframe to store the values w, x, y, z
  result_df <- data.frame(w = w, x = x, y = y, z = z)
  return(result_df)
}
#
master_df <- data.frame(w = numeric(), x = numeric(), y = numeric(), z = numeric())
# Loop through the data frames
for (i in 1:length(AFEPSI_list)) {
  data <- geneswithAFEs(AFEPSI_list[[i]], ALEPSI_list[[i]])
  # Add the result to the master dataframe
  master_df <- rbind(master_df, data)
}
master_df
means <- colMeans(master_df)
means


#numbers used in geometric test below 
#gtex smaples- using readnum10 and conservative parameters, rounded to whole numbers
##total: 14298
##overlap:259 
##first = 1447
##set = 852
##RESULTS
## Rep factor: 3.0 (rep > 1 indicates more overlap than expected)
## p < 1.363e^-65 
  
  
#### ext fig 1b ##


maxPSI <- function(a,b){
  a <- a[a$AFEPSI >=.05,]
  amax <- a %>% group_by(gene) %>% slice(which.max(AFEPSI)) %>% ungroup()
  b <- b[b$ALEPSI >=.05,]
  bmax <- b %>% group_by(gene) %>% slice(which.max(ALEPSI)) %>% ungroup()
  first <- amax[,c(1,10)]
  last <- bmax[,c(1,10)]
  merge <- merge(first, last, by = "gene", all = FALSE)
  # alternativesonly <- merge[merge$AFEPSI !=1 & merge$ALEPSI !=1,]
  return(merge)
}
master_df <- NULL
#for maxPSI, we use a single sample per tissue as representative 
half_AFEPSI_list <- AFEPSI_list[c(1,3,5,7,9,11,13,15,17,19,21)]
half_ALEPSI_list <- ALEPSI_list[c(1,3,5,7,9,11,13,15,17,19,21)]
# Iterate through the lists of dataframes
for (i in seq_along(half_AFEPSI_list)) {
  # Apply your maxPSI function to generate a processed dataframe
  merge_result <- maxPSI(half_AFEPSI_list[[i]], half_ALEPSI_list[[i]])
  
  # Merge the resulting dataframe with the master dataframe
  if (is.null(master_df)) {
    master_df <- merge_result
  } else {
    # Merge the new dataframe with the master dataframe by the "gene" column
    master_df <- merge(master_df, merge_result, by = "gene", all = TRUE)
  }
}

master_df


###### RNA pol II mutant Data ##### 
setwd("/projectnb/encore/carroll8/PITA/bentley_data/hit_index_bam_STAR99/")
bentley_AFEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/bentley_data/hit_index_bam_STAR99/", pattern = "conservative.AFEPSI", recursive = FALSE)
bentley_ALEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/bentley_data/hit_index_bam_STAR99/", pattern = "conservative.ALEPSI", recursive = FALSE)

#load them in alphabetically 
bentley_AFEPSI_list <- lapply(bentley_AFEPSI, read.delim, sep = "\t")
bentley_ALEPSI_list <- lapply(bentley_ALEPSI, read.delim, sep = "\t")
bentley_fast_AFEPSI_list <- bentley_AFEPSI_list[c(1,2)]
bentley_slow_AFEPSI_list <- bentley_AFEPSI_list[c(3,4)]
bentley_wt_AFEPSI_list <- bentley_AFEPSI_list[c(5,6)]
bentley_fast_ALEPSI_list <- bentley_ALEPSI_list[c(1,2)]
bentley_slow_ALEPSI_list <- bentley_ALEPSI_list[c(3,4)]
bentley_wt_ALEPSI_list <- bentley_ALEPSI_list[c(5,6)]


##heatmaps were generated with the same code as in FIGURE1EB  PSIbyGO_exactly3 function above

#################### mouse data ####################

setwd("/projectnb/encore/carroll8/PITA/mouse/correctly_aligned_bams/hit_index/HITindex_aspublished_readnum2_conservative_parameters/")
mouse_AFEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/mouse/correctly_aligned_bams/hit_index/HITindex_aspublished_readnum2_conservative_parameters/", pattern = "AFEPSI", recursive = FALSE)
mouse_ALEPSI <- list.files(path = "/projectnb/encore/carroll8/PITA/mouse/correctly_aligned_bams/hit_index/HITindex_aspublished_readnum2_conservative_parameters/", pattern = "ALEPSI", recursive = FALSE)

mouse_AFEPSI_list <- lapply(mouse_AFEPSI, read.delim, sep = "\t")
mouse_ALEPSI_list <- lapply(mouse_ALEPSI, read.delim, sep = "\t")
mouse_slow_AFEPSI_list <- mouse_AFEPSI_list[c(1,2,3)]
mouse_slow_ALEPSI_list <- mouse_ALEPSI_list[c(1,2,3)]
mouse_wt_AFEPSI_list <- mouse_AFEPSI_list[c(4,5,6)]
mouse_wt_ALEPSI_list <- mouse_ALEPSI_list[c(4,5,6)]

##heatmaps were generated with the same code as in FIGURE1EB PSIbyGO_exactly3 function above 


