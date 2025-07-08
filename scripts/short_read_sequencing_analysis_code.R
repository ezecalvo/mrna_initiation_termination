library(tidyverse)

# ----------------------- SHORT READ ANALYSIS CODE ----------------------- ## 

#Load in all the AFEPSI and ALEPSI outputs from hit-index pipeline for the gTEX samples as two lists of dataframes 
AFEPSI_list <- readRDS("/path/to/list/of/all/AFEPSI_files")
ALEPSI_list <- readRDS("/path/to/list/of/all/ALEPSI_files")

#load in list of metadata used in analyses- file paths, tissue names, sample names, and file types for ease of tracking specific samples
AFEPSI_meta <- readRDS("/path/to/list/of/all/AFEPSI_meta_data")
ALEPSI_meta <- readRDS("/path/to/list/of/all/ALEPSI_meta_data")


## ------ total number of annotated first and last exons per gene (1b) ------ ## 

#extract all first exons observations across all samples with gene-exon information
all_AFE_files <- list.files(path = gtex_folder, pattern = "\\.AFEPSI$", recursive = TRUE, full.names = TRUE) #load in all AFEPSI files 
all_first_exons <- data.frame(gene = character(),exon = character(),stringsAsFactors = FALSE)
for (AFEPSI_file in all_AFE_files) {
  df <- tryCatch(read.csv(AFEPSI_file, header = TRUE), error = function(e) NULL)
  if (is.null(df)) next
  if (!"exon" %in% colnames(df)) {
    message(paste("Skipping", AFEPSI_file, "- no 'exon' column"))
    next
  }
  all_first_exons <- unique(df[, c("gene", "exon")])
  all_first_exons <- unique(rbind(all_first_exons, all_first_exons))
}
#extract all last exon observations across all samples with gene-exon information
all_ALE_files <- list.files(path = gtex_folder, pattern = "\\.ALEPSI$", recursive = TRUE, full.names = TRUE) #load in all ALEPSI files 
all_last_exons <- data.frame(gene = character(),exon = character(),stringsAsFactors = FALSE)
for (ALEPSI_file in all_ALE_files) {
  df <- tryCatch(read.csv(ALEPSI_file, header = TRUE), error = function(e) NULL)
  if (is.null(df)) next
  if (!"exon" %in% colnames(df)) {
    message(paste("Skipping", ALEPSI_file, "- no 'exon' column"))
    next
  }
  all_last_exons <- unique(df[, c("gene", "exon")])
  all_last_exons <- unique(rbind(all_last_exons, all_last_exons))
}
# number of genes with x first exons (1-7, 8 or more)
all_first_exons_FEpergene <- as.data.frame(table(all_first_exons$gene))
all_first_exons_genes_FEpergene <- as.data.frame(table(all_first_exons_FEpergene$Freq))
all_first_exons_genes_FEpergene_upto7 <- all_first_exons_genes_FEpergene[c(1:7),]
all_first_exons_genes_FEpergene_8andup <- all_first_exons_genes_FEpergene[7:nrow(all_first_exons_genes_FEpergene),]
FEpergene_8andup_sum <- sum(all_first_exons_genes_FEpergene_8andup$Freq)
FErow8andup <- c(8, FEpergene_8andup_sum)
all_first_exons_genes_data <- rbind(all_first_exons_genes_FEpergene_upto7, FErow8andup)
# number of genes with x last exons (1-7, 8 or more)
all_last_exons_LEpergene <- as.data.frame(table(all_last_exons$gene))
all_last_exons_genes_LEpergene <- as.data.frame(table(all_last_exons_LEpergene$Freq))
all_last_exons_genes_LEpergene_upto7 <- all_last_exons_genes_LEpergene[c(1:7),]
all_last_exons_genes_LEpergene_8andup <- all_last_exons_genes_LEpergene[7:nrow(all_last_exons_genes_LEpergene),]
LEpergene_8andup_sum <- sum(all_last_exons_genes_LEpergene_8andup$Freq)
LErow8andup <- c(8, LEpergene_8andup_sum)
all_last_exons_genes_data <- rbind(all_last_exons_genes_LEpergene_upto7, LErow8andup) 
## number of genes with x FE and x LE dataframes
all_first_exons_genes_data$type <- rep("FE")
all_last_exons_genes_data$type <- rep("LE")
all_exons_expressed <- rbind(all_first_exons_genes_data, all_last_exons_genes_data)
figure1bdata <- all_exons_expressed


## ------ Distribution of pearson's r valuess of # FE~LEs across all samples (figure 1c) ------ ##

#FUNCTION AFEvsALE: outputs # of FE and LEs per gene, # of observatoins, Pearson's R, and p value #
AFEvsALE <- function(a, b) {
  # Remove exons with PSI < 0.05
  a <- a[a$AFEPSI >= 0.05, ] #threshold for expression
  b <- b[b$ALEPSI >= 0.05, ] #threshold for expression
  # Count FEs and LEs per gene
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  # Merge on common genes
  AFEALEgenes <- merge(a_table, b_table, by = "Var1", all = FALSE)
  colnames(AFEALEgenes) <- c("gene", "Freq.x", "Freq.y")
  # Filter for genes with at least 2 FEs and LEs
  AFEALEgenes <- AFEALEgenes[AFEALEgenes$Freq.x >= 2 & AFEALEgenes$Freq.y >= 2, ]
  # Cap Freq.x and Freq.y at 5
  AFEALEgenes$Freq.x <- pmin(AFEALEgenes$Freq.x, 5)
  AFEALEgenes$Freq.y <- pmin(AFEALEgenes$Freq.y, 5)
  # Count gene combinations
  AFEALEgenes_count <- AFEALEgenes %>% group_by(Freq.x, Freq.y) %>% tally()
  # Pearson correlation + p-values
  correlation_test <- cor.test(AFEALEgenes$Freq.x, AFEALEgenes$Freq.y, method = "pearson")
  rounded_cor <- round(correlation_test$estimate, 3)
  p_value <- correlation_test$p.value
  print(paste("r =", rounded_cor))
  print(paste("p =", signif(p_value, 3)))
  return(list(AFEALEgenes, AFEALEgenes_count, rounded_cor, p_value))
}
# apply function to all AFEPSI and ALEPSI samples in list
allsamples_number_AFE_vs_ALE <- as.data.frame(mapply(AFEvsALE, AFEPSI_list, ALEPSI_list)) 
cor_values <- sapply(allsamples_number_AFE_vs_ALE, function(x) x[[3]][1]) #extract pearson's (3rd element)
cor_vec <- as.numeric(cor_values)
cor_clean <- cor_values[is.finite(cor_values)] #remove non finite values 
figure1cdata <- cor_clean #data to be plotted in density plot 
figure1c_r_mean <- mean(cor_clean, na.rm = TRUE) #extract mean cor value


## ------ Number of first and last exons per gene including single FE/LE (figure 1c inset) ------ ##

AFEvsALE_including1 <- function(a, b) {
  # Remove exons with PSI < 0.05
  a <- a[a$AFEPSI >= 0.05, ]
  b <- b[b$ALEPSI >= 0.05, ]
  # Count FEs and LEs per gene
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  # Merge on common genes
  AFEALEgenes <- merge(a_table, b_table, by = "Var1", all = FALSE)
  colnames(AFEALEgenes) <- c("gene", "Freq.x", "Freq.y")
  # Filter for genes with at least 2 FEs and LEs
  AFEALEgenes <- AFEALEgenes[AFEALEgenes$Freq.x >= 1,]
  # Cap Freq.x and Freq.y at 5
  AFEALEgenes$Freq.x <- pmin(AFEALEgenes$Freq.x, 5)
  AFEALEgenes$Freq.y <- pmin(AFEALEgenes$Freq.y, 5)
  # Count gene combinations
  AFEALEgenes_count <- AFEALEgenes %>% group_by(Freq.x, Freq.y) %>% tally()
  # Pearson correlation + p-values
  correlation_test <- cor.test(AFEALEgenes$Freq.x, AFEALEgenes$Freq.y, method = "pearson")
  rounded_cor <- round(correlation_test$estimate, 3) 
  p_value <- correlation_test$p.value
  print(paste("r =", rounded_cor))
  print(paste("p =", signif(p_value, 3)))
  return(list(AFEALEgenes, AFEALEgenes_count, rounded_cor, p_value))
}

# reviewer asked to plot ALL genes including those with 1 FE or LE- adjust the AFEvsALE code for this output: 
plotting_afevsale_1ormore <- as.data.frame(mapply(AFEvsALE_including1, AFEPSI_list, ALEPSI_list)) 
# use the correlations from figure1cdata to show the correlation between number of FEs and LEs in multi FE/LE genes.
#Below we need to extract the data, correlation value, and p value of one representative samples to be plotted
fig1c_inset_correlation <- figure1cdata[[12363]][[3]]
fig1c_inset_pval <- figure1cdata[[12363]][[4]]
fig1c_inset_data <- figure1cdata[[12363]][[2]]


## ------ Correlation between AFEPSI~ALEPSI based on ordinal position in all genes with 2+ AFE/ALEs (figure1da) ------  ##

#FUNCTION PSIbyGO_2to5_only takes AFEPSI and ALEPSI outputs to compare the AFEPSI~ALEPSI per gene by genomic order#

PSIbyGO_2to5_only <- function(a, b) {
  # Filter for usable PSI values
  a <- a[a$AFEPSI >= 0.05, ]
  b <- b[b$ALEPSI >= 0.05, ]
  # Keep genes with 2 to 5 AFEs and ALEs
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  a_filtered_genes <- a_table[a_table$Freq >= 2 & a_table$Freq <= 5, ]
  b_filtered_genes <- b_table[b_table$Freq >= 2 & b_table$Freq <= 5, ]
  genes <- merge(a_filtered_genes, b_filtered_genes, by = "Var1")
  colnames(genes) <- c("gene", "AFE_count", "ALE_count")
  # Filter AFE/ALE data to include only those genes
  a <- a %>% filter(gene %in% genes$gene)
  b <- b %>% filter(gene %in% genes$gene)
  # Split exon string into genomic coordinates
  a <- separate(a, exon, into = c("chromosome", "start", "end"))
  b <- separate(b, exon, into = c("chromosome", "start", "end"))
  a$start <- as.numeric(a$start)
  b$start <- as.numeric(b$start)
  # Strand-specific ordering
  a_plus  <- a %>% filter(strand == "+") %>% arrange(gene, start)
  a_minus <- a %>% filter(strand == "-") %>% arrange(gene, desc(start))
  b_plus  <- b %>% filter(strand == "+") %>% arrange(gene, start)
  b_minus <- b %>% filter(strand == "-") %>% arrange(gene, desc(start))
  # Assign genomic order per gene
  assign_order <- function(df) {
    df$Order <- ave(seq_len(nrow(df)), df$gene, FUN = seq_along)
    return(df)
  }
  a_ordered <- bind_rows(assign_order(a_plus), assign_order(a_minus))
  b_ordered <- bind_rows(assign_order(b_plus), assign_order(b_minus))
  # Split by order index (1–5)
  go_AFE_list <- split(a_ordered, a_ordered$Order)
  go_ALE_list <- split(b_ordered, b_ordered$Order)
  # Combine all pairwise combinations
  result_list <- list()
  for (i in 1:5) {
    for (j in 1:5) {
      if (!is.null(go_AFE_list[[i]]) && !is.null(go_ALE_list[[j]])) {
        merged_df <- merge(go_AFE_list[[i]], go_ALE_list[[j]], by = "gene")
        if (nrow(merged_df) > 0) {
          merged_df <- merged_df %>%
            left_join(genes, by = "gene")
          result_list[[paste0("x_", i, "v", j)]] <- merged_df
        }
      }
    }
  }
  combined_df <- bind_rows(result_list)
  return(combined_df)
}
# loop through all pairs of AFEPSI/ALEPSI and apply to all gtex samples #
alldata <- data.frame()
start_id <- 0
for (i in seq_along(AFEPSI_list)) {
  tryCatch({
    message(sprintf("Processing sample %d of %d", i, length(AFEPSI_list)))
    data <- PSIbyGO_2to5_only(AFEPSI_list[[i]], ALEPSI_list[[i]])
    if (nrow(data) > 0) {
      # Assign unique gene IDs per sample to avoid overwriting
      start_id <- start_id + length(unique(data$gene))
      data$gene <- as.numeric(factor(data$gene, levels = unique(data$gene))) + start_id
      alldata <- bind_rows(alldata, data)
    } else {
      message(sprintf("No data returned for sample %d", i))
    }
  }, error = function(e) {
    message(sprintf("Error in sample %d: %s", i, e$message))
  })
}
# compute correlations
results_cor_list <- vector("list", length = 25)
for (i in 1:5) {
  for (j in 1:5) {
    index <- (i - 1) * 5 + j
    message(sprintf("Correlating AFE order %d vs ALE order %d (%d/25)", i, j, index))
    subset_df <- alldata[alldata$Order.x == i & alldata$Order.y == j, ]
    if (nrow(subset_df) >= 2) {
      results_cor_list[[index]] <- cor(subset_df$AFEPSI, subset_df$ALEPSI)
    } else {
      results_cor_list[[index]] <- NA
    }
  }
}
figure1da_data <- results_cor_list

## ------- Correlation between the AFEPSI~ALEPSI in genes with EXACTLY 3 FE and 3 LE  based on ordinal position (figure1db) ------- ##

# FunctionPSIbyGO_exactly 3 
PSIbyGO_exactly3 <- function(a, b) {
  a <- a[a$AFEPSI >= 0.05, ]
  b <- b[b$ALEPSI >= 0.05, ]
  a_table <- as.data.frame(table(a$gene))
  b_table <- as.data.frame(table(b$gene))
  a_3 <- a_table[a_table$Freq == 3, ]
  b_3 <- b_table[b_table$Freq == 3, ]
  genes <- merge(a_3, b_3, by = "Var1")
  a_3 <- a %>% filter(gene %in% genes$Var1)
  b_3 <- b %>% filter(gene %in% genes$Var1)
  a_3 <- separate(a_3, exon, into = c("chromosome", "start", "end"))
  b_3 <- separate(b_3, exon, into = c("chromosome", "start", "end"))
  a_3$start <- as.numeric(a_3$start)
  b_3$start <- as.numeric(b_3$start)
  a_3_plus  <- a_3 %>% filter(strand == "+") %>% arrange(start)
  a_3_minus <- a_3 %>% filter(strand == "-") %>% arrange(desc(start))
  b_3_plus  <- b_3 %>% filter(strand == "+") %>% arrange(start)
  b_3_minus <- b_3 %>% filter(strand == "-") %>% arrange(desc(start))
  a_3_plus$Order  <- ave(seq_len(nrow(a_3_plus)), a_3_plus$gene, FUN = seq_along)
  a_3_minus$Order <- ave(seq_len(nrow(a_3_minus)), a_3_minus$gene, FUN = seq_along)
  b_3_plus$Order  <- ave(seq_len(nrow(b_3_plus)), b_3_plus$gene, FUN = seq_along)
  b_3_minus$Order <- ave(seq_len(nrow(b_3_minus)), b_3_minus$gene, FUN = seq_along)
  AFEs <- bind_rows(a_3_plus, a_3_minus)
  ALEs <- bind_rows(b_3_plus, b_3_minus)
  go_AFE_list <- split(AFEs, AFEs$Order)
  go_ALE_list <- split(ALEs, ALEs$Order)
  pairwise_list <- list()
  for (i in 1:3) {
    for (j in 1:3) {
      pairwise_list[[paste0("x_", i, "_v_", j)]] <- merge(go_AFE_list[[i]], go_ALE_list[[j]], by = "gene")
    }
  }
  combined_df <- bind_rows(pairwise_list)
  return(combined_df)
}

# loop through all pairs of AFEPSI/ALEPSI per gtex sample and apply above function to all gtex samples  #
results_list <- vector("list", length = length(AFEPSI_list))

for (i in seq_along(AFEPSI_list)) {
  tryCatch({
    results_list[[i]] <- PSIbyGO_exactly3(AFEPSI_list[[i]], ALEPSI_list[[i]])
    message(sprintf("looped through iteration %d", i))
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
  })
}
alldata <- bind_rows(results_list)
# for this analyses- we do PER GENE mean, so that not a single highly expressed gene that uses 3 and 3 is overrepresented- want all genes represented 
gene_means <- alldata %>%
  group_by(gene, Order.x, Order.y) %>%
  summarize(
    AFEPSI = mean(AFEPSI, na.rm = TRUE),
    ALEPSI = mean(ALEPSI, na.rm = TRUE),
    .groups = "drop"
  )
# correlation analysis 
results_cor_list <- vector("list", length = 9)
cor_matrix <- matrix(nrow = 3, ncol = 3)
for (i in 1:3) {
  for (j in 1:3) {
    subset_df <- gene_means %>%
      filter(Order.x == i, Order.y == j)
    cor_val <- cor(subset_df$AFEPSI, subset_df$ALEPSI, use = "pairwise.complete.obs")
    results_cor_list[[(i - 1) * 3 + j]] <- cor_val
    cor_matrix[i, j] <- cor_val
  }
}
rownames(cor_matrix) <- paste0("AFE", 1:3)
colnames(cor_matrix) <- paste0("ALE", 1:3)
figure1db_data <- as.vector(cor_matrix)


##this code was also used for pol ii mutant analyses- mutant cell line and mouse data in manuscript

## ------ Average # of genes with AFES, genes with ALEs, and genes with both (figure S1A) ----- ##

#get the number of genes with AFEs, ALEs, BOTH AFEs and ALEs, and total number of genes
geneswithAFEsandorALEs <- function(a, b) {
  a <- a[a$AFEPSI >= 0.05, ]
  atable <- as.data.frame(table(a$gene))
  alternative_first <- atable[atable$Freq >= 2, ]
  b <- b[b$ALEPSI >= 0.05, ]
  btable <- as.data.frame(table(b$gene))
  alternative_last <- btable[btable$Freq >= 2, ]
  # Find genes with both AFE and ALE events
  both <- merge(alternative_last, alternative_first, by = "Var1", all = FALSE)
  altboth <- both[both$Freq.x >= 2 & both$Freq.y >= 2, ]
  # All genes present in either table
  allgenes <- merge(atable, btable, by = "Var1", all = TRUE)
  # Extract counts
  w <- nrow(allgenes)         # total genes
  x <- nrow(altboth)          # genes with both AFE and ALE
  y <- nrow(alternative_first) # genes with AFE
  z <- nrow(alternative_last)  # genes with ALE
  return(data.frame(w = w, x = x, y = y, z = z))
}
# Initialize dataframe to store all samples' results
master_df <- data.frame(w = numeric(), x = numeric(), y = numeric(), z = numeric())
# Loop through all samples
for (i in 1:length(AFEPSI_list)) {
  result <- geneswithAFEsandorALEs(AFEPSI_list[[i]], ALEPSI_list[[i]])
  master_df <- rbind(master_df, result)
}
# Calculate mean values across all samples
means <- colMeans(master_df)
figures1adata <- means


## -------- Mean pearson's R by tissue category (figureS1B) -------- ##

#load into the data from figure1c- corr between AFEPSI~ALEPSI 
figure1cdata 
# Step 2: Extract Pearson correlations (3rd element of each list)
cor_values <- sapply(figure1cdata , `[[`, 3)
# Step 3: Add correlation values to metadata
AFEPSI_meta <- AFEPSI_meta %>%
  mutate(
    pearson_correlation = cor_values,
    # Step 4: Assign compact tissue categories- categorize tissues into groups 
    category = case_when(
      str_detect(tissue, "adipose") ~ "adipose",
      str_detect(tissue, "blood_whole|ebv") ~ "blood",
      str_detect(tissue, "blood-vessel") ~ "blood_vessel",
      str_detect(tissue, "^brain|spinal-cord") ~ "brain",
      str_detect(tissue, "breast") ~ "breast",
      str_detect(tissue, "colon|esophagus|stomach|small-intestine") ~ "digestive",
      str_detect(tissue, "adrenal|pancreas|pituitary|thyroid") ~ "endocrine",
      str_detect(tissue, "ovary|uterus|vagina|cervix|fallopian") ~ "female_reproductive",
      str_detect(tissue, "prostate|testis") ~ "male_reproductive",
      str_detect(tissue, "heart") ~ "heart",
      str_detect(tissue, "kidney") ~ "kidney",
      str_detect(tissue, "liver") ~ "liver",
      str_detect(tissue, "lung") ~ "lung",
      str_detect(tissue, "muscle") ~ "muscle",
      str_detect(tissue, "nerve") ~ "nerve",
      str_detect(tissue, "salivary") ~ "salivary_gland",
      str_detect(tissue, "skin") ~ "skin",
      str_detect(tissue, "spleen|bladder") ~ "bladder/spleen",
      TRUE ~ "other"
    )
  )
# Step 5: Summarize mean correlation per category
mean_corr_by_category <- AFEPSI_meta %>%
  group_by(category) %>%
  summarise(
    mean_correlation = mean(pearson_correlation, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_correlation))
# Step 6: Save outputs
figures1bdata <- mean_corr_by_category

## -------- Max AFEPSI and Max ALEPSI per gene (figure S1C) -------- ##

#Function maxPSI to get the max AFEPSI and ALEPSI per gene 
maxPSI <- function(a, b, sample_id = NA) {
  a <- a[a$AFEPSI >= 0.05, ]
  amax <- a %>% group_by(gene) %>% slice(which.max(AFEPSI)) %>% ungroup()
  b <- b[b$ALEPSI >= 0.05, ]
  bmax <- b %>% group_by(gene) %>% slice(which.max(ALEPSI)) %>% ungroup()
  first <- amax[, c("gene", "AFEPSI")]
  last  <- bmax[, c("gene", "ALEPSI")]
  merged <- merge(first, last, by = "gene", all = FALSE)
  merged$sample_id <- sample_id
  return(merged)
}
#do across all samples 
all_results <- list()
for (i in seq_along(AFE_list)) {
  result <- maxPSI(AFE_list[[i]], ALE_list[[i]], sample_id = i)
  all_results[[i]] <- result
}
# Combine all results into one big data frame, mean maxPSI per gene 
master_df <- bind_rows(all_results)
mean_maxPSI_per_gene <- master_df %>%
  group_by(gene) %>%
  summarise(
    mean_max_AFEPSI = mean(AFEPSI, na.rm = TRUE),
    mean_max_ALEPSI = mean(ALEPSI, na.rm = TRUE),
    n_samples = n()
  ) %>%
  filter(n_samples >= 1000)  #keep genes that are expressed in over 1000 samples
# Reshape to long format
heatmap_data <- mean_maxPSI_per_gene %>%
  select(gene, mean_max_AFEPSI, mean_max_ALEPSI) %>%
  pivot_longer(cols = c(mean_max_AFEPSI, mean_max_ALEPSI),
               names_to = "exon_type", values_to = "PSI")
# reshape for plotting
psi_matrix <- heatmap_data %>%
  pivot_wider(names_from = exon_type, values_from = PSI) %>%
  column_to_rownames("gene") %>%
  as.matrix()

figures1cdata <- psi_matrix


## ---------- promoter sequence analysis  ---------- ##

#load in the AFEPSI hitindex outputs- replicates. 
rep1 <- read.delim("/path/to/K652/or/relevant/afepsi/replicate1.AFEPSI")
rep2 <- read.delim("/path/to/K652/or/relevant/afepsi/replicate2.AFEPSI")
#assign genomic orders to the AFES (here acting as proxies for promoters)

genomic_order <- function(a) {
  a <- a[a$AFEPSI >= 0.05, ]
  exon_counts <- as.data.frame(table(a$gene))
  genes_with_2plus <- exon_counts$Var1[exon_counts$Freq >= 2]
  a_2plus <- a[a$gene %in% genes_with_2plus, ]
  a_2plus <- separate(a_2plus, exon, into = c("chromosome", "start", "end"))
  a_2plus$start <- as.numeric(a_2plus$start)
  a_plus  <- a_2plus %>% filter(strand == "+") %>% arrange(start)
  a_minus <- a_2plus %>% filter(strand == "-") %>% arrange(desc(start))
  a_plus$Order  <- ave(seq_len(nrow(a_plus)), a_plus$gene, FUN = seq_along)
  a_minus$Order <- ave(seq_len(nrow(a_minus)), a_minus$gene, FUN = seq_along)
  AFEs <- bind_rows(a_plus, a_minus)
  # Collapse any Order >= 2 to "2"- so we are comparing most upstream promoter to downstream promoters
  AFEs$Order <- ifelse(AFEs$Order >= 2, 2, 1)
  AFEs$Order <- as.character(AFEs$Order)
  return(AFEs)
}

ordered_rep1 <- genomic_order(rep1)
ordered_rep2 <- genomic_order(rep2)
#load in the classifications of K652 genes
pita_classifications <- read.csv("../pita_classif_restrictive.csv") 
ordered_rep1 <- ordered_rep1 %>%
  mutate(classification = pita_classifications$classification[match(gene, pita_classifications$gene_id)])
ordered_rep2 <- ordered_rep2 %>%
  mutate(classification = pita_classifications$classification[match(gene, pita_classifications$gene_id)])
ordered_rep1_PITA <- ordered_rep1[ordered_rep1$classification == "PITA",]
ordered_rep1_noPITA <- ordered_rep1[ordered_rep1$classification == "no PITA",]
ordered_rep2_PITA <- ordered_rep2[ordered_rep2$classification == "PITA",]
ordered_rep2_noPITA <- ordered_rep2[ordered_rep2$classification == "no PITA",]
#extract promoter sequences of PITA and non PITA genes from both replicates
pos1_vs_downstream_PITA_rep1_bed <- ordered_rep1_PITA %>%
  mutate(
    tss = start,
    promoter_start = ifelse(strand == "+", pmax(0, tss - 275), tss - 50),
    promoter_end   = ifelse(strand == "+", tss + 50, tss + 275),
    bed_start = pmin(promoter_start, promoter_end),
    bed_end   = pmax(promoter_start, promoter_end),
    name = paste0(gene, "_", Order)
  ) %>%
  transmute(
    chrom = as.character(chromosome),  # REMOVE chr prefix
    bed_start, bed_end, name, score = "-", strand
  ) %>%
  arrange(name)
write.table(pos1_vs_downstream_PITA_rep1_bed, "pos1_vs_downstream_rep1_promoter_PITA_coordinates.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#do this for all 4- both reps, both types of genes
#use bedtofasta to get fasta files
#complete iProEP analysis and CpG island analysis 


