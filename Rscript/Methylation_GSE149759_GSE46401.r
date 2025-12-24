# Load required libraries
library(GEOquery)
library(minfi)
library(limma)
library(dplyr)
library(sva)
library(ggplot2)
library(ggpubr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(Gviz)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Set working directory
setwd('/data/lyr/data/cooperative_project/AS_project/methylation_merge/')

# Download and process GEO datasets
gset2_1 <- getGEO("GSE149759", GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
gset2_2 <- getGEO("GSE149759", GSEMatrix = TRUE, AnnotGPL = FALSE)[[2]]
gset3 <- getGEO("GSE46401", GSEMatrix = TRUE, AnnotGPL = FALSE)[[2]]

# Remove rows with missing values
gset2_1 <- gset2_1[complete.cases(exprs(gset2_1)), ]
gset2_2 <- gset2_2[complete.cases(exprs(gset2_2)), ]
gset3 <- gset3[complete.cases(exprs(gset3)), ]

# Find common CpG sites across all datasets
common_cpg_ids <- intersect(intersect(gset2_1@featureData@data$ID, 
                                      gset2_2@featureData@data$ID), 
                            gset3@featureData@data$ID)

# Define sample groups for each dataset
group_codes_gse149759_1 <- "111000"
sample_groups_gse149759_1 <- ifelse(strsplit(group_codes_gse149759_1, split = "")[[1]] == 0, "AS", "WT")

group_codes_gse149759_2 <- "100000000000"
sample_groups_gse149759_2 <- ifelse(strsplit(group_codes_gse149759_2, split = "")[[1]] == 0, "AS", "WT")

group_codes_gse46401 <- "0000000000000001111111111111110000000000000000000"
sample_groups_gse46401 <- ifelse(strsplit(group_codes_gse46401, split = "")[[1]] == 0, "AS", "WT")

# Subset datasets to common CpG sites
gset2_1 <- gset2_1[common_cpg_ids, ]
gset2_2 <- gset2_2[common_cpg_ids, ]
gset3 <- gset3[common_cpg_ids, ]

# Combine expression matrices from all datasets
combined_expr <- cbind(exprs(gset2_1)[common_cpg_ids, ],
                       exprs(gset2_2)[common_cpg_ids, ],
                       exprs(gset3)[common_cpg_ids, ])

# Create metadata dataframe
metadata <- data.frame(
  sample = colnames(combined_expr),
  group = c(sample_groups_gse149759_1, 
            sample_groups_gse149759_2, 
            sample_groups_gse46401),
  batch = c(rep('GSE149759_1', length(sample_groups_gse149759_1)),
            rep('GSE149759_2', length(sample_groups_gse149759_2)),
            rep('GSE46401', length(sample_groups_gse46401)))
)

# Create design matrix for linear model
design <- model.matrix(~ 0 + group, data = metadata)

# Perform batch correction using ComBat
corrected_beta <- ComBat(
  dat = combined_expr,
  batch = metadata$batch,
  mod = model.matrix(~ metadata$group)
)

# Fit linear model for differential methylation analysis
fit <- lmFit(corrected_beta, design)

# Define contrast (AS vs WT)
contrast_text <- 'groupAS - groupWT'
contrast_matrix <- makeContrasts(contrasts = contrast_text, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, 0.01)

# Extract top differentially methylated CpG sites
differential_methylation_results <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

# Load annotation for Illumina 450k methylation array
annotation_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_df <- data.frame(annotation_450k)

# Extract KIF13B CpG probes
kif13b_cpg_probes <- intersect(rownames(differential_methylation_results),
                               rownames(annotation_df)[annotation_df$UCSC_RefGene_Name == 'KIF13B'])

# Get results for KIF13B probes
kif13b_results <- differential_methylation_results[kif13b_cpg_probes, ]

# Extract beta values for KIF13B probes
kif13b_beta_values <- corrected_beta[kif13b_cpg_probes, ]

# Sort metadata by group (WT first, then AS)
metadata <- arrange(metadata, desc(group))

# Get annotation information for KIF13B probes
annotation_kif13b <- annotation_df[annotation_df$UCSC_RefGene_Name == 'KIF13B', ]
annotation_kif13b_subset <- annotation_kif13b[kif13b_cpg_probes, ]

# Reorder beta values according to sorted metadata
reordered_kif13b_beta <- kif13b_beta_values[, metadata$sample]

# Rename columns for better visualization
colnames(reordered_kif13b_beta) <- c(paste(rep('WT', 19), seq(1, 19), sep = '_'),
                                     paste(rep('AS', 48), seq(1, 48), sep = '_'))

# Create GRanges object for methylation data visualization
methylation_granges <- GRanges(
  seqnames = annotation_kif13b_subset$chr,
  ranges = IRanges(start = annotation_kif13b_subset$pos, end = annotation_kif13b_subset$pos),
  beta = reordered_kif13b_beta
)

# Create annotation track for CpG probes
annotation_track <- AnnotationTrack(
  methylation_granges,
  name = 'Probes',
  stacking = 'hide',
  background.title = "#98FB98",
  cex.title = 0.5,
  rotation.title = 0
)

# Create genome axis track
genome_axis_track <- GenomeAxisTrack()

# Create data track for methylation beta values
methylation_data_track <- DataTrack(
  methylation_granges,
  name = "Methylation Beta",
  ylim = c(0, 1),
  showAxis = TRUE,
  background.title = "darkblue",
  cex.title = 0.8,
  col = c('#CC6602', '#099396')
)
# Process GTF file to create gene model (if not already done)
# This section is for creating the hg19_model.csv file used earlier
library(rtracklayer)

gtf_path <- "/data/pubres/genomes/GRCh37/GRCh37.annotation.gtf"
gtf <- import(gtf_path)
gtf_df <- as.data.frame(gtf)

# Extract gene information
genes <- gtf_df %>%
  filter(type == "gene") %>%
  select(gene_id, gene_name, gene_type) %>%
  distinct()

# Extract transcript information
transcripts <- gtf_df %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id, transcript_type, transcript_name) %>%
  distinct()

# Extract exon and UTR information
exons_utr <- gtf_df %>%
  filter(type %in% c("exon", "five_prime_utr", "three_prime_utr")) %>%
  mutate(
    feature = case_when(
      type == "five_prime_utr" ~ "utr5",
      type == "three_prime_utr" ~ "utr3",
      TRUE ~ gene_type
    )
  ) %>%
  select(seqnames, start, end, strand, feature, gene_id, exon_id, transcript_id)

# Merge all information
gene_model <- exons_utr %>%
  left_join(genes, by = "gene_id") %>%
  left_join(transcripts, by = c("gene_id", "transcript_id")) %>%
  mutate(
    width = end - start + 1,
    symbol = gene_name
  ) %>%
  select(
    chromosome = seqnames,
    start,
    end,
    width,
    strand,
    feature,
    gene = gene_id,
    exon = exon_id,
    transcript = transcript_id,
    symbol
  ) %>%
  arrange(start)

# Clean up IDs (remove version numbers)
gene_model$gene <- gsub("\\..*", "", gene_model$gene)
gene_model$exon <- gsub("\\..*", "", gene_model$exon)
gene_model$transcript <- gsub("\\..*", "", gene_model$transcript)

# Save gene model
write.csv(gene_model, 'hg19_model.csv', row.names = FALSE)


# Load gene model for KIF13B
gene_model <- read.csv('hg19_model.csv')
kif13b_model <- gene_model[gene_model$symbol == 'KIF13B', ]

# Define genome and chromosome for visualization
genome <- 'hg19'
chromosome <- 'chr8'

# Create gene region track for KIF13B
gene_region_track <- GeneRegionTrack(
  kif13b_model,
  genome = genome,
  chromosome = chromosome,
  name = "Gene Model",
  transcriptAnnotation = "symbol",
  background.title = '#F0E68C',
  background.panel = "#dbeeff",
  cex.title = 0.8,
  col = NULL
)

# Plot tracks
plotTracks(
  list(genome_axis_track, gene_region_track, annotation_track, methylation_data_track),
  groups = c(rep('WT', 19), rep('AS', 48)),
  type = c('smooth'),
  legend = TRUE,
  aggregateGroups = FALSE,
  aggregation = "median",
  title.width = 0.35,
  extend.left = 0.001,
  extend.right = 10000,
  from = 28870000
)

# Save KIF13B beta values
write.csv(kif13b_beta_values, './raw_data/GSE149759_GSE46401_Kif13b_probes.csv', quote = FALSE)

# Calculate median beta values for KIF13B across samples
kif13b_median_beta <- apply(kif13b_beta_values, 2, median)
kif13b_median_beta <- kif13b_median_beta[metadata$sample]

# Create boxplot data
boxplot_data <- cbind(metadata[, c('sample', 'group')], kif13b_median_beta)
colnames(boxplot_data)[3] <- 'Kif13b_median_beta'
boxplot_data$group <- factor(boxplot_data$group, levels = c("WT", "AS"))

# Perform Wilcoxon test
wilcox_test <- wilcox.test(
  boxplot_data$Kif13b_median_beta[boxplot_data$group == "WT"],
  boxplot_data$Kif13b_median_beta[boxplot_data$group == "AS"],
  exact = FALSE,
  alternative = 'less'
)

# Create boxplot
ggplot(boxplot_data, aes(x = group, y = Kif13b_median_beta, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("WT" = "#099396", "AS" = "#CC6602")) +
  labs(
    x = "Group",
    y = "Median KIF13B methylation level (beta value)",
    title = "Median Methylation Level"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  geom_signif(
    comparisons = list(c("WT", "AS")),
    annotations = ifelse(wilcox_test$p.value < 0.001,
                         "p < 0.001",
                         sprintf("p = %.3f", wilcox_test$p.value)),
    y_position = max(boxplot_data$Kif13b_median_beta) * 1.08,
    tip_length = 0.01,
    textsize = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Save median beta values
write.csv(kif13b_median_beta, './raw_data/GSE149759_GSE46401_Kif13b_median_probes.csv', quote = FALSE)

# Load libraries for heatmap visualization
library(ComplexHeatmap)
library(circlize)

# Define group colors
group_colors <- c("WT" = "#099396", "AS" = "#CC6602")

# Prepare metadata for heatmap
metadata_heatmap <- arrange(metadata, desc(group))

# Create heatmap annotation
heatmap_annotation <- HeatmapAnnotation(
  Group = metadata_heatmap$group,
  col = list(Group = group_colors)
)

# Reorder beta values for heatmap
kif13b_beta_heatmap <- kif13b_beta_values[, metadata_heatmap$sample]

# Calculate WT mean for centering
wt_samples <- metadata_heatmap$sample[metadata_heatmap$group == 'WT']
wt_mean_beta <- rowMeans(kif13b_beta_heatmap[, wt_samples], na.rm = TRUE)

# Center beta values relative to WT mean
kif13b_centered_beta <- kif13b_beta_heatmap - wt_mean_beta

# Define color palette for heatmap
color_function <- colorRamp2(c(-0.2, 0, 0.2), c("#099396", "white", "#CC6602"))

# Create heatmap
Heatmap(kif13b_centered_beta,
        name = "Beta value\n(relative to WT)",
        col = color_function,
        top_annotation = heatmap_annotation,
        show_row_names = TRUE,
        row_title = "KIF13B Region Probes",
        column_title = "Methylation Level of KIF13B",
        cluster_columns = FALSE,
        show_column_names = FALSE)

