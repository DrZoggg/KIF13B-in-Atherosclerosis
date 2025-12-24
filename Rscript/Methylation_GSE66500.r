# Load required libraries
library(GEOquery)
library(minfi)
library(limma)
library(sva)
library(ggplot2)
library(ggpubr)
library(tidyr)

# Set working directory
setwd('/data/lyr/data/cooperative_project/AS_project/methylation_merge/')

# Download and process GEO dataset GSE66500
gset <- getGEO("GSE66500", GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]

# Remove rows with missing values
gset <- gset[complete.cases(exprs(gset)), ]

# Define sample groups
group_codes <- "00000000000000000001111111111111111111"
sample_groups <- ifelse(strsplit(group_codes, split = "")[[1]] == 0, "Asymptomatic", "Symptomatic")

# Extract expression matrix
expression_matrix <- exprs(gset)

# Create metadata dataframe
metadata <- data.frame(
  sample = colnames(expression_matrix),
  group = sample_groups
)

# Create design matrix for linear model
design <- model.matrix(~ 0 + group, data = metadata)

# Fit linear model for differential methylation analysis
fit <- lmFit(expression_matrix, design)

# Define contrast (Symptomatic vs Asymptomatic)
contrast_text <- 'groupSymptomatic - groupAsymptomatic'
contrast_matrix <- makeContrasts(contrasts = contrast_text, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, 0.01)

# Extract top differentially methylated CpG sites
differential_methylation_results <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)

# Filter out lowly expressed CpG sites
differential_methylation_results <- differential_methylation_results[differential_methylation_results$AveExpr > 0.05, ]

# Load annotation for Illumina 450k methylation array
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_df <- data.frame(annotation_450k)

# Extract KIF13B CpG probes
kif13b_cpg_probes <- intersect(rownames(differential_methylation_results),
                               annotation_df$Name[annotation_df$UCSC_RefGene_Name == 'KIF13B'])

# Get results for KIF13B probes
kif13b_results <- differential_methylation_results[kif13b_cpg_probes, ]

# Extract beta values for KIF13B probes
kif13b_beta_values <- expression_matrix[kif13b_cpg_probes, ]

# Load libraries for visualization
library(Gviz)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Get annotation information for KIF13B probes
annotation_kif13b <- annotation_df[annotation_df$UCSC_RefGene_Name == 'KIF13B', ]
annotation_kif13b_subset <- annotation_kif13b[kif13b_cpg_probes, ]

# Reorder beta values according to metadata order
reordered_kif13b_beta <- kif13b_beta_values[annotation_kif13b_subset$Name, metadata$sample]

# Rename columns for better visualization
colnames(reordered_kif13b_beta) <- c(paste(rep('Asymptomatic', 19), seq(1, 19), sep = '_'),
                                     paste(rep('Symptomatic', 19), seq(1, 19), sep = '_'))

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
  col = c('#b2182c', '#054d7b')
)

# Load gene model for KIF13B
library(rtracklayer)
library(dplyr)

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
  groups = c(rep('WT', 19), rep('AS', 19)),
  type = c('smooth'),
  legend = TRUE,
  aggregateGroups = FALSE,
  aggregation = "median",
  title.width = 0.35,
  extend.left = 0.001,
  extend.right = 10000,
  from = 28880000
)

# Save KIF13B beta values
write.csv(kif13b_beta_values, './raw_data/GSE66500_Kif13b_probes.csv', quote = FALSE)

# Load previously saved methylation data (if needed)
# load('methylation_C.Rdata')

# Calculate median beta values for KIF13B across samples
kif13b_median_beta <- apply(kif13b_beta_values, 2, median)
kif13b_median_beta <- kif13b_median_beta[metadata$sample]

# Create boxplot data
boxplot_data <- cbind(metadata, kif13b_median_beta)
colnames(boxplot_data)[3] <- 'Kif13b_median_beta'
boxplot_data$group <- factor(boxplot_data$group, levels = c("Asymptomatic", "Symptomatic"))

# Perform Wilcoxon test
asymptomatic_data <- boxplot_data$Kif13b_median_beta[boxplot_data$group == "Asymptomatic"]
symptomatic_data <- boxplot_data$Kif13b_median_beta[boxplot_data$group == "Symptomatic"]
wilcox_test <- wilcox.test(asymptomatic_data, symptomatic_data, exact = FALSE, 
                           alternative = 'less', correct = TRUE)

# Create boxplot
ggplot(boxplot_data, aes(x = group, y = Kif13b_median_beta, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Asymptomatic" = "#054d7b", "Symptomatic" = "#b2182c")) +
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
    comparisons = list(c("Asymptomatic", "Symptomatic")),
    annotations = ifelse(wilcox_test$p.value < 0.001,
                         "p < 0.001",
                         sprintf("p = %.3f", wilcox_test$p.value)),
    y_position = max(boxplot_data$Kif13b_median_beta) * 1.08,
    tip_length = 0.01,
    textsize = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Save median beta values
write.csv(kif13b_median_beta, './raw_data/GSE66500_Kif13b_median_probes.csv', quote = FALSE)

# Load libraries for heatmap visualization
library(ComplexHeatmap)
library(circlize)

# Define group colors
group_colors <- c("Asymptomatic" = "#054d7b", "Symptomatic" = "#b2182c")

# Prepare metadata for heatmap
metadata_heatmap <- arrange(metadata, group)

# Create heatmap annotation
heatmap_annotation <- HeatmapAnnotation(
  Group = metadata_heatmap$group,
  col = list(Group = group_colors)
)

# Reorder beta values for heatmap
kif13b_beta_heatmap <- kif13b_beta_values[, metadata_heatmap$sample]

# Calculate Asymptomatic mean for centering
asymptomatic_samples <- metadata_heatmap$sample[metadata_heatmap$group == 'Asymptomatic']
asymptomatic_mean_beta <- rowMeans(kif13b_beta_heatmap[, asymptomatic_samples], na.rm = TRUE)

# Center beta values relative to Asymptomatic mean
kif13b_centered_beta <- kif13b_beta_heatmap - asymptomatic_mean_beta

# Define color palette for heatmap
color_function <- colorRamp2(c(-0.1, 0, 0.1), c("#054d7b", "white", "#b2182c"))

# Create heatmap
Heatmap(kif13b_centered_beta,
        name = "Beta value\n(relative to Asymptomatic)",
        col = color_function,
        top_annotation = heatmap_annotation,
        show_row_names = TRUE,
        row_title = "KIF13B Region Probes",
        column_title = "Methylation Level of KIF13B",
        cluster_columns = FALSE,
        show_column_names = FALSE)