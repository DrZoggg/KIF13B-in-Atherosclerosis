# Load required libraries
library(promor)
library(tidyverse)
library(limma)

# Set working directory
setwd('/data/lyr/data/cooperative_project/AS_project/BulkandProtein/')

# Read protein expression matrix and metadata
protein_mtx <- read_csv('PXD062283_mtx.csv')
protein_metadata <- read_table('metadata_PXD062283.tsv')

# Preprocess metadata
protein_metadata$Sample_ID <- paste('S', protein_metadata$Sample_ID, sep = '') # Rename SampleID to match data matrix
protein_metadata <- protein_metadata[protein_metadata$Group != -1, ] # remove -1 group （Means Data Loss in original dataset）
protein_metadata <- protein_metadata[
  protein_metadata$Sample_ID %in% intersect(colnames(protein_mtx), protein_metadata$Sample_ID), 
] # Remove loss data in original data matrix 
# remaining 86 data in total

# Preprocess protein matrix, choose Protein Names and expression value columns
protein_mtx <- cbind(protein_mtx$PG.Genes, 
                     protein_mtx[, intersect(colnames(protein_mtx), protein_metadata$Sample_ID)])

# Aggregate protein data by gene (keep maximum value)
protein_mtx <- aggregate(protein_mtx, 
                         by = list(protein_mtx$`protein_mtx$PG.Genes`), 
                         FUN = max)
# Set row names and select expression columns
rownames(protein_mtx) <- protein_mtx$Group.1
protein_mtx <- protein_mtx[, 3:ncol(protein_mtx)]
# Now protein_mtx was processed to a matrix of [Protein-Gene,Sample]

# Create stage labels based on Group values
protein_metadata$Stage <- ifelse(protein_metadata$Group == 0, "Stage4",
                                 ifelse(protein_metadata$Group == 1, "Stage5",
                                        ifelse(protein_metadata$Group == 2, "Stage6", NA)))

# Create unique sample names for each stage
protein_metadata <- protein_metadata %>%
  group_by(Stage) %>%
  mutate(
    sample_num = row_number(),
    new_name = paste0(Stage, "_", sample_num)
  ) %>%
  ungroup()

# Arrange metadata by stage and sample number
protein_metadata <- arrange(protein_metadata, Stage, sample_num)

# Reorder protein matrix columns to match metadata order
protein_mtx <- protein_mtx[, protein_metadata$Sample_ID]

# Rename columns with new sample names
colnames(protein_mtx) <- protein_metadata$new_name

# Save raw matrix for reference
raw_mtx <- protein_mtx

# Filter proteins with too many missing values （over 90% NAs）
raw_filtered <- filterbygroup_na(raw_mtx, set_na = 0.9)

# Impute missing values using minimum detection method（minDet）
imp_df_mp <- impute_na(raw_filtered, seed = 42, method = 'minDet')

# Log2 transformation (add 1 to avoid log(0))
imp_df_mp <- log2(imp_df_mp + 1)

# Normalize the data （limma package quantile normalization）
norm_df <- normalize_data(imp_df_mp)

# Prepare stage factor for differential expression analysis
protein_metadata$Stage <- factor(protein_metadata$Stage, 
                                 levels = c('Stage4', 'Stage5', 'Stage6'))
boxplot_df <- data.frame(
  'KIF13B' = norm_df['KIF13B', ],
  'group' = protein_metadata$Stage
)
# Here we used Stage4 and Stage6 for further t-test
t_test_result = t.test(boxplot_df$KIF13B[boxplot_df$group=='Stage4'],boxplot_df$KIF13B[boxplot_df$group=='Stage6'],var.equal = T)
t_test_result
#	Two Sample t-test

#data:  boxplot_df$KIF13B[boxplot_df$group == "Stage4"] and boxplot_df$KIF13B[boxplot_df$group == "Stage6"]
#t = 2.3273, df = 55, p-value = 0.02366
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.4661991 6.2473176
#sample estimates:
#  mean of x mean of y 
#14.15925  10.80249 
# Save boxplot data and normalized matrix
write.csv(boxplot_df, 'Protein_KIF13B.csv')
write.csv(norm_df, 'PXD062283_normed_mtx.csv')
# 
