# Load required libraries
library(dplyr)
library(tidyr)
library(DESeq2)
library(edgeR)

# Set working directory to the location of raw fastq count files
setwd('/data/lyr/data/cooperative_project/AS_project/BulkandProtein/rawdata/fastq_raw/')

# Read and merge count files
count_files <- list.files(path = './count/', pattern = '.count', full.names = TRUE)
expression_matrix <- read.table(count_files[1])

# Merge all count files into a single expression matrix
for (file in count_files[-1]) {
  data <- read.table(file)
  expression_matrix <- merge(expression_matrix, data, by = 'V1')
}

# Clean up column names
file_names <- list.files(path = './count/', pattern = '.count', full.names = FALSE)
colnames(expression_matrix) <- c('ENSG', gsub('.count', '', file_names))

# Remove the first 5 rows (likely non-gene features)
expression_matrix <- expression_matrix[-c(1:5), ]

# Load annotation table to map Ensembl IDs to gene symbols
anno_table <- read.table('/data/lyr/databases/GRCh38.table')
rownames(anno_table) <- anno_table[, 1]

# Map Ensembl IDs to gene symbols
expression_matrix$ENSG <- anno_table[expression_matrix$ENSG, ]$V2

# Aggregate expression values by gene symbol (take mean of duplicate genes)
expression_matrix <- aggregate(expression_matrix, by = list(expression_matrix$ENSG), FUN = mean)

# Set row names to gene symbols
rownames(expression_matrix) <- expression_matrix$Group.1

# Load sample metadata
mapping_table <- read.table('../../meta_table', sep = '\t')
rownames(mapping_table) <- mapping_table$V1
colnames(mapping_table) <- c('sample', 'group')

# Prepare coldata for DESeq2
coldata <- mapping_table
coldata$group <- factor(coldata$group, levels = c('stable', 'unstable'))

# Remove unnecessary columns from expression matrix
expression_matrix <- expression_matrix[, -c(1, 2)]

# Filter lowly expressed genes using edgeR's cpm function
keep <- rowMeans(cpm(expression_matrix)) > 1
expression_matrix_filtered <- expression_matrix[keep, ]

# Convert counts to integer type for DESeq2
expr <- expression_matrix_filtered
gene_names <- rownames(expr)
expr <- apply(expr, 2, as.integer)
rownames(expr) <- gene_names

# Create DESeq2 dataset and run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = expr, colData = coldata, design = ~ group)
dds <- DESeq(dds)

# Perform variance stabilizing transformation
vst_counts <- vst(dds)
vst_data <- vst_counts@assays@data@listData[[1]]

# Extract KIF13B expression data for boxplot and statistical testing
boxplot_data_KIF13B <- data.frame(
  'KIF13B' = vst_data['KIF13B', ],
  'group' = coldata$group
)

# Perform t-test to compare KIF13B expression between groups
t_test_result <- t.test(
  boxplot_data_KIF13B$KIF13B[boxplot_data_KIF13B$group == 'stable'],
  boxplot_data_KIF13B$KIF13B[boxplot_data_KIF13B$group == 'unstable'],
  var.equal = T
)

# Print t-test results
print(t_test_result)
#Two Sample t-test

#data:  boxplot_data_KIF13B$KIF13B[boxplot_data_KIF13B$group == "stable"] and boxplot_data_KIF13B$KIF13B[boxplot_data_KIF13B$group == "unstable"]
#t = 2.6513, df = 6, p-value = 0.03796
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.01892059 0.47194156
#sample estimates:
#  mean of x mean of y 
#11.00533  10.75990 

# Save KIF13B expression data for visualization
write.csv(boxplot_data_KIF13B, 'KIF13B_vst_counts.csv')