#!/bin/bash
suppressPackageStartupMessages({
library(reticulate)
np <- import("numpy")
library(tidyverse)
library(R.utils)
library(data.table)
library(DESeq2)
library(readxl)
library(argparser, quietly=TRUE)
library("BiocParallel")
})

register(MulticoreParam(40))

today <-
  Sys.Date()
today <- format(today, format="%y%m%d")

# Create a parser
p <- arg_parser("Run R deseq with subsetting")

# Add command line arguments
p <- add_argument(p, "index_dir", help='default location to look for all files')
p <- add_argument(p, "meta_file", help="metadata spreadsheet, currently csv format")
p <- add_argument(p, "meta_column", help="column of meta file to use for comparisons")
p <- add_argument(p, "meta_comparison", help="items in meta column for comparison that will be used for this run, in a vector")
p <- add_argument(p, "n_cpus", help="number of cores to specify for biocparallel", default=40)
p <- add_argument(p, "output_dir", help='output directory for published files', default=paste0('/net/seq/data2/projects/encode4/4078_index/deseq/', today,'/'))
p <- add_argument(p, "--binary_file", help="binary matrix, .npy format", default = 'binary.only_autosomes.filtered.matrix.npy')
p <- add_argument(p, "--count_file", help="count matrix, .npy format", default = 'counts.only_autosomes.filtered.matrix.npy')
p <- add_argument(p, "--scale_factors", help="scale factors to apply to count matrix during normalization, .npy format", default = 'normalization/normalized.only_autosomes.filtered.scale_factors.mean_normalized.npy')
#p <- add_argument(p, "--n_cpus", help="number of cores to specify for biocparallel", default=40)
#p <- add_argument(p, "--output_dir", help='output directory for published files', default=paste0('/net/seq/data2/projects/encode4/4025_index/deseq/', today,'/'))

# Parse the command line arguments
argv <- parse_args(p)

print(argv$meta_column)
print(argv$meta_comparison)

# Rscript round.R 3.141 -d 2

sample_order <- read_table(paste0(argv$index_dir, 'samples_order.txt'), col_names=FALSE)
sample_order

# test_meta <- read_excel(argv$meta_file) %>% filter(ag_id %in% sample_order$X1)

test_meta <- read.csv(argv$meta_file, header = TRUE) %>% filter(ag_id %in% sample_order$X1)

meta <- test_meta[ order(match(test_meta$ag_id, sample_order$X1)), ]
# test_meta$combined_annotation <- paste(test_meta$core_annotation, test_meta$extended_annotation, sep = ' ')
# print(test_meta$combined_annotation)

# meta <- test_meta %>% mutate('stripped_cluster' = str_replace_all(combined_annotation, regex("\\W+"), ''))
# print(meta$stripped_cluster)

#meta <- test_meta %>% mutate('stripped_annotation' = str_replace_all( !! rlang::sym(argv$meta_column) , regex("\\W+"), ''))
#print(paste0('the alphanumeric only version of the meta column is: ', meta$stripped_annotation %>% unique()))

comparisons <- argv$meta_comparison %>% str_trim() %>% str_replace_all(regex("\\W+"), '') %>% str_split('_') %>% data.frame()
print(comparisons)

#print(paste0('subsetting annotations in ',str(comparisons)))
meta_sub <- meta %>% filter( !! rlang::sym(argv$meta_column) %in% comparisons[,1])

if(nrow(meta_sub) == 0L) stop('Error in metadata subset, no datasets left')

print(head(meta_sub$ag_id))
print(paste0('number of samples: ',nrow(meta_sub)))
cluster_group <- paste0(str_replace_all(comparisons[1,1], regex("\\W+"), ''), 'vs', str_replace_all(comparisons[2,1], regex("\\W+"), ''))




binary_file = paste0(argv$index_dir,argv$binary_file)
print(paste0('loading binary matrix at ',binary_file))
binary_matrix <- np$load(binary_file, mmap_mode='r')
print('converting to data frame')
binary_matrix = data.frame(binary_matrix %>% t(), row.names = sample_order$X1)
print(paste0('binary file has ',nrow(binary_matrix), ' rows and ', ncol(binary_matrix), ' columns'))
print('subsetting to sample groups')
binary_matrix <- binary_matrix %>% filter(row.names(binary_matrix) %in% meta_sub$ag_id)
print(paste0('binary matrix has ',nrow(binary_matrix), ' rows and ', ncol(binary_matrix), ' columns after removing samples'))
print(colSums(binary_matrix))
binary_sub <- binary_matrix[colSums(binary_matrix) > 0]


count_file = paste0(argv$index_dir,argv$count_file)
print(paste0('loading counts at ',count_file))
counts <- np$load(count_file, mmap_mode='r')
print('converting to data frame')
counts = data.frame(counts %>% t(), row.names = sample_order$X1)
print(paste0('count file has ',nrow(counts), ' rows and ', ncol(counts), ' columns'))
print('subsetting to sample groups')
counts <- counts %>% filter(row.names(counts) %in% meta_sub$ag_id)
print(paste0('counts matrix has ',nrow(counts), ' rows and ', ncol(counts), ' columns after removing samples'))
counts_sub <- counts[colSums(binary_matrix) > 0]
print(counts_sub)
print(paste0('dimensions of count matrix after subsetting is: ', nrow(counts_sub), ' rows and ', ncol(counts_sub), ' columns'))

sf_file = paste0(argv$index_dir, argv$scale_factors)
print(paste0('loading scale factors at ',sf_file))
norm_factors <- np$load(sf_file, mmap_mode='r')
print('converting to data frame')
norm_factors = data.frame(norm_factors %>% t(), row.names = sample_order$X1)
print(paste0('scale factors file has ',nrow(norm_factors), ' rows and ', ncol(norm_factors), ' columns'))
print('subsetting to sample groups')
norm_factors <- norm_factors %>% filter(row.names(norm_factors) %in% meta_sub$ag_id)
print(paste0('scale factors matrix has ',nrow(norm_factors), ' rows and ', ncol(norm_factors), ' columns after removing samples'))
norm_factors <- norm_factors[colSums(binary_matrix) > 0]
print(paste0('dimensions of normalization matrix after subsetting is: ', str(dim(norm_factors))))
print('converting to data matrix')
norm_factors = norm_factors %>% t() %>% data.matrix()

print(formula(paste("~",argv$meta_column)))
print(argv$meta_column)

dds <- DESeqDataSetFromMatrix(countData = counts_sub %>% t(),
                              colData = meta_sub,
                              design= formula(paste("~",argv$meta_column)) )
normalizationFactors(dds) <- norm_factors

dds <- DESeq(dds, parallel=TRUE)
data_file <- paste0(argv$output_dir, cluster_group,'_deseq_wnorm_postdeseq.rds')
saveRDS(dds, data_file)

res <- results(dds, parallel = TRUE)
data.frame(res) %>% write.csv(file = paste0(argv$output_dir, cluster_group, '_deseq_wnorm_results.tsv'),sep='\t')


# apply(X = contrast_files_pivot, MARGIN = 1, FUN = run_deseq)


# data_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/subsets/stratified_samples/'
# results_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/new_normalization/'
# already_done <- dir(results_dir) %>% data.frame() %>% filter(grepl('.tsv',.)) %>%
#     mutate('contrast' = str_extract_all(., '^.+_', simplify=FALSE) %>% paste0() %>% str_replace('_deseq_wnorm_results_', ''))
# contrast_files <- contrast_files %>% filter(contrast %in% number_of_files$contrast) %>% separate(contrast,c('cluster1', 'cluster2'), '_', remove = FALSE) %>%
#     filter(!is.na(cluster2)) %>% subset(!(contrast %in% c(already_done$contrast)))
# print(contrast_files)
# contrast_files_pivot <- contrast_files %>% pivot_wider(id_cols = c(contrast, cluster1, cluster2), names_from = matrix_type, values_from = c(full_path))

