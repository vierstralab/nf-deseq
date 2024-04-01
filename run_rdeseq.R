#!/bin/bash

library(reticulate)
np <- import("numpy")
library(tidyverse)
library(R.utils)
library(data.table)
library(DESeq2)
library(readxl)

today <-
  Sys.Date()
today <- format(today, format="%y%m%d")

library("BiocParallel")
register(MulticoreParam(40))


#data_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/subsets/byclusters/'
index_dir = '/net/seq/data2/projects/sabramov/SuperIndex/clean_index_clustering/dnase_3501/stratified_sampling.downsampled.0326/subset/output/'
data_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/subsets/stratified_samples/'
results_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/new_normalization/'

contrast_files <- dir(data_dir) %>% data.frame() %>% filter(grepl('.npy',.)) %>% mutate('full_path' = paste0(data_dir, .)) %>% 
    mutate('contrast' = str_extract_all(., '^.+_', simplify=FALSE) %>% str_replace('_binary.only.+|_counts.only.+|_normalized.only.+', ''),
          'matrix_type' = str_extract(.,'counts|scale_factor|binary'))
number_of_files  <- contrast_files %>% group_by(contrast) %>% summarize('n' = n()) %>% filter(n == 3)


already_done <- dir(results_dir) %>% data.frame() %>% filter(grepl('.tsv',.)) %>%
    mutate('contrast' = str_extract_all(., '^.+_', simplify=FALSE) %>% paste0() %>% str_replace('_deseq_wnorm_results_', ''))


contrast_files <- contrast_files %>% filter(contrast %in% number_of_files$contrast) %>% separate(contrast,c('cluster1', 'cluster2'), '_', remove = FALSE) %>%
    filter(!is.na(cluster2)) %>% subset(!(contrast %in% c(already_done$contrast)))
print(contrast_files)
contrast_files_pivot <- contrast_files %>% pivot_wider(id_cols = c(contrast, cluster1, cluster2), names_from = matrix_type, values_from = c(full_path))
# contrast_files_pivot

sample_order <- read_table(paste0(index_dir, 'samples_order.txt'), col_names=FALSE)
sample_order


test_meta <- read_excel('/home/jvierstra/proj/ENCODE4/meta+sample_ids+cluster_id.xlsx') %>% filter(ag_id %in% sample_order$X1)

test_meta <- test_meta[ order(match(test_meta$ag_id, sample_order$X1)), ]

test_meta$combined_annotation <- paste(test_meta$core_annotation, test_meta$extended_annotation, sep = ' ')

meta <- test_meta %>% mutate('stripped_cluster' = str_replace_all(combined_annotation, regex("\\W+"), ''))

print(meta %>% group_by('stripped_cluster') %>% summarize('n'=n()))

#meta <- read.csv('/net/seq/data2/projects/sabramov/ENCODE4/dnase-cavsv2/dnase.v3/meta+sample_ids+cluster_id.tsv', sep = '\t') %>%
#    mutate('stripped_cluster' = str_replace_all(cluster, regex("\\W+"), ''))

run_deseq <- function(row){
    meta_sub <- meta %>% filter(stripped_cluster %in% c(row['cluster1'], row['cluster2']))
    print(meta_sub$ag_id)
    print(nrow(meta_sub))
    print(ncol(meta_sub))
    cluster_group <- row['contrast']
    print(cluster_group)
    norm_file <- paste0(data_dir,cluster_group,'_normalized.only_autosomes.filtered.scale_factors.npy')
    print(norm_file)
    count_file <- paste0(data_dir, cluster_group,'_counts.only_autosomes.filtered.matrix.npy')
    print(count_file)
    norm_factors <- np$load(norm_file, mmap_mode='r')
    counts <- np$load(count_file, mmap_mode='r')
    norm_factors = data.matrix(norm_factors)
    dim(norm_factors)
    counts = data.frame(counts %>% t(), row.names = meta_sub$ag_id)
    dim(counts)
    dds <- DESeqDataSetFromMatrix(countData = counts %>% t(),
                                  colData = meta_sub,
                                  design= ~ stripped_cluster)
    # data_file <- paste0(results_dir, cluster_group,'_deseq_dataset.rds')
    # saveRDS(dds, data_file)
    normalizationFactors(dds) <- norm_factors
    # data_file <- paste0(results_dir, cluster_group, '_deseq_wnorm.rds')
    # saveRDS(dds, data_file)
    
    dds <- DESeq(dds, parallel=TRUE)
    data_file <- paste0(results_dir, cluster_group,'_deseq_wnorm_postdeseq_stratified.rds')
    saveRDS(dds, data_file)
    
    res <- results(dds, parallel = TRUE)
    data.frame(res) %>% write.csv(file = paste0(results_dir, cluster_group, '_deseq_wnorm_results_stratified.tsv'),sep='\t')

}

apply(X = contrast_files_pivot, MARGIN = 1, FUN = run_deseq)
