# library(variancePartition)
library(reticulate)
np <- import("numpy")
library(tidyverse)
library(R.utils)
library(data.table)
library(DESeq2)

data_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/index_atleast5samples/'

norm_file <- paste0(data_dir, 'sf_matrix_atleast5samples.npy')
count_file <- paste0(data_dir, 'count_matrix_atleast3samples.npy')

norm_factors <- np$load(norm_file)
counts <- np$load(count_file)


meta <- read.csv('/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/index_atleast3samples/240123_meta_wcluster.csv')

# sample_list = read.csv('/home/acote/projects/encode4/data/metadata/2023-08-07_meta_with_new_cat.csv') %>% select(-X)
# ordered_sample_list = read.table('/net/seq/data2/projects/sabramov/SuperIndex/dnase-0520/output/indivs_order.txt', sep = '\t') %>% 
#     t() %>% data.frame() %>% filter(. %in% sample_list$id)
# ordered_sample_list$.

norm_factors = data.frame(norm_factors, row.names = meta$id)
norm_factors
counts = data.frame(counts, row.names = meta$id)
counts

sample_listo <- sample_list[match(ordered_sample_list$., sample_list$id), ]
head(sample_listo)


dds <- DESeqDataSetFromMatrix(countData = counts %>% t(),
                              colData = meta,
                              design= ~ cluster)

normalizationFactors(dds) <- norm_factors %>% t()
dds

data_file <- paste0(data_dir, '240206_deseq_wnorm.rds')
saveRDS(dds, data_file)

dds <- DESeq(dds)
data_file <- paste0(data_dir, '240206_deseq_wnorm_postdeseq.rds')
saveRDS(dds, data_file)


res <- results(dds)
data_file <- paste0(data_dir, '240206_deseq_wnorm_results.rds')
saveRDS(dds, data_file)

res