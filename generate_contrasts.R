#!/bin/bash
library(tidyverse)
library(R.utils)
library(readxl)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("generate contrasts from provided meta")

# Add command line arguments
p <- add_argument(p, "meta_file", help="metadata spreadsheet, currently excel format")
p <- add_argument(p, "meta_column", help="column of meta file to use for comparisons")

# Parse the command line arguments
argv <- parse_args(p)

# # Rscript round.R 3.141 -d 2
test_meta <- read_excel(argv$meta_file)
test_meta$combined_annotation <- paste(test_meta$core_annotation, test_meta$extended_annotation, sep = ' ')
test_meta <- test_meta %>% mutate('stripped_cluster' = str_replace_all(combined_annotation, regex("\\W+"), ''))
#print(test_meta$stripped_cluster %>% unique())
meta_combinations <- combn(test_meta[[argv$meta_column]] %>% unique(), 2, simplify=FALSE)
#meta_combinations
concat_comb = sapply(meta_combinations, paste, collapse = "_")
print(as.data.frame(concat_comb),row.names=F)
#print(length(meta_combinations))
#for (i in meta_combinations) {
#  print(i)
#}
