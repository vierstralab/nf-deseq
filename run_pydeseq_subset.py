#!/bin/bash
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from os import chdir
import os
import re
import subprocess
from datetime import datetime
import time
import multiprocessing
import scipy
import itertools
from datetime import datetime
date = datetime.now().strftime('%Y-%m-%d')

output_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/index_K562_pluripotent/'
meta = pd.read_csv('/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/index_atleast3samples/240123_meta_wcluster.csv', index_col=0)
index_dir = '/net/seq/data2/projects/acote/projects/encode4/4025_index/subsets/reproducible_DHS/'

norm_matrix = np.load(index_dir + 'sf_matrix.npy')
bin_matrix = np.load(index_dir + 'bin_matrix_wsf.npy')
count_matrix = np.load(index_dir + 'count_matrix_wsf.npy')
sample_order = pd.read_table(index_dir + 'sample_order.txt', header=None).reset_index()
meta = pd.read_csv('/net/seq/data2/projects/acote/projects/encode4/4025_index/deseq/index_atleast3samples/240123_meta_wcluster.csv', index_col=0)

meta_sub = meta[meta['cluster'].isin(['Myeloid leukemia (K562)', 'Pluripotent/Pluripotent-derived'])]
sample_sub = sample_order[sample_order[0].isin(meta_sub.index)]



bin_matrix = bin_matrix[:, sample_sub['index']]
rowsum = bin_matrix.sum(-1)
bin_matrix = bin_matrix[rowsum >= 1, :]
print('subsetted binary matrix shape')
print(bin_matrix.shape)

count_matrix = count_matrix[:, sample_sub['index']]
count_matrix = count_matrix[rowsum >= 5, :]
print('subsetted binary matrix shape')
print(bin_matrix.shape)

norm_matrix = norm_matrix[:, sample_sub['index']]
norm_matrix = norm_matrix[rowsum >= 5, :]
print('subsetted sf matrix shape')
print(norm_matrix.shape) 

np.save(output_dir + 'bin_matrix_sf_atleast1sample.npy', bin_matrix)
np.save(output_dir + 'count_matrix_sf_atleast1sample.npy', count_matrix)
np.save(output_dir + 'sf_matrix_atleast1sample.npy', norm_matrix)