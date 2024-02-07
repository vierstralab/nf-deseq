#!/bin/bash

import os
import pickle as pkl
import itertools
import pandas as pd
import numpy as np
import re
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pydeseq2.utils as deautils
import itertools
from datetime import datetime
from optparse import OptionParser
date = datetime.now().strftime('%Y-%m-%d')


parser = OptionParser()
(options, args) = parser.parse_args()

dds_file = args[0]
coi = args[1]
contrast = args[2]
output_dir = args[3]
output_prefix = args[4]

print(contrast)
contrast_list = contrast.split(',')
print(contrast_list[0])
print(contrast_list[1])
#print([coi, contrast_list[0] + ', ' + contrast[1]])
print([coi, contrast_list[0],  contrast_list[1]])

with open(dds_file, "rb") as f:
    dds = pkl.load(f)

stat_res_B_vs_A = DeseqStats(dds, contrast=[coi, contrast_list[0], contrast_list[1]])
stat_res_B_vs_A.summary()
stat_res_B_vs_A.results_df.to_csv(os.path.join(output_dir, date + '_' + output_prefix + '_' + re.sub(r'\W+', '', contrast_list[0]) + '_vs_' + re.sub(r'\W+', '', contrast_list[1]) + ".csv"))
