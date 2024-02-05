#!/bin/bash
import os
import pickle as pkl

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pydeseq2.utils as deautils

from datetime import datetime
date = datetime.now().strftime('%Y-%m-%d')

from optparse import OptionParser


parser = OptionParser()
parser.add_option("-c", "--counts", dest="count_matrix",
                  help="count matrix for deseq")
parser.add_option("-m", "--meta", dest="meta",
                  help="metadata file")
parser.add_option("-d", "--design", dest="design_factors",
                  help="metadata column for comparisons")
parser.add_option("-n", "--ncpus", dest="n_cpus", default=8,
                  help="number of cpus to use for main deseq run")
parser.add_option("-o", "--output", dest="output_dir",
                  help="output folder for pkl version of deseq dataset")

(options, args) = parser.parse_args()

output_prefix = args[0]

dds = DeseqDataSet(
    counts=options.count_matrix,
    metadata=options.meta,
    design_factors=options.design_factors,
    refit_cooks=True,
    n_cpus=options.n_cpus
)

dds.deseq2()

data_file = os.path.join(options.output_dir, date + '_' + output_prefix + "_dds_cluster.pkl")

with open(data_file, "wb") as f:
    pkl.dump(dds, f)

# print(data_file)