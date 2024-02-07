#!/bin/bash
import os
import pandas as pd
import numpy as np
import itertools
from optparse import OptionParser


parser = OptionParser()
(options, args) = parser.parse_args()

meta = args[0]
coi = args[1]


meta = pd.read_csv(meta, index_col=0)

clusters = meta[coi].drop_duplicates()
clusters_comb = itertools.combinations(clusters, 2)

for i in clusters_comb:
	print(i[0] + ',' + i[1])
