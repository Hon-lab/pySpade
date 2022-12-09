#gene_expr_level.py

import os
import re
import sys
import collections
import argparse
import tables
import itertools
import matplotlib
import numba

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix

from importlib import resources
from pySpade.utils import get_logger

logger = get_logger(logger_name=__name__)

def Calculate_gene_expression(TRANSCRIPTOME_DF, 
                              GENE_FILE, 
                              OUTPUT_FILE):

    if (TRANSCRIPTOME_DF.endswith('pkl')):
        sub_df = pd.read_pickle(TRANSCRIPTOME_DF)
    elif (TRANSCRIPTOME_DF.endswith('h5')):
        sub_df = pd.read_hdf(TRANSCRIPTOME_DF, 'df')

    #load quired genes
    GENE_LIST = []
    with open(GENE_FILE) as ft:
        for line in ft:
            gene = line.strip()
            GENE_LIST.append(gene)
    

    total_seq_reads = sub_df.sum(axis=0)
    logger.info('Finished loading data.')

    #calculate average cpm and write to output file
    output_file=open(OUTPUT_FILE, 'w')
    output_file.write('Gene' + '\t' + 'Average expression (cpm)' + '\t' + 'Median cpm' + '\t' + 'Portion of cell express' + '\n')
    for gene in GENE_LIST:
        if (gene in sub_df.index) == False:
            logger.critical(str(gene) + ' is missing from the transcriptome dataframe.')
            continue
        dist = sub_df.loc[gene ,:] / total_seq_reads * 1000000
        ave_cpm = np.mean(dist)
        median_cpm = np.median(dist)
        perc_cell = np.sum(dist > 0) /len(dist)
        output_file.write(str(gene) + '\t' + str(ave_cpm) + '\t' + str(median_cpm) + '\t' + str(perc_cell) + '\n')

        logger.info('Finish processing ' + str(gene))

    output_file.close()


if __name__ == '__main__':
    pass
