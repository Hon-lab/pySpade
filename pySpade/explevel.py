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
import glob
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

def Calculate_gene_expression(TRANSCRIPTOME_DIR, 
                              GENE_FILE, 
                              OUTPUT_FILE):

    TRANSCRIPTOME_FILE = glob.glob(TRANSCRIPTOME_DIR + 'Singlet_sub_df.*')
    if TRANSCRIPTOME_FILE[0].endswith('h5') == True:
        sub_df = pd.read_hdf(TRANSCRIPTOME_FILE[0], 'df')
    if TRANSCRIPTOME_FILE[0].endswith('pkl') == True:
        sub_df = pd.read_pickle(TRANSCRIPTOME_FILE[0])

    #load quired genes
    GENE_LIST = []
    if GENE_FILE.endswith('.txt'):
        with open(GENE_FILE) as ft:
            for line in ft:
                gene = line.strip()
                GENE_LIST.append(gene)
    else:
        for i in GENE_FILE.split(','):
            GENE_LIST.append(i)

    total_seq_reads = sub_df.sum(axis=0)
    logger.info('Finished loading data.')

    #calculate average cpm and write to output file
    if OUTPUT_FILE != None:
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

    else:
        logger.info('Gene' + '\t' + 'Average expression (cpm)' + '\t' + 'Median cpm' + '\t' + 'Portion of cell express' + '\n')
        for gene in GENE_LIST:
            if (gene in sub_df.index) == False:
                logger.critical(str(gene) + ' is missing from the transcriptome dataframe.')
                continue
            dist = sub_df.loc[gene ,:] / total_seq_reads * 1000000
            ave_cpm = np.mean(dist)
            median_cpm = np.median(dist)
            perc_cell = np.sum(dist > 0) /len(dist)
            logger.info(str(gene) + '\t' + str(ave_cpm) + '\t' + str(median_cpm) + '\t' + str(perc_cell) + '\n')

    

if __name__ == '__main__':
    pass
