#process.py

from asyncio.log import logger
import os
import re
import sys
import collections
import argparse
import tables
import itertools
import matplotlib
import numba

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from pySpade.utils import get_logger, get_matrix_from_h5, process_feature_df, turn_point, filter_umi, get_sgrna_per_cell

logger = get_logger(logger_name=__name__)

def process_singlet(WORK_DIR,
                    SGRNA,
                    CELL_MULTIPLEX,
                    OUTPUT_DIR,
                    COMP='False'):
 
    logger.info('Reading data.')

    if (SGRNA.split('.')[-2] != 'pkl') and (SGRNA.split('.')[-1] != 'pkl') and (SGRNA.split('.')[-2] != 'csv') and (SGRNA.split('.')[-1] != 'csv'):
        logger.critical('sgRNA matrix format incorrect')
        sys.exit(0)

    if (COMP != 'True') and (COMP != 'False'):
        logger.critical("Incorrect compression method.")
        sys.exit(0)

    #calculate HTO and filter
    ab_list = []
    with open (CELL_MULTIPLEX) as f:
        for line in f:
            ab = line.strip()
            ab_list.append(ab)
    logger.info('Number of antibodies: ' + str(len(ab_list)))
    
    Count_list = []
    CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
    filtered_matrix_h5 = WORK_DIR + 'filtered_feature_bc_matrix.h5'
    filtered_feature_bc_matrix = get_matrix_from_h5(filtered_matrix_h5)
    ab_idx = np.where(filtered_feature_bc_matrix.feature_ref['feature_type'] == b'Antibody Capture')[0]
    if len(ab_idx) != len(ab_list):
        logger.critical('Input Antibody list is not consistent with mapping info. Using mapping info for further analysis.')
    
    #slice antibody matrix to filter UMI
    df = pd.DataFrame(data=filtered_feature_bc_matrix.matrix[ab_idx, :].todense(), \
                      columns=filtered_feature_bc_matrix.barcodes.astype(str), \
                      index=filtered_feature_bc_matrix.feature_ref['name'][ab_idx].astype(str))
    HTO_df_adj, _ = filter_umi(df, copy=True)
    del df
    Count_list.append(np.sum(HTO_df_adj > 0, axis=0).values)

    HTO_df_adj_bool = HTO_df_adj > 0
    del HTO_df_adj
    singlet_HTO_df = HTO_df_adj_bool.T[(HTO_df_adj_bool.sum(axis=0).values == 1)]

    singlet_rate = collections.Counter(Count_list[0])[1] / len(filtered_feature_bc_matrix.barcodes)
    logger.info('Finished processing experimental doublets.')

    #Load sgrna df 
    if (SGRNA.endswith('pkl')) | (SGRNA.endswith('pkl.gz')):
        sgRNA_df = pd.read_pickle(SGRNA)
    elif (SGRNA.endswith('csv')) | (SGRNA.endswith('csv.gz')):
        sgRNA_df = pd.read_csv(SGRNA, index_col=0)

    sgRNA_df_bool = sgRNA_df > 0

    #split the cells into different antibodies
    singlet_cell_ID = []
    for i in ab_list:
        cell_ID = set(list(singlet_HTO_df.index[singlet_HTO_df[i] == True])).intersection(set(sgRNA_df_bool.columns))
        singlet_cell_ID.append(cell_ID)

    sgRNA_num = []
    for i in np.arange(len(singlet_cell_ID)):
        sgRNA_num.append(list(np.sum(sgRNA_df_bool[singlet_cell_ID[i]], axis=0).values))
    
    del singlet_cell_ID 
    All_singlet_ID = HTO_df_adj_bool.columns[(HTO_df_adj_bool > 0).sum(axis=0).values == 1].values
    new_columns = set(sgRNA_df.columns.values).intersection(set(All_singlet_ID))

    new_sgRNA_df = sgRNA_df[new_columns]
    del sgRNA_df
    logger.info('Generated sgRNA matrix without experimental doublets.')

    #filter sgRNA outliers 
    CorrSgrnaPerCell, sgrna_mean , sgrna_median = get_sgrna_per_cell(new_sgRNA_df, return_mean=True, return_median=True)
    Q1 = CorrSgrnaPerCell.quantile(0.25)
    Q3 = CorrSgrnaPerCell.quantile(0.75)
    IQR = Q3 - Q1
    
    #generate new singlet df
    no_outlier_sgrna_df = new_sgRNA_df.loc[:,((CorrSgrnaPerCell < (Q1 - 1.5 * IQR)) |(CorrSgrnaPerCell > (Q3 + 1.5 * IQR))).values == False]
    del new_sgRNA_df
    CorrSgrnaPerCell_out, sgrna_mean_out , sgrna_median_out = get_sgrna_per_cell(no_outlier_sgrna_df, return_mean=True, return_median=True)
    
    logger.info('Generated transcriptome and sgRNA matrix without sgRNA outlier cells.')

    #write output stats file
    stats_output = open(OUTPUT_DIR + '/stats.txt', 'w')
    stats_output.write('singlet rate: ' + str(singlet_rate) + '\n')
    stats_output.write('Cells before filter experimental doublets: ' + str(len(filtered_feature_bc_matrix.barcodes)) + '\n')
    for a, i in zip(ab_list, np.arange(len(sgRNA_num))):
        stats_output.write('Antibody: ' + a + ' with sgRNA median: ' + str(np.median(sgRNA_num[i])) + '\n')
    stats_output.write('Cells after filter experimental doublets: ' + str(len(new_columns)) + '\n')
    stats_output.write('sgRNA mean before filter sgRNA outlier: ' + str(sgrna_mean) + '\n')
    stats_output.write('sgRNA median before filter sgRNA outlier: ' + str(sgrna_median) + '\n')
    stats_output.write('sgRNA mean after filter sgRNA outlier: ' + str(sgrna_mean_out) + '\n')
    stats_output.write('sgRNA median after filter sgRNA outlier: ' + str(sgrna_median_out) + '\n')
    stats_output.write('Cells after filter experimental doublets and sgRNA outliers: ' + str(len(no_outlier_sgrna_df.columns)) + '\n')
    stats_output.close()

    #write output matrix
    sub_df = pd.DataFrame(data=filtered_feature_bc_matrix.matrix.todense(), columns=filtered_feature_bc_matrix.barcodes.astype(str), index=filtered_feature_bc_matrix.feature_ref['name'].astype(str))
    np.save(OUTPUT_DIR + 'Trans_genome_seq', sub_df.index)
    if (COMP=='True'):
        no_outlier_sgrna_df.to_hdf(OUTPUT_DIR + '/Singlet_sgRNA_df.h5', key='df', mode='w', complevel=9, complib='zlib')    
        sub_df[no_outlier_sgrna_df.columns].to_hdf(OUTPUT_DIR + '/Singlet_sub_df.h5', key='df', mode='w', complevel=9, complib='zlib')
    elif (COMP=='False'):
        no_outlier_sgrna_df.to_hdf(OUTPUT_DIR + '/Singlet_sgRNA_df.h5', key='df', mode='w')    
        sub_df[no_outlier_sgrna_df.columns].to_hdf(OUTPUT_DIR + '/Singlet_sub_df.h5', key='df', mode='w')
    

if __name__ == '__main__':
    pass
