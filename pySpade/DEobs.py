#DEobs.py

import os
import sys
import re
import collections
import argparse
import tables
import itertools
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse
from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from random import shuffle

from pySpade.utils import get_logger, read_annot_df, cpm_normalization, metacelll_normalization, read_sgrna_dict, find_all_sgrna_cells, perform_DE, hypergeo_test

np.random.seed(0)
logger = get_logger(logger_name=__name__)

def DE_observe_cells(sub_df_file,
                    sgrna_df,
                    sgrnas_file,
                    output_dir,
                    num_processing=1,
                    norm='cpm',
                    ):

    #check the normalization method   
    if (norm != 'cpm') and (norm != 'metacell'):
        logger.critical("Incorrect normalization method. Has to be either 'cpm' or 'metacell'")
        sys.exit(0)

    #read the 10X hdf5 file
    logger.info('Reading transcriptome file')

    if (sub_df_file.endswith('pkl')):
        sub_df = pd.read_pickle(sub_df_file)
    elif (sub_df_file.endswith('h5')):
        sub_df = pd.read_hdf(sub_df_file, 'df')
       
    idx = np.arange(0, len(sub_df.index))

    #normalize the matrix.       
    if (norm == 'cpm'):
        cpm_matrix = cpm_normalization(sub_df)
        
    elif (norm == 'metacell'):
        cpm_matrix = metacelll_normalization(sub_df)
     
    logger.info('Finished transcriptome normalization.')
    

    #load the sgRNA file
    logger.info('Loading sgRNA df.')
    if (sgrna_df.endswith('pkl')):
        sgrna_df_adj_bool = pd.read_pickle(sgrna_df) > 0 
    elif (sgrna_df.endswith('h5')):
        sgrna_df_adj_bool = pd.read_hdf(sgrna_df, 'df') > 0

    [g,c] = sub_df.shape
    assert np.sum(sub_df.columns == sgrna_df_adj_bool.columns) == c #make sure the sequences are the same compare with two df
    
    del sub_df
    #perform hypergeometric test for every single gene in the dataframe
    sgrna_dict  = read_sgrna_dict(sgrnas_file)

    logger.info('Start DE analysis.')
    for k in sgrna_dict:
        logger.info(f'Start analysis of region: {k}')

        #Find if any sgRNA is missing from the sgRNA df
        missing_sgrna_list = list(set(sgrna_dict[k]) - set(sgrna_df_adj_bool.index))
        if len(missing_sgrna_list) == len(sgrna_dict[k]):
            logger.info('Missing all the sgRNA in this region.')
            continue
        if len(missing_sgrna_list) > 0:
            logger.info('Missing ' + str(len(missing_sgrna_list)) + ' sgrna.')
            for i in missing_sgrna_list:
                logger.info('Missing sgrna: ' + str(i))

        #remove the missing sgrna and continue DE analysis 
        if len(missing_sgrna_list) > 0:
            remain_sgrna_list = list(set(sgrna_dict[k]) - (set(missing_sgrna_list)))
            sgrna_dict.update({k:remain_sgrna_list})

        #idx index of cells containing the given sgRNA
        sgrna_idx = find_all_sgrna_cells(sgrna_df_adj_bool, sgrna_dict[k])

        #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
        pval_list_down = np.zeros(g)
        pval_list_up = np.zeros(g)

        fc_list = np.ones(g)
        cpm_list = np.zeros(g)
       
        #perform the differential gene analysis
        num_sgrna_cell, pval_list_up, pval_list_down, fc_list, cpm_list = perform_DE(
            sgrna_idx,
            cpm_matrix,
            idx,
            num_processing,
            pval_list_down,
            pval_list_up,
            fc_list,
            cpm_list
        )

        #save all the output
        io.savemat(
            '%s/%s-up_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_up})

        io.savemat(
            '%s/%s-down_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_down})

        io.savemat('%s/%s-%s-foldchange'%(output_dir, k[0:], num_sgrna_cell),
                   {'matrix':fc_list})
        
        io.savemat(
            '%s/%s-cpm'%(output_dir, k[0:]),
            {'matrix':cpm_list})

        

if __name__ == '__main__':
    pass
