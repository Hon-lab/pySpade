#DErand.py

import os
import sys
import re
import collections
import argparse
import tables
import itertools
import time
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse
import glob

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from random import shuffle

from pySpade.utils import get_logger, read_annot_df, cpm_normalization, metacelll_normalization, perform_DE, hypergeo_test

np.random.seed(0)
logger = get_logger(logger_name=__name__)

def DE_random_cells(sub_df_file,
                    sgrna_df,
                    num_cell,
                    RAND_M,
                    output_dir,
                    iteration=1000,
                    num_processing=1,
                    norm='cpm'):
    
    #check the normalization method
    if (norm != 'cpm') and (norm != 'metacell'):
        logger.critical("Incorrect normalization method. Has to be either 'cpm' or 'metacell'")
        sys.exit(0)

    #check the normalization method
    if (RAND_M != 'equal') and (RAND_M != 'sgrna'):
        logger.critical("Incorrect randomization method. Has to be either 'equal' or 'sgrna'")
        sys.exit(0)

    #read the 10X hdf5 file
    logger.info('Reading transcriptome file')
    if (sub_df_file.endswith('pkl')):
        sub_df = pd.read_pickle(sub_df_file)
    elif (sub_df_file.endswith('h5')):
        sub_df = pd.read_hdf(sub_df_file, 'df')
    
    #read the plotting annotation
    annot_df = read_annot_df()
    
    #filter the genes
    nonzero_idx = np.where(np.sum(sub_df > 0, axis = 1) > 1)[0]
    idx = list(set(nonzero_idx) & set(annot_df.idx))
    del nonzero_idx
    del annot_df
        
    #normalize the matrix.
    if (norm == 'cpm'):
        cpm_matrix = cpm_normalization(sub_df)
        
    elif (norm == 'metacell'):
        cpm_matrix = metacelll_normalization(sub_df)
     
    logger.info('Finished transcriptome normalization.')
    
    #create input ndarray
    input_array = cpm_matrix[idx]
    del cpm_matrix

    if (sgrna_df.endswith('pkl')):
        sgrna_df_adj_bool = pd.read_pickle(sgrna_df) > 0 
    elif (sgrna_df.endswith('h5')):
        sgrna_df_adj_bool = pd.read_hdf(sgrna_df, 'df') > 0

    [g,c] = sgrna_df_adj_bool.shape
    assert np.sum(sub_df.columns == sgrna_df_adj_bool.columns) == c
    if (RAND_M == 'sgrna'):
        total_choose_num = np.sum(np.sum(sgrna_df_adj_bool))     
        sgrna_weight_list = list(np.sum(sgrna_df_adj_bool, axis=0).values / total_choose_num)
    del sgrna_df_adj_bool
    logger.info('Finished processing sgRNA df.')

    combined_up_array = np.empty([iteration, len(sub_df.index)])
    combined_down_array = np.empty([iteration, len(sub_df.index)])
    combined_cpm_array = np.empty([iteration, len(sub_df.index)])
     
        
    cell_ID_dict = defaultdict(list)
    #idx index of cells containing the given sgRNA
    #random choose cells considering sgrna number per cell 
    logger.info('Start randomization Calculation.')
    for i in np.arange(iteration):
        counter = i
        if (counter % 50 == 0):
            logger.info(str(counter))

        if (RAND_M == 'sgrna'):
            sgrna_idx = np.random.choice(np.arange(c), size=num_cell, replace=False, p=sgrna_weight_list) #number (select from idx)
        elif (RAND_M == 'equal'):
            sgrna_idx = np.random.choice(np.arange(c), size=num_cell, replace=False)
                
        cell_ID_dict[i].append(list(sgrna_idx))
        #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
        pval_list_down = np.zeros(len(sub_df.index))
        pval_list_up = np.zeros(len(sub_df.index))
        fc_list = np.ones(len(sub_df.index))
        cpm_list = np.zeros(len(sub_df.index))
            
        #perform the differential gene analysis by using Virtual FACS
        num_sgrna_cell, pval_list_up, pval_list_down, fc_list, cpm_list = perform_DE(
                                                                            sgrna_idx,
                                                                            input_array,
                                                                            idx,
                                                                            num_processing,
                                                                            pval_list_down,
                                                                            pval_list_up,
                                                                            fc_list,
                                                                            cpm_list
                                                                        )

        #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
        combined_up_array[(counter)] = pval_list_up
        combined_down_array[(counter)] = pval_list_down
        combined_cpm_array[(counter)] = cpm_list
        
        #write the selected cell ID into file 
        output_file = open(output_dir + '/' + str(num_cell) + '-rand_cell_ID.txt', 'w')
        for i in cell_ID_dict.keys():
            select_cell_ID = cell_ID_dict[i]
            output_file.write(str(i) + '\t' + ','.join(str(j) for j in select_cell_ID) + '\n')
        output_file.close()

        up_pval_matrix   = csr_matrix(combined_up_array)
        down_pval_matrix = csr_matrix(combined_down_array)
        region_cpm_matrix = csr_matrix(combined_cpm_array)
   
        #save all the output
        io.savemat(
            '%s/%s-up_log-pval'%(output_dir, num_cell),
            {'matrix':up_pval_matrix}
        )

        io.savemat(
            '%s/%s-down_log-pval'%(output_dir, num_cell),
            {'matrix':down_pval_matrix}
        )

        io.savemat(
            '%s/%s-cpm'%(output_dir, num_cell),
            {'matrix':region_cpm_matrix}
        )

    #Process 1000 iteration
    logger.info('Start randomization summary.')
    up_distribution_mean = np.mean(np.array(up_pval_matrix.tocsr().todense()), axis=0)
    up_distribution_std = np.std(np.array(up_pval_matrix.tocsr().todense()), axis=0)
    down_distribution_mean = np.mean(np.array(down_pval_matrix.tocsr().todense()), axis=0)
    down_distribution_std = np.std(np.array(down_pval_matrix.tocsr().todense()), axis=0)
    cpm_mean = np.mean(np.array(region_cpm_matrix.tocsr().todense()), axis=0)
    
    up_01_perc = np.percentile(np.array(up_pval_matrix.tocsr().todense()), 0.1, axis=0)
    down_01_perc = np.percentile(np.array(down_pval_matrix.tocsr().todense()), 0.1, axis=0)

    np.save(output_dir + 'Up_dist_mean-%s'%(num_cell), up_distribution_mean)
    np.save(output_dir + 'Up_dist_std-%s'%(num_cell), up_distribution_std)
    np.save(output_dir + 'Down_dist_mean-%s'%(num_cell), down_distribution_mean)
    np.save(output_dir + 'Down_dist_std-%s'%(num_cell), down_distribution_std)
    np.save(output_dir + 'Cpm_mean-%s'%(num_cell), cpm_mean)

    np.save(output_dir + 'Up_01_perc-%s'%(num_cell), up_01_perc)
    np.save(output_dir + 'Down_01_perc-%s'%(num_cell), down_01_perc)

if __name__ == '__main__':
    pass
