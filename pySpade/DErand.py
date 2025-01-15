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
import warnings
warnings.filterwarnings('ignore')

from pySpade.utils import get_logger, read_annot_df, read_sgrna_dict, cpm_normalization, metacelll_normalization, perform_DE, hypergeo_test, hypergeo_test_NC, perform_DE_NC, get_num_processes, find_all_sgrna_cells

np.random.seed(0)
logger = get_logger(logger_name=__name__)

def DE_random_cells(sub_df_file,
                    sgrna_df,
                    sgrnas_file,
                    num_cell,
                    RAND_M,
                    output_dir,
                    iteration=1000,
                    num_processing=1,
                    norm='cpm',
                    background_cells='complement'):

    num_processing = get_num_processes()    
    logger.info(str(num_processing) + ' cpu for parallel computing.')
    
    #check the normalization method
    if (norm != 'cpm') and (norm != 'metacell'):
        logger.critical("Incorrect normalization method. Has to be either 'cpm' or 'metacell'")
        sys.exit(0)

    #check the normalization method
    if (RAND_M != 'equal') and (RAND_M != 'sgrna'):
        logger.critical("Incorrect randomization method. Has to be either 'equal' or 'sgrna'")
        sys.exit(0)

    #read the 10X hdf5 file
    logger.info('Reading transcriptome file.')
    if (sub_df_file.endswith('pkl')):
        trans_df = pd.read_pickle(sub_df_file)
    elif (sub_df_file.endswith('h5')):
        trans_df = pd.read_hdf(sub_df_file, 'df')
    
    #convert trans_df to sparse
    sub_df = trans_df.astype(pd.SparseDtype('int32', 0))
    del trans_df

    #read the plotting annotation
    annot_df = read_annot_df()
    idx = np.arange(0, len(sub_df.index))
  
    #normalize the matrix.
    if (norm == 'cpm'):
        cpm_matrix =  cpm_normalization(sub_df)
        
    elif (norm == 'metacell'):
        cpm_matrix = metacelll_normalization(sub_df)
     
    logger.info('Finished transcriptome normalization.')
    
    if (sgrna_df.endswith('pkl')):
        sgrna_df_adj = pd.read_pickle(sgrna_df) > 0 
    elif (sgrna_df.endswith('h5')):
        sgrna_df_adj = pd.read_hdf(sgrna_df, 'df') > 0

    sgrna_df_adj_bool = sgrna_df_adj.astype(pd.SparseDtype('int16', 0))
    del sgrna_df_adj

    [g,c] = sgrna_df_adj_bool.shape
    assert np.sum(sub_df.columns == sgrna_df_adj_bool.columns) == c
    del sub_df
    if (RAND_M == 'sgrna'):
        total_choose_num = np.sum(np.sum(sgrna_df_adj_bool))     
        sgrna_weight_list = list(np.sum(sgrna_df_adj_bool, axis=0).values / total_choose_num)

    logger.info('Finished processing sgRNA df.')

    combined_up_array = np.empty([iteration, len(idx)])
    combined_down_array = np.empty([iteration, len(idx)])
    combined_cpm_array = np.empty([iteration, len(idx)])
     
        
    cell_ID_dict = defaultdict(list)
    if background_cells != 'complement':
        sgrna_dict  = read_sgrna_dict(sgrnas_file)
        if background_cells in list(sgrna_dict.keys()) == False:
            logger.critical('Background region is not in sgRNA dictionary text file. Please fix the dictionary file or use default parameter. Job cancels.')
            sys.exit(0)

        #Calculate the Negative control cells idx for later use
        #exclude the dropout NC sgRNAs 
        missing_NC_list = list(set(sgrna_dict[background_cells]) - set(sgrna_df_adj_bool.index))
        remain_NC_list = list(set(sgrna_dict[background_cells]) - (set(missing_NC_list)))
        logger.info(str(len(remain_NC_list)) + ' sgRNAs remain to find background cells.')
        sgrna_dict.update({background_cells:remain_NC_list})

        NC_idx = find_all_sgrna_cells(sgrna_df_adj_bool, sgrna_dict[background_cells])
        
        del sgrna_df_adj_bool
        if len(NC_idx) == 0:
            logger.critical('No cell as background. Job cancels.')
            sys.exit(0)
        logger.info(str(len(NC_idx)) + ' cells as background for differential expression analysis.')

        logger.info('Start randomization Calculation.')
        for i in np.arange(iteration):
            counter = i
            if (counter % 25 == 0):
                logger.info('Finished ' + str(counter) + ' iterations.')
            if (RAND_M == 'sgrna'):
                sgrna_idx = np.random.choice(np.setxor1d(np.arange(c), NC_idx), size=num_cell, replace=False, p=[sgrna_weight_list[i] for i in np.setxor1d(np.arange(c), NC_idx)]) #number (select from idx)
            elif (RAND_M == 'equal'):
                sgrna_idx = np.random.choice(np.setxor1d(np.arange(c), NC_idx), size=num_cell, replace=False)
                    
            cell_ID_dict[i].append(list(sgrna_idx))
            #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
            pval_list_down = np.zeros(len(idx))
            pval_list_up = np.zeros(len(idx))
            fc_list = np.ones(len(idx))
            cpm_list = np.zeros(len(idx))
                
            #perform the differential gene analysis
            num_sgrna_cell, pval_list_up, pval_list_down, fc_list, cpm_list = perform_DE_NC(
                                                                                sgrna_idx,
                                                                                NC_idx,
                                                                                cpm_matrix,
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

    else: 
        del sgrna_df_adj_bool
        logger.info('Start randomization Calculation.')
        for i in np.arange(iteration):
            counter = i
            if (counter % 25 == 0):
                logger.info('Finished ' + str(counter) + ' iterations.')
            if (RAND_M == 'sgrna'):
                sgrna_idx = np.random.choice(np.arange(c), size=num_cell, replace=False, p=sgrna_weight_list) #number (select from idx)
            elif (RAND_M == 'equal'):
                sgrna_idx = np.random.choice(np.arange(c), size=num_cell, replace=False)
                    
            cell_ID_dict[i].append(list(sgrna_idx))
            #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
            pval_list_down = np.zeros(len(idx))
            pval_list_up = np.zeros(len(idx))
            fc_list = np.ones(len(idx))
            cpm_list = np.zeros(len(idx))
                
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
    logger.info('Start modeling randomization distribution.')
    cpm_mean = np.mean(np.array(region_cpm_matrix.tocsr().todense()), axis=0)

    #model distribution with gamma distribution
    A = np.zeros(len(idx))
    B = np.zeros(len(idx))
    C = np.zeros(len(idx))
    D = np.zeros(len(idx))
    E = np.zeros(len(idx))
    F = np.zeros(len(idx))
    
    logger.info('Process up-regulation.')
    for i in idx:
        if (i % 10000 == 0):
            logger.info('Finished ' + str(i) + ' genes.')

        up_array = -np.array(up_pval_matrix.tocsr()[:, i].todense()).flatten()    
        if np.sum(np.isfinite(up_array)) <= 50:
            continue
        if np.all(up_array[np.isfinite(up_array)] == up_array[np.isfinite(up_array)][0]):
            a, b, c = (np.nan, 0.0, 1.0)
        else: 
            a, b, c = stats.gamma.fit(up_array[np.isfinite(up_array)])
        A[i] = a 
        B[i] = b
        C[i] = c

    logger.info('Process down-regulation.')
    for i in idx:
        if (i % 10000 == 0):
            logger.info('Finished ' + str(i) + ' genes.')

        down_array = -np.array(down_pval_matrix.tocsr()[:, i].todense()).flatten()
        if np.sum(np.isfinite(down_array)) <= 50:
            continue
        if np.all(down_array[np.isfinite(down_array)] == down_array[np.isfinite(down_array)][0]):
            d, e, f = (np.nan, 0.0, 1.0)
        else: 
            d, e, f = stats.gamma.fit(down_array[np.isfinite(down_array)])
        D[i] = d
        E[i] = e
        F[i] = f

    #save the parameters of gamma distribution     
    np.save(output_dir + 'Up_dist_gamma-%s-A'%(num_cell), A)
    np.save(output_dir + 'Up_dist_gamma-%s-B'%(num_cell), B)
    np.save(output_dir + 'Up_dist_gamma-%s-C'%(num_cell), C)
    np.save(output_dir + 'Down_dist_gamma-%s-A'%(num_cell), D)
    np.save(output_dir + 'Down_dist_gamma-%s-B'%(num_cell), E)
    np.save(output_dir + 'Down_dist_gamma-%s-C'%(num_cell), F)

    np.save(output_dir + 'Cpm_mean-%s'%(num_cell), cpm_mean)

    logger.info('Job is done.')
    
if __name__ == '__main__':
    pass
