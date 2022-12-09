#global.py
import os
import re
import sys
import collections
import argparse
import tables
import itertools 
import scipy
import matplotlib
import csv
import glob

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse
import scipy.io as sio

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.io as io

from pySpade.utils import get_logger, read_annot_df, read_sgrna_dict, load_data

logger = get_logger(logger_name=__name__)

def global_analysis(data_dir,
                    SGRNA_DICT,
                    DISTRI_DIR,
                    OUTPUT_DF):
    
    logger.info('Loading files.')

    #read the plotting annotation
    annot_df = read_annot_df()

    ##Load sgRNA dict: All regions 
    sgrna_dict  = read_sgrna_dict(SGRNA_DICT)

    #Save the randomized bin with list 
    
    Num = [int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]
    if len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]) != \
        len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Up_dist_mean-*')]):
        logger.critical('Background distribution files error!')
        sys.exit(0)

    logger.info(f'{len(Num)} sets of background distribution for global analysis.')

    #Calcualte pval
    df_column_list = [
        'idx', 'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 
        'region', 'num_cell', 'bin',
        'pval', 'fc', 'padj-Gaussian', 'fc_by_rand_dist_cpm']
    global_hits_df = pd.DataFrame(columns=df_column_list)

    #Generate bin-region dictionary 
    bin_dict = defaultdict(list)
    for region in list(sgrna_dict.keys()):
        fc_files = glob.glob(data_dir + region + '*-foldchange')
        fc_file = fc_files[0]
        cell_num = int(fc_file.split('/')[-1].split('-')[2])
        chosen_bin = Num[np.argmin([np.absolute(n-cell_num) for n in Num])]
        bin_dict[chosen_bin].append(region)
    logger.info('Finish generating bin-region dictionary')
    
    pval_cutoff = -1
    fc_cutoff = 0.1

    for b in Num:
        logger.info(f'Processing bin: {b}')
        logger.info(f'Num of regions in this bin: {len(bin_dict[b])}')

        #Load the background files
        down_mean = np.load(DISTRI_DIR + 'Down_dist_mean-%s.npy'%(str(b)))
        down_std = np.load(DISTRI_DIR + 'Down_dist_std-%s.npy'%(str(b))) 
        up_mean = np.load(DISTRI_DIR + 'Up_dist_mean-%s.npy'%(str(b)))
        up_std = np.load(DISTRI_DIR + 'Up_dist_std-%s.npy'%(str(b)))
        cpm_mean = np.load(DISTRI_DIR + 'Cpm_mean-%s.npy'%(str(b)))

        for region in bin_dict[b]:
            pval_list_up, pval_list_down, cpm, fc, cell_num = load_data(data_dir, region)
            fc_cpm = (cpm + 0.01)/(cpm_mean + 0.01)

            #preprocess pval list 
            pval_list_up[np.argwhere(pval_list_up == 0)] = 1
            pval_list_up[np.isinf(pval_list_up)] = 0
            pval_list_up[np.isnan(pval_list_up)] = 0
            pval_list_down[np.argwhere(pval_list_down == 0)] = 1
            pval_list_down[np.isinf(pval_list_down)] = 0
            pval_list_down[np.isnan(pval_list_down)] = 0

            #split up-regulation and down-regulation
            up_i = np.where(np.array(fc_cpm) > (1 + fc_cutoff))[0]
            up_idx = np.array(list(set(up_i).intersection(set(np.where(pval_list_up < pval_cutoff)[0]))))

            down_i = np.where(np.array(fc_cpm) < (1 - fc_cutoff))[0]
            down_idx = np.array(list(set(down_i).intersection(set(np.where(pval_list_down < pval_cutoff)[0]))))

            #Calculate p-value adj for Gaussian method 
            down_zscore_list = (pval_list_down[down_idx] - down_mean[down_idx]) / down_std[down_idx]
            down_padj_list = scipy.stats.norm.logsf(abs(down_zscore_list))
            down_hit_fc_list = fc_cpm[down_idx]

            up_zscore_list = (pval_list_up[up_idx] - up_mean[up_idx]) / up_std[up_idx]
            up_padj_list = scipy.stats.norm.logsf(abs(up_zscore_list))
            up_hit_fc_list = fc_cpm[up_idx]

            #save to csv file: down-regulation gene 
            global_gene_series = annot_df.set_index('idx').loc[down_idx, :]
            global_gene_series['region'] = region
            global_gene_series['num_cell'] = cell_num
            global_gene_series['bin'] = b
            global_gene_series['pval'] = pval_list_down[down_idx]
            global_gene_series['fc'] = fc[down_idx]
            global_gene_series['padj-Gaussian'] = down_padj_list
            global_gene_series['fc_by_rand_dist_cpm'] = down_hit_fc_list
            global_gene_series['idx'] = global_gene_series.index
            global_hits_df = global_hits_df.append(global_gene_series)

            #save to csv file: up-regulation gene 
            global_gene_series = annot_df.set_index('idx').loc[up_idx, :]
            global_gene_series['region'] = region
            global_gene_series['num_cell'] = cell_num
            global_gene_series['bin'] = b
            global_gene_series['pval'] = pval_list_up[up_idx]
            global_gene_series['fc'] = fc[up_idx]
            global_gene_series['padj-Gaussian'] = up_padj_list
            global_gene_series['fc_by_rand_dist_cpm'] = up_hit_fc_list
            global_gene_series['idx'] = global_gene_series.index
            global_hits_df = global_hits_df.append(global_gene_series)
            
            logger.info(f'Finish region: {region}')
    
    global_hits_df = global_hits_df.reindex(columns=df_column_list)

    if OUTPUT_DF.endswith('.csv'):
        global_hits_df[global_hits_df['padj-Gaussian'] < pval_cutoff].to_csv(OUTPUT_DF)
    else:
        global_hits_df[global_hits_df['padj-Gaussian'] < pval_cutoff].to_csv(OUTPUT_DF + '.csv')

if __name__ == '__main__':
    pass