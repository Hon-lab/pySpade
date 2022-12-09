#local.py

import os
import re
import sys
import collections
import argparse
from threading import local
import tables
import itertools 
import scipy
import csv
import glob

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse
import scipy.io as sio

from collections import defaultdict
from scipy import sparse
from scipy.sparse import csr_matrix

from pySpade.utils import get_logger, read_annot_df, get_neighbor_genes, get_distance

logger = get_logger(logger_name=__name__)

def local_analysis(data_dir,
                    DISTRI_DIR,
                    OUTPUT_DF):
    
    logger.info('Loading files.')

    #read the plotting annotation
    annot_df = read_annot_df()

    #read all the perturbation files
    pval_files = glob.glob(data_dir + '*-down_log-pval')

    df_column_list = [
        'idx', 'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 'region', 'dist', 'num_cell', 'bin', 
        'pval', 'fc', 'padj-Gaussian', 'fc_by_rand_dist_cpm']
    local_gene_df = pd.DataFrame(columns=df_column_list)

    #Read the background distribution file
    Num = [int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]
    if len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]) != \
        len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Up_dist_mean-*')]):
        logger.critical('Background distribution files error!')
        sys.exit(0)

    #start calculating local hits for each perturbation region
    logger.info('Start analysis of ' + str(len(pval_files)) + ' regions.')
    for f in pval_files:         
        region = f.split('/')[-1].split('_')[0].replace('-down', '')  
        logger.info(f'  Processing region: {region}')  
        
        cpm_file = data_dir + region + '-cpm'
        cpm = sio.loadmat(cpm_file)['matrix'][0]        
        fc_files = glob.glob(data_dir + region + '*-foldchange')
        if len(fc_files) == 1:
            fc_file = fc_files[0]
            cell_num = int(fc_file.split('/')[-1].split('-')[-2])
        else:    
            numbers = np.array([int(float(i.split('/')[-1].split('-')[2])) for i in fc_files])
            chosen_num = np.max(numbers)
            cell_num = int(chosen_num)
            fc_file = data_dir + region + '-' + str(chosen_num) + '-foldchange'
        
        fc = sio.loadmat(fc_file)['matrix'][0]  #compare to all the other cells as background

        #get the gene idx within local analysis window and filter with fold change 
        local_gene_idx = get_neighbor_genes(region, 2e6, annot_df)
        fc_cutoff = 0.05
        down_i = np.where(np.array(fc) < (1 - fc_cutoff))[0]
        gene_idx = np.array(list(set(down_i).intersection(set(local_gene_idx))))

        if len(gene_idx) ==0:
            logger.info(f'  No down-regulation genes within local analysis windown. ')
            continue
            
        #read the pval matrix
        pval = sio.loadmat(f)['matrix'][0] #raw hypergeom p value
        
        #Calculate the adjusted pval with closet cell number distribution
        chosen_dist = Num[np.argmin([np.absolute(n - cell_num) for n in Num])]
        
        #Load the file
        down_mean = np.load(DISTRI_DIR + 'Down_dist_mean-%s.npy'%(str(chosen_dist)))
        down_std = np.load(DISTRI_DIR + 'Down_dist_std-%s.npy'%(str(chosen_dist)))
        cpm_mean = np.load(DISTRI_DIR + 'Cpm_mean-%s.npy'%(str(chosen_dist)))
        
        fc_rand = cpm/cpm_mean
        
        #Calculate adjusted p value of Gaussian padj
        down_zscore_list = (pval[gene_idx] - down_mean[gene_idx]) / down_std[gene_idx]
        down_padj_list = scipy.stats.norm.logsf(abs(down_zscore_list))
        hits_fc = fc[gene_idx]
        hits_fc_rand = fc_rand[gene_idx]

        #save to csv file 
        local_gene_series = annot_df.set_index('idx').loc[gene_idx, :]

        local_gene_series['region'] = region

        dist_list = []
        for i in local_gene_series['pos']:
            dist = get_distance(region, i, annot_df)
            dist_list.append(int(dist))
        local_gene_series['distance'] = dist_list
        local_gene_series['num_cell'] = cell_num
        local_gene_series['bin'] = chosen_dist
        local_gene_series['pval'] = pval[gene_idx] 
        local_gene_series['fc'] = hits_fc
        local_gene_series['padj-Gaussian'] = down_padj_list
        local_gene_series['fc_by_rand_dist_cpm'] = hits_fc_rand

        local_gene_df = local_gene_df.append(local_gene_series)        

    local_gene_df = local_gene_df.reindex(columns=df_column_list)

    if OUTPUT_DF.endswith('.csv'):
        local_gene_df.to_csv(OUTPUT_DF)
    else:
        local_gene_df.to_csv(OUTPUT_DF + '.csv')

if __name__ == '__main__':
    pass