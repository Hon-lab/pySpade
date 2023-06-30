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

from pySpade.utils import get_logger, read_annot_df, get_neighbor_genes, get_distance, read_sgrna_dict

logger = get_logger(logger_name=__name__)

def local_analysis(FILE_DIR,
                    OBS_DIR,
                    DISTRI_DIR,
                    SGRNA_DICT,
                    OUTPUT_DF):
    
    logger.info('Loading files.')

    #read the gene sequence file 
    gene_seq = np.load(FILE_DIR + 'Trans_genome_seq.npy', allow_pickle=True)
    if len(gene_seq) != len(set(gene_seq)):
        logger.critical('Duplication of mapping genes.')
    #read the plotting annotation
    annot_df_dup = read_annot_df()
    #There are many non-coding genes duplication in the annot_df, only keep one.
    annot_df = annot_df_dup.drop_duplicates(subset='gene_names', keep='first')
    
    #Load sgRNA dict: All regions 
    sgrna_dict  = read_sgrna_dict(SGRNA_DICT)

    #read all the perturbation files
    pval_files = glob.glob(OBS_DIR + '*-down_log-pval')

    df_column_list = [
        'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 'region', 'distance', 'num_cell', 'bin', 
        'pval', 'fc', 'padj-Gaussian', 'fc_by_rand_dist_cpm', 'cpm_perturb', 'cpm_bg']
    local_gene_df = pd.DataFrame(columns=df_column_list)

    #Read the background distribution file
    Num = [int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]
    if len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_mean-*')]) != \
        len([int(i.split('/')[-1].split('.')[0].split('-')[-1]) for i in glob.glob(DISTRI_DIR+ '/Up_dist_mean-*')]):
        logger.critical('Background distribution files error!')
        sys.exit(0)

    #start calculating local hits for each perturbation region
    logger.info('Start analysis of ' + str(len(pval_files)) + ' regions.')
    for region in list(sgrna_dict.keys()):
        logger.info(f'  Processing region: {region}')
        if region.startswith('chr') == False:   
            logger.info('No chromosome coordination in this region, cannot compute local analysis.')
            continue

        cpm_file = OBS_DIR + region + '-cpm'
        cpm = sio.loadmat(cpm_file)['matrix'][0]        
        fc_files = glob.glob(OBS_DIR + region + '*-foldchange')
        if len(fc_files) == 1:
            fc_file = fc_files[0]
            cell_num = int(fc_file.split('/')[-1].split('-')[-2])
        else:    
            numbers = np.array([int(float(i.split('/')[-1].split('-')[2])) for i in fc_files])
            chosen_num = np.max(numbers)
            cell_num = int(chosen_num)
            fc_file = OBS_DIR + region + '-' + str(chosen_num) + '-foldchange'
        
        fc = sio.loadmat(fc_file)['matrix'][0]  #compare to all the other cells as background

        #get the gene idx within local analysis window and filter with fold change 
        local_gene = get_neighbor_genes(region, 2e6, annot_df)
        fc_cutoff = 0.01
        #Calculate the overlap genes with annot_df, and only save the information on those genes.
        down_idx = np.where(np.array(fc) < (1 - fc_cutoff))[0]
        unique_elements, unique_indices = np.unique(gene_seq, return_index=True)
        down_keep_genes = list((set(annot_df['gene_names']).intersection(set(gene_seq[down_idx]))).intersection(set(local_gene)))
        down_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, down_keep_genes))[0]]))

        if len(down_keep_genes_idx) ==0:
            logger.info(f'  No down-regulation genes within local analysis windown. ')
            continue
            
        #read the pval matrix
        pval = sio.loadmat(OBS_DIR + region + '-down_log-pval')['matrix'][0] #raw hypergeom p value
        
        #Calculate the adjusted pval with closet cell number distribution
        chosen_dist = Num[np.argmin([np.absolute(n - cell_num) for n in Num])]
        
        #Load the file
        down_mean = np.load(DISTRI_DIR + 'Down_dist_mean-%s.npy'%(str(chosen_dist)))
        down_std = np.load(DISTRI_DIR + 'Down_dist_std-%s.npy'%(str(chosen_dist)))
        cpm_mean = np.load(DISTRI_DIR + 'Cpm_mean-%s.npy'%(str(chosen_dist)))
        
        fc_rand = cpm/cpm_mean
        
        #Calculate adjusted p value of Gaussian padj
        down_zscore_list = (pval[down_keep_genes_idx] - down_mean[down_keep_genes_idx]) / down_std[down_keep_genes_idx]
        down_padj_list = scipy.stats.norm.logsf(abs(down_zscore_list))
        hits_fc = fc[down_keep_genes_idx]
        hits_fc_rand = fc_rand[down_keep_genes_idx]

        #save to csv file 
        local_gene_series = annot_df[annot_df['gene_names'].isin(down_keep_genes)].set_index('idx')
        local_gene_series['region'] = region
        dist_list = []
        for i in local_gene_series['pos']:
            dist = get_distance(region, i)
            dist_list.append(dist)
        local_gene_series['distance'] = dist_list
        local_gene_series['num_cell'] = cell_num
        local_gene_series['bin'] = chosen_dist
        local_gene_series['pval'] = pval[down_keep_genes_idx] 
        local_gene_series['fc'] = hits_fc
        local_gene_series['padj-Gaussian'] = down_padj_list
        local_gene_series['fc_by_rand_dist_cpm'] = hits_fc_rand
        local_gene_series['cpm_perturb'] = cpm[down_keep_genes_idx]
        local_gene_series['cpm_bg'] = cpm_mean[down_keep_genes_idx]
        local_gene_df = local_gene_df.append(local_gene_series)        
        local_gene_df = local_gene_df.reindex(columns=df_column_list)

    if OUTPUT_DF.endswith('.csv'):
        local_gene_df.to_csv(OUTPUT_DF)
    else:
        local_gene_df.to_csv(OUTPUT_DF + '.csv')

if __name__ == '__main__':
    pass
