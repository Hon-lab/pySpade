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
    chr_order = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
    ]

    #read the gene sequence file 
    gene_seq = np.load(FILE_DIR + 'Trans_genome_seq.npy', allow_pickle=True)
    if len(gene_seq) != len(set(gene_seq)):
        logger.critical('Duplication of mapping genes. Duplicates are removed in the analysis.')
    unique_elements, counts = np.unique(gene_seq, return_counts=True)
    duplicate_elements = unique_elements[counts > 1]
 
    #read the plotting annotation
    annot_df_dup = read_annot_df()
    #There are many non-coding genes duplication in the annot_df, only keep one.
    annot_df = annot_df_dup.drop_duplicates(subset='gene_names', keep='first')
    
    #Load sgRNA dict: All regions 
    sgrna_dict  = read_sgrna_dict(SGRNA_DICT)

    df_column_list = [
        'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 'region', 'distance', 'num_cell', 'bin', 
        'log(pval)-hypergeom', 'fc', 'Significance_score', 'fc_by_rand_dist_cpm', 'pval-empirical', 'cpm_perturb', 'cpm_bg']
    local_gene_df = pd.DataFrame(columns=df_column_list)

    #Read the background distribution file
    Num = list(set([int(i.split('/')[-1].split('.')[0].split('-')[-2]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_gamma-*')]))

    #start calculating local hits for each perturbation region
    for region in list(sgrna_dict.keys()):
        logger.info(f'Processing region: {region}')
        if region.startswith('chr') == False:   
            logger.info('No chromosome coordination in this region, cannot compute local analysis.')
            continue
        if (region.split(':')[0] in chr_order) == False:   
            logger.info('Wrong chromosome coordination in this region, cannot compute local analysis.')
            continue
        if os.path.isfile(OBS_DIR + region + '-cpm') == False:
            logger.info('Skip missing region: ' + region + '. Potential dropout region.')
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
        pval_cutoff = 0

        #read the pval matrix and process
        pval = sio.loadmat(OBS_DIR + region + '-down_log-pval')['matrix'][0] #raw hypergeom p value
        pval[np.argwhere(pval == 0)] = 1
        pval[np.isinf(pval)] = 0
        pval[np.isnan(pval)] = 0

        pvalup = sio.loadmat(OBS_DIR + region + '-up_log-pval')['matrix'][0] #raw hypergeom p value
        pvalup[np.argwhere(pvalup == 0)] = 1
        pvalup[np.isinf(pvalup)] = 0
        pvalup[np.isnan(pvalup)] = 0

        #Calculate the overlap genes with annot_df, and only save the information on those genes.
        down_i = np.where(np.array(fc) < (1 - fc_cutoff))[0]
        down_idx = np.array(list(set(down_i).intersection(set(np.where(pval < pval_cutoff)[0]))))
        up_i = np.where(np.array(fc) > (1 - fc_cutoff))[0]
        up_idx = np.array(list(set(up_i).intersection(set(np.where(pvalup < pval_cutoff)[0]))))
        #Remove the duplicate value of gene index, return a list with unique gene names. 
        #Unique_elements is sorted by the gene names, unique_indices return the original index from gene_seq.
        unique_elements, unique_indices = np.unique(gene_seq, return_index=True)
        down_keep_genes = list((set(annot_df['gene_names']).intersection(set(gene_seq[down_idx]))).intersection(set(local_gene)))
        down_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, down_keep_genes))[0]]))
        up_keep_genes = list((set(annot_df['gene_names']).intersection(set(gene_seq[up_idx]))).intersection(set(local_gene)))
        up_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, up_keep_genes))[0]]))

        if (len(down_keep_genes_idx) == 0) & (len(up_keep_genes_idx) == 0):
            logger.info(f'No DE genes within local analysis windown. ')
            continue
        
        #Calculate the adjusted pval with closet cell number distribution
        chosen_dist = Num[np.argmin([np.absolute(n - cell_num) for n in Num])]
        
        #Load the file
        down_A = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-A.npy'%(str(chosen_dist)))
        down_B = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-B.npy'%(str(chosen_dist)))
        down_C = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-C.npy'%(str(chosen_dist)))
        up_A = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-A.npy'%(str(chosen_dist)))
        up_B = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-B.npy'%(str(chosen_dist)))
        up_C = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-C.npy'%(str(chosen_dist)))
        cpm_mean = np.load(DISTRI_DIR + 'Cpm_mean-%s.npy'%(str(chosen_dist)))
        fc_rand = cpm/cpm_mean
       
        #Calculate adjusted p value of gamma adj
        down_padj_list = []
        for i in down_keep_genes_idx:
            down_padj = scipy.stats.gamma.logsf(-pval[i], down_A[i], down_B[i], down_C[i])
            down_padj_list.append(down_padj)
        hits_fc = fc[down_keep_genes_idx]
        hits_fc_rand = fc_rand[down_keep_genes_idx]

        up_padj_list = []
        for i in up_keep_genes_idx:
            up_padj = scipy.stats.gamma.logsf(-pvalup[i], up_A[i], up_B[i], up_C[i])
            up_padj_list.append(up_padj)
        uphits_fc = fc[up_keep_genes_idx]
        uphits_fc_rand = fc_rand[up_keep_genes_idx]

        #Load pvalue matrix and calculate empirical p-value
        rand_down_file = sio.loadmat(DISTRI_DIR + '%s-down_log-pval'%(str(chosen_dist)))
        rand_down_matrix = []
        rand_down_matrix = sp_sparse.vstack(rand_down_file['matrix'])
        iter_num, gene_num = rand_down_matrix.shape
        emp_pval = np.sum(np.asarray(rand_down_matrix.tocsr()[:, down_keep_genes_idx].todense()) < pval[down_keep_genes_idx], axis=0) / iter_num
        
        rand_up_file = sio.loadmat(DISTRI_DIR + '%s-up_log-pval'%(str(chosen_dist)))
        rand_up_matrix = []
        rand_up_matrix = sp_sparse.vstack(rand_up_file['matrix'])
        iter_num, gene_num = rand_up_matrix.shape
        emp_pvalup = np.sum(np.asarray(rand_up_matrix.tocsr()[:, up_keep_genes_idx].todense()) < pval[up_keep_genes_idx], axis=0) / iter_num
        
        #Save to csv file 
        local_gene_series = annot_df[annot_df['gene_names'].isin(gene_seq[down_keep_genes_idx])].set_index('idx').sort_index()
        local_gene_series['region'] = region
        dist_list = []
        for i in local_gene_series['pos']:
            dist = get_distance(region, i)
            dist_list.append(dist)
        local_gene_series['distance'] = dist_list
        local_gene_series['num_cell'] = cell_num
        local_gene_series['bin'] = chosen_dist
        local_gene_series['log(pval)-hypergeom'] = pval[down_keep_genes_idx] 
        local_gene_series['fc'] = hits_fc
        local_gene_series['Significance_score'] = down_padj_list
        local_gene_series['fc_by_rand_dist_cpm'] = hits_fc_rand
        local_gene_series['pval-empirical'] = emp_pval
        local_gene_series['cpm_perturb'] = cpm[down_keep_genes_idx]
        local_gene_series['cpm_bg'] = cpm_mean[down_keep_genes_idx]

        local_gene_series = annot_df[annot_df['gene_names'].isin(gene_seq[up_keep_genes_idx])].set_index('idx').sort_index()
        local_gene_series['region'] = region
        dist_list = []
        for i in local_gene_series['pos']:
            dist = get_distance(region, i)
            dist_list.append(dist)
        local_gene_series['distance'] = dist_list
        local_gene_series['num_cell'] = cell_num
        local_gene_series['bin'] = chosen_dist
        local_gene_series['log(pval)-hypergeom'] = pvalup[up_keep_genes_idx] 
        local_gene_series['fc'] = uphits_fc
        local_gene_series['Significance_score'] = up_padj_list
        local_gene_series['fc_by_rand_dist_cpm'] = uphits_fc_rand
        local_gene_series['pval-empirical'] = emp_pvalup
        local_gene_series['cpm_perturb'] = cpm[up_keep_genes_idx]
        local_gene_series['cpm_bg'] = cpm_mean[up_keep_genes_idx]

        local_gene_df = pd.concat([local_gene_df, local_gene_series], ignore_index=True)
        local_gene_df = local_gene_df.reindex(columns=df_column_list)

    if OUTPUT_DF.endswith('.csv'):
        local_gene_df[~local_gene_df['gene_names'].isin(duplicate_elements)].to_csv(OUTPUT_DF)
    else:
        local_gene_df[~local_gene_df['gene_names'].isin(duplicate_elements)].to_csv(OUTPUT_DF + '.csv')
        
    logger.info('Job is done.')

if __name__ == '__main__':
    pass
