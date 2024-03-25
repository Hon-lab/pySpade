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

from pySpade.utils import get_logger, read_annot_df, read_sgrna_dict, load_data

logger = get_logger(logger_name=__name__)

def global_analysis(FILE_DIR,
                    OBS_DIR,
                    SGRNA_DICT,
                    DISTRI_DIR,
                    OUTPUT_DF):
    
    logger.info('Loading files.')
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

    ##Load sgRNA dict: All regions 
    sgrna_dict  = read_sgrna_dict(SGRNA_DICT)

    #Save the randomized bin with list 
    Num = list(set([int(i.split('/')[-1].split('.')[0].split('-')[-2]) for i in glob.glob(DISTRI_DIR+ '/Down_dist_gamma-*')]))
    logger.info(f'{len(Num)} sets of background distribution for global analysis.')

    #Calcualte pval
    df_column_list = [
        'idx', 'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 
        'region', 'num_cell', 'bin',
        'log(pval)-hypergeom', 'fc', 'Significance_score', 'fc_by_rand_dist_cpm', 'pval-empirical', 'cpm_perturb', 'cpm_bg']
    global_hits_df = pd.DataFrame(columns=df_column_list)

    #Generate bin-region dictionary 
    bin_dict = defaultdict(list)
    for region in list(sgrna_dict.keys()):
        fc_files = glob.glob(OBS_DIR + region + '*-foldchange')
        if len(fc_files) == 0:
            logger.critical('Cannot find foldchange file: ' + str(region))
            continue
        fc_file = fc_files[0]
        cell_num = int(fc_file.split('/')[-1].split('-')[-2])
        chosen_bin = Num[np.argmin([np.absolute(n-cell_num) for n in Num])]
        bin_dict[chosen_bin].append(region)
    logger.info('Finish generating bin-region dictionary')
    
    for b in Num:
        logger.info(f'Processing bin: {b}')
        logger.info(f'Num of regions in this bin: {len(bin_dict[b])}')

        #Load the background files
        down_A = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-A.npy'%(str(b)))
        down_B = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-B.npy'%(str(b)))
        down_C = np.load(DISTRI_DIR + 'Down_dist_gamma-%s-C.npy'%(str(b)))
        up_A = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-A.npy'%(str(b)))
        up_B = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-B.npy'%(str(b)))
        up_C = np.load(DISTRI_DIR + 'Up_dist_gamma-%s-C.npy'%(str(b)))
        cpm_mean = np.load(DISTRI_DIR + 'Cpm_mean-%s.npy'%(str(b)))

        for region in bin_dict[b]:
            pval_list_up, pval_list_down, cpm, fc, cell_num = load_data(OBS_DIR, region)
            fc_cpm = (cpm + 0.01)/(cpm_mean + 0.01)

            #preprocess pval list 
            pval_list_up[np.argwhere(pval_list_up == 0)] = 1
            pval_list_up[np.isinf(pval_list_up)] = 0
            pval_list_up[np.isnan(pval_list_up)] = 0
            pval_list_down[np.argwhere(pval_list_down == 0)] = 1
            pval_list_down[np.isinf(pval_list_down)] = 0
            pval_list_down[np.isnan(pval_list_down)] = 0

            #split up-regulation and down-regulation
            pval_cutoff = 0
            fc_cutoff = 0
            up_i = np.where(np.array(fc_cpm) > (1 + fc_cutoff))[0]
            up_idx = np.array(list(set(up_i).intersection(set(np.where(pval_list_up < pval_cutoff)[0]))))
            down_i = np.where(np.array(fc_cpm) < (1 - fc_cutoff))[0]
            down_idx = np.array(list(set(down_i).intersection(set(np.where(pval_list_down < pval_cutoff)[0]))))

            #Calculate the overlap genes with annot_df, and only save the information on those genes.
            unique_elements, unique_indices = np.unique(gene_seq, return_index=True)
            up_keep_genes = list(set(annot_df['gene_names']).intersection(set(gene_seq[up_idx]))) #up-regulated genes found in both annotation file and transcriptome df
            up_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, up_keep_genes))[0]]))
            down_keep_genes = list(set(annot_df['gene_names']).intersection(set(gene_seq[down_idx])))
            down_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, down_keep_genes))[0]]))

            #Calculate p-value adj for gamma distribution
            down_padj_list = []
            for i in down_keep_genes_idx:
                down_padj = scipy.stats.gamma.logsf(-pval_list_down[i], down_A[i], down_B[i], down_C[i])
                down_padj_list.append(down_padj)
            down_hit_fc_list = fc_cpm[down_keep_genes_idx]
            down_cpm = cpm[down_keep_genes_idx]
            down_cpm_bg = cpm_mean[down_keep_genes_idx]

            up_padj_list = []
            for i in up_keep_genes_idx:
                up_padj = scipy.stats.gamma.logsf(-pval_list_up[i], up_A[i], up_B[i], up_C[i])
                up_padj_list.append(up_padj)
            up_hit_fc_list = fc_cpm[up_keep_genes_idx]
            up_cpm = cpm[up_keep_genes_idx]
            up_cpm_bg = cpm_mean[up_keep_genes_idx]

            #Load background distribution and calculate emprical p-value
            rand_down_file = sio.loadmat(DISTRI_DIR + '%s-down_log-pval'%(str(b)))
            rand_down_matrix = []
            rand_down_matrix = sp_sparse.vstack(rand_down_file['matrix'])
            iter_num, gene_num = rand_down_matrix.shape
            emp_pval_down = np.sum(np.asarray(rand_down_matrix.tocsr()[:, down_keep_genes_idx].todense()) < pval_list_down[down_keep_genes_idx], axis=0) / iter_num

            rand_up_file = sio.loadmat(DISTRI_DIR + '%s-up_log-pval'%(str(b)))
            rand_up_matrix = []
            rand_up_matrix = sp_sparse.vstack(rand_up_file['matrix'])
            iter_num, gene_num = rand_up_matrix.shape
            emp_pval_up = np.sum(np.asarray(rand_up_matrix.tocsr()[:, up_keep_genes_idx].todense()) < pval_list_up[up_keep_genes_idx], axis=0) / iter_num

            #save to csv file: down-regulation gene 
            global_gene_series = annot_df[annot_df['gene_names'].isin(gene_seq[down_keep_genes_idx])].set_index('idx').sort_index()
            global_gene_series['region'] = region
            global_gene_series['num_cell'] = cell_num
            global_gene_series['bin'] = b
            global_gene_series['log(pval)-hypergeom'] = pval_list_down[down_keep_genes_idx]
            global_gene_series['fc'] = fc[down_keep_genes_idx]
            global_gene_series['Significance_score'] = down_padj_list
            global_gene_series['fc_by_rand_dist_cpm'] = down_hit_fc_list
            global_gene_series['pval-empirical'] = emp_pval_down
            global_gene_series['cpm_perturb'] = down_cpm
            global_gene_series['cpm_bg'] = down_cpm_bg
            global_gene_series['idx'] = global_gene_series.index
            #global_hits_df = global_hits_df.append(global_gene_series)
            global_hits_df = pd.concat([global_hits_df, global_gene_series], ignore_index=True)

            #save to csv file: up-regulation gene 
            global_gene_series = annot_df[annot_df['gene_names'].isin(gene_seq[up_keep_genes_idx])].set_index('idx').sort_index()
            global_gene_series['region'] = region
            global_gene_series['num_cell'] = cell_num
            global_gene_series['bin'] = b
            global_gene_series['log(pval)-hypergeom'] = pval_list_up[up_keep_genes_idx]
            global_gene_series['fc'] = fc[up_keep_genes_idx]
            global_gene_series['Significance_score'] = up_padj_list
            global_gene_series['fc_by_rand_dist_cpm'] = up_hit_fc_list
            global_gene_series['pval-empirical'] = emp_pval_up
            global_gene_series['cpm_perturb'] = up_cpm
            global_gene_series['cpm_bg'] = up_cpm_bg
            global_gene_series['idx'] = global_gene_series.index
            #global_hits_df = global_hits_df.append(global_gene_series)
            global_hits_df = pd.concat([global_hits_df, global_gene_series], ignore_index=True)
            
            logger.info(f'Finish region: {region}')
    
    global_hits_df = global_hits_df.reindex(columns=df_column_list)
    rem_dup_global_hits_df = global_hits_df[~global_hits_df['gene_names'].isin(duplicate_elements)]

    if OUTPUT_DF.endswith('.csv'):
        rem_dup_global_hits_df[rem_dup_global_hits_df['Significance_score'] < pval_cutoff].to_csv(OUTPUT_DF)
    else:
        rem_dup_global_hits_df[rem_dup_global_hits_df['Significance_score'] < pval_cutoff].to_csv(OUTPUT_DF + '.csv')

    logger.info('Job is done.')
    
if __name__ == '__main__':
    pass
