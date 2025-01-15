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
import warnings
warnings.filterwarnings('ignore')

from pySpade.utils import get_logger, read_annot_df, read_sgrna_dict, load_data, get_num_processes

logger = get_logger(logger_name=__name__)

def global_analysis(FILE_DIR,
                    OBS_DIR,
                    SGRNA_DICT,
                    DISTRI_DIR,
                    OUTPUT):
    
    logger.info('Loading files.')
    #read the gene sequence file 
    gene_seq = np.load(FILE_DIR + 'Trans_genome_seq.npy', allow_pickle=True)
    if len(gene_seq) != len(set(gene_seq)):
        logger.critical('Duplication of mapping genes. Duplicates are removed in the analysis.')
    unique_elements, counts = np.unique(gene_seq, return_counts=True)
    duplicate_elements = unique_elements[counts > 1]
    #read the plotting annotation
    annot_df_dup = read_annot_df()
    #Add the genes that are missing in the annot_df 
    new_genes = list(set(gene_seq).difference(set(annot_df_dup['gene_names'])))
    start_idx = annot_df["idx"].max() + 1
    new_data = {
    "idx": range(start_idx, start_idx + len(new_genes)),
    "gene_names": new_genes,
    "chromosome": ['NA'] * len(new_genes),  # Placeholder values
    "pos": ['NA'] * len(new_genes),      # Placeholder values
    "strand": ['+'] * len(new_genes),     # Placeholder values
    "color_idx": ['0'] * len(new_genes), # Placeholder values
    "chr_idx": ['24'] * len(new_genes),   # Placeholder values
    }
    new_df = pd.DataFrame(new_data)
    result_df = pd.concat([annot_df_dup, new_df], ignore_index=True)
    #There are many non-coding genes duplication in the annot_df, only keep one.
    annot_df = result_df.drop_duplicates(subset='gene_names', keep='first')

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
        fc_files = glob.glob(OBS_DIR + region + '-' + '*-foldchange')
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
            pval_list_up = pval_list_up[:len(cpm_mean)]
            pval_list_down = pval_list_down[:len(cpm_mean)]
            cpm = cpm[:len(cpm_mean)]
            fc_cpm = (cpm + 0.01)/(cpm_mean + 0.01) #modify back!!
            
            if cell_num < 20:
                logger.info('Not enough cells (' + str(cell_num) + ') for analysis.')
                continue

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
            if (len(down_idx) == 0) & (len(up_idx) == 0):
                logger.info(f'No DE genes within local analysis windown. ')
                continue
            #Calculate the overlap genes with annot_df, and only save the information on those genes.
            unique_elements, unique_indices = np.unique(gene_seq, return_index=True)
            if len(up_idx) != 0:
                up_keep_genes = list(set(annot_df['gene_names']).intersection(set(gene_seq[up_idx]))) #up-regulated genes found in both annotation file and transcriptome df
                up_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, up_keep_genes))[0]]))
            if len(down_idx) != 0:
                down_keep_genes = list(set(annot_df['gene_names']).intersection(set(gene_seq[down_idx])))
                down_keep_genes_idx = sorted(list(unique_indices[np.where(np.isin(unique_elements, down_keep_genes))[0]]))
            if (len(down_keep_genes_idx) == 0) & (len(up_keep_genes_idx) == 0):
                logger.info(f'No DE genes within local analysis windown. ')
                continue
            
            #Calculate p-value adj for gamma distribution
            down_padj_list = []
            num_process = get_num_processes()
            with Pool(processes=num_process) as p:
                for down_adj in p.starmap(scipy.stats.gamma.logsf, zip(
                    -pval_list_down[down_keep_genes_idx], 
                    down_A[down_keep_genes_idx], 
                    down_B[down_keep_genes_idx], 
                    down_C[down_keep_genes_idx])
                ):
                    down_padj_list.append(down_adj)
            down_hit_fc_list = fc_cpm[down_keep_genes_idx]
            down_cpm = cpm[down_keep_genes_idx]
            down_cpm_bg = cpm_mean[down_keep_genes_idx]

            up_padj_list = []
            with Pool(processes=num_process) as p:
                for up_adj in p.starmap(scipy.stats.gamma.logsf, zip(
                    -pval_list_up[up_keep_genes_idx], 
                    up_A[up_keep_genes_idx], 
                    up_B[up_keep_genes_idx], 
                    up_C[up_keep_genes_idx])
                ):
                    up_padj_list.append(up_adj)
            up_hit_fc_list = fc_cpm[up_keep_genes_idx]
            up_cpm = cpm[up_keep_genes_idx]
            up_cpm_bg = cpm_mean[up_keep_genes_idx]

            #Load background distribution and calculate emprical p-value
            rand_down_file = sio.loadmat(DISTRI_DIR + '%s-down_log-pval'%(str(b)))
            rand_down_matrix = []
            rand_down_matrix = sp_sparse.vstack([rand_down_file['matrix']])
            iter_num, gene_num = rand_down_matrix.shape
            emp_pval_down = np.sum(np.asarray(rand_down_matrix.tocsr()[:, down_keep_genes_idx].todense()) < pval_list_down[down_keep_genes_idx], axis=0) / iter_num

            rand_up_file = sio.loadmat(DISTRI_DIR + '%s-up_log-pval'%(str(b)))
            rand_up_matrix = []
            rand_up_matrix = sp_sparse.vstack([rand_up_file['matrix']])
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
    
    #save as csv
    rem_dup_global_hits_df[rem_dup_global_hits_df['Significance_score'] < pval_cutoff].to_csv(OUTPUT + '/unfiltered_global_df.csv')

    logger.info('Job is done.')
    
if __name__ == '__main__':
    pass
