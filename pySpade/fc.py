#fc.py

import os
import re
import sys
import collections
import argparse
import tables
import itertools
import numba
import glob

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

from pySpade.utils import get_logger, read_sgrna_dict, find_single_sgrna_cells, find_all_sgrna_cells, set_axis_style, sc_exp_box_plot

logger = get_logger(logger_name=__name__)

def Calculate_fc(TRANSCRIPTOME_DIR, 
                 SGRNA_DICT, 
                 TARGET_FILE, 
                 OUTPUT_FOLDER):

    TRANSCRIPTOME_FILE = glob.glob(TRANSCRIPTOME_DIR + 'Singlet_sub_df.*')
    if (TRANSCRIPTOME_FILE[0].endswith('pkl')) == True:
        sub_df = pd.read_pickle(TRANSCRIPTOME_FILE[0])
    elif (TRANSCRIPTOME_FILE[0].endswith('h5')) == True:
        sub_df = pd.read_hdf(TRANSCRIPTOME_FILE[0], 'df')

    total_seq_reads = sub_df.sum(axis=0)
    SGRNA_FILE = glob.glob(TRANSCRIPTOME_DIR + 'Singlet_sgRNA*')
    if (SGRNA_FILE[0].endswith('pkl')):
        sgrna_df_bool = pd.read_pickle(SGRNA_FILE[0]) > 0 
    elif (SGRNA_FILE[0].endswith('h5')):
        sgrna_df_bool = pd.read_hdf(SGRNA_FILE[0], 'df') > 0

    if list(sub_df.columns) != list(sgrna_df_bool.columns):
        logger.critical('File format error, make sure transcriptome df columns are the same as sgrna df.')

    #load sgRNA dictionary
    sgrna_dict  = read_sgrna_dict(SGRNA_DICT)
    
    #load regions and their target genes
    query_region_list = []
    query_gene_list = []
    with open(TARGET_FILE) as ft:
        for line in ft:
            region, gene = line.strip().split('\t')
            query_region_list.append(region)
            query_gene_list.append(gene)

    logger.info('Finished loading data.')

    #prepare cpm array for each query gene
    All_cpm_array = []
    for gene in query_gene_list:
        cpm_array = sub_df.loc[gene ,:] / total_seq_reads * 1000000
        All_cpm_array.append(cpm_array)

    #calculate repression efficiency for single sgRNA
    All_FC_list = []
    All_sgRNA = []
    All_num_cell = []
    All_pert_cpm_list = []
    All_bg_cpm_list = []
    region_fc_list = []
    region_cell_list = []
    region_pval_list = []
    region_pert_cpm_list = []
    region_bg_cpm_list = []
    All_pval_list = []
    for i , target_gene, non_zero_array in zip(query_region_list, query_gene_list, All_cpm_array):
        logger.info(f'Calculate gene: {target_gene} in region: {i}')
        
        #Find if any sgRNA is missing from the sgRNA df
        missing_sgrna_list = list(set(sgrna_dict[i]) - set(sgrna_df_bool.index))
        if len(missing_sgrna_list) == len(sgrna_dict[i]):
            logger.info('Missing all the sgRNA in this region.')
            continue
        if len(missing_sgrna_list) > 0:
            logger.info('Missing ' + str(len(missing_sgrna_list)) + ' sgrna.')
            for j in missing_sgrna_list:
                logger.info('Missing sgrna: ' + str(j))

        #remove the missing sgrna and continue DE analysis 
        if len(missing_sgrna_list) > 0:
            remain_sgrna_list = list(set(sgrna_dict[i]) - (set(missing_sgrna_list)))
            sgrna_dict.update({i:remain_sgrna_list})

        sgrna_idx = find_all_sgrna_cells(sgrna_df_bool, sgrna_dict[i])
        other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
        region_cell_num = len(sgrna_idx)
        region_cell_list.append(region_cell_num)

        sgrna_cells_expression = non_zero_array[sgrna_idx]
        other_cells_expression = non_zero_array[other_idx]
        fc = (np.mean(sgrna_cells_expression) + 0.01) / (np.mean(other_cells_expression) + 0.01)
        pval = stats.ttest_ind(sgrna_cells_expression, other_cells_expression, equal_var=False).pvalue
        region_fc_list.append(fc)
        region_pval_list.append(pval)
        region_pert_cpm_list.append(str(np.mean(sgrna_cells_expression)))
        region_bg_cpm_list.append(str(np.mean(other_cells_expression)))

        #Box plot for each region and the query gene 
        save_name = i.replace(':', '_')
        box_plot_file = OUTPUT_FOLDER + 'Box_plot-' + save_name + '-' + target_gene + '-10perc.tiff'
        sc_exp_box_plot(sgrna_cells_expression, other_cells_expression, target_gene, i, box_plot_file)

        sgrna_seq = []
        fc_list = []
        num_cell_list = []
        pert_cpm_list = []
        bg_cpm_list = []
        pval_list = []
        #Now we calculate fold change for individual sgrna in the perturbation
        for sgrna in sgrna_dict[i]:
            try:
                sgrna_idx = find_single_sgrna_cells(sgrna_df_bool, sgrna)
                unperturb_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
                perturb_cpm = np.mean(non_zero_array[sgrna_idx]) 
                bg_cpm = np.mean(non_zero_array[unperturb_idx])
                fc = (np.mean(non_zero_array[sgrna_idx]) + 0.01) / (np.mean(non_zero_array[unperturb_idx]) + 0.01)
                pval = stats.ttest_ind(non_zero_array[sgrna_idx], non_zero_array[unperturb_idx], equal_var=False).pvalue
                num_cell = str(len(sgrna_idx))
            except:
                fc = 1
                num_cell = '0'

            pert_cpm_list.append(perturb_cpm) 
            bg_cpm_list.append(bg_cpm)   
            fc_list.append(fc)
            pval_list.append(pval)
            sgrna_seq.append(sgrna)
            num_cell_list.append(num_cell)

        All_FC_list.append(fc_list)
        All_sgRNA.append(sgrna_seq)
        All_num_cell.append(num_cell_list)
        All_pval_list.append(pval_list)
        All_pert_cpm_list.append(pert_cpm_list)
        All_bg_cpm_list.append(bg_cpm_list)

    #write to output file
    output_file=open(OUTPUT_FOLDER + 'fold_change.txt', 'w')
    output_file.write('Region (sgRNA sequence)' + '\t' + 'Target gene' + '\t' + 'Number of cells' + '\t' + 'Fold change' + '\t' + 'P value' + '\t' + 'Perturb cpm' + '\t' + 'Background cpm' + '\n')

    for idx in np.arange(len(query_region_list)):
        i = query_region_list[idx]
        target_gene = query_gene_list[idx]
        output_file.write(str(i) + '\t' + str(target_gene) + '\t' + str(region_cell_list[idx]) + '\t' + str(region_fc_list[idx]) + '\t' + str(region_pval_list[idx]) + '\t' + str(region_pert_cpm_list[idx]) + '\t' + str(region_bg_cpm_list[idx])+ '\n')
        for sg in All_sgRNA[idx]:
            sec_idx = list(All_sgRNA[idx]).index(sg)
            num_cell = All_num_cell[idx][sec_idx]
            repression = All_FC_list[idx][sec_idx]
            pvalue = All_pval_list[idx][sec_idx]
            pert_cpm = All_pert_cpm_list[idx][sec_idx]
            sgbg_cpm = All_bg_cpm_list[idx][sec_idx]
            output_file.write(str(sg) + '\t' + str(target_gene)+ '\t' + str(num_cell) + '\t' + str(repression) + '\t' + str(pvalue) + '\t' + str(pert_cpm) + '\t' + str(sgbg_cpm) + '\n')
    
    output_file.close()

if __name__ == '__main__':
    pass
