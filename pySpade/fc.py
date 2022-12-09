#fc.py

from asyncio.log import logger
import os
import re
import sys
import collections
import argparse
import tables
import itertools
import numba

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix

from pySpade.utils import get_logger, read_sgrna_dict, read_annot_df, find_single_sgrna_cells, find_all_sgrna_cells

logger = get_logger(logger_name=__name__)

def Calculate_fc(TRANSCRIPTOME_DF, 
                 SGRNA, 
                 SGRNA_DICT, 
                 TARGET_FILE, 
                 OUTPUT_FILE):
    
    if (TRANSCRIPTOME_DF.endswith('pkl')):
        sub_df = pd.read_pickle(TRANSCRIPTOME_DF)
    elif (TRANSCRIPTOME_DF.endswith('h5')):
        sub_df = pd.read_hdf(TRANSCRIPTOME_DF, 'df')

    total_seq_reads = sub_df.sum(axis=0)

    if (SGRNA.endswith('pkl')):
        sgrna_df_bool = pd.read_pickle(SGRNA) > 0 
    elif (SGRNA.endswith('h5')):
        sgrna_df_bool = pd.read_hdf(SGRNA, 'df') > 0

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
    
    #load annotation df 
    annot_df = read_annot_df()
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
    region_fc_list = []
    all_cell_idx = list(np.argwhere(sub_df.columns)[0])
    for i , target_gene, non_zero_array in zip(query_region_list, query_gene_list, All_cpm_array):
        logger.info(f'Calculate gene: {target_gene} in region: {i}')        
        sgrna_idx = find_all_sgrna_cells(sgrna_df_bool, sgrna_dict[i])        
        other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))

        sgrna_cells_expression = non_zero_array[sgrna_idx]
        other_cells_expression = non_zero_array[other_idx]
        fc = (np.mean(sgrna_cells_expression) + 0.01) / (np.mean(other_cells_expression) + 0.01)
        region_fc_list.append(fc)

        sgrna_seq = []
        fc_list = []
        num_cell_list = []
        
        #Now we calculate fold change for individual sgrna in the perturbation
        for sgrna in sgrna_dict[i]:
            try:
                sgrna_idx = find_single_sgrna_cells(sgrna_df_bool, sgrna)
                unperturb_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
                fc = (np.mean(non_zero_array[sgrna_idx]) + 0.01) / (np.mean(non_zero_array[unperturb_idx]) + 0.01)
                num_cell = str(len(sgrna_idx))
            except:
                fc = 1
                num_cell = '0'
                
            fc_list.append(fc)
            sgrna_seq.append(sgrna)
            num_cell_list.append(num_cell)

        All_FC_list.append(fc_list)
        All_sgRNA.append(sgrna_seq)
        All_num_cell.append(num_cell_list)
    
    #write to output file
    output_file=open(OUTPUT_FILE, 'w')
    output_file.write('Region: target gene: region fold change' + '\n')
    output_file.write('\t' + 'sgRNA sequence: number of cells: fold change' + '\n')

    for idx in np.arange(len(query_region_list)):
        i = query_region_list[idx]
        target_gene = query_gene_list[idx]
        output_file.write(str(i) + ': ' + str(target_gene) + ': ' + str(region_fc_list[idx]) + '\n')
        for sg in All_sgRNA[idx]:
            sec_idx = list(All_sgRNA[idx]).index(sg)
            num_cell = All_num_cell[idx][sec_idx]
            repression = All_FC_list[idx][sec_idx]
            output_file.write('\t' + str(sg) + ': ' + str(num_cell) + ': ' + str(repression) + '\n')
    
    output_file.close()

if __name__ == '__main__':
    pass
