#utils.py
import logging
import numpy as np
import pandas as pd
from multiprocessing import Pool
import scipy.stats as stats
import scipy.sparse as sp_sparse
import collections
from collections import defaultdict
import itertools 
from importlib import resources
from scipy import sparse, io
from scipy.sparse import csr_matrix
import glob
import re
import tables

def get_logger(logger_name, log_file=False):
    """Creates a custom logger."""

    FORMATTER = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logger = logging.getLogger(logger_name)
    logger.setLevel(level=logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(FORMATTER))
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(filename='pySpade.log')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(FORMATTER))
        logger.addHandler(file_handler)

    return logger


def process_feature_df(mtx, barcodes, feature):
    gene_name = []
    for i in feature:
        gene_name.append(i[1])
    transcriptome_df = pd.DataFrame.sparse.from_spmatrix(data = mtx.tocsr(),
                                                         columns=barcodes, index=gene_name)
    return transcriptome_df

CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode('utf-8')
            feature_ref[key] = getattr(feature_group, key).read()
         
        return CountMatrix(feature_ref, barcodes, matrix)

def turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = np.sum(sgRNA_count > 0)

    #calculate the turning point by using the max derivative
    turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
    
    return(sgRNA_count.iloc[turning_point])

def filter_umi (df, copy=False):
    df = df.copy() if copy else df
    feature_cutoff = [[turn_point(i, df)] for i in list(df.index)]
    
    for i in range(0, len(feature_cutoff)):
        ZERO_HTO = df.iloc[i, :].loc[df.iloc[i, :] <= feature_cutoff[i][0]].index
        df.at[df.index[i], ZERO_HTO] = 0
    return df, feature_cutoff

def get_sgrna_per_cell(df, return_mean=True, return_median=True):
    CorrSgrnaPerCell = np.sum(df > 0, 0)
    sgrna_mean = np.mean(CorrSgrnaPerCell)
    sgrna_median = np.median(CorrSgrnaPerCell)
    
    if return_mean & return_median:
        return CorrSgrnaPerCell, sgrna_mean, sgrna_median
    elif return_mean:
        return CorrSgrnaPerCell, sgrna_mean
    elif return_median:
        return CorrSgrnaPerCell, sgrna_median
    else:
        return CorrSgrnaPerCell


def read_sgrna_dict(SGRNA_FILE):
    sgrna_dict  = {}
    with open(SGRNA_FILE) as f:
        for line in f:
            region_id, sgrna_string = line.strip().split("\t")
            sgrnas = sgrna_string.split(";")
            sgrna_dict.update({region_id : sgrnas})
    return(sgrna_dict)

def read_annot_df() -> pd.DataFrame:
    with resources.path('pySpade', 'plot_annotation.txt') as df:
        annot_df = pd.read_csv(df,
                   header=None,
                   sep='\t',
                   names=['idx', 'gene_names', 'chromosome', 'pos', 'strand', 'color_idx', 'chr_idx'])
        return annot_df

#be sure that when calculating find_sgrna_cells, the column sequences are the same between transcriptome df and sgrna df
#modified, specific for single sgRNA
def find_single_sgrna_cells(sgRNA_df, sgrna):
    cell_bc = sgRNA_df.loc[:, (sgRNA_df.loc[sgrna] > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(sgRNA_df.columns == i)
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def find_all_sgrna_cells(sgRNA_df, sgRNA_dict):
    cell_bc = sgRNA_df.loc[:, (sgRNA_df.loc[sgRNA_dict].sum() > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(sgRNA_df.columns == i)
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def cpm_normalization(trans_df):
    [g,c] = trans_df.shape
    cpm_matrix = np.zeros((g,c))
    cpm_matrix = np.array(trans_df / trans_df.sum(axis = 0) * 1e6)
    return cpm_matrix

def metacelll_normalization(trans_df):
    [g,c] = trans_df.shape
    cpm_matrix = np.zeros((g,c))
    uniq_id = set()
    for x in trans_df.columns:
        uniq_id.add(x[-1])
    for lib in sorted(uniq_id):
        index = [i for i,e in enumerate(trans_df.columns) if e.split('-')[1] == lib]
        one_cell_cpm = np.array((trans_df.iloc[:,index].sum(axis = 1) + 1)\
                        / np.sum(trans_df.iloc[:,index].values) * 1e6).flatten()
        for i in index:
            cell_cpm = trans_df.iloc[:,i] / np.sum(trans_df.iloc[:,i].values) * 1e6
            cpm_matrix[:,i] = cell_cpm / one_cell_cpm  

def hypergeo_test(non_zero_array, sgrna_idx, i):
    #find indecies of cells in which expression of given gene is
    #equal or less than the median of this gene in the whole population
    median_cell_idx  = np.argwhere(non_zero_array <= np.median(non_zero_array))

    #find the same cells subset in the cells with a given sgRNA
    overlap_cell_idx = np.intersect1d(median_cell_idx, sgrna_idx)
    
    #calculate the median fold change
    other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
    
    fc = (np.mean(non_zero_array[sgrna_idx]) + 0.01) / (np.mean(non_zero_array[other_idx]) + 0.01)
    cpm = np.mean(non_zero_array[sgrna_idx])
    
    #perform hypergeometric test, get the upper tail
    k = len(overlap_cell_idx)
    M = len(non_zero_array)
    n = len(median_cell_idx)
    N = len(sgrna_idx)
    try:
        pval_up = stats.hypergeom.logcdf(k, M, n, N).item()
    except:
        pval_up = float('nan')
        
    try:
        pval_down = stats.hypergeom.logsf(k, M, n, N).item()
    except:
        pval_down = float('nan')
    
    return pval_down, pval_up, fc, cpm

def perform_DE(sgrna_idx, input_array, idx, num_processes, pval_list_down, pval_list_up, fc_list, cpm_list):
    nonzero_pval_list_up = []
    nonzero_pval_list_down = []
    nonzero_fc_list = []
    nonzero_cpm_list = []
    
    with Pool(processes=num_processes) as p:
        for pval_down, pval_up, fc, cpm in p.starmap(hypergeo_test, zip(
                input_array,
                itertools.repeat(sgrna_idx),
                idx)
        ):
            nonzero_pval_list_down.append(pval_down)
            nonzero_pval_list_up.append(pval_up)
            nonzero_fc_list.append(fc)
            nonzero_cpm_list.append(cpm)
    for i in idx:
        pval_list_up[i] = nonzero_pval_list_up.pop(0)
        pval_list_down[i] = nonzero_pval_list_down.pop(0)
        fc_list[i] = nonzero_fc_list.pop(0)
        cpm_list[i] = nonzero_cpm_list.pop(0)
    return len(sgrna_idx), pval_list_up, pval_list_down, fc_list, cpm_list

def load_data(data_dir, region):
    up_pval_file   = data_dir + region + '-up_log-pval'
    down_pval_file = data_dir + region + '-down_log-pval'
    cpm_file = data_dir + region + '-cpm'
    fc_files = glob.glob(data_dir + region + '*-foldchange')
    if len(fc_files) == 1:
        fc_file = fc_files[0]
    else:
        
        numbers = np.array([int(float(i.split('/')[-1].split('-')[2])) for i in fc_files])
        chosen_num = np.max(numbers)
        fc_file = data_dir + region + '-' + str(chosen_num) + '*-foldchange'
    
    pval_list_up = io.loadmat(up_pval_file)['matrix'][0]
    pval_list_down = io.loadmat(down_pval_file)['matrix'][0]
    cpm = io.loadmat(cpm_file)['matrix'][0]
    fc = io.loadmat(fc_file)['matrix'][0]
    
    num_sgrna_cell = int(fc_file.split('/')[-1].split('-')[2])
    return pval_list_up, pval_list_down, cpm, fc, num_sgrna_cell

def get_neighbor_genes(region, query_range, annot_df):
    
    length_list = [
        0, 248956422,491149951,689445510,879660065,1061198324,
        1232004303,1391350276,1536488912,1674883629,1808681051,
        1943767673,2077042982,2191407310,2298451028,2400442217,
        2490780562,2574038003,2654411288,2713028904,2777473071,
        2824183054,2875001522,3031029399
    ]
    
    chr_order = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
    ]

    enh_chrom, left, right = re.split('[:|-]', region)
    position = int(left) + length_list[chr_order.index(enh_chrom)]
    gene_idx = annot_df.loc[(annot_df.pos < position + query_range) \
                            & (annot_df.pos > position - query_range)\
                            & (annot_df.chromosome == enh_chrom)].index.values
    return(gene_idx)

def get_distance(region, gene_pos, annot_df):
    length_list = [
        0, 248956422,491149951,689445510,879660065,1061198324,               
        1232004303,1391350276,1536488912,1674883629,1808681051,
        1943767673,2077042982,2191407310,2298451028,2400442217,
        2490780562,2574038003,2654411288,2713028904,2777473071,
        2824183054,2875001522,3031029399
    ]
    
    chr_order = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
    ]
    
    enh_chrom, left, right = re.split('[:|-]', region)
    dist =  int(gene_pos) - int(left) - length_list[chr_order.index(enh_chrom)]
    
    return np.absolute(dist)