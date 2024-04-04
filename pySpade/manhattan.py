#manhattan.py

import os
import re
import sys
import collections
import argparse
import tables
import itertools
import matplotlib
import numba
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix

from importlib import resources
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pySpade.utils import get_logger, outlier_plot, get_fc_range

logger = get_logger(logger_name=__name__)

def manhattan_plots(FILE_DIR, 
                    GLBOAL_HITS, 
                    OUTPUT,
                    CUTOFF_EXP=0.05,
                    CUTOFF_FC=0.2, 
                    CUTOFF_SIG=-5):
    
    length_list = [
        0, 248956422,491149951,689445510,879660065,1061198324,               
        1232004303,1391350276,1536488912,1674883629,1808681051,
        1943767673,2077042982,2191407310,2298451028,2400442217,
        2490780562,2574038003,2654411288,2713028904,2777473071,
        2824183054,2875001522,3031029399, 3088256814]
    
    chr_order = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    
    df_column_list = [
        'idx', 'gene_names', 'chromosome', 'pos', 'strand', 
        'color_idx', 'chr_idx', 
        'region', 'num_cell', 'bin',
        'log(pval)-hypergeom', 'fc', 'Significance_score', 'fc_by_rand_dist_cpm', 'pval-empirical', 'cpm_perturb', 'cpm_bg']
    
    #Load unfilter global hits df 
    global_df = pd.read_csv(GLBOAL_HITS)[df_column_list]
    #Load expressed genes 
    gene_seq = np.load(FILE_DIR + 'Trans_genome_seq.npy', allow_pickle=True)
    express_level = np.load(FILE_DIR + 'Perc_cell_expr.npy')
    express_idx = np.where(express_level > CUTOFF_EXP)[0]

    #filter df 
    express_df = global_df[global_df['gene_names'].isin(gene_seq[express_idx])\
                      &((global_df['fc_by_rand_dist_cpm'] > (1+CUTOFF_FC)) | (global_df['fc_by_rand_dist_cpm'] < (1-CUTOFF_FC)))\
                      &(global_df['Significance_score'] < CUTOFF_SIG)\
                      &(global_df['pval-empirical'] < 0.001)
                      &(global_df['log(pval)-hypergeom'] < -2)]
    #save filter df 
    express_df.to_csv(OUTPUT + '/filtered_df.csv')
    logger.info('filtered df size: ' + str(express_df.shape))

    #Manhattan plots
    for region in list(set(express_df['region'])):
        logger.info('Process: ' + str(region))
        express_subset_df = express_df[express_df['region'] == region]
        hits, _ = express_subset_df.shape
        if hits == 0:
            logger.info('No hits in this region.')
            continue
        
        #Prepare input for Manhattan plots
        up_idx = np.where(express_subset_df['fc_by_rand_dist_cpm'] > 1)[0]
        down_idx = np.where(express_subset_df['fc_by_rand_dist_cpm'] < 1)[0]
        plot_y_val = [0] * (len(express_subset_df['Significance_score'].values))

        for i in up_idx:
            plot_y_val[i] = -express_subset_df['Significance_score'].values[i]

        for i in down_idx:
            plot_y_val[i] = express_subset_df['Significance_score'].values[i]

        plot_y_val=np.array(plot_y_val) 
        plot_x_val=express_subset_df['pos'].values
        plot_y_val[np.isinf(plot_y_val)] = 0

        odd_idx = np.where(express_subset_df.color_idx == 0)
        even_idx= np.where(express_subset_df.color_idx == 1)

        fc = express_subset_df['fc_by_rand_dist_cpm'].values


        num_sgrna_cell = express_subset_df['num_cell'].values[0]
        if len(re.split(r'[:-]+', region)) == 3:
            _, left, right = re.split(r'[:-]+', region)
            if (_ in chr_order) == True: 
                enh_chrom = _
            else:
                logger.info('Not valid chromosome coordinates, plot without perturbation region line.')
                enh_chrom = 'chr1'
                left = 0
        else:
            logger.info('No chromosome coordinates info, plot without perturbation region line.')
            enh_chrom = 'chr1'
            left = 0
            
        up_cutoff = 1
        down_cutoff = 1


        #start plotting 
        fig = plt.figure(figsize=(14,8))
        gs = gridspec.GridSpec(nrows=1, ncols=11)

        #plot all genes
        ax0 = fig.add_subplot(gs[:, 0:9])

        ax0.scatter(plot_x_val[odd_idx],
                plot_y_val[odd_idx],
                s=1,
                color='#4d4d4d',
                marker='.')

        ax0.scatter(plot_x_val[even_idx],
                plot_y_val[even_idx],
                s=1,
                color='#e0e0e0',
                marker='.')

        ax0.set_title('%s (%d cells)'%(region, num_sgrna_cell),
                    fontsize=18)

        ax0.set_ylabel('Significance score',
                    fontsize = 15)

        #configurate the axis
        [ymin, ymax] = ax0.get_ylim()
        max_yval = max([np.absolute(ymin), np.absolute(ymax)])
        ax0.set_ylim([round(-max_yval-1),round(max_yval+1)])
        ax0.set_xlim([-1e8, length_list[-1] + 1e8])
        [ymin, ymax] = ax0.get_ylim()
        ax0.tick_params(direction='in')
        [xmin, xmax] = ax0.get_xlim()

        #use absolute value for the y-axis
        corrected_ylabels = np.array([])
        labels = [np.absolute(int(i)) for i in ax0.get_yticks()]

        ax0.set_yticklabels(labels)

        #change the x-axis labels to chromosome names
        xtick_pos = np.array([])
        for i,e in enumerate(length_list):
            if i == 0:
                continue
            chrom_midpoint = (length_list[i-1] + e) / 2
            xtick_pos = np.append(xtick_pos, chrom_midpoint)

        print_ChrNames = np.array([])
        for i in chr_order:
            print_ChrNames = np.append(print_ChrNames, i[:1].upper() + i[1:])

        ax0.set_xticklabels(print_ChrNames, 
                            rotation='60',
                            va='top',
                            ha='center',
                            style='oblique',
                            family='monospace')

        for i,e in enumerate(length_list):
            if i == 0:
                continue
            if i % 2 == 0:
                ax0.fill_betweenx([ymin, ymax],
                                [length_list[i-1], length_list[i-1]],
                                [e, e],
                                color='#e0e0e0',
                                alpha=0.1)
            if i % 2 == 1:
                ax0.fill_betweenx([ymin, ymax],
                                [length_list[i-1], length_list[i-1]],
                                [e, e],
                                color='#4d4d4d',
                                alpha=0.1)
        #setup the grid
        #[s.set_visible(False) for s in ax0.spines.values()]
        ax0.yaxis.grid(linestyle = '--')
        ax0.set_xticks(xtick_pos)

        #plot a vertical line at the position of enhancer
        ax0.axvline(int(left)+length_list[chr_order.index(enh_chrom)],
                color = '#7A68A6',
                ymin=ymin,
                ymax=ymax,
                linestyle='-.',
                alpha = 0.8)

        #plot the outliers
        outliers = []
        counter = 0
        for i,e in enumerate(plot_y_val):
            if (e < (-1* down_cutoff) or e > up_cutoff):
                outliers.append(counter)
            counter += 1

        if np.any(outliers):
            for j in outliers:
                if (plot_y_val[j] > 1) or (plot_y_val[j] < -1):
                    gene_name = express_subset_df.iloc[j,1]
                    gene_chr = express_subset_df.iloc[j,2]

                    if plot_y_val[j] > 0:
                        ax0.text(plot_x_val[j] + (xmax*0.01), plot_y_val[j] + (ymax*0.01),
                                '%s'%(gene_name),
                                color = '#A60628',
                                fontsize=15)
                    else:
                        ax0.text(plot_x_val[j] + (xmax*0.01), plot_y_val[j] + (ymax*0.01),
                                '%s'%(gene_name),
                                color='#348ABD',
                                fontsize=15)
            outlier_plot(ax0, fc, plot_x_val, plot_y_val, outliers, up_idx, '#A60628')
            outlier_plot(ax0, fc, plot_x_val, plot_y_val, outliers, down_idx, '#348ABD')


        #plot a legend for the circle size
        ax1 = fig.add_subplot(gs[:,10])

        y_len = ymax - ymin
        y_val = []
        for i in range(0,4):
            y_val.append(ymax - y_len * ((i + 1) * 0.05))

        size = [150, 80, 30, 1]
        legend_text = ['>=4-fold', '>=2-fold',
                    '<2-fold', 'not significant']
        for i,size in enumerate(size):
            ax1.scatter(0.5, y_val[i], 
                        color = 'k',
                        s=size,
                        marker='o')
            ax1.text(0.7, y_val[i],
                    '%s'%(legend_text[i]),
                    ha='left',
                    va='center',
                    fontsize = 12)

        ax1.axis('off')
        ax1.set_ylim([ymin, ymax])
        ax1.set_xlim([0.1,1])
        ax1.set_title("Mean CPM FC",
                    ha='left',
                    va='bottom',
                    fontsize = 15)

        fig.savefig(OUTPUT + '/%s-manhattan_plot.pdf'%(region))

if __name__ == '__main__':
    pass
