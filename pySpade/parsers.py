#parsers.py

import sys
import argparse
from pySpade import __version__

def parse_args(args=sys.argv[1:]):

    # create top-level parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=('pySpade \n' + "Version: " +
                     __version__))

    # create sub-level parser
    subparsers = parser.add_subparsers(title='functions',
                                       dest='command',
                                       metavar='')

    add_process_subparser(subparsers)
    add_explevel_subparser(subparsers)
    add_fc_subparser(subparsers)
    add_DEobs_subparser(subparsers)
    add_DErand_subparser(subparsers)
    add_local_subparser(subparsers)
    add_global_subparser(subparsers)
    add_manhattan_subparser(subparsers)

    if len(sys.argv) > 1:
        if (sys.argv[1] == '--version' or sys.argv[1] == '-v'):
            print(" version: %s" % __version__)
            exit()
    else:
        args = parser.parse_args(['-h'])
        exit()

    return parser.parse_args(args)

# process
def add_process_subparser(subparsers):
    parser = subparsers.add_parser(
        'process',
        help='process mapping output and reformat for downstream analysis.',
        description=(
            'Process transcriptome output and sgrna output to remove experimental doublets. '
            'Transcriptome matrix is from Cellranger output (outs folder),'
            'sgrna matrix is from fba output, other accepted format include pkl and csv file.'
            'sgrna matrix column: cell barcodes consistent with transcriptome matrix, rows: sgrna sequence.'
            'The final output format is h5 file.'))

    parser.add_argument('-f', 
                        '--feature_bc', 
                        dest='feature_bc', 
                        required=True,
                        type=str,
                        help='Specify the output folder from cellranger pipeline (outs folder)')

    parser.add_argument('-s', 
                        '--sgrna',
                        dest='input_sgrna',
                        required=True,
                        type=str,
                        help='Specify the sgrna matrix file.'
                        'Please make sure the barcodes are consistent with transcriptome if there are multiple libraries. File format: pkl or csv')

    parser.add_argument('-m', 
                        '--cell_multiplex',
                        dest='cell_multiplex',
                        required=False,
                        default=False,
                        help='specify the used antibody name in txt file separate with new line. Make sure that the antibody is mapped together with transcriptome. Cellranger is recommended.')

    parser.add_argument('-o', 
                        '--output_directory',
                        dest='output_directory', 
                        required=True,
                        help='Specify output directory. The output file format will be h5')

    parser.add_argument('-c', 
                        '--comp', 
                        dest = 'comp', 
                        required=False,
                        type=str,
                        default='False',
                        help='Final h5 file compress or not. "True" or "False", default is False.')

# explevel
def add_explevel_subparser(subparsers):
    parser = subparsers.add_parser(
        'explevel',
        help='check the average expression level of query genes in single cell matrix',
        description=(
            'Check the average expression level of query genes in single cell matrix'
            'Input: processed transcriptome matrix from the process output,'
            'query genes list has to be txt file, genes are seperated with new line.'))

    parser.add_argument('-t', 
                        '--transcriptome_dir', 
                        dest='transcriptome_dir', 
                        required=True,
                        type=str,
                        help='specify the directory from process function.')

    parser.add_argument('-g', 
                        '--gene',
                        dest='gene', 
                        required=True,
                        type=str,
                        help='specify the query genes.')

    parser.add_argument('-o', 
                        '--output_file', 
                        dest='output_file', 
                        required=False,
                        help='specify output file.')

#fc
def add_fc_subparser(subparsers):
    parser = subparsers.add_parser(
        'fc',
        help='check the fold change of genes in given perturbation',
        description=(
            'Check the fold change of perturbed region and individual sgRNA for query region and gene.'
            'Input: processed transcriptome matrix and sgrna matrix from the process output,'
            'sgrna dict file: perturbation region hg38 coordinates and the sgrna sequence targeting that region.'
            'Region and sgrnas separated by tab, and sgrnas separated by comma'))

    parser.add_argument('-t', 
                        '--transcriptome_dir', 
                        dest='transcriptome_dir', 
                        required=True,
                        type=str,
                        help='specify the processed transcriptome matrix file')

    parser.add_argument('-d', 
                        '-dict', 
                        dest='dict', 
                        required=True,
                        type=str, 
                        help='specify the sgRNA annotation file: perturbation coordinates (hg38) and the sgRNA name.')

    parser.add_argument('-r', 
                        '--region', 
                        dest='region', 
                        required=True,
                        type=str,
                        help='specify the query regions and their target genes to calculate repression efficiency.')

    parser.add_argument('-o', 
                        '--output_folder', 
                        dest='output_folder', 
                        required=True,
                        help='specify output folder directory.')

#DEobs
def add_DEobs_subparser(subparsers):
    parser = subparsers.add_parser(
        'DEobs',
        help='perform differential expression analysis of observe cells',
        description=(
            'Perfrom the genome wide differential expression analysis of all the perturbation regions.'
            'Input: processed transcriptome matrix and sgrna matrix from the process output'
            'sgrna dict file: perturbation region hg38 coordinates and the sgrna sequence targeting that region.'
            'region and sgrnas separated by tab, and  sgrnas separated by comma'
            'Output: up regulation p-value, downregulation p-value, fold change(compare with all the otehr cells) and average cpm'))

    parser.add_argument('-t', 
                        '--transcriptome_df', 
                        dest='transcriptome_df', 
                        required=True,
                        type=str,
                        help='specify the processed transcriptome matrix file (.h5)')
    
    parser.add_argument('-s', 
                        '--sgrna', 
                        dest='input_sgrna', 
                        required=True,
                        type=str,
                        help='specify the processed sgrna matrix file. (.h5)')
    
    parser.add_argument('-d', 
                        '-dict', 
                        dest='dict', 
                        required=True,
                        type=str, 
                        help='specify the perturbation coordinates (hg38) and the sgRNA txt file.')

    parser.add_argument('-r', 
                        '--threads', 
                        dest = 'threads', 
                        required=False,
                        type=int,
                        default=1,
                        help='set number of barcode comparison threads. The default is 1')

    parser.add_argument('-n', 
                        '--norm', 
                        dest = 'norm_method', 
                        required=False,
                        type=str,
                        default='cpm',
                        help='choose normalization methods: "cpm" or "metacell".')

    parser.add_argument('-b', 
                        '--bg', 
                        dest = 'bg', 
                        required=False,
                        type=str,
                        default='complement',
                        help='the background cells for comparason. Default is complementary (all the other cells). Specify the key in sgRNA txt file.')
    
    parser.add_argument('-o', 
                        '--output_dir', 
                        dest='output_dir', 
                        required=True,
                        help='specify an output directory.')

#DErand
def add_DErand_subparser(subparsers):
    parser = subparsers.add_parser(
        'DErand',
        help='perform differential expression analysis of random selection background',
        description=(
            'Perfrom the genome wide differential expression analysis of 1000 random selection cells.'
            'There are two options for random selection: all cells with equal probability or probability based on sgrna number in the cells'
            'User should specify the cell number to select randomly.'
            'It is recommended with either exact cell number or bins (large amount of perturbation region).'
            ))

    parser.add_argument('-t', 
                        '--transcriptome_df', 
                        dest='transcriptome_df', 
                        required=True,
                        type=str,
                        help='specify the processed transcriptome matrix file')

    parser.add_argument('-s', 
                        '--sgrna', 
                        dest='input_sgrna', 
                        required=True,
                        type=str,
                        help='specify the processed sgrna matrix file.')
    
    parser.add_argument('-d', 
                        '-dict', 
                        dest='dict', 
                        required=True,
                        type=str, 
                        help='specify the perturbation coordinates (hg38) and the sgRNA txt file.')

    parser.add_argument('-m', 
                        '--num',
                        dest='num',
                        required=True,
                        type=int,
                        help='specify the number of cells to do random iteration.')

    parser.add_argument('-i',
                        '--iteration',
                        dest='iteration',
                        required=False,
                        type=int,
                        default=1000,
                        help='specify the number of iteration to perform.')

    parser.add_argument('-r', 
                        '--threads', 
                        dest = 'threads', 
                        required=False,
                        type=int,
                        default=1,
                        help='set number of barcode comparison threads. The default is 1')
    
    parser.add_argument('-n', 
                        '--norm', 
                        dest = 'norm_method', 
                        required=False,
                        type=str,
                        default='cpm',
                        help='choose normalization methods: "cpm" or "metacell".')

    parser.add_argument('-a', 
                        '--randomization_method', 
                        dest = 'randomization_method', 
                        required=True,
                        type=str,
                        help='choose randomization methods: "equal" or "sgrna".')
    
    parser.add_argument('-b', 
                        '--bg', 
                        dest = 'bg', 
                        required=False,
                        type=str,
                        default='complement',
                        help='the background cells for comparason. Default is complementary (all the other cells). Specify the key in sgRNA txt file.')

    parser.add_argument('-o', 
                        '--output_dir', 
                        dest='output_dir', 
                        required=True,
                        help='specify an output directory.')

#local
def add_local_subparser(subparsers):
    parser = subparsers.add_parser(
        'local',
        help='perform local hit analysis (+-2Mb) with observation data and random background',
        description=('Using the observation p-value and randomization bavckground p-value'
                     'to calculate the adjusted p-value based on gamma distribution approximation'
                     'Local hits calculation includes the genes within plus and minus 2 Mb of the perturbation region.'
                     'The output is a csv file with all hits information.'))
   
    parser.add_argument('-f', 
                        '--file_dir', 
                        dest='file_dir', 
                        required=True,
                        type=str,
                        help='specify the file directory of "process" function output, the Trans_genome_seq.npy file is required at this step.')
   
    parser.add_argument('-d', 
                        '--data_dir', 
                        dest='data_dir', 
                        required=True,
                        type=str,
                        help='specify the p-value matrix directory of observation test.')

    parser.add_argument('-t', 
                        '--distr', 
                        dest='distr', 
                        required=True,
                        type=str,
                        help='specify the random cell mean/std/10_perc file directory.')
    
    parser.add_argument('-s', 
                        '-sgrna_dict', 
                        dest='sgrna_dict', 
                        required=True,
                        type=str, 
                        help='specify the perturbation coordinates (hg38) and the sgRNA txt file.')

    parser.add_argument('-o', 
                        '--output_file', 
                        dest='output_file', 
                        required=True,
                        type=str,
                        help='specify an output file name incluseing the directory, it has to be in csv format.')

#global
def add_global_subparser(subparsers):
    parser = subparsers.add_parser(
        'global',
        help='perform global hit analysis with observation data and random background',
        description=('Using the observation p-value and randomization bavckground p-value'
                     'to calculate the adjusted p-value based on gamma distribution approximation'
                     'The output is a csv file with all hits information.'))
    
    parser.add_argument('-f', 
                        '--file_dir', 
                        dest='file_dir', 
                        required=True,
                        type=str,
                        help='specify the file directory of "process" function output, the Trans_genome_seq.npy file is required at this step.')
    
    parser.add_argument('-d', 
                        '--data_dir', 
                        dest='data_dir', 
                        required=True,
                        type=str,
                        help='specify the p-value matrix directory of observation test. (DEobs output folder)')

    parser.add_argument('-s', 
                        '-sgrna_dict', 
                        dest='sgrna_dict', 
                        required=True,
                        type=str, 
                        help='specify the perturbation coordinates (hg38) and the sgRNA txt file.')
        
    parser.add_argument('-t', 
                        '--distr', 
                        dest='distr', 
                        required=True,
                        type=str,
                        help='specify the random cell file directory. (DErand output folder)')

    parser.add_argument('-o', 
                        '--output_file', 
                        dest='output_file', 
                        required=True,
                        type=str,
                        help='specify an output file name including the directory, it has to be in csv format.')
    
#manhattan
def add_manhattan_subparser(subparsers):
    parser = subparsers.add_parser(
        'manhattan',
        help='generate Manhattan plots for each perturbation region',
        description=('Use the output csv file from global function to generate Manhattan plots for each perturbation region'))
    
    parser.add_argument('-f', 
                        '--file_dir', 
                        dest='file_dir', 
                        required=True,
                        type=str,
                        help='specify the file directory of "process" function output, the Trans_genome_seq.npy file is required at this step.')
    
    parser.add_argument('-g', 
                        '--global_csv', 
                        dest='global_csv', 
                        required=True,
                        type=str,
                        help='specify the csv file directory from the output of global function.')
              
    parser.add_argument('-cx', 
                        '--cutoff_expression', 
                        dest='cutoff_expression', 
                        required=False,
                        type=float,
                        default=0.05,
                        help='specify the cutoff of expressed genes. Default is 0.05 (genes expressed in more than 5 percent of cells)')

    parser.add_argument('-cf', 
                        '--cutoff_fc', 
                        dest='cutoff_fc', 
                        required=False,
                        type=float,
                        default=0.2,
                        help='specify the cutoff of fold change. Default is 0.2 (fold change is more than 20 percent)')
    
    parser.add_argument('-cs', 
                        '--cutoff_significance', 
                        dest='cutoff_significance', 
                        required=False,
                        type=float,
                        default=-5,
                        help='specify the cutoff of Significance_score. Default is -5 (Significance score is smaller than -5)')
    
    parser.add_argument('-o', 
                        '--output_folder', 
                        dest='output_folder', 
                        required=True,
                        type=str,
                        help='specify an output folder directory.')
