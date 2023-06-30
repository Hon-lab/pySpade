#__main__.py

import sys
import importlib
from pySpade import __version__
from pySpade.parsers import parse_args
from pySpade.utils import get_logger

def main():
    args = parse_args()

    logger = get_logger(logger_name=__name__)
    logger.info(f'pySpade version: {__version__}')
    logger.info(f'Python version: {sys.version_info.major}.{sys.version_info.minor}')

    if not sys.version_info.major == 3 and sys.version_info.minor >= 7:
        logger.critical('Please use Python >= 3.7')
        sys.exit(1)


    if (args.command == 'process'):
        logger.info('Running process command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.process_singlet(
            WORK_DIR=args.feature_bc,
            SGRNA=args.input_sgrna,
            CELL_MULTIPLEX=args.cell_multiplex,
            OUTPUT_DIR=args.output_directory,
            COMP=args.comp
        )
        logger.info('Done.')


    elif (args.command == 'explevel'):
        logger.info('Running explevel command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.Calculate_gene_expression(
            TRANSCRIPTOME_DIR=args.transcriptome_dir,
            GENE_FILE=args.gene,
            OUTPUT_FILE=args.output_file
        )
        logger.info('Done.')


    elif (args.command == 'fc'):
        logger.info('Running fc command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.Calculate_fc(
            TRANSCRIPTOME_DIR=args.transcriptome_dir,
            SGRNA_DICT = args.dict,
            TARGET_FILE = args.region,
            OUTPUT_FOLDER = args.output_folder
        )        
        logger.info('Done.')


    elif (args.command == 'DEobs'):
        logger.info('Running DEobs command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.DE_observe_cells(
            sub_df_file = args.transcriptome_df,
            sgrna_df = args.input_sgrna, 
            sgrnas_file = args.dict,
            num_processing = args.threads,
            norm = args.norm_method,
            output_dir = args.output_dir
            )        
        logger.info('Done.')

    elif (args.command == 'DErand'):
        logger.info('Running DErand command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.DE_random_cells(
            sub_df_file = args.transcriptome_df,
            sgrna_df = args.input_sgrna,
            num_cell = args.num,
            iteration = args.iteration,
            num_processing = args.threads,
            norm = args.norm_method,
            RAND_M = args.randomization_method,
            output_dir = args.output_dir
            )        
        logger.info('Done.')

    elif (args.command == 'local'):
        logger.info('Running local command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.local_analysis(
            FILE_DIR = args.file_dir,
            OBS_DIR = args.data_dir,
            DISTRI_DIR = args.distr,
            SGRNA_DICT = args.sgrna_dict,
            OUTPUT_DF = args.output_file
            )        
        logger.info('Done.')

    elif (args.command == 'global'):
        logger.info('Running global command ...')
        m = importlib.import_module(name=f'pySpade.{args.command}')

        _ = m.global_analysis(
            FILE_DIR = args.file_dir,
            OBS_DIR = args.data_dir,
            SGRNA_DICT = args.sgrna_dict,
            DISTRI_DIR = args.distr,
            OUTPUT_DF = args.output_file
            )        
        logger.info('Done.')

if __name__ == "__main__":
    main()