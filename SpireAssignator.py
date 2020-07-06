#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Jose Rodriguez"
__credits__ = ["Jose Rodriguez", "Andrea Laguillo Gómez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es;andrea.laguillo@cnic.es"
__status__ = "Development"

# import modules
import numpy as np
import pandas as pd

def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    #df = pd.read_csv(infile, skiprows=1, sep="\t", float_precision='high')
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    return df

def pickSpires(df, percentage, cometcol, recomcol, outcol):
    '''    
    label rows where the Recom score improves over the Comet score more
    than the threshold percentage.
    '''
    df[outcol] = np.nan
    
    return df

def main(args):
    '''
    Main function
    '''
    # Variables
    percentage = float(config._sections['SpireAssignator']['percentage'])
    cometcol = config._sections['SpireAssignator']['cometcolumn']
    recomcol = config._sections['SpireAssignator']['recomcolumn']
    outcol = config._sections['SpireAssignator']['spirecolumn']
    # Read infile
    logging.info("Reading input file...")
    df = readInfile(Path(args.infile))
    # Create Spire column
    df = pickSpires(df,
                    percentage,
                    cometcol,
                    recomcol,
                    outcol)

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Spire Assignator',
        epilog='''
        Example:
            python SpireAssignator.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakModeller.ini")
    
    parser.add_argument('-i',  '--infile', required=True, help='Input file that contains the peak picking')  
    
    parser.add_argument('-p', '--percentage', help='Threshold for % of improvement in Recom score')
    parser.add_argument('-c', '--cometcolumn', help='Name of column containing Comet score')
    parser.add_argument('-r', '--recomcolumn', help='Name of column containing Recom score')

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.percentage is not None:
        config.set('SpireAssignator', 'percentage', str(args.percentage))
        config.set('Logging', 'create_ini', '1')
    if args.cometcolumn is not None:
        config.set('SpireAssignator', 'cometcolumn', str(args.cometcolumn))
        config.set('Logging', 'create_ini', '1')
    if args.recomcolumn is not None:
        config.set('SpireAssignator', 'recomcolumn', str(args.recomcolumn))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1: #TODO: check that other modules will not overwrite
        with open(os.path.dirname(args.infile) + '/PeakModeller.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_spire_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_spire_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')