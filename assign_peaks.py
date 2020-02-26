#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Jose Rodriguez"
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import sys
import argparse
import logging
import pandas as pd
import concurrent.futures

###################
# Local functions #
###################
def concatInfiles(file):
    '''
    Pre-processing the data: assign target-decoy, correct monoisotopic mass, calculate cXCorr
    '''    
    # read input file
    df = pd.read_csv(file, sep="\t")
    return df

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''    
    logging.info("get the list of files with the inputs")
    with open(args.infile) as f:
        infiles = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    infiles = [x.strip() for x in infiles] 
    logging.debug(infiles)


    logging.info("concat input files")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        ddf = executor.map(concatInfiles, infiles)
    logging.info("concat")
    ddf = pd.concat(ddf)
    
    logging.info("sort by DeltaMax cal")
    ddf.sort_values(by=['expMH'], inplace=True)
    ddf.reset_index(drop=True, inplace=True)
    
    ddf.to_csv("test.tsv", sep="\t")
    

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Assign peaks',
        epilog='''
        Example:
            python assign_peaks.py

        ''')
    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')
    parser.add_argument('-i',  '--infile', required=True, help='Input file with the list of files that contains the peak picking')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # logging debug level. By default, info level
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