#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "andrea.laguillo@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import concurrent.futures
from itertools import repeat

#infile = r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\cXcorr_Len_Rank_Results_TargetData_Calibration.txt"

###################
# Local functions #
###################
def concatInfiles(infile, fwhm_fname):
    '''    
    Concat input files...
    '''
  
    # read input file
    # use high precision with the floats
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    # assign type to categorical columns
    df['filename'] = df['filename'].astype('category')
    df['Label'] = df['Label'].astype('category')
    df['IsotpicJump'] = df['IsotpicJump'].astype('category')
    return df

def smoothing():
    '''
    1
    '''
    return

def generate_histogram(df, bin_width):
    '''
    Group by DeltaMass into bins of the size specified.
    '''
    
    def _decimal_places(x):
        s = str(x)
        if not '.' in s:
            return 0
        return len(s) - s.index('.') - 1
    
    # sort by deltamass
    df.sort_values(by=['Cal_Delta_MH'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # make bins
    bins = list(np.arange(int(round(df['Cal_Delta_MH'][0])),
                          int(round(df['Cal_Delta_MH'].iloc[-1]))+bin_width,
                          bin_width))
    bins = [round(x, _decimal_places(bin_width)) for x in bins]
    df['bin'] = pd.cut(df['Cal_Delta_MH'], bins=bins)
    
    return df

def first_derivative():
    '''
    1
    '''
    return

def second_derivative():
    '''
    1
    '''
    return

def filter_peaks():
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    return

def main(args):
    '''
    Main function
    '''
    
    logging.info("read input file list")
    with open(args.infile) as f:
        infiles = f.readlines()
    infiles = [x.strip() for x in infiles] # remove whitespace
    
    logging.info("concat input files") # TODO: do we do this for each file or all together
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        df = executor.map(concatInfiles, infiles, repeat(args.fwhm_filename))
    df = pd.concat(df)
    df.reset_index(drop=True, inplace=True)
    
    # first pass for smoothing
    # second pass for filtering:
        # make bins
    df = generate_histogram(df, args.bins)
        # calculate derivatives
        # check which bins pass
    # write apex list in txt

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Modeller',
        epilog='''
        Example:
            python peak_modeller.py

        ''')
    parser.add_argument('-i', '--infile', required=True, help='Input file with the peak file(s) to be filtered')
    
    parser.add_argument('-b', '--bins', required=True, help='Width of the bins')
    parser.add_argument('-p', '--points', required=True, help='Number of points (bins) to use for slope calculation')
    parser.add_argument('-s', '--slope', required=True, help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', required=True, help='Threshold for number of PSMs')
    #parser.add_argument('-m', '--mode', required=True, help='0=filter by slope, 1=filter by frequency, 2=filter by both')
    # ALWAYS FILTER BY BOTH

    parser.add_argument('-w', '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
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