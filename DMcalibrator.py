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
pd.options.mode.chained_assignment = None  # default='warn'

#infile = r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\recom.txt"

###################
# Aminoacids #
###################

AAs = {'A': 71.037114,
       'R': 156.101111,
       'N': 114.042927,
       'D': 115.026943,
       'C': 103.009185,
       'E': 129.042593,
       'Q': 128.058578,
       'G': 57.021464,
       'H': 137.058912,
       'I': 113.084064,
       'L': 113.084064,
       'K': 128.094963,
       'M': 131.040485,
       'F': 147.068414,
       'P': 97.052764,
       'S': 87.032028,
       'T': 101.047679,
       'U': 150.953630,
       'W': 186.079313,
       'Y': 163.063320,
       'V': 99.068414,
       'R': 156.101111,
       'Y': 163.063329,
       'W': 186.079313,
       'O': 132.089878}
M_proton = 1.007825
M_oxygen = 15.994915

###################
# Local functions #
###################

# Calibrate mass separately for each raw file.
def readInfile(infile):
    '''    
    Read input file and determine if it is Comet or Recom.
    '''
    df = pd.read_csv(infile, skiprows=1, sep="\t", float_precision='high')
    recom = 0
    if 'Closest_Deltamass' in df.columns:
        recom = 1
    return df, recom

def getTheoMZ(df):
    '''    
    Calculate theoretical MZ using the PSM sequence.
    '''
    df.insert(df.columns.get_loc("exp_mz")+1, 'theo_mz', 0)
    
    def _PSMtoMZ(sequence, charge):
        total_aas = 2*M_proton + M_oxygen
        for aa in sequence:
            if aa in AAs:
                total_aas += AAs[aa]
            #else: # aminoacid not in list (ask for user input?)
                # TODO
        MZ = (total_aas + charge*M_proton) / charge
        return MZ
    
    df['theo_mz'] = df.apply(lambda x: _PSMtoMZ(x['plain_peptide'], x['charge']), axis = 1)
    
    return df


#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    df, recom = readInfile(args.infile)
    
if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Modeller',
        epilog='''
        Example:
            python peak_modeller.py

        ''')
    parser.add_argument('-i', '--infile', required=True, help='Input file from COMET or RECOM')
    
    parser.add_argument('-c', '--cxcorrmin', required=True, help='Minimum cXcorr')
    parser.add_argument('-p', '--ppmmax', required=True, help='Maximum PPM')

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