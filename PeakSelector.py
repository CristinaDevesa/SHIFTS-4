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
import configparser
import logging
import pandas as pd
import numpy as np
import concurrent.futures
from itertools import repeat

def filter_peaks():
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    return

def main(args):
    '''
    Main function
    '''
    
    logging.info("read input DMTable")
    df = pd.read_csv(infile, sep="\t", float_precision='high')
    #filter
    # write apex list in txt

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Modeller',
        epilog='''
        Example:
            python PeakModeller.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakModeller.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='DMTable to be filtered')
    parser.add_argument('-j', '--infile', required=True, help='DMHistogram with PSM and slope information')
    
    parser.add_argument('-s', '--slope', help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', help='Threshold for number of PSMs')
    #parser.add_argument('-m', '--mode', required=True, help='0=filter by slope, 1=filter by frequency, 2=filter by both')
    # ALWAYS FILTER BY BOTH

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.bins is not None:
        config.set('PeakModeller', 'bins', str(args.bins))
        config.set('Logging', 'create_ini', '1')
    if args.points is not None:
        config.set('PeakModeller', 'points', str(args.points))
        config.set('Logging', 'create_ini', '1')
    if args.slope is not None:
        config.set('PeakSelector', 'slope', str(args.slope))
        config.set('Logging', 'create_ini', '1')
    if args.frequency is not None:
        config.set('PeakSelector', 'frequency', str(args.frequency))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/DMcalibrator.ini', 'w') as newconfig:
            config.write(newconfig)
            
    # TODO: add mode to log!

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