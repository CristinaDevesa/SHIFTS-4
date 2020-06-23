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

def filter_peaks(df_hist, slope, frequency):
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    # TODO: allow specify slope and count columns in INI?
    df_hist = df_hist[df_hist['slope1'] >= slope]
    df_hist = df_hist[df_hist['count'] >= frequency]
    return df_hist

def peak_apex(bins_df):
    '''
    Calculate apex for each peak.
    '''
    apex_list = []
    for i in range(1, len(bins_df)):
        if bins_df.loc[i, 'slope2'] is not None:
            i1 = bins_df.loc[i-1, 'bin']
            i2 = bins_df.loc[i, 'bin']
            # Check intervals are consecutive, and there is a change in sign of slope2
            if i1.right == i2.left and bins_df.loc[i, 'slope2'] < 0 and bins_df.loc[i-1, 'slope2'] >= 0:
                peak = pd.Series([bins_df.loc[i-1, 'midpoint'], np.nan, bins_df.loc[i, 'midpoint']],
                                  index=[bins_df.loc[i-1, 'slope2'], 0, bins_df.loc[i, 'slope2']])
                peak = peak.interpolate(method='index')
                apex_list.append(peak[0])
    return apex_list

def main(args):
    '''
    Main function
    '''
    
    logging.info("Reading input file...")
    df_hist = pd.read_csv(args.histogram, sep="\t", float_precision='high')
    logging.info("Filtering...")
    df_hist = filter_peaks(df_hist, #slope, #frequency)
    apex_list = peak_apex(df_hist) #TODO: could there be problems bc of filtering?
    # write apex list in txt
    logging.info("Writing apex list...")
    outfile = args.infile[:-8] + '_ApexList.txt'
    with open(outfile, 'w') as f:
        for apex in apex_list:
            f.write("%s\n" % apex)
    logging.info("Peak Selection finished")

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
    
    parser.add_argument('-i', '--infile', required=True, help='DMHistogram to be filtered')
    
    parser.add_argument('-s', '--slope', help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', help='Threshold for number of PSMs')
    #parser.add_argument('-m', '--mode', required=True, help='0=filter by slope, 1=filter by frequency, 2=filter by both')
    # ALWAYS FILTER BY BOTH

    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
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

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')