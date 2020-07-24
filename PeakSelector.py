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
import re
import pandas as pd
import numpy as np

def readHistogram(infile):
    df_hist = pd.read_csv(args.infile, sep="\t", float_precision='high')
    df_hist = df_hist.dropna() # Remove rows with missing values (will always have some in beginning and end)
    return df_hist

def filterPeaks(df_hist, slope, frequency):
    '''
    Find peaks that are above the thresholds for slope and PSMs.
    '''
    # TODO: allow specify slope and count columns in INI?
    df_hist = df_hist[abs(df_hist['slope1']) >= slope]
    df_hist = df_hist[df_hist['count'] >= frequency]
    return df_hist

def parseInterval(bins_df):
    '''
    Read 'bin' column as an interval.
    '''
    for i in range(0, len(bins_df)):
        to_interval = bins_df.loc[i, 'bin']
        left = float(re.findall(r'-?\d+\.\d+', to_interval)[0]) 
        right = float(re.findall(r'-?\d+\.\d+', to_interval)[1])
        bins_df.loc[i, 'bin'] = pd.Interval(left, right, closed='right')
    return bins_df

def areValid(intervals):
    '''
    Check whether intervals in a list are contiguous, have a change in
    sign of the slope, and the central point is the closest to 0.
    '''
    cont = 0
    zero_bin = min(intervals['slope2'].tolist(), key=abs)
    zero_index = intervals['slope2'].tolist().index(zero_bin)
    if zero_index == len(intervals)//2: # Central point is closest to 0
       if (intervals.loc[zero_index, 'slope2'] >= 0 and intervals.loc[zero_index+1, 'slope2'] < 0) or (intervals.loc[zero_index, 'slope2'] <= 0 and intervals.loc[zero_index-1, 'slope2'] > 0): # Change in sign  
           bin_list = intervals['bin'].tolist()
           cont = 1
           for i in range(1, len(bin_list)):
               if bin_list[i-1].right != bin_list[i].left: # Not contiguous
                   cont = 0
    if cont == 0:
        return False
    else:
        return True
    
def interpolateApex(bin_subset):
    x_list = bin_subset['midpoint'].tolist()
    y_list = bin_subset['count'].tolist()
    sum1, sum2 = 0, 0
    for i in range(len(x_list)):
        sum1 += (x_list[i] - np.mean(x_list)) * (y_list[i] - np.mean(y_list))
        sum2 += (x_list[i] - np.mean(x_list)) ** 2
    working_slope = sum1 / sum2
    intercept = np.mean(y_list) - working_slope*np.mean(x_list)
    apex = -intercept / working_slope # x where y=0
    return apex

def peakApex(bins_df, apex_points):
    '''
    Calculate apex for each peak.
    '''
    apex_list = []
    # for i in range(1, len(bins_df)):
    #     if bins_df.loc[i, 'slope2'] is not None:
    #         i1 = bins_df.loc[i-1, 'bin'] # TODO parse interval
    #         i2 = bins_df.loc[i, 'bin']
    #         # Check intervals are consecutive, and there is a change in sign of slope2
    #         if i1.right == i2.left and bins_df.loc[i, 'slope2'] < 0 and bins_df.loc[i-1, 'slope2'] >= 0:
    #             peak = pd.Series([bins_df.loc[i-1, 'midpoint'], np.nan, bins_df.loc[i, 'midpoint']],
    #                               index=[bins_df.loc[i-1, 'slope2'], 0, bins_df.loc[i, 'slope2']])
    #             peak = peak.interpolate(method='index')
    #             apex_list.append(peak[0])
    for i in range(apex_points//2, len(bins_df)-apex_points//2):
        # Check there is a change of sign
        intervals = []
        for j in range(i-apex_points//2, i+apex_points//2+1):
            intervals.append(bins_df.loc[j])
        intervals = pd.DataFrame(intervals)
        intervals.reset_index(drop=True, inplace=True)
        if areValid(intervals):
            peak = interpolateApex(intervals)
            apex_list.append(peak)
    return apex_list

def main(args):
    '''
    Main function
    '''
    
    # Main variables
    slope = float(config._sections['PeakSelector']['slope'])
    frequency = int(config._sections['PeakSelector']['frequency'])
    apex_points = int(config._sections['PeakSelector']['apex_points'])
    
    # Read DM Histogram
    logging.info("Reading input file...")
    df_hist = readHistogram(args.infile)
    df_hist.reset_index(drop=True, inplace=True)
    
    # Filter by slope and frequency, calculate apexes
    logging.info("Filtering...")
    df_hist = filterPeaks(df_hist, slope, frequency)
    df_hist.reset_index(drop=True, inplace=True)
    df_hist = parseInterval(df_hist)
    apex_list = peakApex(df_hist, apex_points)
    apex_info = str(len(apex_list)) + " peaks"
    logging.info(apex_info)
    
    # Write apex list
    logging.info("Writing apex list...")
    outfile = args.infile[:-15] + 'ApexList.txt'
    with open(outfile, 'w') as f:
        for apex in apex_list:
            f.write("%s\n" % apex)
    logging.info("Peak Selection finished")

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Selector',
        epilog='''
        Example:
            python PeakSelector.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakModeller.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='DMHistogram to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    parser.add_argument('-s', '--slope', help='Threshold for slope of DM peak')
    parser.add_argument('-f', '--frequency', help='Threshold for number of PSMs')
    parser.add_argument('-p', '--apex_points', help='Number of points (bins) to use for apex calculation')

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
    if args.apex_points is not None:
        config.set('PeakSelector', 'apex_points', str(args.apex_points))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/PeakModeller.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-15] + 'ApexList_log.txt'
    log_file_debug = outfile = args.infile[:-15] + 'ApexList_log_debug.txt'
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