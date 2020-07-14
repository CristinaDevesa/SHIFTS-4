#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Andrea Laguillo Gómez"
__email__ = "jmrodriguezc@cnic.es;andrea.laguillo@cnic.es"
__status__ = "Development"

# import modules
import argparse
import configparser
import concurrent.futures
from itertools import repeat
import logging
import numpy as np
import os
import pandas as pd
import sys

# TODO recom parameter: if no recom column specified, proceed without spire FDR

###################
# Local functions #
###################

def get_spire_FDR(df, xcorr_type): #TODO: we don't have xcorr_type, we have recom_data, take out column names
    #This will be for the group of scans in a peak that are contained within 
    #one recom-assignated theoretical deltamass. Then, when we do peak_FDR, we
    #include these as well as the rest of the values in the peak.
    #Cuando el peak y el spire se solapen, esos escanes lo sometería a ambas FDR
    # How we will handle filtering is still to be determined.
    '''
    Calculate spire FDR for each spire in one bin (1 Da)
    '''
    df['SpireFDR'] = -1
    df['Rank'] = -1
    df['Spire_Rank_T'] = -1
    df['Spire_Rank_D'] = -1
    # TODO: Operations
    # identify spires (filter by recom HERE OR IN RECOM?)
    peaks = df[df['Peak'] == 'PEAK'] # filter by Peak
    recom_peaks = peaks[peaks['XcorType'] == 'RECOM'] # TODO: XcorType needs to be created
    #recom-identified scans are not necessarily peaks, what to do?
    grouped_recom_peaks = recom_peaks.groupby(['ClosestPeak']) # group by ClosestPeak
    for group in grouped_recom_peaks:
        group_index = group[1].index.values
        df.loc[group_index] # repeat steps of local_FDR
        # sort bin
        if xcorr_type == 0: # by Comet Xcorr
            df.loc[group_index].sort_values(by=['Xcor', 'Label'], inplace=True)
        else: # by Comet cXcorr
            df.loc[group_index].sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
        # count targets and decoys
        df.loc[group_index]['Rank'] = df.loc[group_index].groupby('Label').cumcount()+1 # This column can be deleted later
        df.loc[group_index]['Spire_Rank_T'] = np.where(df.loc[group_index]['Label']=='Target', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Spire_Rank_T'] = df.loc[group_index]['Spire_Rank_T'].replace(to_replace=0, method='ffill')
        df.loc[group_index]['Spire_Rank_D'] = np.where(df.loc[group_index]['Label'] == 'Decoy', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Spire_Rank_D'] =  df.loc[group_index]['Spire_Rank_D'].replace(to_replace=0, method='ffill')
        # calculate local FDR
        df.loc[group_index]['SpireFDR'] = df.loc[group_index]['Spire_Rank_D']/df.loc[group_index]['Spire_Rank_T']
    # TODO: End Operations
    df.drop(['Rank'], axis = 1, inplace = True)
    return df

def get_peak_FDR(df, recom_data):
    '''
    Calculate peak FDR for each peak in one bin (1 Da)
    '''
    df['PeakFDR'] = -1
    df['Rank'] = -1
    df['Peak_Rank_T'] = -1
    df['Peak_Rank_D'] = -1
    # identify peaks
    peaks = df[df['Peak'] == 'PEAK'] # filter by Peak
    grouped_peaks = peaks.groupby(['ClosestPeak']) # group by ClosestPeak
    # df.get_group("group")
    #grouped_peaks.groups # group info
    for group in grouped_peaks:
        group_index = group[1].index.values
        df.loc[group_index] # repeat steps of local_FDR
        # sort bin
        if recom_data == 0: # by Comet Xcorr
            df.loc[group_index].sort_values(by=['Xcor', 'Label'], inplace=True)
        else: # by Comet cXcorr
            df.loc[group_index].sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
        # count targets and decoys
        df.loc[group_index]['Rank'] = df.loc[group_index].groupby('Label').cumcount()+1 # This column can be deleted later
        df.loc[group_index]['Peak_Rank_T'] = np.where(df.loc[group_index]['Label']=='Target', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Peak_Rank_T'] = df.loc[group_index]['Peak_Rank_T'].replace(to_replace=0, method='ffill')
        df.loc[group_index]['Peak_Rank_D'] = np.where(df.loc[group_index]['Label'] == 'Decoy', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Peak_Rank_D'] =  df.loc[group_index]['Peak_Rank_D'].replace(to_replace=0, method='ffill')
        # calculate local FDR
        df.loc[group_index]['PeakFDR'] = df.loc[group_index]['Peak_Rank_D']/df.loc[group_index]['Peak_Rank_T']
    df.drop(['Rank'], axis = 1, inplace = True)
    return df

def get_local_FDR(df, recom_data):
    '''
    Calculate local FDR for one bin (1 Da)
    '''
    # sort bin
    if recom_data == 0: # by Comet Xcorr
        df.sort_values(by=['Xcor', 'Label'], inplace=True)
    else: # by Comet cXcorr
        df.sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
        
    # count targets and decoys
    df['Rank'] = df.groupby('Label').cumcount()+1 # This column can be deleted later
    df['Local_Rank_T'] = np.where(df['Label']=='Target', df['Rank'], 0)
    df['Local_Rank_T'] = df['Local_Rank_T'].replace(to_replace=0, method='ffill')
    df['Local_Rank_D'] = np.where(df['Label'] == 'Decoy', df['Rank'], 0)
    df['Local_Rank_D'] =  df['Local_Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # calculate local FDR
    df['LocalFDR'] = df['Local_Rank_D']/df['Local_Rank_T']
    return df

def get_global_FDR(df, recom_data):
    '''
    Calculate global FDR
    '''
    # sort by score
    if recom_data == 0: # by Comet Xcorr
        df.sort_values(by=['Xcor', 'Label'], inplace=True)
    else: # by Comet cXcorr
        df.sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
        
    # count targets and decoys
    df['Rank'] = df.groupby('Label').cumcount()+1 # This column can be deleted later
    df['Global_Rank_T'] = np.where(df['Label']=='Target', df['Rank'], 0)
    df['Global_Rank_T'] = df['Global_Rank_T'].replace(to_replace=0, method='ffill')
    df['Global_Rank_D'] = np.where(df['Label'] == 'Decoy', df['Rank'], 0)
    df['Global_Rank_D'] =  df['Global_Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # calculate local FDR
    df['LocalFDR'] = df['Local_Rank_D']/df['Local_Rank_T']
    return df

def filtering(df, fdr_filter, target_filter):
    return df

def bin_operations(df, recom_data, peak_label):
    '''
    Main function that handles the operations by BIN
    '''
    # calculate local FDR
    df = get_local_FDR(df, recom_data)
    
    # calculate peak FDR
    df = get_peak_FDR(df, recom_data)
    
    # calculate spire FDR
    df = get_spire_FDR(df, recom_data)
    
    return df

#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    # Main variables
    score_column = config._sections['PeakFDRer']['score_column']
    fdr_filter = config._sections['PeakFDRer']['fdr_filter']
    target_filter = config._sections['PeakFDRer']['target_filter']
    recom_data = config._sections['PeakFDRer']['recom_data']
    peak_label = config._sections['PeakAssignator']['peak_label']
    col_CalDeltaMH = config._sections['PeakAssignator']['caldeltamh_column']
    
    #Read input file
    df = pd.read_feather(args.infile)
    
    logging.info("parallel the operations by BIN")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:        
        df = executor.map(bin_operations, list(df.groupby("bin")), repeat(recom_data), 
                                                                   repeat(peak_label) # TODO: missing args
                                                                   ) 
    df = pd.concat(df)
    logging.info("Calculate gobal FDR")
    df = get_global_FDR(df, args.xcorr)
    logging.info("Sort by DeltaMax cal")
    df.sort_values(by=[col_CalDeltaMH], inplace=True)
    
    # Filtering
    df = filtering(df, fdr_filter, target_filter)
    df.reset_index(drop=True, inplace=True)

    # d_h = df.head()
    # d_t = df.tail()
    # d_h.to_csv("kk_head.tsv", sep="\t")
    # d_t.to_csv("kk_tail.tsv", sep="\t")
    

    logging.info("print output")
    # https://towardsdatascience.com/the-best-format-to-save-pandas-data-414dca023e0d
    # begin:printHDF5
    # Note: Explote the Memory!!!
    # assign NumExpr for the tables module
    # tables.parameters.MAX_NUMEXPR_THREADS = args.n_workers
    # df.to_hdf('data.h5', key='df', mode='w')
    # end:printHDF5
    # df.to_csv('data.tsv', sep="\t", index=False)
    df.to_feather('data.ftr')
    

    

if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak FDRer',
        epilog='''
        Example:
            python PeakFDRer.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakModeller.ini")
    
    parser.add_argument('-i',  '--infile', required=True, help='Input feather file with the peak assignation')
    
    parser.add_argument('-s',  '--score_column', help='Name of column with score for FDR calculation')
    parser.add_argument('-f',  '--fdr_filter', help='FDR value to filter by')
    parser.add_argument('-t',  '--target_filter', help='Filter targets, 0=no 1=yes')
    parser.add_argument('-r',  '--recom_data', help='Score for FDR calculation: 0=Xcorr, 1=cXcorr (default: %(default)s)')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.score_column is not None:
        config.set('PeakFDRer', 'score_column', str(args.score_column))
        config.set('Logging', 'create_ini', '1')
    if args.fdr_filter is not None:
        config.set('PeakFDRer', 'fdr_filter', str(args.fdr_filter))
        config.set('Logging', 'create_ini', '1')
    if args.target_filter is not None:
        config.set('PeakFDRer', 'target_filter', str(args.target_filter))
        config.set('Logging', 'create_ini', '1')
    if args.recom_data is not None:
        config.set('PeakFDRer', 'recom_data', str(args.recom_data))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/PeakModeller.ini', 'w') as newconfig:
            config.write(newconfig)

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