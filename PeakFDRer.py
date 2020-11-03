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

def read_experiments(experiments_table):
    '''
    Read input file containing groups and filenames in tab-separated format.
    '''
    df = pd.read_csv(experiments_table, sep="\t", names=['Experiment', 'Filename'])
    df['Experiment'] = df['Experiment'].astype('string')
    df['Filename'] = df['Filename'].astype('string')
    if df['Filename'].duplicated().any(): # Check no repeats
        sys.exit("Experiments table contains repeat values in the filename column")
    #exp_groups = exp_df.groupby(by = exp_df.columns[0], axis = 0)
    #for position, exp in exp_groups:
        #TODO: read filepath or everything in folder
    return df

def make_groups(df, groups):
    '''
    Add group column to input file with the peak assignation.
    '''
    def _match_file(groups, filename):
        # if filename in groups['Filename'].unique():
            # group = df.loc[df['Filename'] == filename]['Experiment']
            # group.reset_index(drop=True, inplace=True)
            # group = group[0]
        # if filename in [x for v in group_dict.values() for x in v]:
        if filename in group_dict:
            group = group_dict.get(filename)[0]
        else:
            group = 'N/A'
        return group
    df['Experiment'] = 'N/A'
    #df['Experiment'] = df.apply(lambda x: _match_file(groups, x['Filename']), axis = 1)
    ###
    group_dict = {}
    for x in range(len(groups)):
        currentid = groups.iloc[x,1]
        currentvalue = groups.iloc[x,0]
        group_dict.setdefault(currentid, [])
        group_dict[currentid].append(currentvalue)
    df['Experiment'] = np.vectorize(_match_file)(group_dict, df['Filename'])
    ###
    #df['Experiment'] = _match_file(groups, df['Filename'])
    if 'N/A' in df['Experiment'].unique():
        logging.info('Warning: ' + str(df['Experiment'].value_counts()['N/A']) + ' rows could not be assigned to an experiment!') # They will all be grouped together for FDR calculations
    return df

def get_spire_FDR(df, score_column, col_Peak, xcorr_type): #TODO: we don't have xcorr_type, we have recom_data, take out column names
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

def get_peak_FDR(df, score_column, col_Peak, closestpeak_column, recom_data):
    '''
    Calculate peak FDR for each peak in one bin (1 Da)
    '''
    df['PeakFDR'] = -1
    df['Rank'] = -1
    df['Peak_Rank_T'] = -1
    df['Peak_Rank_D'] = -1
    # identify peaks
    peaks = df[df[col_Peak] == 'PEAK'] # filter by Peak
    grouped_peaks = peaks.groupby([closestpeak_column]) # group by ClosestPeak
    # df.get_group("group")
    #grouped_peaks.groups # group info
    for group in grouped_peaks:
        group_index = group[1].index.values
        df.loc[group_index] # repeat steps of local_FDR
        # sort bin
        # if recom_data == 0: # by Comet Xcorr
        #     df.loc[group_index].sort_values(by=['Xcor', 'Label'], inplace=True)
        # else: # by Comet cXcorr
        #     df.loc[group_index].sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
        df.loc[group_index].sort_values(by=[score_column, 'Label'], inplace=True)
        # count targets and decoys
        df.loc[group_index]['Rank'] = df.loc[group_index].groupby('Label').cumcount()+1 # This column can be deleted later
        df.loc[group_index]['Peak_Rank_T'] = np.where(df.loc[group_index]['Label']=='Target', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Peak_Rank_T'] = df.loc[group_index]['Peak_Rank_T'].replace(to_replace=0, method='ffill')
        df.loc[group_index]['Peak_Rank_D'] = np.where(df.loc[group_index]['Label'] == 'Decoy', df.loc[group_index]['Rank'], 0)
        df.loc[group_index]['Peak_Rank_D'] =  df.loc[group_index]['Peak_Rank_D'].replace(to_replace=0, method='ffill')
        # calculate peak FDR
        df.loc[group_index]['PeakFDR'] = df.loc[group_index]['Peak_Rank_D']/df.loc[group_index]['Peak_Rank_T']
    df.drop(['Rank'], axis = 1, inplace = True)
    return df

def get_local_FDR(df, score_column, recom_data):
    '''
    Calculate local FDR for one bin (1 Da)
    '''
    # sort bin
    #if recom_data == 0: # by Comet Xcorr
        #df.sort_values(by=['Xcor', 'Label'], inplace=True)
    #else: # by Comet cXcorr
        #df.sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
    df.sort_values(by=[score_column, 'Label'], inplace=True)
        
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

def get_global_FDR(df, score_column, recom_data):
    '''
    Calculate global FDR
    '''
    # get the EXPERIMENT value from the input tuple df=(experiment,df)
    (experiment_value, df) = df[0], df[1]
    # sort by score
    # if recom_data == 0: # by Comet Xcorr
    #     df.sort_values(by=['Xcor', 'Label'], inplace=True)
    # else: # by Comet cXcorr
    #     df.sort_values(by=['CorXcor', 'Label'], inplace=True) # TODO: Fix SHIFTS cXcorr
    df.sort_values(by=[score_column, 'Label'], inplace=True)
        
    # count targets and decoys
    df['Rank'] = df.groupby('Label').cumcount()+1 # This column can be deleted later
    df['Global_Rank_T'] = np.where(df['Label']=='Target', df['Rank'], 0)
    df['Global_Rank_T'] = df['Global_Rank_T'].replace(to_replace=0, method='ffill')
    df['Global_Rank_D'] = np.where(df['Label'] == 'Decoy', df['Rank'], 0)
    df['Global_Rank_D'] =  df['Global_Rank_D'].replace(to_replace=0, method='ffill')
    df.drop(['Rank'], axis = 1, inplace = True)
    
    # calculate global FDR
    df['GlobalFDR'] = df['Global_Rank_D']/df['Global_Rank_T']
    return df

def filtering(df, fdr_filter, target_filter): # This goes on a separate module now
    if target_filter: # =! 0
        df[df['Label'] == 'Target']
    if fdr_filter: # =! 0
        df[df['GlobalFDR'] >= fdr_filter]
    return df

def bin_operations(df, score_column, recom_data, peak_label, col_Peak, closestpeak_column):
    '''
    Main function that handles the operations by BIN
    '''
    
    # get the BIN value from the input tuple df=(bin,df)
    (bin_value, df) = df[0], df[1]
    
    # calculate local FDR
    df = get_local_FDR(df, score_column, recom_data)
    
    # calculate peak FDR
    df = get_peak_FDR(df, score_column, col_Peak, closestpeak_column, recom_data)
    
    # calculate spire FDR
    #if recom_data: #recom_data =! 0
    #df = get_spire_FDR(df, score_column, col_Peak, recom_data)
    
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
    recom_data = config._sections['PeakFDRer']['recom_data']
    peak_label = config._sections['PeakAssignator']['peak_label']
    col_Peak = config._sections['PeakAssignator']['peak_column']
    col_CalDeltaMH = config._sections['PeakAssignator']['caldeltamh_column']
    closestpeak_column = config._sections['PeakAssignator']['closestpeak_column']
    # fdr_filter = config._sections['PeakFDRer']['fdr_filter']
    # target_filter = config._sections['PeakFDRer']['target_filter']
    
    # Read input file
    logging.info('Read input file')
    #df = pd.read_feather(args.infile)
    df = pd.read_csv(args.infile, sep="\t", float_precision='high')
    
    # Add groups
    groups = read_experiments(args.experiment_table)
    df = make_groups(df, groups)
    
    logging.info("Parallel the operations by BIN")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:        
        df = executor.map(bin_operations, list(df.groupby('bin')), repeat(score_column),
                                                                   repeat(recom_data), 
                                                                   repeat(peak_label),
                                                                   repeat(col_Peak),
                                                                   repeat(closestpeak_column)) 
    df = pd.concat(df)
    
    logging.info("Calculate gobal FDR")
    # df = get_global_FDR(df, score_column, recom_data)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        df = executor.map(get_global_FDR, list(df.groupby('Experiment')), repeat(score_column),
                                                                          repeat(recom_data))
    df = pd.concat(df)
    
    logging.info("Sort by DeltaMass cal")
    df.sort_values(by=[col_CalDeltaMH], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # TODO: groups?????
    
    # Filtering # This goes on a separate module now
    # df = filtering(df, fdr_filter, target_filter)
    # df.reset_index(drop=True, inplace=True)

    # d_h = df.head()
    # d_t = df.tail()
    # d_h.to_csv("kk_head.tsv", sep="\t")
    # d_t.to_csv("kk_tail.tsv", sep="\t")
    

    logging.info("Write output file")
    # https://towardsdatascience.com/the-best-format-to-save-pandas-data-414dca023e0d
    # begin:printHDF5
    # Note: Explote the Memory!!!
    # assign NumExpr for the tables module
    # tables.parameters.MAX_NUMEXPR_THREADS = args.n_workers
    # df.to_hdf('data.h5', key='df', mode='w')
    # end:printHDF5
    # df.to_csv('data.tsv', sep="\t", index=False)
    
    outfile = args.infile[:-4] + '_FDR.txt'
    #df.to_feather(outfile)
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
    

    

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
    parser.add_argument('-e',  '--experiment_table', required=True, help='Tab-separated file containing experiment names and file paths')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    parser.add_argument('-s',  '--score_column', help='Name of column with score for FDR calculation')
    #parser.add_argument('-f',  '--fdr_filter', help='FDR value to filter by')
    #parser.add_argument('-t',  '--target_filter', help='Filter targets, 0=no 1=yes')
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
    # if args.fdr_filter is not None:
    #     config.set('PeakFDRer', 'fdr_filter', str(args.fdr_filter))
    #     config.set('Logging', 'create_ini', '1')
    # if args.target_filter is not None:
    #     config.set('PeakFDRer', 'target_filter', str(args.target_filter))
    #     config.set('Logging', 'create_ini', '1')
    if args.recom_data is not None:
        config.set('PeakFDRer', 'recom_data', str(args.recom_data))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/PeakModeller.ini', 'w') as newconfig:
            config.write(newconfig)

    # logging debug level. By default, info level
    log_file = args.infile[:-4] + '_FDR_log.txt'
    log_file_debug = args.infile[:-4] + '_FDR_log_debug.txt'
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