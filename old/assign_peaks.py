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
import os
import sys
import argparse
import logging
import re
import pandas as pd
import concurrent.futures
from itertools import repeat
import tables
import numpy as np

#feather_data = pd.read_feather(r"C:\Users\Andrea\Desktop\data.ftr")
#apex_list = _extract_ApexList(r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\bin002_w7_SL4500Target_ApexList_mass.txt")

###################
# Local functions #
###################
def concatInfiles(infile, fwhm_fname):
    '''    
    Concat input files...
    adding Experiment column (dirname of input file), and adding a FWHM columns by Experiment
    '''
    def _extract_FWHM(file):
        with open(file) as f:
            data = f.read()
            m = re.findall(r'FWHM:\s*([^\n]*)', data)
            if m and len(m)>0:
                return m[0]
            else:
                sys.exit("ERROR! FWHM is not defined for {}".format(file))
        
    # read input file
    # use high precision with the floats
    df = pd.read_csv(infile, sep="\t", float_precision='high')    
    # add folder name into column
    foldername = os.path.dirname(infile)
    df['Experiment'] = foldername
    # add fwhm column
    fwhm_file = "{}/{}".format(foldername, fwhm_fname)
    fwhm = _extract_FWHM(fwhm_file)
    df['FWHM'] = float(fwhm)
    # assign type to categorical columns
    df['Experiment'] = df['Experiment'].astype('category')
    df['filename'] = df['filename'].astype('category')
    df['Label'] = df['Label'].astype('category')
    df['IsotpicJump'] = df['IsotpicJump'].astype('category')
    return df

def closest_peak(apex_list, delta_MH):
    '''
    Assign a delta_MH value to the closest apex in a list
    '''
    peak = min(apex_list, key = lambda x : abs(x - delta_MH))
    return peak

    
def find_orphans(nsigma, fwhm, peak, delta_MH):
    '''
    Identify orphans and peaks
    '''
    # window = float(nsigma) * fwhm / 2
    distance = abs(peak - delta_MH)
    max_distance = abs(float(nsigma) * fwhm / 2)
    if distance <= max_distance:
        ID = 'PEAK' # better to use the actual value?
    else:
        ID = 'ORPHAN'
    return ID

def get_spire_FDR(df, xcorr_type):
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
    recom_peaks = peaks[peaks['XcorType'] == 'RECOM'] # XcorType needs to be created
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

def get_peak_FDR(df, xcorr_type):
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
        if xcorr_type == 0: # by Comet Xcorr
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

def get_local_FDR(df, xcorr_type):
    '''
    Calculate local FDR for one bin (1 Da)
    '''
    # sort bin
    if xcorr_type == 0: # by Comet Xcorr
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

def get_global_FDR(df, xcorr_type):
    '''
    Calculate global FDR
    '''
    # sort by score
    if xcorr_type == 0: # by Comet Xcorr
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

def bin_operations(df, apex_list, nsigma, xcorr_type):
    '''
    Main function that handles the operations by BIN
    '''
    # get the BIN value from the input tuple df=(bin,df)
    (bin_value, df) = int(df[0]), df[1]
    
    # assign to peaks
    df['ClosestPeak'] = df.apply(lambda x: closest_peak(apex_list, x['Cal_Delta_MH']), axis = 1)

    # identify orphans
    df['PeakAssignation'] = df.apply(lambda x: find_orphans(nsigma, x['FWHM'], x['ClosestPeak'], x['Cal_Delta_MH']), axis = 1)
    df['PeakAssignation'] = df['PeakAssignation'].astype('category')
    
    # calculate local FDR
    df = get_local_FDR(df, xcorr_type)
    
    # calculate peak FDR
    df = get_peak_FDR(df, xcorr_type)
    
    # calculate spire FDR
    df = get_spire_FDR(df, xcorr_type)
    
    # create deltamass column
    df['deltaMass'] = df.apply(lambda x: x['deltaPeptide'] if (x['PeakAssignation']=='ORPHAN') else x['ClosestPeak'])
    
    # def peak_FDR():
      # for each peak sort by xcorr (comet) # should we separate recom peaks?
      # for each peak separate targets and decoys
      # make new fdr column rank_D/rank/T
    # def local_FDR():
      # for each bin sort by xcorr (comet)
      # for each bin separate targets and decoys
      # make new fdr_column rank_D/rank_T
    #def recom_FDR():
      # for each recom peak sort by xcorr (recom)
    
    # TO CHECK, print by BIN
    # outfile = os.path.join("D:/tmp/kk/", bin+"_kk.tsv")
    # df.to_csv(outfile, sep="\t", index=False)
    return df  
    
#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    # main variables
    col_CalDeltaMH = 'Cal_Delta_MH'

    
    logging.info("get the list of files with the inputs")
    with open(args.infile) as f:
        infiles = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    infiles = [x.strip() for x in infiles] 
    logging.debug(infiles)
    
     # read apex list
    def _extract_ApexList(file):
        with open(file) as f:
            data = f.read().split('\n')
            data = [x for x in data if x.strip()]
            data = np.array(data, dtype=np.float64)
            return data
    
    foldername = os.path.dirname(args.appfile)
    apex_file = "{}/{}".format(foldername, args.appfile)
    apex_list = _extract_ApexList(apex_file)


    logging.info("concat input files...")
    logging.info("adding Experiment column (dirname of input file),")
    logging.info("and adding a FWHM column by Experiment")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
        df = executor.map(concatInfiles, infiles, repeat(args.fwhm_filename))
    df = pd.concat(df)
    df.reset_index(drop=True, inplace=True)
 
    logging.info("create a column with the bin")
    df['bin'] = df[col_CalDeltaMH].astype(str).str.extract(r'^([^\.]*)')


    logging.info("parallel the operations by BIN")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:        
        df = executor.map(bin_operations, list(df.groupby("bin")), repeat(apex_list),
                                                                   repeat(args.nsigma),
                                                                   repeat(args.xcorr))
    df = pd.concat(df)
    logging.info("calculate gobal FDR")
    df = get_global_FDR(df, args.xcorr)
    logging.info("sort by DeltaMax cal")
    df.sort_values(by=[col_CalDeltaMH], inplace=True)
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
        description='Assign peaks',
        epilog='''
        Example:
            python assign_peaks.py

        ''')
    parser.add_argument('-i',  '--infile', required=True, help='Input file with the list of files that contains the peak picking')
    parser.add_argument('-a',  '--appfile', required=True, help='File with the apex list of Mass')
    
    parser.add_argument('-f',  '--fwhm_filename', default='MAD_and_FWHM_calculations.txt', help='File name with the FWHM value (default: %(default)s)')    
    parser.add_argument('-mn', '--minDelta', default=-500, help='Minimum Delta Mass (default: %(default)s)')
    parser.add_argument('-mx', '--maxDelta', default=500, help='Maximum Delta Mass (default: %(default)s)')
    parser.add_argument('-s',  '--nsigma', default=3, help='Coefficient of Sigma (default: %(default)s)')
    parser.add_argument('-x',  '--xcorr', default=1, help='Score for FDR calculation: 0=Xcorr, 1=cXcorr (default: %(default)s)')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
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