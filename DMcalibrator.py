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
import configparser
import argparse
import logging
import math
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

#infile = r"C:\Users\Andrea\Desktop\SHIFTS-4\testing\recom.txt"
# os.chdir(r"C:\Users\Andrea\Desktop\SHIFTS-4")

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
       'O': 132.089878}
M_proton = 1.007825
M_oxygen = 15.994915

###################
# Local functions #
###################

# Calibrate mass separately for each raw file.

def readInfile(infile):
    '''    
    Read input file to dataframe and determine if it is Comet or Recom.
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
    df.insert(df.columns.get_loc('exp_mz')+1, 'theo_mz', np.nan)
    
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

def getErrors(df):
    '''    
    Calculate absolute (in m/z) and relative (in ppm) errors.
    '''
    df.insert(df.columns.get_loc('theo_mz')+1, 'abs_error', np.nan)
    df.insert(df.columns.get_loc('abs_error')+1, 'rel_error', np.nan)
    df['abs_error'] = df['exp_mz'] - df['theo_mz']
    df['rel_error'] = df['abs_error'] / df['theo_mz'] * 1e6
    return df

def filterPeptides(df, recom, scoremin, ppmmax, scorecolumn, chargecolumn, mzcolumn, seqcolumn):
    '''    
    Filter and keep target peptides that match Xcorrmin and PPMmax conditions.
    This high-quality subpopulation will be used for calibration.
    '''
    
    def _correctXcorr(charge, xcorr, length):
        if charge < 3:
            cxcorr = math.log10(xcorr) / math.log10(2*length)
        else:
            cxcorr = math.log10(xcorr/1.22) / math.log10(2*length)
        return cxcorr
    
    if recom == 0: #non-recom input
        #calculate cxcorr
        df.insert(df.columns.get_loc('xcorr')+1, 'xcorr_corr', np.nan)
        df['xcorr_corr'] = df.apply(lambda x: _correctXcorr(x['charge'],
                                                            x['xcorr'],
                                                            len(x['plain_peptide'])),
                                    axis = 1)
        #keep targets
        df_filtered = df[~df['protein'].str.startswith('DECOY')]
        #keep score > scoremin
        df_filtered = df_filtered[df_filtered['xcorr_corr']>=scoremin]
        #keep abs_error <= ppmmax
        df_filtered = df_filtered[df_filtered['abs_error']<=ppmmax]
    else: #recom input
        #make best_cxcorr column
        df.insert(df.columns.get_loc('Best_Xcorr')+1, 'Best_cXcorr', np.nan)
        df['Best_cXcorr'] = df.apply(lambda x: x['xcorr_corr'] if (x['xcorr_corr']>x['Closest_Xcorr_corr']) else x['Closest_Xcorr_corr'], axis = 1)
        
        #keep targets
        df_filtered = df[~df['protein'].str.startswith('DECOY')]
        #keep cxcorr > scoremin
        df_filtered = df_filtered[df_filtered['Best_cXcorr']>=scoremin]
        #keep abs_error <= ppmmax
        df_filtered = df_filtered[df_filtered['abs_error']<=ppmmax]
    return df_filtered

def getSysError(df_filtered):
    '''
    Calculate systematic error and average PPM error.
    '''
    sys_error = df_filtered['abs_error'].median()
    
    phi = math.sqrt(2) * math.erf(-1)
    mad = df_filtered['abs_error'].mad()
    avg_ppm_error = mad / phi
    return sys_error, avg_ppm_error

def rawCorrection(df, sys_error):
    '''
    Correct exp_mz values from infile using the systematic error.
    '''
    df.insert(df.columns.get_loc('exp_mz')+1, 'exp_mz_cal', np.nan)
    df['exp_mz_cal'] = df['exp_mz'] - sys_error
    return df

def getDMcal(df):
    '''
    Calculate calibrated DM values.
    '''
    df.insert(df.columns.get_loc('exp_mz_cal')+1, 'exp_dm_cal', np.nan)
    df['exp_dm_cal'] = df['exp_mz_cal'] - df['theo_mz']
    #if recom == 0: #comet input
        #  se calcularían a partir de ExpMZCal y TeorMZs para comet
    #else: #recom input #actually we can just handle this later
    return df


#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    logging.info("read input file list")
    with open(args.infile) as f:
        infile_list = f.readlines()
    infile_list = [x.strip() for x in infile_list] # remove whitespace
    
    #TODO: parallelize?
    # with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:
    
    log_str = str(len(infile_list)) + " file(s) to calibrate..."
    logging.info(log_str)
    i = 0
    for infile in infile_list: 
        i += 1
        log_str = "calibrating file " + str(i) + " of " + str(len(infile_list))
        logging.info(log_str)
        # Read infile
        df, recom = readInfile(infile)
        # Calculate theoretical MZ
        df = getTheoMZ(df)
        # Calculate errors
        df = getErrors(df)
        # Filter identifications
        df_filtered = filterPeptides(df,
                                     recom,
                                     float(config._sections['Filtering']['score_min']),
                                     float(config._sections['Filtering']['ppm_max']),
                                     int(config._sections['Input']['scorecolumn']),
                                     int(config._sections['Input']['zcolumn']),
                                     int(config._sections['Input']['mzcolumn']),
                                     int(config._sections['Input']['seqcolumn']))
        # Use filtered set to calculate systematic error
        sys_error, avg_ppm_error = getSysError(df_filtered)
        # Use systematic error to correct infile
        df = rawCorrection(df, sys_error)
        # Calculate DMCal 
        df = getDMcal(df)
        #Write to txt file
        outfile = infile[:-4] + '_calibrated.txt'
        df.to_csv(outfile, index=False, encoding='utf-8')

    
if __name__ == '__main__':

    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='DMcalibrator',
        epilog='''
        Example:
            python DMcalibrator.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/DMcalibrator.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='List of input files from COMET or RECOM')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # these will overwrite the config if specified
    parser.add_argument('-s', '--scoremin', default=None, help='Minimum score')
    parser.add_argument('-p', '--ppmmax', default=None, help='Maximum PPM error')
    parser.add_argument('-sc', '--scorecolumn', default=None, help='Position of the column containing the score')
    parser.add_argument('-zc', '--chargecolumn', default=None, help='Position of the column containing the charge')
    parser.add_argument('-mc', '--mzcolumn', default=None, help='Position of the column containing the experimental m/z')

    parser.add_argument('-w', '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.scoremin is not None:
        config._sections['Filtering']['score_min'] = args.scoremin
    if args.ppmmax is not None:
        config._sections['Filtering']['ppm_max'] = args.ppmmax
    if args.scorecolumn is not None:
        config._sections['Input']['scorecolumn'] = args.scorecolumn
    if args.mzcolumn is not None:
        config._sections['Input']['mzcolumn'] = args.mzcolumn
    if args.zcolumn is not None:
        config._sections['Input']['zcolumn'] = args.zcolumn
    if args.seqcolumn is not None:
        config._sections['Input']['seqcolumn'] = args.seqcolumn
    # TODO: if something is changed, write a copy of ini
        

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