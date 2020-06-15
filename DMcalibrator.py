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
from pathlib import Path
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
# Local functions #
###################

# Calibrate mass separately for each raw file.

def readInfile(infile):
    '''    
    Read input file to dataframe and determine if it is Comet or Recom.
    '''
    df = pd.read_csv(infile, skiprows=1, sep="\t", float_precision='high')
    return df

def getTheoMZ(df):
    '''    
    Calculate theoretical MZ using the PSM sequence.
    '''
    AAs = dict(config._sections['Aminoacids'])
    if 'theo_mz' not in df:
        df.insert(df.columns.get_loc(config._sections['Input']['mzcolumn'])+1, 'theo_mz', np.nan)
    
    def _PSMtoMZ(sequence, charge):
        total_aas = 2*config._sections['Masses']['M_proton'] + config._sections['Masses']['M_oxygen']
        for aa in sequence:
            if aa.lower() in AAs:
                total_aas += AAs[aa]
            #else: # aminoacid not in list (ask for user input?)
                # TODO
        MZ = (total_aas + charge*config._sections['Masses']['M_proton']) / charge
        return MZ
    
    df['theo_mz'] = df.apply(lambda x: _PSMtoMZ(x[config._sections['Input']['seqcolumn']], x[config._sections['Input']['zcolumn']]), axis = 1)
    return df

def getErrors(df):
    '''    
    Calculate absolute (in m/z) and relative (in ppm) errors.
    '''
    if 'abs_error' not in df:
        df.insert(df.columns.get_loc('theo_mz')+1, 'abs_error', np.nan)
    if 'rel_error' not in df:
        df.insert(df.columns.get_loc('abs_error')+1, 'rel_error', np.nan)
    df['abs_error'] = df[config._sections['Input']['mzcolumn']] - df['theo_mz']
    df['rel_error'] = df['abs_error'] / df['theo_mz'] * 1e6
    return df

def filterPeptides(df, scoremin, ppmmax, scorecolumn, chargecolumn, mzcolumn, seqcolumn):
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
    
    #keep targets
    df_filtered = df[~df[config._sections['Input']['proteincolumn']]
                     .str.startswith(config._sections['Input']['decoyprefix'])]
    #keep score > scoremin
    df_filtered = df_filtered[df_filtered[config._sections['Input']['scorecolumn']]
                              >=scoremin]
    #keep abs_error <= ppmmax
    df_filtered = df_filtered[df_filtered['abs_error']
                              <=ppmmax]

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
    if 'exp_mz_cal' not in df:
        df.insert(df.columns.get_loc(config._sections['Input']['mzcolumn'])+1, 'exp_mz_cal', np.nan)
    df['exp_mz_cal'] = df[config._sections['Input']['mzcolumn']] - sys_error
    return df

def getDMcal(df, dmcolumn):
    '''
    Calculate calibrated DM values.
    '''
    if 'exp_dm_cal' not in df:
        df.insert(df.columns.get_loc(config._sections['Input']['mzcolumn'])+1,
                  'exp_dm_cal',
                  np.nan)
    df['exp_dm_cal'] = df[config._sections['Input']['mzcolumn']] - df['theo_mz']
    return df


#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
      
    log_str = "Calibrating file: " + str(Path(args.infile))
    logging.info(log_str)
    # Read infile
    df = readInfile(Path(args.infile))
    # Calculate theoretical MZ
    df = getTheoMZ(df)
    # Calculate errors
    df = getErrors(df)
    # Filter identifications
    df_filtered = filterPeptides(df,
                                 config._sections['Filtering']['score_min'],
                                 config._sections['Filtering']['ppm_max'],
                                 config._sections['Input']['scorecolumn'],
                                 config._sections['Input']['zcolumn'],
                                 config._sections['Input']['mzcolumn'],
                                 config._sections['Input']['seqcolumn'])
    # Use filtered set to calculate systematic error
    sys_error, avg_ppm_error = getSysError(df_filtered)
    # Use systematic error to correct infile
    df = rawCorrection(df, sys_error)
    # Calculate DMCal 
    df = getDMcal(df, config._sections['Input']['dmcolumn'])
    #Write to txt file
    logging.info("Writing output file...")
    outfile = args.infile[:-4] + '_calibrated.txt'
    df.to_csv(outfile, index=False, encoding='utf-8')
    logging.info("Calibration finished")

    
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
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
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
        config._sections['Logging']['create_INI'] = 1
    if args.ppmmax is not None:
        config._sections['Filtering']['ppm_max'] = args.ppmmax
        config._sections['Logging']['create_INI'] = 1
    if args.scorecolumn is not None:
        config._sections['Input']['scorecolumn'] = args.scorecolumn
        config._sections['Logging']['create_INI'] = 1
    if args.mzcolumn is not None:
        config._sections['Input']['mzcolumn'] = args.mzcolumn
        config._sections['Logging']['create_INI'] = 1
    if args.zcolumn is not None:
        config._sections['Input']['zcolumn'] = args.zcolumn
        config._sections['Logging']['create_INI'] = 1
    if args.seqcolumn is not None:
        config._sections['Input']['seqcolumn'] = args.seqcolumn
        config._sections['Logging']['create_INI'] = 1
    # if something is changed, write a copy of ini
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