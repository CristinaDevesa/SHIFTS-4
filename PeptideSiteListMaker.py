#!/usr/bin/env python
# coding: utf-8


# Module metadata variables
__author__ = "Cristina Amparo Devesa Arbiol"
__credits__ = ["Cristina Amparo Devesa Arbiol", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "cristinaamparo.devesa@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# Import modules
import pandas as pd
from os import remove
import numpy as np
from pandas import ExcelWriter
from optparse import OptionParser
import configparser
import argparse
import os
import logging
from pathlib import Path
import sys



###################
# Local functions #
###################
def readInfile(infile,cal_Dm_mh_colum_name):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t",                               
    float_precision='high',low_memory=False,dtype={cal_Dm_mh_colum_name:str})
    df[cal_Dm_mh_colum_name].astype("float64").dtypes
    return df





def peptideSiteListMaker(df,seq,counts):
    
    """
    PeptideSiteListMaker returns the frequency of the scan, the amino acid in the first position, cleaned sequence, and DM.
    """

    clean_seq = seq[:seq.find("[")]+seq[seq.find("]")+1:] # Clean sequence is obtained.
    aa = seq[seq.find("[")-1] # Aminoacid first position
    freq=counts.loc[seq] # Frecuency of the scan

    return clean_seq,seq,aa,freq



##################
# Main functions #
##################
def main(file,infile1):
    
    """
    Reading configuration file
    """

    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)
    logging.info("Reading PeptideSiteListMaker configuration file")  
    seq_column_name = config["PeptideSiteListMaker_Parameters"].get("sequence_column_name") # Sequence column name
    DM_column_name = config["PeptideSiteListMaker_Parameters"].get("DM_column_name") # DM column name
    
 
    logging.info("Processing input file")
    df = readInfile(infile1,DM_column_name) # A dataframe is created based on input file
    counts = df[seq_column_name].value_counts() # The number of times that the species appear is saved in the variable counts

    df2 = pd.DataFrame(columns=["EqSequence","Species","DM","AminoAcid","ScanFreq"],dtype=float) # Dataframe 2 is created with the aim of 

    cont = 0
    seqlist=[] # In this list it will be saved the sequences already analyzed
    for index, row in df.iterrows():

        if row[seq_column_name].find("_") == -1 and row[seq_column_name] not in seqlist:
            
            seqlist.append(row[seq_column_name])
            clean_seq,seq,aa,freq=peptideSiteListMaker(df,row[seq_column_name],counts)
            

            df2.loc[cont+1,"EqSequence"] = clean_seq
            df2.loc[cont+1,"Species"] = seq 
            df2.loc[cont+1,"DM"] = float(row[DM_column_name])
            df2.loc[cont+1,"AminoAcid"] = aa
            df2.loc[cont+1,"ScanFreq"] = float(freq)
            cont=cont+1


    logging.info("Writing output file")
    name = infile1[:-4]
    df2.to_csv(name+'__PeptideSiteListMade.txt', index=False, sep='\t', encoding='utf-8')
    
    logging.info('end script')



if __name__ == '__main__':
    
    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PeptideSiteListMaker',
        epilog='''
        Example:
            python PeptideSiteListMaker.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
        

    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '__PeptideSiteListMade_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '__PeptideSiteListMade_log_debug.txt'
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

    infile1=args.infile
    PeptideSiteListMakerini= args.config
    
    
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    
    # configuration files are read 
    try:
        open('Solver.ini',"r")
        solverini ='Solver.ini'
        logging.info("Modified Solver configuration file is going to be use")

    except:
        open("config/Solver.ini","r")
        solverini = "config/Solver.ini"

    main(solverini, infile1)
