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

#imports
from os import remove
import sys
from optparse import OptionParser
import configparser
import pandas as pd
import numpy as np
import argparse
import os
import logging
from pathlib import Path



###################
# Local functions #
###################

def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    
    df = pd.read_csv(infile, sep="\t", float_precision='high',low_memory=False)
    return df



def StickerSolver( Theo_mh,Exp_mh,seq,Error,dic_Mod):
    
    """
    StickerSolver function returs the sequence and label, that match with the smallest error
    according to the dic_Mod (dicctionary with all possible options).
    Input:
         Theo_mh             > Theorethical mh  
         Exp_mh              > Experimental mh
         seq                 > Sequence with DM
         Error               > Error in ppms stablished by the user
         dic_Mod             > Dicctionary  in which label is the key and the corresponding mass the value.

         
         
    The function returns:
        StickerLabel             > Label option that best fits
        StickerLable_ppm         > Error (ppm) of label selection
        StickerLable_description > Label description
    
    """
    

    
    # These variables are assigned for the cases in which there is no match.
    #StickerLabel = ""
    #StickerLable_ppm = ""
    #StickerLable_description = ""

   
    SequenceMassmod = seq[seq.find("[")+1:]
    SequenceMassmod = SequenceMassmod[:SequenceMassmod.find("]")]  # Mass modification of the sequence is obtained.
    seq = seq[:seq.find("[")]+seq[seq.find("]")+1:]  # Clean sequence is obtained.
  
                    
    # The dictionary option that obtains the smallest error is saved.
    minimun_DiffPPM = 2*10**10 
    for Label in dic_Mod:
        DiffPPM =abs(((Theo_mh+dic_Mod[Label][0]-Exp_mh)*1000000)/(Theo_mh+dic_Mod[Label][0]))
        if DiffPPM < minimun_DiffPPM:
            minimun_DiffPPM = DiffPPM 
            StickerLabel = Label
            StickerLabel_ppm = minimun_DiffPPM
            StickerLabel_description = dic_Mod[Label][1]


   
    return StickerLabel,StickerLabel_ppm, StickerLabel_description



##################
# Main functions #
##################

def main(solverconfig,infile,modsfileuser,modsfileUnimod):
    
    """
    Reading configuration file and processing file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(solverconfig) # Reading configuration file
    
    logging.info("Reading Sticker configuration file")
    
    Error = config["Sticker_Parameters"].getfloat("Relative_Error_ppm") # Relative error (ppm)
    Exp_colum_name = config["Sticker_Parameters"].get("exp_mh_colum_name") # Experimental mh  column name
    Theo_colum_name = config["Sticker_Parameters"].get("theo_mh_colum_name")# Theorethical mh colum name 
    Seq_colum_name = config["Sticker_Parameters"].get("Sequence_colum_name") # Sequence colum name

    StickerLabel_User_output_column_name = config["Sticker_Parameters"].get("StickerLabel_User_output_column_name") # Relative error (ppm)# Column name of the output where the chosen label (from  User file) is annotated
    StickerLabel_Unimod_output_column_name =  config["Sticker_Parameters"].get("StickerLabel_Unimod_output_column_name")# Column name of the output where the chosen label (from  Unimod file) is annotated 
    StickerLabel_ppm_User_output_column_name =  config["Sticker_Parameters"].get("StickerLabel_ppm_User_output_column_name")# Column name of the output where the calculated error in ppm for the selected label (from  User file) is annotated 
    StickerLabel_ppm_Unimod_output_column_name =  config["Sticker_Parameters"].get("StickerLabel_ppm_Unimod_output_column_name")# Column name of the output where the calcçulated error in ppm for the selected label (from  Unimod file) is annotated 
    StickerLabel_Description_output_column_name = config["Sticker_Parameters"].get("StickerLabel_Description_output_column_name")	# Column name of the output where the description of the label is annotated
    output_file_suffix = config["Sticker_Parameters"].get("output_file_suffix") # Chosen suffix for output file


    # A dicctionary that save the labels and their corresponding masses is created
    modfileo=open(modsfileuser,"r") # User file
    dic_Mod = {} 
    next(modfileo)
    for line in modfileo:
        line=line.strip("\n")
        fields=line.split("\t")
        dic_Mod[" "+fields[1]] = float(fields[0].replace(",",".")),""
    
    modfileounimod=open(modsfileUnimod,"r") # Unimod file
    dic_Mod_unimod = {} 
    next(modfileounimod)
    for line in modfileounimod:
        line=line.strip("\n")
        fields=line.split("\t")
        dic_Mod_unimod[" "+fields[1]] = float(fields[0].replace(",",".")),fields[2]
            
            
    
    # Input file is read as data frame and colums names are added to the header
    df = readInfile(infile)
    try: 
        df.drop([StickerLabel_User_output_column_name,StickerLabel_ppm_User_output_column_name,StickerLabel_Unimod_output_column_name,StickerLabel_ppm_Unimod_output_column_name,StickerLabel_Description_output_column_name], axis=1)
    except:
        pass
    df[StickerLabel_User_output_column_name] = ""
    df[StickerLabel_ppm_User_output_column_name] = np.nan
    df[StickerLabel_Unimod_output_column_name] = ""
    df[StickerLabel_ppm_Unimod_output_column_name] = np.nan
    df[StickerLabel_Description_output_column_name] = ""

    cont = 0
    logging.info("Processing input file")
    for index, row in df.iterrows():     
        StickerLabel,StickerLabel_ppm, StickerLabel_description= StickerSolver(row[Theo_colum_name],row[Exp_colum_name],row[Seq_colum_name],Error,dic_Mod)
        StickerLabel_unimod,StickerLabel_ppm_unimod, StickerLabel_description_unimod = StickerSolver(row[Theo_colum_name],row[Exp_colum_name],row[Seq_colum_name],Error,dic_Mod_unimod)
       
            
        # New columns are completed  
        df.loc[cont,StickerLabel_User_output_column_name] = StickerLabel
        df.loc[cont,StickerLabel_ppm_User_output_column_name] = StickerLabel_ppm
        df.loc[cont,StickerLabel_Unimod_output_column_name] = StickerLabel_unimod
        df.loc[cont,StickerLabel_ppm_Unimod_output_column_name] = StickerLabel_ppm_unimod
        df.loc[cont,StickerLabel_Description_output_column_name] = StickerLabel_description_unimod
        cont=cont+1

    # write outputfile
    logging.info("Writing output file")

    outfilename = infile[:-4]+output_file_suffix

    outfile= outfilename+".txt"
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')
  

    logging.info('end script')


        





if __name__ == '__main__':
    
    try:
        remove('Solver.ini')
    except:
        None

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Sticker',
        epilog='''
        Example:
            python Sticker.py
        ''')
      
    # Default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-m', '--modsfile', required=True, help='Path to modifications file')
    parser.add_argument('-u', '--unimodmodsfile', required=True, help='Path to Unimod  file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # These will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    
    if args.relerror is not None:
        config.set('DM0Solver_Parameters', 'Relative_Error', str(args.relerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + '_Sticker_log.txt'
    log_file_debug = outfile = args.infile[:-4] + '_Sticker_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level = logging.DEBUG,
                            format = '%(asctime)s - %(levelname)s - %(message)s',
                            datefmt = '%m/%d/%Y %I:%M:%S %p',
                            handlers = [logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level = logging.INFO,
                            format = '%(asctime)s - %(levelname)s - %(message)s',
                            datefmt = '%m/%d/%Y %I:%M:%S %p',
                            handlers = [logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    infile1 = args.infile
    modsfile = args.modsfile
    unimodmodsfile = args.unimodmodsfile

    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    
    # Configuration files are read   
    try:
        open('Solver.ini',"r")
        solverini ='Solver.ini'
        logging.info("Modified Solver configuration file is going to be use")


    except:
        open("config/Solver.ini","r")
        solverini = "config/Solver.ini"
    main(solverini, infile1,modsfile,unimodmodsfile)
