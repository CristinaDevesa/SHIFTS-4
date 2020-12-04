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
import re
import sys
import math
import time
from pathlib import Path
from optparse import OptionParser
import requests
import configparser
from Bio import SeqIO
import numpy as np 
import pandas as pd 
from pandas import ExcelWriter
import argparse
import os
import logging
from pathlib import Path
import tkinter as tk


###################
# Local functions #
###################
def readInfile(infile):
    '''    
    Read input file to dataframe.
    '''
    df = pd.read_csv(infile, sep="\t", float_precision='high', low_memory=False)
    return df


def theoretical_mh_by_hand(subseq,label_mass,dic_mod,dic_aa,selectedaa,Mproton,Hydrogen,O2,decnum):
    
    """
    Theoretical mass is calculated taking into account fix modifications, label and subsequence. This functions returns theoretical
    mass fix modifications positions and the subsequence adding fix modifications.
    """
    decnum=int(decnum.replace("f","").replace(".",""))
    H20 = 2*Hydrogen+O2
    pattern = "[a-z]"
    
    total = 0
    newsequence = []
    c = 0
    mods_position = []
    
    # For each amino acid, depending on which is its condition (fix modifications, label, N-termina), mass is added  
    for aa in subseq:
        c = c+1
        
        if re.search(pattern,aa) != None: # If aa in N-terminal position  
            aa = aa.upper()
            
            if aa in dic_mod.keys(): # If it has a fix modification
                total = total+float(dic_mod[aa][1])
                total = total+label_mass
                total = round(total,decnum)
                newsequence.append(dic_mod[aa][0]+"TMT")
                mods_position.append(str(c)+"_"+aa+"_"+str(label_mass)+"_N")
                mods_position.append(str(c)+"_S_"+str(dic_mod[aa][2]))  
                    
            else: # If it has not a fix modification
          
                total = total+float(dic_aa[aa]) 
                total = total+label_mass
                total = round(total,decnum)
                cosa = float(dic_aa[aa])+float(label_mass)
                newsequence.append(aa+"-TMT")
                mods_position.append(str(c)+"_"+aa+"_"+str(label_mass)+"_N")
                
        else: # If it has not N-terminal position
         
            if aa in dic_mod.keys(): # If it has a fix modification
    
                total = total+float(dic_mod[aa][1])
                total = round(total,decnum)
                newsequence.append(dic_mod[aa][0])
                mods_position.append(str(c)+"_S_"+str(dic_mod[aa][2]))
   
                
            else: # If it has not a fix modification
                
                total = total+dic_aa[aa]
                total = round(total,decnum)
                newsequence.append(aa)
              
                
    mods_position = ",".join(mods_position)
    total = total+H20+Mproton 
    total = round(total,decnum)

    return total,newsequence,mods_position



def tag(seq,subseq):
    
    """
    Tag function returns a variable indicating if the cut is tryptic or not taking into account the 
    sequence plus one extra amino acid at both sides
    """

    Truncation = []
    number = subseq.find(seq)
    number2 = number+len(seq)
    loss = subseq.replace(seq,"")
    left = ""
    right = ""
    
    # Take into account if there is amino acid at both ends or just at one 
    if number2-number <= len(seq) and number != 0:
        left = subseq[0:number]
        right = subseq[number2:] 
    elif number != 0:
        left = loss    
    elif number2-number <= len(seq):
        right = loss
          
    
    # If the cut is tryptic or not is determined at both sides
    if left != "":
        if left[-1] == "K" or left[-1] == "R" and seq[0] != "P":   
            Truncation.append("YeS")
        else:
             Truncation.append("No") 
            
    if right != "":
        if seq[-1] == "K" or seq[-1] == "R" and right[0] != "P":
             Truncation.append("YeS")
        else:
             Truncation.append("No")
           
          
    # if one of the ends have a non tryptic cut means that there is a truncation
    if "No" in Truncation:
        Truncation1="No"
    else:
        Truncation1="YeS"

    return Truncation1



def Obtain_values(seq,MasterProtein_column,dic_fasta):
    
    """
    Taking in to account sequence, master protein and fasta dictionary of the input file this function returns 
    the entire sequence corresponding to the  master protein and all the initial and final positions that match with the 
    selected sequence (seq).
    """
    
    clean_seq = seq[:seq.find("[")]+seq[seq.find("]")+1:].upper() #The clean sequence is obtained.
    MasterProtein=MasterProtein_column.strip("\n").split("_")
    MasterProtein= MasterProtein[0]+"_"+MasterProtein[1] # The id is extracted from the Master Protein name 
    
    # The fasta sequence corresponding to this identifier is saved 
    for iden in dic_fasta:
        if MasterProtein==iden:
            result=str(dic_fasta[iden].seq.upper()).replace("X","L")
            break
   
    pattern=re.compile(clean_seq.replace("L","l").replace("I","[IL]").replace("l","[IL]")) # Problems that may exist with leucine and isoleucine are solved
    
    dic_seqs={}
    pos = 0
    
    # The corresponding fasta sequence is rigorously scrutinized so that no chance is missed  
    while True:
        match = pattern.search(result, pos)
        if not match:
            break
        s = match.start()
        e = match.end()
        s = s+1
        if s-1 == 0:
            pos1 = 0
        else:
            pos1 = s-2
            
        try:
            p2 = result[e+1]
            pos2 = e+1
        except:
            pos2 = e
    
        s=s-1
        
        cut = tag(result[s:e],result[pos1:pos2])
        dic_seqs[str(s)+":"+str(e)]=cut
        pos = e #  Move forward in text for the next search
    
    # If the same sequence is found in the same protein sequence those that have a tryptic cut will be selected     
    if "YeS" in dic_seqs.values():
        for key in list(dic_seqs):
            if dic_seqs[key] == "No":
                dic_seqs.pop(key)  
           
    return dic_seqs,result





def best_combination(subseq,Exp_Mh,cont,j1,j2,Error,label_mass,dic_mod,selectedaa,dic_CombList,dic_aa,decnum,Mproton,Hydrogen,O2):
    """
    Best_combiantions function returns the combinations of "Combination list" that give rise to an error less than or equal to the allowed and 
    varibales that indicate  whether TrunKSolver should stop extending the length of the sequence being analyzed
    """
    # New subsequence mass is calculated by theoretical_mh_by_hand function
    ther,newsequence,mods_position = theoretical_mh_by_hand(subseq,label_mass,dic_mod,dic_aa,selectedaa,Mproton,Hydrogen,O2,decnum) 
    
    # initial parameters are set 
    Error2=100
    ngreater = 0
    minimun_DiffMasa = Error2
    minimun_DiffPPM = Error
    TrunkSequence = ""
    TrunkDM = ""
    TrunkLabel = ""
    Trunk_Label_ppm = ""
    DiffMasa2 = ""
    New_DM = ""
    New_Theo_MH = ""
    
    
    DiffMasa2 = Exp_Mh-ther # Mass difference betwen experimental an theorethical mass 
   
    # All posible options of "Combination List" are examined
    number_option = 0   
    for iden in dic_CombList: 
        number_option = number_option+1     
        mass = dic_CombList[iden]      
        total_value2 = ther+mass
        total_value2 = float(format (total_value2,decnum))
        DiffMasa = (total_value2-Exp_Mh) #Mass difference between experimental and theorical plus one option of "Combination List"     
       
        # Mass difference is calculated in ppms  
        DiffPPM = abs(((total_value2-Exp_Mh)*1000000)/total_value2)
  
        # The lowest DiffPPM value is saved 
        if abs(DiffPPM) <= Error:
            if DiffPPM < minimun_DiffPPM:
                minimun_DiffMasa = DiffMasa2
                minimun_DiffPPM = DiffPPM
                DiffMasa_sequence = format( minimun_DiffMasa,decnum)
                TrunkSequence = (("".join(subseq)).upper()+"_"+str(DiffMasa_sequence))
                TrunkDM = minimun_DiffMasa
                TrunkLabel = iden
                Trunk_Label_ppm = minimun_DiffPPM
                New_DM = (Exp_Mh-total_value2)
                New_Theo_MH = total_value2 

        elif total_value2 > Exp_Mh+Error2:
            ngreater = ngreater+1
      
    # If all possible option of the number options are greater than the experimental mass plus the Error TrunkSolver function 
    # will stop extending the sequence at that end
    if ngreater == number_option:     
        if cont == 2:
            j2 = "subseq2_stop"
        elif cont == 1:
            j1 = "subseq1_stop"

    return minimun_DiffPPM,TrunkSequence,TrunkDM,TrunkLabel,mods_position,Trunk_Label_ppm,j2,j1,New_DM, New_Theo_MH





def TrunkSolver(seq,dic_seqs,Exp_mh,calibrated_delta_MH,result,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,Theo_mh):
    
    addition = "" # Variable cretaed with the aim of controlling if there are more than one possible subsequence  
    
    # All possible sequence found in protein sequence will be analyzed 
    for key in dic_seqs:
        minimun = 2000000 # initial parameter
        fields = key.split(":")
        position = int(fields[0])
        final_position = int(fields[1])
  
        result=result.upper()
        result1 = list(result) 
        deltapeptide = seq[:seq.find("[")]+seq[seq.find("]")+1:]

        i = 1
        k = 1
        j1 = ""
        j2 = ""
        dic_r = {}

        # While the extension of the subsequence can be elongated one aminoacid will be added 
        while j1 != "subseq1_stop" or j2 != "subseq2_stop":
            listseqs = [] # all possible subsequences will be saved in this list

            if j1 != "subseq1_stop":
                if position+i <= len(result1[position:])+position:
                    subseq1 = result1[position:position+i+1]
                    subseq1[0] = subseq1[0].lower()
                    listseqs.append(subseq1)
                    i = i+1
                    cont = "1"
                else:
                    j1 = "subseq1_stop"


            if j2 != "subseq2_stop":
                if final_position-k >= 0:
                    subseq2 = result1[final_position-k:final_position]
                    subseq2[0] = subseq2[0].lower()
                    listseqs.append(subseq2)
                    k = k+1
                    cont = "2"
                else:
                    j2 = "subseq2_stop"

           
            if len(listseqs) > 0:
                
                for subseq in listseqs: # best combination of this subsequence is saved
                    minimun_DiffPPM,TrunkSequence,TrunkDM,TrunkLabel,mods_position,Trunk_Label_ppm,j2,j1,New_DM, New_Theo_MH = best_combination(subseq,Exp_mh,int(cont),j1,j2,Error,float(NT_label),dic_mod,selectedaa,dic_CombList,dic_aa,decnum,Mproton,Hydrogen,O2)
                    
                    if TrunkSequence != "": # the best options of all the subsquences are saved in a dicctionary 
                        
                        # If the PPM difference is greater than the minimun the option will be discarded
                        if abs(minimun_DiffPPM) <= minimun:
                            minimun = abs(minimun_DiffPPM)
                            final_TrunkSequence1 = TrunkSequence
                            final_TrunkSequence2 = final_TrunkSequence1[:final_TrunkSequence1.find("_")]
                            ptag = result[result.find(final_TrunkSequence2)-1:result.find(final_TrunkSequence2)+len(final_TrunkSequence2)+1]
                    
                            if deltapeptide==final_TrunkSequence2: 
                                cut=" "
                            else:
                                cut = tag(final_TrunkSequence2,ptag)
                            dic_r[minimun_DiffPPM] = TrunkDM,TrunkSequence,TrunkLabel,mods_position,Trunk_Label_ppm,cut,New_DM, New_Theo_MH


    # if there is more than one possibility those that have a tryptic digestion will have preference
    if dic_r :
        for key in dic_r:
            if str(dic_r.values()).find("YeS") != -1:
                if dic_r[key][5] == "YeS":
                    minimun = abs(key)
                    final_TrunkDM = dic_r[key][0]
                    final_TrunkSequence = dic_r[key][1]
                    final_TrunkLabel = "TrypticCut;"+dic_r[key][2]
                    final_mods_position = dic_r[key][3]
                    final_Trunk_Label_ppm = dic_r[key][4]
                    final_New_DM = dic_r[key][6]
                    final_New_Theo_MH = dic_r[key][7]
            elif dic_r[key][5] == "No":                
                minimun = abs(key)
                final_TrunkDM = dic_r[key][0]
                final_TrunkSequence = dic_r[key][1]
                final_TrunkLabel = "Truncation;"+dic_r[key][2]
                final_mods_position = dic_r[key][3]
                final_Trunk_Label_ppm = dic_r[key][4]
                final_New_DM = dic_r[key][6]
                final_New_Theo_MH = dic_r[key][7]
            elif dic_r[key][5] == " ":                
                minimun = abs(key)
                final_TrunkDM = dic_r[key][0]
                final_TrunkSequence = dic_r[key][1]
                final_TrunkLabel = dic_r[key][2]
                final_mods_position = dic_r[key][3]
                final_Trunk_Label_ppm = dic_r[key][4]
                final_New_DM = dic_r[key][6]
                final_New_Theo_MH = dic_r[key][7]
          
        
        # if there are more tha one possibility a variable, all of them will be saved 
        c = 0
        for key in dic_r:
            if dic_r[key][0] != final_TrunkDM:
                c = c+1
                if c == 1:
                    addition = addition+dic_r[key][1]+","+str(dic_r[key][4])
                else:
                    addition = addition+"   ;   "+dic_r[key][1]+","+str(dic_r[key][4])
                    
                    
    # If no option is chosen, these variables will acquire their initial value
    else:
        final_TrunkDM = calibrated_delta_MH
        final_TrunkSequence = seq
        final_TrunkLabel = ""
        final_mods_position = ""
        final_Trunk_Label_ppm = ""
        final_New_Theo_MH = Theo_mh 
        final_New_DM = calibrated_delta_MH
        
        
    
       
    return final_TrunkSequence,final_TrunkDM,final_TrunkLabel,final_mods_position,minimun,final_Trunk_Label_ppm,len(dic_r),addition,final_New_DM,final_New_Theo_MH





##################
# Main functions #
##################

def main(file,file1,infile1, infilefasta):
    """
    Reading configuration file
    """
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file) # Read MassMod configuration file

    
    logging.info("Reading MassMod configuration file")
        
    Mproton = config["Masses"].getfloat("m_proton")
    Hydrogen = config["Masses"].getfloat("m_hydrogen")
    O2 = config["Masses"].getfloat("m_oxygen")
        
    # dicctionary of aa at their masses is created
    dic_aa = {}       
    for option, value in config["Aminoacids"].items():
        dic_aa[option.strip("\n").upper()] = float(value.strip("\n").replace(",","."))
        
    # A  dicctionary whith fix mosidifcations indicating their "new name", and the mass is created  
    dic_mod = {}   
    for option, value in config["Fix_modifications"].items():
        option = option.upper()
        new_name = value.split("\t")[1][1:]
        value = float(value.split("\t")[0].strip("\n").replace(",","."))
        if option.upper() == "NT":
            NT_label = value
        
        else:        
            if value == NT_label:
                selectedaa = option.upper() # The aa that has a fix label (ej ; K-TMT) is saved
        
            value2 = value+dic_aa[option]
            dic_mod[option[0]] = new_name,value2,value
            
    
    
    
    
    config.read(file1) # Read Trunksolver configuration file
    
    logging.info("Reading TrunkSolver configuration file")
              
    Error = config["TrunkSolver_Parameters"].getfloat("Relative_Error") # Relative error is saved        
    Exp_mh_column_name = config["TrunkSolver_Parameters"].get("exp_mh_column_name") # Experimental mmh column name
    Theo_mh_column_name = config["TrunkSolver_Parameters"].get("theo_mh_column_name") # Theoretical mh name
    fix_mod_column_name = config["TrunkSolver_Parameters"].get("static_modifications_column_name") # Fix modifications column name  
    MasterProtein_column_name = config["TrunkSolver_Parameters"].get("MasterProtein_column_name") # Master protein column nam e
    Seq_column_name = config["TrunkSolver_Parameters"].get("Sequence_column_name") # Sequence colum name
    Delta_MH_cal_column_name = config["TrunkSolver_Parameters"].get("Calibrated_delta_mh_column_name") # DM column name
    New_Deltamass_column_name = config["TrunkSolver_Parameters"].get("New_Deltamass_column_name") # New deltamass  column name
    New_Theo_mh_column_name = config["TrunkSolver_Parameters"].get("New_Theo_mh_column_name") # New theoretical mh column name
    decnum =  config["TrunkSolver_Parameters"].get("decnum") # Number of decimales
    decnum = "."+str(decnum)+"f"




    # All labels that want to be checked are save in dic_CombList dicctionary
    dic_CombList = {}   
    for option, value in config["TrunkSolver_CombList"].items():
        dic_CombList[option.strip("\n").upper()] = float(value.strip("\n").replace(",","."))
        
    
    

    dic_fasta = SeqIO.index(infilefasta, "fasta")
    df = readInfile(infile1)
    df["TrunkSequence"]=""
    df["TrunkDM"]=np.nan
    df["TrunkLabel"]=""
    df["TrunkLabel_ppm"]=np.nan
    df[New_Theo_mh_column_name]=np.nan
    df[New_Deltamass_column_name]=np.nan  
    df["Static_modifications_position"]=""
    df["match_number"]=""
    df["Possible_option"]=""
    
    logging.info("Processing input file")
    cont=0
    for index, row in df.iterrows():
  
        if row[Seq_column_name].find("_") != -1:
            final_TrunkSequence = row[Seq_column_name]
            final_TrunkDM = row[Delta_MH_cal_column_name]
            final_TrunkLabel = " "
            final_mods_position = " "
            final_Trunk_Label_ppm = " "
            match_number = 0
            final_New_Theo_MH = row[Theo_mh_column_name]
            final_New_DM = row[Delta_MH_cal_column_name]
                

        else:


    
            dic_seqs,result=Obtain_values(row[Seq_column_name],row[MasterProtein_column_name],dic_fasta)

            final_TrunkSequence,final_TrunkDM,final_TrunkLabel,final_mods_position,minimun,final_Trunk_Label_ppm,match_number,addition,final_New_DM,final_New_Theo_MH = TrunkSolver(row[Seq_column_name],dic_seqs,row[Exp_mh_column_name],row[Delta_MH_cal_column_name],result,Error,dic_aa,dic_CombList,dic_mod,NT_label,selectedaa,decnum,Mproton,Hydrogen,O2,row[Theo_mh_column_name])
 
            
        
        df.loc[cont,"TrunkSequence"] = final_TrunkSequence
        df.loc[cont,"TrunkDM"] = final_TrunkDM
        df.loc[cont,"TrunkLabel"] = final_TrunkLabel
        df.loc[cont,"TrunkLabel_ppm"] = final_Trunk_Label_ppm
        df.loc[cont,New_Theo_mh_column_name]= final_New_Theo_MH
        df.loc[cont,New_Deltamass_column_name] = final_New_DM
        df.loc[cont,"Static_modifications_position"] = final_mods_position
        df.loc[cont,"match_number"]= match_number
        df.loc[cont,"Possible_option"]= addition

        cont=cont+1
    

    
    # write outputfile
    logging.info("Writing output file")

    outfilename = infile1[:-4]
    outfile= outfilename+"_TrunkSolved"+".txt"
    df.to_csv(outfile, index=False, sep='\t', encoding='utf-8')


    logging.info('end script')



if __name__ == '__main__':
    try:
        remove('Solver.ini')
    except:
        None


    # parse arguments
    parser = argparse.ArgumentParser(
        description='TrunkSolver',
        epilog='''
        Example:
            python TrunkSolver.py
        ''')
      
    # default TrunkSolver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Path to input file')
    parser.add_argument('-f', '--fastafile', required=True, help='Path to input fastafile')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    
    # these will overwrite the config if specified
    parser.add_argument('-r', '--relerror', help='Maximum allowable relative error (ppm)')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
   
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.relerror is not None:
        config.set('TrunkSolver_Parameters', 'Relative_Error', str(args.relerror))
        config.set('Logging', 'create_ini', '1')
   
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/Solver.ini', 'w') as newconfig:
            config.write(newconfig)
        
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'TrunkSolver_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'TrunkSolver_log_debug.txt'
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
    infilefasta=args.fastafile
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    

    # configuration files are read      
    try:
        open('Solver.ini',"r")
        trunkini='Solver.ini'
        logging.info("Modified Trunkolver configuration file is going to be use")
        
    except:
        open("config/Solver.ini","r")
        trunkini="config/Solver.ini"

    main("config/MassMod.ini",trunkini, infile1,infilefasta)
