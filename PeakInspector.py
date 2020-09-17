#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Rafael Barrero Rodriguez"
__credits__ = ["Rafael Barrero Rodriguez", "Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "rbarreror@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# Import modules
import os
import sys
from pathlib import Path
import argparse
import configparser
import logging
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.ticker
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import SingleIntervalTicker, LinearAxis
from bokeh.layouts import gridplot
import tkinter as tk

import pdb


###################
# Local functions #
###################

class Histogram_Class:
    '''
    DataFrame_object with data from histogram:
        self.df: Pandas dataframe
        self.path: Directory path from which dataframe was taken
    '''

    def __init__(self):
        self.path = None
        self.df = None


    def ReadFromPath(self, path):
        '''
        Method used to read histogram from path and storing it ass a pandas dataframe
        in self.df and the source in self.path
        '''

        log_str = 'Reading file: ' + str(Path(path))
        logging.info(log_str)

        try:
            self.df = pd.read_csv(path, sep="\t", float_precision='high')
            self.path = path
        
        except:
            log_str = 'error: ' + str(path) + ' could not be openned.'
            logging.info(log_str)
            sys.exit(log_str)
        

    def ReadFromGUI(self):
        '''
        Method used to get filename from a dialog and read the histogram selected
        '''

        logging.info("Getting filename from GUI")

        filename = tk.filedialog.askopenfilename()

        self.ReadFromPath(filename)
        

def plot_bottom_graph(main_plot, letter, letter_to_colInfo, df):
    '''
    Represent bottom graph using main_plot as reference. It will use its
    x axis, so they are coupled. The function will return the bottom plot figure
    object.
    '''

    # Extract information from bottom plot
    bottom_plot_name = letter_to_colInfo[letter]['Name']
    bottom_plot_color = letter_to_colInfo[letter]['Color']
    # bottom_plot_index = letter_to_colInfo[letter]['Index']

    bottom_plot_column_name = letter_to_colInfo[letter]['ColumnName']

    try:
        bottom_plot_index = np.where(df.df.columns == bottom_plot_column_name)[0][0]

    except IndexError:
        logging.info(f"Error: {bottom_plot_column_name} column was not found. Is it correctly written?")
        return ""

    #Build bottom plot
    bottom_plot = figure(title=bottom_plot_name + " representation",\
         x_axis_label='Delta mass', y_axis_label=bottom_plot_name,\
         width=1300, height=400, x_range=main_plot.x_range,\
         tools = "pan,xzoom_in,xzoom_out,ywheel_zoom,box_zoom,reset,save,undo,hover", tooltips=[("Name", "$name")])

    try:
        bottom_plot.xaxis.ticker.desired_num_ticks = 30
		
    except AttributeError:
        logging.info(f"AttributeError: {err}")
        print("bokeh package needs to be updated (pip install bokeh -U)")
        sys.exit()
		
    bottom_plot.line(df.df.iloc[:, 1], df.df.iloc[:, bottom_plot_index], line_width=2, color=bottom_plot_color, name=bottom_plot_name)

    return bottom_plot


def plot_threshold(pi, threshold, df):
    '''
    Function used to plot threshold lines
    Input:  pi: Figure in which threshold is plotted
            threshold: List with threshold values
            df: Pandas data frame with all values
    Return: It returns the figure with added threshold
    '''

    min_mz = np.min(df.df.iloc[:, 1])
    max_mz = np.max(df.df.iloc[:, 1])

    x_axis = (min_mz, max_mz)

    for threshold_i in threshold:
        # threshold_i_Y_values = np.ones_like(df.df.iloc[:, 1])*threshold_i
        y_axis = np.ones_like(x_axis)*threshold_i
        pi.line(x_axis, y_axis, line_color='black', line_dash="4 4")
    
    return [pi]


def plot_pleak(pi, peaks_list, column_name, df):
    '''
    Input:
        - pi: Bokeh plot figure in which peaks are plotted
        - peaks_list: List of floats containing theoretical DM given by the user
        - column_name: String with name of the column plotted
        - df: Dataframe object containing the histogram
    Return:
        - [pi]: List containing the modified plot
    '''

    min_value = np.min(df.df.loc[:, column_name])
    min_value = min_value - abs(0.1*min_value)
    max_value = np.max(df.df.loc[:, column_name])*1.1

    # y_axis = np.arange(min_value, max_value)
    y_axis = (min_value, 0, max_value)
    for peak, peak_name in peaks_list:
        x_axis = np.ones_like(y_axis)*peak
        pi.line(x_axis, y_axis, line_color='green', line_width=2, name=peak_name)
    
    return [pi]


def iniMaker(filename, plot_letters, letter_to_colInfo):
    '''

    '''

    config.set('Parameters', 'Infile', str(args.infile))
    config.set('Parameters', 'Peaks', str(args.tDM))

    config.set('Plots', 'frequency', str('A' in plot_letters))
    config.set('Plots', 'smooth_frequency', str('B' in plot_letters))
    config.set('Plots', 'slope1', str('C' in plot_letters))
    config.set('Plots', 'slope2', str('D' in plot_letters))

    if letter_to_colInfo['A']['Threshold']:
        config.set('Thresholds', 'frequency_T', str(letter_to_colInfo['A']['Threshold'][0]))
    else:
        config.set('Thresholds', 'frequency_T', '0')

    if letter_to_colInfo['B']['Threshold']:
        config.set('Thresholds', 'smooth_frequency_T', str(letter_to_colInfo['B']['Threshold'][0]))
    else:
        config.set('Thresholds', 'smooth_frequency_T', '0')

    if letter_to_colInfo['C']['Threshold']:
        config.set('Thresholds', 'slope1_T1', str(letter_to_colInfo['C']['Threshold'][0]))
    else:
        config.set('Thresholds', 'slope1_T1', '0')

    if len(letter_to_colInfo['C']['Threshold']) == 2:
        config.set('Thresholds', 'slope1_T2', str(letter_to_colInfo['C']['Threshold'][1]))
    else:
        config.set('Thresholds', 'slope1_T2', '0')

    if letter_to_colInfo['D']['Threshold']:
        config.set('Thresholds', 'slope2_T1', str(letter_to_colInfo['D']['Threshold'][0]))
    else:
        config.set('Thresholds', 'slope2_T1', '0')

    if len(letter_to_colInfo['D']['Threshold']) == 2:
        config.set('Thresholds', 'slope2_T2', str(letter_to_colInfo['D']['Threshold'][1]))
    else:
        config.set('Thresholds', 'slope2_T2', '0')

    ini_path = os.path.join(os.path.dirname(args.infile), filename[:-5] + '.ini')

    with open(ini_path, 'w') as newconfig:
        logging.info(f"Saving .ini: {ini_path}")
        config.write(newconfig)



def savePlot(plot, plot_letters, letter_to_colInfo):
    '''
    Input:
        - plot: Bokeh plot to be saved
        - plot_letters: List of strings containing the letters associated to the plotted graphs
        - letter_to_colInfo: Dictionary that associates each letter to the plot information
    '''

    # Build filename: 
    filename = '_'.join([letter_to_colInfo[letter]['ColumnName'] + '_' + '_'.join([str(int(thr)) for thr in letter_to_colInfo[letter]['Threshold']]) \
        for letter in plot_letters])
    
    filename += '.html'
    filename = filename.replace('__', '_')
    filename = filename.replace('-', 'm')

    filename = os.path.join(os.path.dirname(args.infile), filename)

    save(plot, filename=filename)
    logging.info(f"Graph saved: {filename}")

    # Create ini associated to this plot
    if config.getint('Logging', 'create_ini'):
        iniMaker(filename, plot_letters, letter_to_colInfo)


def plot_graphs(plot_letters, letter_to_colInfo, df, peaks_list):
    '''
    Function used to make all plots. The first plot will receive a different
    treatment than the others
    Input:  plot_letters: list of letters selected by the user
            letter_to_colInfo: Dictionary with the information of each plot
            df: Dataframe with all data
    '''
    
    log_str = 'Plotting the following graphs: ' + ", ".join([letter_to_colInfo[letter]['Name'] for letter in plot_letters])
    logging.info(log_str)

    output_file(df.path + '_plot.html')

    # Extract information from the first plot
    first_plot_name = letter_to_colInfo[plot_letters[0]]['Name']
    first_plot_color = letter_to_colInfo[plot_letters[0]]['Color']
    #first_plot_index = letter_to_colInfo[plot_letters[0]]['Index']
    
    first_plot_column_name = letter_to_colInfo[plot_letters[0]]['ColumnName']

    try:
        first_plot_index = np.where(df.df.columns == first_plot_column_name)[0][0]

    except IndexError:
        logging.info(f"Error: {first_plot_column_name} column was not found. Is it correctly written?")
        return 0

    
    # Build the first plot
    p1 = figure(title=first_plot_name + " representation",\
         x_axis_label='Delta mass', y_axis_label=first_plot_name,\
         width=1300, height=400, tools = "pan,xzoom_in,xzoom_out,ywheel_zoom,box_zoom,reset,save,undo,hover", tooltips=[('Name', '$name')])

    try:
        p1.xaxis.ticker.desired_num_ticks = 30
	
    except AttributeError as err:
        logging.info(f"AttributeError: {err}")
        print("bokeh package needs to be updated (pip install bokeh -U)")
        sys.exit()

    p1.line(df.df.iloc[:, 1], df.df.iloc[:, first_plot_index], line_width=2, color=first_plot_color, name=first_plot_name)

    # If there are more plots, these are represented below using plot_bottom_graph function
    if len(plot_letters) > 1:

        bottom_graphs_list = [[plot_bottom_graph(p1, letter, letter_to_colInfo, df)] for letter in plot_letters[1:]]
        bottom_graphs_list = [plot for plot in bottom_graphs_list if plot != [""]]

        all_graphs_list = [[p1]] + bottom_graphs_list

        # Plot threshold
        logging.info("Plotting thresholds")
        all_graphs_list = [plot_threshold(pi[0], letter_to_colInfo[letter]['Threshold'], df) \
            for letter, pi in zip(plot_letters, all_graphs_list)]
        
        logging.info("Plotting theoretical DM")
        all_graphs_list = [plot_pleak(pi[0], peaks_list, letter_to_colInfo[letter]['ColumnName'], df) \
            for letter, pi in zip(plot_letters, all_graphs_list)]

        # Show the plot
        plot = gridplot(all_graphs_list)
        show(plot)
        savePlot(plot, plot_letters, letter_to_colInfo)

    else:
        logging.info("Plotting threshold")
        plot_threshold(p1, letter_to_colInfo[plot_letters[0]]['Threshold'], df)

        logging.info("Plotting theoretical DM")
        plot_pleak(p1, peaks_list, letter_to_colInfo[plot_letters[0]]['ColumnName'], df)

        show(p1)
        savePlot(p1, plot_letters, letter_to_colInfo)
    
    
    logging.info('Graphs plotted')

    return 0


def parse_entry_threshold(entry_list):
    '''
    Function used to parse threshold values introduced by the user
    '''

    return [float(entry.get()) for entry in entry_list if entry.get().lower() not in ['', 'none']]


def getPeaks(peaks_str):
    '''
    Input:
        - peaks_str: String containing the theoretical DM given by the user
    Return:
        - peaks_list: List of pairs. The first element of each pair is a float with the theoretical DM and the
        second element is its associated name
    '''

    peaks_str += '\n' if peaks_str[-1] != '\n' else peaks_str

    match = re.search(r'(-?\d+(?:\.\d*)?)(?:\n|(?:\s*,\s*|\s+)([^\n]*)\n)', peaks_str)
    peaks_list = []
    
    while match:
        # List of double containing m/z and name. If there is no name it stores "NA"
        peak_tmp = [float(match.groups()[0])]
        peak_tmp.append("NA") if not match.groups()[1] else peak_tmp.append(match.groups()[1])
        # Add pair to peak_list object
        peaks_list.append(peak_tmp) 

        # recompose peaks_str and apply re.search
        peaks_str = peaks_str[match.span()[1]:]
        match = re.search(r'(-?\d+(?:\.\d*)?)(?:\n|(?:\s*,\s*|\s+)([^\n]*)\n)', peaks_str)

    return peaks_list



def click_plot_buttom(letter_to_colInfo, df, freq, smooth_freq, s1, s2, freq_entry, smooth_freq_entry,\
    s1_entry1, s1_entry2, s2_entry1, s2_entry2, peaks):
    '''
    Function executed when the user press the buttom Plot. It receives all parameters required to plot
    '''

    all_entries = [[freq_entry], [smooth_freq_entry], [s1_entry1, s1_entry2], [s2_entry1, s2_entry2]]

    plot_letters = [letter for letter, bool_plot in zip(['A', 'B', 'C', 'D'], \
        [freq.get(), smooth_freq.get(), s1.get(), s2.get()]) if bool_plot]

    [letter_to_colInfo[letter].update({'Threshold': parse_entry_threshold(entry_list)})\
         for letter, entry_list in zip(['A', 'B', 'C', 'D'], all_entries)]

    # Get theoretical DM given by the user
    peaks_list = getPeaks(peaks)
    
    plot_graphs(plot_letters, letter_to_colInfo, df, peaks_list)



##################
# Main functions #
##################

def main(args):
    '''
    Main function
    '''

    # Read infile
    df_object = Histogram_Class()

    if args.infile:
        df_object.ReadFromPath(args.infile)
    
    # Create a dictionary used to store parameters associated to each plot
    letter_to_colInfo = {'A': {'ColumnName': 'count', 'Index': 2, 'Name': 'Frequency', 'Color': 'navy', 'Threshold': []},\
                         'B': {'ColumnName': 'smooth_count', 'Index': 3, 'Name': 'Smooth Frequency', 'Color': 'steelblue', 'Threshold': []},\
                         'C': {'ColumnName': 'slope1', 'Index': 4, 'Name': 'Slope 1', 'Color': 'firebrick', 'Threshold': []},\
                         'D': {'ColumnName': 'slope2', 'Index': 5, 'Name': 'Slope 2', 'Color': 'salmon', 'Threshold': []}}
    
    
    # Create GUI using tkinter
    window = tk.Tk()
    window.title("PeakInspector")
    window.geometry('400x350')

    # Create Menu to import file
    menu = tk.Menu(window)
    new_item = tk.Menu(menu)
    new_item.add_command(label='Open', command=df_object.ReadFromGUI)
    menu.add_cascade(label='File', menu=new_item)

    window.config(menu=menu)

    # Create first section, that enables plot selection
    lbl = tk.Label(window, text="Select Data to Plot", font=("comicsans", 10))
    lbl.place(x=5, y=5)

    frequency_check_state = tk.BooleanVar()
    frequency_check_state.set(args.f)
    frequency_check = tk.Checkbutton(window, text=' Frequency', var=frequency_check_state, font=("comicsans", 10))
    frequency_check.place(x=25, y=30)

    smooth_frequency_check_state = tk.BooleanVar()
    smooth_frequency_check_state.set(args.sf)
    smooth_frequency_check = tk.Checkbutton(window, text=' Smooth Frequency', var=smooth_frequency_check_state, font=("comicsans", 10))
    smooth_frequency_check.place(x=25, y=60)

    slope1_check_state = tk.BooleanVar()
    slope1_check_state.set(args.s1)
    slope1_check = tk.Checkbutton(window, text=' Slope 1', var=slope1_check_state, font=("comicsans", 10))
    slope1_check.place(x=25, y=90)

    slope2_check_state = tk.BooleanVar()
    slope2_check_state.set(args.s2)
    slope2_check = tk.Checkbutton(window, text=' Slope 2', var=slope2_check_state, font=("comicsans", 10))
    slope2_check.place(x=25, y=120)

    # Create second section with entry spaces to introduce thresholds
    lbl2 = tk.Label(window, text="Define Thresholds (e.g. 5e3)", font=("comicsans", 10))
    lbl2.place(x=5, y=180)

    # Frequency entry
    lbl2_freq = tk.Label(window, text="Frequency", font=("comicsans", 10))
    lbl2_freq.place(x=5, y=210)

    lbl2_freq_entry = tk.Entry(window, width=10)
    lbl2_freq_entry.place(x=5, y=230)

    lbl2_freq_entry.insert(tk.END, args.fT)

    # Smooth Frequency entry
    lbl2_smooth_freq = tk.Label(window, text="Smooth Frequency", font=("comicsans", 10))
    lbl2_smooth_freq.place(x=95, y=210)

    lbl2_smooth_freq_entry = tk.Entry(window, width=10)
    lbl2_smooth_freq_entry.place(x=95, y=230)

    lbl2_smooth_freq_entry.insert(tk.END, args.sfT)

    # Slope 1 entry
    lbl2_slope1 = tk.Label(window, text="Slope 1", font=("comicsans", 10))
    lbl2_slope1.place(x=225, y=210)

    lbl2_slope1_entry1 = tk.Entry(window, width=10)
    lbl2_slope1_entry1.place(x=225, y=230)

    lbl2_slope1_entry1.insert(tk.END, args.s1T1)

    lbl2_slope1_entry2 = tk.Entry(window, width=10)
    lbl2_slope1_entry2.place(x=225, y=250)

    lbl2_slope1_entry2.insert(tk.END, args.s1T2)


    # Slope 2 entry
    lbl2_slope2 = tk.Label(window, text="Slope 2", font=("comicsans", 10))
    lbl2_slope2.place(x=315, y=210)

    lbl2_slope2_entry1 = tk.Entry(window, width=10)
    lbl2_slope2_entry1.place(x=315, y=230)

    lbl2_slope2_entry1.insert(tk.END, args.s2T1)


    lbl2_slope2_entry2 = tk.Entry(window, width=10)
    lbl2_slope2_entry2.place(x=315, y=250)

    lbl2_slope2_entry2.insert(tk.END, args.s2T2)


    # Create Plot buttom
    plot_btn = tk.Button(window, text="Plot", command=lambda: click_plot_buttom(letter_to_colInfo, df_object,\
            frequency_check_state, smooth_frequency_check_state, slope1_check_state, slope2_check_state, \
            lbl2_freq_entry, lbl2_smooth_freq_entry, lbl2_slope1_entry1, lbl2_slope1_entry2, lbl2_slope2_entry1,\
            lbl2_slope2_entry2, args.tDM))
    plot_btn.place(x=5, y=290)


    # If parameters -f, -sf, -s1 or -s2 were used in the execution, the plot is shown. It is as if plot button
    # was pressed
    if any([args.f, args.sf, args.s1, args.s2]) and args.infile:

        click_plot_buttom(letter_to_colInfo, df_object,\
                frequency_check_state, smooth_frequency_check_state, slope1_check_state, slope2_check_state, \
                lbl2_freq_entry, lbl2_smooth_freq_entry, lbl2_slope1_entry1, lbl2_slope1_entry2, lbl2_slope2_entry1,\
                lbl2_slope2_entry2, args.tDM)

    window.mainloop()

    return 0
    
    

if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PeakInspector',
        epilog='''
        Example:
            python PeakInspector.py
        
        '''
    )

    defaultconfig = os.path.join(os.path.dirname(__file__), "config/PeakInspector.ini")

    parser.add_argument('-i', '--infile', help='Path to input file', type=str)
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file', type=str)

    parser.add_argument('-v', dest='verbose', action='store_true', help='Increase output verbosity')

    # Arguments to set which data is plotted
    parser.add_argument('-f', action='store_true', help='Plot frequency data')
    parser.add_argument('-sf', action='store_true', help='Plot smooth frequency data')
    parser.add_argument('-s1', action='store_true', help='Plot slope 1 data')
    parser.add_argument('-s2', action='store_true', help='Plot slope 2 data')

    # Arguments to set thresholds
    parser.add_argument('-fT', default='', help='Define threshold for frequency plot', type=str)
    parser.add_argument('-sfT', default='', help='Define threshold for smooth frequency plot', type=str)
    parser.add_argument('-s1T1', default='', help='Define threshold 1 for slope 1 plot', type=str)
    parser.add_argument('-s1T2', default='', help='Define threshold 2 for slope 1 plot', type=str)
    parser.add_argument('-s2T1', default='', help='Define threshold 1 for slope 2 plot', type=str)
    parser.add_argument('-s2T2', default='', help='Define threshold 2 for slope 2 plot', type=str)

    # Arguments to set theoretical DM
    parser.add_argument('-tDM', default='', help='Set theoretical DM peaks to be plotted (e.g. "0.98, 100")', type=str)

    args = parser.parse_args()

    # If user use config file, parse it
    if args.config:
        config  = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(args.config)
        
        
    if not args.f and config['Plots']['frequency'].lower().strip() == 'true':
        args.f = config['Plots']['frequency']
        config.set('Logging', 'create_ini', '1')
        
    if not args.sf and config['Plots']['smooth_frequency'].lower().strip() == 'true':
        args.sf = config['Plots']['smooth_frequency']
        config.set('Logging', 'create_ini', '1')

    if not args.s1 and config['Plots']['slope1'].lower().strip() == 'true':
        args.s1 = config['Plots']['slope1']
        config.set('Logging', 'create_ini', '1')

    if not args.s2 and config['Plots']['slope2'].lower().strip() == 'true':
        args.s2 = config['Plots']['slope2']
        config.set('Logging', 'create_ini', '1')

    if not args.fT:
        args.fT = config['Thresholds']['frequency_T']
        config.set('Logging', 'create_ini', '1')

    if not args.sfT:
        args.sfT = config['Thresholds']['smooth_frequency_T']
        config.set('Logging', 'create_ini', '1')

    if not args.s1T1:
        args.s1T1 = config['Thresholds']['slope1_T1']
        config.set('Logging', 'create_ini', '1')

    if not args.s1T2:
        args.s1T2 = config['Thresholds']['slope1_T2']
        config.set('Logging', 'create_ini', '1')

    if not args.s2T1:
        args.s2T1 = config['Thresholds']['slope2_T1']
        config.set('Logging', 'create_ini', '1')

    if not args.s2T2:
        args.s2T2 = config['Thresholds']['slope2_T2']
        config.set('Logging', 'create_ini', '1')

    if not args.infile:
        args.infile = config['Parameters']['Infile']
        config.set('Logging', 'create_ini', '1')

    if not args.tDM:
        args.tDM = config['Parameters']['Peaks']
        config.set('Logging', 'create_ini', '1')
        

    # logging debug level. By default, info level
    if args.infile:
        log_file = outfile = args.infile[:-4] + '_log.txt'
        log_file_debug = outfile = args.infile[:-4] + '_log_debug.txt'
    
    else:
        log_file = outfile = 'log.txt'
        log_file_debug = outfile = 'log_debug.txt'


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