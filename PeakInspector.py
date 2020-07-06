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
import matplotlib.pyplot as plt
import matplotlib.ticker
from bokeh.plotting import figure, output_file, show
from bokeh.models import SingleIntervalTicker, LinearAxis
from bokeh.layouts import gridplot
import tkinter as tk


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
        

'''
def readInfile(infile):
    
    Read table in infile given by the user using pandas
    

    try:
        df = pd.read_csv(infile, sep="\t", float_precision='high')
    except:
        log_str = 'error: ' + str(infile) + ' could not be openned.'
        logging.info(log_str)
        sys.exit(log_str)

    return df
'''

def plot_bottom_graph(main_plot, letter, letter_to_colInfo, df):
    '''
    Represent bottom graph using main_plot as reference. It will use its
    x axis, so they are coupled. The function will return the bottom plot figure
    object.
    '''

    # Extract information from bottom plot
    bottom_plot_name = letter_to_colInfo[letter]['Name']
    bottom_plot_index = letter_to_colInfo[letter]['Index']
    bottom_plot_color = letter_to_colInfo[letter]['Color']

    #Build bottom plot
    bottom_plot = figure(title=bottom_plot_name + " representation",\
         x_axis_label='Delta mass', y_axis_label=bottom_plot_name,\
         width=1300, height=400, x_range=main_plot.x_range, tools = "pan,yzoom_in,yzoom_out,wheel_zoom,box_zoom,reset,save,undo")#, y_range=main_plot.y_range)

    try:
        bottom_plot.xaxis.ticker.desired_num_ticks = 30
		
    except AttributeError:
        logging.info(f"AttributeError: {err}")
        print("bokeh package needs to be updated (pip install bokeh -U)")
        sys.exit()
		
    bottom_plot.line(df.df.iloc[:, 1], df.df.iloc[:, bottom_plot_index], line_width=2, color=bottom_plot_color)

    return bottom_plot

def plot_graphs(plot_letters, letter_to_colInfo, df):
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
    first_plot_index = letter_to_colInfo[plot_letters[0]]['Index']
    first_plot_color = letter_to_colInfo[plot_letters[0]]['Color']

    # Build the first plot
    p1 = figure(title=first_plot_name + " representation",\
         x_axis_label='Delta mass', y_axis_label=first_plot_name,\
         width=1300, height=400, tools = "pan,yzoom_in,yzoom_out,wheel_zoom,box_zoom,reset,save,undo")

    try:
        p1.xaxis.ticker.desired_num_ticks = 30
	
    except AttributeError as err:
        logging.info(f"AttributeError: {err}")
        print("bokeh package needs to be updated (pip install bokeh -U)")
        sys.exit()

    p1.line(df.df.iloc[:, 1], df.df.iloc[:, first_plot_index], line_width=2, color=first_plot_color)

    # If there are more plots, these are represented below using plot_bottom_graph function
    if len(plot_letters) > 1:

        bottom_graphs_list = [[plot_bottom_graph(p1, letter, letter_to_colInfo, df)] for letter in plot_letters[1:]]

        all_graphs_list = [[p1]] + bottom_graphs_list

        # Plot threshold
        logging.info("Plotting thresholds")
        all_graphs_list = [plot_threshold(pi[0], letter_to_colInfo[letter]['Threshold'], df) \
            for letter, pi in zip(plot_letters, all_graphs_list)]

        # Show the plot
        plot = gridplot(all_graphs_list)
        show(plot)

    else:
        logging.info("Plotting threshold")
        plot_threshold(p1, letter_to_colInfo[plot_letters[0]]['Threshold'], df)

        show(p1)
    
    
    logging.info('Graphs plotted')

    return 0


def parse_entry_threshold(entry_list):
    '''
    Function used to parse threshold values introduced by the user
    '''

    return [float(entry.get()) for entry in entry_list if entry.get().lower() not in ['', 'none']]


def plot_threshold(pi, threshold, df):
    '''
    Function used to plot threshold lines
    Input:  pi: Figure in which threshold is plotted
            threshold: List with threshold values
            df: Pandas data frame with all values
    Return: It returns the figure with added threshold
    '''

    for threshold_i in threshold:
        threshold_i_Y_values = np.ones_like(df.df.iloc[:, 1])*threshold_i
        pi.line(df.df.iloc[:, 1], threshold_i_Y_values, line_color='black', line_dash="4 4")
    
    return [pi]


def click_plot_buttom(letter_to_colInfo, df, freq, smooth_freq, s1, s2, freq_entry, smooth_freq_entry,\
    s1_entry1, s1_entry2, s2_entry1, s2_entry2):
    '''
    Function executed hen the user press the buttom Plot. It receives all parameters required to plot
    '''

    all_entries = [[freq_entry], [smooth_freq_entry], [s1_entry1, s1_entry2], [s2_entry1, s2_entry2]]

    plot_letters = [letter for letter, bool_plot in zip(['A', 'B', 'C', 'D'], \
        [freq.get(), smooth_freq.get(), s1.get(), s2.get()]) if bool_plot]

    [letter_to_colInfo[letter].update({'Threshold': parse_entry_threshold(entry_list)})\
         for letter, entry_list in zip(['A', 'B', 'C', 'D'], all_entries)]
    
    plot_graphs(plot_letters, letter_to_colInfo, df)



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
    letter_to_colInfo = {'A': {'Index': 2, 'Name': 'Frequency', 'Color': 'navy', 'Threshold': []},\
                         'B': {'Index': 3, 'Name': 'Smooth Frequency', 'Color': 'steelblue', 'Threshold': []},\
                         'C': {'Index': 4, 'Name': 'Slope 1', 'Color': 'firebrick', 'Threshold': []},\
                         'D': {'Index': 5, 'Name': 'Slope 2', 'Color': 'salmon', 'Threshold': []}}
    
    
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
            lbl2_slope2_entry2))
    plot_btn.place(x=5, y=290)


    # If parameters -f, -sf, -s1 or -s2 were used in the execution, the plot is shown. It is as if plot button
    # was pressed
    if any([args.f, args.sf, args.s1, args.s2]) and args.infile:

        click_plot_buttom(letter_to_colInfo, df_object,\
                frequency_check_state, smooth_frequency_check_state, slope1_check_state, slope2_check_state, \
                lbl2_freq_entry, lbl2_smooth_freq_entry, lbl2_slope1_entry1, lbl2_slope1_entry2, lbl2_slope2_entry1,\
                lbl2_slope2_entry2)

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

    parser.add_argument('-i', '--infile', help='Path to input file', type=str)
    parser.add_argument('-c', '--config', help='Path to custom config.ini file', type=str)

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

    args = parser.parse_args()

    # If user use config file, parse it
    if args.config:
        config  = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(args.config)

        # Assign to parse variables which data is plotted
        args.f, args.sf, args.s1, args.s2 = config['Plots']['frequency'], config['Plots']['smooth_frequency'],\
            config['Plots']['slope1'], config['Plots']['slope2']

        # Assign to parse variables the thresholds defined
        args.fT, args.sfT, args.s1T1, args.s1T2, args.s2T1, args.s2T2 = config['Thresholds']['frequency_T'],\
            config['Thresholds']['smooth_frequency_T'], config['Thresholds']['slope1_T1'], \
            config['Thresholds']['slope1_T2'], config['Thresholds']['slope2_T1'], config['Thresholds']['slope2_T2']


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