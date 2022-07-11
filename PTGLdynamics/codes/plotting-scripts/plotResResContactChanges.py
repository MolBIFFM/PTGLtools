########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "1.0.0"  

########### built-in imports ###########

import sys
import os
import argparse
import logging
import traceback
import pathlib
import re
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import platform
import subprocess
import numpy as np


########### functions ###########

def log(message, level=""):
    """Prints the message according to level of severity and output settings"""
    if (level == "c"):
        logging.critical(message)
    if (level == "e"):
        logging.error(message)
    elif (level == "w"):
        logging.warning(message)
    elif (level == "i"):
        logging.info(message)
    elif (level == "d"):
        logging.debug(message)
    else:
        if (args.outputfile == ""):
            print(message)
        else:
            global output_file
            output_file.write(message + "\n")
            

def check_dir_args(argument):
    """Checks directory arguments and returns its value if the directory exists or exits the program otherwise."""
    
    if( argument != ""):
        if(os.path.isdir(argument)):
            return argument
        else:
            logging.error("Specified directory '%s' does not exist. Exiting now.", argument)
            sys.exit(1)
    else:
        return os.getcwd()
        
def check_input_files(inputfile):
    """Returns the argument if the file is readable. Otherwise raises an error and exits."""
    if(inputfile != ""):
        if(os.access(inputfile, os.R_OK)):
            return inputfile
        else:
            logging.error("Specified input file '%s' is not readable. Exiting now.", inputfile)
            sys.exit(1)
    else:
        return ''

def read_in_file(file):
    """Returns a nested dictionary with each residue´s change and residues ordered according to their chain id. """
    changes = {}
    with open(file, 'r') as f:
        csv_lines = f.read().split("\n")
    header = csv_lines[0].split(',')
    header = header[0]
    for line in csv_lines[1:]:
        if line != '':
            columns = line.split(',')
            keys = columns[0].split('|')
            chain_id = keys[1]

            if chain_id in changes:
                changes[chain_id] [columns[0]] = columns[1]
                
            else:
                changes[chain_id] = {columns[0] : columns[1]}

    return changes, header

def split(elem):
    """Returns the residues ID from the 3-tuple."""
    key = elem[0]
    keys = key.split('|')
    residue = keys[0]
    residue = residue[2:-1]
    residue = int(residue)
    return residue
        
        
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
            
########### command line parser ###########
"""Arguments with hyphen like '--input-dir' are called with an underscore within the programm: args.input_dir
However, in the command line call the hyphen is used: --input-dir <path> """

## create the parser
cl_parser = argparse.ArgumentParser(description="Plots the changes in residue residue contacts. One or two csv files computed by calculateChanges can be given, for each chain the absolute change for each residue is plotted.",
                                    fromfile_prefix_chars="@")

## add arguments

# add mutually exclusive group for silent / verbose
loudness = cl_parser.add_mutually_exclusive_group()
loudness.add_argument('-s',
                       '--silent',
                       action='store_true',
                       help='only print results to stdout, e.g., no status messages and meta information')

loudness.add_argument('-v',
                       '--verbose',
                       action='store_true',
                       help='print more to stdout, e.g., status messages and meta information')

loudness.add_argument('-d',
                      '--debug',
                      action='store_true',
                      help='print everything including debug information')

# add command line arguments
cl_parser.add_argument('--version',
                       action='version',
                       version='%(prog)s ' + version)

cl_parser.add_argument('file_one',
                       metavar = 'file-one',
                       help = 'Enter the path to a csv file computed by calculateChanges containing the absolute changes.')

cl_parser.add_argument('-file-two',
                       '--file-two',
                       default = '',
                       metavar = 'file-two',
                       help='Enter the path to a second csv file computed by calculateChanges containing the absolute changes.')                                            

cl_parser.add_argument('-p',
                       '--output-dir',
                       metavar = 'output-directory',
                       default = '',
                       help = 'specify a path to your output directory. Otherwise the current folder is used.')
                       
cl_parser.add_argument('-n',
                       '--name-pdf-plots',
                       metavar = 'name-pdf-plots',
                       default = '',
                       help = 'specify the name of the output pdf where all plots are stored')    
                       
cl_parser.add_argument('--show-plots',
                       action='store_true',
                       help = 'open the matplot plots.')                                           

args = cl_parser.parse_args()


########### check arguments ###########

# assign log level
log_level = logging.WARNING
if (args.debug):
    log_level = logging.DEBUG
elif (args.verbose):
    log_level = logging.INFO
logging.getLogger().setLevel(log_level)

# inputfiles
file_one = check_input_files(args.file_one)
file_two = check_input_files(args.file_two)    

# output directory
output_dir = check_dir_args(args.output_dir)

# name of pdf output
if(args.name_pdf_plots != ""):
    if(args.name_pdf_plots.endswith('.pdf')):
        name_pdf_plots = args.name_pdf_plots
    else:
        name_pdf_plots = args.name_pdf_plots + '.pdf'
else:
    name_pdf_plots = "plots_res_res_contact_changes.pdf"

########### vamos ###########
no_U_dict = {"1":"Nqo1","2":"Nqo2","3":"Nqo3","4":"Nqo4","5":"Nqo5","6":"Nqo6","9":"Nqo9","7":"Nqo15","J":"Nqo10","K":"Nqo11","L":"Nqo12","M":"Nqo13","N":"Nqo14","H":"Nqo8","W":"Nqo16","A":"Nqo7"}
_start_time = time.time()

figure_count = 0
pp = PdfPages(name_pdf_plots)

changes_1, header_1 = read_in_file(file_one)
filename_1 = os.path.basename(file_one).split('.')

changes_2 = ''
header_2 = ''
if file_two != '':
    changes_2, header_2 = read_in_file(file_two)
    filename_2 = os.path.basename(file_two).split('.')
    if header_2 != header_1:
        log("The files aren´t in the same format. Exiting.", "e")
        exit()

for key in changes_1:
    changes_1[key] = dict(sorted(changes_1[key].items(), key=split))
           
for key in changes_1:
    if len(changes_1.get(key)) < 2:
        pass
    else:
        x = []
        y = []
        z = []

        for residue in changes_1[key]:
            residue_parts = residue.split('|')
            residue_id = residue_parts[0].replace('(', '')
            residue_id = residue_id.replace("'", "")
            x.append(residue_id)
            y.append(int(changes_1[key].get(residue)))
            if changes_2 != '':
                if residue in changes_2[key]:
                    value = int(changes_2[key].get(residue))
                    del changes_2[key][residue]
                else:
                    value = 0
                z.append(value)
        if changes_2 != '':
            if changes_2[key] != '':
                 for residue in changes_2[key]:
                    x.append(residue)
                    y.append(0)
                    z.append(int(changes_2[key].get(residue)))
        key_2 = key.replace("'", "")
        key_nq = no_U_dict[key_2]            
        counter_parts = 0
        parts = round(len(x) / 100)
        x = np.array_split(x, parts)
        y = np.array_split(y, parts)
        z = np.array_split(z, parts)
        while counter_parts < parts:                      
            plt.figure(figure_count) 
            figure_count +=1
            plt.gca().xaxis.grid(True)
            plt.xlabel("residues in " + str(key_nq) + " part " + str(counter_parts + 1) + " of " + str(parts), fontsize=4)
            plt.ylabel("sum of changes over all time steps")                   
            plt.xticks(rotation=90, fontsize=3)
            plt.plot(x[counter_parts],y[counter_parts], '-o', color="royalblue", alpha=0.65, label="NoQ", linewidth=1) 
            if file_two != '':
                plt.plot(x[counter_parts],z[counter_parts],'-o', color="red", alpha=0.65, label="Qox", linewidth=1)   
            plt.legend(loc="upper left", fontsize=6, framealpha=0.25)
            plt.suptitle('Inter- and intrachain contacts in ' + str(key_nq) + " part " + str(counter_parts + 1) + " of " + str(parts))
            pp.savefig()
            counter_parts += 1
            if (args.show_plots == False):
                plt.close('all')  
                                 
pp.close() 
if (args.show_plots):
    log("Finished calculating. Showing plots.",'i')
    plt.show()
else:
    log("Finished calculating. Opening the pdf with the plots.", 'i')
    if platform.system() == 'Linux':
        subprocess.call(('xdg-open', name_pdf_plots)) 


log('Created plots. Exiting plotResResContactChanges.py.', 'i')
log("-- %s seconds ---"% (time.time()- _start_time), 'i')

