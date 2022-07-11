 
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
import fnmatch

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
    """Returns a dictionary with all changes. """
    changes = {}
    with open(file, 'r') as f:
        csv_lines = f.read().split("\n")
    header = csv_lines[0].split(',')
    header = header[0]
    for line in csv_lines[1:]:
        if line != '':
            columns = line.split(',')
            changes[columns[0]] = columns[1]
    return changes, header
        
        
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
            
########### command line parser ###########
"""Arguments with hyphen like '--input-dir' are called with an underscore within the programm: args.input_dir
However, in the command line call the hyphen is used: --input-dir <path> """

## create the parser
cl_parser = argparse.ArgumentParser(description="Takes two csv files computed by calculateChanges with the absolute changes and calculates the difference.",
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

cl_parser.add_argument('file_two',
                       metavar = 'file-two',
                       help='Enter the path to a second csv file computed by calculateChanges containing the absolute changes.')                                           
                       
cl_parser.add_argument('--exclude-coloring',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains that should not be considered in coloring, but in calculation.')
                       
cl_parser.add_argument('--exclude-calculation-chains',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains that should not be considered in calculation and coloring.')                         

cl_parser.add_argument('-p',
                       '--output-dir',
                       metavar = 'output-directory',
                       default = '',
                       help = 'specify a path to your output directory. Otherwise the current folder is used.')

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

exclude_chains_coloring = args.exclude_coloring
exclude_calculation_chains = args.exclude_calculation_chains

log(f"Excluded chains from calcutation: {exclude_calculation_chains}", "d")

########### vamos ###########

changes_one, header_one =  read_in_file(file_one) 
changes_two, header_two =  read_in_file(file_two)

if header_one != header_two:
    log("Non matching csv files. Enter two csv files in the same format.", "e")
    exit()

log("Changes calculated for the first file :" + str(changes_one), 'i')
log("Changes calculated for the second file :" + str(changes_two), 'i')

changes = changes_one

for key in changes_two:
    chainID = key.split("|")[1].replace("'", "")
    #log(chainID, "d")
    if chainID in exclude_calculation_chains:
        continue
    old_value = changes.get(key)
    if old_value != None:
        new_value = int(old_value) - int(changes_two[key])
    else:
        new_value = 0 - int(changes_two[key])
    changes[key] = new_value
    
log("Changes calculated for both files: " + str(changes), 'i')
        
# Create outputfile
comp_datasets = open(output_dir + '/' + 'compared_datasets.csv','w')
comp_datasets.write(header_one + "," +  "change in percent" + '\n')

for key in changes:
    comp_datasets.write(key + ',' + str(changes[key]) + '\n')

comp_datasets.close()

all_changes = [int(v) for v in set(changes.values())]
value_max = int(max(all_changes))
value_min = int(min(all_changes))
change = value_max - value_min

for key in changes:
    value = changes.get(key)
    percent = (int(value) - value_min) / change
    percent = round(percent, 2)
    changes[key] = percent
    
percents = open(output_dir + '/' + 'compared_datasets_changes_in_percent.csv','w')
percents.write(header_one + "," +  "change in percent" + '\n')

for key in changes:
    percents.write(key + ',' + str(changes[key]) + '\n')

percents.close()

log('finished comparing the datasets and creating the PyMOL script', 'i')



