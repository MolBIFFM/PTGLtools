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

def get_data_from_csv_files(file):
    """ Reads in the data from a csv file and returns the included chains and the change in edgeweights for each chain. """
    chains = []
    changes = {}
    with open(file, "r") as f:
        lines = f.read().split("\n")
        del lines[0]
    
        for line in lines:
            if line != '':
                columns = line.split(',')            
                nodes = columns[0].split(' ')
                node_1 = nodes[0]            
                node_2 = nodes[1]            
                node_1 = node_1.replace(' ', '')
                node_1 = node_1.replace('{', '')
                node_2 = node_2.replace(' ', '')
                node_2 = node_2.replace('}', '')
    
                chains.append(node_1)
                chains.append(node_2)
            
                # Exclude chains from calculation
                if (node_1 in exclude_calculation_chains) or (node_2 in exclude_calculation_chains):
                    pass
                    
                else:
                    old_value_node_1 = changes.get(node_1)
                    old_value_node_2 = changes.get(node_2)
                    if old_value_node_1 == None:
                        old_value_node_1 = 0
                    if old_value_node_2 == None:
                        old_value_node_2 = 0
                    new_value_node_1 = old_value_node_1 + int(columns[1])
                    new_value_node_2 = old_value_node_2 + int(columns[1])
    
                    changes[node_1] = new_value_node_1
                    changes[node_2] = new_value_node_2  
    f.close()
    return chains, changes
    
    
        
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
            
########### command line parser ###########
"""Arguments with hyphen like '--input-dir' are called with an underscore within the programm: args.input_dir
However, in the command line call the hyphen is used: --input-dir <path> """

## create the parser
cl_parser = argparse.ArgumentParser(description="Takes the two csv files 'change_each_edge' created by 'createPymolScript' and creates a new Pymol script to display in which dataset more changes occured. If the change of edgeweights for a chain is higher in the first file than in the second file, the chain will be colored red, otherwise blue.",
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
                       help = 'Enter the path to csv file 1 in "changes_each_edge" format.')

cl_parser.add_argument('file_two',
                       metavar = 'file-two',
                       help='Enter the path to the second csv file in "changes_each_edge" format.')
                       
cl_parser.add_argument('inputfile',
                       metavar = 'inputfile',
                       default = '',
                       help = 'specify a path to a pdb or mmCIF file to load into PyMOL.')                                              
                       
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
inputfile = check_input_files(args.inputfile)  

# output directory
output_dir = check_dir_args(args.output_dir)

exclude_chains_coloring = args.exclude_coloring
exclude_calculation_chains = args.exclude_calculation_chains

########### vamos ###########

chains_one, changes_one =  get_data_from_csv_files(file_one) 
chains_two, changes_two =  get_data_from_csv_files(file_two) 

log("Changes calculated for the first file :" + str(changes_one), 'i')
log("Changes calculated for the second file :" + str(changes_two), 'i')

changes = changes_one

for key in changes_two:
    old_value = changes.get(key)
    if old_value != None:
        new_value = old_value - changes_two[key]
    else:
        new_value = 0 - changes_two[key]
    changes[key] = new_value
    
log("Changes calculated for both files: " + str(changes), 'i')
        
# Create outputfile
comp_datasets = open(output_dir + '/' + 'compared_datasets.csv','w')
comp_datasets.write("chain" + "," +  "change" + '\n')

for key in changes:
    comp_datasets.write(key + ',' + str(changes[key]) + '\n')

comp_datasets.close()

# Create PyMol script.
pymol_script = open(output_dir + '/' + 'PyMOL_script_compared_datasets.py', 'w')
pymol_script.write('"""' + '\n' + 'This script shows the variation in edgeweights of chains of two different datasets. Chains are colored in different shades of red, if the first entered file contains higher edgeweight changes than the second one. Blue colores symbolize the other way around. Run this script in PyMOL using the command line with the following command:' + '\n' + 'run ' + output_dir + 'PyMOL_script_compared_datasets.py' + '\n' + '"""' + '\n')
pymol_script.write("cmd.load('" + inputfile + "')" + "\n")


not_in_both_molecules = list(set(chains_one) - set(chains_two)) + list (set(chains_two) - set(chains_one))
for chain in not_in_both_molecules:
    if chain in changes:
        del changes[chain]
   
chains = list(set(chains_one + chains_two)) 


key_max = max(changes.keys(), key = (lambda k: changes[k]))
key_min = min(changes.keys(), key = (lambda k: changes[k]))
value_max = changes[key_max]
value_min = changes[key_min]


for key in changes:
    chain = key
    chains.remove(chain)
    value = changes.get(key)
    if value == 0:
        pymol_script.write("cmd.color('white', 'chain " + chain + "')" + "\n") 
        
    elif (value > 0):
        percent = round((value / value_max), 2)
        R = 255
        G = int(255 * (1.0  - percent))
        B = int(255 * (1.0 - percent))
        pymol_script.write("cmd.set_color('color" + chain + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + chain + "', 'chain " + chain + "')" + "\n")
        
    elif (value < 0):
        percent = round((value / value_min), 2)
        R = int(255 * (1.0 - percent))
        G = int(255 * (1.0 - percent))
        B = 255
        pymol_script.write("cmd.set_color('color" + chain + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + chain + "', 'chain " + chain + "')" + "\n")
    
# Remaining uncolored chains
for chain in chains:               
    pymol_script.write("cmd.color('grey', 'chain " + chain + "')" + "\n") 
        

pymol_script.close()

log('finished comparing the datasets and creating the PyMOL script', 'i')


