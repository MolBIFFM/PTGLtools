########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "2.0.0"  


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


def get_working_dir(new_dir):
    """Changes the working directory to the given path and returns the new working directory"""
    os.chdir(new_dir)
    return os.getcwd() + '/'

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.
    creates a list for each file consisting of the different int and string parts of the name
    afterwards the file list is sorted considering those changed names only
 
    Required arguments:
    l -- The iterable to be sorted.
    
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)
            
            
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
           
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="Compare the complex graphs of consecutive timesteps and sum up the edge weights to create a csv file and a PyMOL script with a heatmap visualisation to display changes in chains edge weights.",
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
                       
cl_parser.add_argument('inputpdb',
                       metavar = 'path',
                       default = '',
                       help = 'specify a path to a pdb or mmCIF file to load into PyMOL.')

cl_parser.add_argument('inputcsv',
                       metavar = 'path',
                       default = '',
                       help = 'specify a path to a csv file containing the changes for each chain or residue.')
                       
cl_parser.add_argument('--exclude-chains',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains that should not be considered in coloring.')
                       
cl_parser.add_argument('--exclude-residues',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify residues that should not be considered in coloring. Enter residues in the following format: "(PDB_ID|Chain|iCode)".')                                                                 

cl_parser.add_argument('-p',
                       '--outputdirectory',
                       metavar = 'path',
                       default = '',
                       help = 'specify a path to your output files. Otherwise the current folder is used.')

args = cl_parser.parse_args()


########### check arguments ###########

# assign log level
log_level = logging.WARNING
if (args.debug):
    log_level = logging.DEBUG
elif (args.verbose):
    log_level = logging.INFO
logging.getLogger().setLevel(log_level)

# input pdb
inputfile = check_input_files(args.inputpdb)

# input csv
csv_file = check_input_files(args.inputcsv)
      
# output directory
output_dir = check_dir_args(args.outputdirectory)

exclude_chains = args.exclude_chains
exclude_residues = args.exclude_residues

########### vamos ###########

log("Version " + version, "i")

changes = {}

with open(csv_file, 'r') as c:
    csv_lines = c.read().split("\n")
header = csv_lines[0].split(',')
if header[0] == "residue (PDB_ID|chain|iCode)":
    header[0] = "residue"
    
for line in csv_lines[1:]:
    if line != '':
        columns = line.split(',')
        changes[columns[0]] = columns[1]
        
# Create PyMol script.
pymol_script = open(output_dir + '/' + 'PyMol_script_Heatmap_visualisation_of_each_' + header[0] + '.py', 'w')
pymol_script.write('"""' + '\n' + 'This script shows the given molecule in a heatmap visualisation in PyMOL: A ' + '\n' + header[0] + ' coloured in blue has fewer changes than a ' + header[0] + 'coloured in red.' + '\n' + 'The calculation is based on the sum of ' + header[1] + '.' + '\n' + 'Run this script in PyMOL using the command line with the following command:' + '\n' + 'run ' + output_dir + '/' +'PyMol_script_Heatmap_visualisation_of_each_' + header[0] + '.py' + '\n' + '"""' + '\n')
pymol_script.write("cmd.load('" + inputfile + "')" + "\n")
pymol_script.write("cmd.color('grey', 'all')" + "\n")


all_changes = set(changes.values())
value_max = int(max(all_changes))
value_min = int(min(all_changes))
change = value_max - value_min

# Exlude specified chains and residues from coloring.
if (exclude_chains != []):
    for chain in exclude_chains:   
        if (chain in changes):
            del changes[chain]
if (exclude_residues != []):
    for residue in exclude_residues:
        if (residue in changes):
            del changes[residue]
            
for key in changes:
    if header[0] == "residue":
        parts = key.replace('(', '').replace(')','').replace("'", "").split('|')
    value = changes.get(key)
    percent = (int(value) - value_min) / change
    percent = round(percent, 2)
    if percent == 0.5:
        if header[0] == "chain":
            pymol_script.write("cmd.color('white', 'chain " + key + "')" + "\n")
        elif header[0] == "residue (PDB_ID|chain|iCode)":
            pymol_script.write("cmd.color('white', 'chain " + str(parts[1]) + " and resi " + str(parts[0]) + str(parts[2]) + "')" + "\n") 
            
    else:
        if percent < 0.5:
            R = int(255 * percent * 2)
            G = int(255 * percent * 2)
            B = 255

        else:
            R = 255
            G = int(255 * (1.0 - percent) * 2)
            B = int(255 * (1.0 - percent) * 2)
            
        if header[0] == "chain":
            pymol_script.write("cmd.set_color('color" + key + "', [" + str(R) + "," + str(G) + "," + str(B) + "])" + "\n")
            pymol_script.write("cmd.color('color" + key + "', 'chain " + key + "')" + "\n")
            
        elif header[0] == "residue":
            pymol_script.write("cmd.set_color('color" + str(parts[1]) + str(parts[0]) + str(parts[2]) + "', [" + str(R) + "," + str(G) + "," + str(B) + "])" + "\n")
            pymol_script.write("cmd.color('color" + str(parts[1]) + str(parts[0]) + str(parts[2]) + "', 'chain " + str(parts[1]) + " and resi " + str(parts[0]) + str(parts[2]) + "')" + "\n")
      

pymol_script.close()

# Run the script in PyMOL.
run_script = 'pymol -q ' + output_dir + '/' + 'PyMol_script_Heatmap_visualisation_of_each_' + header[0] + '.py'
os.system(run_script)


log("Finished calculations. Created PyMOL script.", 'i') 
