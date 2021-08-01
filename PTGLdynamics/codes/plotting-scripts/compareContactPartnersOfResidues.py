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
            
            

def check_input_files(inputfile):
    """Returns the argument if the file is readable. Otherwise raises an error and exits."""
    if(inputfile != ""):
        if(os.access(inputfile, os.R_OK)):
            return inputfile
        else:
            logging.error("Specified input file '%s' is not readable. Exiting now.", inputfile)
            sys.exit(1)
    else:
        logging.error("This program needs two input files to work. Exiting now.")
        
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
        
def read_file_in(file):
    """Reads a csv file in and returns sets of residue contacts for each residue in the file."""
    with open(file, "r") as f:
        csv_lines = f.read().split("\n")
    
    residues_contacts = {}    
    for line in csv_lines[1:]:
        if line != '':
            columns = line.split(',')            
            residues = set(columns[2].split('#'))
            residues_contacts[columns[0]] = residues

        else:
            break
         
    return residues_contacts
          

########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
           
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="Compares two csv files containing all contact partners of a residue. Returns the number of changes in contact partners in a csv file.",
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

cl_parser.add_argument('-file1',
                       '--inputfile1',
                       metavar = 'path',
                       help='the first csv file containing the residue contacts for a residue.')

cl_parser.add_argument('-file2',
                       '--inputfile2',
                       metavar = 'path',
                       help='the second csv file containing the residue contacts for a residue.')

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

# input file 1
file1 = check_input_files(args.inputfile1)
        
# input file 2
file2= check_input_files(args.inputfile2)
        
# output directory
output_dir = check_dir_args(args.outputdirectory)

########### vamos ###########

log("Version " + version, "i")

residues_contacts_file1 = read_file_in(file1)
log("Residue contacts in file 1" + str(residues_contacts_file1), 'i')
residues_contacts_file2 = read_file_in(file2)
log("Residue contacts in file 2" + str(residues_contacts_file2), 'i')

symmetric_difference = {}

for key in residues_contacts_file1:
    if key in residues_contacts_file2:
        residues_file1 = residues_contacts_file1.get(key)
        residues_file2 = residues_contacts_file2.get(key)
        difference_residues = residues_file1.symmetric_difference(residues_file2)
        difference = len(difference_residues)
        symmetric_difference[key] = difference
        del residues_contacts_file2[key]
    else:
        symmetric_difference[key] = len(residues_contacts_file1.get(key))
        
for key in residues_contacts_file2:
    symmetric_difference[key] = len(residues_contacts_file2.get(key))
        
log("Change in contacts " + str(symmetric_difference), 'i')

name_file1 = os.path.basename(file1).split(".")
name_file1 = name_file1[0].split("contact_list_for_each_residue_")

name_file2 = os.path.basename(file2).split(".")
name_file2 = name_file2[0].split("contact_list_for_each_residue_")
               
change_in_contacts = open(output_dir + '/' + 'contact_lists_for_each_residue_compared_' + name_file1[1] + '_' + name_file2[1] + '.csv','w')
change_in_contacts.write("Residue (PDB_ID | chain)" + "," +  "# Contact changes" + '\n')

for key in symmetric_difference:
    change_in_contacts.write(key + ',' + str(symmetric_difference.get(key)) + '\n')

change_in_contacts.close() 

log("Finished comparing of two csv files.", 'i')

