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
cl_parser = argparse.ArgumentParser(description="Takes several csv files with residue-residue contacts and computes a list of residue contacts for each residue that has contacts.",
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
                       
cl_parser.add_argument('inputdir',
                       metavar = 'inputdir',
                       default = '',
                       help = 'specify a path to a folder containing csv files with all residue-residue contacts between chains.')

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

# input directory
input_dir = check_dir_args(args.inputdir)
      
# output directory
output_dir = check_dir_args(args.outputdirectory)

########### vamos ###########

log("Version " + version, "i")

work_dir = get_working_dir(input_dir)
list_work_dir = os.listdir(work_dir)
list_work_dir = sorted_nicely(list_work_dir)


files = 0
for file in list_work_dir:
    if file.endswith(".csv"):
        files += 1
        matches = {}
        with open(file, "r") as f:
            csv_lines = f.read().split("\n")
        
        for line in csv_lines[1:]:
            if line != '':
                columns = line.split(',')
                if columns[0] != '' and columns[4] !='':
                    key = (columns[0], columns[3])
                    value = matches.get(key)
                    if value == None: 
                        matches[key] = {(columns[4], columns[7])}
                    else:
                        value.update([(columns[4], columns[7])])
                        matches[key] = value       
                if columns[4] != '' and columns[7] !='':
                    key = (columns[4], columns[7])
                    value = matches.get(key)
                    if value == None:                        
                        matches[key] = {(columns[0], columns[3])}
                    else:
                        value.update([(columns[0], columns[3])])
                        matches[key] = value

            else:
                break
        log("Contact list for each residue: " + str(matches), 'i')
        
        name = os.path.basename(file).split(".")
        name = name[0].split("contacts_")
               
        contact_lists= open(output_dir + '/' + 'contact_list_for_each_residue_' + name[1] + '.csv','w')
        contact_lists.write("Residue (PDB_ID | chain)" + "," +  "# Contacts" + "," + "Contact with residues" + '\n')

        for key in matches:
            key_replaced = str(key).replace(', ', '|')
            contacts = 0
            matches_replaced = matches[key]
            matched_residues = ''
            for match in matches_replaced:
                contacts += 1
                matched_residues = matched_residues + '#' + str(match).replace(', ', '|')
            contact_lists.write(key_replaced + ',' + str(contacts) + ',' + matched_residues[1:] + '\n')

        contact_lists.close()      
           
log(str(files) + " files containing the contacts with other residues for each residue computed.", 'i')

