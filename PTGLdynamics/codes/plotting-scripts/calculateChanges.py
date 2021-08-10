 
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
    
def divide_by_chainlength(changes, csv_file):
    with open(csv_file, 'r') as f:
        lines_chainlength = f.read().split("\n")

    for chain in lines_chainlength[1:]:
        if chain != '':
            columns = chain.split(',')
            if (len(columns) == 2):
                changes[columns[0]] = int(changes.get(columns[0])) / int(columns[1])  
                    
    return changes   
            
            
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
           
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="Calculates the changes for each chain or residue in a molecule. Chain change calculations can be done based on CG´s or residue-residue contacts.",
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
                       metavar = 'input-directory',
                       default = '',
                       help = 'specify a path to the directory where you store input csv files comparing two timesteps.')

cl_parser.add_argument('calculation',
                       metavar = 'tuple',
                       choices = ["(chain, res)", "(chain, CG)", "(res, res)"],
                       default = '',
                       help = 'Choose which changes should be calculated: Change for each chain based on CG´s: "(chain, CG)", change for each chain based on residue-residue contacts: "(chain, res)", change for each residue based on residue-residue contacts: "(res, res)".')

cl_parser.add_argument('-f',
                       '--first-timestep',
                       type = int,
                       default = 1,
                       help='The first complex graph to be compared')

cl_parser.add_argument('-l',
                       '--last-timestep',
                       type = int,
                       default = 5,
                       help='The last complex graph to be compared')

cl_parser.add_argument('--exclude-chains',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains that should not be considered in calculation and coloring. Only used in (chain, CG) calculation.')  
                       
cl_parser.add_argument('--exclude-edges',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify edges that should not be considered in calculation and coloring. Enter edges as "{edge1  edge2}". Only used in (chain, CG) calculation.')
                       
cl_parser.add_argument('--divide-by-chainlength',
                       metavar = 'path',
                       default = '',
                       help = 'specify a path to a "number_of_residues_in_each_chain" csv file to divide the change of each chain by its chain length. Only used in (chain, CG).')                       

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

# calculation
calculation = args.calculation
calculate = str(calculation.split(',')[0].replace('(', ''))
based_on = str(calculation.split(',')[1].replace(')', '').replace(' ', ''))
      
# output directory
output_dir = check_dir_args(args.outputdirectory)

# additional arguments
exclude_chains = args.exclude_chains
exclude_edges = args.exclude_edges
divide_chainlength = check_input_files(args.divide_by_chainlength)
print(divide_chainlength)

if (args.first_timestep >= args.last_timestep):
    log("The given first frame number is greater than the number given as the last timestep.",'e')
    exit()

########### vamos ###########
log("Version " + version, "i")

work_dir = get_working_dir(input_dir)
list_work_dir = os.listdir(work_dir)
list_work_dir = sorted_nicely(list_work_dir)

frames = list(range(args.first_timestep, args.last_timestep))
files_counter = 0
files = []
    
for frame in frames:
    pattern = ''
    if based_on == 'CG':
        pattern = "compared*frame"+ str(frame) + "_complex_chains_albelig_CG*frame" + str((frame + 1)) + "_complex_chains_albelig_CG.csv"
    elif based_on == 'res':
        pattern = "change_contacts_for_each_residue_compared*frame"+ str(frame) + "*frame" + str((frame + 1)) + ".csv"
    matching = fnmatch.filter(list_work_dir, pattern)
    files = files + matching
    files_counter += 1

log(files, 'i')

if len(files) < files_counter:
    log("There are fewer csv files than expected.", 'i')
if len(files) == 0:
    log("No input data. Exiting.", 'i')
    exit()

changes = {}
if based_on == 'res':
    for file in files:
        with open(file, "r") as f:
            csv_lines = f.read().split('\n')
        
        for line in csv_lines[1:]:
            if line != '':
                columns = line.split(',')
                old_value = changes.get(columns[0])
                if (old_value == None):
                    changes[columns[0]] = int(columns[1])
                else:
                    changes[columns[0]] = int(old_value) + int(columns[1])
            else:
                break
            
    if calculate == 'res':
        log("Changes in residue-residue contacts: " + str(changes), 'i')
    elif calculate == 'chain':
        changes_chains = {}
        for key in changes:
            parts = key.replace("(", "").replace(")", "").replace("'", "").split("|")
            value = changes.get(key)

            if parts[1] in changes_chains:
                changes_chains[parts[1]] = changes_chains.get(parts[1]) + value
            else:
                changes_chains[parts[1]] = value
                
        changes = changes_chains

        log("Changes each chain based on residue-residue contacts: " + str(changes), 'i')

if calculation == '(chain, CG)':
    # for each edge: sum of changes
    changes_edges = {}

    for file in files:
        with open(file, "r") as f:
            csv_lines = f.read().split("\n")

        for line in csv_lines[1:]:
            columns = line.split(',')
            if columns[0] != '' and columns[3] !='':
                old_value = changes_edges.get(columns[0])
                if old_value == None:
                    old_value = 0
                new_value = int(old_value) + int(columns[3])
                changes_edges[columns[0]] = new_value
            else:
                break

    # Exclude edges from calculation.
    if (exclude_edges != []):
        for edge in exclude_edges:
            if edge in changes_edges:
                del changes_edges[edge]   
             
            else:
                log("You entered non valid edge names.", 'i')

    log("Changes for each edge: " + str(changes_edges), 'i')

    # Calculate change for each chain
    for key in changes_edges:
        nodes = key.split(' ')
        node_1 = nodes[0].replace(' ', '').replace('{', '')
        node_2 = nodes[1].replace(' ', '').replace('}', '')

    
        # Exclude chains from calculation
        if (node_1 in exclude_chains) or (node_2 in exclude_chains):
            pass
        
        else:
            old_value_node_1 = changes.get(node_1)
            old_value_node_2 = changes.get(node_2)
            if old_value_node_1 == None:
                old_value_node_1 = 0
            if old_value_node_2 == None:
                old_value_node_2 = 0
            new_value_node_1 = int(old_value_node_1) + changes_edges[key]
            new_value_node_2 = int(old_value_node_2) + changes_edges[key]
    
            changes[node_1] = new_value_node_1
            changes[node_2] = new_value_node_2
            
    # Divide by chainlength
    if(divide_chainlength != ''):
        changes_divided = divide_by_chainlength(changes, divide_chainlength)
        changes = changes_divided        

    log("Changes for each chain: " + str(changes), 'i')
  
# Create output csv file.
changes_out = open(output_dir + '/' + 'changes_each_' + calculate + '_based_on_' + based_on + '_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.csv','w')

if calculation == "(chain, res)":
    changes_out.write("chain" + "," +  "change in res-res contacts" + '\n')
elif calculation == "(res, res)":
    changes_out.write("residue (PDB_ID|chain|iCode)" + "," +  "change in res-res contacts" + '\n')
elif calculation == "(chain, CG)":
    changes_out.write("chain" + "," +  "change in CG´s" + '\n')

for key in changes:
    changes_out.write(key + ',' + str(changes[key]) + '\n')

changes_out.close()  

log("'calculate_changes.py' computations done.", 'i')
