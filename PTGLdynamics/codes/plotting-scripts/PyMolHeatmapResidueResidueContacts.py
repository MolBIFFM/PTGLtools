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
cl_parser = argparse.ArgumentParser(description="Runs a PyMOL script with a heatmap visualisation. Residues are coloured according to their number of changes in residue-residue contacts with residues of other chains. Residues having a high number of contact changes are coloured red, residues with few changes blue.",
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
                       
cl_parser.add_argument('inputfile',
                       metavar = 'inputfile',
                       default = '',
                       help = 'specify a path to a pdb or mmCIF file to load into PyMOL.')

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
                       
cl_parser.add_argument('--exclude-coloring',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains whose residues should not be considered in coloring, but in calculation.')                       

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

# input file
inputfile = check_input_files(args.inputfile)
      
# output directory
output_dir = check_dir_args(args.outputdirectory)

exclude_res_coloring = args.exclude_coloring

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
    pattern = "contact_lists_for_each_residue_compared*frame"+ str(frame) + "*frame" + str((frame + 1)) + ".csv"
    matching = fnmatch.filter(list_work_dir, pattern)
    files = files + matching
    files_counter += 1

log(files, 'i')

if len(files) < files_counter:
    log("There are fewer csv files than expected.", 'i')
if len(files) == 0:
    log("No input data. Exiting.", 'i')
    exit()
    
    
contact_changes = {}
for file in files:
    with open(file, "r") as f:
        csv_lines = f.read().split('\n')
        
    for line in csv_lines[1:]:
        if line != '':
            columns = line.split(',')
            old_value = contact_changes.get(columns[0])
            if (old_value == None):
                contact_changes[columns[0]] = int(columns[1])
            else:
                contact_changes[columns[0]] = int(old_value) + int(columns[1])
        else:
            break

log("Changes in residue-residue contacts: " + str(contact_changes), 'i')  

# Create output csv file.
res_res_contact_changes = open(output_dir + '/' + 'changes_in_res_res_contacts_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.csv','w')
res_res_contact_changes.write("residue (PDB_ID | Chain_ID)" + "," +  "change in res-res contacts" + '\n')

for key in contact_changes:
    res_res_contact_changes.write(key + ',' + str(contact_changes[key]) + '\n')

res_res_contact_changes.close()  

# Exlude residues from specified chains from coloring.
if (exclude_res_coloring != []):
    for element in exclude_res_coloring:
        element = str(element)
    to_delete = []
    for key in contact_changes:
        parts_key = key.split('|')
        chain_id = parts_key[1].replace(')','')
        chain_id = chain_id.replace("'", "")
        if chain_id in exclude_res_coloring:
            to_delete.append(key)

for item in to_delete:
    del contact_changes[item]

pymol_script = open(output_dir + '/' + 'PyMol_script_res_res_contacts_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.py', 'w')   

pymol_script.write('"""' + '\n' + 'This script shows the given molecule in a heatmap visualisation in PyMOL. Residues are coloured according to their number of changes in residue-residue contacts with residues of other chains. Residues having a high number of contact changes are coloured red, residues with few changes blue. Residues having no contacts with residues of other chains are coloured grey. Run this script in PyMOL using the command line with the following command:' + '\n' + 'run ' + output_dir +  'PyMol_script_res_res_contacts_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.py' + '\n' + '"""' + '\n') 

pymol_script.write("cmd.load('" + inputfile + "')" + "\n")   
pymol_script.write("cmd.color('grey', 'all')" + "\n")

all_changes = set(contact_changes.values())
value_max = max(all_changes)
value_min = min(all_changes)
change = value_max - value_min

for residue in contact_changes:
    res = residue.replace('(', '')
    res = res.replace(')', '')
    parts = res.split('|')
    pdb_id = parts[0].replace("'", "")
    chain_id = parts[1].replace("'", "")
    value = contact_changes.get(residue)
    percent = (value - value_min) / change
    percent = round(percent, 2)
    if percent == 0.5:
        pymol_script.write("cmd.color('white', 'chain " + chain_id + " and resi " + pdb_id + "')" + "\n")       

    #left, blue part of gradient
    elif percent < 0.5:
        R = int(255 * percent * 2)
        G = int(255 * percent * 2)
        B = 255
        pymol_script.write("cmd.set_color('color" + str(chain_id) + str(pdb_id) + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + str(chain_id) + str(pdb_id) + "', 'chain " + chain_id + " and resi " + pdb_id + "')" + "\n")

        
    #right, red part of gradient
    elif percent > 0.5:
        R = 255
        G = int(255 * (1.0 - percent) * 2)
        B = int(255 * (1.0 - percent) * 2)
        pymol_script.write("cmd.set_color('color" + str(chain_id) + str(pdb_id) + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + str(chain_id) + str(pdb_id) + "', 'chain " + chain_id + " and resi " + pdb_id + "')" + "\n")
 
pymol_script.close()

#run_script = 'pymol -q ' + output_dir + '/' + 'PyMol_script_res_res_contacts_frame' + str(args.first_timestep) + '_to_frame' + #str(args.last_timestep) + '.py'
#os.system(run_script)



log("Finished calculations. Created PyMOL script.", 'i')

