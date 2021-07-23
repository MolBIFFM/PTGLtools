########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "1.0.1"  


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
                       help='Specify chains that should not be considered in coloring, but in calculation.')
                       
cl_parser.add_argument('--exclude-calculation-chains',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify chains that should not be considered in calculation and coloring.')  
                       
cl_parser.add_argument('--exclude-calculation-edges',
                       type = str,
                       nargs = '+',
                       default = [],
                       help='Specify edges that should not be considered in calculation and coloring. Enter edges as "{edge1  edge2}".')                                            

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


if (args.first_timestep >= args.last_timestep):
    log("The given first frame number is greater than the number given as the last timestep.",'e')
    exit()
    
exclude_chains_coloring = args.exclude_coloring
exclude_calculation_chains = args.exclude_calculation_chains
exclude_calculation_edges = args.exclude_calculation_edges

########### vamos ###########

log("Version " + version, "i")

work_dir = get_working_dir(input_dir)
list_work_dir = os.listdir(work_dir)
list_work_dir = sorted_nicely(list_work_dir)

frames = list(range(args.first_timestep, args.last_timestep))
files_counter = 0
files = []
for frame in frames:
    pattern = "compared*frame"+ str(frame) + "_complex_chains_albelig_CG*frame" + str((frame + 1)) + "_complex_chains_albelig_CG.csv"
    matching = fnmatch.filter(list_work_dir, pattern)
    files = files + matching
    files_counter += 1

log(files, 'i')

if len(files) < files_counter:
    log("There are fewer csv files than expected.", 'i')
if len(files) == 0:
    log("No input data. Exiting.", 'i')
    exit()

chains = []
    
# for each edge: sum of changes
changes_edges = {}

for file in files:
    with open(file, "r") as f:
        csv_lines = f.read().split("\n")

    del csv_lines[0]

    for line in csv_lines:
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
if (exclude_calculation_edges != ''):
    for edge in exclude_calculation_edges:
        if edge in changes_edges:
            nodes = edge.split(' ')
            node_1 = nodes[0]
            node_2 = nodes[1]
            node_1 = node_1.replace(' ', '')
            node_1 = node_1.replace('{', '')
            node_2 = node_2.replace(' ', '')
            node_2 = node_2.replace('}', '')
    
            chains.append(node_1)
            chains.append(node_2)
            
            del changes_edges[edge]   
             
        else:
            log("You entered non valid edge names.", 'i') 

log(changes_edges, 'i')

# Output changes_edges as a csv file
change_each_edge = open(output_dir + '/' + 'change_each_edge_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.csv','w')
change_each_edge.write("edge" + "," +  "change" + '\n')

for key in changes_edges:
    change_each_edge.write(key + ',' + str(changes_edges[key]) + '\n')

change_each_edge.close()

# for each node: sum of changes
changes_nodes = {}
for key in changes_edges:
    nodes = key.split(' ')
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
        old_value_node_1 = changes_nodes.get(node_1)
        old_value_node_2 = changes_nodes.get(node_2)
        if old_value_node_1 == None:
            old_value_node_1 = 0
        if old_value_node_2 == None:
            old_value_node_2 = 0
        new_value_node_1 = int(old_value_node_1) + changes_edges[key]
        new_value_node_2 = int(old_value_node_2) + changes_edges[key]
    
        changes_nodes[node_1] = new_value_node_1
        changes_nodes[node_2] = new_value_node_2
    
log(changes_nodes, 'i')

# Save dictionary as a csv file.
change_each_chain = open(output_dir + '/' + 'change_each_chain_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.csv','w')
change_each_chain.write("chain" + "," +  "change" + '\n')

for key in changes_nodes:
    change_each_chain.write(key + ',' + str(changes_nodes[key]) + '\n')

change_each_chain.close()


# Create PyMol script.
pymol_script = open(output_dir + '/' + 'PyMol_script_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.py', 'w')
pymol_script.write('"""' + '\n' + 'This script shows the given molecule in a heatmap visualisation in PyMOL. Chains colored in blue have lower edge weight changes than chains colored in red. Run this script in PyMOL using the command line with the following command:' + '\n' + 'run ' + output_dir + 'PyMol_script_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.py' + '\n' + '"""' + '\n')
pymol_script.write("cmd.load('" + inputfile + "')" + "\n")


key_max = max(changes_nodes.keys(), key = (lambda k: changes_nodes[k]))
key_min = min(changes_nodes.keys(), key = (lambda k: changes_nodes[k]))
value_max = changes_nodes[key_max]
value_min = changes_nodes[key_min]
change = value_max - value_min

# Exlude specified chains for coloring.
if (exclude_chains_coloring != []):
    for chain in exclude_chains_coloring:   
        if (chain in changes_nodes):
            del changes_nodes[chain]

chains = list(set(chains))

for key in changes_nodes:
    chain = key
    chains.remove(chain)
    value = changes_nodes.get(key)
    percent = (value - value_min) / change
    percent = round(percent, 2)
    if percent == 0.5:
        pymol_script.write("cmd.color('white', 'chain " + chain + "')" + "\n")       

    #left, blue part of gradient
    elif percent < 0.5:
        R = int(255 * percent * 2)
        G = int(255 * percent * 2)
        B = 255
        pymol_script.write("cmd.set_color('color" + chain + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + chain + "', 'chain " + chain + "')" + "\n")

        
    #right, red part of gradient
    elif percent > 0.5:
        R = 255
        G = int(255 * (1.0 - percent) * 2)
        B = int(255 * (1.0 - percent) * 2)
        pymol_script.write("cmd.set_color('color" + chain + "', [" + str(R) + ","
                           + str(G) + "," + str(B) + "])" + "\n")
        pymol_script.write("cmd.color('color" + chain + "', 'chain " + chain + "')" + "\n")

# Remaining uncolored chains
for chain in chains:               
    pymol_script.write("cmd.color('grey', 'chain " + chain + "')" + "\n") 
        

pymol_script.close()

# Run the script in PyMOL.
run_script = 'pymol -q ' + output_dir + '/' + 'PyMol_script_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.py'
os.system(run_script)


log("Finished calculations. Created PyMOL script.", 'i')

