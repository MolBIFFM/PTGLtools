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
import shutil


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

def new_directory(name):
    """Creates a new directory 'name' if it does not exist yet. Returns the absolute path of that directory"""
    if not os.path.isdir(name):
        os.makedirs(name)
    return os.path.abspath(name)

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

cl_parser.add_argument('--divide-by-chainlength',
                       metavar = 'divide-by-chainlength',
                       default = '',
                       help = 'specify a path to a "number_of_residues_in_each_chain" csv file to divide the change of each chain by its chain length.')

# Args for the computation of the csv file containing the number of residues in each chain. 
cl_parser.add_argument('--pdb-or-cif-file',
                       metavar = 'pdb-or-cif-file',
                       default = '',
                       help = 'specify a path to a pseudo pdb or cif file. Enter a corresponding dssp file and if you entered a pseudo pdb file a header file as well to create a csv file with the number of residues in each chain.')

cl_parser.add_argument('--dssp-file',
                       metavar = 'dssp-file',
                       default = '',
                       help = 'specify a path to a dssp file corresponding to your cif or pseudo pdb file.')

cl_parser.add_argument('--headerfile',
                       metavar = 'headerfile',
                       default = '',
                       help = 'specify a path to a headerfile for the computation of a mmcif file if you have entered a file in pseudo pdb format.')

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

os.chdir(output_dir)
input_dir_csv_file = new_directory('files_for_GraCom_computation') + '/'
os.chdir(input_dir)

# pdb of cif file
if args.pdb_or_cif_file != '':
    if(os.access(args.pdb_or_cif_file, os.R_OK)) and (args.pdb_or_cif_file.endswith('.pdb') or args.pdb_or_cif_file.endswith('.cif')):        
        shutil.copy(args.pdb_or_cif_file, input_dir_csv_file + os.path.basename(args.pdb_or_cif_file))
        pdb_or_cif_file = input_dir_csv_file + os.path.basename(args.pdb_or_cif_file)
    else:
        log("Can´t read the given pdb or cif file or wrong file format. Continuing without computation of a csv file.", 'e')
        pdb_or_cif_file = ''
else:
    pdb_or_cif_file = ''
    log("No pseudo pdb or mmCIF file given. Can´t create csv file with number of residues in each chain.", 'i')

# dssp file
if args.dssp_file != '':
    if(os.access(args.dssp_file, os.R_OK)) and (args.dssp_file.endswith('.dssp')):
        shutil.copy(args.dssp_file, input_dir_csv_file + os.path.basename(args.dssp_file))
        dssp_file = input_dir_csv_file + os.path.basename(args.dssp_file)
    else:
        log("Can´t read the given dssp file or wrong file format. Continuing without computation of a csv file.", 'e')
        dssp_file = ''
else:
    dssp_file = ''
    log("No dssp file found or given. Can´t create csv file with number of residues in each chain.", 'i')

# headerfile
if args.headerfile != '':
    if(os.access(args.headerfile, os.R_OK)):
        headerfile = args.headerfile
    else:
        log("Can´t read the headerfile. Continuing without computation of a csv file if you entered a pdb file.", 'e')
        headerfile = ''
else:
    headerfile = ''
    log("No headerfile given. Continuing without computation of a csv file if you entered a pdb file.", 'i')

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

# Create csv file with number of residues in each chain.
path_PyMolHeatmap = os.path.dirname(__file__)
path_toMmCif = ((path_PyMolHeatmap.split('plotting-scripts'))[0]) + 'toMmCIF.py'
path_graCom = ((path_PyMolHeatmap.split('PTGLdynamics'))[0]) + 'PTGLgraphComputation/dist/PTGLgraphComputation.jar'

cif_file = ''
if pdb_or_cif_file.endswith('.pdb') and dssp_file != '' and headerfile != '':
    createMmCifFile = 'python3 ' + path_toMmCif + ' -i ' +  input_dir_csv_file + ' -p ' + input_dir_csv_file + ' --headerfile ' + headerfile
    os.chdir(input_dir_csv_file)
    os.system(createMmCifFile)
    os.chdir(work_dir)
    cif_file = os.path.abspath(input_dir_csv_file) + '/' + (os.path.basename(pdb_or_cif_file).split('.pdb'))[0] + '.cif'
else:
    cif_file = pdb_or_cif_file

if cif_file != '' and dssp_file != '':
    createCsvFile = 'java -jar ' + path_graCom + ' none -p ' + cif_file + ' -d ' + dssp_file + ' -o ' + output_dir + ' -I --set "PTGLgraphComputation_B_csv_number_residues_chains" "true" --set "PTGLgraphComputation_B_debug_only_parse" "true"' 
    os.chdir(output_dir)
    os.system(createCsvFile)  
    os.chdir(work_dir)
    log("Csv file with the number of residues in each chain created.", 'i')
    
    

# Read in file with the number of residues in each chain.
if(args.divide_by_chainlength != ""):
    file_chainlength = args.divide_by_chainlength
    if(os.access(file_chainlength, os.R_OK)):
        with open(file_chainlength, 'r') as fc:
            lines_chainlength = fc.read().split("\n")

        for chain in lines_chainlength[1:]:
            if chain != '':
                columns = chain.split(',')
                if (len(columns) == 2):
                    changes_nodes[columns[0]] = int(changes_nodes.get(columns[0])) / int(columns[1])

        log("Changes divided by chain lengths " + str(changes_nodes), 'i')
            
    else:
        logging.error("Specified input file '%s' is not readable. Continuing without it.", file_chainlength)

change_each_chain_divided_chainlength = open(output_dir + '/' + 'change_each_chain_divided_by_chain_length_frame' + str(args.first_timestep) + '_to_frame' + str(args.last_timestep) + '.csv','w')
change_each_chain_divided_chainlength.write("chain" + "," +  "change divided by chain length" + '\n')

for key in changes_nodes:
    change_each_chain_divided_chainlength.write(key + ',' + str(changes_nodes[key]) + '\n')

change_each_chain_divided_chainlength.close()


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

