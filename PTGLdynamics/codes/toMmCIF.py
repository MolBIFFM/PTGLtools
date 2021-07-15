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
import time
import fileinput
import re
import string


########### functions ###########


def check_file_writable(fp):
    """Checks if the given filepath is writable"""
    if os.path.exists(fp):
        if os.path.isfile(fp):
            return os.access(fp, os.W_OK)
        else:
            return False
    # target does not exist, check perms on parent dir
    parent_dir = os.path.dirname(fp)
    if not parent_dir: parent_dir = '.'
    return os.access(parent_dir, os.W_OK)


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
    
    

########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
            
########### command line parser ###########
"""Arguments with hyphen like '--input-dir' are called with an underscore within the programm: args.input_dir
However, in the command line call the hyphen is used: --input-dir <path> """

## create the parser
cl_parser = argparse.ArgumentParser(description="Creates a mmCIF file out of a pseudo legacyPDB file. ",
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

cl_parser.add_argument('-f',
                       '--starting-line',
                       metavar = 'starting line',
                       type = int,
                       default = 1,
                       help = 'Enter the starting line for the _atom_site loop.')

cl_parser.add_argument('--headerfile',
                       metavar = 'headerfile',
                       default = '',
                       help='The headerfile is used as a header for the mmCIF files.')

cl_parser.add_argument('-i',
                       '--input-dir',
                       metavar = 'input-directory',
                       default = '',
                       help = 'specify a path to your input directory. Otherwise the current folder is used.')

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

# headerfile
if (args.headerfile != ""):
    if(os.access(args.headerfile, os.R_OK)):
        headerfile = args.headerfile
    else:
        logging.error("Specified headerfile '%s' is not readable. Continuing without headerfile.", args.headerfile)
        headerfile = ''
else:
    headerfile = ''

# input directory
input_dir = check_dir_args(args.input_dir)
    

# output directory
output_dir = check_dir_args(args.output_dir)


########### vamos ###########

log("Version " + version, "i")

_start_time = time.time()
output_dir = os.path.abspath(output_dir) + '/'
os.chdir(input_dir)
input_dir = os.getcwd() + '/'


#reading .txt file which is the header of a file
header = ''
if headerfile != '':
    f = open(headerfile, "r")
    header = f.read()
    header = header 

#amino acids (residues)in column 18-20 
aaCode = [ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

loop_mmCIF = ("loop_" + '\n' +
              "_atom_site.group_PDB" + '\n' +
              "_atom_site.id" + '\n' +
              "_atom_site.type_symbol" +'\n'+
              "_atom_site.label_atom_id" +'\n'+
              "_atom_site.label_alt_id" +'\n'+
              "_atom_site.label_comp_id" +'\n'+
              "_atom_site.label_asym_id" +'\n'+
              "_atom_site.label_entity_id" +'\n'+
              "_atom_site.label_seq_id" +'\n'+
              "_atom_site.pdbx_PDB_ins_code" +'\n'+
              "_atom_site.Cartn_x" +'\n'+
              "_atom_site.Cartn_y" +'\n'+
              "_atom_site.Cartn_z" +'\n'+
              "_atom_site.occupancy" +'\n'+
              "_atom_site.B_iso_or_equiv" +'\n'+
              "_atom_site.pdbx_formal_charge" +'\n'+
              "_atom_site.auth_seq_id" +'\n'+
              "_atom_site.auth_comp_id" +'\n'+
              "_atom_site.auth_asym_id" +'\n'+
              "_atom_site.auth_atom_id" +'\n'+
              "_atom_site.pdbx_PDB_model_num" +'\n')


log(__file__, 'i')

#get number of files in os.listdir(input_dir) to avoid going through the modified ones several times
number_files = len(os.listdir(input_dir))

list_input_dir = os.listdir(input_dir)
list_input_dir = sorted_nicely(list_input_dir)

for file in list_input_dir:
    if file.endswith(".pdb"):
        data = []
        with open(file) as f:
            for line in f:
                data.append(line)
        
        os.chdir(output_dir)
        output = open(file, "w")

        del data[0:args.starting_line]

        line = str(header) + str(loop_mmCIF) 
        output.write(line)
        

        currChainID = '' 
        prevChainID = ''
        chain_counter = 0

        for line in data:
                
            if line[0:3] == 'END':
                line = '#'
                output.write(line)
                
            if 'ATOM' == line[0:4]:
            
                if line[72:75] == "CHA":
                    chainID = line[75:76]
                else: 
                    chainID = line[72]                 

                currChainID = chainID
                    
                if (prevChainID != currChainID):
                    chain_counter += 1
                               

                # Add column label_alt_id
                line = line[:17] + '. ' + line[17::]

                # Add pdbx_PDB_ins_code
                line = line[:29] + '?' + line[29::]

                # blanks between x,y,z
                line = line[:41] + ' ' + line[41:49]+ ' ' + line[49::]

                
                # Add pdbx_formal_charge
                line = line[:72] + '?' + line[72::]

                
                # Add auth_seq_id
                auth_seq_id = line[25:29]
                line = line[:75] + auth_seq_id +  line[75::]
                
                # auth_comp_id
                auth_comp_id = line[19:23]
                line = line[:79] + auth_comp_id + line[79::]
                
               
                # Add pdbx_PDB_model_num
                line = line[:92] + '  1' + line[92::]
                

                #Add label_asym_id
                label_asym_id = dict(zip(range(1,27), string.ascii_uppercase))                
                line = line[:23] + '  ' + label_asym_id[chain_counter] + ' ' + line[24::]                              

                # auth_atom_id
                type_symbol = line[13:16]
                atom = line[13]
                if (line[13:15] == 'FE'):
                    atom = line[13:15]                    
                line = line[:93] + ' ' + type_symbol + ' ' + line[96::]
                
                #Delete CHA
                if line[89:92] == 'CHA':
                    line = line[:89] + line[92::]
                else:
                    line = line[:90] + line[92::]

                # entity_id
                asym_to_entity_id = {'1':'1', '2':'2', '3':'3', '4':'4', '5':'5', '6':'6', '9':'7', 'A':'10', 'J':'11', 'K':'12', 'L':'13', 'M':'14', 'N':'15', 'H':'16', '7':'8', 'W':'9', 'F':'18', 'U': '20'}
                line = line[:28] + ' ' + str(asym_to_entity_id[chainID]) + '  ' + line[28::]
                
                # Add column type_symbol
                line = line[:12] + ' ' + atom + ' ' + line[12::]                
              
                if line[22:25] not in aaCode:
                    line = 'HETATM '+ line[6:]

                output.write(line)
                prevChainID = currChainID

        output.close()
        old_file = os.path.abspath(file)
        new_file = old_file[:-3]+ 'cif'
        os.rename(old_file, new_file)
        os.chdir(input_dir)

log('finish pdb file overwritting', 'i')
log("-- %s seconds ---"% (time.time()- _start_time), 'i')

