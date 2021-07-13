########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "1.0.2"


########### built-in imports ###########

import sys
import os
import argparse
import logging
import traceback


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


########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
        
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="This script returns all occurences of a specified argument in a given gml file in csv format.",
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

# add general command line arguments
cl_parser.add_argument('--version',
                       action='version',
                       version='%(prog)s ' + version)

# add positional arguments
cl_parser.add_argument('inputfile',
                       metavar = 'inputfile',
                       default = '',
                       help = 'gml file to process.')

cl_parser.add_argument('attribute_name',
                       metavar = 'attribute_name',
                       default = '',
                       help = 'search the gml for all values of this argument')

# add optional arguments
cl_parser.add_argument('-o',
                       '--outputfile',
                       metavar = 'path',
                       default = '',
                       help = 'save results to file')

cl_parser.add_argument('-t',
                       '--type-to-consider',
                       metavar = 'type-to-consider',
                       default = 'b',
                       help = "only consider arguments from 'n' = nodes, 'e' = edges, 'b' = default = both")

# add mutually exclusive group for  configuration of header
cl_parser_group_header_configuration = cl_parser.add_mutually_exclusive_group()
cl_parser_group_header_configuration.add_argument('-c',
                                                  '--complete-output',
                                                  action='store_true',
                                                  help = "Print attribute names and not only values")

cl_parser_group_header_configuration.add_argument('--omit-header',
                                                  action = 'store_true',
                                                  help = 'Omit printing any header, i.e. only data is printed')

args = cl_parser.parse_args()

########### check arguments ###########

# assign log level
log_level = logging.WARNING
if (args.debug):
    log_level = logging.DEBUG
elif (args.verbose):
    log_level = logging.INFO
logging.getLogger().setLevel(log_level)

# input file
if(os.access(args.inputfile, os.R_OK)):
    gml_input = args.inputfile
else:
    logging.error("Specified input file '%s' is not readable. Exiting now.", args.inputfile)
    sys.exit(1)

# output file
if (args.outputfile != ""):
    if(check_file_writable(args.outputfile)):
        output_file = open(args.outputfile, "w")
    else:
        logging.error("Specified output file '%s' is not writable. Exiting now.", args.outputfile)
        sys.exit(1)

if (args.type_to_consider != "b"):
    if(args.type_to_consider in ["n","e","b"]):
        type_to_consider = args.type_to_consider
    else:
        logging.error("The specified type '%s' is not valid. Exiting now.", args.type_to_consider)
else:
    type_to_consider = "b"
        
########### vamos ###########

log("Version " + version, "i")


### read file in

log("Reading file in", 'i')
with open(gml_input, "r") as f:
    gml_lines = f.read().split("\n")
    
    
### parse information

attribute_name = args.attribute_name

node_section = False
edge_section = False
list_edges = [] # contains element of form [ source, target, attribute ]
list_nodes = [] # contains element of form [  , attribute ]

tmp_id = ""
tmp_source = ""
tmp_target = ""
tmp_attribute = ""

for line in gml_lines:
    if "node [" in line:
        node_section = True
    if "edge [" in line:
        edge_section = True
    elif line.replace("\t", "") == "]":
        if node_section:
            if(tmp_attribute != ""):
                list_nodes.append([tmp_id, "", tmp_attribute])
            node_section = False
            tmp_id = ""
            tmp_attribute = ""
        elif edge_section:
            if(tmp_attribute != ""):
                list_edges.append([tmp_source, tmp_target, tmp_attribute])
            edge_section = False
            tmp_source = ""
            tmp_target = ""
            tmp_attribute = ""
    else:
        # no control statements
        tmp_fields = line.replace("\t","").split(" ", 1)
        if node_section:
            if tmp_fields[0] == "id":
                tmp_id = tmp_fields[1]
            if tmp_fields[0] == attribute_name:
                tmp_attribute = tmp_fields[1].lstrip()
        elif edge_section:
            if tmp_fields[0] == "source":
                tmp_source = tmp_fields[1].lstrip()
            if tmp_fields[0] == "target":
                tmp_target = tmp_fields[1].lstrip()
            if tmp_fields[0] == attribute_name:
                tmp_attribute = tmp_fields[1].lstrip()


### print results

log("Creating result file", 'i')
log("output_file: '" + args.outputfile + "'", "d")
log("type_to_consider: '" + type_to_consider + "'", 'd')

if (not args.omit_header):
    if(args.complete_output):
        log('id/source,target,' + attribute_name)
    else:
        log(attribute_name)
    
if(type_to_consider in ['n','b']):
    for i in range(len(list_nodes)):
        if(args.complete_output):
            log(list_nodes[i][0] + ',' + list_nodes[i][1] + ',' + list_nodes[i][2])
        else:
            log(list_nodes[i][2])
if(type_to_consider in['e','b']):
    for i in range(len(list_edges)):
        if(args.complete_output):
            log(list_edges[i][0] + ',' + list_edges[i][1] + ',' + list_edges[i][2])
        else:
            log(list_edges[i][2])


### tidy up
if (args.outputfile != ""):
    output_file.close()

