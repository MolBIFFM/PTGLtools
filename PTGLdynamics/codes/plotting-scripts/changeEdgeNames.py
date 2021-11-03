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
                
########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########

            
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="Changes the chain ID to a label ID given by the gml file and creates a csv file. You can change the names in the third column to the names of your choice and run the script again with the csv file given as an argument (-c)",
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

cl_parser.add_argument('inputfile',
                       metavar = 'inputfile',
                       default = '',
                       help = 'get data as a csv file created by "evalEdgeWeights.py". Specify the exact input format with -e or -a. Default is -a.')

cl_parser.add_argument('gmlfile',
                       metavar = 'gml-inputfile',
                       default = '',
                       help = 'specify a path to a gml file created by PTGLgraphComputation, using it for chain_id and label_id.')     
                       
cl_parser.add_argument('-C',
                       '--chainnames-csv-file',
                       metavar = 'chainnames-inputfile',
                       default = '',
                       help = "get chain names in csv file. The third column is used for renaming the chains.")                    

cl_parser.add_argument('-p',
                       '--output-dir',
                       metavar = 'output-directory',
                       default = '',
                       help = 'specify a path to your output files. Otherwise the current folder is used.')


# add mutually exclusive group for  configuration of header
cl_parser_group_input_format = cl_parser.add_mutually_exclusive_group()
cl_parser_group_input_format.add_argument('-e',
                                          '--eval-edge-weights-input',
                                          action='store_true',
                                          help = "Set when the input files are in the original eval_edge_weights.py output format.")

cl_parser_group_input_format.add_argument('-a',
                                          '--all-edges-input',
                                          action='store_true',
                                          help = "Set when the input files are in the 'all_edges' output format from eval_edge_weights.py")



cl_parser.add_argument('-o',
                       '--outputfile',
                       metavar = 'path',
                       default = '',
                       help = 'save results to file')


args = cl_parser.parse_args()

if (args.eval_edge_weights_input == False) and (args.all_edges_input == False):
    args.all_edges_input = True

########### check arguments ###########

# assign log level
log_level = logging.WARNING
if (args.debug):
    log_level = logging.DEBUG
elif (args.verbose):
    log_level = logging.INFO
logging.getLogger().setLevel(log_level)

# output file
if (args.outputfile != ""):
    if(check_file_writable(args.outputfile)):
        output_file = open(args.outputfile, "w")
    else:
        logging.error("Specified output file '%s' is not writable. Exiting now.", args.outputfile)
        sys.exit(1)

# input file
csv = check_input_files(args.inputfile)
    
# csv-file with chainnames
chainnames_csv_file = check_input_files(args.chainnames_csv_file)
       
# gml-file
gml_file = check_input_files(args.gmlfile)

# output directory
output_dir = check_dir_args(args.output_dir)
    
########### vamos ###########

log("Version " + version, "i")

id_label_dict = {}

with open(gml_file, "r") as g:
    node_section = False
    for line in g: 
        if "node [" in line:
            node_section = True
        if "id" in line and node_section == True:
           idline = line.split('id ')            
           new_key = idline[1]
           new_key = new_key[:-1]
        if "label" in line and node_section == True:
           labelline = line.split('"')
           id_label_dict[new_key] = labelline[1]
           node_section = False
                
                
# csv file with chainnames automatically generated, get input from gml-files
auto_chainnames = output_dir + '/' + 'chainnames.csv'
with open(auto_chainnames, 'w') as ac:
    ac.write("gml_id" + "," +  "chain_id" + "," + "given chainname" + '\n')
    for key in id_label_dict:
        ac.write(key + ',' + id_label_dict[key] + ',' + id_label_dict[key] + '\n')   

# use possibly given csv file for naming the chains. 
if (args.chainnames_csv_file != ''):  
    with open(chainnames_csv_file, "r") as c:
        chainnames_lines = c.read().split("\n")

    for chain in chainnames_lines[1:]:
        if chain != '':
            columns = chain.split(',')
            if (len(columns) == 3):
                id_label_dict[columns[0]] = columns[2]
            else:
                log("The given csv file does not consist of three columns. Using the label IDs from the gml file instead.", 'e')
                break 

with open(csv, "r") as f:
    csv_lines = f.read().split("\n")

    
if(args.all_edges_input):
    edge_names = csv_lines[0] # sth like: {0 | 1},{0 | 14},{0 | 16},{0 | 2}
    
    edge_names = edge_names.replace('{','')
    edge_names = edge_names.replace('}','')
    edge_names = edge_names.replace('|','')   
    edge_names = edge_names.split(',') 

    new_line = ""
    
    if (edge_names == ['']):
        log('File contains no edge names. Exiting the script.', 'w')
        exit() 
        
    for j in range(len(edge_names)):
        #new_line += "{"
        edge_names[j] = edge_names[j].split('  ')
        
        for k in range(len(edge_names[j])):
            try:
                edge_names[j][k] = id_label_dict[edge_names[j][k]]
            except KeyError:
                log("Seems like you did not use the corresponding file format and argument. Exiting.", "e")
                sys.exit()
        
        new_edge_name = "{"+edge_names[j][0] + " | " + edge_names[j][1] + "}"

        new_line += new_edge_name
        
        if(j < len(edge_names)-1):
            new_line += ","
        else:
            csv_lines[0] = new_line
            
    # write result to output file
    #csv_name = csv.replace(".csv", "")
    #with open(csv_name+'_new.csv', "w") as f:
    for elem in csv_lines:
        log(elem)

           
elif(args.eval_edge_weights_input):
    for j in range(1, len(csv_lines)):
        if(csv_lines[j] == ''):
            continue
        
        line = csv_lines[j].split(',')

        source = line[0]
        target = line[1]

        log("Source " + str(source) , 'i')
        log("Target " + str(target) , 'i')

        try:
            source = id_label_dict[source]
            target = id_label_dict[target]
        except KeyError:
                log("Seems like you did not use the corresponding file format and argument. Exiting.", "e")
                sys.exit()

        # write result back to csv_lines
        new_line = source + ',' + target
        for k in range(2, len(line)):
            new_line += ',' + line[k]
        csv_lines[j] = new_line
   
    # write result to output file
    #csv_name = csv.replace(".csv", "")
    #with open(csv_name+'_new.csv', "w") as f:
    for elem in csv_lines:
        log(elem)
    

log("All done, exiting changeEdgeNames", 'i')
        
                       

# tidy up
if (args.outputfile != ""):
    output_file.close()
