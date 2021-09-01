 
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
import math
import copy
import decimal
import time



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
            
            
def new_directory(name):
    """Creates a new directory 'name' if it does not exist yet. Returns the absolute path of that directory"""
    if not os.path.isdir(name):
        os.makedirs(name)
    return os.path.abspath(name)

            
def get_working_dir(new_dir):
    """Changes the working directory to the given path and returns the new working directory"""
    os.chdir(new_dir)
    return os.getcwd() + '/'


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

## create the parser
cl_parser = argparse.ArgumentParser(description="This script reads multiple csv files containing source, target and attribute values. The values of entries with similar source and target entries are evaluated and the minimum, maximum, mean and medium values are calculated. Optional a csv file with all values of the same source and target entries in the same column can be created.",
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

cl_parser.add_argument('-o',
                       '--outputfile',
                       metavar = 'path',
                       default = '',
                       help = 'save results to file')

cl_parser.add_argument('inputdir',
                       metavar = 'input-directory',
                       #nargs = '*',
                       default = '',
                       help = 'specify a path to the directory where you store your input csv files.')

cl_parser.add_argument('-p',
                       '--output-dir',
                       metavar = 'output-directory',
                       default = '',
                       help = 'specify a path to your output files. Otherwise the current folder is used.')


cl_parser.add_argument('-e',
                       '--all-edge-value-output',
                       action='store_true',
                       help = 'store a csv with all edge values of the attribute')

args = cl_parser.parse_args()


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
        
# input directory
input_dir = check_dir_args(args.inputdir)

# output directory
output_dir = check_dir_args(args.output_dir)

    
########### vamos ###########

log("Version " + version, "i")

_start_time = time.time()

output_dir = os.path.abspath(output_dir) + '/'

################# iterate over csv files #####################

counter = 0
edge_list = [] # stores edges as: [source, target, min, max, edge_values_sum -> calculate mean later, [edge_values] -> calculate median later]

work_dir = get_working_dir(input_dir)

for csv in os.listdir(work_dir):
    
    if(csv.endswith(".csv") == False):
        continue
    
    #read in file if readable
    if(os.access(csv, os.R_OK) == False):
        logging.error("Specified file '%s' is not readable. Continuing without it.", csv)
        continue
    
    
    with open(csv, "r") as f:
        csv_lines = f.read().split("\n")
    
    
    #clone edge_list:
    cloned = copy.deepcopy(edge_list)
    
    for line in csv_lines[1:]:
        line = line.split(",")
        
        if(len(line) < 2):
            break
        
        edge_found = False

        
        for i in range(len(edge_list)): #iterate through edges from previous files

            
            for j in range(len(edge_list[i])): # compare one edge to line

                if(line[j] == edge_list[i][j]):
                    
                    if(j == 1): #similar source and target
                        
                        edge_found = True
                        line[2] = decimal.Decimal(line[2])
                        if(line[2] < edge_list[i][2]):
                            cloned[i][2] = line[2]
                        if(line[2] > edge_list[i][3]):
                            cloned[i][3] = line[2]
                        cloned[i][4] += line[2]
                        cloned[i][5].append(line[2])

                        break
                    
                else: # line and edge are different
                    break
        
        if(edge_found == False):
            line[2] = decimal.Decimal(line[2])
            if(counter == 0):
                cloned.append([line[0],line[1],line[2],line[2], line[2], [line[2]]])
            else:
                cloned.append([line[0],line[1],decimal.Decimal(0.0),line[2], line[2], []])
                for i in range(counter-1):
                    cloned[-1][5].append(decimal.Decimal(0.0))
                cloned[-1][5].append(line[2])
    
    
    for elem in cloned:
        if(len(elem[5]) <= counter):
            elem[5].append(decimal.Decimal(0.0))
            elem[2] = decimal.Decimal(0.0)
    

    edge_list = cloned
    
    counter += 1


# calculate mean and median for each edge

median_list_item = ((counter-1) / 2)
median_list_item = math.modf(median_list_item) # stores a number as ( decimal place , integer place )
log(median_list_item, 'd')
log(edge_list, 'd')

edges = []

edge_list = sorted(edge_list, key = lambda l:l[0]+l[1])
log('Median','d')
for elem in edge_list:
    
    edges.append(elem[5]) # get all edge weights in an extra list for another csv
    
    
    elem[4] = decimal.Decimal(elem[4] / (counter))
    
    elem[5] = sorted(elem[5])

    if(median_list_item[0] != 0.0):
        median_1 = elem[5][int(median_list_item[1]+1)]
        median_2 = elem[5][int(median_list_item[1])]
        log(median_1, 'd')
        log(median_2,'d')
        median = decimal.Decimal((median_1 + median_2) / 2)
    else:
        median = elem[5][int(median_list_item[1])]
    elem.pop()
    elem.append(median)

log(edges,'d')
log("Printing results", 'i')

if(args.all_edge_value_output):
    edges_out = open(output_dir + 'all_edge_values.csv','w')

log('source,target,min,max,mean,median')
for i in range(len(edge_list)):
    log(edge_list[i][0] + ',' + edge_list[i][1] + ',' + str(edge_list[i][2]) + ',' + str(edge_list[i][3]) + ',' + str(edge_list[i][4]) + ',' + str(edge_list[i][5]))
    
    if(args.all_edge_value_output):
        edges_out.write("{" + edge_list[i][0] + " | " + edge_list[i][1] + "}")
        if(i < len(edge_list)-1):
                edges_out.write(",")
        else:
            edges_out.write("\n")

if(args.all_edge_value_output):
    for i in range(counter):
        for j in range(len(edges)):
            
            edges_out.write(str(edges[j][i]))
            if(j < len(edges)-1):
                edges_out.write(",")
            else:
                edges_out.write("\n")


    edges_out.close()

log("-- %s seconds ---"% (time.time()- _start_time), 'i')
log("All done, exiting eval_edges", 'i')


# tidy up
if (args.outputfile != ""):
    output_file.close()


