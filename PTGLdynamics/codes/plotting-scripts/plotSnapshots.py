########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "1.0.0"  # TODO version of this template, change this to 1.0.0 for a new script or 2.0.0 if you upgrade another script to this template's architecture


########### built-in imports ###########

import sys
import os
import argparse
import logging
import traceback
import time
import matplotlib.pyplot as plt
import decimal
from matplotlib.backends.backend_pdf import PdfPages
import platform
import subprocess


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


def check_for_U(name):
    """checks the file name wether it belongs to the dataset with or without U"""
    if('withU' in name):
        ending = '-with_U'
    elif('noU' in name):
        ending = '-noU'
    else:
        ending = ''
    return ending

########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########
          
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="Reads csv with edge weights in columns and creates plot with line for each edge as well as plots for each edge seperately.",
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

cl_parser.add_argument('inputfiles',
                       metavar = 'inputfile',
                       nargs = "+",
                       default = '',
                       help = 'get data as csv file input')

cl_parser.add_argument('-o',
                       '--outputfile',
                       metavar = 'path',
                       default = '',
                       help = 'save results to file')

cl_parser.add_argument('-n',
                       '--name-pdf-plots',
                       metavar = 'name-pdf-plots',
                       default = '',
                       help = 'specify the name of the output pdf where all plots are stored')

"""
########### DOES NOT WORK ##############
cl_parser.add_argument('-l',
                       '--log-scale',
                       action='store_true',
                       help = 'a log scale is used for the y-axis')
"""

cl_parser.add_argument('-m',
                       '--mean-interval',
                       metavar = 'mean-interval',
                       default = 50,
                       type = int,
                       help = 'specify the interval length for the mean-graph. The default value is 50.')

cl_parser.add_argument('-r',
                       '--y-axis-range',
                       metavar = 'y-axis-range',
                       type = float,
                       nargs = 2,
                       default = [],
                       help = 'specify the range of the y-axis in floats.')
                       
cl_parser.add_argument('--show-plots',
                       action='store_true',
                       help = 'open the matplot plots.')                       

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

# input file
csv = []
for elem in args.inputfiles:
    if(os.access(elem, os.R_OK) == True):
        csv.append(elem)
    else:
        logging.error("Specified file '%s' is not readable. Exiting.", elem)
        exit()
if(csv == []):
    logging.error("No readable files found. Exiting now.", args.inputfiles)
    sys.exit(1)

    
# name of pdf output
if(args.name_pdf_plots != ""):
    if(args.name_pdf_plots.endswith('.pdf')):
        name_pdf_plots = args.name_pdf_plots
    else:
        name_pdf_plots = args.name_pdf_plots + '.pdf'
else:
    name_pdf_plots = "multipage.pdf"

# interval length
if(args.mean_interval !=0):
    interval = args.mean_interval

# y axis range
y_axis_range = []
if(args.y_axis_range != [] and len(args.y_axis_range) == 2):
    y_axis_range = args.y_axis_range

########### vamos ###########

log("Version " + version, "i")

_start_time = time.time()

edges = [] # list with list for each file containing list for each edge with all its weights [file][edge][weights]
edge_keys = [] # list with list for each file containing lists with edge labels [file][edge labels]
fig_count = 0

mean = [] # list that stores the mean values calculated for all edge values [file][edge][mean for current interval]


pp = PdfPages(name_pdf_plots)



########## get data from CSVs ############

for i in range(len(csv)): # [file]

    with open(csv[i], "r") as f:
        csv_lines = f.read().split("\n")
    
    if(csv_lines != ""): #create a list in edges and edge_keys for each file
        edges.append([])
        edge_keys.append([])
        mean.append([])
        
    empty_lines_count = 0
    
    for j in range(len(csv_lines)): # [weights]
        if(csv_lines[j].replace('\t','') == ''): # count empty lines at the ending of a file - they are substracted later
            log('empty line','i')
            empty_lines_count +=1
            continue
        
        csv_lines[j] = csv_lines[j].split(",")
        for k in range(len(csv_lines[j])): # [edge]
            if(j==0):
                edge_keys[-1].append([csv_lines[j][k]]) #get the edge names
            elif(j == 1):
                edges[-1].append([decimal.Decimal(csv_lines[j][k])]) #create a new list for each edge in the file
                mean[-1].append([decimal.Decimal(csv_lines[j][k])])
            else:
                edges[-1][k].append(decimal.Decimal(csv_lines[j][k])) #append all entries to the corresponding edge (in the csv each column represents one edge)

                mean[-1][k].append(decimal.Decimal(csv_lines[j][k]))
                if(j<interval):
                    for h in range(len(mean[-1][k])-1): #go through all values for this edge and add the current value (the new appended entry is spared)
                        mean[-1][k][h] += decimal.Decimal(csv_lines[j][k]) # calculate the sum of the first #interval values
                else:
                    for h in range(j-interval,j-1): #add the current value to all values in the interval
                        mean[-1][k][h] += decimal.Decimal(csv_lines[j][k])
                
    #plot each file with all its edges

    x = list(range(len(csv_lines)-empty_lines_count-2))
    w = list(range(int(interval/2), len(csv_lines) - empty_lines_count -int(interval/2)-2))

    plt.figure(fig_count)
    fig_count +=1
    for j in range(len(edges[-1])):
        y = [decimal.Decimal(edges[-1][j][y]) for y in x]

        dataset = check_for_U(csv[i])
        plt.plot(x,y, label=str(edge_keys[-1][j]) + dataset , alpha=0.5)
        
    
    for j in range(len(edges[-1])):

        z =  [decimal.Decimal(mean[-1][j][z-int(interval/2)]/interval) for z in w] #calculate mean for this point

        dataset = check_for_U(csv[i])
        plt.plot(w,z, label=str(edge_keys[-1][j]) + dataset + '-mean' , alpha=0.25)

    """
    ############# DOES NOT WORK ##############
    if(args.log_scale):
        plt.yscale('log')
    """

    plt.xlabel("timestep")
    plt.ylabel("edge weight")

    if(y_axis_range !=[]):
        plt.ylim(y_axis_range)
    #plt.set_title(csv)
    plt.legend(loc="upper left", fontsize=2)
    #plt.show()
    pp.savefig()
    if (args.show_plots == False):
        plt.close('all')
    

########### get edges #############

all_keys = []

for elem in edge_keys:
    for edge in elem:
        if edge in all_keys:
            continue
        else:
            all_keys.append(edge)
        

########## plot data from different files #########


#plot each edge with all file appearances

for i in range(len(all_keys)):

    plt.figure(fig_count)
    fig_count +=1
    plt.xlabel("timestep")
    plt.ylabel("edge weight")
    for j in range(len(edge_keys)):

        if all_keys[i] in edge_keys[j]:
            #log('plotting ' + str(fig_count),'i')
            
            k = edge_keys[j].index(all_keys[i])
            x = list(range(len(edges[j][k])))

            y = [decimal.Decimal(edges[j][k][y]) for y in x]

            dataset = check_for_U(csv[j])
            plt.plot(x,y, label=str(all_keys[i]) + dataset , alpha=0.25)
            #if(args.log_scale):
            #    plt.yscale('log')
            
    for j in range(len(edge_keys)):

        if all_keys[i] in edge_keys[j]:
            #log('plotting ' + str(fig_count),'i')
            
            k = edge_keys[j].index(all_keys[i])
            w = list(range(int(interval/2),len(edges[j][k])-int(interval/2)))

            z = [decimal.Decimal(mean[j][k][z-int(interval/2)]/interval) for z in w] #calculate mean for this point

            dataset = check_for_U(csv[j])
            plt.plot(w,z, label=str(all_keys[i])+ dataset + '-mean' , alpha=0.25)


    """
    ############ DOES NOT WORK #############
    if(args.log_scale):
        plt.yscale('log')
    """
    plt.legend(loc="upper left", fontsize=8)
    if(y_axis_range !=[]):
        plt.ylim(y_axis_range)
    pp.savefig()
    if (args.show_plots == False):
        plt.close('all')

pp.close()
if (args.show_plots):
    log("Finished calculating. Showing plots.",'i')
    plt.show()
else:
    log("Finished calculating. Opening the pdf with the plots.", 'i')
    if platform.system() == 'Linux':
        subprocess.call(('xdg-open', name_pdf_plots))



log("-- %s seconds ---"% (time.time()- _start_time), 'i')
log("All done, exiting plot_snapshots", 'i')

# tidy up
if (args.outputfile != ""):
    output_file.close()

