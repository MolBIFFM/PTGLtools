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
import matplotlib.pyplot as plt
import decimal
from matplotlib.backends.backend_pdf import PdfPages

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
cl_parser = argparse.ArgumentParser(description="The script calculates the sum of all entries in the same row of a csv as well as a mean graph for the sum graph and outputs the result with matplotlib",
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
if(os.access(args.inputfile, os.R_OK) == True):
    csv = args.inputfile
else:
    logging.error("Specified file '%s' is not readable.", args.inputfile)
    sys.exit(1)

# name of pdf output
if(args.name_pdf_plots != ""):
    if(args.name_pdf_plots.endswith('.pdf')):
        name_pdf_plots = args.name_pdf_plots
        legend = args.name_pdf_plots.replace('.pdf','')
    else:
        name_pdf_plots = args.name_pdf_plots + '.pdf'
        legend = args.name_pdf_plots
else:
    name_pdf_plots = "edge_sum.pdf"
    legend = "edge_sum"

# interval length
if(args.mean_interval !=0):
    interval = args.mean_interval

# y axis range
y_axis_range = []
if(args.y_axis_range != [] and len(args.y_axis_range) == 2):
    y_axis_range = args.y_axis_range

########### vamos ###########

log("Version " + version, "i")
edge_weights = [] # [line][0==timestep ; 1==sum]
mean = []

with open(csv, "r") as f:
    csv_lines = f.read().split("\n")

for j in range(len(csv_lines)):
    if(j == 0):
        edge_weights.append([])
        edge_weights[j].append("timestep")
        edge_weights[j].append("sum")
    else:
        if(csv_lines[j] != ''):
            edge_weights.append([])
            weights = csv_lines[j].split(',')
            log("weights list",'i')
            log(weights,'i')
            edge_weights[j].append(j)
            edge_weights[j].append(0)
            for k in range(len(weights)): #edge weights at that time step
                #log(type(elem),'i')
                #log(elem,'i')
                elem = decimal.Decimal(weights[k])
                edge_weights[j][1] += elem
            mean.append(edge_weights[j][1])
            log('mean:','i')
            log(mean,'i')
            if(j<interval):
                for h in range(len(mean)-1): #go through all values and add the current value (the new appended entry is spared)
                    mean[h] += edge_weights[j][1] # calculate the sum of the first #interval values
            else:
                for h in range(j-interval,j-1): #add the current value to all values in the interval
                    mean[h] += edge_weights[j][1]
            
            

# plotting
pp = PdfPages(name_pdf_plots)

plt.figure(figsize = (20,10))

x = list(range(len(edge_weights)-1))
y = [edge_weights[y+1][1] for y in x]

plt.plot(x,y, label=legend , alpha=0.5)

w = list(range(int(interval/2), len(edge_weights) -int(interval/2)-1))
z =  [decimal.Decimal(mean[z-int(interval/2)]/interval) for z in w] #calculate mean for this point

plt.plot(w,z, label=legend+'-mean')

plt.xlabel("timestep")
plt.ylabel("sum of contacts")
plt.legend(fontsize = 'small')
if(y_axis_range !=[]):
    plt.ylim(y_axis_range)
pp.savefig()
plt.close('all')

pp.close()
#log("Finished calculating. Showing plots.",'i')
#plt.show()

# write result to output file
for elem in edge_weights:
    log(str(elem[0]) + ',' + str(elem[1]))

# tidy up
if (args.outputfile != ""):
    output_file.close()

