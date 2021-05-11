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
import pathlib
import shutil
import re
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
    
def check_arguments_args(argument):
    """Checks if the given argument is empty. Returns the given value."""
    if(argument != ""):
        return argument
    else:
        return ""
        
        
def execute_getAttribute(elem, plotting_dir, args, filepath, output_dir, working_dir, ending_of_files):
    """Creates the command line command for a given filepath and executes getAttribute."""
    filename_gml = os.path.basename(filepath)
    filename_csv = filename_gml[:filename_gml.index(ending_of_files)] 
    get_attribute = 'python3 ' + plotting_dir + elem + ' ' + filepath + ' ' + args + ' -o ' + output_dir + filename_csv + 'attribute_extracted.csv'
    log(get_attribute, 'd')
    os.chdir(output_dir)  
    os.system(get_attribute)
    os.chdir(working_dir)          

########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########

            
########### command line parser ###########
"""Arguments with hyphen like '--input-dir' are called with an underscore within the programm: args.input_dir
However, in the command line call the hyphen is used: --input-dir <path> """

## create the parser
cl_parser = argparse.ArgumentParser(description="This is a pipeline that reads several snapshots in pdb-format. " +
                                    "The pipeline runs PTGLgraphComputation on each of them and compares the differences in between their Complex Graph outputs. "+
                                    "If not specified otherwise the results are written in the current folder in csv-format. " +
                                    "Use: 'python3 <path>/ptglDynamics.py <arguments>' to run the programm.",
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

cl_parser.add_argument('-i',
                       '--input-dir',
                       metavar = 'input-directory',
                       default = '',
                       help = 'specify a path to your input files. Otherwise the current folder is used.')

cl_parser.add_argument('-p',
                       '--output-dir',
                       metavar = 'output-directory',
                       default = '',
                       help = 'specify a path to your output files. Otherwise the current folder is used.')

cl_parser.add_argument('-c',
                       '--compoundfile',
                       metavar = 'compoundfile',
                       default = '',
                       help = 'to integrate a header in pdb files specify the path of your compound file.')

cl_parser.add_argument('-a',
                       '--applications',
                       metavar = 'applications',
                       nargs = "*",
                       type = str,
                       default = ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py', 'getAttributeDataFromGml.py'],
                       help = "to execute only the specified scripts. The scripts must be part of the PTGLdynamics set which contains: 'toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py', 'getAttributeDataFromGml.py'")

cl_parser.add_argument('-m',
                       '--dssp-input-dir',
                       metavar = 'dssp-input-dir',
                       default = '',
                       help = 'to specify a different input directory for dssp files if only applications after the dssp and pdb modification are executed.')

cl_parser.add_argument('-u',
                       '--sub-dir-structure',
                       action='store_true',
                       help='display the results in sub directories in the output directory.')

cl_parser.add_argument('--PTGLgraphComputation-path',
                       default = (os.path.dirname(__file__) + '/PTGLgraphComputation/dist/PTGLgraphComputation.jar'),
                       help = 'Absolute path to a custom PTGLgraphComputation JAR file. Otherwise assuming built version of PTGLtools.')

cl_parser.add_argument('-k', 
                       '--PTGLgraphComputation-args',
                       metavar = 'PTGLgraphComputation-args',
                       type = str,
                       default = '',
                       help = 'a string with the PTGLgraphComputation arguments you want to use and its values to execute PTGLgraphComputation in different ways using PTGLgraphComputations command line arguments. Insert arguments like this: -a="<arguments and their inputs>" ')

cl_parser.add_argument('-b',
                       '--toLegacyPDB-args',
                       metavar = 'toLegacyPDB-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for toLegacyPDB you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -b="<arguments and their inputs>" ')

cl_parser.add_argument('-e',
                       '--dsspcmbi-args',
                       metavar = 'dsspcmbi-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for dsspcmbi you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -e="<arguments and their inputs>" ')

cl_parser.add_argument('-f',
                       '--postProcessDssp-args',
                       metavar = 'postProcessDssp-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for postProcessDssp you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -f="<arguments and their inputs>" ')

cl_parser.add_argument('-g',
                       '--gmlCompareEdgeWeightsAndSubsets-args',
                       metavar = 'gmlCompareEdgeWeightsAndSubsets-args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for gmlCompareEdgeWeightsAndSubsets you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -g="<arguments and their inputs>" ')

cl_parser.add_argument('-j',
                       '--different-dssp-folders',
                       action = 'store_true',
                       help = 'saves the dssp from dsspcmbi and the post processed dssp in different folders if the sub directory structure is activated.')

cl_parser.add_argument('-t',
                       '--getAttributeDataFromGml-args',
                       metavar = 'getAttributeDataFromGml-args',
                       default = 'numAllResResContacts',
                       help = 'a string with the arguments for getAttributeDataFromGml you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -t="<arguments and their inputs>". Put the positional attribute argument first, followed by optional arguments. Default attribute argument is numAllResResContacts. ')

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
input_dir = check_dir_args(args.input_dir)
    

# output directory
original_output_dir = check_dir_args(args.output_dir)
    
# compoundfile directory
if (args.compoundfile != ""):
    if(os.access(args.compoundfile, os.R_OK)):
        compound = os.path.abspath(args.compoundfile)
        cmd_compound = ' -c ' + compound
    else:
        logging.error("Specified compound file '%s' is not readable. Continuing without compound file.", args.compoundfile)
        compound = ''
        cmd_compound = ''
else:
    compound = ''
    cmd_compound = ''

# list of applications
programm_list = []
if (args.applications != []):
    for programm in args.applications:
        if programm in ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py', 'getAttributeDataFromGml.py']:
            programm_list.append(programm)
        else:
            logging.error("Specified programm '%s' is not part of the ptglDynamics pipeline. Continuing without it.", programm)
else:
    programm_list = ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py', 'getAttributeDataFromGml.py']

# dssp directory
dssp_input_dir = check_dir_args(args.dssp_input_dir)

# PTGLgraphComputation arguments
add_PTGLgraphComputation_args = check_arguments_args(args.PTGLgraphComputation_args)

# toLegacyPDB arguments
add_toLegacyPDB_args = check_arguments_args(args.toLegacyPDB_args)

# postProcessDssp arguments
add_postProcessDssp_args = check_arguments_args(args.postProcessDssp_args)
    
# gmlCompareEdgeWeightsAndSubsets arguments
add_gml_comparison_args = check_arguments_args(args.gmlCompareEdgeWeightsAndSubsets_args)

# dsspcmbi arguments
add_dsspcmbi_args = check_arguments_args(args.dsspcmbi_args)

# getAttributeDataFromGml arguments
add_getAttributeDataFromGml_args = check_arguments_args(args.getAttributeDataFromGml_args)

    
# different dssp folders
if (args.different_dssp_folders):
    dir_names = {'toLegacyPDB.py':'legacyPDB', 'dsspcmbi':'oldDssp', 'postProcessDssp.py':'newDssp', 'PTGLgraphComputation':'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py': 'gml', 'getAttributeDataFromGml.py': 'csv'}
else:
    dir_names = {'toLegacyPDB.py':'legacyPDB', 'dsspcmbi':'dssp', 'postProcessDssp.py':'dssp', 'PTGLgraphComputation':'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py': 'csv', 'getAttributeDataFromGml.py': 'csv'}

########### vamos ###########

log("Version " + version, "i")


############### Declaration of variables ###############
_start_time = time.time()

ptglDynamics_path = os.path.dirname(__file__)
PTGLgraphComputation_path = args.PTGLgraphComputation_path

cmd_start = 'python3 ' + ptglDynamics_path + '/codes/'

input_dir = os.path.abspath(input_dir) + '/'
original_output_dir = os.path.abspath(original_output_dir) + '/'

pdb_dir = input_dir
plotting_dir = ptglDynamics_path + '/codes/plotting-scripts/'

# define dssp_dir considering that it can be specified via command line
if(dssp_input_dir == os.getcwd()):
    dssp_dir = input_dir
else:
    dssp_dir = os.path.abspath(dssp_input_dir) + '/'

PTGLgraphComputation_dir = input_dir
gml_dir = input_dir

work_dir = get_working_dir(input_dir)
list_work_dir = []

log("work_dir: ", 'd')
log(work_dir, 'd')
log("programm list: ", 'd')
log(programm_list, 'd')


################## Go through the given applications and execute them #####################
for elem in programm_list:
    log("elem: " + elem, 'd')
    
    # Get the output directory
    if (args.sub_dir_structure) and (dir_names[elem] != ''):
        os.chdir(original_output_dir)
        out_dir = new_directory(dir_names[elem]) + '/'
        os.chdir(work_dir)
        
        if(elem == "gmlCompareEdgeWeightsAndSubsets.py"):
            os.chdir(original_output_dir)
            gml_dir = new_directory("gml") + '/'
            os.chdir(work_dir)
        
    else:
        out_dir = original_output_dir
        
    #execute different scripts:
    if (elem == 'toLegacyPDB.py'):

        work_dir = get_working_dir(pdb_dir)
        exec_string = cmd_start + elem + ' ' + add_toLegacyPDB_args + ' -i ' + work_dir + ' -p ' + out_dir + cmd_compound
        log('exec_string ' + exec_string, 'd')
        os.chdir(out_dir)
        os.system(exec_string)
        pdb_dir = os.path.abspath(out_dir) + '/'
        
        log('toLegacyPDB computations are done.', 'i')
        
      
    elif (elem == 'dsspcmbi'):

        work_dir = get_working_dir(pdb_dir)
        
        list_work_dir = os.listdir(work_dir)
        list_work_dir = sorted_nicely(list_work_dir)
        
        for pdb in list_work_dir:
            if pdb.endswith(".pdb"):
                dssp = out_dir + pathlib.Path(pdb).stem + '.dssp'
                exec_string = ptglDynamics_path + '/codes/' + elem + ' ' + work_dir + pdb + ' ' + dssp  + ' ' + add_dsspcmbi_args
                log('exec_string ' + exec_string, 'd')
                os.chdir(out_dir)
                os.system(exec_string)
                os.chdir(work_dir)
        dssp_dir = os.path.abspath(out_dir) + '/'
        
        log('dsspcmbi computations are done.', 'i')

    
    elif (elem == 'postProcessDssp.py'):

        work_dir = get_working_dir(dssp_dir)
        exec_string = cmd_start + elem + ' ' + add_postProcessDssp_args + ' -i ' + work_dir + ' -p ' + out_dir
        log('exec_string ' + exec_string, 'd')
        os.chdir(out_dir)
        os.system(exec_string)
        dssp_dir = os.path.abspath(out_dir) + '/'
        
        log('postProcessDssp computations are done.', 'i')

    
    elif (elem == 'PTGLgraphComputation'):

        work_dir = get_working_dir(pdb_dir)
        
        list_work_dir = os.listdir(work_dir)
        list_work_dir = sorted_nicely(list_work_dir)
        len_work_dir = len(list_work_dir)
        
        counter = 0;
        
        for pdb in list_work_dir:
            if (pdb.endswith(".pdb")):
                
                pdb_id = pathlib.Path(pdb).stem
                pdb_id_folder = ''
                if(args.sub_dir_structure):
                    os.chdir(out_dir)
                    pdb_id_folder = new_directory(pdb_id)
                    os.chdir(work_dir)
                else:
                    pdb_id_folder = out_dir

                dssp = dssp_dir + pathlib.Path(pdb).stem + '.dssp'

                PTGLgraphComputation = 'java -jar ' + PTGLgraphComputation_path + ' ' + pdb_id + ' ' + add_PTGLgraphComputation_args + ' -p ' + work_dir + pdb + ' -d ' + dssp + ' -o ' + pdb_id_folder

                log(PTGLgraphComputation,'d') 
                os.chdir(out_dir)
                os.system(PTGLgraphComputation)
                os.chdir(work_dir)
                
                counter += 1
                if(counter % 5 == 0):
                    log(str(counter) + ' / ' + str(len_work_dir) + ' file computations with PTGLgraphComputation are done.', 'i')
        
        PTGLgraphComputation_dir = os.path.abspath(out_dir) + '/'
        
        
    elif (elem == 'gmlCompareEdgeWeightsAndSubsets.py'):

        work_dir = get_working_dir(PTGLgraphComputation_dir)
        log('PTGLgraphComputation_dir: ' + PTGLgraphComputation_dir,'d')
        prevGml = ''
        list_work_dir = os.listdir(work_dir)
        list_work_dir = sorted_nicely(list_work_dir)
        for entry in list_work_dir:
            log(entry, 'd')
            
            if (os.path.isdir(entry)):

                temp_work_dir = get_working_dir(entry)
                list_temp_work_dir = os.listdir(temp_work_dir)
                log(list_temp_work_dir,'d')
                
                list_temp_work_dir = sorted_nicely(list_temp_work_dir);
                for gml in list_temp_work_dir:
                    log('in dir: ' + gml, 'd')
                    if gml.endswith("complex_chains_albelig_CG.gml"):
                        log(gml_dir,'d')
                        try:
                            shutil.copy(gml, gml_dir + gml)
                        except shutil.SameFileError:
                            log("Source and destination represents the same file.", 'i')
                        
                        log(prevGml,'d')
                        if prevGml != '':
                            os.chdir(gml_dir)
                            gmlComparison = cmd_start + elem +' ' + add_gml_comparison_args + ' -a ' + prevGml + ' -b ' + gml + ' -p ' + out_dir
                            log(gmlComparison,'d')
                            os.system(gmlComparison)
                            os.chdir(temp_work_dir)

                        prevGml = gml
                os.chdir('../')
            else:
                log('not in dir: ' + entry, 'd')
                if entry.endswith("complex_chains_albelig_CG.gml"):
                    try:
                        shutil.copy(entry, gml_dir + entry)
                    except shutil.SameFileError:
                        log("Source and destination represents the same file.", 'i')


                    if prevGml != '':
                        os.chdir(gml_dir)
                        gml_comparison = cmd_start + elem +' ' + add_gml_comparison_args + ' -a ' + prevGml + ' -b ' + entry + ' -p ' + out_dir
                        log(gml_comparison,'d')
                        os.system(gml_comparison)
                        os.chdir(work_dir)
                    prevGml = entry
        
        log('gmlCompareEdgeWeightsAndSubsets computations are done.', 'i')


    elif (elem == 'getAttributeDataFromGml.py'):
        ending_of_files = 'complex_chains_albelig_CG.gml'

        # use gml-files computed by PTGLGraphComputation.    
        work_dir = get_working_dir(PTGLgraphComputation_dir)        
        list_work_dir = os.listdir(work_dir)
        list_work_dir = sorted_nicely(list_work_dir)  

        for data in list_work_dir:

            if (os.path.isdir(data)):   
                temp_work_dir = get_working_dir(data)
                list_temp_work_dir = os.listdir(temp_work_dir)
                list_temp_work_dir = sorted_nicely(list_temp_work_dir)
 
                for file in list_temp_work_dir:
                    if file.endswith(ending_of_files):
                        filepath = temp_work_dir + file
                        execute_getAttribute(elem, plotting_dir, add_getAttributeDataFromGml_args, filepath, out_dir, temp_work_dir, ending_of_files)

                os.chdir('../')

            elif (data.endswith(ending_of_files)):  
                datapath = work_dir + data
                execute_getAttribute(elem, plotting_dir, add_getAttributeDataFromGml_args, datapath, out_dir, work_dir, ending_of_files)
    
        getAttributeDataFromGml_dir = os.path.abspath(out_dir) + '/'
        log('getAttributeDataFromGml computations are done.', 'i')
    

log("-- %s seconds ---"% (time.time()- _start_time), 'i')
log("All done, exiting ptglDynamics.", 'i')

# tidy up
if (args.outputfile != ""):
    output_file.close()

