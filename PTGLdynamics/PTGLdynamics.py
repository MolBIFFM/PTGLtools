########### settings ###########

# This script's version as MAJOR.MINOR.PATCH
#   major: big (re-)implementation, new user interface, new overall architecture
#   minor: new functions, git merge
#   patch: fixes, small changes
#   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
# -> only increment with commit / push / merge not while programming
version = "1.1.0"  # TODO version of this template, change this to 1.0.0 for a new script or 2.0.0 if you upgrade another script to this template's architecture  


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

# Add new subscripts here
programs = ['toLegacyPDB.py', 'toMmCIF.py', 'dsspcmbi', 'postProcessDssp.py', 'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py', 'getAttributeDataFromGml.py', 'evalEdgesWeights.py', 'changeEdgeNames.py', 'sumEdgeWeights.py', 'plotSnapshots.py', 'PyMolHeatmapVisualisation.py']

default_path_graCom = os.path.dirname(__file__)
parts = default_path_graCom.split('/PTGLdynamics')
default_path_graCom = parts[0] + '/PTGLgraphComputation/dist/PTGLgraphComputation.jar'


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

cl_parser.add_argument('-H',
                       '--headerfile',
                       metavar = 'headerfile',
                       default = '',
                       help = 'to integrate a header in your mmCif files specify the path of your header file.')
                       
cl_parser.add_argument('-C',
                       '--compoundfile',
                       metavar = 'compoundfile',
                       default = '',
                       help = 'to integrate a header in your pdb files containing 3D structural data specify the path of your header file.')                       

cl_parser.add_argument('-a',
                       '--applications',
                       metavar = 'applications',
                       nargs = "*",
                       type = str,
                       default = list(programs),
                       help = "to execute only the specified scripts. The scripts must be part of the PTGLdynamics set which contains: " + str(programs))

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
                       default = default_path_graCom,
                       help = 'Absolute path to a custom PTGLgraphComputation JAR file. Otherwise assuming built version of PTGLtools.')

cl_parser.add_argument('--PTGLgraphComputation-args',
                       metavar = 'PTGLgraphComputation-args',
                       type = str,
                       default = '',
                       help = 'a string with the PTGLgraphComputation arguments you want to use and its values to execute PTGLgraphComputation in different ways using PTGLgraphComputations command line arguments. Insert arguments like this: -a="<arguments and their inputs>" ')

cl_parser.add_argument('--toLegacyPDB-args',
                       metavar = 'toLegacyPDB-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for toLegacyPDB you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -b="<arguments and their inputs>" ')

cl_parser.add_argument('--dsspcmbi-args',
                       metavar = 'dsspcmbi-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for dsspcmbi you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -e="<arguments and their inputs>" ')

cl_parser.add_argument('--postProcessDssp-args',
                       metavar = 'postProcessDssp-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for postProcessDssp you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -f="<arguments and their inputs>" ')

cl_parser.add_argument('--gmlCompareEdgeWeightsAndSubsets-args',
                       metavar = 'gmlCompareEdgeWeightsAndSubsets-args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for gmlCompareEdgeWeightsAndSubsets you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -g="<arguments and their inputs>" ')

cl_parser.add_argument('-j',
                       '--different-dssp-folders',
                       action = 'store_true',
                       help = 'saves the dssp from dsspcmbi and the post processed dssp in different folders if the sub directory structure is activated.')
                       
cl_parser.add_argument('-I',
                       '--toMmCIF',
                       action='store_true',
                       help = 'converts non legacy PDB files to mmCIF format.')
                       
cl_parser.add_argument('--toMmCIF-args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for toMmCIF you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: --toMmCIF-args="<arguments and their inputs>" ')

cl_parser.add_argument('--getAttributeDataFromGml-args',
                       metavar = 'getAttributeDataFromGml-args',
                       default = 'numAllResResContacts -c',
                       help = 'a string with the arguments for getAttributeDataFromGml you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -t="<arguments and their inputs>". Put the positional attribute argument first, followed by optional arguments. Default attribute argument is numAllResResContacts. ')
                       
cl_parser.add_argument('--evalEdgesWeights-args',
                       metavar = 'evalEdgesWeights-args',
                       default = '',   
                       help = 'a string with the arguments for evalEdgesWeights you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -z="<arguments and their inputs>". For default are two csv files in your output folder created, all_edge_values.csv and edges_weights.csv. If the second file name is altered, changeEdgeNames.py has to be executed manually on that file.')
                       
cl_parser.add_argument('--changeEdgeNames-args',
                       metavar = 'changeEdgeNames-args',
                       default = '',   
                       help = 'a string with the arguments for changeEdgeNames you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -y="<arguments and their inputs>". Using the files from the inputfolder as default. Executing the script alone needs at least all positional arguments, only one csv at a time can be executed.')
                       
cl_parser.add_argument('--sumEdgeWeights-args',
                       metavar = 'sumEdgeWeights-args',
                       default = '',   
                       help = 'a string with the arguments for sumEdgeWeights you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -x="<arguments and their inputs>".')

cl_parser.add_argument('--plotSnapshots-args',
                       metavar = 'plotSnapshots-args',
                       default = '',   
                       help = 'a string with the arguments for plotSnapshots you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -n="<arguments and their inputs>".')
                       
cl_parser.add_argument('--pyMolHeatmapVisualisation-args',
                       metavar = 'pyMolHeatmapVisualisation-arguments',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for createPymolScript you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: --createPymolScript-args="<arguments and their inputs>", including the positional arguments.')                       

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
    
# headerfile directory
if (args.headerfile != ""):
    if(os.access(args.headerfile, os.R_OK)):
        header_mmcif = os.path.abspath(args.headerfile)
        cmd_header_mmcif = ' --headerfile ' + header_mmcif
    else:
        logging.error("Specified header file '%s' is not readable. Continuing without header file. PTGLgraphComputation will not be working properly though.", args.headerfile)
        header_mmcif = ''
        cmd_header_mmcif = ''
else:
    header_mmcif = ''
    cmd_header_mmcif = ''
    if (args.toMmCIF):
        logging.error("No header file was specified. PTGLgraphComputation will not be working properly.") 
    
# compoundfile directory
if (args.compoundfile != ""):
    if(os.access(args.compoundfile, os.R_OK)):
        header = os.path.abspath(args.compoundfile)
        cmd_header = ' -c ' + header
    else:
        logging.error("Specified header file '%s' is not readable. Continuing without header file.", args.compoundfile)
        header = ''
        cmd_header = ''
else:
    header = ''
    cmd_header = ''    

# list of applications
program_list = []
if (args.applications != []):
    for program in args.applications:
        if program in programs:            
            program_list.append(program)
        else:
            logging.error("Specified program '%s' is not part of the ptglDynamics pipeline. Continuing without it.", program)
else:
    program_list = list(programs)
    
if (args.toMmCIF):
    program_list.remove('toLegacyPDB.py')
else:
    program_list.remove('toMmCIF.py')


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

# toMmCIF arguments
add_toMmCIF_args = check_arguments_args(args.toMmCIF_args)

# getAttributeDataFromGml arguments
add_getAttributeDataFromGml_args = check_arguments_args(args.getAttributeDataFromGml_args)
if ('-c' not in add_getAttributeDataFromGml_args):
    log('Not using -c as an argument for getAttributeDataFromGml may cause errors in later scripts.', 'w')

# evalEdgesWeights arguments
add_evalEdgesWeights_args = check_arguments_args(args.evalEdgesWeights_args)

# changeEdgeNames arguments
add_changeEdgeNames_args = check_arguments_args(args.changeEdgeNames_args) 

# sumEdgeWeights arguments
add_sumEdgeWeights_args = check_arguments_args(args.sumEdgeWeights_args)

# plotSnapshots arguments
add_plotSnapshots_args = check_arguments_args(args.plotSnapshots_args)

# createPymolScript arguments
add_pyMolHeatmapVisualisation_args = check_arguments_args(args.pyMolHeatmapVisualisation_args)
    
# different dssp folders
if (args.different_dssp_folders):
    dir_names = {'toLegacyPDB.py':'legacyPDB', 'toMmCIF.py':'mmCIF', 'dsspcmbi':'oldDssp', 'postProcessDssp.py':'newDssp', 'PTGLgraphComputation':'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py': 'gml', 'getAttributeDataFromGml.py': 'csv', 'evalEdgesWeights.py':'csv', 'changeEdgeNames.py':'csv', 'sumEdgeWeights.py':'pdf', 'plotSnapshots.py':'pdf', 'PyMolHeatmapVisualisation.py':'PyMOL'}
else:
    dir_names = {'toLegacyPDB.py':'legacyPDB', 'toMmCIF.py':'mmCIF', 'dsspcmbi':'dssp', 'postProcessDssp.py':'dssp', 'PTGLgraphComputation':'PTGLgraphComputation', 'gmlCompareEdgeWeightsAndSubsets.py': 'csv', 'getAttributeDataFromGml.py': 'csv', 'evalEdgesWeights.py':'csv', 'changeEdgeNames.py':'csv', 'sumEdgeWeights.py':'pdf', 'plotSnapshots.py':'pdf', 'PyMolHeatmapVisualisation.py':'PyMOL'}

########### vamos ###########

log("Version " + version, "i")


############### Declaration of variables ###############
_start_time = time.time()

ptglDynamics_path = os.path.dirname(__file__)
PTGLgraphComputation_path = args.PTGLgraphComputation_path

cmd_start = 'python3 ' + ptglDynamics_path + '/codes/'

input_dir = os.path.abspath(input_dir) + '/'
original_output_dir = os.path.abspath(original_output_dir) + '/'

file_dir = input_dir
plotting_dir = ptglDynamics_path + '/codes/plotting-scripts/'

# define dssp_dir considering that it can be specified via command line
if(dssp_input_dir == os.getcwd()):
    dssp_dir = input_dir
else:
    dssp_dir = os.path.abspath(dssp_input_dir) + '/'

PTGLgraphComputation_dir = input_dir
gml_dir = input_dir
getAttributeDataFromGml_dir = input_dir
evalEdgesWeights_dir = input_dir
changeEdgeNames_dir = input_dir
compareSubsets_dir = input_dir

work_dir = get_working_dir(input_dir)
list_work_dir = []

log("work_dir: ", 'd')
log(work_dir, 'd')
log("program list: ", 'd')
log(program_list, 'd')


################## Go through the given applications and execute them #####################
for elem in program_list:
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

        work_dir = get_working_dir(file_dir)
        exec_string = cmd_start + elem + ' ' + add_toLegacyPDB_args + ' -i ' + work_dir + ' -p ' + out_dir + cmd_header
        log('exec_string ' + exec_string, 'd')
        os.chdir(out_dir)
        os.system(exec_string)
        file_dir = os.path.abspath(out_dir) + '/'
        
        log('toLegacyPDB computations are done.', 'i')
        
    elif (elem == 'toMmCIF.py'):

        work_dir = get_working_dir(file_dir)
        exec_string = cmd_start + elem + ' ' + add_toMmCIF_args + ' -i ' + work_dir + ' -p ' + out_dir + cmd_header_mmcif
        log('exec_string ' + exec_string, 'd')
        os.chdir(out_dir)
        os.system(exec_string)
        file_dir = os.path.abspath(out_dir) + '/'
        
        log('toMmCIF computations are done.', 'i')
        
      
    elif (elem == 'dsspcmbi'):

        work_dir = get_working_dir(file_dir)
        
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

        work_dir = get_working_dir(file_dir)
        
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
                    
        compareSubsets_dir = os.path.abspath(out_dir) + '/'
        
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
        
    elif (elem == 'evalEdgesWeights.py'):
        ending_of_files = 'attribute_extracted.csv'
        work_dir = get_working_dir(getAttributeDataFromGml_dir)
        os.chdir(getAttributeDataFromGml_dir)
        input_dir_evalEdges = new_directory('attribute_csv_files') + '/' 
        os.chdir(work_dir)
        list_work_dir = os.listdir(work_dir)
        list_work_dir = sorted_nicely(list_work_dir) 
         
        for file in list_work_dir:
            if file.endswith(ending_of_files):
                try:
                    shutil.copy(file, input_dir_evalEdges + file)
                except shutil.SameFileError:
                        log("Source and destination represents the same file.", 'i')
                        
        if (len(os.listdir(input_dir_evalEdges)) != 0):
             
            if (add_evalEdgesWeights_args == ''): 
                add_evalEdgesWeights_args = '-e -o ' + out_dir + 'edges_weights.csv'
                
            evalEdgesWeights = 'python3 ' + plotting_dir + elem + ' ' + input_dir_evalEdges + ' ' + add_evalEdgesWeights_args  
            log(evalEdgesWeights, 'd')
            os.chdir(out_dir) 
            os.system(evalEdgesWeights)
            os.chdir(work_dir)
            
        else:
            os.rmdir(input_dir_evalEdges)
        
        evalEdgesWeights_dir = os.path.abspath(out_dir) + '/'
        
        log('evalEdgesWeights computations are done.', 'i')
    
    
    elif (elem == 'changeEdgeNames.py'): 
        if (add_changeEdgeNames_args == ''):  
            # get gml-file
            gml_path = get_working_dir(gml_dir)
            gml_files = os.listdir(gml_dir)
            gml_files = sorted_nicely(gml_files)
            gml_file = ''
            for gml in gml_files:
                    if gml.endswith("complex_chains_albelig_CG.gml"):
                        gml_file = os.path.abspath(gml)
                        break
            work_dir = get_working_dir(evalEdgesWeights_dir)
            log(work_dir, 'd')
            list_work_dir = os.listdir(work_dir)
            list_work_dir = sorted_nicely(list_work_dir) 
            for file in list_work_dir: 
                if (file == 'all_edge_values.csv'):
                    changeEdgeNames = 'python3 ' + plotting_dir + elem + ' ' + work_dir + file + ' ' + gml_file + ' -a -o ' + work_dir + 'all_edges_values_biological_names.csv'  
                    log(changeEdgeNames, 'd') 
                    os.chdir(out_dir) 
                    os.system(changeEdgeNames)
                    os.chdir(work_dir)
                elif (file.endswith('edges_weights.csv')):
                    changeEdgeNames = 'python3 ' + plotting_dir + elem + ' ' + work_dir + file + ' ' + gml_file + ' -e -o ' + work_dir + 'edges_weights_biological_names.csv'
                    log(changeEdgeNames, 'd')
                    os.chdir(out_dir) 
                    os.system(changeEdgeNames)
                    os.chdir(work_dir)
        
        elif (add_changeEdgeNames_args !=''):
            work_dir = get_working_dir(evalEdgesWeights_dir)
            changeEdgeNames = 'python3 ' + plotting_dir + elem + ' ' + add_changeEdgeNames_args
            os.chdir(out_dir) 
            os.system(changeEdgeNames)
            os.chdir(work_dir)
            
        changeEdgeNames_dir = os.path.abspath(out_dir) + '/'

        log('changeEdgeNames computations are done.', 'i')
        
    elif (elem == 'sumEdgeWeights.py'):
        work_dir = get_working_dir(changeEdgeNames_dir)
        log(work_dir, 'd')
        
        if (add_sumEdgeWeights_args == ''):  
            list_work_dir = os.listdir(work_dir)
            list_work_dir = sorted_nicely(list_work_dir)     
            for file in list_work_dir:      
                if ('all_edges_values_biological_names.csv' in file):
                    sumEdgeWeights = 'python3 ' + plotting_dir + elem + ' ' + work_dir + file + ' -n ' + out_dir + 'edge_sum.pdf'
                    log(sumEdgeWeights, 'd')
                    os.chdir(out_dir) 
                    os.system(sumEdgeWeights)
                    os.chdir(work_dir)
                
        elif (add_sumEdgeWeights_args != ''):
            sumEdgeWeights = 'python3 ' + plotting_dir + elem + ' ' + add_sumEdgeWeights_args
            log(sumEdgeWeights, 'd')
            os.chdir(out_dir) 
            os.system(sumEdgeWeights)
            os.chdir(work_dir)
                        
        log('sumEdgeWeights computations are done.', 'i')


    elif (elem == 'plotSnapshots.py'):
        work_dir = get_working_dir(changeEdgeNames_dir)
        log(work_dir, 'd')

        if (add_plotSnapshots_args == ''):
            list_work_dir = os.listdir(work_dir)
            list_work_dir = sorted_nicely(list_work_dir)     
            for file in list_work_dir:      
                if ('all_edges_values_biological_names.csv' in file):
                    plotSnapshots = 'python3 ' + plotting_dir + elem + ' ' + work_dir + file + ' -n ' + out_dir + 'edge_plots.pdf'
                    log(plotSnapshots, 'd')
                    os.chdir(out_dir) 
                    os.system(plotSnapshots)
                    os.chdir(work_dir)
                    
        elif (add_plotSnapshots_args != ''):
            plotSnapshots = 'python3 ' + plotting_dir + elem + ' ' + add_plotSnapshots_args
            log(plotSnapshots, 'd')
            os.chdir(out_dir) 
            os.system(plotSnapshots)
            os.chdir(work_dir)

        log('plotSnapshots computations are done.', 'i')
        

    elif (elem == 'PyMolHeatmapVisualisation.py'):
        # create csv file with number of residues in each chain
        files_dir = get_working_dir(file_dir)
        list_file_dir = os.listdir(files_dir)
        list_file_dir = sorted_nicely(list_file_dir)
        
        os.chdir(out_dir)
        input_dir_csv_file = new_directory('files_for_GraCom_computation') + '/' 
        os.chdir(files_dir)
        
        cif_or_pdb_file = ''
        for file in list_file_dir:
            if file.endswith('.cif'):
                shutil.copy(file, input_dir_csv_file + file)
                cif_or_pdb_file = input_dir_csv_file + os.path.basename(file)
                break
                                          
        if (cif_or_pdb_file == ''):
            list_input_dir = os.listdir(input_dir)
            list_input_dir = sorted_nicely(list_file_dir)
            for file in list_input_dir:
                if file.endswith('.pdb'):
                    shutil.copy(file, input_dir_csv_file + file)
                    cif_or_pdb_file = input_dir_csv_file + os.path.basename(file)
                    break
        
        # get corresponding dssp file
        base = os.path.basename(cif_or_pdb_file)
        name = os.path.splitext(base)[0]
        dssp_path = dssp_dir + name + '.dssp'
        dssp_file = ''
        if os.path.exists(dssp_path):
            dssp_file = dssp_path
            shutil.copy(dssp_file, input_dir_csv_file + name + '.dssp')
            dssp_file = input_dir_csv_file + name + '.dssp'
        else:
            log("There exists no dssp file " + dssp_path + ". Can´t create csv file with number of residues in each chain.", 'e')   
           
        if dssp_file != '' and cmd_header_mmcif != '': 
            cif_file = cif_or_pdb_file       
            if cif_or_pdb_file.endswith('.pdb'):
                createMmCifFile = cmd_start + 'toMmCIF.py -i ' + os.path.dirname(cif_or_pdb_file) + ' -p ' + out_dir + cmd_header_mmcif
                os.chdir(out_dir)
                os.system(createMmCifFile)
                cif_file = os.path.abspath(out_dir) + '/' + name + '.cif'
                os.chdir(files_dir)
            
            createCsvFile = 'java -jar ' + default_path_graCom + ' none -p ' + cif_file + ' -d ' + dssp_file + ' -o ' + input_dir_csv_file + ' -I -N --silent --dont-write-images --set "PTGLgraphComputation_B_csv_number_residues_chains" "true"' 
            log("Csv file with the number of residues in each chain created.", 'i')
            os.chdir(out_dir)
            os.system(createCsvFile)  
            os.chdir(files_dir) 
            
        else:
            log("There was no header file given for a mmcif file or no dssp file found. As PTGLgraphComputation will not compute output without, there can not be created the csv file containing the number of residues for each chain.", 'e')
                
        # Create PyMOL script
        if (add_pyMolHeatmapVisualisation_args == ''):            
            work_dir = get_working_dir(compareSubsets_dir)
            log(work_dir, 'd')
            
            createPymolScript = 'python3 ' + plotting_dir + elem + ' ' + work_dir + ' ' + cif_or_pdb_file + ' -p ' + out_dir 
            log(createPymolScript, 'd')
            os.chdir(out_dir) 
            os.system(createPymolScript)
            os.chdir(work_dir)
                    
        elif (add_pyMolHeatmapVisualisation_args != ''):
            createPymolScript = 'python3 ' + plotting_dir + elem + ' ' + add_pyMolHeatmapVisualisation_args
            log(createPymolScript, 'd')
            os.chdir(out_dir) 
            os.system(createPymolScript)
            os.chdir(work_dir)

        log('PyMolHeatmapVisualisation computations are done.', 'i')
            
        
           

log("-- %s seconds ---"% (time.time()- _start_time), 'i')
log("All done, exiting ptglDynamics.", 'i')

# tidy up
if (args.outputfile != ""):
    output_file.close()

