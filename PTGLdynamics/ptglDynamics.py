#################################################
#################### HOW TO #####################
############### USE THIS TEMPLATE ###############
"""
- argparse
    - describe what the script does (see TODO)
    - add command line arguments with it
    - argument file can be specified with @filepath
- logging
    - use the function 'log' for each and every print!
        -> d(ebug) for debug messages
        -> i(nfo) for status and meta information
        -> w(arning), e(rror), c(ritical) for more severe messages
        -> without a level to print results and if required obligatory non-result prints
            -> command line option 'silent' should suppress those obligatory prints so that only results are printed
- pipelining
    - write the results with -o to an output file to use as next input or
    - use 'silent' mode (if obligatory prints exist) and pipe the output of stdout directly to next script
        -> all debug, info, warning etc. are printed to stderr
- imports from not built-in modules in the respective section with error handling and useful prints where to find the module and if possible how to install it
- change version to your script's version (see TODO)
        
(Remove this HOW TO and have fun programming!)
"""


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
    if not os.path.isdir(name):
        os.makedirs(name)
    return os.path.realpath(name)


########### configure logger ###########

logging.basicConfig(format = "[%(levelname)s] %(message)s")
            
            
########### not built-in imports ###########

# EXAMPLE - can be removed
"""
    from pygmlparser.Parser import Parser
except ModuleNotFoundError as exception:
    log(exception,"e")
    traceback.print_exc()
    log("  Module from https://github.com/hasii2011/PyGMLParser seems to be missing.", "e")
    log("  Try installing it with 'pip3 install PyGMLParser'. Exiting now.", "e")
    sys.exit(1)
"""

            
########### command line parser ###########

## create the parser
cl_parser = argparse.ArgumentParser(description="This is a pipeline that reads several snapshots in pdb-format. The pipeline runs plcc on each of them and compares the differences in between their Complex Graph outputs. If not specified otherwise the results are written in the current folder in csv-format.",
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

cl_parser.add_argument('-id',
                       '--inputdirectory',
                       default = '',
                       help = 'specify a path to your input files. Otherwise the current folder is used.')

cl_parser.add_argument('-od',
                       '--outputdirectory',
                       default = '',
                       help = 'specify a path to your output files. Otherwise the current folder is used.')

cl_parser.add_argument('-c',
                       '--compoundfile',
                       default = '',
                       help = 'to integrate a header in pdb files specify the path of your compound file.')

cl_parser.add_argument('-p',
                       '--programms',
                       nargs = "*",
                       type = str,
                       default = ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'plcc', 'gmlCompareEdgeWeightsAndSubsets.py'],
                       help = 'to execute only the specified scripts.')

cl_parser.add_argument('-u',
                       '--subdirectorystructure',
                       action='store_true',
                       help='display the results in sub-directories in the output directory.')

cl_parser.add_argument('-a',
                       '--plcc_args',
                       type = str,
                       default = '',
                       help = 'a string with the plcc arguments you want to use and its values to execute plcc in different ways using plccs command line arguments. Insert arguments like this: -a="<arguments and their inputs>" ')

cl_parser.add_argument('-b',
                       '--toLegacyPDB_args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for toLegacyPDB you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -b="<arguments and their inputs>" ')

cl_parser.add_argument('-e',
                       '--postProcessDssp_args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for postProcessDssp you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -e="<arguments and their inputs>" ')

cl_parser.add_argument('-f',
                       '--gmlCompareEdgeWeightsAndSubsets_args',
                       type = str,
                       default = '',
                       help = 'a string with the arguments for gmlCompareEdgeWeightsAndSubsets you want to use and its values to execute the script in different ways using your command line arguments. Insert arguments like this: -f="<arguments and their inputs>" ')


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
if (args.inputdirectory != ""):
    if(os.path.isdir(args.inputdirectory)):
        i_dir = args.inputdirectory
    else:
        logging.error("Specified input directory '%s' is not readable. Exiting now.", args.inputdirectory)
        sys.exit(1)
else:
    i_dir = os.getcwd() + '/'
    

# output directory
if (args.outputdirectory != ""):
    if(os.path.isdir(args.outputdirectory)):
        o_dir = args.outputdirectory
    else:
        logging.error("Specified output directory '%s' is not writable. Exiting now.", args.outputdirectory)
        sys.exit(1)
else:
    o_dir = os.getcwd() + '/'
    
    
# output directory
if (args.compoundfile != ""):
    if(os.access(args.compoundfile, os.R_OK)):
        compound = args.compoundfile
        cmd_compound = ' -c ' + compound
    else:
        logging.error("Specified compound file '%s' is not readable. Continuing without compound file.", args.compoundfile)
        compound = ''
        cmd_compound = ''
else:
    compound = ''
    cmd_compound = ''

# list of programms
programm_list = []
if (args.programms != []):
    for programm in args.programms:
        if programm in ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'plcc', 'gmlCompareEdgeWeightsAndSubsets.py']:
            programm_list.append(programm)
        else:
            logging.error("Specified programm '%s' is not part of the ptglDynamics pipeline. Continuing without it.", programm)
else:
    programm_list = ['toLegacyPDB.py', 'dsspcmbi', 'postProcessDssp.py', 'plcc', 'gmlCompareEdgeWeightsAndSubsets.py']
    
# plcc arguments
if (args.plcc_args != ''):
    add_plcc_args = args.plcc_args
else:
    add_plcc_args = ''

# toLegacyPDB arguments
if (args.toLegacyPDB_args != []):
    add_toLegacyPDB_args = args.toLegacyPDB_args
else:
    add_toLegacyPDB_args = ''

# postProcessDssp arguments
if (args.postProcessDssp_args != ''):
    add_postProcessDssp_args = args.postProcessDssp_args
else:
    add_postProcessDssp_args = ''

# gmlCompareEdgeWeightsAndSubsets arguments
if (args.gmlCompareEdgeWeightsAndSubsets_args != ''):
    add_gmlComparison_args = args.gmlCompareEdgeWeightsAndSubsets_args
else:
    add_gmlComparison_args = ''

########### vamos ###########

log("Version " + version, "i")

# TODO add your code here



#general variables
dir_names = {'toLegacyPDB.py':'legacyPDB', 'dsspcmbi':'dssp', 'postProcessDssp.py':'dssp', 'plcc':'plcc', 'gmlCompareEdgeWeightsAndSubsets.py': 'gml'}
curr_path = os.path.dirname(__file__)
curr_wd = os.getcwd()
plccJarPath = os.path.dirname(curr_path) + '/plcc/dist/plcc.jar'

#Variables of Command Line Arguments
cmd_start = 'python3 ' + curr_path + '/codes/'


pdb_o_dir = i_dir
dssp_o_dir = i_dir
#plcc_o_dir = i_dir
plcc_o_dir = curr_wd
gml_dir = i_dir
log("i_dir: " + i_dir, 'd')
log("programm list: ", 'd')
log(programm_list, 'd')

for elem in programm_list:
    log("elem: " + elem, 'd')
    if (args.subdirectorystructure) and (dir_names[elem] != ''):
        out_path = new_directory(o_dir + dir_names[elem]) + '/'
    else:
        out_path = o_dir
        
    #execute different scripts:
    if (elem == 'toLegacyPDB.py'):
        #exec_string = cmd_start + elem + ' -id ' + i_dir + ' -od ' + out_path + cmd_compound + ' -d &> toLegacyPDB_out.txt'
        exec_string = cmd_start + elem + ' ' + add_toLegacyPDB_args + ' -id ' + i_dir + ' -od ' + out_path + cmd_compound
        log('exec_string ' + exec_string, 'd')
        os.system(exec_string)
        pdb_o_dir = out_path
        
      
    elif (elem == 'dsspcmbi'):
        
        for pdb in os.listdir(pdb_o_dir):
            if pdb.endswith(".pdb"):
                dssp = out_path + pathlib.Path(pdb).stem + '.dssp'
                exec_string = curr_path + '/codes/' + elem + ' ' + pdb_o_dir + pdb + ' ' + dssp
                log('exec_string ' + exec_string, 'd')
                os.system(exec_string)
        dssp_o_dir = out_path
    
    elif (elem == 'postProcessDssp.py'):
        #exec_string = cmd_start + elem + ' -id ' + dssp_o_dir + ' -od ' + out_path + ' -d &> post_process_out.txt'
        exec_string = cmd_start + elem + ' ' + add_postProcessDssp_args + ' -id ' + dssp_o_dir + ' -od ' + out_path
        log('exec_string ' + exec_string, 'd')
        os.system(exec_string)
        dssp_o_dir = out_path
       
    
    elif (elem == 'plcc'):
        for pdb in os.listdir(pdb_o_dir):
            if (pdb.endswith(".pdb")):
                pdb_id = pathlib.Path(pdb).stem
                pdb_id_folder = ''
                if(args.subdirectorystructure):
                    pdb_id_folder = new_directory(out_path + pdb_id)
                else:
                    pdb_id_folder = out_path
                
                pdb = pdb_o_dir + pdb
                dssp = dssp_o_dir + pathlib.Path(pdb).stem + '.dssp'
                
                #plcc = 'java -jar ' + plccJarPath + ' ' + pdb_id + ' -p ' + pdb + ' -d ' + dssp + ' -o ' + pdb_id_folder + ' -O g' 
                
                plcc = 'java -jar ' + plccJarPath + ' ' + pdb_id + ' ' + add_plcc_args + ' -p ' + pdb + ' -d ' + dssp + ' -o ' + pdb_id_folder
                log(plcc,'d')        
                os.system(plcc)
        
        #plcc_o_dir = out_path
        
        
    elif (elem == 'gmlCompareEdgeWeightsAndSubsets.py'):
        prevGml = ''
        for entry in os.listdir(plcc_o_dir):
            if (os.path.isdir(entry)):
                for gml in os.listdir(entry):
                    if gml.endswith("complex_chains_albelig_CG.gml"):
                        try:
                            shutil.copy(gml, out_path + gml)
                        except shutil.SameFileError:
                            log("Source and destination represents the same file.", 'i')
                            pass

                        if prevGml != '':
                            gmlComparison = cmd_start + elem +' ' + add_gmlComparison_args + ' -i1 ' + out_path + prevGml + ' -i2 ' + out_path + gml + ' -od ' + o_dir
                            log(gmlComparison,'d')
                            os.system(gmlComparison)
                        #log(prevGml +' '+gml,'i')
                        prevGml = gml
            else:
                if entry.endswith("complex_chains_albelig_CG.gml"):
                    try:
                        shutil.copy(entry, out_path + entry)
                    except shutil.SameFileError:
                        log("Source and destination represents the same file.", 'i')
                        pass

                    if prevGml != '':
                        gmlComparison = cmd_start + elem + ' ' + add_gmlComparison_args + ' -i1 ' + out_path + prevGml + ' -i2 ' + out_path + entry + ' -od ' + o_dir
                        log(gmlComparison,'d')
                        os.system(gmlComparison)
                    #log(prevGml +' '+ entry, 'i')
                    prevGml = entry
    
        gml_dir = out_path



# tidy up
if (args.outputfile != ""):
    output_file.close()
