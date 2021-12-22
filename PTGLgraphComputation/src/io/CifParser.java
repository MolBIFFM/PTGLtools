/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;

// imports
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import settings.Settings;
import proteinstructure.AminoAcid;
import proteinstructure.Atom;
import proteinstructure.Chain;
import proteinstructure.Model;
import proteinstructure.Molecule;
import proteinstructure.ProtMetaInfo;
import proteinstructure.Residue;
import proteinstructure.RNA; 
import proteinstructure.Ligand;
import tools.DP;


/**
 *
 * @author niclas
 */
class CifParser {
    
    // declare class vars
    // control vars
    static boolean dataInitDone;
    static boolean silent;
    
    // data structures
    protected static HashMap<String, String> metaData;
    private static ArrayList<ProtMetaInfo> allProteinMetaInfos = new ArrayList<>();
    private static int lastIndexProtMetaInfos = 0;  // used to traverse allProteinMetaInfos quicker for sucessive request, e.g., in the same order they were added
    private static String pdbID;
    private static HashMap<String, HashMap<String, String>> entityInformation = new HashMap<>();  // <entity ID, <column head, data>>
    protected static HashMap<String, String> chainIdentity = new HashMap<>();       // matches chain ID with its molecule type
    protected static HashMap<String, HashMap<String, String>> chemicalComponents = new HashMap<>();   // stores all information on chemical components
    
    // - - - vars for parsing - - -
    private static Boolean dataBlockFound = false;  // for now only parse the first data block (stop if seeing 2nd block)
    private static Integer numLine = 0;
    private static Boolean inLoop = false;
    private static Boolean inString = false;
    private static String[] lineData;  // array which holds the data items from one line (seperated by spaces)
    private static String interruptedLine = "";  // used to build up a line that is interrupted by a multi-line string
    
    // - - variables for one loop (reset when hitting new loop) - -
    private static String currentCategory = null;
    // Key: name of column; Val: position in list
    //  -> if auth columns from atom_site not present they will be mapped to the PDB columns
    //     therefore always use auth columns unless you explicitly want the PDB ones
    private static HashMap<String,Integer> colHeaderPosMap = new HashMap<>();
    private static Boolean columnsChecked = false;
    private static ArrayList<String> currentLineValues = new ArrayList<String>();    //stores values of non-loop-blocks so they can be parsed as "fake-loops"
    
    // - - atom_site - -
    private static int ligandsTreatedNum = 0;
    private static int RnaTreatedNum = 0;
    private static int freeResTreatedNum = 0;
    private static int numberAtoms = 0;
    
    // - variables for successive matching atom -> residue/RNA : Molecule -> chain -
    private static Model m = null;
    private static Molecule lastMol = null;    // starts as first residue and is always the actual one
    private static Molecule tmpMol = null;      // used to save lastMol if getResidue returns null
    private static Chain tmpChain = null;
    private static Residue res = null;
    private static Ligand lig = null;
    private static RNA rna = null;

    // - variables per (atom) line -
    private static Integer atomSerialNumber, coordX, coordY, coordZ, molNumPDB, entityID;
    private static String atomRecordName, atomName, chainID, altChainID, chemSym, altLoc, iCode, molNamePDB;
    private static Double oCoordX, oCoordY, oCoordZ;            // the original coordinates in Angstroem (coordX are 10th part Angstroem)
    private static Float oCoordXf, oCoordYf, oCoordZf;
    private static int lastLigandNumPDB = 0; // used to determine if atom belongs to new ligand residue
    private static String lastChainID = ""; // s.a.
    private static int lastRnaNumPDB = 0;
    private static String[] tmpLineData;
    private static String tmpModelID;
    private static String chainNum = null;  // identifier of the current chain e.g. A
    private static String chainType = null; // type of the current chain e.g. RNA
    private static String nameOrgCommon = "";
    private static String nameOrgScientific = "";
    private static String nameOrgCommonSource = "";  // line that an orgName comes from
    private static String nameOrgScientificSource = "";  // line that an orgName comes from

    // - variables for already printed warnings -
    private static Boolean furtherModelWarningPrinted = false;
    
    // - - - - - - 
    
    // constants
    static final String[] NO_VALUE_PLACEHOLDERS = {".", "?"};  // characters that are used in mmCIFs as placeholder when no value exists
    static final String SINGLE_LINE_STRING_MARKER = "'";
    
    
    /**
     * Calls hidden FileParser method initVariables and inits additional CIF Parser variables.
     * @param pf PDB file path
     */
    private static void initVariables(String pf) {
        metaData = new HashMap<>();
    }
    
    /**
     * Like initData but for mmCIF data.
     * @param pf Path to a PDB file. Does NOT test whether it exist, do that earlier.
     * @param df Path to a DSSP file. Does NOT test whether it exist, do that earlier.
     * @return 
     */
    protected static Boolean initData(String pf) {
        
        initVariables(pf);
        
        silent = FileParser.settingSilent();
        
        if(parseData()) {
            dataInitDone = true;
            return(true);
        }
        else {
            System.err.println("ERROR: Could not parse dssp and pdb data.");
            dataInitDone = false;
            System.exit(1);
            return(false);          // for the IDE ;)
        }
        
    }
    
    
    /**
     * Goes through all lines of mmCIF PDB and DSSP file and applies appropriate function to handle each line.
     * Goes through mmCIF file exactly once, saves only what is needed and matches atom<->molecule<->chain<->model on the fly.
     * @return ignore (?)
     */
    private static Boolean parseData() {
        
        createResidues();
        lastMol = new Residue(); // create artificial molecule to fill so there is no NullPointerException, it will be overwritten once atoms are parsed
        
        try {
            BufferedReader in = new BufferedReader(new FileReader(FileParser.pdbFile));
            String line;
            while ((line = in.readLine()) != null) {
                numLine ++;
                                                             
                //
                // - - control sequences - -
                //
                // check for beginning of loop (table-like lines)
                if (line.startsWith("loop_")) {
                    // loops must not be nested by mmCIF definition
                    if (inLoop) {
                        DP.getInstance().e("FP_CIF", "Found a nested loop starting in line " + numLine + " which is forbidden by the mmCIF definition. Exiting now.");
                        System.exit(1);
                    }
                    // from now on: inLoop was false
                    
                    inLoop = true;
                    continue;  // nothing else to parse here
                }
                
                // check if line is a comment ending a loop
                if (line.startsWith("#")) {
                    inLoop = false;
                    colHeaderPosMap.clear();
                    columnsChecked = false;
                    currentCategory = null;
                    inString = false;
                    interruptedLine = "";
                    currentLineValues.clear();
                    continue;  // nothing else to do / parse here
                }
                // from now on: line does not start with '#' (is no comment)
                
                // - - handle multi-line strings - -
                // check if line encloses a string
                if (line.startsWith(";")) {
                    inString = !inString;
                    
                    if (!inString) {
                        interruptedLine += "\n;";  // add newline char to mark end of multi-line string & parse combined line
                        line = "";  // so we dont add it multiple times and with space
                    }
                }
                
                // check if line is part of a multi-line string
                if (inString) {
                        interruptedLine += line;  // keep the semicolon to hide special characters inside. reset is done after switch/case treating a (combined) line
                    continue;  // only parse combined line, nothing else to parse here
                }
                // from now on: cannot be inString and line should always hold all data items (not splitted), despite when new line and single quotation string
                
                //
                // - - parse data - -
                //
                
                // 1) check for data block - do this FIRST b/c line is too short for a 'normal' line
                if (line.startsWith("data_")) {
                    if (! handleDataLine(line)) {
                        // encountered second data block -> stop here
                        DP.getInstance().w("FP_CIF", " Parsing of first data block ended at line " + numLine.toString()
                            + " as right now only the first data block is parsed.");
                        break;
                    }
                    continue;  // nothing else to parse here
                }
                // from now on: not in data block line
                
                // 2) check for category - do this SECOND b/c column header are too short for a 'normal' line
                if (line.startsWith("_")) {
                    currentCategory = line.split("\\.")[0];
                    if (inLoop) {
                        addColumnHeaderToMap(line);
                        continue;  // nothing else to parse here
                    }
                    else {
                        String[] elementsSeperatedByDotAndTab = line.split("\\.| ");
                        colHeaderPosMap.put(elementsSeperatedByDotAndTab[1].trim(), colHeaderPosMap.size());  // add header to colHeaderPosMap
                        String value = (line.substring(line.indexOf(' ')+1)).trim();
                        currentLineValues.add(value);  // save value in currentLineValues
                    }
                }
                // from now on: not in loop column header definition
                
                // get the data items
                if (! interruptedLine.equals("")) {
                    // we want to combine the line with previous line(s)
                    if (inLoop){
                        line = interruptedLine + " " + line;
                    }
                    else {
                        // since a value is always added to currentLineValues (even if it is empty) when filling the colHeaderPosMap, 
                        // the empty value has to be deleted to be replaced by the interruptedLine
                        currentLineValues.remove(currentLineValues.size() -1);
                        currentLineValues.add(interruptedLine);
                    }
                    interruptedLine = "";
                }
                
                if (inLoop) {
                    lineData = lineToArray(line);
                } else {
                    lineData = trimSpecialChars(currentLineValues);
                }
                
                // check for minimum length of lineData (2 for non-loop and #header for loop) and merge with coming line if not
                if ((inLoop && lineData.length < colHeaderPosMap.keySet().size()) || (!inLoop && lineData.length < 1)) {
                        interruptedLine += line;
                        continue;
                }
                
                if (lineData.length > colHeaderPosMap.keySet().size()) {
                    DP.getInstance().w("FP_CIF", "Line " + numLine + " (together with previous if combined) seems to be too long. Trying to ignore it, but may miss data.");
                }
                // from now on: line contains expected lenght of data items

                switch(currentCategory) {
                case "_exptl":
                    // check for experimental method
                    handleExptlLine();
                    break;
                case "_reflns":
                    handleReflnsLine();
                    break;
                case "_refine":
                    // check for resolution (meta data)
                    // TODO: really the best attribute(s) to do this?
                    handleResolutionLine();
                    break;
                case "_entity":
                    // parse entity information to retrieve mol name later
                    handleEntityLine();
                    break;
                case "_entity_poly":
                    // check for homologues chains
                    handleEntityPolyLine();
                    break;
                case "_atom_site":
                    // check for atom coordinate data
                    handleAtomSiteLine();
                    break;
                case "_chem_comp":
                    // check for information on chemical components
                    handleChemComp();
                    break;
                case "_entity_src_gen":
                    handleEntitySrcGen();
                    break;
                case "_entity_src_nat":
                    handleEntitySrcNat();
                    break;
                case "_pdbx_entity_src_syn":
                    handlePdbxEntitySrcSyn();
                    break;
                case "_struct":
                    handleStructLine();
                    break;
                case "_struct_keywords":
                    handleStructKeywords();
                    break;
                case "_pdbx_database_status":
                    handlePdbxDatabaseStatus();
                    break;
                }
                
                // reset here, b/c we only get here when a (combined) line was treated
                interruptedLine = "";


            } // end of reading in lines           
	} catch (IOException e) {
            System.err.println("ERROR: Could not parse PDB file.");
            System.err.println("ERROR: Message: " + e.getMessage());
            System.exit(1);
	}
        
        if (! (silent || FileParser.essentialOutputOnly)) {
            System.out.println("  PDB: Found in total " + FileParser.s_chains.size() + " chains.");
        }
        
        // alt loc treatment copy&pasted from old parser
        if(! (FileParser.silent || FileParser.essentialOutputOnly)) {
            System.out.println("    PDB: Hit end of PDB file at line " + numLine + ".");

            // remove duplicate atoms from altLoc here
            System.out.println("    PDB: Selecting alternative locations for atoms of all residues.");
        }
        ArrayList<Atom> deletedAtoms;
        int numAtomsDeletedAltLoc = 0;
        int numResiduesAffected = 0;
        Molecule m;
        for(int i = 0; i < FileParser.s_molecules.size(); i++) {
            m = FileParser.s_molecules.get(i);
            deletedAtoms = m.chooseYourAltLoc();


            if(deletedAtoms.size() > 0) {
                numResiduesAffected++;
            }

            //delete atoms from global atom list as well
            for(Atom a : deletedAtoms) {
                if(FileParser.s_atoms.remove(a)) {
                    numAtomsDeletedAltLoc++;
                } else {
                    DP.getInstance().w("Atom requested to be removed from global list does not exist in there.");
                }
            }
        }
        
        // create all protein meta data from entityInformation
        fillProteinMetaData();
        
        // add empty Strings to metadata to avoid SQL null errors
        fillMetadataEmptyStrings();
        if (FileParser.s_chains.size() > 62 || numberAtoms > 99999) {
            metaData.put("isLarge", "true");
        } else {
            metaData.put("isLarge", "false");
        }
        
        // all lines have been read
        return(true);
    }
    
    
    /**
     * Calls DSSP parser to create all residues from DSSP file.
     */
    private static void createResidues() {
        if(! silent) {
            System.out.println("  Creating all Molecules...");
        }

        DsspParser.createAllResiduesFromDsspData(true);

        // If there is no data part at all in the DSSP file, the function readDsspToData() will catch
        //  this error and exit, this code will never be reached in that case.
        if(FileParser.s_molecules.size() < 1) {
            DP.getInstance().e("FP_CIF", "DSSP file contains no residues (maybe the PDB file only holds DNA/RNA data). Exiting.");
            System.exit(2);
        }
    }
    
    
    /**
     * Handles a line starting with 'data' defining a data block. Parses the PDB ID.
     * @return if this is the first data block encountered
     */
    private static boolean handleDataLine(String line) {
        if (dataBlockFound) {
            return false;
        }
        else {
            dataBlockFound = true;
            if (line.length() > 5) {
                pdbID = line.substring(5, line.length()).toLowerCase();
                if (! silent) {
                    System.out.println("  PDB: Found the first data block named: " + pdbID);
                }
            }
            else {
                if (! silent) {
                    DP.getInstance().w("FP_CIF", "Expected first data block to be named after PDB ID, but found no name. "
                            + "Protein meta information for the chains will contain to PDB ID because of this.");
                    pdbID = "";
                }
            }
            return true;
        }
    }
    
    
    /**
     * Handles a line starting with '_exptl.' defining the experimental method.
     * Saves the experiment method in metaData.
     */
    private static void handleExptlLine() {
        if (colHeaderPosMap.get("method") != null ){
            if (! lineData[colHeaderPosMap.get("method")].equals("?")){
                metaData.put("experiment", lineData[colHeaderPosMap.get("method")]);
            }
        }
        else{
            metaData.put("experiment", "");
        }
    }
    
    
    /**
     * Handles lines starting with '_struct.'.
     * Saves the title of the PDB File in metaData.
     */
    private static void handleStructLine(){
        if (colHeaderPosMap.get("title") != null ){
            if (! lineData[colHeaderPosMap.get("title")].equals("?")){
                metaData.put("title", lineData[colHeaderPosMap.get("title")]);
            }
        }
        else{
            metaData.put("title", "");
        }
    }
    
    
    /**
     * Handles lines starting with '_struct_keywords'.
     * Saves possible keywords in metaData.
     */
    private static void handleStructKeywords(){
        if (colHeaderPosMap.get("text") != null ){
            if (! lineData[colHeaderPosMap.get("text")].equals("?")){
                metaData.put("keywords", lineData[colHeaderPosMap.get("text")]);
            }
        }
        else{
            metaData.put("keywords", "");
        }
        if (colHeaderPosMap.get("pdbx_keywords") != null ){
            if (! lineData[colHeaderPosMap.get("pdbx_keywords")].equals("?")){
                metaData.put("header", lineData[colHeaderPosMap.get("pdbx_keywords")]);
            }
        }
        else{
            metaData.put("header", "");
        }
    }
    
    
    /**
     * Handles lines starting with '_pdbx_database_status'.
     * Saves the date in metaData.
     */
    private static void handlePdbxDatabaseStatus(){
        if (colHeaderPosMap.get("recvd_initial_deposition_date") != null ){
            if (! lineData[colHeaderPosMap.get("recvd_initial_deposition_date")].equals("?")){
                metaData.put("date", lineData[colHeaderPosMap.get("recvd_initial_deposition_date")]);
            }
        }
        else{
            metaData.put("date", "");
        }
    }
    
    
    /**
     * Handles a line starting with '_refine.'.
     * Defining the resolution. Saves the value in metaData. Is probably not the best way to extract the resolution.
     * Resolution can be parsed different ways, this is the preferred one.
     */
    private static void handleResolutionLine() {
        if (colHeaderPosMap.get("ls_d_res_high") != null ){
            if (! lineData[colHeaderPosMap.get("ls_d_res_high")].equals("?")){
                metaData.put("resolution", lineData[colHeaderPosMap.get("ls_d_res_high")]);
            }
        }
        else{
            if (! metaData.containsKey("resolution")) {
                metaData.put("resolution", "");
            }
        }
    }
    
    /**
     * Handles a line starting with '_reflns.'.
     * Alternative way to parse resolution.
     */
    private static void handleReflnsLine() {
        if (! metaData.containsKey("resolution") && colHeaderPosMap.get("d_resolution_high") != null) {
            if (! lineData[colHeaderPosMap.get("d_resolution_high")].equals("?")){
                metaData.put("resolution", lineData[colHeaderPosMap.get("d_resolution_high")]);
            }
        }
        else{
            if (! metaData.containsKey("resolution")){
                metaData.put("resolution", "");
            }
        }
    }
    
    
    /**
     * Handles a line starting with '_entity.'.
     * Parses information on which entities (chains/ligands) exist in the molecule. Information is saved in entityInformation map.
     */
    private static void handleEntityLine() {
        String tmpValue;

        String tmpEntityID = lineData[0];
        entityInformation.put(tmpEntityID, new HashMap<String, String>());

        // get all available information per entity, despite ID (saved as superior key)
        for (String colHeader : colHeaderPosMap.keySet()) {
            if (! colHeader.equals("id")) {
                tmpValue = (valueIsAssigned(lineData[colHeaderPosMap.get(colHeader)]) ? lineData[colHeaderPosMap.get(colHeader)] : null);  // assign value or null
                entityInformation.get(tmpEntityID).put(colHeader, (valueIsAssigned(tmpValue) ? tmpValue : ""));
            }

        }            
    }
    
    
    /**
     * Handle a line starting with '_entity_poly.' holding entity information of polymers.
     * Fills, for example, the homologuesMap and the chainIdentity Map which maps chain ID to its molecule type.
     */
    private static void handleEntityPolyLine() {
            if (colHeaderPosMap.get("pdbx_strand_id") != null ) {
                FileParser.fillHomologuesMapFromChainIdList(lineData[colHeaderPosMap.get("pdbx_strand_id")].split(","));
                
                // Multiple chains can be listed under "pdbx_strand_id" so they have to be separated and individually added to chainIdentity
                String[] chainList = lineData[colHeaderPosMap.get("pdbx_strand_id")].split(",");
                for (String s : chainList) {
                    chainIdentity.put(s, lineData[colHeaderPosMap.get("type")]);
                }
            }
        }
    
    
    /**
     * Handle a line starting with '_atom_site.' holding atom coordinates. Creates the atoms, ligands, chains and models.
     * Matches atoms <-> residues/rna/ligands <-> chains <-> models.
     */
    private static void handleAtomSiteLine() {
        
        Chain lastChain = new Chain(" ");
        
        // atom coordinates should always be within a loop      
        if (! inLoop) {
            DP.getInstance().e("FP_CIF", "Parsing line " + numLine + ". Atom coordinates seem not be within a loop. Is the file broken? Exiting now.");
            System.exit(2);
        }

        // we are in the row section (data!)

        numberAtoms++;

        // check once if required column headers are present
        if (! columnsChecked) {
            ArrayList<String> missingCols = checkColumns(currentCategory, new ArrayList<>(colHeaderPosMap.keySet()));
            if (missingCols.size() > 0) {
                DP.getInstance().e("FP_CIF", "Missing following columns in " + currentCategory + 
                        ": " + missingCols);
                DP.getInstance().e("FP_CIF", " Exiting now.");
                System.exit(1);
            }

            // if auth columns not present map them to PDB ones
            // pdbx_PDB_model_num not checked as no equivalent existing (just use default model 1 if not existing)
            String[] authCols = {"auth_atom_id", "auth_asym_id", "auth_comp_id", "auth_seq_id"};
            // Matching equivalents to author columns
            String[] pdbCols = {"label_atom_id", "label_asym_id", "label_comp_id", "label_seq_id"};
            for (int i = 0; i < authCols.length; i++) {
                if (colHeaderPosMap.get(authCols[i]) == null) {
                    colHeaderPosMap.put(authCols[i], colHeaderPosMap.get(pdbCols[i]));
                    if (! silent) {
                        System.out.println("   Using " + pdbCols[i] + " instead of "+ 
                                "missing column " + authCols[i]);
                    }
                }
            }

            columnsChecked = true;
        }

        // - - model - -
        // Look if model numbers are included
        if (colHeaderPosMap.get("pdbx_PDB_model_num") != null) {
            tmpModelID = lineData[colHeaderPosMap.get("pdbx_PDB_model_num")];

            // save modelID for print later
            if (! FileParser.s_allModelIDsFromWholePDBFile.contains(tmpModelID)) {
                FileParser.s_allModelIDsFromWholePDBFile.add(tmpModelID);
            }

            if (m == null) {
                // use first model
                m = new Model(tmpModelID);
                FileParser.s_models.add(m);
                if(! (FileParser.silent || FileParser.essentialOutputOnly)) {
                    System.out.println("   PDB: New model '" + m.getModelID() + "' found");
                }
            } else {
                // same model as before?
                if (! m.getModelID().equals(tmpModelID)) {
                    if (metaData.get("experiment").contains("NMR")) {
                        if (! furtherModelWarningPrinted) {
                            System.out.println("   PDB: Found further models. Ignoring them.");
                            furtherModelWarningPrinted = true;
                        }
                        // skip this line if the structure is NMR
                        return;
                    }
                }
            }
        } else {
            // create default model instead
            m = new Model("1");
            FileParser.s_models.add(m);
            System.out.println("   PDB: No model column. Creating default model '1'");
        }

        // - - chain - -
        // check for a new chain (always hold the current 
        // get chain ID
        if (colHeaderPosMap.get("auth_asym_id") != null && lineData.length >= colHeaderPosMap.get("auth_asym_id") + 1) {
                String tmp_cID = lineData[colHeaderPosMap.get("auth_asym_id")];
                
                // get macromolID
                String tmpMolId;
                if (colHeaderPosMap.get("label_entity_id") != null && lineData.length >= colHeaderPosMap.get("label_entity_id") + 1) {
                    tmpMolId = lineData[colHeaderPosMap.get("label_entity_id")];
                } else {
                    tmpMolId = "";
                }
                
                if (tmpChain == null) {
                    tmpChain = getOrCreateChain(tmp_cID, m, tmpMolId);
                } else 
                    if (! (tmpChain.getPdbChainID().equals(tmp_cID))) {
                        tmpChain = getOrCreateChain(tmp_cID, m, tmpMolId);
                    }
        }

        // - - atom - -
        // reset variables
        atomSerialNumber = molNumPDB = coordX =  coordY = coordZ = null;
        atomRecordName = atomName = molNamePDB = chainID = chemSym = altLoc = null;
        iCode = " "; // if column does not exist or ? || . is assigned use 1 blank (compare old parser)
        oCoordX = oCoordY = oCoordZ = null;            // the original coordinates in Angstroem (coordX are 10th part Angstroem)
        oCoordXf = oCoordYf = oCoordZf = null;

        // chain name
        chainID = lineData[colHeaderPosMap.get("auth_asym_id")];      // chain ID as set by author --> preferably use this one for all kinds of tasks
        altChainID = lineData[colHeaderPosMap.get("label_asym_id")];  // chain ID computed by the PDB

        // PDBx field alias atom record name
        if (colHeaderPosMap.get("group_PDB") != null) {
            if (colHeaderPosMap.get("group_PDB") < 0) {
                atomRecordName = lineData[colHeaderPosMap.get("group_PDB")];
            }
        } else {
            if( ! Settings.getBoolean("PTGLgraphComputation_B_no_parse_warn")) {
                DP.getInstance().w("FP_CIF", "Seems like _atom_site.group_PDB is missing. Trying to ignore it.");  
                colHeaderPosMap.put("group_PDB", -1);  // save that warning has been printed
            } 
        }

        // atom id alias serial number
        atomSerialNumber = Integer.valueOf(lineData[colHeaderPosMap.get("id")]); // there should be no need to trim as whitespaces should be ignored earlier

        // detailed atom name
        // old PDB files used spacing to differentiate between atoms
        // e.g. " CA " = C alpha, how to deal with this? mmCIF has no spacings
        // for now workaround for probable C alpha
        if (lineData[colHeaderPosMap.get("label_atom_id")].equals("CA")) {
            atomName = " " + lineData[colHeaderPosMap.get("label_atom_id")] + " ";
        } else {
            atomName = lineData[colHeaderPosMap.get("label_atom_id")];
        }

        // alternative location
        if (colHeaderPosMap.get("label_alt_id") != null) {
            altLoc = lineData[colHeaderPosMap.get("label_alt_id")];
        }

        // residue name or rna name 
         //resNamePDB = lineData[colHeaderPosMap.get("label_comp_id")];
        molNamePDB = lineData[colHeaderPosMap.get("label_comp_id")];

        // residue number or rna number 

        // use auth_seq_id > label_seq_id (hope DSSP does so too)
        // resNumPDB = Integer.valueOf(lineData[colHeaderPosMap.get("auth_seq_id")]);
        molNumPDB = Integer.valueOf(lineData[colHeaderPosMap.get("auth_seq_id")]);
        
        // entity ID (mostly used for classification of ligands)
        entityID = Integer.valueOf(lineData[colHeaderPosMap.get("label_entity_id")]);

        // insertion code
        // only update if column and value exist, otherwise stick to blank ""
        if (colHeaderPosMap.get("pdbx_PDB_ins_code") != null) {
            if (! (lineData[colHeaderPosMap.get("pdbx_PDB_ins_code")].equals("?") || lineData[colHeaderPosMap.get("pdbx_PDB_ins_code")].equals("."))) {
                iCode = lineData[colHeaderPosMap.get("pdbx_PDB_ins_code")];
            }
        }

        // coordX
        // for information on difference between ptgl and PTGLgraphComputation style look in old parser
        if (Settings.getBoolean("PTGLgraphComputation_B_round_coordinates")) {
            oCoordXf = Float.valueOf(lineData[colHeaderPosMap.get("Cartn_x")]) * 10;
            coordX = Math.round(oCoordXf);
        } else {
            oCoordX = Double.valueOf(lineData[colHeaderPosMap.get("Cartn_x")]) * 10.0;
            coordX = oCoordX.intValue();
         }


        // coordY
        if (Settings.getBoolean("PTGLgraphComputation_B_round_coordinates")) {
            oCoordYf = Float.valueOf(lineData[colHeaderPosMap.get("Cartn_y")]) * 10;
            coordY = Math.round(oCoordYf);
        } else {
            oCoordY = Double.valueOf(lineData[colHeaderPosMap.get("Cartn_y")]) * 10.0;
            coordY = oCoordY.intValue();
        }

        // coordZ
        if (Settings.getBoolean("PTGLgraphComputation_B_round_coordinates")) {
            oCoordZf = Float.valueOf(lineData[colHeaderPosMap.get("Cartn_z")]) * 10;
            coordZ = Math.round(oCoordZf);
        } else {
            oCoordZ = Double.valueOf(lineData[colHeaderPosMap.get("Cartn_z")]) * 10.0;
            coordZ = oCoordZ.intValue();            
        }

        // chemical symbol
        chemSym = lineData[colHeaderPosMap.get("type_symbol")];

        // standard AAs and (some) non-standard, atm: UNK, MSE
        //   -> may be changed below if it is free (treat as ligand then)

        // TODO: possible to ignore alt loc atoms right now?

        if(FileParser.isDNAresidueName(FileParser.leftInsertSpaces(molNamePDB, 3))) {
            if( ! Settings.getBoolean("PTGLgraphComputation_B_no_parse_warn")) {
                DP.getInstance().w("Atom #" + atomSerialNumber + " in PDB file belongs to DNA residue (residue 3-letter code is '" + molNamePDB + "'), skipping.");
            }
            return;  // atom is not used
        }

        if( ! Settings.getBoolean("PTGLgraphComputation_B_include_rna")) {
            if(checkType(Molecule.RESIDUE_TYPE_RNA)) {
                if( ! Settings.getBoolean("PTGLgraphComputation_B_no_parse_warn")) {
                    DP.getInstance().w("Atom #" + atomSerialNumber + " in PDB file belongs to RNA residue (residue 3-letter code is '" + molNamePDB + "'), skipping.");
                }
                return;  // atom is not used
            }
        }

        // >> AA <<
        // update lastMol if the atom in the current line belongs to a new molecule than the previous line 
        //     -> enables getting DsspResNum for atom from res
        // match res <-> chain here
        // load new Residue into lastMol if we approached next Residue, otherwise only add new atom
        
        
        // entity stuff wird noch zu Problememn fürhern, weil in dssp_cif steht das glaub ich erst am Ende
        
        // todo chain ID prüfen
        // wird das schon von - - chain - - gemacht?
        if (!chainID.equals(lastMol.getChainID())){
            // new chain
            // lastChainID = chainID ?
        }
        
        if (!Objects.equals(molNumPDB, lastMol.getPdbNum())){
            // new residue
            res = new Residue();
            
            res.setPdbNum(molNumPDB);
            res.setType(Molecule.RESIDUE_TYPE_AA);  // kann man hier glaube ich noch nicht machen
            // DsspNum sollte nicht mehr benötigt werden
            res.setChainID(chainID);
            res.setiCode(iCode);
            res.setName3(molNamePDB);
            res.setAAName1(Residue.getAAName1fromAAName3(molNamePDB));
            res.setChain(FileParser.getChainByPdbChainID(chainID));
            res.setModelID(m.getModelID());
            // can't do that yet res.setSSEString(); Dazu muss ich glaube ich in den DSSP ranges schauen
            
            FileParser.s_molecules.add(res);
            
            lastMol = res;
            lastMol.setModelID(m.getModelID());
            lastMol.setChain(tmpChain);
            tmpChain.addMolecule(lastMol);
            
            // assign PDB res name (which differs in case of modifed residues)
            lastMol.setName3(molNamePDB);
            lastMol.setEntityID(entityID);
            
        }
        
        if (! (Objects.equals(molNumPDB, lastMol.getPdbNum()) && chainID.equals(lastMol.getChainID()) && iCode.equals(lastMol.getiCode()))) {
            tmpMol = FileParser.getResidueFromList(molNumPDB, chainID, iCode);  // null if not in DSSP data -> rna/ligand/free AA
            // check that a peptide residue could be found                   
            if (checkType(tmpMol.RESIDUE_TYPE_LIGAND) || entityInformation.get(String.valueOf(entityID)).get("type").equals("non-polymer")) {
                // residue is not in DSSP file and is not part of a chain -> must be free (modified) amino acid, ligand or RNA
                if (! silent) {
                    // print note only once
                    if (! molNumPDB.equals(lastLigandNumPDB))

                        if(Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 1) {
                            System.out.println("   PDB: Found a ligand, RNA or free (modified) amino acid at PDB# " + molNumPDB + ". Free amino acids are treated as ligands.");
                        }
                }

            } else {
                
                if(Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                    System.out.println("    [DEBUG LV 2] Found an amino acid at PDB# " + molNumPDB + " that is not listed in the DSSP file (might be at a chain break). Parsing it as part of a chain.");
                }
                
                // sometimes residues are missing from the dssp file if they are incomplete (mostly at chain breaks)
                // in this case, they have to be parsed here
                if (tmpMol == null) {
                    res = new Residue();
                                        
                    res.setPdbNum(molNumPDB);
                    res.setType(Molecule.RESIDUE_TYPE_AA);
                    
                    // assign fake dssp number taking into account other elements that have been given a dssp number
                    freeResTreatedNum++;
                    res.setDsspNum(assignDsspNum());

                    res.setChainID(chainID);
                    res.setiCode(iCode);
                    res.setName3(molNamePDB);
                    res.setAAName1(Residue.getAAName1fromAAName3(molNamePDB));
                    res.setChain(FileParser.getChainByPdbChainID(chainID));
                    res.setModelID(m.getModelID());
                    res.setSSEString("C");
                    
                    lastChainID = chainID;
                    FileParser.s_molecules.add(res);
                    lastMol = res;
                } else {
                    lastMol = tmpMol;
                }

                lastMol.setModelID(m.getModelID());
                lastMol.setChain(tmpChain);
                tmpChain.addMolecule(lastMol);

                // assign PDB res name (which differs in case of modifed residues)
                lastMol.setName3(molNamePDB);
                lastMol.setEntityID(entityID);
            }
        }

        Atom a = new Atom();

        if (FileParser.isIgnoredAtom(chemSym)) {
            if( ! (Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H"))) {
                if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
                    System.out.println("DEBUG Ignored atom line " + numLine.toString() + 
                            " as it is either in ignored list or handle_hydrogens turned off.");
                }
                return;
            }
        }

        // only ATOMs, not HETATMs, have a DSSP entry
        if((Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H"))) {
            a.setDsspResNum(null);
        }
        else {
            a.setDsspResNum(lastMol.getDsspNum());
        }
        
        if ((checkType(Molecule.RESIDUE_TYPE_RNA) && entityInformation.get(String.valueOf(entityID)).get("type").equals("polymer"))) {  // Nucleotides are only parsed as RNA if they are polymers, not if they act as single ligands.
            // >> RNA <<
            // if the line we are currently in belongs to the same molecule as the previous one, we only create a new atom for this line.
            // otherwise, a new RNA molecule is created
            if( ! ( molNumPDB.equals(lastRnaNumPDB) && chainID.equals(lastChainID) ) ) {
                rna = new RNA();
                rna.setPdbNum(molNumPDB);
                rna.setType(Molecule.RESIDUE_TYPE_RNA);                
                
                RnaTreatedNum++;
                rna.setDsspNum(assignDsspNum());
                
                rna.setChainID(chainID);
                rna.setiCode(iCode);
                rna.setName3(molNamePDB);
                rna.setAAName1(molNamePDB);
                rna.setChain(FileParser.getChainByPdbChainID(chainID));
                rna.setModelID(m.getModelID());
                rna.setSSEString("PTGLgraphComputation_S_rnaSseCode");
                rna.setEntityID(entityID);
                
                lastMol = rna;
                                
                if(FileParser.isIgnoredLigRes(molNamePDB)) {
                    // RNA chains can contain ligands that should be ignored such as water
                    // In this case, no new atom needs to be saved
                    RnaTreatedNum--;
                    a.setAtomtype(Atom.ATOMTYPE_IGNORED_LIGAND);
                    FileParser.ignoredLigands += 1;
                    
                    if(Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0){
                        DP.getInstance().w("Ignored RNA-Ligand found at PDB line " + molNumPDB + ". Name: " + molNamePDB);
                    }
                    return;
                    
                } else {                    
                    lastRnaNumPDB = molNumPDB;
                    lastChainID = chainID;
                    FileParser.s_molecules.add(rna);

                    Integer rnaIndex = FileParser.s_molecules.size() - 1;
                    FileParser.s_rnaIndices.add(rnaIndex);

                    FileParser.getChainByPdbChainID(chainID).addMolecule(rna);
                    
                    if(Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
                        if(! silent) {
                            DP.getInstance().d("New RNA molecule named " + molNamePDB + ", DSSPNumber " + rna.getDsspNum() + ", added in PDB line " + molNumPDB + " to chain " + chainID + ".");
                        }
                    }
                }
            }       
        }
        
        // If a molecule is not parsed at this point, it has to be a ligand
        else if (checkType(Molecule.RESIDUE_TYPE_LIGAND) ||                                          // if the molecule is categorized as a ligand through chem_comp map
                (entityInformation.get(String.valueOf(entityID)).get("type").equals("non-polymer"))) // if entity is defined as 'non-polymer' (e.g. free AA/RNA)
            
        {
            // >> LIG <<
                     
            // idea: add always residue (for consistency) but atom only if it is not an ignored ligand

            // check if we have created ligand residue for s_residue
            if( ! ( molNumPDB.equals(lastLigandNumPDB) && chainID.equals(lastChainID) ) ) {

                // create new Residue from info, we'll have to see whether we really add it below though
                lig = new Ligand();

                lig.setPdbNum(molNumPDB);
                lig.setType(Molecule.RESIDUE_TYPE_LIGAND);

                // assign fake DSSP Num increasing with each seen ligand
                ligandsTreatedNum ++;
                lig.setDsspNum(assignDsspNum());
                
                lig.setChainID(chainID);
                lig.setiCode(iCode);
                lig.setName3(molNamePDB);
                lig.setAAName1(AminoAcid.getLigandName1());
                lig.setChain(FileParser.getChainByPdbChainID(chainID));
                // still just assigning default model 1
                lig.setModelID(m.getModelID());
                lig.setSSEString(Settings.get("PTGLgraphComputation_S_ligSSECode"));
                lig.setEntityID(entityID);
                
                lastMol = lig;
                
                // add ligand to list of residues if it not on the ignore list
                if(FileParser.isIgnoredLigRes(molNamePDB)) {
                    ligandsTreatedNum--;    // We had to increment before to determine the fake DSSP res number, but this ligand won't be stored so decrement to previous value.
                    FileParser.ignoredLigands += 1;

                    if(Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0){
                        System.out.println("    PDB: Ignored ligand '" + molNamePDB + "-" + molNumPDB + "' at PDB line " + molNumPDB + ".");
                    }
                } else {

                    lig.setLigName((chemicalComponents.get(molNamePDB)).get("name"));
                    lig.setLigFormula((chemicalComponents.get(molNamePDB)).get("formula"));
                    lig.setLigSynonyms((chemicalComponents.get(molNamePDB)).get("pdbx_synonyms"));
                    
                    lastLigandNumPDB = molNumPDB;
                    lastChainID = chainID;

                    FileParser.s_molecules.add(lig);
                    Integer ligandIndex = FileParser.s_molecules.size() - 1;
                    FileParser.s_ligandIndices.add(ligandIndex);

                    FileParser.getChainByPdbChainID(chainID).addMolecule(lig);

                    // do we need this?
                    //resIndex = s_residues.size() - 1;
                    //resIndexDSSP[resNumDSSP] = resIndex;
                    //resIndexPDB[molNumPDB] = resIndex;      // This will crash because some PDB files contain negative residue numbers so fuck it.
                    if(! (FileParser.silent || FileParser.essentialOutputOnly)) {
                        System.out.println("   PDB: Added ligand monomer '" +  molNamePDB + "-" + molNumPDB + "', chain " + chainID + " (line " + numLine + ", ligand #" + ligandsTreatedNum + ", Fake DSSP #" + lig.getDsspNum() + ").");
                        System.out.println("   PDB:   => Ligand name = '" + lig.getLigName() + "', formula = '" + lig.getLigFormula() + "', synonyms = '" + lig.getLigSynonyms() + "'.");
                    }

                }
            }

            if(FileParser.isIgnoredLigRes(molNamePDB)) {
                a.setAtomtype(Atom.ATOMTYPE_IGNORED_LIGAND);       // invalid ligand (ignored)

                // We do not need these atoms and they may lead to trouble later on, so
                //  just return without adding the new Atom to any Residue here so this line
                //  is skipped and the next line can be handled.
                //  If people want all ligands they have to change the isIgnoredLigRes() function.

                // DEBUG
                // DP.getInstance().w("FP_CIF", " Ignored ligand atom of '" +  resNamePDB + "-" + molNumPDB + "', chain " + chainID + " (line " + numLine + ", ligand #" + ligandsTreatedNum + ", Fake DSSP #" + lig.getDsspResNum().toString() + ").");
                return; // can we do this here? Does it cut off other important stuff? -> added to res but not atom (s.a.)
            }
            else {
                a.setAtomtype(Atom.ATOMTYPE_LIGAND);       // valid ligand
                //a.setDsspResNum(getDsspResNumForPdbFields(molNumPDB, chainID, iCode));  // We can't do this because the fake DSSP residue number has not yet been assigned
            }
        }

        // now create the new Atom

        // lastMol may be NULL
        // Note that the command above may have returned NULL, we care for that below

        a.setPdbAtomNum(atomSerialNumber);
        a.setAtomName(atomName);
        a.setAltLoc(altLoc);
        a.setMolecule(lastMol);
        a.setChainID(chainID);        
        a.setChain(FileParser.getChainByPdbChainID(chainID));
        a.setPdbResNum(molNumPDB);
        // we cant get the DSSP res num easily here and have to do it later (I guess)
        // old parser seems to assign 0 here whatsoever so we just leave the default value there
        // a.setDsspResNum(resNumDSSP);
        a.setCoordX(coordX);
        a.setCoordY(coordY);
        a.setCoordZ(coordZ);
        a.setChemSym(chemSym);

        // from old parser, not working with models right now
        /*
        if(curModelID != null) {
            a.setModelID(curModelID);
            a.setModel(getModelByModelID(curModelID));
        }
        */
        
            if (lastMol.getType() == Molecule.RESIDUE_TYPE_AA || lastMol.getType() == Molecule.RESIDUE_TYPE_RNA){
                if (lastMol == null) {
                    DP.getInstance().w("Molecule with PDB # " + molNumPDB + " of chain '" + chainID + "' with iCode '" + iCode + "' not listed in CIF data, skipping atom " + atomSerialNumber + " belonging to that residue (PDB line " + numLine.toString() + ").");
                    return;
                } else {
                    if(Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H")) {
                        lastMol.addHydrogenAtom(a);
                    }
                    else {
                        // add Atom to list of atoms of current molecule as well as list of all atoms
                        FileParser.s_atoms.add(a);
                        if (checkType(Molecule.RESIDUE_TYPE_AA)){
                            if(Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                                System.out.println("    [DEBUG LV 2] New AA atom added: " + a.toString());
                            }
                            a.setAtomtype(Atom.ATOMTYPE_AA);
                            lastMol.addAtom(a);
                        }
                        if (checkType(Molecule.RESIDUE_TYPE_RNA)){
                            if(Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                                System.out.println("    [DEBUG LV 2] New RNA atom added: " + a.toString());
                            }
                            a.setAtomtype(Atom.ATOMTYPE_RNA);
                            a.setMolecule(rna);
                            rna.addAtom(a);
                        }
                    }
                }
            }
            else {
                if (! (lig == null)){
                    if(Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                       System.out.println("    [DEBUG LV 2] New ligand atom added: " + a.toString());
                    }
                    lig.addAtom(a);
                    a.setMolecule(lig);
                    FileParser.s_atoms.add(a);
                }

            }
    }
            
    
    /**
     * Handles lines starting with _chem_comp by filling the chemicalComponents HashMap.
     * Entries can look like this: MET={name=METHIONINE, pdbx_synonyms=?, formula=C5 H11 N O2 S, id=MET, type=L-peptide linking, formula_weight=149.211}}
     */
    private static void handleChemComp(){
        String category = null;              // categories such as type, name
        String value = null;              // values such as peptide, Methionine
        HashMap<String, String> tmpComponent = new HashMap<>();     // stores categories and values for current component
        for (String cat : colHeaderPosMap.keySet()){
            category = cat;
            value = lineData[colHeaderPosMap.get(cat)];
            tmpComponent.put(category, value);
            chemicalComponents.put(lineData[colHeaderPosMap.get("id")], tmpComponent);      // matches one component with all its category/value pairings
        }
    }
    
    
    /**
     * Handles lines starting with '_entity_src_gen'.
     * Assigns the common and the scientific name of the molecule to variables that are later transferred to ProtMetaInfo.
     */
    private static void handleEntitySrcGen(){
        if (colHeaderPosMap.get("gene_src_common_name") != null ){
            if (! lineData[colHeaderPosMap.get("gene_src_common_name")].equals("?")){
                nameOrgCommon = lineData[colHeaderPosMap.get("gene_src_common_name")];
                nameOrgCommonSource = "_entity_src_gen";
            }
        }
        if (colHeaderPosMap.get("pdbx_gene_src_scientific_name") != null ){
            if (! lineData[colHeaderPosMap.get("pdbx_gene_src_scientific_name")].equals("?")){
                nameOrgScientific = lineData[colHeaderPosMap.get("pdbx_gene_src_scientific_name")];
                nameOrgScientificSource = "_entity_src_gen";
            }
        }
    }
    
    
    /**
     * Handles Lines starting with '_entity_src_nat.'.
     * Alternative way to parse organism name.
     */
    private static void handleEntitySrcNat(){
        if (colHeaderPosMap.get("common_name") != null && ! nameOrgCommonSource.equals("_entity_src_gen")){
            if (! lineData[colHeaderPosMap.get("common_name")].equals("?")){
                nameOrgCommon = lineData[colHeaderPosMap.get("common_name")];
                nameOrgCommonSource = "_entity_src_nat";
            }
        }
        if (colHeaderPosMap.get("pdbx_organism_scientific") != null && ! nameOrgScientificSource.equals("_entity_src_gen")){
            if (! lineData[colHeaderPosMap.get("pdbx_organism_scientific")].equals("?")){
                nameOrgScientific = lineData[colHeaderPosMap.get("pdbx_organism_scientific")];
                nameOrgScientificSource = "_entity_src_nat";
            }
        }
    }
    
    
    /**
     * Handles lines starting with '_pdbx_entity_src_syn'.
     * Alternative way to parse organism name.
     */
    private static void handlePdbxEntitySrcSyn(){
        if (colHeaderPosMap.get("organism_common_name") != null && ! nameOrgCommon.equals("")){
            if (! lineData[colHeaderPosMap.get("organism_common_name")].equals("?")){
                nameOrgCommon = lineData[colHeaderPosMap.get("organism_common_name")];
                nameOrgCommonSource = "_pdbx_entity_src_syn";
            }
        }
        if (colHeaderPosMap.get("organism_scientific") != null && ! nameOrgScientific.equals("")){
            if (! lineData[colHeaderPosMap.get("organism_scientific")].equals("?")){
                nameOrgScientific = lineData[colHeaderPosMap.get("organism_scientific")];
                nameOrgScientificSource = "_pdbx_entity_src_syn";
            }
        }
    }

    
    /**
     * Checks the string value for being not null, empty or an mmCIF placeholder.
     * @param value
     * @return 
     */
    private static boolean valueIsAssigned(String value) {
        if (value == null) {
            return false;
        } else {
            if (value.length() == 0) {
                return false;
            }
        }
        
        for (String placeholder : NO_VALUE_PLACEHOLDERS) {
            if (placeholder.equals(value)) {
                return false;
            }
        }
        return true;
    }
    
    
    /**
     * Trims ' or ; when they exist as first an last char, then removes \n.
     * 
     */
    private static String[] trimSpecialChars(ArrayList<String> stringList) {
        ArrayList<String> charsToDelete = new ArrayList<>(Arrays.asList((";"), ("\'")));
        String[] returnList = new String[stringList.size()];
        String element;
        for (int i = 0; i < stringList.size(); i++) {
            element = stringList.get(i);  // start with respective element
            if (stringList.get(i).length() > 0) {
                // First check if there are special chars around the line
                if (charsToDelete.contains(String.valueOf(stringList.get(i).charAt(0))) && charsToDelete.contains(String.valueOf(stringList.get(i).charAt((stringList.get(i).length()) - 1)))){
                    element = stringList.get(i).substring(1, stringList.get(i).length() - 1);
                    if (element.endsWith("\n")){
                        // Secondly, check if there is still a newline char at the end of the line
                        element = String.valueOf(element).replaceAll("\n", "");
                    }
                }
            }
            returnList[i] = element;
        }
        return returnList;
    }
    
    
    /**
     * Returns an array of 'words' separated by an arbitrary amount of spaces. Considers in-line strings.
     * @param handleSpecialCharacters whether single and double quotation should be considered.
     *   Should be false whenever a combined multi-line string is provided and true otherwise.
     * @param line
     * @return
     */
    private static String[] lineToArray(String line) {
        ArrayList<String> dataItemList = new ArrayList<>();
        String tmpItem = "";
        boolean inLineString = false;  // marked by single-quoation mark
        boolean inLineQuote = false;
        boolean inMultiString = false;  // surrounded by semi-colon and hides everything
        
        // used to  determine whether single-quotation starts in-line string
        Character prev_char;
        Character next_char;
        
        for (int i = 0; i < line.length(); i++) {
            
            char tmpChar = line.charAt(i);
            
            if (inMultiString) {
                // accept everything but ending newline-character+semi-colon or newline-character
                
                // end of multi-line string
                if (tmpChar == ';' && line.charAt(i - 1 ) == '\n') {
                    // newline can only appear before multi-line-ending semi-colon, so it should be okay to assume there is always a char before
                    inMultiString = false;
                    continue;
                }
                
                // newline-character only used to mark multi-line ending semi-colon, so exclude it from data
                if (tmpChar == '\n') { continue; }
            } else {
            
                if (inLineQuote) {
                    // accept everything but ending quote
                    if (tmpChar == '"') {
                        inLineQuote = false;
                        continue;
                    }
                } else {
                    
                    // if there is previous / next character assign it, else assign ' ', i.e., handle single quoations on start / end as in-line string
                    prev_char = (i > 0 ? line.charAt(i - 1) : ' ');
                    next_char = (i < line.length() - 1 ? line.charAt(i + 1) : ' ');
                    
                    if (inLineString) {
                        // accept everything but ending single-quote
                        if (tmpChar == '\'' && next_char == ' ') {
                            inLineString = false;
                            continue;
                        }
                    } else {
                        // we are not in multi-line, single-line string or quotation

                        switch (tmpChar) {
                            case ' ':
                                if (tmpItem.length() > 0) {
                                    dataItemList.add(tmpItem);                   
                                    tmpItem = "";
                                }
                                break;
                            case '"':
                                inLineQuote = true;
                                break;
                            case ';':
                                if (prev_char == ' ' || next_char == ' '){  // in case a ';' is embedded in between two characters it is assumed it is part of a string
                                    inMultiString = true;
                                }
                                break;
                            case '\'':
                                if (prev_char == ' ') {
                                    inLineString = true;
                                } else {
                                    // in-word single-quotation, e.g., 5'-O
                                    tmpItem += tmpChar;
                                }
                                break;
                            default:
                                tmpItem += tmpChar;
                        }  // end of switch-case
                        continue;
                    }  // end of in-lineString
                }  // end of if in-quote
            }  // end of if in-multiString
            // only reached when in some flagged part and end not reached, so simply add character
            tmpItem += tmpChar;
        }  // end of for-loop over line's characters
        
        // lines appear to end with a whitespace, still care for when they do not
        if (tmpItem.length() > 0) {
            dataItemList.add(tmpItem);
        }
        
        if (inLineString) {
            DP.getInstance().w("FP_CIF", "In line " + numLine + " an not ending inline-string (single quotation mark) was found. Is the file correct? Trying to ignore it.");
        }

        // return Array instead of list
        return dataItemList.toArray(new String[dataItemList.size()]);
    }
    
    
    /**
     * Fills important metaData fields that are not existing with empty Strings.
     */
    private static void fillMetadataEmptyStrings() {
        List<String> mdFields = Arrays.asList("title", "keywords", "experiment", "resolution", "date", "header");
        for (String field : mdFields) {
            if (! metaData.containsKey(field)) {
                metaData.put(field, "");
            }
        }
    }
    
    
    /**
     * Fills allProteinMetaInfos from entityInformation. PMIs already have to be created (e.g. during parsing and creation of chains).
     */
    private static void fillProteinMetaData() {
        String tmpValue;
        ArrayList<String> tmpArrayList;
        HashMap <String, String> tmpEntityInfo;

        for (ProtMetaInfo pmi : allProteinMetaInfos) {
            tmpEntityInfo = entityInformation.get(pmi.getMacromolID());

            // pdbx_description -> molName
            tmpValue = tmpEntityInfo.get("pdbx_description");
            if (tmpValue != null) {
                pmi.setMolName(tmpValue);
            } 
            
            // pdbx_ec -> ecnumber
            tmpValue = tmpEntityInfo.get("pdbx_ec");
            if (tmpValue != null) {
                pmi.setECNumber(tmpValue);
            }
            
            pmi.setOrgCommon(nameOrgCommon);
            pmi.setOrgScientific(nameOrgScientific);
            
            // FileParser.homologuesMap -> allMolChains
            //   seems like this is never used, but cant hurt to fill it, since we have the information
            tmpArrayList = FileParser.homologuesMap.get(pmi.getChainid());
            if (tmpArrayList != null) {
                
                tmpValue = tmpArrayList.toString().replace("[", "").replace("]", "");
                
                tmpValue = (tmpValue.length() > 0 ? tmpValue + ", " + pmi.getChainid() : pmi.getChainid());  // include your own ID as it says '>all', but only append if there are other chains              
                pmi.setAllMolChains(tmpValue);
            }   
        }       
    }
    
    /**
     * Checks for the presence of predefined (hard coded) column headers.
     * @param categoryName The name of the category
     * @param columnHeaders The column headers found in the file
     * @return A list of all required columns that were missing
     */
    private static ArrayList<String> checkColumns(String categoryName, ArrayList<String> columnHeaders) {
        // different columns are expected depending on category
        // define them here:
        String[] reqColumns;
        switch (categoryName) {
            case "_atom_site":
                reqColumns = new String[] {"id", "type_symbol", "label_atom_id", "label_comp_id", 
                    "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z"};
                break;
            default:
                if (! Settings.getBoolean("PTGLgraphComputation_B_no_warn")) {
                    DP.getInstance().w("FP_CIF", "Tried to check table of category " + categoryName +
                            " for presence of important columns, but function is not defined for that " +
                            "category. Ignoring the check and moving on.");
                }
                reqColumns = new String[0];
        }
        
        ArrayList<String> missingColumns = new ArrayList<>();

        Boolean found;
        for (String reqColumn : reqColumns) {
            found = false;
            for (String colHeader: columnHeaders) {
                if (colHeader.equals(reqColumn)) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                    missingColumns.add(reqColumn);
            }
        }
        return missingColumns;
    }
    
    
    /**
     * Handles line defining loop's definition and adds it to colHeaderPosMap.
     * @param line 
     */
    private static void addColumnHeaderToMap(String line) {
        if (inLoop){
            String[] elementsSeperatedByDot = line.split("\\.");

            if (elementsSeperatedByDot.length < 2) {
                DP.getInstance().w("FP_CIF", " Expected table definition in line " + 
                    numLine.toString() + " but couldnt parse it. Skip it (may miss important data!).");
                return;
            }
            // from now on elementsSeperatedByDot has at least two elements

            colHeaderPosMap.put(elementsSeperatedByDot[1].trim(), colHeaderPosMap.size());
        }
    }

    
    /**
     * Creates from a list all aub lists excluding one different element.
     * @param targetArray
     * @return 
     */
    private static ArrayList<ArrayList<String>> arrayToSubListsWithoutOne(String[] targetArray) {
        ArrayList<ArrayList<String>> subLists = new ArrayList<>();
        ArrayList<String> tmpList = new ArrayList<>();
        
        for (int i = 0; i < targetArray.length; i++) {
            tmpList = new ArrayList<>();
            for (int j = 0; j < targetArray.length; j++) {
                if (i != j) {
                    // sub list shall exclude one element
                    tmpList.add(targetArray[j]);
                    System.out.println("tmpList: " + tmpList);
                }
            }
            subLists.add(tmpList);
        }
        
        return subLists;
    }
    
    
    /**
     * Gets chain by ID if existing otherwise creates it. Also creates the ProtMetaInfo for this chain if necessary.
     * @param cID chain ID as String
     * @param m Model to which the chain belongs
     * @param entityID the entity_id from a cif file describing the macromolecular ID from legacy file
     * @return 
     */
    private static Chain getOrCreateChain(String cID, Model m, String entityID) {
        for (Chain existing_c : FileParser.s_chains) {
            if (existing_c.getPdbChainID().equals(cID)) {
                return existing_c;
            }
        }
        
        // reaching this code only if chain didnt exist
        //   that means that also not ProtMetaInfo exists, so create one here
        Chain c = new Chain(cID);      
        ProtMetaInfo pmi = new ProtMetaInfo(pdbID, cID);
        
        c.setModel(m);
        c.setModelID(m.getModelID());
        m.addChain(c);
        if (chainIdentity.get(cID) == null) {
            c.setMoleculeType("non-polymer");
        } else {
            c.setMoleculeType(chainIdentity.get(cID));
        }
        c.setAltChainID(altChainID);
                
        c.setHomologues(FileParser.homologuesMap.get(cID));
        //c.setMacromolID(entityID);
        pmi.setMacromolID(entityID);
        
        FileParser.s_chains.add(c);
        if (! (FileParser.silent || FileParser.essentialOutputOnly)) {
            System.out.println("   PDB: New chain named " + cID + " found.");
        }
        
        allProteinMetaInfos.add(pmi);
        return c;
    }
        

    protected static ProtMetaInfo getProteinMetaInfo(String pdbID, String chainID) {
        Boolean foundPMI = false;
        ProtMetaInfo pmi = null;
        Integer currentIndex = lastIndexProtMetaInfos;

        // iterate up to allProteinMetaInfos.size() times
        for(Integer i = 0; i < allProteinMetaInfos.size(); i++) {
                       
            // start at last occurence
            currentIndex = (lastIndexProtMetaInfos + i) % allProteinMetaInfos.size();
            
            pmi = allProteinMetaInfos.get(currentIndex);
                      
            if(pmi.getPdbid().equals(pdbID)  && pmi.getChainid().equals(chainID)) {
                foundPMI = true;
                break;
            }
        }
        
        if (! foundPMI) {
            pmi = new ProtMetaInfo(pdbID, chainID);
            // if no matching pmi was found (should not happen), return empty pmi and throw warning
            DP.getInstance().w("No protein chain meta information for PDB ID: " + pdbID + " and chain ID " + chainID + " found."
                + " Returning empty informtation instead.");
        }
        lastIndexProtMetaInfos = currentIndex;
        return pmi;
    }
    
    
    /**
     * Returns whether the molecule type matches the requested one.
     * Checks whether the the molecule type in the chemical components map corresponds to the requested one.
     * Careful with this function! There are amino acid ligands that are not always reliably categorized as ligands but instead as amino acids.
     * @param requestedType integer of the type one is looking for (0 for AA, 1 for ligand, 3 for RNA)
     */
    protected static Boolean checkType(Integer requestedType){
        if (! (requestedType == 0 || requestedType == 1 || requestedType == 2 || requestedType == 3)){
            DP.getInstance().w("Tried to check molecule type, but requested type was not recognized. Trying to move on without checking the type.");
            return false;
        }
        Integer actualType = 1;
        String chemType = ((chemicalComponents.get(molNamePDB)).get("type"));
        String chemTypeLowerCase = chemType.toLowerCase();
        int intIndex = chemTypeLowerCase.indexOf("rna");
        if (intIndex != -1){
            actualType = 3;
        }
        else{
            intIndex = chemTypeLowerCase.indexOf("peptide linking");  // if it says "peptide" only, it might be a ligand and not linking
            if (intIndex != -1){
                actualType = 0;
            }
        }
        return (actualType == requestedType) ? true : false;
    }
    
    
    /**
     * Assigns a fake DSSP-Number based on all objects that have received DSSP-Numbers.
     * Only AAs in DSSP file get a real DSSP number, all other objects (RNA, ligands, free AAs) have to be assigned one.
     */
    protected static Integer assignDsspNum () {
        return DsspParser.lastUsedDsspNum + RnaTreatedNum + ligandsTreatedNum + freeResTreatedNum;
    }
      
}
