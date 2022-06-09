/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;

// imports
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.TreeMap;
import settings.Settings;
import proteingraphs.MolContactInfo;
import proteingraphs.ProtGraph;
import proteinstructure.AminoAcid;
import proteinstructure.Atom;
import proteinstructure.BindingSite;
import proteinstructure.Chain;
import proteinstructure.Model;
import proteinstructure.ProtMetaInfo;
import proteinstructure.Residue;
import proteinstructure.Ligand;
import resultcontainers.ProteinResults;
import tools.DP;


/**
 * File parser for legacy PDB file format.
 * 
 * @author jnw
 */
class LegacyParser {
    static Integer maxResidues = 11000;  // could be useless, but was included in pre-dissolved file parser
    static String firstModelName = "1";  // the model ID that identifies the first model in a PDB file
    static String defaultModelName = firstModelName;
    static boolean dataInitDone;
    static boolean silent;
    static boolean essentialOutputOnly;
    
    static Integer curLineNumPDB = null;
    static String curLinePDB = null;
    static String curModelID = null;
    static String oldModelID = null;
    static String curChainID = null;
    static String oldChainID = null;
    static Integer curResNumPDB = null;
    static Residue curRes = null;
    
    
    /**
     * @param df Path to a DSSP file. Does NOT test whether it exist, do that earlier.
     * @param pf Path to a PBD file. Does NOT test whether it exist, do that earlier.
     */
    protected static Boolean initData(String pf) {
        
        silent = FileParser.settingSilent();
        essentialOutputOnly = FileParser.essentialOutputOnly;

        // resIndexPDB = new Integer[maxResidues];      // Removed because some PDB files have negative residue numbers, they break this. :/ So we
        //                                              //  have to go through the whole list (which is slow and stupid, bah).
        FileParser.resIndexDSSP = new Integer[maxResidues];

        
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
     * Parses all data, goes through all the lines and calls the appropriate function to handle each line.
     * @return ignore
     */
    private static Boolean parseData() {

        if(! silent) {
            System.out.println("  Parsing pdb and dssp file lines...");
        }

        // init class vars
        curLineNumPDB = 0;
        curLinePDB = "";
        curModelID = defaultModelName; 
        curChainID = " ";
        oldChainID = " ";
        oldModelID = "";

        if(! silent) {
            System.out.println("  Scanning whole PDB file for models...");
        }
        createAllModelIDsFromWholePdbFile();   // fills s_allModelIDsFromWholePDBFile

        // create all models, chains and residues first. parse atoms afterwards.
        if(! silent) {
            System.out.println("  Creating all Models from handled PDB lines...");
        }
        createAllModelsFromHandledPdbLines();   // fills s_models

        if(FileParser.s_models.size() > 1) {
            System.out.println("ERROR: Found > 1 model (" + FileParser.s_models.size() + " to be precise) models in the parsed PDB file lines, something went wrong. Exiting.");
            System.exit(1);
        }
        
        if(! silent) {
            System.out.println("  Setting chain homologues...");
        }
        setHomologuesOfChains();   // fills homologuesMap     

        if(! silent) {
            System.out.println("  Creating all Chains...");
        }
        createAllChainsFromPdbData();   // fills s_chains
        
        if(Settings.getBoolean("PTGLgraphComputation_B_parse_binding_sites")) {
            if(! silent) {
                System.out.println("  Creating binding sites...");
            }
                createAllBindingSitesFromPdbData(); // fills s_bindingsites
            if(! (silent || essentialOutputOnly)) {
                for(BindingSite s : FileParser.s_sites) {
                    System.out.println("    PDB: " + s.toString());
                }
                System.out.println("    Found " + FileParser.s_sites.size() + " binding sites.");
            }
        }
        

        if(! silent) {
            System.out.println("  Creating all Residues...");
        }
        
        DsspParser.createAllResiduesFromDsspData(false);    // fills s_residues

        // If there is no data part at all in the DSSP file, the function readDsspToData() will catch
        //  this error and exit, this code will never be reached in that case.
        if(FileParser.s_molecules.size() < 1) {
            System.err.println("ERROR: DSSP file contains no residues (maybe the PDB file only holds DNA/RNA data). Exiting.");
            System.exit(2);
        }
        if(! silent) {
            System.out.println("  Creating all Ligand Residues...");
        }
        createAllLigandResiduesFromPdbData();       // adds stuff to s_residues

        if(! silent) {
            System.out.println("  Creating all SSEs according to DSSP definition...");
        }
        //s_dsspSSEs = createAllSSEsFromResidueList();     // fills s_dsspSSEs

        //System.out.println("   Done with all DSSP SSEs, ligand SSEs not created yet:");
        //for(Integer i = 0; i < s_dsspSSEs.size(); i++) {
        //    System.out.println("    SSE #" + i + ": " + s_dsspSSEs.get(i) + ".");
        //}

        //System.out.println("  Creating all ligand SSEs...");
        //createAllLigandSSEsFromResidueList();       // adds stuff to s_dsspSSEs

        //System.out.println("  Creating modified SSE list according to PTGL definition...");
        //createAllPtglSSEsFromSSEList();         // fills s_ptglSSEs


        if(! silent) {
            System.out.println("  Creating all Atoms...");
        }
        Boolean ignoreRestOfFile = false;
       
        // parse all PDB lines based on format definitions at http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {

            curLineNumPDB = i + 1;
            curLinePDB = FileParser.pdbLines.get(i);

            //System.out.println("  PDB line " + currentLineNumPDB + ": " + currentLinePDB + "");


            // handle MODEL lines
            if(curLinePDB.startsWith("MODEL ")) {
                handlePdbLineMODEL();
            }

            // handle any ATOM lines (ATOM, HETATM)
            else if(curLinePDB.startsWith("ATOM  ") || curLinePDB.startsWith("HETATM")) {
                handlePdbLineANYATOM();
            }

            //  handle TER lines
            else if(curLinePDB.startsWith("TER   ")) {
                handlePdbLineTER();
            }

            //  handle all other lines
            else {
                handlePdbLineOTHER();
            }

        }

        if(! (silent || essentialOutputOnly)) {
            System.out.println("    PDB: Hit end of PDB file at line " + curLineNumPDB + ".");

            // remove duplicate atoms from altLoc here
            System.out.println("    PDB: Selecting alternative locations for atoms of all residues.");
        }
        ArrayList<Atom> deletedAtoms;
        int numAtomsDeletedAltLoc = 0;
        int numResiduesAffected = 0;
        Residue r;
        for(int i = 0; i < FileParser.s_molecules.size(); i++) {
            if (FileParser.s_molecules.get(i) instanceof Residue) {
                r = (Residue) FileParser.s_molecules.get(i);
                deletedAtoms = r.chooseYourAltLoc();

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
        }
        
        // assign residues to binding sites
        if(Settings.getBoolean("PTGLgraphComputation_B_parse_binding_sites")) {
            Integer siteResPdbResNum; String siteResChainID; String siteResName;
            int numResAssigned = 0;
            int numWaterResIgnored = 0;
            for(BindingSite s : FileParser.s_sites) {
                List<String[]> siteResidueInfos = s.getResidueInfos();
                // try-catch for coping with the Int/Str-ParseError if insertion codes
                //    are used in the SITE fields
                try {
                    for(String[] resInfo : siteResidueInfos) {

                        siteResPdbResNum = Integer.parseInt(resInfo[2]);

                        siteResChainID = resInfo[1];
                        siteResName = resInfo[0];
                        Residue res = FileParser.getResByPdbFields(siteResPdbResNum, siteResChainID, null);
                        if(res != null) {
                            res.addPartOfBindingSite(s);
                            numResAssigned++;
                        } else {
                            if(siteResName.toUpperCase().equals("HOH")) {
                                numWaterResIgnored++;
                            }
                            else {
                                DP.getInstance().w("FileParser", "Could not assign residue " + siteResName + " #" + siteResPdbResNum + " of chain " + siteResChainID + " to binding site, no such residue found in residue list.");
                            }
                        }
                    }
                } catch (NumberFormatException e) {
                    DP.getInstance().w("FP", "NumberFormatException while parsing the PDB-File.");
                    DP.getInstance().w("FP", "Trying to go on now. Set PTGLgraphComputation_B_parse_binding_sites=false "
                            + "if insertion codes are used in the SITE fields "
                            + "(skips the detection of binding sites then).");
                }                
            }
            if(! (silent || essentialOutputOnly)) {
                System.out.println("    PDB: Assigned " + numResAssigned + " residues to be part of a binding site (ignored " + numWaterResIgnored + " waters).");
            }
        }
        
        
        if(! (silent || essentialOutputOnly)) {
            System.out.println("    PDB: Deleted " + numAtomsDeletedAltLoc + " duplicate atoms from " + numResiduesAffected + " residues which had several alternative locations.");

            // report statistics
            System.out.println("  All data parsed. Found " + FileParser.s_models.size() + " models, " +
                                                           FileParser.s_chains.size() + " chains, " +
                                                           FileParser.s_molecules.size() + " molecules, " +
                                                           FileParser.s_atoms.size() + " atoms.");
        }

        return(true);
    }
    
    
    /**
     * Handles PDB MODEL lines.
     * @return ignored
     */
    private static boolean handlePdbLineMODEL() {
        // Now handled by function createAllModelsFromPdbData() because the models have to exist
        //  before we create the other stuff (atoms, residues, etc.) so we can assign the models to them.
        // This function still updates the global 'curModelID' variable that's used by other functions (the one that
        // creates that atoms) though, so it is still required.
        String mID = "";

        try {
            mID = (curLinePDB.substring(10, 16)).trim();
        } catch(Exception e) {
            System.err.println("ERROR: Hit MODEL line at PDB line number " + curLinePDB + " but parsing the line failed.");
            e.printStackTrace();
            System.exit(1);
        }

        curModelID = mID;
        return(true);
    }

    /**
     * Handle PDB ATOM and HETATM lines. Currently this function creates all the atoms. It tracks chains etc via global variables. It would be better if this was 
     * turned into a function createAllAtomsFromPdbData() that works like the createAll* functions.
     * @return true if the line could be parsed, false otherwise
     */
     private static boolean handlePdbLineANYATOM() {
        

        Integer atomSerialNumber, resNumPDB, resNumDSSP, rnaNumPDB;
        Integer coordX, coordY, coordZ;
        atomSerialNumber = resNumPDB = resNumDSSP = coordX = coordY = coordZ = 0;
        String atomRecordName, atomName, resNamePDB, chainID, chemSym, altLoc, iCode;
        atomRecordName = atomName = resNamePDB = chainID = chemSym = altLoc = iCode = "";
        Double oCoordX, oCoordY, oCoordZ;            // the original coordinates in Angstroem (coordX are 10th part Angstroem)
        oCoordX = oCoordY = oCoordZ = 0.0;
        Float oCoordXf, oCoordYf, oCoordZf;
        oCoordXf = oCoordYf = oCoordZf = 0.0f;

        try {
            atomRecordName = curLinePDB.substring(0, 6).trim();
            atomSerialNumber = Integer.valueOf((curLinePDB.substring(6, 11)).trim());
            // 11 is ignored: blank
            atomName = curLinePDB.substring(12, 16);
            altLoc = curLinePDB.substring(16, 17);
            resNamePDB = curLinePDB.substring(17, 20);
            // 20 is ignored: blank
            chainID = curLinePDB.substring(21, 22);
            resNumPDB = Integer.valueOf((curLinePDB.substring(22, 26)).trim());
            iCode = curLinePDB.substring(26, 27);       // don't trim this!
            // 27 - 29 are ignored: blanks
            
            
  
            if(Settings.getBoolean("PTGLgraphComputation_B_round_coordinates")) {
            // PTGLgraphComputation style: round the coordinates
                oCoordXf = Float.valueOf((curLinePDB.substring(30, 38)).trim()) * 10;
                oCoordYf = Float.valueOf((curLinePDB.substring(38, 46)).trim()) * 10;
                oCoordZf = Float.valueOf((curLinePDB.substring(46, 54)).trim()) * 10;
                coordX = Integer.valueOf(Math.round(oCoordXf));
                coordY = Integer.valueOf(Math.round(oCoordYf));
                coordZ = Integer.valueOf(Math.round(oCoordZf));
                
                // now, is has become 88
                //System.out.println("[atom#" + atomSerialNumber + "] oCoords=(" + oCoordXf + "," + oCoordYf + "," + oCoordZf + "), Coords=(" + coordX + "," + coordY + "," + coordZ + ")");
                
            }
            else {
                // PTGL style: always round them down, i.e., simply ignore the last 2 digits
                oCoordX = Double.valueOf((curLinePDB.substring(30, 38)).trim()) * 10.0;
                oCoordY = Double.valueOf((curLinePDB.substring(38, 46)).trim()) * 10.0;
                oCoordZ = Double.valueOf((curLinePDB.substring(46, 54)).trim()) * 10.0;
                
                // now, an example PDB file value of "8.796" was transformed into -87.96f                
                
                coordX = oCoordX.intValue();
                coordY = oCoordY.intValue();
                coordZ = oCoordZ.intValue();                                
                
                //System.out.println("[atom#" + atomSerialNumber + "] oCoords=(" + oCoordX + "," + oCoordY + "," + oCoordZ + "), Coords=(" + coordX + "," + coordY + "," + coordZ + ")");
                
            }                        
 
            // 54 - 59 are ignored: occupancy
            // 60 - 65 are ignored: temp factor
            // 66 - 72 are ignored: blanks
            // 72 - 75 are ignored: segment identifier
            chemSym = curLinePDB.substring(76, 78);
            
            //if(chemSym.trim().equals("H")) {
            //    System.out.println("Found hydrogen line '" + curLinePDB + "'.");
            //}

            // 78 - 79 are ignored: atom charge
        } catch(Exception e) {
            System.err.println("ERROR: Hit ATOM/HETATM line at PDB line number " + curLineNumPDB + " but parsing the line failed (length " + curLinePDB.length() + "): '" + e.getMessage() + "'.");
            return false;
            //System.exit(1);
        }

        // Ignore alternate location identifiers (only use the 1st one. It is identified by an altLoc field that is empty, "A" or "1").
        //if( ! (altLoc.equals(" ") || altLoc.equals("A") || altLoc.equals("1") )) {
        //    //System.out.println("INFO: Ignored alternate atom location identifier '" + altLoc + "' for atom #" + atomSerialNumber + " in PDB file.");
        //    return(false);
        //}

        Atom a = new Atom();

        // handle stuff that's different between ATOMs and HETATMs
        if(atomRecordName.equals("ATOM")) {

            if(FileParser.isIgnoredAtom(chemSym)) {
                //a.setAtomtype(3);
                
                if( ! (Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H"))) {
                    return(false);
                } 
            }

            // the N terminus is in an ATOM line, never in a HETATM line
            // TODO: How to determine the N terminus? I can't see anything special about the N terminus lines, the code below finds all N.
            //if(atomName.equals(" N  ")) {
            //    System.out.println("  Found N terminus at PDB line number " + curLineNumPDB + ".");
            //}

            // set atom type
            a.setAtomtype(Atom.ATOMTYPE_AA);

            // only ATOMs, not HETATMs, have a DSSP entry
            //a.setDsspResNum(getDsspResNumForPdbResNum(resNumPDB));
            if((Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H"))) {
                a.setDsspResNum(null);
            }
            else {
                a.setDsspResNum(DsspParser.getDsspResNumForPdbFields(resNumPDB, chainID, iCode));
            }
        }
        else {          // HETATM

            if(FileParser.isIgnoredLigRes(resNamePDB)) {
                a.setAtomtype(Atom.ATOMTYPE_IGNORED_LIGAND);       // invalid ligand (ignored)
                FileParser.ignoredLigands += 1;

                // We do not need these atoms and they may lead to trouble later on, so
                //  just return without adding the new Atom to any Residue here so this line
                //  is skipped and the next line can be handled.
                //  If people want all ligands they have to change the isIgnoredLigRes() function.
                return(false);
            }
            else {
                a.setAtomtype(Atom.ATOMTYPE_LIGAND);       // valid ligand
                //a.setDsspResNum(getDsspResNumForPdbFields(resNumPDB, chainID, iCode));  // We can't do this because the fake DSSP residue number has not yet been assigned
            }
        }

        // now create the new Atom        
        Residue tmpRes = FileParser.getResidueFromList(resNumPDB, chainID, iCode);
        // Note that the command above may have returned NULL, we care for that below
       
        a.setPdbAtomNum(atomSerialNumber);
        a.setAtomName(atomName);
        a.setAltLoc(altLoc);
        a.setMolecule(tmpRes);
        a.setChainID(chainID);        
        a.setChain(FileParser.getChainByPdbChainID(chainID));
        a.setPdbResNum(resNumPDB);
        a.setDsspResNum(resNumDSSP); // always 0 atm (default value and changed later?)
        a.setCoordX(coordX);
        a.setCoordY(coordY);
        a.setCoordZ(coordZ);
        a.setChemSym(chemSym);

        if(curModelID != null) {
            a.setModelID(curModelID);
            a.setModel(FileParser.getModelByModelID(curModelID));
        }
        
        if(tmpRes == null) {
            DP.getInstance().w("Residue with PDB # " + resNumPDB + " of chain '" + chainID + "' with iCode '" + iCode + "' not listed in DSSP data, skipping atom " + atomSerialNumber + " belonging to that residue.");
            return(false);
        } else {            
            
            if(Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce") && chemSym.trim().equals("H")) {
                tmpRes.addHydrogenAtom(a);
            }
            else {
                tmpRes.addAtom(a);
                FileParser.s_atoms.add(a);
            }
        }

        return(true);
    }
     
     
    /**
    * Parses a TER record line of a PDB file and writes info to STDOUT. 
    * @return true
    */
    private static boolean handlePdbLineTER() {

        
        Integer cTerminusResNumPDB = -1;
        String cID = "";
        String iCode = " ";

        try {
            
            // Ugly hack for the non-standard PDB files produced by the 'reduce' hydrogen program follows.
            // Reduce rewrites the PDB file and does NOT fill the empty lines with spaces, so we need to cope with shorter line lengths here.
            if(curLinePDB.length() < 27) {
                cID = curLinePDB.substring(21, 22);
                cTerminusResNumPDB = Integer.valueOf((FileParser.getLineContentsToEndFrom(22, curLinePDB).trim()));                                        
                iCode = ""; // icode parsing is not supported with reduce-made PDB files, I'm not sure what they do with it and how we should discriminate it from the last digit of the residue number
                DP.getInstance().w("Main", " Non-standard PDB 'TER' line of length " + curLinePDB.length() + " encountered. Parsing of iCode not supported for these lines. Assuming terminating residue number '" + cTerminusResNumPDB + "' and empty iCode for chain '" + cID + "'.");
            }
            else {
                // this is the method for standard PDB files
                cTerminusResNumPDB = Integer.valueOf((curLinePDB.substring(22, 26)).trim());
                cID = curLinePDB.substring(21, 22);                        
                iCode = curLinePDB.substring(26, 27);
            }
            
        } catch (Exception e) {
            System.err.println("WARNING: Hit TER line at PDB line number " + curLineNumPDB + " but parsing the line failed: '" + e.getMessage() + "', ignoring.");
            //e.printStackTrace();
            //System.exit(1);
        }
        

        if(! (silent || essentialOutputOnly)) {
            System.out.println("    PDB: Found C Terminus of chain " + cID + " at PDB line " + curLineNumPDB + ", PDB residue number " + cTerminusResNumPDB + " iCode '" + iCode + "'.");
        }
        //TODO: mark the proper residue as C terminus

        return(true);
    }


    /**
     * Does nothing, other lines are ignored atm.
     * @return always returns true
     */
    private static boolean handlePdbLineOTHER() {
        //System.out.println("  Ignored PDB line " + currentLineNumPDB + ", starts with '" + currentLinePDB.substring(0, 7) + "'.");
        return(true);
    }
    
    
    private static void createAllModelsFromHandledPdbLines() {

        Integer pLineNum = 0;
        String mID = "";
        Model m;
        String pLine = "";

        // create default model
        FileParser.s_models.add(new Model(defaultModelName));

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {

            pLineNum = i + 1;
            pLine = FileParser.pdbLines.get(i);

            if(pLine.startsWith("MODEL ")) {

                try {
                    mID = (pLine.substring(10, 16)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Hit MODEL line at PDB line number " + pLineNum + " but parsing the line failed.");
                    e.printStackTrace();
                    System.exit(-1);
                }

                // Model found
                if( ! FileParser.modelExistsWithModelID(mID)) {
                    m = new Model(mID);
                    if(! (silent || essentialOutputOnly)) {
                        System.out.println("    PDB: New PDB Model (model ID '" + mID + "') starts at PDB line " + pLineNum + ".");
                    }
                    FileParser.s_models.add(m);
                }
            }

        }

        // create the default model if this is a non-NMR file (crystal data) that contains no models
        if(FileParser.s_models.size() < 1) {
            if(! (silent || essentialOutputOnly)) {
                System.out.println("    PDB: No models found in handled PDB lines. This most likely is a crystal data (non-NMR) file. Adding default model '" + defaultModelName + "'.");
            }
            FileParser.s_models.add(new Model(defaultModelName));
            
        }

        if(! (silent || essentialOutputOnly)) {
            System.out.println("    PDB: Handled PDB lines contain data from " + FileParser.s_models.size() + " model(s).");
        }

    }


    // Counts all models. Not really needed for the program itself since it only uses the first model, no matter how
    // many others follow.
    private static void createAllModelIDsFromWholePdbFile() {

        ArrayList<String> allPDBLines = FileParser.slurpFile(FileParser.pdbFile);
        Integer pLineNum = 0;
        String mID = "";
        String pLine = "";
        Integer numModels = 0;

        if(! silent) {
            System.out.println("  Counting total number of models in the whole PDB file '" + FileParser.pdbFile + "' (" + allPDBLines.size() + " lines)...");
        }

        for(Integer i = 0; i < allPDBLines.size(); i++) {

            pLineNum = i + 1;
            pLine = allPDBLines.get(i);

            if(pLine.startsWith("MODEL ")) {

                try {
                    mID = (pLine.substring(10, 16)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Hit MODEL line at PDB line number " + pLineNum + " but parsing the line failed.");
                    e.printStackTrace();
                    System.exit(-1);
                }

                // Model found
                //System.out.println("    PDB: New PDB Model (model ID '" + mID + "') starts at PDB line " + pLineNum + ".");
                numModels++;
                FileParser.s_allModelIDsFromWholePDBFile.add(mID);
            }

        }

        if(! silent) {
            System.out.println("    PDB: Scanned whole PDB file for Models, found " + numModels + ".");

            if(numModels < 1) {
                System.out.println("    PDB: This most likely is crystal data (a non-NMR file without models). A default model will be created.");
            }
        }

        if(numModels > 1) {
            if(! silent) {
                System.out.println("    PDB: Multiple models found, all but the default model will be ignored.");
            }

            if(Settings.getBoolean("PTGLgraphComputation_B_split_dsspfile_warning")) {
                System.err.println("*************************************** WARNING ******************************************");
                DP.getInstance().w("Multiple models detected in PDB file. I'm fine with that but unless you did split");
                DP.getInstance().w(" the PDB file into separate models for DSSP, the current DSSP file is broken.");
                DP.getInstance().w(" I parse that file and rely on it. You know the deal: garbage in, garbage out.");
                DP.getInstance().w(" I'll continue but you have been warned. ;)");
                DP.getInstance().w(" (Set 'PTGLgraphComputation_B_split_dsspfile_warning' to 'false' in the config file to suppress this message.)");
                System.err.println("*************************************** WARNING ******************************************");
            }
        }
        
    }
    
    
    private static void setHomologuesOfChains() {

        Integer pLineNum = 0;
        String pLine = "";
        String chainLine = "";
        String[] curChains;
        Boolean nextLineContinuesChainList = false;
          
        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {
          
            pLineNum = i + 1;
            pLine = FileParser.pdbLines.get(i);

                
            if (pLine.startsWith("COMPND  ")) {
                try {
                    if (pLine.substring(11, 17).equals("CHAIN:") || nextLineContinuesChainList) {
                       
                        if (nextLineContinuesChainList) {                      
                            chainLine = chainLine + pLine.substring(11, pLine.length());
                            nextLineContinuesChainList = false;
                        } else {                          
                            chainLine = pLine.substring(18, pLine.length());                           
                        }
                        
                        //when end of chain list is not reached
                        if (!chainLine.contains(";")) {
                            nextLineContinuesChainList = true;
                        } 

                        // do this not before chain list is complete
                        if (!nextLineContinuesChainList) {
                            chainLine = chainLine.replaceAll(" ", "").replace(";", "");
                            curChains = chainLine.split(",");
                            
                            FileParser.fillHomologuesMapFromChainIdList(curChains);
                        }
                    }
                } catch (Exception e) {
                    System.err.println("Hit errenous COMPND line at line #" + pLineNum + "', could not add homologues entry.");
                }
            }             
        }
    }
    
    
    /**
     * Parses binding site data from the PDB data.
     */
    private static void createAllBindingSitesFromPdbData() {
        Integer pLineNum = 0;
        String pLine = "";
        Boolean haveBeenInSiteSectionAlready = Boolean.FALSE;
        String lastSiteInRemark = null;
     
        
        String modelID = null;
        String lastModelID = null;
        String siteName;
        String siteLineNumString;
        Integer siteLineNum;
        BindingSite curSite = null;
        String siteNumResDeclaredString;
        Integer siteNumResDeclared;

        // parse the REMARK lines
        String siteIdentifier = null;
        String siteEvidenceCode;
        String siteDescription;            
        Map<String, BindingSite> sites = new HashMap<>();
            
        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {

            pLineNum = i + 1;
            pLine = FileParser.pdbLines.get(i);
                        
            if(pLine.startsWith("REMARK 800")) {                
                int lineLength = pLine.length();
                if(pLine.startsWith("REMARK 800 SITE_IDENTIFIER:")) {
                    siteIdentifier = pLine.substring(27, (lineLength - 1)).trim();
                    lastSiteInRemark = siteIdentifier;
                    if(! sites.containsKey(siteIdentifier)) { 
                        sites.put(siteIdentifier, new BindingSite(siteIdentifier));
                        //System.out.println("    Adding new binding site '" + siteIdentifier + "' from remark line.");
                    }
                }
                
                if(pLine.startsWith("REMARK 800 EVIDENCE_CODE:")) {
                    siteEvidenceCode = pLine.substring(25, (lineLength - 1)).trim();
                    if(lastSiteInRemark != null) {
                        sites.get(lastSiteInRemark).setEvidenceCode(siteEvidenceCode);
                        if(lastSiteInRemark == null) {
                            DP.getInstance().w("FileParser", "No site encountered before hitting site evidence code line.");
                        }
                    }
                }
                
                if(pLine.startsWith("REMARK 800 SITE_DESCRIPTION:")) {
                    siteDescription = pLine.substring(28, (lineLength - 1)).trim();
                    if(lastSiteInRemark != null) {
                        sites.get(lastSiteInRemark).setDescription(siteDescription);
                        if(lastSiteInRemark == null) {
                            DP.getInstance().w("FileParser", "No site encountered before hitting site description line.");
                        }
                    }
                }
            }
                        
            
            // keep track of current model
            if(pLine.startsWith("MODEL  ")) {
                try {
                    modelID = (pLine.substring(10, 16)).trim();
                    lastModelID = modelID;
                } catch(Exception e) {
                    DP.getInstance().e("FileParser", " Hit broken MODEL line at PDB line number " + pLineNum + " while looking for SITEs.");
                    System.exit(1);
                }
            }
            
            // parse the SITE lines
            if(pLine.startsWith("SITE  ")) {
                
                haveBeenInSiteSectionAlready = true;
                siteLineNumString = pLine.substring(4, 10).trim();
                
                if(siteLineNumString.equals("***")) {
                    DP.getInstance().w("FileParser", " Binding site PDB line is marked with '***', ignoring it. This could be the result of using a non-standard format PDB file, e.g., from the REDUCE software.");
                    continue;
                }
                
                siteLineNum = 0;
                try {
                    siteLineNum = Integer.parseInt(siteLineNumString);
                } catch(Exception e) {
                    DP.getInstance().w("FileParser", "Parsing site internal line num failed for line '" + pLine + "', assuming 0.");
                }
                siteNumResDeclaredString = pLine.substring(14, 17).trim();
                siteNumResDeclared = Integer.parseInt(siteNumResDeclaredString);
                siteName = pLine.substring(11, 14).trim();
                curSite = sites.get(siteName);
                if(curSite == null) {
                    DP.getInstance().w("FileParser", "Site not found in map of sites when hitting SITE line in PDB file, should have been created from data in REMARK section already. Creating it now.");
                    sites.put(siteName, new BindingSite(siteName));
                    curSite = sites.get(siteName);
                }
                
                // set model and number of declared residues on the first line -- we could also do this on any other (or each) line, it would make no difference
                if(siteLineNum.equals(1)) {
                    curSite.setNumResDeclared(siteNumResDeclared);
                    curSite.setModelID(lastModelID);
                }
                
                // parse the residue info from this line and add it
                List<String[]> residueInfos = parseBindingSiteResidueInfos(pLine);
                //System.out.println("Added " + residueInfos.size() + " residues to new site " + siteName + " in site line number " +  siteLineNum + "");
                curSite.addResidueInfos(residueInfos);                                        
            }
        }

        for(String key : sites.keySet()) {
            FileParser.s_sites.add(sites.get(key));
        }
        
    }
    
    /**
     * Parses residue info from a PDB SITE line, which looks like 'SITE     1 AC1 15 ASN A  10  LYS A  12  HIS A  95  GLU A 165                    
SITE     2 AC1 15 ALA A 169  ILE A 170  GLY A 171  SER A 211                    
SITE     3 AC1 15 LEU A 230  GLY A 232  GLY A 233  HOH A 620                    
SITE     4 AC1 15 HOH A 621  HOH A 622  HOH A 623                               
'
     * @param line the input SITE line, like 'SITE     1 AC1 15 ASN A  10  LYS A  12  HIS A  95  GLU A 165'
     * @return a list of residue infos, consisting of residue name, chain, and PDB residue number (e.g., ["ASN", "A", "170"]).
     */
    private static List<String[]> parseBindingSiteResidueInfos(String line) {
        List<String[]> residueInfos = new ArrayList<>();
        
        // parse first res info
        try {
            String[] resInfo = parseResInfo(line.substring(18, 28).trim());
            if(resInfo != null) { residueInfos.add(resInfo); }
        } catch(Exception e) { DP.getInstance().w("FileParser", "parseBindingSiteResidueInfos: First residue in SITE line not found, but should have at least 1 residue."); }
        
        // parse 2nd res info, if available
        try {
            String[] resInfo = parseResInfo(line.substring(29, 39).trim());
            if(resInfo != null) { residueInfos.add(resInfo); }
        } catch(Exception e) { /* It's fine if there is only 1 residue per line */ }
        
        // parse 3rd res info, if available
        try {
            String[] resInfo = parseResInfo(line.substring(40, 50).trim());
            if(resInfo != null) { residueInfos.add(resInfo); }
        } catch(Exception e) { /* It's fine if there is only 1 residue per line */ }
        
        // parse 4th res info, if available
        try {
            String[] resInfo = parseResInfo(line.substring(51, 61).trim());
            if(resInfo != null) { residueInfos.add(resInfo); }
        } catch(Exception e) { /* It's fine if there is only 1 residue per line */ }
        
        
        return residueInfos;
    }
    
    /**
     * Parses residue info from a SITE line substring
     * @param s input string like 'ASN A  10'
     * @return the tokenized info or null if the string did not contain info
     */
    private static String[] parseResInfo(String s) {
        if(s == null) { return null; }
        if(s.isEmpty()) { return null; }
        StringTokenizer t = new StringTokenizer(s, " ");
        String[] data = new String[3];
        if(t.countTokens() == 3) {
            for(int i = 0; i < 3; i++) {                
                data[i] = t.nextToken();
            }
            //System.out.println("Found 3 tokens: " + data[0] + "," + data[1] + "," + data[2]);
            return data;
        }
        else {
            //System.err.println("parseResInfo: Found " + t.countTokens() + " tokens instead of 3.");
        }
        return null;
    }

    private static void createAllChainsFromPdbData() {

        Integer pLineNum = 0;
        String cID = "";
        Chain c, d;
        String pLine = "";
        Boolean cAlreadyExists;
        String modelID = null;

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {

            pLineNum = i + 1;
            pLine = FileParser.pdbLines.get(i);
            cAlreadyExists = false;

            // keep track of current model
            if(pLine.startsWith("MODEL  ")) {
                try {
                    modelID = (pLine.substring(10, 16)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Hit broken MODEL line at PDB line number " + pLineNum + " while looking for Chains.");
                    e.printStackTrace();
                    System.exit(1);
                }
            }

            // parse the chain info from the ATOM and HETATM lines
            if(pLine.startsWith("ATOM  ") || pLine.startsWith("HETATM")) {

                try {
                    cID = pLine.substring(21, 22);
                } catch(Exception e) {
                    System.err.println("ERROR: Could not parse Chain ID colum of PDB line " + pLineNum + ".");
                    e.printStackTrace();
                    System.exit(1);
                }

                if(cID.equals(" ")) {
                    cID = "_";
                    System.out.println("    PDB: Found chain with empty name, renamed to '" + cID + "'.");
                }

                // Chain found, add it if it is a new one
                for(Integer j = 0; j < FileParser.s_chains.size(); j++) {      // s_chains is small so this should be ok
                    d = FileParser.s_chains.get(j);
                    if(d.getPdbChainID().equals(cID)) {
                        cAlreadyExists = true;
                        break;
                    }
                }
                
                if(! cAlreadyExists) {
                    c = new Chain(cID);
                    
                    
                    if(modelID != null) {
                        c.setModelID(modelID);       // left at null otherwise, which is fine
                    }
                    
                    c.setHomologues(FileParser.homologuesMap.get(cID));

                    FileParser.s_chains.add(c);

                    if(! (silent || essentialOutputOnly)) {
                        System.out.println("    PDB: New PDB Chain (chain ID '" + cID + "') starts at PDB line " + pLineNum + ".");
                    }
                }
            }

        }

        if(! silent) {
            System.out.println("    PDB: Scanned PDB file for Chains, found " + FileParser.s_chains.size() + ".");
        }

    }
    
    
    /**
     * Parses the PDB data and creates the ligand list from it.
     */
    private static Integer createAllLigandResiduesFromPdbData() {

        Integer curLigNum = 0;

        Integer atomSerialNumber, resNumPDB, resNumDSSP, lastLigandNumPDB, resIndex;
        atomSerialNumber = resNumPDB = resNumDSSP = lastLigandNumPDB = resIndex = 0;

        String atomName, resNamePDB, molNamePDB, chainID, chemSym, modelID, lastChainID, iCode;
        atomName = molNamePDB =  resNamePDB = chainID = chemSym = modelID = lastChainID = iCode = "";

        String lf, ln, ls;      // temp for lig formula, lig name, lig synonyms
        lf = ln = ls = "";

        modelID = "1";          // default model
        Ligand lig;
        Integer pLineNum = 0;
        String pLine = "";

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {

            pLine = FileParser.pdbLines.get(i);
            pLineNum = i + 1;

            // keep track of current model
            if(pLine.startsWith("MODEL ")) {
                try {
                    modelID = (pLine.substring(10, 16)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Hit broken MODEL line at PDB line number " + pLineNum + " while looking for Chains.");
                    e.printStackTrace();
                    System.exit(1);
                }
            }

            // keep track of last Residue (resNumPDB will be overwritten in next step)
            //lastLigandNumPDB = resNumPDB;

            // we have to create the ligand residues from the HETATM lines
            if(pLine.startsWith("HETATM")) {

                try {

                    atomSerialNumber = Integer.valueOf((pLine.substring(6, 11)).trim());
                    atomName = pLine.substring(12, 16);
                    resNamePDB = pLine.substring(17, 20);
                    chainID = pLine.substring(21, 22);
                    resNumPDB = Integer.valueOf((pLine.substring(22, 26)).trim());
                    iCode = pLine.substring(26, 27);
                    chemSym = pLine.substring(76, 78);
                } catch(Exception e) {
                    System.err.println("ERROR: Hit HETATM line at PDB line number " + pLineNum + " but parsing the line failed (length " + pLine.length() + ").");
                    System.err.println("ERROR: (continued) Message was '" + e.getMessage() + "'.");
                    //e.printStackTrace();
                    System.exit(1);
                    //continue;
                }
                
                // Create new ligand residue here if this Atom does not belong to the same residue as the last one did
                if( ! ( resNumPDB.equals(lastLigandNumPDB) && chainID.equals(lastChainID) ) ) {

                    curLigNum++;

                    // create new Residue from info, we'll have to see whether we really add it below though
                    lig = new Ligand();
                    lig.setPdbNum(resNumPDB);
                    lig.setType(Residue.RESIDUE_TYPE_LIGAND);
                    resNumDSSP = DsspParser.getLastUsedDsspResNumOfDsspFile() + curLigNum; // assign an unused fake DSSP residue number
                    lig.setDsspNum(resNumDSSP);
                    lig.setChainID(chainID);
                    lig.setiCode(iCode);
                    lig.setName3(resNamePDB);
                    lig.setAAName1(AminoAcid.getLigandName1());
                    lig.setChain(FileParser.getChainByPdbChainID(chainID));
                    lig.setModelID(modelID);
                    lig.setSSEString(Settings.get("PTGLgraphComputation_S_ligSSECode"));
                                      

                    // add ligand to list of residues if it not on the ignore list
                    if(FileParser.isIgnoredLigRes(resNamePDB)) {
                        curLigNum--;    // We had to increment before to determine the fake DSSP res number, but
                                        //  this ligand won't be stored so decrement to previous value.
                        //System.out.println("    PDB: Ignored ligand '" + resNamePDB + "-" + resNumPDB + "' at PDB line " + pLineNum + ".");
                    }
                    else {
                        // add info from PDB HET fields (HET, HETNAM, HETSYN, FORMUL)
                        // Note: we now use prepared statements so any strange chars do no longer lead to irritations or security trouble
                        Boolean removeStuff = Settings.getBoolean("PTGLgraphComputation_B_uglySQLhacks");
                        lf = getLigFormula(resNamePDB);
                        if(removeStuff) {
                            lf = lf.replaceAll("\\s", "");               // remove all whitespace
                            lf = lf.replaceAll("~", "");                 // remove tilde char (it causes SQL trouble during DB insert)
                            lf = lf.replaceAll("\\\\", "");              // remove all backslashes (it causes SQL WARNING 'nonstandard use of escape in a string literal' during DB insert)
                            lf = lf.replaceAll("'", "");              // remove all ticks (obviously SQL trouble)
                        }

                        ln = getLigName(resNamePDB);
                        if(removeStuff) {
                            ln = ln.replaceAll(" ", "_");                // replace spaces with underscores
                            ln = ln.replaceAll("\\s", "");               // remove all other whitespace
                            ln = ln.replaceAll("~", "");
                            ln = ln.replaceAll("\\\\", "");
                            ln = ln.replaceAll("'", "");
                        }

                        ls = getLigSynonyms(resNamePDB);
                        if(removeStuff) {
                            ls = ls.replaceAll(" ", "_");
                            ls = ls.replaceAll("\\s", "");
                            ls = ls.replaceAll("~", "");
                            ls = ls.replaceAll(";", ".");
                            ls = ls.replaceAll("\\\\", "");
                            ls = ls.replaceAll("'", "");
                        }

                        lig.setLigFormula(lf);
                        lig.setLigName(ln);
                        lig.setLigSynonyms(ls);

                        lastLigandNumPDB = resNumPDB;
                        lastChainID = chainID;
                        
                        //TODO: Add a check for the molecular weight of the ligand here and only add ligands
                        //      which are within the range (range should be defined in cfg file).
                        //      Problem atm is that the mol weight is not in the PDB file. Idea: count the atoms
                        //      instead, use a range over number of atoms.
                        
                        FileParser.s_molecules.add(lig);
                        
                        FileParser.getChainByPdbChainID(chainID).addMolecule(lig);

                        resIndex = FileParser.s_molecules.size() - 1;
                        FileParser.resIndexDSSP[resNumDSSP] = resIndex;
                        //resIndexPDB[resNumPDB] = resIndex;      // This will crash because some PDB files contain negative residue numbers so fuck it.
                        if(! (silent || essentialOutputOnly)) {
                            System.out.println("    PDB: Added ligand '" +  resNamePDB + "-" + resNumPDB + "', chain " + chainID + " (line " + pLineNum + ", ligand #" + curLigNum + ", DSSP #" + resNumDSSP + ").");
                            System.out.println("    PDB:   => Ligand name = '" + lig.getLigName() + "', formula = '" + lig.getLigFormula() + "', synonyms = '" + lig.getLigSynonyms() + "'.");
                        }
                        
                    }

                }
                else {
                    //System.out.println("    PDB: Atom " + atomSerialNumber + " of reside " + resNumPDB + " belongs to old ligand " +  lastLigandNumPDB + " we already added.");
                }
                
            } // end of: if HETATM
        } // end of: iteration through all PDB lines


        return(curLigNum);
    }
    
    
    // parses the lines of a PDB file to find the formula of a hetero residue in the FORMUL records
    protected static String getLigFormula(String ligName3) {
        String hetID, continuation, asterisk, formula, line;
        hetID = continuation = asterisk = formula = line = "";
        Integer compNum = 0;

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {
            line = FileParser.pdbLines.get(i);

            if(line.startsWith("FORMUL")) {

                // see http://www.wwpdb.org/documentation/format32/sect4.html#FORMUL
                try {

                    compNum = Integer.valueOf((line.substring(8, 10)).trim());
                    hetID = line.substring(12, 15);                         // Can't trim this as spaces are important!
                    continuation = (line.substring(16, 18)).trim();
                    asterisk = (line.substring(18, 19)).trim();
                    
                    // Ugly hack for the non-standard PDB files produced by the 'reduce' hydrogen program follows.
                    // Reduce rewrites the PDB file and does NOT fill the empty lines with spaces, so we need to cope with shorter line lengths here.
                    if(line.length() < 71) {
                        formula = FileParser.getLineContentsToEndFrom(19, line);
                    }
                    else {
                        // this is the method for standard PDB files
                        formula = (line.substring(19, 71)).trim();
                    }
                } catch(NumberFormatException e) {
                    System.err.println("ERROR: Parsing FORMUL line at PDB line number " + i + " failed.");
                    //e.printStackTrace();
                    continue;
                }
                catch(StringIndexOutOfBoundsException s) {
                    System.err.println("ERROR: Parsing FORMUL line at PDB line number " + i + " failed.");
                    continue;
                }

                if(continuation.equals("")) {
                    // We only handle the first line of very long formulas that span multiple lines, the rest
                    //  is truncated. This is a quick hack and should be fixed later. No time atm.
                    if(hetID.equals(ligName3)) {
                        // We found the right line.
                        return(formula);
                    }
                }
            }
            
        }

        // This should never be reached.
        DP.getInstance().w("WARNING: getLigFormula(): No formula found for hetero group " + ligName3 + ".");     
        return("");
    }

    // parses the lines of a PDB file to find the chemical name of a hetero residue in the HETNAM records
    protected static String getLigName(String ligName3) {
        String hetID, continuation, name, line;
        hetID = continuation = name = line = "";

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {
            line = FileParser.pdbLines.get(i);

            if(line.startsWith("HETNAM")) {

                // see http://www.wwpdb.org/documentation/format32/sect4.html#HETNAM
                try {
                    continuation = (line.substring(8, 10)).trim();
                    hetID = line.substring(11, 14);                         // Can't trim this as spaces are important!
                    
                    // Ugly hack for the non-standard PDB files produced by the 'reduce' hydrogen program follows.
                    // Reduce rewrites the PDB file and does NOT fill the empty lines with spaces, so we need to cope with shorter line lengths here.
                    if(line.length() < 71) {
                        name = FileParser.getLineContentsToEndFrom(15, line);
                    }
                    else {
                        name = (line.substring(15, 71)).trim();
                    }
                } catch(Exception e) {
                    System.err.println("ERROR: Parsing HETNAM line at PDB line number " + i + " failed.");
                    //e.printStackTrace();
                    //System.exit(-1);
                    continue;
                }

                if(continuation.equals("")) {
                    // We only handle the first line of very long names that span multiple lines, the rest
                    //  is truncated. This is a quick hack and should be fixed later. No time atm.
                    if(hetID.equals(ligName3)) {
                        // We found the right line.
                        return(name);
                    }
                }
            }

        }


        // Water, deuterated water and methanol have no HETNAM records when used as solvent, so need to warn if none was found for them.
        if(ligName3.equals("HOH") || ligName3.equals("DOD") || ligName3.equals("MOH")) {
            if(ligName3.equals("HOH")) { return("WATER"); }
            if(ligName3.equals("DOD")) { return("DEUTERATED_WATER"); }
            if(ligName3.equals("MOH")) { return("METHANOL"); }
            return("SOLVENT");              // never reached, just for the IDE.
        }
        else {
            // All others should have a HETNAM record so this should never be reached.
            System.out.println("WARNING: getLigName(): No name found for hetero group " + ligName3 + ".");
            return("");
        }
    }

    // parses the lines of a PDB file to find the synonyms for a hetero residue in the HETSYN records
    protected static String getLigSynonyms(String ligName3) {
        String hetID, continuation, synonyms, line;
        hetID = continuation = synonyms = line = "";

        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {
            line = FileParser.pdbLines.get(i);

            if(line.startsWith("HETSYN")) {

                // see http://www.wwpdb.org/documentation/format32/sect4.html#HETSYN
                try {
                    continuation = (line.substring(8, 10)).trim();
                    hetID = line.substring(11, 14);                         // Can't trim this as spaces are important!
                    synonyms = (line.substring(15, 71)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Parsing HETSYN line at PDB line number " + i + " failed.");
                    //e.printStackTrace();
                    //System.exit(-1);
                    continue;
                }

                if(continuation.equals("")) {
                    // We only handle the first line of very long synonyms that span multiple lines, the rest
                    //  is truncated. This is a quick hack and should be fixed later. No time atm.
                    if(hetID.equals(ligName3)) {
                        // We found the right line.
                        return(synonyms);
                    }
                }
            }

        }

        // This could be reached for various valid entries I guess (or do they have a HETSYN line with empty syn field?).
        //System.out.println("WARNING: getLigSynonyms(): No synonym found for hetero group " + ligName3 + ".");
        return("");
    }
    
    
    /**
     * Parses all lines of the PDB file and extracts the source organism of the macromolecular structure with chain ID 'chain'.
     * @param pdbid the PDB ID
     * @param chainid the chain ID
     * @return  the PMI
     */
    protected static ProtMetaInfo getMetaInfo(String pdbid, String chainid) {

        ProtMetaInfo pmi = new ProtMetaInfo(pdbid, chainid);
        String mol_id = pmi.setYourMolID(FileParser.pdbLines);

        if(pmi.isReady()) {
            //System.out.println("    Extracted MOL_ID '" + mol_id + "' for chain '" + chainid + "' from PDB header.");
            if(pmi.parseAllMetaData(FileParser.pdbLines)) {
                //if(! silent) {
                //    System.out.println("    Retrieved all meta data for chain '" + chainid + "' from PDB header.");
                //}
            }
        }
        else {
            DP.getInstance().w("getMetaInfo(): Could not extract MOL_ID for chain '" + chainid + "' from PDB header, meta data unknown.");
        }
            
        return(pmi);
    }
    
    
    /**
     * Retrieves meta data that applies to the whole PDB file from its lines, e.g. experiment type and resolution.
     * Only use this function for the old parser (legacy PDB files)! Use getMetaData() for the new parser (mmCIF).
     * @return A HashMap of (String, String) pairs (key, value) with information on the PDB file. The following
     * Strings are set: 'resolution' (which may be cast to Double), 'experiment', 'keywords', 'header', 'title', 'isLarge' as false.
     */
    protected static HashMap<String, String> getPDBMetaData() {

        HashMap<String, String> md = new HashMap<String, String>();
        Double resolution = -1.0;
        String experiment = "NA";
        String keywords = "";
        String header = "NA";
        String title = "";
        String date = "NA";

        String line;
        for(Integer i = 0; i < FileParser.pdbLines.size(); i++) {
            line = FileParser.pdbLines.get(i);
            
            if(line.startsWith("EXPDTA")) {
                // line looks like this: 'EXPDTA    SOLUTION NMR'
                experiment = line.substring("EXPDATA".length()).trim();
            }
            else if(line.startsWith("REMARK   2 RESOLUTION.")) {
                resolution = getResFromREMARK2Line(line.trim());
            }
            else if(line.startsWith("KEYWDS")) {
                // Multiple KEYWDS lines may exist, e.g.:
                //    KEYWDS    PROTEIN RNA COMPLEX, ANTITERMINATOR COMPLEX, RNA HAIRPIN,
                //    KEYWDS   2 TRANSCRIPTION/RNA COMPLEX
                // but the data part always starts in column 10 (if the first is 0)
                keywords += line.substring(10).trim();
            }
            else if(line.startsWith("HEADER")) {
                try {
                    header = line.substring(10, 50).trim();
                    date = line.substring(50, 59).trim();
                } catch(Exception e) {
                    // Just leave it at "NA" then.
                }
            }
            else if(line.startsWith("TITLE")) {
                // Multiple TITLE lines may exist, but the data part always starts in column 10 (if the first is 0)
                title += line.substring(10).trim();
            }
            else {
                // we aren't interested in other lines atm
            }
        }

        // All data has been parsed and should contain stuff.
        if(keywords.isEmpty()) { keywords = "NA"; }
        if(title.isEmpty()) { title = "NA"; }

        md.put("resolution", resolution.toString());
        md.put("experiment", experiment);
        md.put("keywords", keywords);
        md.put("header", header);
        md.put("title", title);
        md.put("date", date);
        md.put("isLarge", "false");  // only old parser uses this function, therefore always return false (old FP cant treat large structures)
        
        // register md with results
        ProteinResults pr = ProteinResults.getInstance();
        pr.addProteinMetaData("Resolution", resolution.toString());
        pr.addProteinMetaData("Experiment", experiment);
        pr.addProteinMetaData("Keywords", keywords);
        pr.addProteinMetaData("Header", header);
        pr.addProteinMetaData("Title", title);
        pr.addProteinMetaData("Date", date);

        return(md);
    }
    
    
    /**
     * Parses the 'REMARK 2 RESOLUTION' record of a PDB file and returns the resolution in Angstroem, or -1.0 if it is 'NOT APPLICABLE' or could not be determined.
     * Lines may look like this 'REMARK   2 RESOLUTION.    2.24 ANGSTROMS.' or this 'REMARK   2 RESOLUTION. NOT APPLICABLE.'
     */
    private static Double getResFromREMARK2Line(String l) {
        Double res = -1.0;

        String prefix = "REMARK   2 RESOLUTION.";

        // Should have been checked already, but you never know.
        if(l.indexOf(prefix) < 0) {
            return(res);
        }

        try {
            Integer startIndex = l.indexOf(prefix) + prefix.length();       // the column where the data starts (cut off the prefix)
            String data = l.substring(startIndex, l.length());              // the data part

            Integer suffixStart = data.indexOf("ANGSTROMS");
            if(suffixStart >= 0) {
                // There should be a resolution in this data part, it looks similar to this: '    2.24 ANGSTROMS.'
                res = Double.valueOf(data.substring(0, suffixStart).trim());
            }
            else {
                // The data part may look similar to this: ' NOT APPLICABLE.' - or be total rubbish, we don't care.
                //System.out.println("Data part contains no 'ANGSTROMS', assuming 'NOT APPLICABLE'.");
                res = -1.0;
            }

        } catch(Exception e) {
            DP.getInstance().w("Could not parse resolution from 'REMARK 2 RESOLUTION' record, assuming 'NOT APPLICABLE'.");
            res = -1.0;
        }
        return(res);
    }
    
    /**
     * Parses a PDB file with multiple MODELS, converts those models to separate chains, and saves the result in a new PDB file.
     * If no or only one model was found, the method stops because there is nothing to convert.
     * @param pdbFile the input file
     * @param outFile the output path for the new file
     * @return true if converting was successful or if there was nothing to convert (no or only one model), false otherwise.
     */
    protected static Boolean convertPdbModelsToChains(String pdbFile, String outFile) {
        // Read the PDB file
        ArrayList<String> file = FileParser.slurpFile(pdbFile);
        
        Integer pLineNum = 0;
        String mID = "";
        // Contains the all the models as <k,v> = <line number were models starts, model ID>
        TreeMap<Integer, String> models = new TreeMap<Integer, String>();
        String pLine = "";
      
        // Iterate through the PDB file and look at each line
        for(Integer i = 0; i < file.size(); i++) {
            pLineNum = i + 1;
            pLine = file.get(i);
            
            if(pLine.startsWith("MODEL ")) {
                // If the line is marked as MODEL, try to get the model ID
                try {
                    mID = (pLine.substring(10, 16)).trim();
                } catch(Exception e) {
                    System.err.println("ERROR: Hit MODEL line at PDB line number " + pLineNum + " but parsing the line failed.");
                    e.printStackTrace();
                    return false;
                }

                // Model found
                if(!models.containsValue(mID)) {
                    System.out.println("    PDB: New PDB Model (model ID '" + mID + "') starts at PDB line " + pLineNum + ".");
                    models.put(pLineNum, mID);
                }
                else {
                    System.err.println("ERROR: Found models with the same ID.");
                    return false;
                }
            }
        }
        
        // Only if there is more than two models we will have to convert something, otherwise stop here
        if(models.size() < 2) {
            //System.out.println(models.size());
            // Initate a string builder that will store the output string to write the new PDB file
            StringBuilder sb = new StringBuilder();
            String lineSep = System.lineSeparator();
            for(Integer i = 0; i < file.size(); i++) {
                pLine = file.get(i);

                if(!(pLine.startsWith("MODEL ") || pLine.startsWith("ENDMDL"))) {
                    sb.append(pLine);
                    sb.append(lineSep);
                }
            }
            
            // Save the converted PDB file
            File convertedPdbFile = new File(outFile);
            try {
                FileWriter fw = new FileWriter(convertedPdbFile.getAbsoluteFile());
                BufferedWriter bw = new BufferedWriter(fw);
                bw.append(sb);
                bw.close();
            } catch (IOException ex) {
                DP.getInstance().e("FileParser", "Could not write converted PDB file '" + outFile + "': '" + ex.getMessage() + "'.");
                return false;
            }
            
            System.out.println("No different models/only one model found. Nothing to convert here but the MODEL and ENDMDL lines have been removed.");
            return true;
        }
        
        // Initate a string builder that will store the output string to write the new PDB file
        StringBuilder sb = new StringBuilder();
        String lineSep = System.lineSeparator();
        
        // Parse the header of the file until the first MODEL is reached
        for(Integer i = 0; i < file.size(); i++) {
            pLineNum = i + 1;
            pLine = file.get(i);
            
            if(pLine.startsWith("MODEL ")) {
                break;
            }
            else {
                sb.append(pLine);
                sb.append(lineSep);
            }
        }
        
             
        // Alphabet that makes up all possible chain IDs
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        String newChainID = null;
        Integer atomID = 0;
        Integer residueID = 0;
        Integer newAtomID = null;
        Integer terCount = 0;
        
        if(models.size() > alphabet.length()) {
            DP.getInstance().w("FileParser", "There are more models than possible chain IDs (" + alphabet.length() + "). Aborting the conversion now.");
            return false;
        }

        // Go through all the found models by starting at the line the models begin
        for(Integer lineNum : models.keySet()) {
            String chainID = "-1";
            
            // Get the first new chain ID from the alphabet that can be used to rename chains
            if(alphabet.length() > 0) {
                newChainID = alphabet.substring(0, 1);
                alphabet = alphabet.replace(newChainID, "");
            }
            else {
                DP.getInstance().w("FileParser", "Not enough letters to set unique chain IDs. Aborting.");
                return false;
            }
            
            while(!file.get(lineNum).startsWith("ENDMDL")) {
                
                // Count all the TER records in the file. This is later used to write the MASTER record
                if(file.get(lineNum).startsWith("TER")) {
                    terCount++;
                }
                
                // Go through the model until its end is reached and extract the chain ID from each entry
                String newLine = file.get(lineNum);
                chainID = newLine.substring(21,22);
                
                // Get the IDs from the current line which will be the old IDs that get overwritten
                Integer oldAtomID = Integer.parseInt(newLine.substring(6, 11).trim());
                Integer oldresidueID = Integer.parseInt(newLine.substring(22, 26).trim());
                // Now add the highest IDs from the previous model to the current IDs to get the new IDs
                // If we are still processing the first model, zero will be added, so nothing changes
                newAtomID = oldAtomID + atomID;
                Integer newResidueID = oldresidueID + residueID;
                
                if(newAtomID > 99999) {
                    DP.getInstance().w("FileParser", "Trying to assign a new atom ID greater than 99999. The PDB file format does not support this.\n"
                            + "Aborting the conversion now.");
                    return false;
                }
                
                if(newResidueID > 9999) {
                    DP.getInstance().w("FileParser", "Trying to assign a new residue ID greater than 9999. The PDB file format does not support this.\n"
                                        + "Aborting the conversion now.");
                    return false;
                }
                
                // Creates a string of spaces. To keep the PDB file format intact it is necessary to add the correct amount of spaces
                // in front of the later inserted string. To do this, we take the maximum length the entry may have (5 for atom IDs, 4 for residue IDs)
                // and substract the length of the the new ID. With this we know how much spaces need to be added as prefix.
                String spacesAtomID = new String(new char[5 - newAtomID.toString().length()]).replace("\0", " ");
                String spacesResidueID = new String(new char[4 - newResidueID.toString().length()]).replace("\0", " ");
                

                lineNum++;
                
                // Now look at the next line to see what chain ID this line has
                String nextLine = file.get(lineNum);
                String nextChainID = nextLine.substring(21, 22);
                
                // If both chain IDs are different and we have not reached the end of the model yet, 
                // get a new chain ID and replace the old chain ID with this new one
                if(!chainID.equals(nextChainID) && !nextLine.startsWith("ENDMDL")) {
                    newChainID = alphabet.substring(0, 1);
                    // Delete the character from the alphabet as we have used it now
                    alphabet = alphabet.replace(newChainID, "");
                    
                    newLine = new StringBuilder(newLine).replace(21, 22, newChainID).toString();
                    // Add the new IDs
                    newLine = new StringBuilder(newLine).replace(6, 11, spacesAtomID + newAtomID.toString()).toString();
                    newLine = new StringBuilder(newLine).replace(22, 26, spacesResidueID + newResidueID.toString()).toString();
                    sb.append(newLine);
                    sb.append(lineSep);
                }
                else {
                    // If the next line is the same as the one before, we do not have to assing a new chain ID
                    // but we still have to rename the old chain ID with the currently used new chain ID
                    newLine = new StringBuilder(newLine).replace(21, 22, newChainID).toString();
                    // Add the new IDs
                    newLine = new StringBuilder(newLine).replace(6, 11, spacesAtomID + newAtomID.toString()).toString();
                    newLine = new StringBuilder(newLine).replace(22, 26, spacesResidueID + newResidueID.toString()).toString();
                    sb.append(newLine);
                    sb.append(lineSep);
                    
                }
            }
            // Get the last IDs before a new model is reached and the IDs would therefore start all over again
            // We will use those last/highest IDs to add them to the IDs from the next model; in this way the
            // IDs will be consecutive
            atomID = Integer.parseInt(file.get(lineNum - 1).substring(6, 11).trim());
            residueID = Integer.parseInt(file.get(lineNum - 1).substring(22, 26).trim());
        }
        
        // Get the missing spaces for the MASTER record entries
        Integer masterNumAtm = newAtomID - terCount;
        String spacesMasterNumAtm = new String(new char[5 - masterNumAtm.toString().length()]).replace("\0", " ");
        String spacesMasterNumTer = new String(new char[5 - terCount.toString().length()]).replace("\0", " ");
        
        // Find the MASTER record and keep all its entries except the number of ATOM/HETATM records and the number of TER records.
        // Those two entries are updated with the new correct numbers of the converted file. Since all the other numbers to not get
        // changed during the converting, they are not touched.
        for(Integer x = file.size() - 1; x > 0; x--) {
            if(file.get(x).startsWith("MASTER ")) {
                String newMaster = file.get(x);
                newMaster = new StringBuilder(newMaster).replace(50, 55, spacesMasterNumAtm + masterNumAtm.toString()).toString();
                newMaster = new StringBuilder(newMaster).replace(55, 60, spacesMasterNumTer + terCount.toString()).toString();
                sb.append(newMaster);
                sb.append(lineSep);
            }
        }
        
        sb.append("END                                                                             ");
        sb.append(lineSep);
        
        // Save the converted PDB file
        File convertedPdbFile = new File(outFile);
        try {
            FileWriter fw = new FileWriter(convertedPdbFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            bw.append(sb);
            bw.close();
        } catch (IOException ex) {
            DP.getInstance().e("FileParser", "Could not write converted PDB file '" + outFile + "': '" + ex.getMessage() + "'.");
            return false;
        }
        //System.out.println(sb.toString());
        //System.out.println(chains.toString());
        //System.out.println(models.toString());
        return true;
    }
    
    
    /**
     * 
     * Debug function for PTGL compatibility. Parses residue contact data from a file in <PDBID>.geo format (or <PDBID>.geolig if isGeoLig is true). Files in this
     * format can be generated by this program itself or by the PTGL program 'geom_neo.c'.
     * 
     * Call './geom_neo <PDBID>' to generate the input file.
     * 
     * The <pdbid.geo> file format:
     * 
     * All of the following fields are separated by a single space ' '.
     * 
     * 
     * 1      1 2   16 23 12 42  38   0 0    23  2 1    0  0 0    0  0 0    0  0 0   8
     * 1438 413 414  3 33 13 52  37   0 0    24  2 1    0  0 0   32  2 8    0  0 0   8
     * -----------------------------------------------------------------------------------
     * contact number (left justified, length 4)
     *        dssp residue number of 1st residue (right justified, length 3)
     *          dssp residue number of 2nd residue (left justified, length 3)
     *              AA code (as a number, internal geom_neo representation, see num() function) of 1st residue (right justified, lenth 2)
     *                 center sphere radius (around C alpha) of 1st residue, truncated to 99 if it was larger (left justified, lenth 2)
     *                    AA code (as a number, internal geom_neo representation, see num() function) of 2nd residue (right justified, lenth 2)
     *                       center sphere radius (around C alpha) of 1st residue, truncated to 99 if it was larger (left justified, lenth 2)
     *                           C alpha distance of these two residues, set to 999 if it was larger (right justified, length 3)
     * 
     *                                possible H bridge no. 1 contact distance (right justified, length 3)
     *                                  possible H bridge no. 2 contact distance (left justified, length 3)
     * 
     *                                       possible backbone-backbone contact atom distance; 0 if no BB contact (right justified, length 3)
     *                                           atom index* of the backbone atom of 1st residue that forms this BB contact, internal geom_neo number representation (right justified, length 2)
     *                                             atom index* of the backbone atom of the 2nd residue that forms this BB contact, internal geom_neo number representation (left justified, length 2)
     * 
     *                                                  possible chain-backbone contact atom distance; 0 if no CB contact (right justified, length 3)
     *                                                     atom index* of sidechain atom of 1st residue (right justified, length 2)
     *                                                       atom index* of backbone atom of 2nd residue (left justified, length 2)
     * 
     *                                                            possible backbone-chain contact atom distance; 0 if no BC contact (right justified, length 3)
     *                                                               atom index* of backbone atom of 1st residue (right justified, length 2)
     *                                                                 atom index* of sidechain atom of 2nd residue (left justified, length 2)
     * 
     *                                                                      possible chain-chain contact atom distance; 0 if no CC contact (right justified, length 3)
     *                                                                         atom index* of sidechain atom of 1st residue (right justified, length 2)
     *                                                                           atom index* of sidechain atom of 2nd residue (left justified, length 2)
     * 
     *                                                                               number of contacts for this residue pair (right justified, length 2)
     * -----------------------------------------------------------------------------------* 
     * 1      1 2   16 23 12 42  38   0 0    23  2 1    0  0 0    0  0 0    0  0 0   8
     * 1438 413 414  3 33 13 52  37   0 0    24  2 1    0  0 0   32  2 8    0  0 0   8
     * 
     * 
     * NOTE: atom indices of geom_neo start at '1', not at '0'!
     * 
     * 
     * @param filePath the path to the PDBID.geo or .geolig file
     * @param isGeoLig whether the file is in geo (false) or geolig (true) format
     * @param ourContacts an ArrayList containing the contact data computed by this application
     * 
     */
    protected static void compareResContactsWithPdbidDotGeoFile(String filePath, Boolean isGeoLig, ArrayList<MolContactInfo> ourContacts) {
        
        Integer radDif1, radDif2, maxRadDif, distDif, maxDistDif, sumDif, sumDifAbs, common;  // difference logging
        
        Integer expLineLength = 79;     // expected line length (fixed, see .geo format documentation)
        
        if(isGeoLig) {
            
            expLineLength = 143;    // expected line length (fixed, see .geolig format documentation)
            
            // not implemented yet
            System.err.println("ERROR: getResContactInfoFromGeoFile(): Ligand data parsing not implemented yet.");
            System.exit(1);
        }
               
        ArrayList<String> lines = FileParser.slurpFile(filePath);        
        String line;
        
            
        // ------------------------------------ contact arrays init done, init line vars -------------------------------
        Integer undef = -1;
        Integer contNum, resADsspNum, resBDsspNum, resAInternalCode, resBInternalCode, resAColSphereRadius, resBColSphereRadius;
        contNum = resADsspNum = resBDsspNum = resAInternalCode = resBInternalCode = resAColSphereRadius = resBColSphereRadius = undef;
        Integer resDist, hb1Dist, hb2Dist, BBdist, BBatomResAIndex, BBatomResBIndex;
        resDist = hb1Dist = hb2Dist = BBdist = BBatomResAIndex = BBatomResBIndex = undef;
        Integer CBdist, CBatomResAIndex, CBatomResBIndex;
        Integer BCdist, BCatomResAIndex, BCatomResBIndex;
        Integer CCdist, CCatomResAIndex, CCatomResBIndex, totalCont;
        
        CBdist = CBatomResAIndex = CBatomResBIndex = BCdist = BCatomResAIndex = BCatomResBIndex = CCdist = CCatomResAIndex = CCatomResBIndex = totalCont = undef;
        
        
        Integer numParsed = 0;      // number of lines successfully parsed from the pdbid.geo file (geom_neo)
        Integer numIgnored = 0;     // number of lines that could not be parsed from the pdbid.geo file (geom_neo)
        Integer numGeomNeoContactsMissingInPlcc = 0;
        Integer numPlccContactsMissingInGeomNeo = 0;
        Integer numLigContactsIgnored = 0;
        
        radDif1 = radDif2 = maxRadDif = distDif = maxDistDif = sumDif = sumDifAbs = common = 0;
        
        // stuff for comparing the other way around
        Integer numRes = FileParser.s_molecules.size();
        System.out.println("Assuming that " + numRes + " residues exist, preparing matrix.");
        Integer numResDSSP = numRes + 1;       // since DSSP residue numbers start with 1, not 0
        Integer [][] geom_neo_contact_exists = new Integer[numResDSSP][numResDSSP];
        for(Integer i = 0; i < numResDSSP; i++) {
            for(Integer j = 0; j < numResDSSP; j++) {
                geom_neo_contact_exists[i][j] = 0;
            }
        }
        
        
        for(Integer i = 0; i < lines.size(); i++) {
            
            line = lines.get(i);
            
            
            /**
            if(line.length() != expLineLength) {
                DP.getInstance().w("Skipping geo line " + (i + 1) + " of file " + filePath + ": length " + line.length() + ", expected length " + expLineLength + ".");
                numIgnored++;
                continue;
            }
            
            
            // ----------------------------------- parse data -------------------------------
            
            try {                
                // residue properties
                contNum = Integer.parseInt(line.substring(0, 4).trim());
                resADsspNum = Integer.parseInt(line.substring(5, 8).trim());
                resBDsspNum = Integer.parseInt(line.substring(9, 12).trim());
                resAInternalCode = Integer.parseInt(line.substring(13, 15).trim());
                resAColSphereRadius = Integer.parseInt(line.substring(16, 18).trim());
                resBInternalCode = Integer.parseInt(line.substring(19, 21).trim());
                resBColSphereRadius = Integer.parseInt(line.substring(22, 24).trim());
                
                // contacts and distances
                resDist = Integer.parseInt(line.substring(25, 28).trim());
                
                hb1Dist = Integer.parseInt(line.substring(29, 32).trim());
                hb2Dist = Integer.parseInt(line.substring(33, 36).trim());
                
                BBdist = Integer.parseInt(line.substring(37, 40).trim());
                BBatomResAIndex = Integer.parseInt(line.substring(41, 43).trim());
                BBatomResBIndex = Integer.parseInt(line.substring(44, 46).trim());
                
                CBdist = Integer.parseInt(line.substring(47, 50).trim());
                CBatomResAIndex = Integer.parseInt(line.substring(51, 53).trim());
                CBatomResBIndex = Integer.parseInt(line.substring(54, 56).trim());
                
                BCdist = Integer.parseInt(line.substring(57, 60).trim());
                BCatomResAIndex = Integer.parseInt(line.substring(61, 63).trim());
                BCatomResBIndex = Integer.parseInt(line.substring(64, 66).trim());
                
                CCdist = Integer.parseInt(line.substring(67, 70).trim());
                CCatomResAIndex = Integer.parseInt(line.substring(71, 73).trim());
                CCatomResBIndex = Integer.parseInt(line.substring(74, 76).trim());
                
                totalCont = Integer.parseInt(line.substring(77, 79).trim());
                
            } catch(Exception e) {
                System.err.println("ERROR: Hit malformed contact line in residue contact file '" + filePath + "'. Skipping line " + (i + 1) + ".");
                numIgnored++;
                continue;
            }
             */
            
            Scanner scanner = new Scanner(line);
            
            try {                
                // residue properties
                contNum = scanner.nextInt();
                resADsspNum = scanner.nextInt();
                resBDsspNum = scanner.nextInt();
                resAInternalCode = scanner.nextInt();
                resAColSphereRadius = scanner.nextInt();
                resBInternalCode = scanner.nextInt();
                resBColSphereRadius = scanner.nextInt();
                
                // contacts and distances
                resDist = scanner.nextInt();
                
                hb1Dist = scanner.nextInt();
                hb2Dist = scanner.nextInt();
                
                BBdist = scanner.nextInt();
                BBatomResAIndex = scanner.nextInt();
                BBatomResBIndex = scanner.nextInt();
                
                CBdist = scanner.nextInt();
                CBatomResAIndex = scanner.nextInt();
                CBatomResBIndex = scanner.nextInt();
                
                BCdist = scanner.nextInt();
                BCatomResAIndex = scanner.nextInt();
                BCatomResBIndex = scanner.nextInt();
                
                CCdist = scanner.nextInt();
                CCatomResAIndex = scanner.nextInt();
                CCatomResBIndex = scanner.nextInt();
                
                totalCont = scanner.nextInt();
                
            } catch(Exception e) {
                DP.getInstance().w("Hit malformed contact line in residue contact file '" + filePath + "'. Skipping line " + (i + 1) + ".");
                numIgnored++;
                continue;
            }
                    //scanner.nextInt();

            
            
            numParsed++;
            geom_neo_contact_exists[resADsspNum][resBDsspNum] = 1;
            
            // ---------------------------------------------
            
            System.out.println("Comparing external ResContact #" + contNum + ": DSSP pair " + resADsspNum + ", " + resBDsspNum + " with radii " + resAColSphereRadius + "," + resBColSphereRadius +  " in CA dist " + resDist + " ...");
                                    
            // compare with our data, i.e., check the number of geom_neo contacts which were not accepted by PTGLgraphComputation
            
            Boolean found = false;
            radDif1 = radDif2 = distDif = 0;
            for(MolContactInfo rc : ourContacts) {
                
                //System.out.println("Comparing with our contact pair " + rc.getDsspNumA() + "," + rc.getDsspNumB() + ".");
                
                if((rc.getDsspNumA().equals(resADsspNum) && rc.getDsspNumB().equals(resBDsspNum)) || (rc.getDsspNumA().equals(resBDsspNum) && rc.getDsspNumB().equals(resADsspNum))) {                    
                    found = true;
                    common++;
                    
                    // log differences
                    sumDif += rc.getCenterSphereRadiusResA() - resAColSphereRadius;
                    radDif1 = Math.abs(rc.getCenterSphereRadiusResA() - resAColSphereRadius);
                    sumDifAbs += radDif1;
                    if(radDif1 > maxRadDif) { maxRadDif = radDif1; }
                    
                    sumDif += rc.getCenterSphereRadiusResB() - resBColSphereRadius;
                    radDif2 = Math.abs(rc.getCenterSphereRadiusResB() - resBColSphereRadius);
                    sumDifAbs += radDif2;
                    if(radDif2 > maxRadDif) { maxRadDif = radDif2; }
                    
                    sumDif += rc.getMolPairDist() - resDist;
                    distDif = Math.abs(rc.getMolPairDist() - resDist);
                    sumDifAbs += distDif;
                    if(distDif > maxDistDif) { maxDistDif = distDif; }                    
                    
                    System.out.println("  Found. We claim radii " + rc.getCenterSphereRadiusResA() + "," + rc.getCenterSphereRadiusResB() +  " in CA dist " + rc.getMolPairDist() + ". (Differences: " + radDif1 + "," + radDif2 + "," + distDif + ")");
                    break;
                }
                else {
                    // ourContacts are ordered...
                    if(rc.getDsspNumA() > resADsspNum) {
                        // ...so checking the rest is useless.
                        found = false;
                        numGeomNeoContactsMissingInPlcc++;
                        break;
                    }
                    
                }
            }
            if(!found) {
                System.out.println("  NOT found by PTGLgraphComputation: DSSP contact pair " + resADsspNum + ", " + resBDsspNum + " not in PTGLgraphComputation internal contact list.");
            }                                 
        }
        
        // We should also check for contacts which occur in our data but not in the comparison file, i.e.,
        //  the number of PTGLgraphComputation contacts which were not accepted by geom_neo 
        Integer dsspResA, dsspResB;
        for(MolContactInfo rc : ourContacts) {
            dsspResA = rc.getDsspNumA();
            dsspResB = rc.getDsspNumB();
            
            // We need to skip the ligand contacts, of course! They are not computed by geom_neo and they
            // would lead to an array out of bounds exception because they are not included in num_res_DSSP, the 
            // number of DSSP residues that determines the matrix size.
            if(rc.getMolA().isLigand() || rc.getMolB().isLigand()) {
                numLigContactsIgnored++;
                continue;
            }

            if(geom_neo_contact_exists[dsspResA][dsspResB] != 1) {
                numPlccContactsMissingInGeomNeo++;
                System.out.println("  NOT found by geom_neo: DSSP contact pair " + dsspResA + ", " + dsspResB + " not in geom_neo contact list file.");

            }
        }
        
        Integer avgDif, avgDifAbs;
        avgDif = sumDif / common;
        avgDifAbs = sumDifAbs / common;
        
        System.out.println("Parsed " + numParsed + " of " + lines.size() + " lines of file '" + filePath + "', ignored " + numIgnored + " (1st line is always empty so 1 is expected here).");
        System.out.println("File contains " + numParsed + " residue contacts, we found " + (ourContacts.size() - numLigContactsIgnored) + " (and " + numLigContactsIgnored + " ligand contacts).");
        System.out.println("Maximum differences: residue collision sphere radius=" + maxRadDif + ", distance=" + maxDistDif + ". (sumDif=" + sumDif + ", sumDifAbs=" + sumDifAbs + ", avgDif=" + avgDif + ", avgDifAbs=" + avgDifAbs + ")");
        System.out.println("Number of geom_neo contacts which were not detected by PTGLgraphComputation: " + numGeomNeoContactsMissingInPlcc + ".");
        System.out.println("Number of PTGLgraphComputation contacts which were not detected by geom_neo: " + numPlccContactsMissingInGeomNeo + ".");
        System.out.println("Number of ligand contacts ignored in PTGLgraphComputation list: " + numLigContactsIgnored + ".");
    }
    
    
    
    /**
     * Debug function for PTGL compatibility. Parses SSE contact data from a file in geo.dat format. Files in this
     * format can be generated by the PTGL program 'bet_neo.c'.
     * 
     * Call './bet_neo <PDBID><CHAINID>' to generate the input file.
     * 
     * A geo.dat line consists of 9 fields separated by spaces:
     * <pdbid><chain> <SSE#1> <SSE#2> <int1> <int2> <int3> <int4> <spatial_relation> <double_difference>
     * 
     * The fields are:
     * 
     * <pdbid><chain>      : the PDBID and CHAIN of the protein (e.g.: 8icdA)
     * <SSE#1>             : sequential number of the 1st contact SSE
     * <SSE#2>             : sequential number of the 2nd contact SSE
     * <int1>              : number of BB contacts (between residues of both SSEs, counting max 1 contact per residue contact pair)
     * <int2>              : number of CB contacts ...
     * <int3>              : number of BC contacts ...
     * <int4>              : number of CC contacts ...
     * <spatial_relation>  : the spatial relation between the 2 SSEs: p=parallel, a=antiparallel, m=mixed
     * <double_difference> : the double difference value that was computed to determine <spatial_relation>
     * 
     * Example line:
     * 3kmfA 1 4 0 5 0 3 m -4
     * 
     * @param compareSSEContactsFile the input file to parse (the geo.dat output file of geom_neo)
     * @param pg the current protein graph. Its contacts should be compared to those in the input file
     */
    protected static void compareSSEContactsWithGeoDatFile(String compareSSEContactsFile, ProtGraph pg) {
        
        ArrayList<String> lines = FileParser.slurpFile(compareSSEContactsFile);        
        String line;
        
            
        // ------------------------------------ contact arrays init done, init line vars -------------------------------
        Integer undef = -1;
        Integer sseANum, sseAIndex, sseBNum, sseBIndex, numBB, numCB, numBC, numCC, doubleDif;
        String pdbidAndChain, chain, pdbid, spatRel;
        pdbidAndChain = chain = pdbid = spatRel = "?";
        sseANum = sseAIndex = sseBNum = sseBIndex = numBB = numCB = numBC = numCC = doubleDif = -1;
        
        Integer numParsed = 0;      // number of lines successfully parsed from the geo.dat file
        Integer numIgnored = 0;     // number of lines that could not be parsed from the geo.dat file
        Integer numBetNeoContactsMissingInPlcc = 0;
        Integer numPlccContactsMissingInBetNeo = 0;
        Integer numLigContactsIgnored = 0;
        Integer common = 0;
        Boolean found = false;
              
        // stuff for comparing the other way around
        Integer numSSEs = pg.numVertices();
                
        Integer numSSEIDs = numSSEs + 1; // the +1 is required because SSE IDs start with 1 (not 0) in the file
        //System.out.println("DEBUG: Assuming that " + numSSEs + " SSEs exist, preparing matrix.");
        Integer [][] bet_neo_contact_exists = new Integer[numSSEIDs][numSSEIDs];    
        for(Integer i = 0; i < numSSEIDs; i++) {
            for(Integer j = 0; j < numSSEIDs; j++) {
                bet_neo_contact_exists[i][j] = 0;
            }
        }        
        
        for(Integer i = 0; i < lines.size(); i++) {
            
            line = lines.get(i);                               
            found = false;
            
            Scanner scanner = new Scanner(line);
            
            try {                
                // scan SSE properties from the line
                pdbidAndChain = scanner.next();
                pdbid = pdbidAndChain.substring(0, 4);
                chain = pdbidAndChain.substring(4, 5);
                sseANum = scanner.nextInt();
                sseAIndex = sseANum - 1;
                sseBNum = scanner.nextInt();
                sseBIndex = sseBNum - 1;
                numBB = scanner.nextInt();
                numCB = scanner.nextInt();
                numBC = scanner.nextInt();
                numCC = scanner.nextInt();
                spatRel = scanner.next();
                //spatRel = scanner.findInLine("[mapl]");         // 'm' or 'a' or 'p' or 'l' for the spatial SSE contact types
                doubleDif = scanner.nextInt();
                
            } catch(Exception e) {
                DP.getInstance().w("Hit malformed SSE contact line in SSE contact file '" + compareSSEContactsFile + "'. Skipping line " + (i + 1) + ".");
                //DP.getInstance().w("Parsed data (pdbidchain, sseA, sseB, numBB, numCB, numBC, numCC, sr, dd): " + pdbidAndChain + ", " + sseA + ", " + sseB + ", " + numBB + ", " + numBC + ", " + numCB + ", " + numCC + ", " + spatRel + ", " + doubleDif + ".");
                numIgnored++;
                continue;
            }
            
            //System.out.println("INFO:    Parsed data (pdbidchain, sseA, sseB, numBB, numCB, numBC, numCC, sr, dd): " + pdbidAndChain + ", " + sseA + ", " + sseB + ", " + numBB + ", " + numBC + ", " + numCB + ", " + numCC + ", " + spatRel + ", " + doubleDif + ".");
            
            // abort this function if the input file we are parsing is not for the correct chain or PDB ID
            if((! pdbid.equals(pg.getPdbid())) || (! chain.equals(pg.getChainid()))) {
                System.out.println("Skipping SSE comparison for this chain, internal data is for chain '" + pg.getChainid() + "' of pdbid '" + pg.getPdbid() + "', file data for '" + chain + "' of '" + pdbid + "'.");
                return;
            }

            
            
            numParsed++;
            bet_neo_contact_exists[sseANum][sseBNum] = 1;     // the -1 is required because PTGLgraphComputation starts indices at 0, not 1
            //bet_neo_contact_exists[sseB - 1][sseA - 1] = 1;       // Don't do this, they are ordered in there and i<j holds for each contact pair (i, j)
            
            // -----------------------------------------------
            
            if(pg.sseContactExistsPos(sseAIndex, sseBIndex)) {
                found = true;
                common++;
            }
            else {
                found = false;
                System.out.println("  NOT found by PTGLgraphComputation: potential SSE contact num pair (" + sseANum + ", " + sseBNum + ") not in internal accepted contact list.");
                numBetNeoContactsMissingInPlcc++;                
            }
        }
        
        // check the other way around
        
        for(Integer i = 0; i < pg.numVertices(); i++) {         // i and j are now SSE indices, not SSE numbers
            for(Integer j = i + 1; j < pg.numVertices(); j++) {
                
                sseAIndex = i;
                sseBIndex = j;
                sseANum = i + 1;
                sseBNum = j + 1;
                
                if(pg.getSSEBySeqPosition(sseAIndex).isLigandSSE() || pg.getSSEBySeqPosition(sseAIndex).isOtherSSE() || pg.getSSEBySeqPosition(sseBIndex).isLigandSSE() || pg.getSSEBySeqPosition(sseBIndex).isOtherSSE()) {
                    numLigContactsIgnored++;
                    continue;
                }
                
                // If PTGLgraphComputation list a contact here...
                if(pg.sseContactExistsPos(sseAIndex, sseBIndex)) {
                    
                    // ...check whether it also exists in bet_neo
                    if(bet_neo_contact_exists[sseANum][sseBNum] != 1) {
                        numPlccContactsMissingInBetNeo++;
                        System.out.println("  NOT found by bet_neo: SSE contact number pair (" + sseANum + ", " + sseBNum + ") not in bet_neo contact list file.");
                    }
                }                
            }
        }
        
        // print results
        System.out.println("Parsed " + numParsed + " of " + lines.size() + " lines of file '" + compareSSEContactsFile + "', ignored " + numIgnored + " (1st line is always empty so 1 is expected here).");
        System.out.println("File contains " + numParsed + " potential SSE contacts, we found " + (pg.numSSEContacts() - numLigContactsIgnored) + " accepted contacts (and " + numLigContactsIgnored + " ligand contacts).");
        System.out.println("Note that the geom_neo list contains all potential SSE contacts (not filtered by rule set yet), while the PTGLgraphComputation list contains only accepted contacts.");
        System.out.println("Number of bet_neo contacts which were not detected by PTGLgraphComputation: " + numBetNeoContactsMissingInPlcc + " (but read above).");
        System.out.println("Number of PTGLgraphComputation contacts which were not detected by bet_neo: " + numPlccContactsMissingInBetNeo + ".");
        System.out.println("Number of ligand contacts ignored in PTGLgraphComputation list: " + numLigContactsIgnored + ".");                
    }
     
     
    
}
    