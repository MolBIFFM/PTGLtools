/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;


// imports
import proteinstructure.Model;
import proteinstructure.Residue;
import proteinstructure.Chain;
import proteinstructure.AminoAcid;
import proteinstructure.Atom;
import proteinstructure.SSE;
import proteinstructure.Molecule;
import tools.DP;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.*;
import java.util.Arrays;
import java.util.List;
import settings.Settings;
import proteingraphs.MolContactInfo;
import proteingraphs.ProtGraph;
import proteinstructure.BindingSite;
import proteinstructure.ProtMetaInfo;


/**
 * A file parser, used to extract information from PDB and DSSP files.
 * 
 * @author ts
 */
public class FileParser {

    // declare class vars
    public static Boolean silent = false;
    public static Boolean essentialOutputOnly = false;
    static String fs = System.getProperty("file.separator");

    static String pdbFile = null;
    static ArrayList<String> pdbLines = null;
    static ArrayList <Molecule> s_molecules = null;
    static ArrayList<Model> s_models = null;
    static ArrayList<String> s_allModelIDsFromWholePDBFile = null;
    static ArrayList<Chain> s_chains = null;
    static ArrayList<Integer> s_residueIndices = null;
    static ArrayList<Integer> s_rnaIndices = null;
    static ArrayList<Integer> s_ligandIndices = null; // a list of all non-ignored ligand indices
    static Integer ignoredLigands = 0;
    static ArrayList<Atom> s_atoms = null;
    static ArrayList<SSE> s_ptglSSEs = null;                // the modified SSE list the PTGL uses
    static HashMap<String, ArrayList<String>> homologuesMap = null;
    static List<BindingSite> s_sites;
    
    static Integer lastIndexGetRes = 0;  // last index from s_residue where getResidueFromList found target
    
    // The list of sulfur bridges (aka disulfide bridges) from the DSSP file. The key is the DSSP sulfur bridge
    // id (an arbitrary character, starting with 'a' for the first bridge usually). The list in the value part contains
    // the DSSP residue IDs of all residues which are part of the bridge (and thus should have length 2).
    static HashMap<Character, ArrayList<Integer>> s_sulfurBridges;
    
    // The list of interchain sulfur bridges from the DSSP file. The list structure is the same as for s_sulfurBridges.
    static HashMap<Character, ArrayList<Integer>> s_interchainSulfurBridges;
    // This list is needed to track if sulfur bridges span between to independent chains.
    static HashMap<Character, String> s_interchainSulfurBridgesChainID;
    
    protected static Boolean dataInitDone = false;

    static Integer[] resIndexPDB = null;        // for a DSSP res number, holds the index of that residue in s_residues
    static Integer[] resIndexDSSP = null;       // for a PDB res number, holds the index of that residue in s_residues


      
    /**
     * Parses the contents of a line from the given startIndex to the end of the line.
     * @param startIndexInclusive the start index, inclusive
     * @param line the whole line as a string
     * @return the contents of the line from start index to end -- may be empty. if the line is too short, the string is also empty (even if it is shorter than start index).
     */
    public static String getLineContentsToEndFrom(int startIndexInclusive, String line) {
        if(line == null) {
            return "";
        }
        int length = line.length();
        if(length <= startIndexInclusive) {
            return "";
        }
        String substr = line.substring(startIndexInclusive, length);
        return substr;
    }
    
    
    /**
     * Inits class variables.
     * @param pf PDB file path
     */
    protected static void initVariables(String pf) {
        silent = Settings.getBoolean("PTGLgraphComputation_B_silent");
        
        pdbFile = pf;
        
        // read all lines of the files into lists
        if(! FileParser.silent) {
            System.out.println("  Reading PDB file...");
        }
        pdbLines = new ArrayList<String>();
        pdbLines = slurpPDBFileToModel(pdbFile, "2");
        if(! FileParser.silent) {
            System.out.println("    Read " + pdbLines.size() + " lines of file '" + pdbFile + "'.");
        }

    
        s_models = new ArrayList<Model>();
        s_allModelIDsFromWholePDBFile = new ArrayList<String>();
        s_chains = new ArrayList<Chain>();
        s_molecules = new ArrayList <Molecule> ();
        s_residueIndices = new ArrayList<>();
        s_rnaIndices = new ArrayList<>();
        s_ligandIndices = new ArrayList<>();
        s_atoms = new ArrayList<Atom>();
        s_ptglSSEs = new ArrayList<SSE>();
        s_sites = new ArrayList<>();
        s_sulfurBridges = new HashMap<Character, ArrayList<Integer>>();
        s_interchainSulfurBridges = new HashMap<Character, ArrayList<Integer>>();
        s_interchainSulfurBridgesChainID = new  HashMap<Character, String>();
        homologuesMap = new HashMap<>();
    }

    
    /**
     * Reads the target text file and returns the data in it.
     * @param file Path to a readable text file. Does NOT test whether it exist, do that earlier.
     * @param bool_dssp boolean vararg. If given: expects one bool which determines wether file is a dssp file.
     * @return all lines of the file as an ArrayList
     */
    public static ArrayList<String> slurpFile(String file, boolean... bool_dssp) {
        
        boolean local_warning = false; // ensures that warning is only printed once
        boolean local_skip_blanks = false;
        
        if (bool_dssp.length > 0) {
            local_skip_blanks = bool_dssp[0];
        }

        ArrayList<String> lines = new ArrayList<String>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = null;
            while ((line = in.readLine()) != null) {
                if (! local_skip_blanks) {
                    lines.add(line);
                }
                else {
                    if (! line.isEmpty()) { // skips empty lines and prints warning
                        lines.add(line);
                    }
                    else {
                        if (! local_warning) {
                            DP.getInstance().w("FP", "DSSP File contains empty lines. Parser skips them and tries to go on.");
                            DP.getInstance().w("FP", " ... You should check if the empty line occurs at the end of the file "
                                    + "(b/c of your DSSP-database or text editor you pasted to)");
                            DP.getInstance().w("FP", " ... or if this indicates a severe error. Garbage in -> garbage out.");
                            local_warning = true;
                        }
                    }
                }
            }
	} catch (IOException e) {
            System.err.println("ERROR: Could not read text file '" + file + "'.");
            e.printStackTrace();
            System.exit(1);
	}

        return(lines);
    }
    
    
    /**
     * Reads the target text file and returns the data in it.
     * @param file Path to a readable text file. Does NOT test whether it exist, do that earlier.
     * @return all lines of the file as a single string
     */
    public static String slurpFileToString(String file) throws IOException {

        String lines = "";

        //try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = null;
            while ((line = in.readLine()) != null) {
                lines += line + "\n";
            }
	//} 
        /*catch (IOException e) {
            System.err.println("ERROR: Could not read text file '" + file + "':" + e.getMessage() + ".");
            //e.printStackTrace();
            System.exit(1);
	}
        */

        return(lines);
    }
    
    
    /**
     * Reads the target text file and returns the data in it.
     * @param file Path to a readable text file. Does NOT test whether it exist, do that earlier.
     * @return the data as a multi line string 
     */
    public static String slurpFileSingString(String file) {

        String lines = "";
        Integer numLines = 0;

        try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = null;
            while ((line = in.readLine()) != null) {
                numLines++;
                lines += line + "\n";
            }
	} catch (IOException e) {
            System.err.println("ERROR: Could not read text file '" + file + "'.");
            e.printStackTrace();
            System.exit(1);
	}

        return(lines);
    }


    // Reads a PDB file until it hit the model with the ID specified in the argument 'stopAtModelID'.
    // This means that if you set the latter to "3" it will stop reading after the last line of model "2", ignoring
    // the rest of the file.
    protected static ArrayList<String> slurpPDBFileToModel(String file, String stopAtModelID) {

        ArrayList<String> lines = new ArrayList<String>();
        Integer numLines = 0;
        String mID = "";

        Boolean error = false;
        Boolean endModelAlreadyHit = false;
        String errMsg = "";

        try {
            BufferedReader in = new BufferedReader(new FileReader(file));
            String line = null;
            while ((line = in.readLine()) != null) {

                // make sure we abort before the model
                if(line.startsWith("MODEL ")) {

                    try {
                        mID = (line.substring(10, 16)).trim();
                    } catch(Exception e) {
                        System.err.println("ERROR: Hit MODEL line at PDB line number " + (numLines + 1) + " but parsing the line failed.");
                        error = true;
                        errMsg = "ERROR: Hit MODEL line at PDB line number " + (numLines + 1) + " but parsing the line failed.";
                        e.printStackTrace();
                        System.exit(-1);
                    }

                    if( mID.equals(stopAtModelID) || endModelAlreadyHit ) {
                        endModelAlreadyHit = true;
                        if(! FileParser.silent) {
                            System.out.println("  Aborting the parsing of PDB file at line " + (numLines + 1) + " as requested because model " + mID + " starts here.");
                        }
                        break;
                    }
                    else {

                        if(! error) {
                            numLines++;
                            lines.add(line);
                        }
                    }

                }
                else {
                    // just add all other lines
                    numLines++;
                    lines.add(line);
                }
                
            }
	} catch (IOException e) {
            System.err.println("ERROR: Could not read text file '" + file + "'.");
            if(error) {
                System.err.println("ERROR: Parsing line failed: '" + errMsg + "'.");
            }
            e.printStackTrace();
            System.exit(1);
	}

        return(lines);
    }


    // getters
    public static ArrayList<Model> getModels() { return(s_models); }
    public static ArrayList<String> getAllModelIDsFromWholePdbFile() { return(s_allModelIDsFromWholePDBFile); }
    public static ArrayList<Chain> getChains() { return(s_chains); }
    public static ArrayList<Molecule>getMolecule(){ return(s_molecules);}
    public static ArrayList<Integer> getResidueIndices() { return(s_residueIndices); }
    public static ArrayList<Integer> getRnaIndices() {return (s_rnaIndices);}
    public static ArrayList<Atom> getAtoms() { return(s_atoms); }
    public static ArrayList<SSE> getPtglSSEs() { return(s_ptglSSEs); }


    public static HashMap<Character, ArrayList<Integer>> getSulfurBridges() {
        return s_sulfurBridges;
    }

    public static HashMap<Character, ArrayList<Integer>> getInterchainSulfurBridges() {
        return s_interchainSulfurBridges;
    }
    
     
    /**
    * Inserts spaces on left end until targetStr has length newLenght.
     * If newLength less than / equals length of targetStr then targetStr is returned
     * @param targetStr
     * @param newLength 
     * @return 
     */
     public static String leftInsertSpaces(String targetStr, int newLength) {
         if (newLength <= targetStr.length()) {
             return targetStr;
         }
         return new String(new char[newLength - targetStr.length()]).replace("\0", " ") + targetStr;
     }
    
     
    /**
     * Returns true if the residue name stands for a 
     * @param resNamePDB
     * @return whether the residue name marks a DNA / RNA residue
     */
    public static boolean isRNAresidueName(String molNamePDB) {
        
        // standard ribonucleotides
        ArrayList<String> RNANames = new ArrayList<>(Arrays.asList("  G", "  U", "  A", "  T", "  C", "  I", "G  ", "U  ", "A  ", "T  ", "C  ", "I  ", "G", "U", "A", "T", "C", "I"));
        if (RNANames.contains(molNamePDB)){
            return true;
        }
        
        return false;
    }
    
    
    public static boolean isDNAresidueName(String molNamePDB) {
        
        // standard deoxribonucleotides
        if(molNamePDB.equals(" DG") || molNamePDB.equals(" DU") || molNamePDB.equals(" DA") || molNamePDB.equals(" DT") || molNamePDB.equals(" DC") || molNamePDB.equals(" DI")) {
            return true;
        }
        
        
        return false;
    }
    
    
    /**
     * Returns true if AAName is standard aminoacid name (3-letter code).
     * @param AAName Amino acid name, 3-letter code, capitalized
     * @param includeNonStandard allows additional to 20 standard amino acids: UNK
     * @return 
     */
    protected static boolean isAminoacidName(String AAName, Boolean includeNonStandard) {
        ArrayList<String> AANames = new ArrayList<>(Arrays.asList("ALA", "ARG", "ASN", "ASP", "CYS", 
            "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", 
            "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"));
        
        if (includeNonStandard) {
            AANames.add("UNK");  // UNKnown, 'X' or any amino acid (often only backbone atoms are present
            AANames.add("MSE");  // Selenomethionine: Methionine with selen instead of sulfur (helps in crystallography)
        }
        
        if (AANames.contains(AAName)) {
            return true;
        } else {
            return false;
        }
    }
    
    
    /**
     * Returns true if chainIdentity equals polyribonucleotide.
     * @param chainID Name of the chain i.e. A
     * @return 
     */
    protected static boolean isRNAchain(String chainID) {
        return (CifParser.chainIdentity.get(chainID).equals("polyribonucleotide"));
    }
    

    /**
     * Calls getResidueFromList and returns its value with printing of error message if null returned.
     * Only residues, no ligands or RNA can be found.
     * @param resNumPDB
     * @param chainID
     * @param iCode
     * @param atomSerialNumber
     * @param numLine current line of pdb file.
     * @return 
     */
    //ändern zu getMolFromListWithErrMsg, return tmpMol
    private static Residue getResFromListWithErrMsg(Integer resNumPDB, String chainID, String iCode, Integer atomSerialNumber, Integer numLine) {
        Residue tmpRes = getResidueFromList(resNumPDB, chainID, iCode);
        if (tmpRes == null) {
            DP.getInstance().w("FP_CIF", " Residue with PDB # " + resNumPDB + 
                    " of chain '" + chainID + "' with iCode '" + iCode + 
                    "' not listed in DSSP data, skipping atom " + atomSerialNumber.toString() + 
                    " belonging to that residue (PDB line " + numLine.toString() + ").");
        }
        return tmpRes;
    }

    
    /**
     * Tries to get the residue with the given PDB residue number, chain ID and insertion code from the internal list of all residues.
     * Always starts at index of last hit (beginning with 0), i.e. rotates through list.
     * Only residues, no ligands or RNA can be found.
     * @param resNumPDB the PDB residue number
     * @param chainID the chain ID of the residue
     * @param iCode the insertion code of the residue
     * @return the residue if such a residue exists, null if no such residue exists.
     */
    protected static Residue getResidueFromList(Integer resNumPDB, String chainID, String iCode) {
       
        Residue tmp;
        Residue found = null;

        // iterate up to s_molecules.size() times
        for(Integer i = 0; i < s_molecules.size(); i++) {
            
            // start at last occurence
            Integer currentIndex = (lastIndexGetRes + i) % s_molecules.size();
            
            if (s_molecules.get(currentIndex) instanceof Residue) {
                tmp = (Residue) s_molecules.get(currentIndex);
                              
                if(tmp.getPdbNum().equals(resNumPDB)) {

                    if(tmp.getChainID().equals(chainID)) {

                        if(tmp.getiCode().equals(iCode)) {
                            found = tmp;
                            lastIndexGetRes = currentIndex;

                            // break here and return found to increase speed
                            return(found);
                        }
                    }
                }
            }
        }

        // Not found in the whole list, something went wrong
        //System.err.println("ERROR: Residue with PDB residue number '" + resNum + "' and chain ID '" + cID + "' and iCode '" + ic + "' does not exist in residue list.");
        return(found);
    } 
    

    /**
     * Finds the residue with PDB residue number 'p' of chain 'chID' with iCode 'ic' in the residue list.
     * Only residues, no ligands or RNA can be found.
     * @param p pdb residue number
     * @param chID chain id
     * @param ic insert code, may be null if you don't care
     * @return The Residue object if found, null otherwise.
     */
    public static Residue getResByPdbFields(Integer p, String chID, String ic) {
        Residue r;
        
        for(Integer i = 0; i < s_molecules.size(); i++) {
            if (s_molecules.get(i) instanceof Residue) {
                r = (Residue) s_molecules.get(i);
                if((r.getPdbNum()).equals(p)) {

                    if((r.getChainID()).equals(chID)) {

                        if( ic == null || (r.getiCode()).equals(ic)) {
                            return(r);
                        }
                    }
                }
            }
        }

        // only reached if not found
        // System.out.println("WARNING: Could not find Residue with PDB number " + p + " and chain ID '" + chID + "' and iCode '" + ic + "'.");
        return(null);
    }


    protected static Boolean isIgnoredLigRes(String r) {
        if(r.equals("DOD") || r.equals("HOH")) {
            return(true);
        } else {
            return(false);
        }
    }


    // determines whether a model with a certain model ID already exists in s_models
    protected static Boolean modelExistsWithModelID(String mID) {

        Model m;

        for(Integer i = 0; i < s_models.size(); i++) {

            m = s_models.get(i);

            if(m.getModelID().equals(mID)) {
                return(true);
            }
        }

        return false;
    }


    public static Chain getChainByPdbChainID(String cID) {
        Chain resultChain = null;
        Chain tmpChain = null;

        for(Integer i = 0; i < s_chains.size(); i++) {
            tmpChain = s_chains.get(i);

            if(cID.equals(tmpChain.getPdbChainID())) {
                resultChain = tmpChain;
                break;      // no need to look further
            }
        }

        if(resultChain == null) {
            DP.getInstance().e("FP", "Could not find Chain '" + cID + "' in s_chains. Maybe using cif-type "
                    + "DSSP file instead of legacy type? Exiting now.");
            System.exit(1);
            return(null);   // for the IDE
        }

        return(resultChain);
    }

    
    public static Model getModelByModelID(String mID) {
        Model resultModel = null;
        Model tmpModel = null;

        for(Integer i = 0; i < s_models.size(); i++) {
            tmpModel = s_models.get(i);

            if(mID.equals(tmpModel.getModelID())) {
                resultModel = tmpModel;
                break;      // no need to look further
            }
        }

        if(resultModel == null) {
            System.err.println("ERROR: Could not find Model '" + mID + "' in s_models.");
            System.exit(-1);
            return(null);   // for the IDE
        }

        return(resultModel);
    }


    // Determines if an ATOM should be ignored via its chemSym string
    protected static Boolean isIgnoredAtom(String cSym) {

        // Ignored atoms: H=hydrogen, Q=pseudo atom used in NMR refinement, D=deuterium (^2H)
        if(cSym.trim().equals("H") || cSym.trim().equals("Q") || cSym.trim().equals("D")) {
            return(true);
        }

        return(false);
    }

    public static Boolean convertPdbModelsToChains(String convertModelsToChainsInputFile, String convertModelsToChainsOutputFile) {
        if (settingCif()) {
            DP.getInstance().w("FP", "Converting mmCIF file's models to chains is not implemented. "
                    + "Turn use of mmCIF parser off and retry with legacy PDB file.");
            return false;
        } else {
            return LegacyParser.convertPdbModelsToChains(convertModelsToChainsInputFile, convertModelsToChainsOutputFile);
        }
    }
    
    public static void initData(String pdbFile, String dsspFile, String outputDir) {
        DsspParser.initVariables(dsspFile);
        FileParser.initVariables(pdbFile);
        
        if (settingCif()) {
            CifParser.initData(pdbFile);
        } else {
            DP.getInstance().w("Legacy-Parser chosen: This parser is no longer kept up to date, results might not be correct anymore. If possible, use Cif-Parser in the future.");
            LegacyParser.initData(pdbFile);
        }
        
        if (Settings.getBoolean("PTGLgraphComputation_B_csv_number_residues_chains")) {
            File numberResInChains = new File(outputDir +  "number_of_residues_in_each_chain.csv");
            if (!numberResInChains.exists()) {
                try {
                    numberResInChains.createNewFile();
                } catch (IOException ex) {
                    DP.getInstance().e("ERROR: Could not create file '" + numberResInChains.getAbsolutePath() + ".");
                }
            }
            String resInChain = "Chain_ID" + "," + "chain length" + "\n";
            FileParser.getChains();
            for(Chain c: s_chains){               
                resInChain += c.getPdbChainID() + "," + c.getChainLength() + "\n";
            }
            IO.writeStringToFile(resInChain, numberResInChains.getAbsolutePath(), true);  
            DP.getInstance().w("Csv file with the number of residues in each chain created.");
            }
        }
    
    
    public static void compareResContactsWithPdbidDotGeoFile(String compareResContactsFile, boolean b, ArrayList<MolContactInfo> cInfo) {
        LegacyParser.compareResContactsWithPdbidDotGeoFile(compareResContactsFile, b, cInfo);
    }
    
    public static void compareSSEContactsWithGeoDatFile(String get, ProtGraph pg) {
        LegacyParser.compareSSEContactsWithGeoDatFile(pdbFile, pg);
    }
    
    
    /**
     * Fills homologuesMap from list of chain IDs.
     * @param ChainIDs
     */
    protected static void fillHomologuesMapFromChainIdList(String[] ChainIDs) {
        for (String chain : ChainIDs) {
            ArrayList<String> homologueChains = new ArrayList<>();
            for (String hChain : ChainIDs) {
                if (!chain.equals(hChain)) {
                    homologueChains.add(hChain);
                }
            }
            homologuesMap.put(chain, homologueChains);
        }
    }
    
    
    public static HashMap<String, String> getMetaData() {
        if (settingCif()) {
            // meta Data in CIF files is created while parsing the file once
            //     -> no need to call a function to create it, just get it!
            return CifParser.metaData;
        } else {
            return LegacyParser.getPDBMetaData();
        }
    }
    
    public static ProtMetaInfo getMetaInfo(String pdbid, String chain) {
        if (settingCif()) {
            return CifParser.getProteinMetaInfo(pdbid, chain);
        } else {
            return LegacyParser.getMetaInfo(pdbid, chain);
        }
    }

    // setting methods
    private static boolean settingCif() {
        return Settings.getBoolean("PTGLgraphComputation_B_use_mmCIF_parser");
    }
    
    // protected so that other parsers can use
    protected static boolean settingSilent() {
        return Settings.getBoolean("PTGLgraphComputation_B_silent");
    }
    
    // protected so that other parsers can use
    protected static boolean settingEssentialOutput() {
        return Settings.getBoolean("PTGLgraphComputation_B_only_essential_output");
    }


}