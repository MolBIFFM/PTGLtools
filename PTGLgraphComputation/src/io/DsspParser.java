/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;

// imports
import java.util.ArrayList;
import java.util.HashMap;
import settings.Settings;
import proteinstructure.AminoAcid;
import proteinstructure.Residue;
import proteinstructure.SSE;
import tools.DP;



/**
 *
 * @author niclas
 */
public class DsspParser {
    
    // declare class vars
    static String dsspFile = null;
    static Integer dsspDataStartLine = null;
    static ArrayList<String> dsspLines = null;
    static ArrayList<SSE> s_dsspSSEs = null;             // all SSEs according to DSSP definition
    
    static Integer curLineNumDSSP = null;
    static String curLineDSSP = null;
    static Integer curResNumDSSP = null;
    static Integer lastUsedDsspNum = null;
    static Integer maxUsedDsspResNumInDsspFile = null;  // used to determine fake DSSP numbers for ligands
    
    static boolean dataInitDone = false;
    static boolean silent;
    public static ArrayList<Integer> chainBreakDssps = new ArrayList<Integer>();
    
    
    /**
     * Inits class variables.
     * @param df DSSP file path
     */
    protected static void initVariables(String df) {
        dsspFile = df;
        silent = Settings.getBoolean("PTGLgraphComputation_B_silent");
        
        // read all lines of the files into lists
        if(! silent) {
            System.out.println("  Reading DSSP file...");
        }
        
        dsspLines = new ArrayList<String>();
        dsspLines = FileParser.slurpFile(dsspFile, true); // vararg tells the function that this is a dssp file
        if(! silent) {
            System.out.println("    Read all " + dsspLines.size() + " lines of file '" + dsspFile + "'.");
        }
        
        dsspDataStartLine = readDsspToData();

        s_dsspSSEs = new ArrayList<SSE>();
    }
    

    // reads the dssp lines (starting at first line) till it arrives at the first data line (skips header stuff)
    private static Integer readDsspToData() {

        curLineNumDSSP = 0;
        curLineDSSP = "";
        Boolean hitDsspData = false;

        for(Integer i = 0; i < dsspLines.size(); i++) {
            curLineNumDSSP++;
            curLineDSSP = dsspLines.get(i);

            if(curLineDSSP.startsWith("  #")) {     // We found the last line of the header. The next line is the first one we are interested in.

                if(curLineNumDSSP >= dsspLines.size() - 1) {        // this may already be the last line if the DSSP file is broken, we shouldn't go to the next line in that case ;)
                    System.err.println("ERROR: DSSP file '" + dsspFile + "' ends after last header line (line " + curLineNumDSSP + ".");
                    System.exit(1);
                }
                else {                                          // yay, data found!
                    curLineNumDSSP++;
                    curLineDSSP = dsspLines.get(i + 1);
                    hitDsspData = true;                    
                    break;                                      // look no further!
                }
            }
        }

        if(hitDsspData) {
            if(! silent) {
                System.out.println("    DSSP: Found start of DSSP data in line " + curLineNumDSSP + ".");
            }
            return(curLineNumDSSP);
        }
        else {
            System.err.println("ERROR: Hit the end of DSSP file at line " + curLineNumDSSP + " without finding the end of the header.");
            System.err.println("ERROR: DSSP file contains no data. Maybe the PDB file contains only DNA/RNA, check COMPND lines. Exiting.");
            System.exit(2);
            return(null);
        }

    }
    
    
    public static Integer getDsspResNumForPdbFields(Integer prn, String chainID, String iCode) {
        Integer foundDsspResNum = null;
        Integer foundPdbResNum = null;
        String foundPdbICode = null;
        String foundChain = null;
        Integer resultDsspResNum = null;
        String dline = null;
        String tmpPdbResNum = null;
        Character lastChar = null;

        //System.out.println("    Starting search at DSSP line number '" + (dsspDataStartLine - 1) + "'.");
        for(Integer i = dsspDataStartLine - 1; i < dsspLines.size(); i++) {
            dline = dsspLines.get(i);
            foundDsspResNum = getDsspNumFromLine(dline);

            //System.out.println("Looking for PDB# " + prn + " at line "  + (i + 1) + ": '" + dline + "'");

            // skip chain brakes
            if(isChainBreakLine(dline)) {
                continue;
            }
            
            // The PDB residue number string may still contain the iCode (e.g. '123A', thus we cannot simply try to cast it to Integer
            try{
                tmpPdbResNum = dline.substring(7, 11).trim();
//                foundChain = dline.substring(11, 12).trim();
                foundChain = dline.substring(159, 163).trim(); // AUTHCHAIN column 160-163. md: we now create the residues from cif file and use auth-info. Therefore we need the authchain here

                if(tmpPdbResNum.length() == 0) { // The next command (lastChar) would be out of bounds anyway if this code is executed
                    DP.getInstance().w("Length of PDB number at line " + (i + 1) + " of DSSP file is 0." );
                }

                lastChar = tmpPdbResNum.charAt(tmpPdbResNum.length() - 1);
                if(Character.isDigit(lastChar)) {
                    // The last char is a digit, so the iCode field is empty.
                    foundPdbICode = " ";
                    foundPdbResNum = Integer.valueOf(tmpPdbResNum);
                }
                else {
                    // The last character is NO digit so it is the iCode
                    if(tmpPdbResNum.length() <= 1) {
                        DP.getInstance().w("Length of PDB number at line " + (i + 1) + " of DSSP file is 1 and this is not a digit." );
                    }
                    else {
                        foundPdbResNum = Integer.valueOf(tmpPdbResNum.substring(0, tmpPdbResNum.length() - 1));
                        foundPdbICode = tmpPdbResNum.substring(tmpPdbResNum.length() - 2, tmpPdbResNum.length() - 1);
                    }
                }
            } catch(Exception e) {
                DP.getInstance().w("Something went wrong with parsing PDB number at line " + (i + 1) + " of DSSP file, ignoring." );
                //System.exit(-1);
            }

            //System.out.println("    Found PDB residue number '" + foundPdbResNum + "', looking for '" + prn + "'. DSSP # is '" + foundDsspResNum + "' here.");

            if(foundPdbResNum.equals(prn)) {
                
                if(foundChain.equals(chainID)) {
                    
                    if(foundPdbICode.equals(iCode)) {
                        resultDsspResNum = foundDsspResNum;
                        break;
                    }
                }                
            }            
        }

        //if(resultDsspResNum == null) {
        //    System.out.println("WARNING: getDsspResNumForPdbFields(): Could not find DSSP residue number for residue with PDB number " + prn + " and chainID '" + chainID + "' and iCode '" + iCode + "'.");
        //}
        
        return(resultDsspResNum);
    }
    
    
    /**
     * Determines whether a line from a DSSP file is a chainbreak line.
     * @return True if it is, false otherwise.
     */
    public static Boolean isChainBreakLine(String dsspLine) {
        if((dsspLine.substring(13, 14)).equals("!")) {
            return(true);
        }
        else {
            return(false);
        }
    }
    
    
    /**
     * Parses the DSSP data and creates the residue list from it.
     * @param isCIF true if using mmCIF parser and mmCIF pdb file as chain IDs may be 4 character long then
     */
    protected static void createAllResiduesFromDsspData(Boolean isCIF) {
        
        String dLine;
        Integer dLineNum, dsspResNum, pdbResNum, resIndex, acc;
        dLineNum = dsspResNum = pdbResNum = 0;
        String dsspChainID, resName1Letter, sseString, iCode;
        dsspChainID = resName1Letter = iCode = "";
        sseString = "UNDEF";
        Character lastChar = null;
        Float phi = 0.0f;
        Float psi = 0.0f;
        int offset;  // if > 99,999 residues everything is shifted
        
        // moved here and was previously executed >before< calling this function at all.
        // md: I moved it to init variables again so I could use the new function grabSseStringForDsspNum to get the sse String for a residue,
        // without creating the residues from this dssp file. This is for retrofitting, so a user could still use a cif file and an old .dssp file
        // dsspDataStartLine = readDsspToData();

        Residue r;

        for(Integer i = dsspDataStartLine - 1; i < dsspLines.size(); i++) {
            dLine = dsspLines.get(i);
            dLineNum = i + 1;
            
            // jnw_2020: This was introduced for file with where DSSP number exceeded 99999 and caused a shift of the whole line
            //   due to changes to the DSSP file this seems not to be needed anymore. Instead we encounter "--->" and the DSSP num is found in later column
            // offset = Math.max(dLine.split(" ")[0].length() - 5, 0);  // typically first 5 columns code for res num, but if exceeded add offset

            if(dLine.substring(13, 14).equals("!")) {       // chain brake
                chainBreakDssps.add(getDsspNumFromLine(dLine));
                if(! silent) {
                    if (! Settings.getBoolean("PTGLgraphComputation_B_no_chain_break_info")) {
                        System.out.println("    DSSP: Found chain brake at DSSP line " + dLineNum + ".");
                    }
                }
            }
            else {          // parse the residue line

                try {
                    dsspResNum = getDsspNumFromLine(dLine);

                    
                    // last used DSSP res num is later needed for ligands
                    // (and we only wand to go through dssp file once)
                    if (isCIF) {
                        // lets hope dssp res num always increases
                        lastUsedDsspNum = dsspResNum;
                    }
                    
                    // 5 is ignored: blank
                    pdbResNum = Integer.valueOf(dLine.substring(6, 10).trim());
                    iCode = dLine.substring(10, 11);                    
                    // with PDB mmCIF files things got more difficult: 4-character chain ids
                    //     prioritize AUTHCHAIN > CHAIN
                    if (dLine.length() >= 163){
                        dsspChainID = dLine.substring(159, 163).trim(); // AUTHCHAIN column 160-163
                    } else {
                        dsspChainID = dLine.substring(149, 153).trim();  // PDBCHAIN column
                    }
                    
                    // 12 is ignored: blank
                    resName1Letter = dLine.substring(13, 14);
                    // 14+15 are ignored: blank
                    sseString = dLine.substring(16, 17);
                    // lots of stuff is ignored here
                    phi = Float.valueOf(dLine.substring(103, 108).trim());  // phi backbone angle
                    psi = Float.valueOf(dLine.substring(109, 114).trim());  // psi backbone angle
                    acc = Integer.valueOf(dLine.substring(35, 38).trim());  // solvent accessible surface
                    // rest is ignored: not needed
                } catch(Exception e) {
                    System.err.println("ERROR: Parsing of DSSP line " + dLineNum + " failed: '" + e.getMessage() + "'. DSSP file looks broken.");
                    //e.printStackTrace();
                    continue;
                }

                // successfully parsed, create Residue and add all the info we know atm
                r = new Residue(pdbResNum, dsspResNum);
                //System.err.println("DEBUG: Created residue with PDB# " + pdbResNum + ", DSSP# " + dsspResNum + ".");
                //System.err.println("DEBUG: Residue: " + r.getFancyName() + ".");
                r.setChainID(dsspChainID);
                // don't set the chain here, we'll set it later with the pdbChainID while creating the atoms
                // EDIT: Why not?! I guess I thought that DSSP and PDB files used different chain IDs for residues when I wrote that.
                // EDIT jnw: do the matching in CIF files later so that each file is only read once
                if (! isCIF) {
                    FileParser.getChainByPdbChainID(dsspChainID).addMolecule(r);
                }


                r.setType(Residue.RESIDUE_TYPE_AA);                   // DSSP files only contain protein residues, ligands are ignored
                r.setPhi(phi);
                r.setPsi(psi);
                r.setAcc(acc);

                // Cysteins which form a sulfur bridge are arbitrary lowercase letters in dssp. A residue
                // with name 'a' means it is a cysteine, but it has forms a sulfur bridge with a second residue in the file that is
                // also labeled 'a'.
                try {
                    char c = (resName1Letter.toCharArray())[0];
                    if( ! Character.isUpperCase(c)) {
                        resName1Letter = "C";   // change residue code cysteine

                        // now go save the sulfur bridge
                        if(FileParser.s_sulfurBridges.containsKey(c)) {
                            // the sulfur bridge partner is already in there
                            (FileParser.s_sulfurBridges.get(c)).add(dsspResNum);

                            // Check whether its a interhcain sulfur bridge (different chain IDs)
                            if(! FileParser.s_interchainSulfurBridgesChainID.get(c).equals(dsspChainID)) {
                                ArrayList<Integer> tmpInterchain = new ArrayList<Integer>();

                                // Get the dsspResNum from the first residue of this interchain sulfur bridge
                                tmpInterchain.add(FileParser.s_sulfurBridges.get(c).get(0));
                                // Also add the dsspResNum of the current (second) residue.
                                tmpInterchain.add(dsspResNum);
                                FileParser.s_interchainSulfurBridges.put(c, tmpInterchain);
                            }


                        } else {
                            // this is the first residue of the sulfur bridge pair
                            ArrayList<Integer> tmp = new ArrayList<Integer>();
                            tmp.add(dsspResNum);
                            FileParser.s_sulfurBridges.put(c, tmp);

                            // If its the first residue of the sulfur bridge save the chain this residue belongs to.
                            // This will be used to check if the second residue is on the same change or not.
                            FileParser.s_interchainSulfurBridgesChainID.put(c, dsspChainID);
                        }
                    }
                } catch (Exception e) {
                    DP.getInstance().w("Something is fishy with residue in line " + dLineNum + " of the DSSP file: '" + e.getMessage() + "'. Sulfur bridge trouble?");
                }


                // Fix DSSP file like the one for 2ZW3.pdb which list cysteine residues as 'o' instead of 'c'.
                if(resName1Letter.equals("O") || resName1Letter.equals("U")) {
                    DP.getInstance().w("Turned AA with 1-letter code '" + resName1Letter + "' into 'C'.");
                    resName1Letter = "C";
                }

                // Define "asparagine or aspartic acid" as asparagine
                if(resName1Letter.equals("B")) {
                    DP.getInstance().w("Turned AA with 1-letter code '" + resName1Letter + "' into 'N'.");
                    resName1Letter = "N";
                }

                r.setiCode(iCode);
                r.setAAName1(resName1Letter);
                r.setName3(AminoAcid.name1ToName3(resName1Letter));
                r.setSSEString(sseString);
                // Note: the SSE itself will be set by the getAllSSEs*() functions, it is still NULL atm

                // set model id to '1'. This is always true because DSSP can't handle models and the input files have to be split for it, thus every PDB file it reads has only 1 model
                r.setModelID("1");


                // DEBUG 3eoj
                //System.err.println("DEBUG --- residue state ---");
                //if(r.getPdbResNum() == 8 || r.getPdbResNum() == 7 || r.getPdbResNum() == 213) {
                //    System.err.println("DEBUG: RESIDUE: " + r.getFancyName() + ".");
                //}

                // add to list of Residues
                FileParser.s_molecules.add(r);
                resIndex = FileParser.s_molecules.size() - 1;
                // add index to s_residueIndices
                FileParser.s_residueIndices.add(resIndex);
                // produces null ponter exception in CIF parser and I dont see where we need it (maybe I'll understand later)
                if (! isCIF) {
                    FileParser.resIndexDSSP[dsspResNum] = resIndex;
                }
                //resIndexPDB[pdbResNum] = resIndex;  // This will crash because some PDB files use negative residue numbers, omfg.
                //System.out.println("    DSSP: Added residue PDB # " +  pdbResNum + ", DSSP # " + dsspResNum + " to s_residues at index " + resIndex + ".");

                
            }
        }  // end of iterating over DSSP lines
        
        if (! silent) {
            System.out.println("    DSSP: Found " + FileParser.s_molecules.size() + " residues.");
        }
        
    }

    
    /**
     * Get the DSSP residue number (ID) from a DSSP line. Correctly treats large structures.
     * @param dsspLine line from a DSSP file (not header)
     * @return DSSP residue number
     */
    private static Integer getDsspNumFromLine(String dsspLine) {
        if (dsspLine.length() > 4) {
            if (!dsspLine.substring(0, 5).equals("---->")) {
                return Integer.valueOf(dsspLine.substring(0, 5).trim());
            } else {
                // larger number than fits into five digits -> DSSP stores them in columns to the right
                if (dsspLine.length() > 173) {
                    return Integer.valueOf(dsspLine.substring(168, 174).trim());
                }
                else{
                    DP.getInstance().e("dsspLine indicates large number to the right, but the line isn't long enough. Exiting.");
                    System.exit(1);
                }
            }
        }
        return null;
    }
    
    // md: doc
    // no optimization because it's only for retrofitting
    protected static String grabSseStringForDsspNum(Integer requestedDsspNum){
        String dline;
        
        for (Integer i = dsspDataStartLine - 1; i < dsspLines.size(); i++) {
            dline = dsspLines.get(i);
            if (requestedDsspNum.equals(getDsspNumFromLine(dline))){
                return dline.substring(16, 17);
            }
        }
        return null;
    }

    public static ArrayList<String> getDsspLines() {
        if(dataInitDone) {
            return(dsspLines);
        }
        else {
            System.err.println("ERROR: FileParser.getDsspLines(): Request for data before initData() was called.");
            System.exit(1);
            return(null);
        }
    }
    
        
    protected static Integer getLastUsedDsspResNumOfDsspFile() {

        Integer dLineNum, dsspResNum;
        dLineNum = dsspResNum = 0;
        String dLine;
        dLine = "";

        if(maxUsedDsspResNumInDsspFile == null) {
            // value is not known yet, we have to determine it from the DSSP file
            for(Integer i = dsspDataStartLine - 1; i < dsspLines.size(); i++) {
                dLine = dsspLines.get(i);
                dLineNum = i + 1;
                
                if( (dLine.length() >= 14)) { // Warning has been printed before
                    if(! dLine.substring(13, 14).equals("!")) {       // if this is NOT a chain brake line

                        try {
                            dsspResNum = Integer.valueOf(dLine.substring(1, 5).trim());
                        } catch (Exception e) {
                            System.err.println("ERROR: Parsing of DSSP line " + dLineNum + " failed. DSSP file broken.");
                            e.printStackTrace();
                            System.exit(1);
                        }
                    }
                }
            }

            // DSSP residue numbers are ordered in the file so after iterating through
            //  all lines of the file, this holds the last used DSSP residue number.

            maxUsedDsspResNumInDsspFile = dsspResNum;
            //System.out.println("    DSSP: Last used DSSP residue number is " + maxUsedDsspResNumInDsspFile + ", found in line " + dLineNum + ".");


        }

        // value is known now and saved in class global var, just return it.
        return(maxUsedDsspResNumInDsspFile);
    }
    
}
