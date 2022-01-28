/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package proteinstructure;

// imports
import proteinstructure.SSE;
import tools.DP;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import plcc.Main;
import settings.Settings;


/**
 * Represents a Residue. Based on information in the PDB and DSSP file.
 * 
 * @author ts
 */
public class Residue extends Molecule implements java.io.Serializable {
   
    
    // declare class vars
                             // 3 letter name
                              // 1 letter name, AAs only
       // PDB insertion code
    
    /** The binding sites, if any, that this residue is part of (for protein residues which are part of the pocket, NOT for ligands which dock into a pocket).  */
    private List<BindingSite> partOfBindingSites;
       
    public static HashMap<String,String> AANames = new HashMap<String, String>(){{
       put("ARG","R");
       put("HIS","H");
       put("LYS","K");
       put("ASP","D");
       put("GLU","E");
       put("SER","S");
       put("THR","T");
       put("ASN","N");
       put("GLN","Q");
       put("CYS","C");
       put("SEC","U");
       put("GLY","G");
       put("PRO","P");
       put("ALA","A");
       put("ILE","I");
       put("LEU","L");
       put("MET","M");
       put("PHE","F");
       put("TRP","W");
       put("TYR","Y");
       put("VAL","V");
    }};
    
    
    /**
     * Returns the one-letter code for amino acids from AAName Map.
     * @param AAName3: Three letter code
     * @return One letter code
     */
    public static String getAAName1fromAAName3(String AAName3) {
        return AANames.getOrDefault(AAName3, "X");
    }
    
    /**
     * Determines all AltLoc identifiers which have an alpha carbon atom assigned to them for this residue.
     * @return the list of all altLocs which include alpha carbons
     */ 
    protected ArrayList<String> getAltLocsWithAlphaCarbonAtoms() {
        
        if(! this.isAA()) {
            DP.getInstance().w("Non-AA residue cannot contain alpha carbon atom but asked to look for one.");
        }
        
        ArrayList<String> altLocsWithAlphaCarbon = new ArrayList<String>();
        
        String candidateAltLoc;
        for(Atom a : this.atoms) {
            if(a.isCalphaAtom()) {
                // Be careful not to add an AltLoc which has several C-alpha 
                // atoms listed twice (should never happen and makes no sense, but
                // with PDB files you never know).
                candidateAltLoc = a.getAltLoc();
                if( ! altLocsWithAlphaCarbon.contains(candidateAltLoc)) {
                    altLocsWithAlphaCarbon.add(candidateAltLoc);
                }
            }
        }
        
        return altLocsWithAlphaCarbon;
    }
    
    
    /**
     * Tells this molecule to choose its alternate location PDB identifier and delete all its
     * atoms which have other AltLocs. For amino acids, this will choose the AltLoc which has a C-alpha atom
     * and maximizes the number of atoms in the molecule. For non-AAs, this function will always choose the AltLoc identifier which
     * maximizes the number of atoms in this molecule. The list of deleted atoms is returned so they can be deleted
     * from the global atom list as well.
     * @return the list of the atoms that were deleted
     */
    @Override public ArrayList<Atom> chooseYourAltLoc() {
        
        if (! this.isAA()) {
            return super.chooseYourAltLoc();
        }
        
        int numAtomsBefore = this.atoms.size();
        if(numAtomsBefore < 1) {
            if(! Settings.getBoolean("PTGLgraphComputation_B_no_parse_warn")) {
                DP.getInstance().w("Molecule " + this.getFancyName() + " of chain " + this.getChainID() + " has NO atoms at all *before* choosing alternative location PDB field and deleting others.");
            }
            return new ArrayList<Atom>();
        }                
        ArrayList<Atom> deletedAtoms;
    
        String chosenAltLoc= getAltLocWithAlphaCarbonAndMostAtoms();
    
        int numAtomsWithChosenAltLoc = this.getNumAtomsWithAltLoc(chosenAltLoc);
    
        if((numAtomsWithChosenAltLoc < 1 || numAtomsWithChosenAltLoc > Main.MAX_ATOMS_PER_AA) && this.isAA()) {
            if( ! Settings.getBoolean("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce")) {
                DP.getInstance().w("Chosen altLoc '" + chosenAltLoc + "' leads to " + numAtomsWithChosenAltLoc + " atoms for AA molecule " + this.getFancyName() + ".");
            }
        }
              
        deletedAtoms = this.deleteAtomsWithAltLocDifferentFrom(chosenAltLoc);
                
        if(this.atoms.size() < 1) {
            DP.getInstance().w("Molecule  " + this.getFancyName() + " of chain " + this.getChainID() + " has no atoms after choosing alternative location PDB field (had " + numAtomsBefore + " before).");
        }
        
        return deletedAtoms;
    }
    
        /**
     * Returns the alternate location identifier that most atoms of this molecule share.
     * @return the most common alternate location identifier of this molecule. If this molecule
     * has no atoms, a space (" ", the PDB default altLoc) will be returned.
     */
       public String getAltLocWithAlphaCarbonAndMostAtoms() {
        
        ArrayList<String> altLocsCA = this.getAltLocsWithAlphaCarbonAtoms();
        if(altLocsCA.size() < 1) {
            DP.getInstance().w("Molecule  " + this.getFancyName() + " has no AltLoc which includes an alpha carbon atom.");
            return this.getAltLocWithMostAtoms();
        }
        
        if(altLocsCA.size() == 1) {
            return altLocsCA.get(0);
        }
        
        String maxAltLoc = " ";
        
        // add all altLocs
        String altLoc;
        Integer newCount;
        Integer maxCount = 0;
        HashMap<String, Integer> atomCountsByAltLoc = new HashMap<String, Integer>();
        for(Atom a : this.atoms) {
            altLoc = a.getAltLoc();
            
            // skip AltLocs which do not contain C alpha atoms
            if( ! altLocsCA.contains(altLoc)) {
                continue;
            }
            
            if(atomCountsByAltLoc.containsKey(altLoc)) {
                // altLoc already listed, so just increase count
                newCount = atomCountsByAltLoc.get(altLoc) + 1;
            } else {
                // altLoc not listed yet so this is the first atom with this type, put it in there
                newCount = 1;                
            }
            atomCountsByAltLoc.put(altLoc, newCount);
            
            // keep track of maximum
            if(newCount > maxCount) {
                maxAltLoc = altLoc;
                maxCount = newCount;
            }
        }
        
        //DEBUG
        /*
        if(this.pdbResNum == 209) {
            System.out.println("DEBUG: " + atomCountsByAltLoc.keySet().size() + " Alternate locations and their atoms counts for " + this.toString() + ": ");
            for(String al : atomCountsByAltLoc.keySet()) {
                System.out.println(" '" + al + "' has " +  atomCountsByAltLoc.get(al) + " atoms (of " + this.atoms.size() + " total).");
            }
            System.out.println("DEBUG: MaxAltLoc=" + maxAltLoc + ", maxCount=" + maxCount + ".");
        }
        */
                
        return maxAltLoc;
    }

    public List<BindingSite> getPartOfBindingSites() {
        return partOfBindingSites;
    }

    public void setPartOfBindingSites(List<BindingSite> partOfBindingSites) {
        this.partOfBindingSites = partOfBindingSites;
    }
    
    /** Assigns a new binding site the residue is part of.
     * @param s the binding site the residue s part of
     */
    public void addPartOfBindingSite(BindingSite s) {
        this.partOfBindingSites.add(s);
    }
    
    private Float phi = null;       // phi backbone angle
    private Float psi = null;       // psi backbone angle
    private Integer acc = null;     // solvent accessible surface, from DSSP file

    
    /**
     * Determines the chemical properties of this AA, according to the 5 types system.
     * @return the chemical property, which is one of the constants at AminoAcid.CHEMPROP5_AA_INT_* (Example: ALA => AminoAcid.CHEMPROP5_AA_INT_SMALL_APOLAR).
     */
    public Integer getChemicalProperty5Type() {
        if(this.isAA()) {            
            return AminoAcid.getChemProp5OfAAByName3(this.Name3);
        }
        else {
            return AminoAcid.CHEMPROP5_AA_INT_UNKNOWN;            
        }
    }
    
    /**
     * Determines the chemical properties of this AA, according to the 3 types system.
     * @return the chemical property, which is one of the constants at AminoAcid.CHEMPROP3_AA_INT_* (Example: ALA => AminoAcid.CHEMPROP3_AA_INT_HYDROPHOBIC).
     */
    public Integer getChemicalProperty3Type() {
        if(this.isAA()) {            
            return AminoAcid.getChemProp3OfAAByName3(this.Name3);
        }
        else {
            return AminoAcid.CHEMPROP3_AA_INT_UNKNOWN;            
        }
    }
    
    /**
     * Determines the chemical properties of this AA in the 5 types system and returns the string.
     * @return the chemical property string, which is one of the constants at AminoAcid.CHEMPROP5_AA_STRING_* (Example: ALA => AminoAcid.CHEMPROP5_AA_SMALL_APOLAR).
     */
    public String getChemicalProperty5OneLetterString() {
        return AminoAcid.getChemProp5_OneLetterString(AminoAcid.getChemProp5OfAAByName3(this.Name3));
    }
    
    /**
     * Determines the chemical properties of this AA in the 3 types system and returns the string.
     * @return the chemical property string, which is one of the constants at AminoAcid.CHEMPROP3_AA_STRING_* (Example: ALA => CHEMPROP3_AA_STRING_HYDROPHOBIC).
     */
    public String getChemicalProperty3OneLetterString() {
        return AminoAcid.getChemProp3_OneLetterString(AminoAcid.getChemProp3OfAAByName3(this.Name3));
    }
    
    // constructor
    public Residue() { this.atoms = new ArrayList<>(); this.hydrogenatoms = new ArrayList<>(); this.partOfBindingSites = new ArrayList<>(); }

    public Residue(Molecule Mol) {
        this.atoms = Mol.getAtoms();
        this.chain = Mol.getChain();
        this.chainID = Mol.getChainID();
        this.modelID = Mol.getModelID();
        this.iCode = Mol.getiCode();
        this.pdbNum = Mol.getPdbNum();
        this.dsspNum = Mol.getDsspNum();
        this.Name3 = Mol.getName3();
        this.AAName1 = Mol.getAAName1();
        this.sseString = Mol.getSSEString();
        this.sse = Mol.getSSE();
        this.sseStringDssp = Mol.getSSEStringDssp();
        this.isPartOfDsspSse = Mol.getDsspSseState();
    }
    
    /**
     * Constructs a new residue with PDB residue number 'prn' and DSSP residue number 'drn'.
     * @param residueNumberPDB the PDB residue number
     * @param residueNumberDSSP the DSSP residue number
     */
    public Residue(Integer residueNumberPDB, Integer residueNumberDSSP) {
        atoms = new ArrayList<>();
        hydrogenatoms = new ArrayList<>();
        pdbNum = residueNumberPDB;
        dsspNum = residueNumberDSSP;
        this.partOfBindingSites = new ArrayList<>();
    }
    
    
    /**
     * Determines whether this residue is part of any binding site.
     * @return whether this residue is part of any binding site
     */
    public Boolean isBindingSiteResidue() {
        return(this.partOfBindingSites.size() > 0);
    }
    

    /**
     * Returns the C alpha atoms of this residue or null if it has none.
     * @return the alpha carbon or null 
     */
    public Atom getAlphaCarbonAtom() {
        if(this.isAA()) {
            for(Atom a : this.atoms) {
                if(a.isCalphaAtom()) {
                    return a;
                }
            }
        }
        return null;
    }
    
    
    /**
     * Sets the center atom of this molecule to be the Calpha, and also sets the center sphere radius for the molecule.
     * If the molecule is not an AA or no Calpha can be found, the original getCenterAtom() function from class Molecule is used.
     * @return the center atom
     */
    @Override public Atom getCenterAtom() {
        
        // For non-AAs and AAs without Calpha use the function in Molecule class
        if((! this.isAA()) || this.getAlphaCarbonAtom() == null) {
            return super.getCenterAtom();
        }
        
    	Atom a, b, center;
        a = b = center = null;
        Integer maxDistForAtom, dist = 0;      // just assign a small start value
//      //Integer atomRadius = Settings.getInteger("PTGLgraphComputation_I_aa_atom_radius");

        if(atoms.size() < 1) {
            if( ! Settings.getBoolean("PTGLgraphComputation_B_no_parse_warn")) {
                DP.getInstance().w("getCenterAtom(): PDB molecule  " + this.pdbNum + " chain " + this.getChainID() + " of type " + getName3() + " has " + atoms.size() + " atoms in default location, returning null.");
            }
            return(null);
        }

            center = this.getAlphaCarbonAtom();
             
             // For the return value of this function alone we would be done, but we want to set the C alpha sphere
             //  radius as well. It is the maximal distance of the center atom to any other atom of this molecule.
            maxDistForAtom = 0;
            for(Integer j = 0; j < this.atoms.size(); j++) {

                b = this.atoms.get(j);
                 
                if(center.equalsAtom(b)) {
                   continue; 
                }
                 
                dist = center.distToAtom(b);

                if(dist > maxDistForAtom) {
                    maxDistForAtom = dist;
                }

            }

            // The distance to any other atom cannot be smaller than 2 * atomRadius, otherwise the vdW radii of
            //  the atoms would overlap. If this occurs in a PDB file, something most likely is wrong with the file.
            //  Note though that this only applies to AA molecules because a ligand could consist of a single atom
            //  and in that case the center sphere radius can be smaller than 2 * atomRadius.
            //if(maxDistForAtom < (2 * atomRadius)) {
            //    maxDistForAtom = 2 * atomRadius;
        

            // The maximal distance of the Atom we chose to any other Atom of this molecule is the C alpha/center sphere radius
        this.centerSphereRadius = maxDistForAtom;
         
        return (center);
    }
    
    
    @Override public String toString() {
        return("[Residue] PDB# " + pdbNum + ", DSSP# " + dsspNum + ", Type " + type + ", AA1 " + AAName1 + ", AA3 " + Name3 + ", Chain " + chainID + ", Model " + modelID + ", # of Atoms " + atoms.size());
    }
    
    public Float getPhi() { if(this.isAA()) { return(phi); } else { return(0.0f); } }
    public Float getPsi() { if(this.isAA()) { return(psi); } else { return(0.0f); } }
    public Integer getAcc() { if(this.isAA()) { return(acc); } else { return(0); } }
    

 
    /**
     * Determines whether this residue lies at the protein surface.
     * @return whether this residue lies at the protein surface.
     */
    public Boolean isSurfaceResidue() {
        if(this.acc == null) {
            return true;    // if we do not have SAS data, assume yes
        }
        return(this.acc >= 15);
    }
  
    //setters
    public void setPhi(Float f) { phi = f; }
    public void setPsi(Float f) { psi = f; }
    public void setAcc(Integer f) { acc = f; }
    
}
    
    
