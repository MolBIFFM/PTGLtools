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
import io.FileParser;
import java.util.ArrayList;
import java.util.Arrays;
import io.IO;
import settings.Settings;
import tools.DP;
import plcc.Main;

/**
 * Represents a protein chain in a PDB file.
 * @author ts
 */
public class Chain implements java.io.Serializable {

    // declare class vars
    private String pdbChainID = null;                // author chain ID from PDB file
    private String altChainID = null;               // alternative chain ID from DSSP file
    private ArrayList<Molecule> molecules = new ArrayList<>();      // a list of all Molecules of the Chain
    private String macromolID = null;                // the macromolecule ID of the chain in the PDB file, defines chains forming a single macromolecule
    private String macromolName = null;              // the macromol name from the PDB file
    private String modelID = null;
    private Model model = null;                                                                                                                                                                                                                             // the Model of this Chain
    private ArrayList<String> homologues = null;     // a list of homologous chains (defined by PDB COMPND)
    private final Integer[] chainCentroid = new Integer[3];  // X-/Y-/Z-coordinates as 10th of Angström of the center of all non-H atoms
    private Integer radiusFromCentroid = null;         // distance from center to farthest non-H atom. -1 if no protein-atoms
    private Integer maxSeqNeighborAADist = null;    // largest distance between sequential residue neighbors excluding ligands (center to center)
    private String moleculeType = null;             // type of molecules that make up the chain, i.e. RNA or Residue
    
    // constructor
    public Chain(String ci) { pdbChainID = ci; molecules = new ArrayList<>();}


    // getters
    /** Returns the chain ID (chain name) from PDB file. Something like "A". */
    public String getPdbChainID() { return(pdbChainID); }
    public String getAltChainID() { return(altChainID); }
    public String getModelID() { return(modelID); }
    public String getMacromolID() { return(macromolID); }
    public String getMacromolName() { return(macromolName); }
    public Model getModel() { return(model); }
    public String getMoleculeType() { return (moleculeType); }
    public ArrayList<Molecule> getMolecules() { return(molecules); }
    public ArrayList<String> getHomologues() { return(homologues); }
    
    public Boolean isRnaChain() { return(this.getMoleculeType().contains("polyribonucleotide")); }
    
    /**
     * Retrieves (and computes if called 1st time) the coordinates of the chain centroid.
     * @return 
     */
    public Integer[] getChainCentroid() {
        if (chainCentroid == null) {
            this.computeChainCentroidAndRadius();
        }
        return(chainCentroid); 
    }
    
    
    /**
     * Retrieves (and computes if called 1st time) the distance from chain center to farthest atom.
     * @return radius as 10th of Angström
     */
    public Integer getRadiusFromCentroid() {
        if (radiusFromCentroid == null) {
            this.computeChainCentroidAndRadius();
        }
        return(radiusFromCentroid);
    }
    
        
    /**
     * Returns a list of all ligand residues in this chain.
     * @return a list of all ligand residues in this chain
     * For this method we need to convert the molecule object into a ligand object 
     */
    public ArrayList<Ligand> getAllLigandResidues() {
        ArrayList<Ligand> ligands = new ArrayList<>();     
        Ligand l = new Ligand();
        for(Molecule m : this.molecules) 
            if(m.isLigand()) {
                l = (Ligand) m;
                ligands.add(l);
            }
        return ligands;
    }
    
    
    /**
     * Returns all molecules of class Residue belonging to this chain.
     * @return ArrayList of residues
     */
    public ArrayList<Residue> getAllAAResidues() {
        ArrayList<Residue> thisResidues = new ArrayList<>();
        for (Molecule m : this.molecules) {
            if (m instanceof Residue) {
                thisResidues.add((Residue) m);
            }
        }
        return thisResidues;
    }
    
    
       /**
     * Returns all molecules of class Residue belonging to this chain as output type Molecule.
     * @return ArrayList of residues which are type Molecule
     */
    public ArrayList<Molecule> getAllAAResiduesAsTypeMolecule() {
        ArrayList<Molecule> thisResidues = new ArrayList<>();
        for (Molecule m : this.molecules) {
            if (m instanceof Residue) {
                thisResidues.add(m);
            }
        }
        return thisResidues;
    }
      
    
    /**
     * Returns a list of all RNA residues in this chain.
     * 
     */
    public ArrayList<RNA> getAllRnaResidues() {
        ArrayList<RNA> allRna = new ArrayList<>();     
        RNA r = new RNA();
        for(Molecule m : this.molecules) 
            if(m.isRNA()) {
                r = (RNA) m;
                allRna.add(r);
            }
        return allRna;
    }
    
    public Integer getChainLength(){
        Integer chainLength = this.getAllRnaResidues().size() + this.getAllAAResidues().size() + this.getAllLigandResidues().size();
        return chainLength;
    }

    public Integer getLength(boolean includeLigands){
        Integer chainLength = this.getAllAAResidues().size() + this.getAllRnaResidues().size();
        if (includeLigands){
            chainLength += this.getAllLigandResidues().size();
        }
        return chainLength;
    }
    
    /**
     * Retrieves (and calculates if called 1st time) the maximum sequence neighbor distance between amino acids.
     * @return the maximum sequence neighbor distance between amino acids in 10th of Angström
     */
    public Integer getMaxSeqNeighborAADist() {
        if (maxSeqNeighborAADist == null) {
            // Would make sense to have the function here in Chain, but the Main relies on it for the old (no chain-sphere speedup, no centroid) calculation.
            //   Would be possible to move the function here and make it static, though.
            maxSeqNeighborAADist = Main.getGlobalMaxSeqNeighborResDist(this.getAllAAResidues());
        }
        return maxSeqNeighborAADist;
    }
    

    // setters
    public void addMolecule(Molecule mol){molecules.add(mol);}
    public void setPdbChainID(String s) { pdbChainID = s; }
    public void setAltChainID(String s) { altChainID = s; }
    public void setMacromolID(String s) { macromolID = s; }
    public void setMacromolName(String s) { macromolName = s; }
    public void setModelID(String s) { modelID = s; }
    public void setModel(Model m) { model = m; }
    public void setMoleculeType (String s) { moleculeType = s; }
    public void setHomologues(ArrayList<String> h) { homologues = h; }
    public void addHomologue(String s) {
        if(!homologues.contains(s)){
            homologues.add(s);
        }
    }

    /**
     * Returns the chemical property string for this chain according to the 5 types system, i.e., a concatenation of all chemProps of all the chain residues.
     * @return the chemical property string, 5 types system
     */
    public String getChainChemProps5StringAllResidues() {
        StringBuilder sb = new StringBuilder();   
        if (this.molecules.size() > 0) {
            if ((this.molecules.get(0)) instanceof RNA || this.molecules.get(0) instanceof Ligand) {
                return null;
            }
        }
        Residue r;
        for(Molecule m : this.molecules) {
            r= (Residue) m;
            sb.append(r.getChemicalProperty5OneLetterString());            
        }        
        return sb.toString();
    }
    
    /**
     * Returns the chemical property string for this chain according to the 3 types system, i.e., a concatenation of all chemProps of all the chain residues.
     * @return the chemical property string, 3 types system
     */
    public String getChainChemProps3StringAllResidues() {
        StringBuilder sb = new StringBuilder();  
        if (this.molecules.size() > 0) {
            if ((this.molecules.get(0)) instanceof RNA || this.molecules.get(0) instanceof Ligand) {
                return null;
            }
        }
        
        Residue r;
        for(Molecule m : this.molecules) {
            r = (Residue)m;
            sb.append(r.getChemicalProperty3OneLetterString());            
        }        
        return sb.toString();
    }
    
    /**
     * Returns the chemical property string for this chain in the 5 types system, but considers only AA residues which are part of an SSE (ignores ligands and AAs in coiled regions).
     * @param sseSeparator the String (most likely a single character) to add between two SSEs. You can specify the empty String if you do not want anything inserted, of course.
     * @return the chemical property string, 5 types system
     */
    public String[] getChainChemPropsStringSSEResiduesOnly(String sseSeparator) {
        StringBuilder sbChemProp = new StringBuilder();        
        StringBuilder sbSSE = new StringBuilder();
        SSE s;
        String sseString = null;
        String sseStringlast = null;
        int numSeps = 0;
        
        // some basic stat stuff, just a quick test
        int numSSETypes = 2; int H = 0; int E = 1; // helix and strand
        int numResH = 0; int numResE = 0;
        Integer[][] chemPropsBySSEType = new Integer[numSSETypes][AminoAcid.ALL_CHEM_PROPS5.length];
        for(int i = 0; i < numSSETypes; i++) {
            Arrays.fill(chemPropsBySSEType[i], 0);
        }
        if (this.molecules.size() > 0) {
            if ((this.molecules.get(0)) instanceof RNA) {
                return null;
            }
        }
        
        Residue r;
        for(Molecule m : this.molecules) {
            r = (Residue )m;
            if(r.isAA()) {
                s = r.getSSE();
                
                if(s != null) {
                    if(s.isOtherSSE()) {
                        continue;
                    }
                    sseString = s.getSSEClass();
                    if(sseStringlast != null && (sseString == null ? sseStringlast != null : !sseString.equals(sseStringlast))) {
                        sbChemProp.append(sseSeparator);
                        sbSSE.append(sseSeparator);
                        numSeps++;
                    }
                    sbChemProp.append(r.getChemicalProperty5OneLetterString());
                    if(sseString.equals("H")) {
                        numResH++;
                        chemPropsBySSEType[H][r.getChemicalProperty5Type()]++;
                    }
                    if(sseString.equals("E")) {
                        numResE++;
                        chemPropsBySSEType[E][r.getChemicalProperty5Type()]++;
                    }
                    sbSSE.append(s.getSSEClass());
                }
                sseStringlast = sseString;
            }            
        }
        
        // stat stuff
        Boolean doStats = false;
        if(doStats) {
            float[] sharesH = new float[AminoAcid.ALL_CHEM_PROPS5.length];
            float[] sharesE = new float[AminoAcid.ALL_CHEM_PROPS5.length];
            for(int i = 0; i < AminoAcid.ALL_CHEM_PROPS5.length; i++) {
                sharesH[i] = (float)numResH / (float)chemPropsBySSEType[H][i];
                sharesE[i] = (float)numResE / (float)chemPropsBySSEType[E][i];
            }
            System.out.println("Chain: chemType by SSEtype[H]: totalH=" + numResH + ", shares: " + IO.floatArrayToString(sharesH));
            System.out.println("Chain: chemType by SSEtype[E]: totalE=" + numResE + ", shares: " + IO.floatArrayToString(sharesE));
        }
        
        //System.out.println("Added " + numSeps + " separators.");
        return new String[] { sbChemProp.toString(), sbSSE.toString() };
    }
    
    
    /**
     * Computes the geometrical center of all atoms and the largest distance from center to an atom (=radius).
     */
    private void computeChainCentroidAndRadius() {
        // compute center
        Integer[] tmpCenter = new Integer[3];
        tmpCenter[0] = tmpCenter[1] = tmpCenter[2] = 0;
        int tmpAtomNumber = 0;
        for (Molecule mol : molecules) {
            for (Atom a : mol.getAtoms()) {
                tmpCenter[0] += a.getCoordX();
                tmpCenter[1] += a.getCoordY();
                tmpCenter[2] += a.getCoordZ();
                tmpAtomNumber += 1;
            }
        }
        
        // there could be the case that only RNA/DNA atoms are in a chain
        // in this case radius is set to -1
        
        if (tmpAtomNumber > 0) {
        
            chainCentroid[0] = (int) (Math.round((double) tmpCenter[0] / tmpAtomNumber));
            chainCentroid[1] = (int) (Math.round((double) tmpCenter[1] / tmpAtomNumber));
            chainCentroid[2] = (int) (Math.round((double) tmpCenter[2] / tmpAtomNumber));

            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
                System.out.println("[DEBUG] Center of chain " + pdbChainID + " is at " + Arrays.toString(chainCentroid));
            }

            // compute radius
            int tmpBiggestDist = 0;
            int tmpCurrentDist;
            for (Molecule mol : molecules) {
                for (Atom a : mol.getAtoms()) {
                    tmpCurrentDist = a.distToPoint(chainCentroid[0], chainCentroid[1], chainCentroid[2]);
                    // System.out.println("[DEBUG] Distance to center from atom " + a.toString() + " is " + String.valueOf(tmpCurrentDist));
                    if (tmpCurrentDist > tmpBiggestDist) {
                        tmpBiggestDist = tmpCurrentDist;
                    }
                }
            }

            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
                System.out.println("[DEBUG] Radius of chain " + pdbChainID + " is " + String.valueOf(tmpBiggestDist));
            }

            radiusFromCentroid = tmpBiggestDist;
                
        } else {
            DP.getInstance().w("Chain " + this.pdbChainID + " seems not to hold protein atoms. No center can be detected.", 2);
            radiusFromCentroid = -1;
        }        
    }
    
    /**
     * This function determines whether we need to look at the residues to check for contacts between
     * this chain and another one. If the center spheres don't overlap, there cannot exist any atom contacts.
     * @param c Chain: the other chain
     * @return Bool: if spheres overlap
     */
    public Boolean contactPossibleWithChain(Chain c) {
        
        // no contact possible if either chain contains no atoms
        if ((! this.containsAtoms()) || (! c.containsAtoms())) {
            return false;
        }

        Integer dist, tmpSum;
        
        tmpSum = 0;
        tmpSum += (this.getChainCentroid()[0] - c.getChainCentroid()[0]) * (this.getChainCentroid()[0] - c.getChainCentroid()[0]);
        tmpSum += (this.getChainCentroid()[1] - c.getChainCentroid()[1]) * (this.getChainCentroid()[1] - c.getChainCentroid()[1]);
        tmpSum += (this.getChainCentroid()[2] - c.getChainCentroid()[2]) * (this.getChainCentroid()[2] - c.getChainCentroid()[2]);
        
        dist = (int)Math.round(Math.sqrt(tmpSum));
        
        // presume there are ligands in the chain
        Integer atomRadius = Settings.getInteger("PTGLgraphComputation_I_lig_atom_radius");
        
        Integer justToBeSure = 4;   // account for small errors due to rounding
        Integer summedSpheres = this.getRadiusFromCentroid() + c.getRadiusFromCentroid() + (atomRadius * 2) + justToBeSure;

        //System.out.println("    Center sphere radius for PDB residue " + this.getPdbResNum() + " = " + this.getCenterSphereRadius() + ", for " + r.getPdbResNum() + " = " + r.getCenterSphereRadius() + ", atom radius is " + atomRadius + ".");
        //System.out.println("    DSSP Res distance " + this.getDsspResNum() + "/" + r.getDsspResNum() + " is " + dist + " (no contacts possible above distance " + maxDistForContact + ").");

        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 0) {
            System.out.println("[DEBUG][CHAIN] Chain " + this.pdbChainID + " and " + c.pdbChainID);
            System.out.println(" ... mid points: " + this.chainCentroid[0] + "|" + this.chainCentroid[1] + "|" + this.chainCentroid[2]);
            System.out.println(" ... mid points: " + c.chainCentroid[0] + "|" + c.chainCentroid[1] + "|" + c.chainCentroid[2]);
            System.out.println(" ... radii: " + this.radiusFromCentroid + " and " + c.radiusFromCentroid);
            System.out.println(" ... distance: " + dist);
            System.out.println(" ... summedSpheres: " + summedSpheres);
        }
        
        if(dist <= summedSpheres) {
            return(true);
        }
        else {
            return(false);
        }
    }
    
    public Boolean containsAtoms() {
        for (Molecule mol : this.getMolecules()) {
            for (Atom a : mol.getAtoms()) {
                return true;
            }
        }
        return false;
    }
}
