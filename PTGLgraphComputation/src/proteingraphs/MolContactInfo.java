/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package proteingraphs;

import java.util.ArrayList;
import proteinstructure.Residue;
import java.util.Arrays;
import plcc.Main;
import proteinstructure.Atom;
import proteinstructure.Molecule;


/**
 * This class implements a molecule level contact, i.e., all information required about the contact between a pair (a, b) of Molecule objects.
 * @author ts
 */
public class MolContactInfo {

    // declare class vars

    // The contents of all the fields of these arrays is explained in Main.calculateAtomContactsBetweenResidues(),
    //  the function that creates these arrays.
    //  Just to be sure a copy of these lines can be found below. The length of all of these is 12.
    
    private Integer[] numPairContacts;
    // The positions in the numPairContacts array hold the number of contacts of each type for a pair of molecule:
    // Some cheap vars to make things easier to understand (a poor replacement for #define):
    public static final Integer TT = 0;         //  0 = total number of contacts            (all molecule type combinations)
    public static final Integer BB = 1;         //  1 = # of backbone-backbone contacts     (protein - protein only)
    public static final Integer CB = 2;         //  2 = # of sidechain-backbone contacts    (protein - protein only)
    public static final Integer BC = 3;         //  3 = # of backbone-sidechain contacts    (protein - protein only)
    public static final Integer CC = 4;         //  4 = # of sidechain-sidechain contacts   (protein - protein only)
    public static final Integer HB = 5;         //  5 = # of H-bridge contacts 1, N=>0      (protein - protein only)
    public static final Integer BH = 6;         //  6 = # of H-bridge contacts 2, 0=>N      (protein - protein only)
    public static final Integer BL = 7;         //  7 = # of backbone-ligand contacts       (protein - ligand only)
    public static final Integer LB = 8;         //  8 = # of ligand-backbone contacts       (protein - ligand only)
    public static final Integer CL = 9;         //  9 = # of sidechain-ligand contacts      (protein - ligand only)
    public static final Integer LC = 10;        // 10 = # of ligand-sidechain contacts      (protein - ligand only)
    public static final Integer LL = 11;        // 11 = # of ligand-ligand contacts         (ligand - ligand only)
    public static final Integer RX = 12;        // 12 = # of RNA-X contacs (X can be protein or ligand)
    public static final Integer XR = 13;        // 13 = # of X-RNA contacs (X can be protein or ligand)
    public static final Integer RR = 14;        // 14 = # of RNA-RNA contacts
    public static final Integer CTL = 15;        // 15 = # of transitive (chain-ligand-chain) contacts by ligands (protein - protein only)
    public static final Integer DISULFIDE = 16; // 16 = # of disulfide bridges
    public static final Integer BBHB = 17;      // 17 = # of interchain H-bridge contacts, backbone-backbone, donor-acceptor (N=>O)
    public static final Integer BBBH = 18;      // 18 = # of interchain H-bridge contacts, backbone-backbone, acceptor-donor (O=>N)
    public static final Integer IVDW = 19;      // 19 = # of interchain van der Waals interactions
    public static final Integer ISS = 20;       // 20 = # of interchain disulfide bridges
    public static final Integer IPI = 21;       // 21 = # of interchain pi-effects
    public static final Integer ISB = 22;       // 22 = # of interchain salt bridges
    public static final Integer BCHB = 23;      // 23 = # of interchain H-bridge contacts, backbone-sidechain, donor-acceptor
    public static final Integer BCBH = 24;      // 24 = # of interhcain H-bridge contacts, backbone-sidechain, acceptor-donor
    public static final Integer CBHB = 25;      // 25 = # of interchain H-bridge contacts, sidechain-backbone, donor-acceptor
    public static final Integer CBBH = 26;      // 26 = # of interchain H-bridge contacts, sidechain-backbone, acceptor-donor
    public static final Integer CCHB = 27;      // 27 = # of interchain H-bridge contacts, sidechain-sidechain, donor-acceptor
    public static final Integer CCBH = 28;      // 28 = # of interchain H-bridge contacts, sidechain-sidechain, acceptor-donor
    public static final Integer NHPI = 29;      // 29 = # of interchain N-H...Pi contacts, backbone-sidechain, donor-acceptor (non-canonical interaction)
    public static final Integer PINH = 30;      // 30 = # of interchain N-H...Pi contacts, sidechain-backbone, acceptor-donor (non-canonical interaction)
    public static final Integer CAHPI = 31;     // 31 = # of interchain CA-H...Pi contacts, backbone-sidechain, donor-acceptor (non-canonical interaction)
    public static final Integer PICAH = 32;     // 32 = # of interchain CA-H...Pi contacts, sidechain-backbone, acceptor-donor (non-canonical interaction) 
    public static final Integer CNHPI = 33;     // 33 = # of interchain N-H...Pi contacts, sidechain-sidechain, donor-acceptor (non-canonical interaction)
    public static final Integer PICNH = 34;     // 34 = # of interchain N-H...Pi contacts, sidechain-sidechain, acceptor-donor (non-canonical interaction)
    public static final Integer SHPI = 35;      // 35 = # of interchain S-H...Pi contacts, sidechain-sidechain, donor-acceptor (non-canonical interaction)
    public static final Integer PISH = 36;      // 36 = # of interchain S-H...Pi contacts, sidechain-sidechain, acceptor-donor (non-canonical interaction)
    public static final Integer XOHPI = 37;     // 37 = # of interchain X-O-H...Pi contacts witch X = Ser || Thr || Tyr, sidechain-sidechain, donor-accepor (non-canonical interaction
    public static final Integer PIXOH = 38;     // 38 = # of interchain X-O-H...Pi contacts witch X = Ser || Thr || Tyr, sidechain-sidechain, acceptor-donor(non-canonical interaction
    public static final Integer PROCDHPI = 39;  // 39 = # of interchain Pro-CD-H...Pi contacts, sidechain-sidechain, donor-acceptor (non-canonical)
    public static final Integer PIPROCDH = 40;  // 40 = # of interchain Pro-CD-H...Pi contacts, sidechain-sidechain, acceptor-donor (non-canonical)
    public static final Integer CCAHCO = 41;    // 41 = # of interchain CA-H...O=C contacts, mainchain-sidechain, donor-acceptor (non-canonical)
    public static final Integer CCOCAH = 42;    // 42 = # of interchain CA-H...O=C contacts, mainchain-sidechain, acceptor-donor (non-canonical)
    public static final Integer BCAHCO = 43;    // 43 = # of interchain CA-H...O=C contacts, mainchain-mainchain, donor-acceptor (non-canonical)
    public static final Integer BCOCAH = 44;    // 44 = # of interchain CA-H...O=C contacts, mainchain-mainchain, acceptor-donor (non-canonical)
    
    
    
    // The different values in 'XY' vs 'YX' are only produced by the the sequential order of the
    //  molecules (e.g., the ligand came fist: this was treated as an 'LB' contact, otherwise it
    //  had been treated as an 'BL' contact (note that molecule pairs are compared only once by the
    //  loop in Main.calculateAllContacts()). Saving this separately seems strange but we have to
    //  do it because geom_neo does it and we want our output to be compatible with it.

    public Integer numTotalLigContactsPair;
    public Integer numTotalRnaContactsPair;

    public Integer[] minContactDistances;
    // Holds the minimal distances of contacts of the appropriate type (see numPairContacts, index 0 is unused)

    public Integer[] contactAtomIndexInMoleculeA;
    // Holds the number Atom x has in its molecule a for the contact with minimal distance of that type.
    // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious)

    public Integer[] contactAtomIndexInMoleculeB;
    // Holds the number Atom y has in its molecule b for the contact with minimal distance of that type.
    // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused because the atom is
    //  fixed for these (HB1: backbone N, HB2: backbone O)

    public Molecule molA, molB;
    public Integer dist;
    
    public ArrayList<String> atomAtomContactType;
 //   private Atom[] atomAtomContact;
    public ArrayList<Atom[]> atomAtomContacts;

    /**
     * Constructor for a molecule pair contact between the molecules (a, b).
     * 
     * @param npcs array holding the number of atom contacts of all different types (e.g., BB, BC, CB, ...) of this molecule pair.
     * @param mcds array holding the minimum contact distances of all different types (BB, ...)
     * @param can_a array holding the atom indices of the contact atoms of the different types (BB, ...) in molecule a
     * @param can_b array holding the atom indices of the contact atoms of the different types (BB, ...) in molecule b
     * @param a the first molecule of (a, b)
     * @param b the second molecule of (a, b)
     * @param d the distance of this molecule pair (molecule center to molecule center)
     * @param nlc the number of total ligand contacts of this molecule pair 
     */    
    public MolContactInfo(Integer[] npcs, Integer[] mcds, Integer[] can_a, Integer[] can_b, Molecule a, Molecule b, Integer d, Integer nlc, Integer nrc) {

        numPairContacts = npcs;
        minContactDistances = mcds;
        contactAtomIndexInMoleculeA = can_a;
        contactAtomIndexInMoleculeB = can_b;
        molA = a;
        molB = b;
        dist = d;
        numTotalLigContactsPair = nlc;
        numTotalRnaContactsPair = nrc;
    }
    
    
    /**
     * Constructor for a molecule pair contact between the molecules (a, b).
     * 
     * @param npcs array holding the number of atom contacts of all different types (e.g., BB, BC, CB, ...) of this molecule pair.
     * @param mcds array holding the minimum contact distances of all different types (BB, ...)
     * @param can_a array holding the atom indices of the contact atoms of the different types (BB, ...) in molecule a
     * @param can_b array holding the atom indices of the contact atoms of the different types (BB, ...) in molecule b
     * @param a the first molecule of (a, b)
     * @param b the second molecule of (a, b)
     * @param d the distance of this molecule pair (molecule center to molecule center)
     * @param nlc the number of total ligand contacts of this molecule pair 
     */    
    public MolContactInfo(Integer[] npcs, Integer[] mcds, Integer[] can_a, Integer[] can_b, Molecule a, Molecule b, Integer d, Integer nlc, Integer nrc, ArrayList<String> aact, ArrayList<Atom[]> aac) {

        numPairContacts = npcs;
        minContactDistances = mcds;
        contactAtomIndexInMoleculeA = can_a;
        contactAtomIndexInMoleculeB = can_b;
        molA = a;
        molB = b;
        dist = d;
        numTotalLigContactsPair = nlc;
        numTotalRnaContactsPair = nrc;
        atomAtomContactType = aact;
        atomAtomContacts = aac;
    }
    
    
    
    /**
     * Creates an empty MolContactInfo between the two molecules and fills in their distance. All other values
     * are set to zero or not initialized.
     * @param a the first molecule
     * @param b the second molecule
     */
    public MolContactInfo(Molecule a, Molecule b) {
        molA = a;
        molB = b;
        dist = a.distTo(b);  
        numTotalLigContactsPair = 0;
        
        
        numPairContacts = new Integer[Main.NUM_MOLECULE_PAIR_CONTACT_TYPES];
        minContactDistances = new Integer[Main.NUM_MOLECULE_PAIR_CONTACT_TYPES];
        contactAtomIndexInMoleculeA = new Integer[Main.NUM_MOLECULE_PAIR_CONTACT_TYPES];
        contactAtomIndexInMoleculeB = new Integer[Main.NUM_MOLECULE_PAIR_CONTACT_TYPES];
        
        Arrays.fill(numPairContacts, 0);
        Arrays.fill(minContactDistances, Integer.MAX_VALUE);
        Arrays.fill(contactAtomIndexInMoleculeA, -1);
        Arrays.fill(contactAtomIndexInMoleculeB, -1);
    }

    /*
     * @return A short String representation of the contact (molA<->molB).
     */
    public String shortStringRep() {
        return(this.getMolA().getUniqueString() + "<-->" + this.getMolB().getUniqueString());
    }

    @Override public String toString() {
        return("[MCI] Types: " + getResTypeStringA() + "/" + getResTypeStringB() + ", DsspNums"+ getDsspNumA() + "/" + getDsspNumB() + 
                ", chainIDs" + molA.getChainID() + "/" + molB.getChainID()) + " ResNums: " + molA.getPdbNum() + "/" + molB.getPdbNum() +
                ", iCodes: '" + molA.iCode + "'/'" + molB.iCode + "'";
    }
    
    
    /**
     * Whether one or both molecules are ligands.
     * @return 
     */
    public Boolean isLigandContact() {
        return molA.isLigand() || molB.isLigand();
    }
    
        /**
     * Whether one or both molecules are RNA.
     * @return 
     */
    public Boolean isRnaContact() {
        return molA.isRNA() || molB.isRNA();
    }


    /** Returns a string indicating whether molecule A is a protein residue ("PRT"), a ligand ("LIG"), an RNA residue ("RNA"), or something else ("OTH").
     * @return  a string indicating whether molecule B is a protein residue ("PRT"), a ligand ("LIG"), an RNA residue ("RNA"), or something else ("OTH").  
     */
    public String getResTypeStringA() {  if(this.molA.isAA()) { return("PRT"); } else if(this.molA.isLigand()) { return("LIG"); } else if(this.molA.isRNA()) { return("RNA"); } else { return("OTH"); } }
    /** Returns a string indicating whether molecule B is a protein residue ("PRT"), a ligand ("LIG"), an RNA residue ("RNA"), or something else ("OTH").
     * @return  a string indicating whether molecule B is a protein residue ("PRT"), a ligand ("LIG"), an RNA residue ("RNA"), or something else ("OTH").
     */
    public String getResTypeStringB() {  if(this.molB.isAA()) { return("PRT"); } else if(this.molB.isLigand()) { return("LIG"); } else if(this.molB.isRNA()) { return("RNA"); } else { return("OTH"); } }


    // getters required for writing the <pdbid>.geo file
    public Integer getDsspNumA() { return(this.molA.getDsspNum()); }
    public Integer getDsspNumB() { return(this.molB.getDsspNum()); }
    public Integer getPdbNumA() { return(this.molA.getPdbNum()); }
    public Integer getPdbNumB() { return(this.molB.getPdbNum()); }
    public String getName3A() { return(this.molA.getName3()); }
    public String getName3B() { return(this.molB.getName3()); }
    public String getName1A() { return(this.molA.getAAName1()); }
    public String getName1B() { return(this.molB.getAAName1()); }
    public Integer getAAIDResA() { return(this.molA.getInternalAAID()); }
    public Integer getAAIDResB() { return(this.molB.getInternalAAID()); }

    public Integer getCenterSphereRadiusResA() { return(this.molA.getSphereRadius()); }
    public Integer getCenterSphereRadiusResB() { return(this.molB.getSphereRadius()); }
    public Integer getMolPairDist() { return(dist); }

    public Integer getHB1Dist() { return(minContactDistances[HB]); }
    public Integer getHB2Dist() { return(minContactDistances[BH]); }

    public Integer getBBContactDist() { return(minContactDistances[BB]); }
    public Integer getBBContactAtomNumA() { return(contactAtomIndexInMoleculeA[BB] + 1); }   // The '+1' is required because geom_neo starts the atom index at '1' but we start it at '0'.
    public Integer getBBContactAtomNumB() { return(contactAtomIndexInMoleculeB[BB] + 1); }   //  This also means that we HAVE to init the atom index arrays in Main.calculateAtomContactsBetweenResidues()
                                                                                            //  with '-1'.  This return statement will add 1, and '0' means 'no contact' for geom_neo. :)
    public Integer getCBContactDist() { return(minContactDistances[CB]); }
    public Integer getCBContactAtomNumA() { return(contactAtomIndexInMoleculeA[CB] + 1); }
    public Integer getCBContactAtomNumB() { return(contactAtomIndexInMoleculeB[CB] + 1); }

    public Integer getBCContactDist() { return(minContactDistances[BC]); }
    public Integer getBCContactAtomNumA() { return(contactAtomIndexInMoleculeA[BC] + 1); }
    public Integer getBCContactAtomNumB() { return(contactAtomIndexInMoleculeB[BC] + 1); }

    public Integer getCCContactDist() { return(minContactDistances[CC]); }
    public Integer getCCContactAtomNumA() { return(contactAtomIndexInMoleculeA[CC] + 1); }
    public Integer getCCContactAtomNumB() { return(contactAtomIndexInMoleculeB[CC] + 1); }

    public Integer getBLContactDist() { return(minContactDistances[BL]); }
    public Integer getBLContactAtomNumA() { return(contactAtomIndexInMoleculeA[BL] + 1); }
    public Integer getBLContactAtomNumB() { return(contactAtomIndexInMoleculeB[BL] + 1); }

    public Integer getLBContactDist() { return(minContactDistances[LB]); }
    public Integer getLBContactAtomNumA() { return(contactAtomIndexInMoleculeA[LB] + 1); }
    public Integer getLBContactAtomNumB() { return(contactAtomIndexInMoleculeB[LB] + 1); }

    public Integer getCLContactDist() { return(minContactDistances[CL]); }
    public Integer getCLContactAtomNumA() { return(contactAtomIndexInMoleculeA[CL] + 1); }
    public Integer getCLContactAtomNumB() { return(contactAtomIndexInMoleculeB[CL] + 1); }

    public Integer getLCContactDist() { return(minContactDistances[LC]); }
    public Integer getLCContactAtomNumA() { return(contactAtomIndexInMoleculeA[LC] + 1); }
    public Integer getLCContactAtomNumB() { return(contactAtomIndexInMoleculeB[LC] + 1); }

    public Integer getLLContactDist() { return(minContactDistances[LL]); }
    public Integer getLLContactAtomNumA() { return(contactAtomIndexInMoleculeA[LL] + 1); }
    public Integer getLLContactAtomNumB() { return(contactAtomIndexInMoleculeB[LL] + 1); }
    
    public Integer getIVDWContactDist() { return(minContactDistances[IVDW]); }
    public Integer getIVDWContactAtomNumA() { return(contactAtomIndexInMoleculeA[IVDW] + 1); }
    public Integer getIVDWContactAtomNumB() { return(contactAtomIndexInMoleculeB[IVDW] + 1); }
    
    public Integer getBBHBContactDist() { return(minContactDistances[BBHB]); }
    public Integer getBBHBContactAtomNumA() { return(contactAtomIndexInMoleculeA[BBHB] + 1); }
    public Integer getBBHBContactAtomNumB() { return(contactAtomIndexInMoleculeB[BBHB] + 1); }
    
    public Integer getBBBHContactDist() { return(minContactDistances[BBBH]); }
    public Integer getBBBHContactAtomNumA() { return(contactAtomIndexInMoleculeA[BBBH] + 1); }
    public Integer getBBBHContactAtomNumB() { return(contactAtomIndexInMoleculeB[BBBH] + 1); }
    
    public Integer getBCHBContactDist() { return(minContactDistances[BCHB]); }
    public Integer getBCHBContactAtomNumA() { return(contactAtomIndexInMoleculeA[BCHB] + 1); }
    public Integer getBCHBContactAtomNumB() { return(contactAtomIndexInMoleculeB[BCHB] + 1); }
    
    public Integer getBCBHContactDist() { return(minContactDistances[BCBH]); }
    public Integer getBCBHContactAtomNumA() { return(contactAtomIndexInMoleculeA[BCBH] + 1); }
    public Integer getBCBHContactAtomNumB() { return(contactAtomIndexInMoleculeB[BCBH] + 1); }
    
    public Integer getCBHBContactDist() { return(minContactDistances[CBHB]); }
    public Integer getCBHBContactAtomNumA() { return(contactAtomIndexInMoleculeA[CBHB] + 1); }
    public Integer getCBHBContactAtomNumB() { return(contactAtomIndexInMoleculeB[CBHB] + 1); }
    
    public Integer getCBBHContactDist() { return(minContactDistances[CBBH]); }
    public Integer getCBBHContactAtomNumA() { return(contactAtomIndexInMoleculeA[CBBH] + 1); }
    public Integer getCBBHContactAtomNumB() { return(contactAtomIndexInMoleculeB[CBBH] + 1); }
    
    public Integer getCCHBContactDist() { return(minContactDistances[CCHB]); }
    public Integer getCCHBContactAtomNumA() { return(contactAtomIndexInMoleculeA[CCHB] + 1); }
    public Integer getCCHBContactAtomNumB() { return(contactAtomIndexInMoleculeB[CCHB] + 1); }
    
    public Integer getCCBHContactDist() { return(minContactDistances[CCBH]); }
    public Integer getCCBHContactAtomNumA() { return(contactAtomIndexInMoleculeA[CCBH] + 1); }
    public Integer getCCBHContactAtomNumB() { return(contactAtomIndexInMoleculeB[CCBH] + 1); }
    
    public Integer getNHPIContactDist() { return(minContactDistances[NHPI]); }
    public Integer getNHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[NHPI] + 1); }
    public Integer getNHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[NHPI] + 1); }
    
    public Integer getPINHContactDist() { return(minContactDistances[PINH]); }
    public Integer getPINHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PINH] + 1); }
    public Integer getPINHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PINH] + 1); }
    
    public Integer getCAHPIContactDist() { return(minContactDistances[CAHPI]); }
    public Integer getCAHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[CAHPI] + 1); }
    public Integer getCAHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[CAHPI] + 1); }
    
    public Integer getPICAHContactDist() { return(minContactDistances[PICAH]); }
    public Integer getPICAHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PICAH] + 1); }
    public Integer getPICAHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PICAH] + 1); }
    
    public Integer getCNHPIContactDist() { return(minContactDistances[CNHPI]); }
    public Integer getCNHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[CNHPI] + 1); }
    public Integer getCNHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[CNHPI] + 1); }
    
    public Integer getPICNHContactDist() { return(minContactDistances[PICNH]); }
    public Integer getPICNHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PICNH] + 1); }
    public Integer getPICNHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PICNH] + 1); }
    
    public Integer getSHPIContactDist() { return(minContactDistances[SHPI]); }
    public Integer getSHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[SHPI] + 1); }
    public Integer getSHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[SHPI] + 1); }
    
    public Integer getPISHContactDist() { return(minContactDistances[PISH]); }
    public Integer getPISHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PISH] + 1); }
    public Integer getPISHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PISH] + 1); }
    
    public Integer getXOHPIContactDist() { return(minContactDistances[XOHPI]); }
    public Integer getXOHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[XOHPI] + 1); }
    public Integer getXOHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[XOHPI] + 1); }
    
    public Integer getPIXOHContactDist() { return(minContactDistances[PIXOH]); }
    public Integer getPIXOHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PIXOH] + 1); }
    public Integer getPIXOHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PIXOH] + 1); }
    
    public Integer getPROCDHPIContactDist() { return(minContactDistances[PROCDHPI]); }
    public Integer getPROCDHPIContactAtomNumA() { return(contactAtomIndexInMoleculeA[PROCDHPI] + 1); }
    public Integer getPROCDHPIContactAtomNumB() { return(contactAtomIndexInMoleculeB[PROCDHPI] + 1); }
    
    public Integer getPIPROCDHContactDist() { return(minContactDistances[PIPROCDH]); }
    public Integer getPIPROCDHContactAtomNumA() { return(contactAtomIndexInMoleculeA[PIPROCDH] + 1); }
    public Integer getPIPROCDHContactAtomNumB() { return(contactAtomIndexInMoleculeB[PIPROCDH] + 1); }   
    
    public Integer getCCAHCOContactDist() { return(minContactDistances[CCAHCO]); }
    public Integer getCCAHCOContactAtomNumA() { return(contactAtomIndexInMoleculeA[CCAHCO] + 1); }
    public Integer getCCAHCOContactAtomNumB() { return(contactAtomIndexInMoleculeB[CCAHCO] + 1); }
    
    public Integer getCCOCAHContactDist() { return(minContactDistances[CCOCAH]); }
    public Integer getCCOCAHContactAtomNumA() { return(contactAtomIndexInMoleculeA[CCOCAH] + 1); }
    public Integer getCCOCAHContactAtomNumB() { return(contactAtomIndexInMoleculeB[CCOCAH] + 1); }
    
    public Integer getBCAHCOContactDist() { return(minContactDistances[BCAHCO]); }
    public Integer getBCAHCOContactAtomNumA() { return(contactAtomIndexInMoleculeA[BCAHCO] + 1); }
    public Integer getBCAHCOContactAtomNumB() { return(contactAtomIndexInMoleculeB[BCAHCO] + 1); }
    
    public Integer getBCOCAHContactDist() { return(minContactDistances[BCOCAH]); }
    public Integer getBCOCAHContactAtomNumA() { return(contactAtomIndexInMoleculeA[BCOCAH] + 1); }
    public Integer getBCOCAHContactAtomNumB() { return(contactAtomIndexInMoleculeB[BCOCAH] + 1); }
    
    public Integer getCTLContactDist() { return(minContactDistances[CTL]); }
    public Integer getCTLContactAtomNumA() { return(contactAtomIndexInMoleculeA[CTL] + 1); }
    public Integer getCTLContactAtomNumB() { return(contactAtomIndexInMoleculeB[CTL] + 1); }
    
    
// the getters for statistics follow
    public Integer getNumContactsTotal() { return(numPairContacts[TT]); }
    public Integer getNumLigContactsTotal() { return(numTotalLigContactsPair); }
    public Integer getNumRnaContactsTotal() { return(numTotalRnaContactsPair); }
    public Integer getNumContactsBB() { return(numPairContacts[BB]); }
    public Integer getNumContactsBC() { return(numPairContacts[BC]); }
    public Integer getNumContactsCB() { return(numPairContacts[CB]); }
    public Integer getNumContactsCC() { return(numPairContacts[CC]); }
    public Integer getNumContactsHB1() { return(numPairContacts[HB]); }
    public Integer getNumContactsHB2() { return(numPairContacts[BH]); }
    public Integer getNumContactsBL() { return(numPairContacts[BL]); }
    public Integer getNumContactsLB() { return(numPairContacts[LB]); }
    public Integer getNumContactsCL() { return(numPairContacts[CL]); }
    public Integer getNumContactsLC() { return(numPairContacts[LC]); }
    public Integer getNumContactsLL() { return(numPairContacts[LL]); }
    public Integer getNumContactsDisulfide() { return(numPairContacts[DISULFIDE]); }
    public Integer getNumContactsBBHB() { return(numPairContacts[BBHB]); }
    public Integer getNumContactsBBBH() { return(numPairContacts[BBBH]); }
    public Integer getNumContactsIVDW() { return(numPairContacts[IVDW]); }
    public Integer getNumContactsISS() { return(numPairContacts[ISS]); }
    public Integer getNumContactsIPI() { return(numPairContacts[IPI]); }
    public Integer getNumContactsISB() { return(numPairContacts[ISB]); }
    public Integer getNumContactsBCHB() { return(numPairContacts[BCHB]); }
    public Integer getNumContactsBCBH() { return(numPairContacts[BCBH]); }
    public Integer getNumContactsCBHB() { return(numPairContacts[CBHB]); }
    public Integer getNumContactsCBBH() { return(numPairContacts[CBBH]); }
    public Integer getNumContactsCCHB() { return(numPairContacts[CCHB]); }
    public Integer getNumContactsCCBH() { return(numPairContacts[CCBH]); }
    public Integer getNumContactsNHPI() { return(numPairContacts[NHPI]); }
    public Integer getNumContactsPINH() { return(numPairContacts[PINH]); }
    public Integer getNumContactsCAHPI() { return(numPairContacts[CAHPI]); }
    public Integer getNumContactsPICAH() { return(numPairContacts[PICAH]); }
    public Integer getNumContactsCNHPI() { return(numPairContacts[CNHPI]); }
    public Integer getNumContactsPICNH() { return(numPairContacts[PICNH]); }
    public Integer getNumContactsSHPI() { return(numPairContacts[SHPI]); }
    public Integer getNumContactsPISH() { return(numPairContacts[PISH]); }
    public Integer getNumContactsXOHPI() { return(numPairContacts[XOHPI]); }
    public Integer getNumContactsPIXOH() { return(numPairContacts[PIXOH]); }
    public Integer getNumContactsPROCDHPI() { return(numPairContacts[PROCDHPI]); }
    public Integer getNumContactsPIPROCDH() { return(numPairContacts[PIPROCDH]); }
    public Integer getNumContactsCCACOH() { return(numPairContacts[CCAHCO]); }
    public Integer getNumContactsCCOCAH() { return(numPairContacts[CCOCAH]); }
    public Integer getNumContactsBCACOH() { return(numPairContacts[BCAHCO]); }
    public Integer getNumContactsBCOCAH() { return(numPairContacts[BCOCAH]); }
    public Integer getNumContactsCTL() { return(numPairContacts[CTL]); }
    
    public ArrayList<String> getAtomAtomContactTypes() { return atomAtomContactType; }
    public ArrayList<Atom[]> getAtomAtomContacts() { return atomAtomContacts; }
    
    public Integer increaseContact(Integer contact, Integer amount) {
        // Setter increases every contact type by a specific amount
        // Should be used with contact letter
        // TT currently needs to be increased separately
        numPairContacts[contact] = numPairContacts[contact] + amount;  
        return numPairContacts[contact];
        
    }
    
    // DEBUG only
    public Molecule getMolA() { return(molA); }
    public Molecule getMolB() { return(molB); }
    
    /**
     * Determines whether the molecule pair described by this RCI is in any contact, i.e., has contacts of any type.
     * @return true if this RCI contains any kind of contact
     */
    public boolean describesAnyContact() {
        return(this.getNumContactsTotal() > 0 || this.getNumLigContactsTotal() > 0);
    }

    /**
     * Determines whether the molecule pair described by this RCI has enough contacts on atom levels to count as PPI.
     * For the pair to be in contact, at least two contacts of any kind or one disulphide bridge contact must present.
     * @return true if this RCI contains a PPI contact
     */
    public boolean describesPPIContact() {
         return(this.getNumContactsTotal() > 1 || this.getNumLigContactsTotal() > 1 || this.getNumContactsISS() > 0);
    }

    // no setters needed, this class is only a nice wrapper for the ugly arrays

}
