/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package proteinstructure;

// imports

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import proteingraphs.Position3D;
import settings.Settings;


/**
 * Represents an Atom, based on information in a PDB file.
 * 
 * @author ts
 */
public class Atom implements java.io.Serializable {

    // declare class vars
    private String name = null;                 // atom name from PDB file, spaces included (they are important and differentiate atoms, e.g. 'CA  '= C alpha from ' CA '=calcium)
    private String chemSym = null;              // chemical symbol of the atom (periodic table notation, extracted from pdb file)
    private String chainID = null;
    private String modelID = null;
    private Model model = null;
    private Integer pdbAtomNumber = null;       // atom number from pdb file
    public Molecule molecule = null;             // Molecule this Atom belongs to
    private Integer type = null;               // atom type:  0=AA, 1=Ligand, 2=Ignored HETATM (e.g. 'DOD'-molecule atoms) 3=Ignored ATOM (e.g. H, Q)
    private Integer pdbResNum = null;
    private Integer dsspResNum = null;
    private Integer coordX = null;              // 3D coordinate X from pdb file, converted to 10th part Angstroem
    private Integer coordY = null;
    private Integer coordZ = null;
    private Float ocoordX = null;              // 3D coordinate X from pdb file, original
    private Float ocoordY = null;
    private Float ocoordZ = null;
    private Integer pdbLineNum = null;
    private Chain chain = null;
    private String altLoc = null;
    
    public static final ArrayList<String> BACKBONE_ATOM_NAMES = new ArrayList<>(Arrays.asList("N", "CA", "C", "O"));  // PDB names of atoms that are backbone atoms
    
    public static final Integer ATOMTYPE_AA = 0;
    public static final Integer ATOMTYPE_LIGAND = 1;
    public static final Integer ATOMTYPE_IGNORED_LIGAND = 2;
    public static final Integer ATOMTYPE_RNA = 3;
    public static final Integer ATOMTYPE_IGNORED_ATOM = 4;

    /**
     * Getter for PDB alternate location identifier.
     * @return the altLoc string
     */
    public String getAltLoc() {
        return altLoc;
    }

    /**
     * Sets the alternate location identifier (from PDB file).
     * @param altLoc the altLoc string (one character, usually " ")
     */
    public void setAltLoc(String altLoc) {
        this.altLoc = altLoc;
    }

    public Boolean isLigandAtom() { return this.type == ATOMTYPE_LIGAND; }
    public Boolean isProteinAtom() { return this.type == ATOMTYPE_AA; }
    public Boolean isOtherAtom() { return this.type == ATOMTYPE_IGNORED_LIGAND; }
    public Boolean isRnaAtom() { return this.type == ATOMTYPE_RNA; }
    
    
    /**
     * Returns the distance from this atom to atom 'a'. It uses the atom centers to calculate
     * the distance, so you have to take care of the collision sphere size yourself.
     * @param a the other Atom
     * @return the euclidian distance, rounded to an Integer
     */
    public Integer distToAtom(Atom a) {
        Integer di;
        Integer deltaOriginal;
        Integer deltaRounded;
        //di = distToPoint(a.getCoordX(), a.getCoordY(), a.getCoordZ());
        // fjg: use the original coordinates at least to calculate the distance and then round it
        di = distToPointFloat(a.getoCoordX(), a.getoCoordY(), a.getoCoordZ());
        deltaOriginal = di;
        deltaRounded = distToPoint(a.getCoordX(), a.getCoordY(), a.getCoordZ());
        if((deltaRounded == 20 && deltaOriginal == 19) || (deltaRounded == 19 && deltaOriginal == 20 ) || (deltaRounded == 21 && deltaOriginal == 20) || (deltaRounded == 20 && deltaOriginal == 21)) {
            System.out.println("Distance between atoms: " + deltaOriginal + " " + deltaRounded);};
        
        if(Settings.getBoolean("PTGLgraphComputation_B_contact_debug_dysfunct")) {
            if(this.isCalphaAtom() && a.isCalphaAtom()) {
                System.out.println("Distance between C-alpha atoms " + this.pdbAtomNumber + " of " + this.getPdbResNum() + " and " + a.pdbAtomNumber + " of " + a.getPdbResNum() + " is " + di + " (-- before sqrt -> due to change of function not given).");
                System.out.println(this.getCoordString() + "/" + a.getCoordString());
            }            
        }
        
        return(di);
    }
    
    /**
     * Returns the distance from this atom to a point in 3D.
     * @param dx X coordinate as 10th of Angström in int
     * @param dy Y coordinate as 10th of Angström in int
     * @param dz Z coordinate as 10th of Angström in int
     * @return the euclidian distance, rounded to an int
     */
    public Integer distToPoint(int dx, int dy, int dz) {
        Double dd = 0.0;
        Integer di;
        
        dd += (coordX - dx) * (coordX - dx);
        dd += (coordY - dy) * (coordY - dy);
        dd += (coordZ - dz) * (coordZ - dz);

        // di = (int)Math.sqrt(dd);
        // jnw: lets round instead of truncate the result
        di = (int)Math.round(Math.sqrt(dd));
        
        return(di);
    }
    
        
    /**
     * Returns the distance from this atom to a point in 3D.
     * @param dx X coordinate as 10th of Angström in float
     * @param dy Y coordinate as 10th of Angström in float
     * @param dz Z coordinate as 10th of Angström in float
     * @return the euclidian distance, rounded to an int
     */
    public Integer distToPointFloat(float dx, float dy, float dz) {
        Double dd = 0.0;
        Integer di;
        
        dd += (coordX - dx) * (coordX - dx);
        dd += (coordY - dy) * (coordY - dy);
        dd += (coordZ - dz) * (coordZ - dz);

        // di = (int)Math.sqrt(dd);
        // jnw: lets round instead of truncate the result
        di = (int)Math.round(Math.sqrt(dd));
        
        return(di);
    }
    
    
    /**
     * Returns the distance from this atom to atom 'a'. It uses the atom centers to calculate
     * the distance, so you have to take care of the collision sphere size yourself.
     * @param a the other Atom
     * @return the euclidian distance, rounded to an Integer
     */
    @Deprecated public Integer distToAtomOld(Atom a) {
        Double distDouble = .0;
        Integer distInt, dx, dy, dz, dd = 0;

        dx = this.getCoordX() - a.getCoordX();
        dy = this.getCoordY() - a.getCoordY();
        dz = this.getCoordZ() - a.getCoordZ();

        dd = dx * dx + dy * dy + dz * dz;
        distDouble = Math.sqrt(dd);
        distInt = Integer.valueOf((int)Math.round(distDouble));
        
        if(Settings.getBoolean("PTGLgraphComputation_B_contact_debug_dysfunct")) {
            if(this.isCalphaAtom() && a.isCalphaAtom()) {
                System.out.println("Distance between C-alpha atoms " + this.pdbAtomNumber + " of " + this.getPdbResNum() + " and " + a.pdbAtomNumber + " of " + a.getPdbResNum() + " is " + distInt + " (" + dd + " before sqrt, " + distDouble + " as Double).");
            }            
        }
        
        return(distInt);
    }
    
    
    /**
     * Determines whether this is a C alpha atom. This is determined from the ATOM NAME entry of the PDB file.
     * @return true if this is a C alpha atom, false otherwise.
     */
    public boolean isCalphaAtom() {
        return(this.isProteinAtom() && this.name.equals(" CA "));
    }



    /**
     * Checks whether a contact (vdW radius overlap) exists to another atom. Returns the distance
     * to the atom (a positive Integer) if a contact exists, -1 otherwise. (This behavior saves
     * us from having to call the function twice during atom contact calculation, which is important
     * since we don't want to do useless Math.sqrt() operations (for performance reasons).
     */
    public Boolean atomContactTo(Atom a) {


        Integer atomRadiusThis;
        Integer atomRadiusOther;
        

        Integer radProt = Settings.getInteger("PTGLgraphComputation_I_aa_atom_radius");
        Integer radLig = Settings.getInteger("PTGLgraphComputation_I_lig_atom_radius");
        Integer radRna = Settings.getInteger("PTGLgraphComputation_I_rna_atom_radius");


        if(this.isProteinAtom()) {
            atomRadiusThis = radProt;
        }
        else if (this.isRnaAtom()) {
            atomRadiusThis = radRna;
        }
        else {
            atomRadiusThis = radLig;
        }
        
        if(a.isProteinAtom()) {
            atomRadiusOther = radProt;
        }
        else if (a.isRnaAtom()) {
            atomRadiusOther = radRna;
        }
        else {
            atomRadiusOther = radLig;
        }


        Integer dist = this.distToAtom(a);
        Integer maxDist = atomRadiusThis + atomRadiusOther;

        //if(dist < 0) {
        //    System.err.println("ERROR: Distance of atoms " + this.getPdbAtomNum() + " and " + a.getPdbAtomNum() + " is " + dist + ", but should be > 0.");
        //    System.exit(1);
        //}
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
            if (dist < maxDist) {
                System.out.println("   [DEBUG LV 2] Atom " + this.getPdbAtomNum() + " " +
                    this.getCoordString() + " and atom " + 
                    a.getPdbAtomNum() + " " + a.getCoordString() + " have distance " + dist + " => contact!");
            } else {
                System.out.println("   [DEBUG LV 2] Atom " + this.getPdbAtomNum() + " " +
                    this.getCoordString() + " and atom " + 
                    a.getPdbAtomNum() + " " + a.getCoordString() + " have distance " + dist);
            }
        }
            

        if( dist < maxDist) {
            // Contact!
            //System.out.println("        ++++ CONTACT between atoms " + this.pdbAtomNumber + " and " + a.getPdbAtomNum() + " in dist " + dist + " (maxDist is " + maxDist + ").");
            return(true);
        }
        else {
            // No contact
            //System.out.println("        ---- No contact between atoms " + this.pdbAtomNumber + " and " + a.getPdbAtomNum() + " in dist " + dist + " (maxDist is " + maxDist + ").");
            return(false);
        }
    }

    @Override public String toString() {
        return("[Atom] #" + this.pdbAtomNumber + " NAME='" + this.name + "' CS='" + this.chemSym + "' Type=" + this.type + " Chain=" + this.chainID + " ResDssp=" + this.getMolecule().getDsspNum() + " ResPDB=" + this.getMolecule().getUniquePDBName() + " Coords=" + getCoordString() + " AltLoc='" + this.altLoc + "'.");
    }


    /**
     * Checks whether a vdW contact exists to another atom (vdW radii overlap).
     * In contrast to atomContactTo() this method considers different vdW radii
     * for different atoms. This method is used for the calculation of the
     * alternative contact model.
     * @param a is the atom to which a contact is checked.
     * @return True if contact exists, otherwise False.
     */
    public Boolean vdwAtomContactTo (Atom a) {
        HashMap<String, Double> vdwRadii = new HashMap<>();
        vdwRadii.put("H", 12.0);
        vdwRadii.put("HE", 14.0);
        vdwRadii.put("LI", 18.1);
        vdwRadii.put("BE", 19.8);
        vdwRadii.put("B", 19.1);
        vdwRadii.put("C", 17.0);
        vdwRadii.put("N", 15.5);
        vdwRadii.put("O", 15.2);
        vdwRadii.put("F", 14.7);
        vdwRadii.put("NE", 15.4);
        vdwRadii.put("NA", 22.7);
        vdwRadii.put("MG", 17.3);
        vdwRadii.put("AL", 22.5);
        vdwRadii.put("SI", 22.2);
        vdwRadii.put("P", 18.0);
        vdwRadii.put("S", 18.0);
        vdwRadii.put("CL", 17.5);
        vdwRadii.put("AR", 17.6);
        vdwRadii.put("K", 27.5);
        vdwRadii.put("CA", 26.2);
        vdwRadii.put("SC", 25.8);
        vdwRadii.put("TI", 24.6);
        vdwRadii.put("V", 24.2);
        vdwRadii.put("CR", 24.5);
        vdwRadii.put("MN", 24.5);
        vdwRadii.put("FE", 24.4);
        vdwRadii.put("CO", 24.0);
        vdwRadii.put("NI", 16.3);
        vdwRadii.put("CU", 14.0);
        vdwRadii.put("ZN", 13.9);
        vdwRadii.put("GA", 18.7);
        vdwRadii.put("GE", 22.9);
        vdwRadii.put("AS", 18.5);
        vdwRadii.put("SE", 19.0);
        vdwRadii.put("BR", 18.3);
        vdwRadii.put("KR", 20.2);
        vdwRadii.put("RB", 32.1);
        vdwRadii.put("SR", 28.4);
        vdwRadii.put("Y", 27.5);
        vdwRadii.put("ZR", 25.2);
        vdwRadii.put("NB", 25.6);
        vdwRadii.put("MO", 24.5);
        vdwRadii.put("TC", 24.4);
        vdwRadii.put("RU", 24.6);
        vdwRadii.put("RH", 24.4);
        vdwRadii.put("PD", 16.3);
        vdwRadii.put("AG", 17.2);
        vdwRadii.put("CD", 16.2);
        vdwRadii.put("IN", 19.3);
        vdwRadii.put("SN", 21.7);
        vdwRadii.put("SB", 22.0);
        vdwRadii.put("TE", 20.0);
        vdwRadii.put("I", 19.8);
        vdwRadii.put("XE", 21.6);
        vdwRadii.put("CS", 34.8);
        vdwRadii.put("BA", 30.3);
        vdwRadii.put("LA", 29.8);
        vdwRadii.put("CE", 28.8);
        vdwRadii.put("PR", 29.2);
        vdwRadii.put("ND", 29.5);
        vdwRadii.put("SM", 29.0);
        vdwRadii.put("EU", 28.7);
        vdwRadii.put("GD", 28.3);
        vdwRadii.put("TB", 27.9);
        vdwRadii.put("DY", 28.7);
        vdwRadii.put("HO", 28.1);
        vdwRadii.put("ER", 28.3);
        vdwRadii.put("TM", 27.9);
        vdwRadii.put("YB", 28.0);
        vdwRadii.put("LU", 27.4);
        vdwRadii.put("HF", 26.3);
        vdwRadii.put("TA", 25.3);
        vdwRadii.put("W", 25.7);
        vdwRadii.put("RE", 24.9);
        vdwRadii.put("OS", 24.8);
        vdwRadii.put("IR", 24.1);
        vdwRadii.put("PT", 17.2);
        vdwRadii.put("AU", 16.6);
        vdwRadii.put("HG", 17.0);
        vdwRadii.put("TL", 19.6);
        vdwRadii.put("PB", 20.2);
        vdwRadii.put("BI", 23.0);
        vdwRadii.put("AC", 28.0);
        vdwRadii.put("TH", 29.3);
        vdwRadii.put("PA", 28.8);
        vdwRadii.put("U", 18.6);
        vdwRadii.put("NP", 28.2);
        vdwRadii.put("PU", 28.1);
        vdwRadii.put("AM", 28.3);
        vdwRadii.put("CM", 30.5);
        vdwRadii.put("BK", 34.0);
        vdwRadii.put("CF", 30.5);
        vdwRadii.put("ES", 27.0);
        
        
        Double atomRadiusThis;
        Double atomRadiusOther;
        
        //Double radLig = Settings.getInteger("PTGLgraphComputation_I_lig_atom_radius").doubleValue();

        if(this.isLigandAtom()) {
            atomRadiusThis = vdwRadii.get(this.chemSym.replaceAll("\\s+",""));
            if(atomRadiusThis == null) {
                atomRadiusThis = 12.0;
            }
        }
        else {
            atomRadiusThis = vdwRadii.get(this.chemSym.replaceAll("\\s+","")); //replaceAll is needed to delete whitespace in front of the chemSym; otherwise you cannot look it up in the hashmap
            if (atomRadiusThis == null) {
                atomRadiusThis = 12.0;
            }
        }
        
        if(a.isLigandAtom()) {
            atomRadiusOther = vdwRadii.get(a.chemSym.replaceAll("\\s+",""));
            if(atomRadiusOther == null) {
                atomRadiusOther = 12.0;
            }
        }
        else {
            atomRadiusOther = vdwRadii.get(a.chemSym.replaceAll("\\s+",""));
            if (atomRadiusOther == null) {
                atomRadiusOther = 12.0;
            }
        }


        Double dist = this.distToAtom(a).doubleValue();
        Double maxDist = atomRadiusThis + atomRadiusOther;

        //TODO: - check whether dist is <0? Why has it been removed in atomContactTo()?
        if( dist < maxDist) {
            // Contact!
            //System.out.println("        ++++ CONTACT between atoms " + this.pdbAtomNumber + " and " + a.getPdbAtomNum() + " in dist " + dist + " (maxDist is " + maxDist + ").");
            return(true);
        }
        else {
            // No contact
            //System.out.println("        ---- No contact between atoms " + this.pdbAtomNumber + " and " + a.getPdbAtomNum() + " in dist " + dist + " (maxDist is " + maxDist + ").");
            return(false);
        }
    }

    
    /**
     * Checks whether the angle between acceptor-donor and acceptor-acceptor antecedent
     * in a possible hydrogen bond is satisfied.
     * This is the case when the angle is greater than 90 degrees.
     * For more information see: doi:10.1006/jmbi.1994.1334
     * @param a is the acceptor atom.
     * @param aa is the acceptor antecedent atom.
     * @return True if the angle is greater than 90 degrees, otherwise False.
     */
    public Boolean hbondAtomAngleBetween(Atom a, Atom aa) {
    
        Integer d1x = this.getCoordX() - a.getCoordX();
        Integer d1y = this.getCoordY() - a.getCoordY();
        Integer d1z = this.getCoordZ() - a.getCoordZ();

        Integer d2x = aa.getCoordX() - a.getCoordX();
        Integer d2y = aa.getCoordY() - a.getCoordY();
        Integer d2z = aa.getCoordZ() - a.getCoordZ();

        Double angle = 0.0;
        angle = Math.acos((d1x * d2x + d1y * d2y + d1z * d2z) / (Math.sqrt(Math.pow(d1x,2) + Math.pow(d1y, 2) + Math.pow(d1z, 2)) * Math.sqrt(Math.pow(d2x, 2) + Math.pow(d2y, 2) + Math.pow(d2z, 2)))) * (180/Math.PI);
        
        if (angle > 90) {
            return(true);
        }

        return(false);
    }
    
    
    /**
     * Compares two Atoms via their PDB atom number.
     * @param other the other atom
     * @return  true if they are the same, false otherwise
     */
    public Boolean equalsAtom(Atom other) {
        return(this.pdbAtomNumber == other.pdbAtomNumber);
    }

    // getters
    public String getAtomName() { return(name); }
    public String getAtomShortName() { return(name.trim()); }
    public String getChemSym() { return(chemSym); }
    public String getChainID() { return(chainID); }
    public String getModelID() { return(modelID); }
    public Chain getChain() { return(chain); }
    public Model getModel() { return(model); }
    public Integer getPdbAtomNum() { return(pdbAtomNumber); }
    public Integer getCoordX() { return(coordX); }
    public Integer getCoordY() { return(coordY); }
    public Integer getCoordZ() { return(coordZ); }
    public Float getoCoordX() { return(ocoordX); }
    public Float getoCoordY() { return(ocoordY); }
    public Float getoCoordZ() { return(ocoordZ); }
    public Integer getPdbLineNum() { return(pdbLineNum); }
    public Molecule getMolecule() { return(molecule); }
    public Integer getPdbResNum() { return(pdbResNum); }
    public Integer getDsspResNum() { return(dsspResNum); }
    public Integer getAtomType() { return(type); }
    public String getCoordString() { return("(" + coordX + "," + coordY + "," + coordZ + ")"); }

    // setters
    public void setAtomName(String s) { name = s; }
    public void setChemSym(String s) { chemSym = s; }
    public void setChainID(String s) { chainID = s; }
    public void setModelID(String s) { modelID = s; }
    public void setPdbAtomNum(Integer i) { pdbAtomNumber = i; }
    public void setCoordX(Integer i) { coordX = i; }
    public void setCoordY(Integer i) { coordY = i; }
    public void setCoordZ(Integer i) { coordZ = i; }
    public void setoCoordX(Float i) { ocoordX = i; }
    public void setoCoordY(Float i) { ocoordY = i; }
    public void setoCoordZ(Float i) { ocoordZ = i; }
    public void setPdbLineNum(Integer i) { pdbLineNum = i; }
    public void setMolecule(Molecule m) { molecule = m; }
    public void setPdbResNum(Integer i) { pdbResNum = i; }
    public void setDsspResNum(Integer i) { dsspResNum = i; }
    public void setAtomtype(Integer i) { type = i; }
    public void setChain(Chain c) { chain = c; }
    public void setModel(Model m) { model = m; }
    
    public Position3D getPosition3D() {
        return new Position3D(coordX / 10.0f, coordY  / 10.0f, coordZ / 10.0f);
    }

}
