/*
 * This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.
 *
 * Copyright Tim Schäfer 2016. VPLG is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */


package proteinstructure;

import java.util.ArrayList;
import java.util.List;

/**
 * Models a binding site in a certain protein structure. Binding site information is parsed from the SITE entries of a PDB file.
 * @author spirit
 */
public class BindingSite {
    
    public BindingSite() {
        this.siteResidues = new ArrayList<>();
        this.ligandResidues = new ArrayList<>();
        this.siteName = "";
        this.modelID = null;
        this.numResDeclared = null;
        this.residueInfos = new ArrayList<>();
    }
    
    public BindingSite(String siteName) {
        this.siteResidues = new ArrayList<>();
        this.ligandResidues = new ArrayList<>();
        this.siteName = siteName;
        this.modelID = null;
        this.numResDeclared = null;
        this.residueInfos = new ArrayList<>();
    }
    
    protected List<Residue> siteResidues;
    protected List<Residue> ligandResidues;
    protected String modelID;
    /** The number of residues declared in the PDB SITE line for this site. */
    protected Integer numResDeclared;
    
    /** Site description from REMARK 800 line of PDB file */
    protected String description;
    /** Site evidence code from REMARK 800 line of PDB file */
    protected String evidenceCode;
    protected List<String[]> residueInfos;

    public List<String[]> getResidueInfos() {
        return residueInfos;
    }

    public void setResidueInfos(List<String[]> residueInfos) {
        this.residueInfos = residueInfos;
    }
    
    public void addResidueInfos(List<String[]> residueInfos) {
        for(String[] r : residueInfos) {
            this.residueInfos.add(r);
        }
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getEvidenceCode() {
        return evidenceCode;
    }

    public void setEvidenceCode(String evidenceCode) {
        this.evidenceCode = evidenceCode;
    }
    
    /** The name assigned in the SITE line of the PDB file. */
    protected String siteName;

    public List<Residue> getSiteResidues() {
        return siteResidues;
    }
    
    public String getSiteName() {
        return this.siteName;
    }
    
    public void addSiteResidue(Residue r) {
        this.siteResidues.add(r);
    }
    
    public void setNumResDeclared(Integer i) {
        this.numResDeclared = i;
    }
    
    public void addLigandResidue(Residue r) {
        this.ligandResidues.add(r);
    }

    public void setSiteResidues(List<Residue> siteResidues) {
        this.siteResidues = siteResidues;
    }

    public List<Residue> getLigandResidues() {
        return ligandResidues;
    }

    public void setLigandResidues(List<Residue> ligandResidues) {
        this.ligandResidues = ligandResidues;
    }
    
    public void setModelID(String modelID) {
        this.modelID = modelID;
    }
    
    public Integer getNumResiduesDeclared() {
        return this.numResDeclared;
    }
    
    @Override
    public String toString() {
        return this.siteName + " (" + this.getDescription() + "), size " + this.getResidueInfos().size() + "";
    }
}