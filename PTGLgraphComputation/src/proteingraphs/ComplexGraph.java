/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Marcus Kessler 2013. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author Marcus Kessler
 * modified by TS
 * modified by fjg
 */
package proteingraphs;

import algorithms.AgglomerativeClustering;
import algorithms.ConsensusTree;
import datastructures.ClusteringResult;
import graphdrawing.PageLayout;
import graphdrawing.DrawTools;
import graphdrawing.DrawResult;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.RenderingHints;
import java.awt.geom.Arc2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.AffineTransform;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.spargel.datastructures.UAdjListGraph;
import net.sourceforge.spargel.writers.GMLWriter;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import io.DBManager;
import io.IO;
import io.FileParser;
import java.io.ByteArrayOutputStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Objects;
import java.util.Set;
import settings.Settings;
import proteinstructure.Chain;
import proteinstructure.Molecule;
import proteinstructure.Residue;
import proteinstructure.SSE;
import tools.DP;
import tools.TextTools;
import proteingraphs.ComplexGraphEdgeWeightTypes.EdgeWeightType;

/**
 *
 * @author marcus
 */
public class ComplexGraph extends UAdjListGraph {

    private final int numberChains;  // Number of chains
    private static List<Chain> allChains; // list of all chains
    public Map<Edge, String[]> chainNamesInEdge;
    public Map<Edge, Integer> numHelixHelixInteractionsMap;
    public Map<Edge, Integer> numHelixStrandInteractionsMap;
    public Map<Edge, Integer> numHelixCoilInteractionsMap;
    public Map<Edge, Integer> numHelixLigandInteractionsMap;
    public Map<Edge, Integer> numStrandStrandInteractionsMap;
    public Map<Edge, Integer> numStrandCoilInteractionsMap;
    public Map<Edge, Integer> numStrandLigandInteractionsMap;
    public Map<Edge, Integer> numCoilCoilInteractionsMap;
    public Map<Edge, Integer> numCoilLigandInteractionsMap;
    public Map<Edge, Integer> numLigandLigandInteractionsMap;
    public Map<Edge, Integer> numDisulfidesMap;
    
    private Map<EdgeWeightType, Map<Edge, BigDecimal>> mapWeightNamesToMapEdgeValues = new HashMap<>();  // <name of weight / normalization, <edge, value>>
    private BigDecimal minimumMultiplicativeLengthNormalizedEdgeWeight;  // the smallest normalized edge weight: used for the lucid normalized edge weights

    public Map<Vertex, String> proteinNodeMap;
    private Map<Vertex, String> labelNodeMap;  // for now, used for labels in interactive assembly prediction and output newick string
    public Map<Vertex, String> molMap;  // contains for each vertex (= protein chain) the corresponding molecule name
    public Map<Vertex, Double> chainLengthMap;  // used for the GML file output
    public Map<Vertex, Double> chainGyradiusMap;  // used for the GML file output
    public Map<Vertex, String> chainTypeMap;
    private Map<Integer, Vertex> mapVertexIdVertex = new HashMap<>();
    public Map<List<Integer>, Integer> numSSEContacts;
    public Map<List<Integer>, List<String>> numSSEContactChainNames;
    
    private final Boolean createContactInfo;
    private ArrayList<String> contactInfo;

    private Integer[][] numChainInteractions;
    private Integer[][] numChainInteractionsCTL;
    private ArrayList<Integer> chainIDList;
    private BigDecimal settingFactor;
    private HashMap<String, Integer> mapChainIdToLength;  // used for the normalized edge weights
    private HashMap<String, Double> mapChainIdToGyradius;  // used for the normalized edge weights
    private Integer[][] homologueChains;
    private final Set<String> molIDs;  // Contains the mol IDs for all chains. Used to get number of mol IDs (= size)
    private final String[] chainResAASeq;
    public Integer neglectedEdges;
    public Boolean includeRna = Settings.getBoolean("PTGLgraphComputation_B_include_rna");
    

    /**
     * The RCSB PDB id this graph is based on.
     */
    private final String pdbid;
    
    private final static String CLASS_TAG = "CG";
    
    /**
     * Constructor.
     * @param pdbid RSCB PDB ID
     * @param chains
     * @param resContacts
     * @param createConInfo Whether contact info for writing of GML should be created. Usually pass PTGLgraphComputation_B_writeComplexContactCSV setting.
     */
    public ComplexGraph(String pdbid, List<Chain> chains, List<MolContactInfo> resContacts, Boolean createConInfo) {
        this.pdbid = pdbid;
        allChains = chains;
        numberChains = chains.size();
        createContactInfo = createConInfo;
        
        numHelixHelixInteractionsMap = createEdgeMap();
        numHelixStrandInteractionsMap = createEdgeMap();
        numHelixCoilInteractionsMap = createEdgeMap();
        numHelixLigandInteractionsMap = createEdgeMap();
        numStrandStrandInteractionsMap = createEdgeMap();
        numStrandCoilInteractionsMap = createEdgeMap();
        numStrandLigandInteractionsMap = createEdgeMap();
        numCoilCoilInteractionsMap = createEdgeMap();
        numCoilLigandInteractionsMap = createEdgeMap();
        numLigandLigandInteractionsMap = createEdgeMap();
        
        for (EdgeWeightType type : EdgeWeightType.values()) {
            mapWeightNamesToMapEdgeValues.put(type, createEdgeMap());
        }
        
        numDisulfidesMap = createEdgeMap();
        chainNamesInEdge = createEdgeMap();
        
        proteinNodeMap = createVertexMap();
        labelNodeMap = createVertexMap();
        molMap = createVertexMap();
        chainLengthMap = createVertexMap();
        chainGyradiusMap = createVertexMap();
        chainTypeMap = createVertexMap();
        
        molIDs = new HashSet<>();
        numSSEContacts = new HashMap<>();
        numSSEContactChainNames = new HashMap<>();
        mapChainIdToLength = new HashMap<>();
	mapChainIdToGyradius = new HashMap<>();
        

        neglectedEdges = 0;
        
        chainResAASeq = new String[numberChains];
        
        // preprocess chains if required
        List<Chain> preprocessedChains = new ArrayList<>();
        if (Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")) {
            for (Chain tmpChain : allChains) {
                Chain newChain = new Chain(tmpChain.getPdbChainID());
                if (tmpChain.getMoleculeType().equals("polyribonucleotide") && includeRna) {
                    for (Molecule tmpMol : tmpChain.getAllRnaResidues()) {
                        // lig off, rna on
                        newChain.addMolecule(tmpMol);
                    }
                }
                for (Molecule tmpMol : tmpChain.getAllAAResidues()) {
                    // lig off, rna off/on
                    newChain.addMolecule(tmpMol);
                }
                preprocessedChains.add(newChain);
            }
        } else {
            if (includeRna) {
                // lig on, rna on
                preprocessedChains = allChains;
            } else {
                for (Chain tmpChain : allChains) {
                    Chain newChain = new Chain(tmpChain.getPdbChainID());
                    if (tmpChain.getMoleculeType().equals("polyribonucleotide")) {
                        // lig on, rna off
                        continue;
                    } else {
                        for (Molecule tmpMol : tmpChain.getMolecules()) {
                            // lig on, rna off
                            newChain.addMolecule(tmpMol);
                        }
                    }
                    preprocessedChains.add(newChain);
                }
            }
        }

        createVertices(preprocessedChains);
        createHomologueChainsMatrix(preprocessedChains);
        
        // preprocess res contacts if required
        List<MolContactInfo> preprocessedResContacts;
        if (Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")) {
            preprocessedResContacts = new ArrayList<>();
            for (MolContactInfo tmpMci : resContacts) {
                if (tmpMci.isRnaContact() && ! includeRna) {
                    // lig off, rna off
                    continue;
                } else {
                    if (! tmpMci.isLigandContact()) {
                        // lig off, rna on
                        preprocessedResContacts.add(tmpMci);
                    }
                }  
            }
        } else {
            if (includeRna) {
                // lig on, rna on
                preprocessedResContacts = resContacts;
            } else {
                preprocessedResContacts = new ArrayList<>();
                for (MolContactInfo tmpMci : resContacts) {
                    if (tmpMci.isRnaContact()) {
                        // lig on, rna off
                        continue;
                    } else {
                        preprocessedResContacts.add(tmpMci);
                    }
                }
            }
        }
        
        calculateNumChainInteractions(preprocessedResContacts);
        createEdges(preprocessedResContacts);
        
        assemblyPrediction();
    }
    
    
    /**
     * Creates the vertices, fills corresponding maps and fills amino acid sequence.
     * @param chains All chains of this complex graph.
     */
    private void createVertices(List<Chain> chains) {
        Vertex v;
        int curNum = 0;
        String labelListAsString = Settings.get("PTGLgraphComputation_S_assembly_prediction_leaf_labels");
        String[] labelList = (labelListAsString.length() == 0 ? new String[0] : labelListAsString.split(","));  // split would return [""] for empty String so care for that
        
        for(Integer i = 0; i < chains.size(); i++) {
            Chain tmpChain = chains.get(i);
            for (Integer j = 0; j < allChains.size(); j++) {    // loop through chains to find matches
                if (tmpChain.getPdbChainID() == allChains.get(j).getPdbChainID()) {     // if the chains match, check their molecule type
                    if (allChains.get(j).getMoleculeType().contains("polyribonucleotide") && ! includeRna) {  // hint: if inclusion of RNA is turned off, RNA/DNA hybrids are not included either
                        continue;
                    }
                    else {
                        v = createVertex();
             
                        proteinNodeMap.put(v, tmpChain.getPdbChainID());
                        labelNodeMap.put(v, labelList.length > curNum ? labelList[curNum] : tmpChain.getPdbChainID());
                        molMap.put(v, FileParser.getMetaInfo(pdbid, tmpChain.getPdbChainID()).getMolName());  // get the mol name from the ProtMetaInfo
                        chainLengthMap.put(v, (double) tmpChain.getLength(! Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")));
                        chainGyradiusMap.put(v, tmpChain.getChainRadiusOfGyration("mass", ! Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")));

                        chainTypeMap.put(v, allChains.get(j).getMoleculeType());
                        mapVertexIdVertex.put(Integer.parseInt(v.toString()), v);

                        molIDs.add(FileParser.getMetaInfo(pdbid, tmpChain.getPdbChainID()).getMolName());

                        int tmpLength = tmpChain.getLength(! Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands"));
                        mapChainIdToLength.put(tmpChain.getPdbChainID(), tmpLength);
			double tmpGyradius = tmpChain.getChainRadiusOfGyration("mass", ! Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands"));
                        mapChainIdToGyradius.put(tmpChain.getPdbChainID(), tmpGyradius);

                        // get AA sequence string for each chainName
                        for(Residue resi : tmpChain.getAllAAResidues()) {

                            if ( ! Settings.get("PTGLgraphComputation_S_ligAACode").equals(resi.getAAName1())) {  // Skip ligands to preserve sequence identity. What to do with "_B_", "_Z_", "_X_" (B,Z,X)?
                                if (chainResAASeq[i] != null) {
                                    chainResAASeq[i] = chainResAASeq[i] + resi.getAAName1();
                                } else {
                                    chainResAASeq[i] = resi.getAAName1();
                                }
                            }
                        }
                        curNum++;
                    }
                }
            } 
        }
    }
    
    
    /**
     * Creates and fills homologueChains matrix.
     * @param chains All chains of this complex graph
     */
    private void createHomologueChainsMatrix(List<Chain> chains) {
    if(! Settings.getBoolean("PTGLgraphComputation_B_silent")) {
            System.out.println("  Computing CG contacts.");
        }
    
        homologueChains = new Integer[chains.size()][chains.size()];
        if (chains.size() > 1) {
            for (Integer i = 0; i < chains.size(); i++) {
                
                for (Integer j = 0; j < chains.size(); j++) {

                    String compareChainID = chains.get(j).getPdbChainID();
                    // make sure no chainName is matched with itself 
                    if (chains.get(i).getHomologues() != null) {
                        if ((chains.get(i).getHomologues().contains(compareChainID)) && (!Objects.equals(i, j))) {
                            homologueChains[i][j] = 1;
                            // fill bottom-left and top-right triangle of matrix (required since chains are unordered b/c of speedup)
                            homologueChains[j][i] = 1;
                        } else {
                            homologueChains[i][j] = 0;
                        }
                    } else {
                        homologueChains[i][j] = 0;
                    }
                }
            }
        }
    }
    
    
    /**
     * Creates and fills number of chain interactions matrix.
     * @param resContacts Residue contacts
     */
    private void calculateNumChainInteractions(List<MolContactInfo> resContacts) {
        numChainInteractions = new Integer[numberChains][numberChains];
        numChainInteractionsCTL = new Integer[numberChains][numberChains];
        ArrayList<Integer> chainIDList = new ArrayList<>();
        Integer roundedCTL;
        
        for(Integer i = 0; i < resContacts.size(); i++){             
            ComplexGraph.Vertex chainA = getVertexFromChain(resContacts.get(i).getMolA().getChainID());
            ComplexGraph.Vertex chainB = getVertexFromChain(resContacts.get(i).getMolB().getChainID());
                      
            Integer chainAint = Integer.parseInt(chainA.toString());
            Integer chainBint = Integer.parseInt(chainB.toString());
            
            if(!chainIDList.contains(chainAint)){
                chainIDList.add(chainAint);
            }
            if(!chainIDList.contains(chainBint)){
                chainIDList.add(chainBint);
            }
            
            // We only want interchain contacts
            if (!chainA.equals(chainB)){
                if(! Settings.getBoolean("PTGLgraphComputation_CTL")){
                    if(numChainInteractions[chainAint][chainBint] == null){
                        numChainInteractions[chainAint][chainBint] = 1;
                        numChainInteractions[chainBint][chainAint] = 1;
                    } else {
                        numChainInteractions[chainAint][chainBint]++;
                        numChainInteractions[chainBint][chainAint]++;
                    }
                }
                else{
                    if(resContacts.get(i).getNumContactsCTL().equals(resContacts.get(i).getNumContactsTotal())){
                        if(numChainInteractionsCTL[chainAint][chainBint] == null){
                            numChainInteractionsCTL[chainAint][chainBint] = 1;
                            numChainInteractionsCTL[chainBint][chainAint] = 1;
                        } else {
                            numChainInteractionsCTL[chainAint][chainBint]++;
                            numChainInteractionsCTL[chainBint][chainAint]++;
                        }
                    }
                    else{
                        if(numChainInteractions[chainAint][chainBint] == null){
                            numChainInteractions[chainAint][chainBint] = 1;
                            numChainInteractions[chainBint][chainAint] = 1;
                        } else {
                            numChainInteractions[chainAint][chainBint]++;
                            numChainInteractions[chainBint][chainAint]++;
                        }
                    }
                }
            }
            
            if(i == resContacts.size()-1){
                if(Settings.getBoolean("PTGLgraphComputation_CTL")){
                    for(Integer k = 0; k < numberChains; k++){
                        for(Integer l = 0; l < numberChains; l++){
                            //weight and round CTL 
                            Integer mergeChainA = chainIDList.get(k);
                            Integer mergeChainB = chainIDList.get(l);
                            if (!mergeChainA.equals(mergeChainB)){
                                if(numChainInteractionsCTL[mergeChainA][mergeChainB] == null){
                                    roundedCTL = 0;
                                } else{
                                    roundedCTL = (int) Math.round(Settings.getFloat("PTGLgraphComputation_CTL_edge_weight")*numChainInteractionsCTL[mergeChainA][mergeChainB]);
                                }
                                
                                if(numChainInteractions[mergeChainA][mergeChainB] == null){
                                    numChainInteractions[mergeChainA][mergeChainB] = roundedCTL;
                                } else{
                                    numChainInteractions[mergeChainA][mergeChainB] = numChainInteractions[mergeChainA][mergeChainB] + roundedCTL;
                                }
                            }   
                        }
                    }
                }
            }
        }      
    }
    
    /**
     * Creates edges and if required contact info for csv.
     * @param resContacts Residue contacts
     */
    private void createEdges(List<MolContactInfo> resContacts) {
        if (createContactInfo) {
            contactInfo = new ArrayList<>();
            contactInfo.add("ChainA;ChainB;ResNameA;ResNameB;resTypeA;resTypeB;BB;BC;BL;CB;CL;CC;HB1;HB2;LB;LC;LL;"
                  + "BBDist;BCDist;BLDist;CBDist;CLDist;CCDist;HB1Dist;HB2Dist;LBDist;LCDist;LLDist");
        }
        
        // inform here if edge threshold is >1 (below it would result in multiple prints)
        if (Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold") > 1) {    
            if (! Settings.getBoolean("PTGLgraphComputation_B_silent")) {
                System.out.println("  Complex graph contact threshold for edges is set to "
                        + Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold").toString()
                        + ". Resulting graphs may differ from default setting '1' where all "
                        + "edges are drawn.");
            }
        }
        
        //get CTL Setting factor
        if(Settings.getBoolean("PTGLgraphComputation_CTL")){
            settingFactor = BigDecimal.valueOf(Settings.getFloat("PTGLgraphComputation_CTL_edge_weight"));
        }
        
        // create edges for all contacts
        for(Integer i = 0; i < resContacts.size(); i++) {
            ComplexGraph.Vertex chainA = getVertexFromChain(resContacts.get(i).getMolA().getChainID());
            ComplexGraph.Vertex chainB = getVertexFromChain(resContacts.get(i).getMolB().getChainID());
            
            Integer chainAint = Integer.parseInt(chainA.toString());
            Integer chainBint = Integer.parseInt(chainB.toString());
            
            // We only want interchain contacts with a certain threshold of contacts
            if (chainsHaveEnoughContacts(chainAint, chainBint)){
                
                MolContactInfo curResCon = resContacts.get(i);
                
                // Die Datenkrake
                String chainAString = curResCon.getMolA().getChainID().toString();
                String chainBString = curResCon.getMolB().getChainID().toString();
                String resNameA = curResCon.getName3A();
                String resNameB = curResCon.getName3B();
                String resTypeA = curResCon.getMolA().getSSETypePlcc();
                String resTypeB = curResCon.getMolB().getSSETypePlcc();
                
                String BBDist = curResCon.getBBContactDist().toString();
                String BCDist = curResCon.getBCContactDist().toString();
                String BLDist = curResCon.getBLContactDist().toString();
                String CBDist = curResCon.getCBContactDist().toString();
                String CLDist = curResCon.getCLContactDist().toString();
                String CCDist = curResCon.getCCContactDist().toString();
                String HB1Dist = curResCon.getHB1Dist().toString();
                String HB2Dist = curResCon.getHB2Dist().toString();
                String LBDist = curResCon.getLBContactDist().toString();
                String LCDist = curResCon.getLCContactDist().toString();
                String LLDist = curResCon.getLLContactDist().toString();
                                
                String numBB = curResCon.getNumContactsBB().toString();
                String numBC = curResCon.getNumContactsBC().toString();
                String numBL = curResCon.getNumContactsBL().toString();
                String numCB = curResCon.getNumContactsCB().toString();
                String numCL = curResCon.getNumContactsCL().toString();
                String numCC = curResCon.getNumContactsCC().toString();
                String numHB1 = curResCon.getNumContactsHB1().toString();
                String numHB2 = curResCon.getNumContactsHB2().toString();
                String numLB = curResCon.getNumContactsLB().toString();
                String numLC = curResCon.getNumContactsLC().toString();
                String numLL = curResCon.getNumContactsLL().toString();
                
                
                // Only if both residues belong to a SSE..
                if (curResCon.getMolA().getSSE() != null && curResCon.getMolB().getSSE() != null) {
                    //.. and are of type helix or strand
                    if ((!Objects.equals(curResCon.getMolA().getSSE().getSSETypeInt(), SSE.SSECLASS_NONE)
                            && !Objects.equals(curResCon.getMolA().getSSE().getSSETypeInt(), SSE.SSECLASS_OTHER))
                            && (!Objects.equals(curResCon.getMolB().getSSE().getSSETypeInt(), SSE.SSECLASS_NONE)
                            && !Objects.equals(curResCon.getMolB().getSSE().getSSETypeInt(), SSE.SSECLASS_OTHER))) {

                        Integer ResASseDsspNum = curResCon.getMolA().getSSE().getStartDsspNum();
                        Integer ResBSseDsspNum = curResCon.getMolB().getSSE().getStartDsspNum();

                        
                        Integer tmp;
                        String tmpName;
                        if (ResASseDsspNum > ResBSseDsspNum) {
                            tmp = ResBSseDsspNum;
                            tmpName = chainBString;

                            ResBSseDsspNum = ResASseDsspNum;
                            chainBString = chainAString;

                            ResASseDsspNum = tmp;
                            chainAString = tmpName;
                        }

                        List<Integer> SSEPair = Arrays.asList(ResASseDsspNum, ResBSseDsspNum);
                        List<String> NamePair = Arrays.asList(chainAString, chainBString);

                        numSSEContactChainNames.put(SSEPair, NamePair);
                        if (numSSEContacts.get(SSEPair) == null) {
                            numSSEContacts.put(SSEPair, 1);
                        } else {
                            numSSEContacts.put(SSEPair, numSSEContacts.get(SSEPair) + 1);
                        }
                    }

                }
                      
                                
                //This is for CSV output only. Like the most of the code above.
                if (createContactInfo) {
                contactInfo.add(chainAString + ";" + chainBString + ";" + resNameA + ";" + resNameB + ";" + 
                            resTypeA + ";" + resTypeB + ";" + numBB + ";" + numBC + ";" + numBL + ";" + 
                            numCB + ";" + numCL + ";" + numCC + ";" + numHB1 + ";" + numHB2 + ";" + 
                            numLB + ";" + numLC + ";" + numLL + ";" + BBDist + ";" + BCDist + ";" + 
                            BLDist + ";" + CBDist + ";" + CLDist + ";" + CCDist + ";" + HB1Dist + ";" + 
                            HB2Dist + ";" + LBDist + ";" + LCDist + ";" + LLDist);
                }

                String[] chainPair = {chainAString, chainBString};
                
                if (getEdge(chainA, chainB) == null){
                    // We don't have an edge yet, but need one, so create an edge
                    ComplexGraph.Edge e1 = createEdge(chainA, chainB);
                    chainNamesInEdge.put(e1, chainPair);
                    if(Settings.getBoolean("PTGLgraphComputation_CTL") && (resContacts.get(i).getNumContactsCTL().equals(resContacts.get(i).getNumContactsTotal()))){
                      mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).put(e1, settingFactor); // rather weird: 1 is added here, therefore the first contact is skipped later. We will have to change this and sum up the others in the end, (which of them depending on the graph type)
                    }
                    else{
                    mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).put(e1, BigDecimal.ONE); // rather weird: 1 is added here, therefore the first contact is skipped later. We will have to change this and sum up the others in the end, (which of them depending on the graph type)    
                    }
                    numHelixHelixInteractionsMap.put(e1, 0);
                    numHelixStrandInteractionsMap.put(e1, 0);
                    numHelixCoilInteractionsMap.put(e1, 0);
                    numHelixLigandInteractionsMap.put(e1, 0);
                    numStrandStrandInteractionsMap.put(e1, 0);
                    numStrandCoilInteractionsMap.put(e1, 0);
                    numStrandLigandInteractionsMap.put(e1, 0);
                    numCoilCoilInteractionsMap.put(e1, 0);
                    numCoilLigandInteractionsMap.put(e1, 0);
                    numLigandLigandInteractionsMap.put(e1, 0);
                    
                    
                    if (resContacts.get(i).getMolA().getSSE()!=null){
                        // the 1st residue of this contact belongs to a valid PTGL SSE
                        int firstSSEClass = resContacts.get(i).getMolA().getSSE().getSSETypeInt();
                        switch (firstSSEClass){                            
                            case 1: // SSECLASS_HELIX
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixHelixInteractionsMap.put(e1, 1);  // we are creating a new edge, so this is the first contact
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numHelixStrandInteractionsMap.put(e1, 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            //System.out.println("Ligand Contact");
                                            numHelixLigandInteractionsMap.put(e1, 1);
                                            break;
                                        case 4:
                                            numHelixCoilInteractionsMap.put(e1, 1);
                                            break;
                                    }
                                }
                                else{
                                    numHelixCoilInteractionsMap.put(e1, 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 2: // SSECLASS_BETASTRAND
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixStrandInteractionsMap.put(e1, 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandStrandInteractionsMap.put(e1, 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            numStrandLigandInteractionsMap.put(e1, 1);
                                            //System.out.println("Ligand Contact");
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numStrandCoilInteractionsMap.put(e1, 1);
                                            break;
                                    }
                                }
                                else{
                                    numStrandCoilInteractionsMap.put(e1, 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 3: // SSECLASS_LIGAND
                                //System.out.println("Ligand Contact");
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixLigandInteractionsMap.put(e1, 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandLigandInteractionsMap.put(e1, 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            numLigandLigandInteractionsMap.put(e1, 1);
                                            //System.out.println("Ligand Contact");
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numCoilLigandInteractionsMap.put(e1, 1);
                                            break;
                                    }
                                }
                                else{
                                    numCoilLigandInteractionsMap.put(e1, 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 4:
                                //System.out.println("Other Contact");
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixCoilInteractionsMap.put(e1, 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandCoilInteractionsMap.put(e1, 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            numCoilLigandInteractionsMap.put(e1, 1);
                                            //System.out.println("Ligand Contact");
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numCoilCoilInteractionsMap.put(e1, 1);
                                            break;
                                    }
                                }
                                else{
                                    numCoilCoilInteractionsMap.put(e1, 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                        }
                    }
                    else{
                        // the first residue of this contact does NOT belong to a valid PTGL SSE, i.e., it is a coil
                        if (resContacts.get(i).getMolB().getSSE()!=null){
                            int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                            switch (secondSSE){
                                case 1: // SSECLASS_HELIX
                                    numHelixCoilInteractionsMap.put(e1, 1);
                                    break;
                                case 2: // SSECLASS_BETASTRAND
                                    numStrandCoilInteractionsMap.put(e1, 1);
                                    break;
                                case 3: // SSECLASS_LIGAND
                                    //System.out.println("Ligand Contact");
                                    numCoilLigandInteractionsMap.put(e1, 1);
                                    break;
                                case 4:
                                    //System.out.println("Other Contact");
                                    numCoilCoilInteractionsMap.put(e1, 1);
                                    break;
                            }

                        }
                        else{
                            numCoilCoilInteractionsMap.put(e1, 1);
                            //System.out.println("Loop-loop Contact");
                        }
                    }
                    //System.out.println("Contact found between chainName " + resContacts.get(i).getMolA().getChainID() + " and chainName " + resContacts.get(i).getMolB().getChainID());
                }
                else{
                    // We already have an edge, just adjust values
                    if(Settings.getBoolean("PTGLgraphComputation_CTL") && (resContacts.get(i).getNumContactsCTL().equals(resContacts.get(i).getNumContactsTotal()))){
                      mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).put(getEdge(chainA, chainB), mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(getEdge(chainA, chainB)).add(settingFactor)); 
                    }
                    else{
                    mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).put(getEdge(chainA, chainB), mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(getEdge(chainA, chainB)).add(BigDecimal.ONE)); 
                    }

                    if (resContacts.get(i).getMolA().getSSE()!=null){
                        // first residue of contact belongs to valid PTGL SSE, i.e., is NOT a coil
                        int firstSSE = resContacts.get(i).getMolA().getSSE().getSSETypeInt();
                        switch (firstSSE){
                            case 1: // SSECLASS_HELIX
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixHelixInteractionsMap.put(getEdge(chainA, chainB), numHelixHelixInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numHelixStrandInteractionsMap.put(getEdge(chainA, chainB), numHelixStrandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            //System.out.println("Ligand Contact");
                                            numHelixLigandInteractionsMap.put(getEdge(chainA, chainB), numHelixLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numHelixCoilInteractionsMap.put(getEdge(chainA, chainB), numHelixCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                    }
                                }
                                else{
                                    numHelixCoilInteractionsMap.put(getEdge(chainA, chainB), numHelixCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 2: // SSECLASS_BETASTRAND
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixStrandInteractionsMap.put(getEdge(chainA, chainB), numHelixStrandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandStrandInteractionsMap.put(getEdge(chainA, chainB), numStrandStrandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            //System.out.println("Ligand Contact");
                                            numStrandLigandInteractionsMap.put(getEdge(chainA, chainB), numStrandLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numStrandCoilInteractionsMap.put(getEdge(chainA, chainB), numStrandCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                    }
                                }
                                else{
                                    numStrandCoilInteractionsMap.put(getEdge(chainA, chainB), numStrandCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 3: // SSECLASS_LIGAND
                                //System.out.println("Ligand Contact");
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixLigandInteractionsMap.put(getEdge(chainA, chainB), numHelixLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandLigandInteractionsMap.put(getEdge(chainA, chainB), numStrandLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            //System.out.println("Ligand Contact");
                                            numLigandLigandInteractionsMap.put(getEdge(chainA, chainB), numLigandLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numCoilLigandInteractionsMap.put(getEdge(chainA, chainB), numCoilLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                    }
                                }
                                else{
                                    numCoilLigandInteractionsMap.put(getEdge(chainA, chainB), numCoilLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                            case 4:
                                //System.out.println("Other Contact");
                                if (resContacts.get(i).getMolB().getSSE()!=null){
                                    int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                                    switch (secondSSE){
                                        case 1: // SSECLASS_HELIX
                                            numHelixCoilInteractionsMap.put(getEdge(chainA, chainB), numHelixCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 2: // SSECLASS_BETASTRAND
                                            numStrandCoilInteractionsMap.put(getEdge(chainA, chainB), numStrandCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 3: // SSECLASS_LIGAND
                                            //System.out.println("Ligand Contact");
                                            numCoilLigandInteractionsMap.put(getEdge(chainA, chainB), numCoilLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                        case 4:
                                            //System.out.println("Other Contact");
                                            numCoilCoilInteractionsMap.put(getEdge(chainA, chainB), numCoilCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                            break;
                                    }
                                }
                                else{
                                    numCoilCoilInteractionsMap.put(getEdge(chainA, chainB), numCoilCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    //System.out.println("Loop Contact");
                                }
                                break;
                        }
                    }
                    else{   // first residue of contact does NOT belong to valid PTGL SSE, i.e., is a coil
                        if (resContacts.get(i).getMolB().getSSE()!=null){
                            int secondSSE = resContacts.get(i).getMolB().getSSE().getSSETypeInt();
                            switch (secondSSE){
                                case 1: // SSECLASS_HELIX
                                    numHelixCoilInteractionsMap.put(getEdge(chainA, chainB), numHelixCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    break;
                                case 2: // SSECLASS_BETASTRAND
                                    numStrandCoilInteractionsMap.put(getEdge(chainA, chainB), numStrandCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    break;
                                case 3: // SSECLASS_LIGAND
                                    //System.out.println("Ligand Contact");
                                    numCoilLigandInteractionsMap.put(getEdge(chainA, chainB), numCoilLigandInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    break;
                                case 4:
                                    //System.out.println("Other Contact");
                                    numCoilCoilInteractionsMap.put(getEdge(chainA, chainB), numCoilCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                                    break;
                            }
                        }
                        else{
                            numCoilCoilInteractionsMap.put(getEdge(chainA, chainB), numCoilCoilInteractionsMap.get(getEdge(chainA, chainB)) + 1);
                            //System.out.println("Loop Contact");
                        }
                    }
                }
                
                // TODO: Test by Tim: maybe we should delete the edge if it has no contacts:
                
                if(mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(getEdge(chainA, chainB)) == BigDecimal.ZERO) {
                    removeEdge(getEdge(chainA, chainB));
                }
            }
        } // end of loop over all res contacts
        computeNormalizedEdgeWeights();  // do this here instead of in loop, so we need to compute it only once
    }
    
    
    private void computeNormalizedEdgeWeights() {
        BigDecimal curMinimumMultNormEdgeWeight = BigDecimal.ONE;
        
        for (Edge e : mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).keySet()) {
            
            BigDecimal tmpContacts = mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(e);
            int tmpNumRes1 = mapChainIdToLength.get(chainNamesInEdge.get(e)[0]);
            int tmpNumRes2 = mapChainIdToLength.get(chainNamesInEdge.get(e)[1]);
	    double tmpGyradius1 = mapChainIdToGyradius.get(chainNamesInEdge.get(e)[0]);
            double tmpGyradius2 = mapChainIdToGyradius.get(chainNamesInEdge.get(e)[1]);
            
            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 3) {
                System.out.println("[DEBUG LV 4] AdditiveLengthNormalization: tmpContacts | tmpNumRes1 | tmpNumRes2");
                System.out.println("  [DEBUG LV 4] " + tmpContacts);
                System.out.println("  [DEBUG LV 4] " + tmpNumRes1);
                System.out.println("  [DEBUG LV 4] " + tmpNumRes2);
            }
            
            BigDecimal tmpAddNormWeight = ComplexGraphEdgeWeightTypes.computeAdditiveLengthNormlization(tmpContacts, tmpNumRes1, tmpNumRes2);
            BigDecimal tmpMultNormWeight = ComplexGraphEdgeWeightTypes.computeMultiplicativeLengthNormlization(tmpContacts, tmpNumRes1, tmpNumRes2);
	    BigDecimal tmpMultIndNormWeight = ComplexGraphEdgeWeightTypes.computeIndividualMultiplicativeLengthNormlization(tmpContacts, tmpNumRes1, tmpNumRes2);
            BigDecimal tmpAddIndNormWeight = ComplexGraphEdgeWeightTypes.computeIndividualAdditiveLengthNormlization(tmpContacts, tmpNumRes1, tmpNumRes2);
            BigDecimal tmpAddGyrNormWeight = ComplexGraphEdgeWeightTypes.computeAdditiveGyrationNormlization(tmpContacts, tmpGyradius1, tmpGyradius2);
            BigDecimal tmpAddGyrIndNormWeight = ComplexGraphEdgeWeightTypes.computeIndividualAdditiveGyrationNormlization(tmpContacts, tmpGyradius1, tmpGyradius2);
            
            curMinimumMultNormEdgeWeight = curMinimumMultNormEdgeWeight.min(tmpMultNormWeight);  // update min if necessary
            
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ADDITIVE_LENGTH_NORMALIZATION).put(e, tmpAddNormWeight);
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.MULTIPLICATIVE_LENGTH_NORMALIZATION).put(e, tmpMultNormWeight);
	    mapWeightNamesToMapEdgeValues.get(EdgeWeightType.MULTIPLICATIVE_INDIVIDUAL_NORMALIZATION).put(e, tmpMultIndNormWeight);
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ADDITIVE_INDIVIDUAL_NORMALIZATION).put(e, tmpAddIndNormWeight);
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ADDITIVE_GYRATION_NORMALIZATION).put(e, tmpAddGyrNormWeight);
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ADDITIVE_GYRATION_INDIVIDUAL_NORMALIZATION).put(e, tmpAddGyrIndNormWeight);
        }
        
        minimumMultiplicativeLengthNormalizedEdgeWeight = curMinimumMultNormEdgeWeight;
        
        // now e.x. minimum multiplicative normalized edge weight -> can compute lucid edge weight
        computeLucidLengthNormalizedEdgeWeight();
    }
    
    
    private void computeLucidLengthNormalizedEdgeWeight() {
        mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).keySet().forEach(e -> {
            BigDecimal curNormWeight = mapWeightNamesToMapEdgeValues.get(EdgeWeightType.MULTIPLICATIVE_LENGTH_NORMALIZATION).get(e);
            BigDecimal lucidWeight = ComplexGraphEdgeWeightTypes.computeLucidMultiplicativeLengthNormlization(curNormWeight, minimumMultiplicativeLengthNormalizedEdgeWeight);
            mapWeightNamesToMapEdgeValues.get(EdgeWeightType.LUCID_MULTIPLICATIVE_NORMALIZATION).put(e, lucidWeight);
        });     
    }
    
    
    private void assemblyPrediction() {
        Boolean silent = Settings.getBoolean("PTGLgraphComputation_B_silent");
        if (!silent) {
            System.out.println("  Computing assembly prediction based on different edge weight types...");
        }
        String[] desiredWeightTypes = Settings.getList("PTGLgraphComputation_L_AP_weight_types");
        for (EdgeWeightType weightType : EdgeWeightType.values()) {
            if (! Arrays.stream(desiredWeightTypes).anyMatch(weightType.shortName::equals)) {
                continue;
            }
            System.out.println("    " + weightType.name);
            
            AgglomerativeClustering clustering;
            if(weightType.shortName.equals("addGyrInd") || weightType.shortName.equals("addGyr")){
                // Maybe a different constructor is needed to deal with the gyradius instead of the length, but the best would be to leave everything behind the same, except for some Double/Integer things
                // provide the chain instances, which can then be accessed by allChains.get(i)
                clustering = 
                        new AgglomerativeClustering(getEdgesAsArray(),  vertexMapToVertexIdMap(chainGyradiusMap), weightType, vertexMapToVertexIdMap(labelNodeMap), allChains);
            }
            else{
                clustering = 
                        new AgglomerativeClustering(getEdgesAsArray(), vertexMapToVertexIdMap(chainLengthMap), weightType, vertexMapToVertexIdMap(labelNodeMap), allChains);
            }
            
            ClusteringResult clusteringResult;
            if(Settings.getInteger("PTGLgraphComputation_I_type_assembly_prediction") == 2){
                int iterations = Settings.getInteger("PTGLgraphComputation_I_stochastic_iterations");
                ConsensusTree consensus = new ConsensusTree(clustering, iterations, clustering.getlabelMap());
                clusteringResult = consensus.computeConsensus();
                consensus.writeStatistics(this.pdbid+"_"+weightType.toString());
            }
            else{
                clusteringResult = clustering.chainInformationClustering();
            }
            
            System.out.println("      " + clusteringResult.toNewickString(vertexMapToVertexIdMap(labelNodeMap)));
            
            //System.out.println("      Large Interface Score: " + clusteringResult.getScore().toString());
            for (String line : clusteringResult.getScoreLines()) {
                System.out.println("      Large Interface Score for " + line);
            }
        }
    }
    

    public Vertex getVertexFromChain(String chainID) {
        Iterator<Vertex> vertIter = this.getVertices().iterator();
        Vertex nextVert;
        while (vertIter.hasNext()) {
            nextVert = vertIter.next();
            if (ComplexGraph.this.proteinNodeMap.get(nextVert).equals(chainID)) {
                return nextVert;
            }
        }
        return null;
    }

    public String getPDBID() {
        return this.pdbid;
    }

    public Boolean chainsHaveEnoughContacts(Integer A, Integer B) {
        if (this.numChainInteractions[A][B] != null) {
            if (this.numChainInteractions[A][B] >= Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold"))
                return true;
            else {
                if (this.numChainInteractions[A][B] > 0) {
                    // just to be sure we do not count anything wrong
                    neglectedEdges++;  // this is WRONG b/c chainsHaveEnoughContacts is called for every mol contact
                    return false;
                }
            }
        } else {
            return false;
        }
        return false;
    }


    public void writeSSEComplexContactInfoToDB(String pdbid) {

    //DBManager.writeSSEComplexContactToDB()
        String chainA;
        String chainB;

        Integer numHelixHelixInteractions;
        Integer numHelixStrandInteractions;
        Integer numHelixLoopInteractions;
        Integer numSSInteractions;
        Integer numStrandLoopInteractions;
        Integer numLoopLoopInteractions;
        Integer numAllInteractions;
        Integer numDisulfides;

        int countInsert = 0;
        int countFail = 0;
        Boolean retVal = false;
        for (Map.Entry<List<Integer>, Integer> pair : this.numSSEContacts.entrySet()) {

            List<Integer> curSSEs = (List<Integer>) pair.getKey();
            chainA = this.numSSEContactChainNames.get(curSSEs).get(0);
            chainB = this.numSSEContactChainNames.get(curSSEs).get(1);
            
            Integer sse1_dssp_start = (Integer)curSSEs.get(0);
            Integer sse2_dssp_start = (Integer)curSSEs.get(1);
            Integer contactCount = (Integer)pair.getValue();
            
          
            Boolean res = false;
            
            
            // this action could result in an error due to the definition of a PTGL SSE
            // e.g. the SSE is too short and is merged to another SSE or not defined in the DB
            try {
                res = DBManager.writeSSEComplexContactToDB(pdbid, chainA, chainB, sse1_dssp_start, sse2_dssp_start, contactCount);
                //it.remove(); // avoids a ConcurrentModificationException
                if(res) {
                    countInsert++;
                    retVal = true;
                }
                else {
                    countFail++;
                    retVal = false;
                }
            } catch (SQLException ex) {
                DP.getInstance().e("ComplexGraph", "writeSSEComplexContactInfoToDB: SQL exception: '" + ex.getMessage() + "'.");
                Logger.getLogger(ComplexGraph.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        if(! Settings.getBoolean("PTGLgraphComputation_B_silent")) {
            System.out.println("    SSE Contacts written to DB: " + countInsert + " inserted, " + countFail + " skipped (contact involved coils).");
        }
    }

    /**
     * Separates interaction maps and associates them with a chain name then
     * calls DBManager.writeComplexContactToDB()
     *
     * @return true if DBManager is succesfull
     */
    public boolean writeChainComplexContactInfoToDB() {

        String chainA;
        String chainB;

        Integer numHelixHelixInteractions;
        Integer numHelixStrandInteractions;
        Integer numHelixLoopInteractions;
        Integer numStrandStrandInteractions;
        Integer numStrandLoopInteractions;
        Integer numLoopLoopInteractions;
        Integer numResResContacts;
        Integer numDisulfides;

        for (Map.Entry pair : this.chainNamesInEdge.entrySet()) {

            Edge curEdge = (Edge) pair.getKey();
            String[] chainPair = (String[]) pair.getValue();

            chainA = chainPair[0];
            chainB = chainPair[1];

            // interactions with ligands are NOT yet written to the database!
            numHelixHelixInteractions = this.numHelixHelixInteractionsMap.get(curEdge);
            numHelixStrandInteractions = this.numHelixStrandInteractionsMap.get(curEdge);
            numHelixLoopInteractions = this.numHelixCoilInteractionsMap.get(curEdge);
            numStrandStrandInteractions = this.numStrandStrandInteractionsMap.get(curEdge);
            numStrandLoopInteractions = this.numStrandCoilInteractionsMap.get(curEdge);
            numLoopLoopInteractions = this.numCoilCoilInteractionsMap.get(curEdge);
            numDisulfides = this.numDisulfidesMap.get(curEdge);
            numResResContacts = mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(curEdge).toBigInteger().intValue();

            Integer[] interactionNums = {numHelixHelixInteractions, numHelixStrandInteractions, numHelixLoopInteractions,
                numStrandStrandInteractions, numStrandLoopInteractions, numLoopLoopInteractions,
                numDisulfides, numResResContacts};

            // make sure no entry is null or something shitty
            for (int i = 0; i < interactionNums.length; i++) {
                if (interactionNums[i] == null) {
                    interactionNums[i] = 0;
                }
            }

            try {
                DBManager.writeChainComplexContactToDB(pdbid, chainA, chainB, interactionNums);
                //it.remove(); // avoids a ConcurrentModificationException
            } catch (SQLException ex) {
                Logger.getLogger(ComplexGraph.class.getName()).log(Level.SEVERE, null, ex);
                return false;
            }
        }
        return true;
    }

    // Old color function: was insensitive to the case of many vertices and few MolIDs (Homologues)
    /*
    private float getUniqueColor(Integer numVertices) {
        float step = 360 / (numVertices + 1); // +1 to avoid double red
        float hue = (lastColorStep + step) / 360;
        lastColorStep += step;      
        return hue;
    }
    */
    
    
    /**
     * Returns a hue based on the molID in relation to the total number of molIDs. Requires molIDs to be filled.
     * All homologous chains get the same color and the colors of all homologues are distinct with equal distance to each other. 
     * @param molID
     * @param numberHomologues total number of colors needed
     * @return 
     */
    public float getUniqueHue(Integer molID, Integer numberHomologues) {
        //return (360f / (molIDs.size()) * molID) / 360f;
        //return 1f / molIDs.size() * (molID - 1);
        return 1f / numberHomologues * (molID - 1);
    }
    
    /**
     * Returns an ArrayList containing [hue,saturation,brightness] for the HSB-color-coding.
     * Iterates over hue and adds iteration over saturation and brightness after a certain number of distinct vertices.
     * since only homologous chains are colored the number of colors is equal to the number of homologous chains
     * @param molIdName molID of the current vertex
     * @param colors contains current HSB-color-coding
     * @param Homologues contains only molIDs that have homologous chains
     * @return
     */
    public ArrayList<Float> getColorCode(String molIdName, ArrayList<Float> colors, ArrayList<String> Homologues){
        
        Integer molID = Homologues.indexOf(molIdName) + 1; //since only homologous chains are colored the molID needs to be set according to the number of homologues
        Integer numberHomologues = Homologues.size();
        
        //Divides the vertices into groups of size x. x can't be larger than 10 because we prefer more groups to larger groups.
        Integer pieces = 10; 
        Integer x = (int)(numberHomologues / pieces);
        if(x > pieces){
            pieces = x;
            x = 10;
        }
        x = x-1;
        
        Integer maxSat = 2; //Maximum of different saturation values.
        
        
        if(x <= 0){ //For graphs with up to 19 vertices only the hue is iterated.
            colors.set(0, getUniqueHue(molID, numberHomologues));
            colors.set(1, 1f);
            colors.set(2, 0.8f);
            
        }
        else if(x <= 2){ // For graphs with up to 39 vertices hue and brightness are iterated.
            if((molID-1) % (x+1) == 0){ //the hue with default saturation and brightness is used first.
                
                colors.set(0, getUniqueHue(molID, numberHomologues));
                colors.set(1, 1f);
                colors.set(2, 0.8f); 
            }
            else{ //The brightness is iterated.
                
                colors.set(0, getUniqueHue(molID, numberHomologues));
                colors.set(1, 1f);
                colors.set(2, 0.6f + (0.4f / x * ((molID-1) % x))); //All values under 0.6 produce dark colors that are harder to distinguish
            }
        }
        else{ // Brightness and saturation are iterated
            if((molID-1) % (x+1) == 0){
               
               colors.set(0, getUniqueHue(molID, numberHomologues));
               colors.set(1, 1f);
               colors.set(2, 0.8f);
               
            }
            else{
                
                colors.set(0, getUniqueHue(molID, numberHomologues));
                Integer y = (int)Math.ceil(x/maxSat) +1; // number of changes the brightness has for one group of x
                Integer counter = ((molID-1)%(x+1)); // calculates the current position of molID in its group of x vertices
                
                //the brightness and saturation are calculated according to the counter
                if((counter == 1)|| (counter % (maxSat+1)) == 1){
                    colors.set(2, 0.6f + (0.4f / y * (counter % y)));
                }
                
                colors.set(1, 0.6f + (0.4f / maxSat * ((counter-1) % (maxSat+1))));
                
            }
        }       
        return colors;
    }
    

    /**
     * Draws a complex graph
     * @param nonProteinGraph whether the graph is a non-protein graph and should be drawn black and white
     * @param cg the complex graph
     * @param molInfoForChains info mapping chain IDs (like "A") to their macromolecule (MOL_ID in PDB file, e.g., "1"). Give an empty one if you dont know
     * @return a draw result
     */
    private static DrawResult drawChainLevelComplexGraphG2D(Boolean nonProteinGraph, ComplexGraph cg, Map<String, String> molInfoForChains) {
        
        Boolean hasRna = false;
        Integer numVerts = cg.getVertices().size();

        Boolean bw = nonProteinGraph;
        
        ArrayList<String> molNames = new ArrayList<String>(cg.molMap.values());

    // All these values are in pixels
        // page setup
        PageLayout pl = new PageLayout(numVerts, molNames);
        Position2D vertStart = pl.getVertStart();

    // ------------------------- Prepare stuff -------------------------
        // TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed into integer pixels
        //BufferedImage bi = new BufferedImage(pl.getPageWidth(), pl.getPageHeight(), BufferedImage.TYPE_INT_ARGB);
        SVGGraphics2D ig2;

        //if(Settings.get("PTGLgraphComputation_S_img_output_format").equals("SVG")) {                    
        // Apache Batik SVG library, using W3C DOM tree implementation
        // Get a DOMImplementation.
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        // Create an instance of the SVG Generator.
        ig2 = new SVGGraphics2D(document);
        //ig2.getRoot(document.getDocumentElement());
        // }
        //else {
        //    ig2 = (SVGGraphics2D)bi.createGraphics();
        //}

        ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // make background white
        ig2.setPaint(Color.WHITE);
        ig2.fillRect(0, 0, pl.getPageWidth(), pl.getPageHeight());
        ig2.setPaint(Color.BLACK);

    //pl.drawAreaOutlines(ig2);
        // prepare font
        Font font = pl.getStandardFont();
        ig2.setFont(font);
        FontMetrics fontMetrics = ig2.getFontMetrics();

    // ------------------------- Draw header -------------------------
        // check width of header string
        String proteinHeader = "Chain Complex Graph of PDB entry " + cg.pdbid + " [V=" + cg.getVertices().size() + ", E=" + cg.getEdges().size() + "]";
        proteinHeader += "\n" + Settings.getImagesCreatorVersionLine();
        if (Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold") > 1) {
            proteinHeader += "\n(Interchain contact threshold is set to " + Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold") + ")"; // could include neglected edges if they were correctly computed
        }
        
        //Integer stringWidth = fontMetrics.stringWidth(proteinHeader);       // Should be around 300px for the text above
        Integer stringHeight = fontMetrics.getAscent();
        // above line should be:
             // Integer stringHeight = fontMetrics.getAscent() * proteinHeader.split("\n").length;
        //   but vertices and edges do not react to it and so labels would not fit. PageLayout.getVertStart would need to respect stringHeight of header
        
        String chainName;    // the SSE number in the primary structure, N to C terminus
        String chainNumber;  // the SSE number in this graph, 1..(this.size)

        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_header")) {
            //ig2.drawString(proteinHeader, pl.headerStart.x, pl.headerStart.y);
            DrawTools.drawStringLineBreaks(ig2, proteinHeader, pl.headerStart.x, pl.headerStart.y);
            //i2.drawString(addInfo, pl.headerStart.x, pl.headerStart.y + pl.textLineHeight);((1f/x*(molID-2)) - (int)(1f/x*(molID-2)))
        }

    // ------------------------- Draw the graph -------------------------
        // Draw the edges as arcs
        java.awt.Shape shape;
        Arc2D.Double arc;
        ig2.setStroke(new BasicStroke(2));  // thin edges
        Integer leftVert, rightVert, leftVertPosX, rightVertPosX, arcWidth, arcHeight, arcTopLeftX, arcTopLeftY, spacerX, spacerY;
        Integer labelPosX, labelPosY;

        String edges = cg.getEdges().toString();        
        for (Integer i = 0; i < cg.getVertices().size(); i++) {
            for (Integer j = i + 1; j < cg.getVertices().size(); j++) {

                String tmpEdgeString = "(" + i + ", " + j + ")";

                // If there is a contact...
                if (edges.contains(tmpEdgeString)) {

                // determine edge type and the resulting color
                    //edgeType = cg.getContactType(i, j);
                    
                    ig2.setPaint(new Color(0.3f, 0.3f, 0.3f, 0.3f)); //sets Color to Gray and Transparency. 0.0f = fully transparent, 1.0f = not transparent
                    
                    
                    if (bw) {
                        ig2.setPaint(Color.LIGHT_GRAY);
                    }      // for non-protein graphs

                // ----- complex graph specific stuff -----
                    // determine chain of SSEs
                    /*
                     for(Integer x = 0; x < cg.chainEnd.size(); x++){
                     if(i < cg.chainEnd.get(x)) {iChainID = x; break;}
                     }
                     for(Integer x = 0; x < cg.chainEnd.size(); x++){
                     if(j < cg.chainEnd.get(x)) {jChainID = x; break;}
                     }
                     if (!Objects.equals(iChainID, jChainID)) {ig2.setPaint(Color.PINK);}
                
                     */
                // ----- end complex graph specific stuff -----

                    // determine the center of the arc and the width of its rectangle bounding box
                    if (i < j) {
                        leftVert = i;
                        rightVert = j;
                    } else {
                        leftVert = j;
                        rightVert = i;
                    }

                    leftVertPosX = pl.getVertStart().x + (leftVert * pl.vertDist);
                    rightVertPosX = pl.getVertStart().x + (rightVert * pl.vertDist);

                    arcWidth = rightVertPosX - leftVertPosX;
                    arcHeight = arcWidth / 2;

                    arcTopLeftX = leftVertPosX;
                    arcTopLeftY = pl.getVertStart().y - arcHeight / 2;

                    spacerX = pl.vertRadius;
                    spacerY = 0;
                    
                    
                    // draw it                                                
                    arc = new Arc2D.Double(arcTopLeftX + spacerX, arcTopLeftY + spacerY, arcWidth, arcHeight, 0, 180, Arc2D.OPEN);
                    shape = ig2.getStroke().createStrokedShape(arc);
                    ig2.fill(shape);
                    ig2.setPaint(Color.GRAY);

                }
            }
        }
        
        for (Edge e : cg.getEdges()) {

            // recreate environment for loop over (i,j)
            int i = Integer.parseInt(e.toString().split(",")[0].replace("(", "").strip());  // num of first chain
            int j = Integer.parseInt(e.toString().split(",")[1].replace(")", "").strip());  // num of second chain
            
            if (i < j) {
                leftVert = i;
                rightVert = j;
            } else {
                leftVert = j;
                rightVert = i;
            }

            // TODO: is it clever to calculate everything again?
            leftVertPosX = pl.getVertStart().x + (leftVert * pl.vertDist);
            rightVertPosX = pl.getVertStart().x + (rightVert * pl.vertDist);

            arcWidth = rightVertPosX - leftVertPosX;
            arcHeight = arcWidth / 2;

            arcTopLeftX = leftVertPosX;
            arcTopLeftY = pl.getVertStart().y - arcHeight / 2;

            spacerX = pl.vertRadius;
            spacerY = 0;
            //calculate label positions
            labelPosX = leftVertPosX + arcWidth / 2 + 2;
            labelPosY = arcTopLeftY + spacerY - 5;
            //draw labels on arcs
            Font labelfont = new Font(Settings.get("PTGLgraphComputation_S_img_default_font"), Font.PLAIN, Settings.getInteger("PTGLgraphComputation_I_img_default_font_size") - 5);
            ig2.setFont(labelfont);
            ig2.setPaint(Color.BLACK);
            
            String cInteractionsString;
            // cInteractionsString = cg.lucidMultiplicativeLengthNormalizedEdgeWeight.get(e).setScale(0).toString();  // label as rounded factor
            cInteractionsString = cg.numChainInteractions[i][j].toString();

            if(cInteractionsString != null) {
                ig2.drawString(cInteractionsString, labelPosX, labelPosY + (stringHeight / 4));
            } else {
                // This may happen if not the top-right part of the numChainInteractions is filled (due to some different ordering in the chains)
                DP.getInstance().w("ComplexGraph", "Tried to read out a null entry from numChainInteractions while writing the edge weights in the chain-level CG. "
                        + "This is probably a programming error, so please inform the developers.");
            }
        }
        
        
        ig2.setFont(font);
        Rectangle2D.Double rect;
        ig2.setStroke(new BasicStroke(2));
        
        Iterator<Vertex> vertIterator = cg.getVertices().iterator();  // this and next line used to iterate vertices to get the number of Homologues
        Vertex curVertice;
        
        ArrayList<Float> colorCode = new ArrayList<Float>(); //ArrayList that contains the color codings for HSB-Colors
        colorCode.add(0.5f);
        colorCode.add(1f);
        colorCode.add(0.8f);
        
        ArrayList<String> Homologues = new ArrayList<String>(); //contains molIDs of homologue vertices
        Integer homologues;
        
        //Array Homologues is filled with homologue molIDs
        for(Integer i = 0; i < cg.getVertices().size(); i++){
            curVertice = vertIterator.next();
            String molID = molInfoForChains.get(cg.proteinNodeMap.get(curVertice));
            if(!Homologues.contains(molID) && molID != ""){
                homologues = Collections.frequency((molInfoForChains.values()), molID);
                if(homologues > 1){
                    Homologues.add(molID);
                }
            }
        }
        
        
        Iterator<Vertex> vertIter = cg.getVertices().iterator();  // this and next line used to iterate vertices to get current vertex number for colorCode
        Vertex curVert;    
        for (Integer i = 0; i < cg.getVertices().size(); i++) {
            curVert = vertIter.next();

            // standard color: saturation and brightness will not be changed
            // set hue, saturation, brighness
            float h = (float) 0.5;
            float s = (float) 1.0; // change this for saturation (higher = more saturated)
            float b = (float) 0.8; // change this for brightness (0.0 -> Dark/Black)
            
            
            String molID = molInfoForChains.get(cg.proteinNodeMap.get(curVert));
                        
            if (! bw) {
                //check if molID has homologues
                //only chains with homologues are colored
                if(Homologues.contains(molID)){
                    colorCode = cg.getColorCode(molID, colorCode, Homologues);
                    ig2.setPaint(Color.getHSBColor(colorCode.get(0), colorCode.get(1), colorCode.get(2)));
                }
                else{
                    ig2.setPaint(Color.GRAY);
                }
            } else {
                ig2.setPaint(Color.GRAY);
            }
            
            
            // following code produces a color palette for m graphs containing n = 1..m vertices
            //   it is best to set min_width and min_height to 1200 and 900
            //   and turn header and footer off
            //   you can switch between 'super old' getUniqueColor and 'old' getUniqueHue and 'new' getColorCode
            //   the result is saved as chain CG
            //   Important to note: the 'old' function is insensitive to the case of many nodes and few molIDs
            
            /*cg.molIDs.clear();
            for (Integer m = 1; m <= 200; m++) {
                System.out.println("m: " + m);
                //cg.lastColorStep = 0;
                cg.molIDs.add(m.toString());
                for (int n = 1; n <= m; n++) {
                    s = 1.0f;
                    rect = new Rectangle2D.Double(0 + n * 10, 0 + m * 10, 10, 10);
                    
                    //ig2.setPaint(Color.getHSBColor(cg.getUniqueHue(n), s, b));
                    
                    //ig2.setPaint(Color.getHSBColor(cg.getUniqueColor(m), s, b));
            
                    //new coding:
                    colorCode = cg.getColorCode(n, colorCode);
                    ig2.setPaint(Color.getHSBColor(colorCode.get(0), colorCode.get(1), colorCode.get(2)));
            
                    ig2.fill(rect);
                }
            }
            i = 1000;  // no edges and vertices
            */
            
            
        // pick color depending on SSE type

            // draw a shape based on SSE type
            
            Chain curChain = getChainByMacromolId(molID);
            if (curChain.isRnaChain()) {
                hasRna = true;
                // create triangle shape for RNA nodes
                int [] xpoints = {vertStart.x + (i * pl.vertDist) + 10, vertStart.x + (i * pl.vertDist) - 4, vertStart.x + (i * pl.vertDist) + 24};
                int [] ypoints = {vertStart.y - pl.getVertDiameter() + 12, vertStart.y - pl.getVertDiameter() + 36, vertStart.y - pl.getVertDiameter() + 36};
                int nPoints = 3;
                ig2.fillPolygon(xpoints, ypoints, nPoints);
            } else {
                // create rectangular shape for residue nodes
                rect = new Rectangle2D.Double(vertStart.x + (i * pl.vertDist) + pl.getVertDiameter() / 2, vertStart.y - pl.getVertDiameter() / 2, pl.getVertDiameter(), pl.getVertDiameter());
                AffineTransform rot_45deg = new AffineTransform();
                rot_45deg.rotate(0.785, vertStart.x + (i * pl.vertDist) + pl.getVertDiameter() / 2, vertStart.y - pl.getVertDiameter() / 2); // rotation around center of vertex
                ig2.fill(rot_45deg.createTransformedShape(rect));
            }
        }
        
        
        // Draw the markers for the N-terminus and C-terminus if there are any vertices in this graph            
        ig2.setStroke(new BasicStroke(2));
        ig2.setPaint(Color.BLACK);

        /*
         if( ! bw) {
         if(cg.getVertices().size() > 0) {                    
         ig2.drawString("N", vertStart.x - pl.vertDist, vertStart.y + 20);    // N terminus label
         ig2.drawString("C", vertStart.x + cg.getVertices().size() * pl.vertDist, vertStart.y + 20);  // C terminus label
         }
         }
         */
    // ************************************* footer **************************************
        if (Settings.getBoolean("PTGLgraphComputation_B_graphimg_footer")) {

        // Draw the vertex numbering into the footer
        
            // old version where in CGs > 99 vertices only each second footer information is printed
            /*
            // Determine the dist between vertices that will have their vertex number printed below them in the footer field
            Integer printNth = 1;
            if (cg.getVertices().size() > 9) {
                printNth = 1;
            }
            if (cg.getVertices().size() > 99) {
                printNth = 2;
            }
            */

            // line markers: S for sequence order, G for graph order
            Integer lineHeight = pl.textLineHeight;
            if (cg.getVertices().size() > 0) {
                ig2.drawString("C#", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (stringHeight / 4));
                ig2.drawString("CN", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + lineHeight + (stringHeight / 4));
                ig2.drawString("ML", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (lineHeight *2) + (stringHeight / 4));
                ig2.drawString("MN", pl.getFooterStart().x - pl.vertDist, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4));
            } else {
                ig2.drawString("(Graph has no vertices.)", pl.getFooterStart().x, pl.getFooterStart().y);
            }
            String[] vertexNameAssignment = cg.proteinNodeMap.toString().replace("{", "").replace("}", "").replace(" ", "").split(",");  // produces array of "x=a" where x is number of vertex and a chain name
            
            Iterator<Vertex> vertIter2 = cg.getVertices().iterator();  // this and next line used to iterate vertices to get current vertex number for MolNames
            Vertex curVert2;
            
            for (Integer i = 0; i < cg.getVertices().size(); i++) {
                curVert2 = vertIter2.next();
                // Draw each label until 999 and from then on only even ones
                if (i >= 999) {
                    if (i % 2 == 0) {
                        continue;
                    }
                }
                chainNumber = "" + (i + 1);
                //sseNumberSeq = "" + (cg.proteinNodeMap.get(i));

                // old version of getting chain name
                /*
                Integer foundIndex = edgesString.indexOf(i.toString() + "=");
                String chainId;
                if (i < 10) {
                    chainId = edgesString.substring(foundIndex + 2, foundIndex + 3);
                } else {
                    chainId = edgesString.substring(foundIndex + 3, foundIndex + 4);
                }
                chainName = "" + chainId; Vertex curVert2;
                */

                chainName = vertexNameAssignment[i].split("=")[1];
                

                //stringWidth = fontMetrics.stringWidth(sseNumberSeq);
                stringHeight = fontMetrics.getAscent();
                
                String molName = cg.molMap.get(curVert2);
                AffineTransform rotateMN = new AffineTransform();
                rotateMN.rotate(0.785d,0,0); // rotation around center of vertex
                Font rotatedFont = font.deriveFont(rotateMN);
                
               
                //Font font = new Font(null, Font.PLAIN, 10);
                

                ig2.drawString(chainNumber, pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (stringHeight / 4));
                ig2.drawString(chainName, pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 1) + (stringHeight / 4));
                ig2.drawString(molInfoForChains.get(chainName), pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));
                
               
                ig2.setFont(rotatedFont);

                ig2.drawString(molName, pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 3) + (stringHeight / 4));    
                
                ig2.setFont(font);
                
// determine chain of SSEs
                /*for(Integer x = 0; x < cg.getVertices().size(); x++){
                 if(i < cg.chainEnd.get(x)) {iChainID = x; break;}
                 }
                 */
                //if(iChainID != -1) {ig2.drawString(cg.allChains.get(iChainID).getPdbChainID(), pl.getFooterStart().x + (i * pl.vertDist) + pl.vertRadius / 2, pl.getFooterStart().y + (lineHeight * 2) + (stringHeight / 4));}
            }
            
            
            // footer
            final Rectangle2D key;
            String footer = "Homologous chains have the same color. Non-homologous chains are gray.";
            
            if (hasRna) {
                String footerAddon = "\nPeptide chains are represented as rectangles, RNA chains are represented as triangles.";
                footer += footerAddon;
                key = ig2.getFontMetrics().getStringBounds(footerAddon, ig2);
            } else {
                key = ig2.getFontMetrics().getStringBounds(footer, ig2);
            }

            int border = 10;
            
            // draw footer
            DrawTools.drawStringLineBreaks(ig2, footer, pl.getFooterStart().x - pl.vertDist - border, pl.getFooterStart().y + (pl.footerHeight - 40) - border);
            
            int numberOfLines = 1 + footer.length() - footer.replace("\n", "").length();
            
            // draw the rectangle around the footer
            ig2.drawRect(pl.getFooterStart().x - pl.vertDist - border - 5, pl.getFooterStart().y + (pl.footerHeight - 55 ) - border, (int) key.getWidth() + border, (int) key.getHeight() + numberOfLines * border);
            
            
            /*
             if(Settings.getBoolean("PTGLgraphComputation_B_graphimg_legend")) {
             if(iChainID != -1){
             SSEGraph.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + lineHeight * 3 + (stringHeight / 4)), pl, pg);
             }
             else{
             SSEGraph.drawLegend(ig2, new Position2D(pl.getFooterStart().x, pl.getFooterStart().y + lineHeight * 2 + (stringHeight / 4)), pl, pg);
             }
             }
             */
        }

    // all done, write the image to disk
        //if(Settings.get("PTGLgraphComputation_S_img_output_format").equals("SVG")) {
    //boolean useCSS = true;
        //FileOutputStream fos = new FileOutputStream(new File("/tmp/mySVG.svg"));
        //Writer out = new OutputStreamWriter(fos, "UTF-8");
        //ig2.stream(out, useCSS); 
        Rectangle2D roi = new Rectangle2D.Double(0, 0, pl.getPageWidth(), pl.getPageHeight());

        DrawResult drawRes = new DrawResult(ig2, roi);
        return drawRes;

    }

    /**
     * Draw a complex graph to an image file in the requested format.
     *
     * @param baseFilePathNoExt the img base file name, no file extension
     * @param drawBlackAndWhite whether to draw in grayscale only
     * @param formats a list of img formats to write
     * @param cg the complex graph to draw
     * @param molInfoForChains info mapping chain IDs (like "A") to their macromolecule (MOL_ID in PDB file, e.g., "1"). Give an empty one if you dont know
     * @return a list of file names that were written to disk, (as a map of
     * formats to file names)
     */
    public static HashMap<DrawTools.IMAGEFORMAT, String> drawComplexGraph(String baseFilePathNoExt, Boolean drawBlackAndWhite, DrawTools.IMAGEFORMAT[] formats, ComplexGraph cg, Map<String, String> molInfoForChains) {

        DrawResult drawRes = ComplexGraph.drawChainLevelComplexGraphG2D(drawBlackAndWhite, cg, molInfoForChains);

        //System.out.println("drawProteinGraph: Basefilepath is '" + baseFilePathNoExt + "'.");
        String svgFilePath = baseFilePathNoExt + ".svg";
        HashMap<DrawTools.IMAGEFORMAT, String> resultFilesByFormat = new HashMap<DrawTools.IMAGEFORMAT, String>();
        try {
            DrawTools.writeG2dToSVGFile(svgFilePath, drawRes);
            resultFilesByFormat.put(DrawTools.IMAGEFORMAT.SVG, svgFilePath);
            resultFilesByFormat.putAll(DrawTools.convertSVGFileToOtherFormats(svgFilePath, baseFilePathNoExt, drawRes, formats));
        } catch (IOException ex) {
            DP.getInstance().e("Could not write protein graph file : '" + ex.getMessage() + "'.");
        }

        if (!Settings.getBoolean("PTGLgraphComputation_B_silent")) {
            StringBuilder sb = new StringBuilder();
            sb.append("      Output complex graph files: ");
            for (DrawTools.IMAGEFORMAT format : resultFilesByFormat.keySet()) {
                String ffile = new File(resultFilesByFormat.get(format)).getName();
                sb.append("(").append(format.toString()).append(" => ").append(ffile).append(") ");
            }
            System.out.println(sb.toString());
        }
        return resultFilesByFormat;
    }

    /**
     * Writes this complex graph to the file 'file' in GML format. Note that
     * this function will overwrite the file if it exists.
     *
     * @param file the target file. Has to be writable.
     * @return true if the file was written, false otherwise
     */
    public boolean writeToFileGML(File file) {
        
        Boolean snakeCase = Settings.getBoolean("PTGLgraphComputation_B_gml_snake_case");
        
        // - - - Checks - - -

        if (!file.exists()) {
            try {
                file.createNewFile();
            } catch (IOException ex) {
                System.err.println("ERROR: Could not create file '" + file.getAbsolutePath() + "': " + ex.getMessage() + ".");
                return false;
            }
        }
        
        
        // - - - Create GML String with SPARGEL - - -
        
        // - - Vertex attributes - -

        GMLWriter<ComplexGraph.Vertex, ComplexGraph.Edge> gw = new GMLWriter<ComplexGraph.Vertex, ComplexGraph.Edge>(this);        
        gw.addVertexAttrWriter(new GMLWriter.AttrWriter<Vertex>() {
            @Override
            public String getAttribute() {
                return "label";
            }

            @Override
            public boolean hasValue(Vertex o) {
                return proteinNodeMap.containsKey(o);
            }

            @Override
            public String write(Vertex o) {
                return '"' + proteinNodeMap.get(o).toString() + '"';
            }
        });
        
        // add a new line per node holding the chain's molecule name
        gw.addVertexAttrWriter(new GMLWriter.AttrWriter<Vertex>() {
            @Override
            public String getAttribute() {
                return TextTools.formatAsCaseStyle(Arrays.asList("mol", "name"), snakeCase);
            }

            @Override
            public boolean hasValue(Vertex o) {
                return molMap.containsKey(o);
            }

            @Override
            public String write(Vertex o) {
                return '"' + molMap.get(o) + '"';
            }
        });
        
        // add a new line per node holding the chain's length
        gw.addVertexAttrWriter(new GMLWriter.AttrWriter<Vertex>() {
            @Override
            public String getAttribute() {
                return TextTools.formatAsCaseStyle(Arrays.asList("chain", "length"), snakeCase);
            }

            @Override
            public boolean hasValue(Vertex o) {
                return chainLengthMap.containsKey(o);
            }

            @Override
            public String write(Vertex o) {
                return chainLengthMap.get(o).toString();
            }
        });
        
        // add a new line per node holding the chain's radius of gyration
        gw.addVertexAttrWriter(new GMLWriter.AttrWriter<Vertex>() {
            @Override
            public String getAttribute() {
                return TextTools.formatAsCaseStyle(Arrays.asList("chain", "gyradius"), snakeCase);
            }

            @Override
            public boolean hasValue(Vertex o) {
                return chainGyradiusMap.containsKey(o);
            }

            @Override
            public String write(Vertex o) {
                return chainGyradiusMap.get(o).toString();
            }
        });
        
        // add a new line per node holding the chain's chemical type (RNA/peptide/non-polymer)
        gw.addVertexAttrWriter(new GMLWriter.AttrWriter<Vertex>() {
            @Override
            public String getAttribute() {
                return TextTools.formatAsCaseStyle(Arrays.asList("chain", "type"), snakeCase);
            }

            @Override
            public boolean hasValue(Vertex o) {
                return chainTypeMap.containsKey(o);
            }

            @Override
            public String write(Vertex o) {
                return '"' + chainTypeMap.get(o).toString() + '"';
            }
        });
        
        
        // - - Edge attributes - -

        //gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>());
        /**
         * Overwrite the edge attribute writer for the labels. This has to be
         * done because the labels need to be enclosed in closing quotation
         * marks to mark them as strings (due to the GML format definition).
         */
        gw.addEdgeAttrWriter(new GMLWriter.AttrWriter<Edge>() {
            @Override
            public String getAttribute() {
                return "label";
            }

            @Override
            public boolean hasValue(Edge e) {
                return mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).containsKey(e);
            }

            @Override
            public String write(Edge e) {
                return '"' + mapWeightNamesToMapEdgeValues.get(EdgeWeightType.ABSOLUTE_WEIGHT).get(e).toString() + '"';
            }

        });
       
        for (EdgeWeightType type : EdgeWeightType.values()) {
            gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList(type.name.split(" ")), snakeCase), mapWeightNamesToMapEdgeValues.get(type)));
        }
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "helix", "helix", "contacts"), snakeCase), numHelixHelixInteractionsMap));
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "helix", "strand", "contacts"), snakeCase), numHelixStrandInteractionsMap));
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "helix", "coil", "contacts"), snakeCase), numHelixCoilInteractionsMap));
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "strand", "strand", "contacts"), snakeCase), numStrandStrandInteractionsMap));
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "strand", "coil", "contacts"), snakeCase), numStrandCoilInteractionsMap));
        gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "coil", "coil", "contacts"), snakeCase), numCoilCoilInteractionsMap));
        if (! Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands")) { 
            gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "helix", "ligand", "contacts"), snakeCase), numHelixLigandInteractionsMap));
            gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "strand", "ligand", "contacts"), snakeCase), numStrandLigandInteractionsMap));
            gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "coil", "ligand", "contacts"), snakeCase), numCoilLigandInteractionsMap));
            gw.addEdgeAttrWriter(new GMLWriter.MapAttrWriter<>(TextTools.formatAsCaseStyle(Arrays.asList("num", "ligand", "ligand", "contacts"), snakeCase), numLigandLigandInteractionsMap));
        }
        
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try {
            gw.write(baos);
        } catch (IOException ex) {
            DP.getInstance().e(CLASS_TAG, "Writing the GML to a ByteArrayOutputStream failed. Please inform a developer. Trying to ignore it and going on.");
            Logger.getLogger(ComplexGraph.class.getName()).log(Level.SEVERE, null, ex);
        }
        String GmlString = new String(baos.toByteArray());
        
        
        // - - - Add graph meta data / graph attributes
        
        // SPARGEL's GMLWriter does not allow to add graph meta data / graph attributes, so we post process the string
        
        String[] GmlLines = GmlString.split("\n");
        String postprocessedGmlString = "";
        String lastLine = "";
        
        LinkedHashMap<String, String> graphAttributes = new LinkedHashMap<>();
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("creator"), snakeCase), "\"PTGLgraphComputation\"");
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("version"), snakeCase), "\"" + Settings.getVersion() + "\"");
        if (FileParser.getWhichSseInfo().equals("author")){
            graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("SseNote \"Created with author classification of SSEs. Results can't be compared to graphs of other proteins!\""), snakeCase), "");
        }
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("ignore", "ligands"), snakeCase), (Settings.getBoolean("PTGLgraphComputation_B_CG_ignore_ligands") ? "1" : "0"));  // whether ligands were ignored
        
        // documentation of edge weight types
        String innerPart = "[\n";
        for (EdgeWeightType type : EdgeWeightType.values()) {
            innerPart += "\t\t" + TextTools.formatAsCaseStyle(Arrays.asList(type.name.split(" ")), snakeCase) + " \"" + 
                    type.description + "\"\n";
        }
        innerPart += "\t]";
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("description", "edge", "weight", "types"), snakeCase), innerPart);
        
        // factor for lucid edge weight
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("min", "contacts", "for", "edge"), snakeCase), Settings.getInteger("PTGLgraphComputation_I_CG_contact_threshold").toString());  // contact threshold
        graphAttributes.put(TextTools.formatAsCaseStyle(Arrays.asList("factor", 
                TextTools.formatAsCaseStyle(Arrays.asList(EdgeWeightType.LUCID_MULTIPLICATIVE_NORMALIZATION.name.split(" ")), snakeCase)), snakeCase), 
                minimumMultiplicativeLengthNormalizedEdgeWeight.toString());  // factor to reconstruct normalized edge weight

        for (String GmlLine : GmlLines) {
            if (lastLine.equals("graph [")) {
                postprocessedGmlString += IO.mapToKeyValueString(graphAttributes, "\t", " ");
            }
            postprocessedGmlString += GmlLine + "\n";
            lastLine = GmlLine;
        }
        
        
        // - - - Write String to file - - -
        
        return IO.writeStringToFile(postprocessedGmlString, file.getAbsolutePath(), true);
    }
    
    
    /**
     * Retrieves the specified vertex ID from an edge. Method relies on the string representation of Edge.
     * @param edge
     * @param vertexNumber (positive) of the desired vertex. Most of the times 0 or 1.
     * @return vertex ID or null if invalid vertex number.
     */
    private Integer getVertexIdFromEdgeString(Edge edge, int vertexNumber) {
        if (vertexNumber < 0) {
            return null;
        }
        int[] vertexIDs = Arrays.stream(edge.toString()
                .replace("(", "").replace(")", "")
                .replace(" ", "")
                .split(","))
                .mapToInt(Integer::parseInt)
                .toArray();
        if (vertexNumber < vertexIDs.length) {
            return vertexIDs[vertexNumber];
        } else {
            return null;
        }
    }
    

    /**
     * Returns an integer array of simple format containing all edges with absolute weights.
     * @return Array of [index v1, index v2, weight]
     */
    private int[][] getEdgesAsArray() {
        EdgeWeightType weightType = EdgeWeightType.ABSOLUTE_WEIGHT;  // We use integer array, so other weight types make no sense
        int numEdges = mapWeightNamesToMapEdgeValues.get(weightType).size();
        int[][] result = new int[numEdges][3];  // v1, v2, weight
        int curEdgeIndex = 0;
        for (Edge edge : mapWeightNamesToMapEdgeValues.get(weightType).keySet()) {           
            result[curEdgeIndex] = new int[]{
                getVertexIdFromEdgeString(edge, 0),
                getVertexIdFromEdgeString(edge, 1),
                mapWeightNamesToMapEdgeValues.get(weightType).get(edge).intValue()
            };
            curEdgeIndex ++;
        }
        return result;
    }
    
    
    /**
     * 
     * @param vertexMap
     * @return 
     */
    private <T> Map<Integer, T> vertexMapToVertexIdMap(Map<Vertex, T> vertexMap) {
        Map<Integer, T> vertexIdMap = new HashMap<>();
        vertexMap.keySet().forEach(vertex -> {
            vertexIdMap.put(Integer.parseInt(vertex.toString()), vertexMap.get(vertex));
        });
        return vertexIdMap;
    }
    
    
    public ArrayList<String> getContactInfo() {
        if (createContactInfo) {
            return contactInfo;
        }
        else {
            DP.getInstance().w("Tried to get CG's contact info despite setting PTGLgraphComputation_B_writeComplexContactCSV was off. Returning null.");
            return null;
        }
    }
     
     
    public static Chain getChainByMacromolId(String inputId) {;
        for (Integer i = 0; i < allChains.size(); i++) {
            Chain curChain = allChains.get(i);
             if (curChain.getMacromolID().equals(inputId)) {
                return(curChain);
            }
        }
        return null;
    }
}
