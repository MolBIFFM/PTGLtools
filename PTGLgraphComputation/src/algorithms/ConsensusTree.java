/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package algorithms;

import java.util.ArrayList;
import datastructures.ClusteringResult;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import settings.Settings;

/**
 * Computes a consensus tree using multiple stochastic clusterings. 
 * Agglomerative clusters should be created using stochastic approach 
 * (PTGLgraphComputation_I_type_assembly_prediction set to 2)
 * @author christian
 */
public class ConsensusTree {
    
    private final ArrayList<ClusteringResult> results;
    private final ArrayList<Map<String, Integer>> resultsCladeCount;
    private final Map<String, Integer> cladeSupport;
    private final Map<Integer, String> labelMap;
    
    /**
     * Constructor for the ConsensusTree - starts clustering directly. 
     * @param aggloClust AgglomerativeClustering - be sure to select the stochastic clustering
     * @param iterations Number of trees which should be generated
     * @param labelMap Labelmapping used for Newick-strings (to retain biological labels)
     */
    public ConsensusTree(AgglomerativeClustering aggloClust, int iterations, Map<Integer, String> labelMap){
        this.results = new ArrayList<>();
        
        //stores individual clade support for XML file
        this.resultsCladeCount = new ArrayList<>(); 
        this.cladeSupport = new HashMap<>();
        this.labelMap = labelMap;
        
        for (int i = 0; i < iterations; i++){
            this.results.add(clustering(aggloClust));
        }   
    }
    
    /**
     * Generates several clusterings based on copies of the given one. 
     * @param aggloClust clustering procedure - used as parameter to not store it separatly
     * @return ClusteringResults
     */
    private ClusteringResult clustering(AgglomerativeClustering aggloClust){
        AgglomerativeClustering cluster = AgglomerativeClustering.newInstance(aggloClust);
        return cluster.chainLengthClustering();
    }
    
    /**
     * Procedure to add nodes to a clade by creating a new ArrayList.
     * @param vertex Vertix which should be added
     * @param subvertices List of vertices belonging to a clade.
     */
    private void addNodeToClade(String vertex, ArrayList<ArrayList<String>> subvertices){
        if(vertex.length() > 0){
            ArrayList<String> temp = new ArrayList<>();
            temp.add(vertex);
            subvertices.add(temp);
        }
    }
    
    /**
     * Extracts clades from a newick string using recursion. Each clade will 
     * be sorted lexicographically. (each part of the split, as also each subgroup itself)
     * 
     * Method should work for non-binary trees - NOT TESTED JET
     * 
     * @param newick Newick string which should be parsed
     * @param position Position where to start (used for new clades)
     * @param clades List of clades which will be filled. 
     * @return Nodes which belong to the calculated clade.
     */
    private ArrayList<String> newClade(String newick, int position, ArrayList<String> clades){
        ArrayList<ArrayList<String>> subvertices = new ArrayList<>(); 
        String vertex = new String();
        int lastSign = 0;
        
        for (int i=position; i < newick.length(); i++){
            Character act = newick.charAt(i);
            
            // needs to be outside of switch to interrupt the for loop
            if (act == ')' && lastSign == 0){
                addNodeToClade(vertex, subvertices);
                break;
            }
          
            if (lastSign == 0){
                switch(act){
                    case '(':   subvertices.add(newClade(newick, i+1, clades));
                                lastSign++;
                    case ';':   break;
                
                    case ',':   addNodeToClade(vertex, subvertices);
                                vertex = new String();
                                break;
                    
                    default:    vertex += act;
                                break;  
                }
            }
            else{
                switch(act){
                    case '(':   lastSign++;
                                break;

                    case ')':   if(lastSign != 0){lastSign--;} break;

                    default:    break;
                }
            }
        }
        
        ArrayList<String> cladeVertices = new ArrayList<>();
        ArrayList<String> vertices = new ArrayList<>();
        for (ArrayList<String> list: subvertices){
            Collections.sort(list);
            cladeVertices.add(String.join(",", list));
            vertices.addAll(list);
        }
        
        // each ArrayList in cladeVertices descripes one subset 
        Collections.sort(cladeVertices);
        String clade = String.join("-", cladeVertices);
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2){
            System.out.println("[DEBUG LV 3] Edge: " + clade);
        }
        
        clades.add(clade);
        return vertices;
    }
    
    /**
     * Counts the clades within a newick string.
     * Only tested for binary trees. 
     * @param newick newick-string of a tree
     * @return HashMap of clade counts (keys = clades as string)
     */
    private Map<String, Integer> countCladesFromNewick(String newick){
        Map<String, Integer> localEdgeCount = new HashMap<>();
        ArrayList<String> clades = new ArrayList<>();
        
        // staring recursion 1 since outer bracket is already considered at this point
        newClade(newick, 1, clades);
        
        for(String clade: clades){
            // each clade will just be seen ones - no need for checking collisions
            localEdgeCount.put(clade, 1); 
        }
        
        return localEdgeCount;
    }
    
    /**
     * Starts counting clades from generated clusterings.
     * Fills cladeSupport and resultsCladeCount with the appropriate values.
     */
    private void countClades(){
        for (ClusteringResult result: results){
            Map<String, Integer> count = countCladesFromNewick(result.toNewickString(labelMap));
            
            // combining counted hashmaps
            for(String key: count.keySet()){
                cladeSupport.merge(key, count.get(key), (oldVal, newVal) -> oldVal+newVal);
            }
            
            // store counted edge for one ClusteringResult to reduce computations
            resultsCladeCount.add(count);
        }
    }
    
    /**
     * Calculates consensus tree based on Clade Support (CS) 
     * @return ClusteringResult from consensus
     */
    public ClusteringResult computeConsensus(){
        countClades();
        ClusteringResult consensus = null;
        int repValue = -1; 
        
        // evaluate
        for (int i=0; i < results.size(); i++){
            Map<String, Integer> resCount = resultsCladeCount.get(i);
            int resRepValue = 0; 
            for (String key: resCount.keySet()){
                resRepValue += cladeSupport.get(key);    
            }
            
            if (resRepValue > repValue){
                consensus = results.get(i); 
                repValue = resRepValue;
            }
            
            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1){
                System.out.println("[DEBUG LV 2] Result:" + results.get(i).toNewickString());
                System.out.println("[DEBUG LV 2] Reresentation value: " + resRepValue);
            }
        }

        return consensus; 
    }
    
    /**
     * Creates XML-files to store statistics in the output directory.
     * Consensus needs to be calculated first, otherwise the file will be empty. 
     * Includes whole CS-map of the whole forest, as also each clustering result
     * in particular. 
     * @param indicator filename extension (e.g. edge-weight type)
     */
    public void writeStatistics(String indicator){
        XMLOutputFactory xof = XMLOutputFactory.newInstance();
        XMLStreamWriter xsw = null;
        String sep = System.getProperty("file.separator");
                
        String savepath = Settings.get("PTGLgraphComputation_S_output_dir") + sep + indicator+"_clustStats.xml";
        try{
            xsw = xof.createXMLStreamWriter(new FileWriter(savepath));
            xsw.writeStartDocument();
            xsw.writeCharacters("\n");
            xsw.writeStartElement("stochastic_clustering");
            // 1) write edgeCount
            xsw.writeCharacters("\n\t");
            xsw.writeStartElement("count_table");
            xsw.writeCharacters("\n");
            
            for(String key: cladeSupport.keySet()){
                xsw.writeCharacters("\t\t ");
                xsw.writeStartElement("entry");
                xsw.writeAttribute("name", key);
                xsw.writeAttribute("clade_support", Integer.toString(cladeSupport.get(key)));
                xsw.writeEndElement();
                xsw.writeCharacters("\n");
            }
            xsw.writeCharacters("\t");
            xsw.writeEndElement();
            xsw.writeCharacters("\n\t");
            
            // 2) write resulting trees and edge list
            xsw.writeStartElement("clustering_results");
            xsw.writeAttribute("number", Integer.toString(results.size()));
            xsw.writeCharacters("\n\t");
            for(int i = 0; i < results.size(); i++){
                xsw.writeCharacters("\t");
                xsw.writeStartElement("result");
                xsw.writeAttribute("index", Integer.toString(i));
                xsw.writeCharacters("\n");
                results.get(i).writeXMLPresentation(labelMap, xsw);        
                xsw.writeCharacters("\n\t\t\t");
                xsw.writeStartElement("edge_list");
                for(String key: resultsCladeCount.get(i).keySet()){
                    xsw.writeCharacters("\n\t\t\t\t");
                    xsw.writeStartElement("entry");
                    xsw.writeAttribute("name", key);
                    xsw.writeAttribute("clade_support", Integer.toString(cladeSupport.get(key)));
                    xsw.writeEndElement();
                }
                xsw.writeCharacters("\n\t\t\t");
                xsw.writeEndElement();
                xsw.writeCharacters("\n\t\t");
                xsw.writeEndElement();
                xsw.writeCharacters("\n\t");
            }
            xsw.writeEndElement();
            xsw.writeCharacters("\n");
            xsw.writeEndElement();
            xsw.writeEndDocument();
        
        }
        catch (IOException | XMLStreamException e){
                System.err.println("Unable to write clustering statistics: " + e.getMessage());
        }
        finally{
            try{
                if (xsw != null){
                        xsw.close();
                }
            }
            catch (XMLStreamException e){
                    System.err.println("Unable to close clustering statistics: " + e.getMessage());
            }
        }
       
    }
    
    // Test
    public static void main(String args[]){
        Settings.set("PTGLgraphComputation_I_debug_level", "3");
        Settings.set("PTGLgraphComputation_I_type_assembly_prediction", "3");
        
        int[][] edgeList = {
            {0, 1, 5},
            {0, 2, 4},
            {0, 3, 3},
            {0, 4, 10},
            {0, 5, 8},
            {0, 6, 7},
            {1, 2, 15},
            {2, 3, 14},
            {5, 7, 18},
            {5, 8, 21},
            {6, 7, 12}
        };
        
        Map<Integer, Integer> chainLengths = new HashMap<>();
        chainLengths.put(0, 40);
        chainLengths.put(1, 8);
        chainLengths.put(2, 10);
        chainLengths.put(3, 20);
        chainLengths.put(4, 5);
        chainLengths.put(5, 12);
        chainLengths.put(6, 15);
        chainLengths.put(7, 11);
        chainLengths.put(8, 7);
        
        Map<Integer, String> labelMap = new HashMap<>();
        for (Integer vertexId : chainLengths.keySet()) {
            labelMap.put(vertexId, String.valueOf(vertexId));
        }
        
        AgglomerativeClustering clustering = new AgglomerativeClustering(edgeList, chainLengths, ComplexGraphEdgeWeightTypes.EdgeWeightType.ADDITIVE_LENGTH_NORMALIZATION, labelMap);
        ConsensusTree ct = new ConsensusTree(clustering, 1, labelMap);
        ClusteringResult cr = ct.computeConsensus();
        ct.writeStatistics("test");
        System.out.println(cr.toNewickString());
    }    
}
