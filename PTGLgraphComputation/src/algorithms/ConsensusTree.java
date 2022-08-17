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
 * Computes a consensus tree using bootstrapping based on several agglomerative clusters. 
 * Agglomerative clusters should be created using stochastic approach. 
 * @author christian
 */
public class ConsensusTree {
    
    private final ArrayList<ClusteringResult> results;
    private final ArrayList<Map<String, Integer>> resultsEdgeCounts;
    private final Map<String, Integer> edgeCount;
    
    public ConsensusTree(AgglomerativeClustering aggloClust, int iterations){
        this.results = new ArrayList<>();
        this.resultsEdgeCounts = new ArrayList<>();
        this.edgeCount = new HashMap<>();
        
        for (int i = 0; i < iterations; i++){
            this.results.add(clustering(aggloClust));
        }   
    }
    
    private ClusteringResult clustering(AgglomerativeClustering aggloClust){
        AgglomerativeClustering cluster = AgglomerativeClustering.newInstance(aggloClust);
        return cluster.chainLengthClustering();
    }
    
    private void addNewEdgeNode(String vertex, ArrayList<ArrayList<String>> subvertices){
        if(vertex.length() > 0){
            ArrayList<String> temp = new ArrayList<>();
            temp.add(vertex);
            subvertices.add(temp);
        }
    }
    
    private ArrayList<String> newEdge(String newick, int position, ArrayList<String> edges){
        ArrayList<ArrayList<String>> subvertices = new ArrayList<>(); 
        String vertex = new String();
        int lastSign = 0;
        
        for (int i=position; i < newick.length(); i++){
            Character act = newick.charAt(i);
            
            // needs to be outside of switch to interrupt the for loop
            if (act == ')' && lastSign == 0){
                addNewEdgeNode(vertex, subvertices);
                break;
            }
          
            if (lastSign == 0){
                switch(act){
                    case '(':   subvertices.add(newEdge(newick, i+1, edges));
                                lastSign++;
                    case ';':   break;
                
                    case ',':   addNewEdgeNode(vertex, subvertices);
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
        
        ArrayList<String> edgeVertices = new ArrayList<>();
        ArrayList<String> vertices = new ArrayList<>();
        for (ArrayList<String> list: subvertices){
            Collections.sort(list);
            edgeVertices.add(String.join(",", list));
            vertices.addAll(list);
        }
        
        Collections.sort(edgeVertices);
        String edge = String.join("-", edgeVertices);
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2){
            System.out.println("[DEBUG LV 3] Edge: " + edge);
        }
        
        edges.add(edge);
        return vertices;
    }
    
    private Map<String, Integer> countEdgesFromNewick(String newick){
        Map<String, Integer> localEdgeCount = new HashMap<>();
        ArrayList<String> edges = new ArrayList<>();
        
        // 1 since outer bracket isn't describing an edge
        newEdge(newick, 1, edges);
        
        for(String edge: edges){
            localEdgeCount.put(edge, 1);
        }
        
        return localEdgeCount;
    }
    
    private void countEdges(){
        for (ClusteringResult result: results){
            Map<String, Integer> count = countEdgesFromNewick(result.toNewickString());
            
            // combining counted hashmaps
            for(String key: count.keySet()){
                edgeCount.merge(key, count.get(key), (oldVal, newVal) -> oldVal+newVal);
            }
            
            // store counted edge for one ClusteringResult to reduce computations
            resultsEdgeCounts.add(count);
        }
    }
    
    public ClusteringResult computeConsensus(){
        countEdges();
        ClusteringResult consensus = null;
        int repValue = -1; 
        
        // evaluate
        for (int i=0; i < results.size(); i++){
            Map<String, Integer> resCount = resultsEdgeCounts.get(i);
            int resRepValue = 0; 
            for (String key: resCount.keySet()){
                resRepValue += edgeCount.get(key);    
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
    
    public void writeStatistics(String indicator){
        XMLOutputFactory xof = XMLOutputFactory.newInstance();
        XMLStreamWriter xsw = null;
        String sep = System.getProperty("file.separator");
                
        String savepath = Settings.get("PTGLgraphComputation_S_output_dir") + sep + indicator+"_clustStats.xml";
        try{
            xsw = xof.createXMLStreamWriter(new FileWriter(savepath));
            xsw.writeStartDocument();
            
            // 1) write edgeCount
            xsw.writeCharacters("\n");
            xsw.writeStartElement("count_table");
            xsw.writeCharacters("\n");
            
            for(String key: edgeCount.keySet()){
                xsw.writeCharacters("\t");
                xsw.writeStartElement("entry");
                xsw.writeAttribute("name", key);
                xsw.writeAttribute("bootstrap", Integer.toString(edgeCount.get(key)));
                xsw.writeEndElement();
                xsw.writeCharacters("\n");
            }
            xsw.writeEndElement();
            xsw.writeCharacters("\n");
            
            // 2) write resulting trees and edge list
            xsw.writeStartElement("clustering_results");
            xsw.writeAttribute("number", Integer.toString(results.size()));
            xsw.writeCharacters("\n");
            for(int i = 0; i < results.size(); i++){
                xsw.writeCharacters("\t");
                xsw.writeStartElement("result");
                xsw.writeAttribute("index", Integer.toString(i));
                xsw.writeAttribute("newick", results.get(i).toNewickString());
                xsw.writeCharacters("\n\t\t");
                xsw.writeStartElement("edge_list");
                for(String key: resultsEdgeCounts.get(i).keySet()){
                    xsw.writeCharacters("\n\t\t\t");
                    xsw.writeStartElement("entry");
                    xsw.writeAttribute("name", key);
                    xsw.writeAttribute("bootstrap", Integer.toString(edgeCount.get(key)));
                    xsw.writeEndElement();
                }
                xsw.writeCharacters("\n\t\t");
                xsw.writeEndElement();
                xsw.writeCharacters("\n\t");
                xsw.writeEndElement();
                xsw.writeCharacters("\n");
            }
            
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
        Settings.set("PTGLgraphComputation_B_type_assembly_prediction", "2");
        
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
        ConsensusTree ct = new ConsensusTree(clustering, 2);
        ClusteringResult cr = ct.computeConsensus();
        ct.writeStatistics("test");
        System.out.println(cr.toNewickString());
    }    
}
