/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package algorithms;

import java.util.ArrayList;
import datastructures.ClusteringResult;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import settings.Settings;

/**
 *
 * @author christian
 */
public class ConsensusTree {
    
    private final ArrayList<ClusteringResult> results;
    private final Map<Integer, ArrayList<String>> computedEdges;
    
    public ConsensusTree(AgglomerativeClustering aggloClust, int iterations){
        this.results = new ArrayList<>();
        this.computedEdges = new HashMap<>();
        
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
        Map<String, Integer> edgeCount = new HashMap<>();
        ArrayList<String> edges = new ArrayList<>();
        
        // 1 since outer bracket isn't describing an edge
        newEdge(newick, 1, edges);
        
        for(String edge: edges){
            edgeCount.put(edge, 1);
        }
        
        return edgeCount;
    }
    
    private Map<String, Integer> countEdges(ArrayList<Map<String, Integer>> resultEdgeCounts){
        Map<String, Integer> edgeCount = new HashMap<>();
        for (ClusteringResult result: results){
            Map<String, Integer> count = countEdgesFromNewick(result.toNewickString());
            
            // combining counted hashmaps
            for(String key: count.keySet()){
                edgeCount.merge(key, count.get(key), (oldVal, newVal) -> oldVal+newVal);
            }
            
            // store counted edge for one ClusteringResult to reduce computations
            resultEdgeCounts.add(count);
        }
        
        return edgeCount;
    }
    
    public ClusteringResult computeConsensus(){
        ArrayList<Map<String, Integer>> resultEdgeCounts = new ArrayList<>();
        Map<String, Integer> edgeCount = countEdges(resultEdgeCounts);  
        ClusteringResult consensus = null;
        int repValue = -1; 
        
        // evaluate
        for (int i=0; i < results.size(); i++){
            Map<String, Integer> resCount = resultEdgeCounts.get(i);
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
        ConsensusTree ct = new ConsensusTree(clustering, 1);
        ClusteringResult cr = ct.computeConsensus();
        System.out.println(cr.toNewickString());
    }    
}
