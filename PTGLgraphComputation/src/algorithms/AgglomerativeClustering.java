/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package algorithms;

import datastructures.ClusteringResult;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;
import java.util.Scanner;
import java.util.Random; 
import static java.util.stream.Collectors.toList;
import proteingraphs.ComplexGraphEdgeWeightTypes.EdgeWeightType;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import tools.DP;
import settings.Settings;

/**
 *
 * @author jnw
 * modified by Christian Ickes
 */
public class AgglomerativeClustering {
    private final static String CLASS_TAG = "CG";
    private static final Random RANDOM = new Random(); 
    
    private ArrayList<Edge> edges;
    private final Map<Integer, Integer> chainLengths;
    private final EdgeWeightType weightType;
    private final Map<Integer, String> labelMap;
    private int curStepNum;
    
    class Edge implements Comparable<Edge> {
        int v1;  // always lower ID
        int v2;  // always greater ID
        private int absoluteWeight;
        private BigDecimal normalizedWeight;
        private final EdgeWeightType normalizationType;
        
        Edge (int v1, int v2, int absoluteWeight, EdgeWeightType normalizationType, BigDecimal normalizedWeight) {
            this.v1 = Math.min(v1, v2);
            this.v2 = Math.max(v1, v2);
            this.absoluteWeight = absoluteWeight;
            this.normalizationType = normalizationType;
            this.normalizedWeight = normalizedWeight;            
        }
        
        public Edge mergeEdgeIntoThis(Edge otherEdge, int lengthChainA, int lengthChainB) {
            if (this.normalizationType != otherEdge.normalizationType) {
                DP.getInstance().e(CLASS_TAG, "Cannot merge edges of different normalization. Returning first edge and trying to got on. "
                        + "Please report this bug to the developers.");
            }
            absoluteWeight += otherEdge.absoluteWeight;
            normalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(absoluteWeight, lengthChainA, lengthChainB, normalizationType);
            return this;
        }
        
        @Override
        public int compareTo(Edge otherEdge) {
            return normalizedWeight.compareTo(otherEdge.normalizedWeight);
        }
        
        @Override
        public int hashCode() {
            return Objects.hash(v1, v2);
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final Edge other = (Edge) obj;
            if (this.v1 != other.v1) {
                return false;
            }
            return this.v2 == other.v2;
        }
        
        @Override
        public String toString() {
            return "Edge {" + v1 + "," + v2 + "}, absolute weight: " + absoluteWeight + " normalization: " + normalizationType.name + ", normalized weight: " + normalizedWeight.toString();
        }
        
        public int[] getVertices() {
            int[] vertices = {v1, v2};
            return vertices;
        }
        
        public int getAbsoluteWeight() {
            return absoluteWeight;
        }
    }
    
    class EdgeComparatorNormalizedWeight implements Comparator<Edge> {
        @Override
        public int compare(Edge e1, Edge e2) {
            return e1.compareTo(e2);
        }
    }
    
    /**
     * Creates a new instance of the clustering to copy it completly. 
     * @param aggloClust
     * @return complete new instance 
     */
    public static AgglomerativeClustering newInstance(AgglomerativeClustering aggloClust){
        return new AgglomerativeClustering(aggloClust);
    }
    
    /**
     * Constructor of an agglomerative clustering.Class cannot be static so that it can hold instances of nested class Edge.
     * @param edgeList Array of [vertex ID 1, vertex ID 2, edge weight]
     * @param chainLengths
     * @param weightType
     * @param labelMap
     */   
    public AgglomerativeClustering(int[][] edgeList, Map<Integer, Integer> chainLengths, EdgeWeightType weightType, Map<Integer, String> labelMap) {
        this.curStepNum = 0;
        
        this.edges = new ArrayList<>();
        for (int[] edge : edgeList) {
            if (edge.length < 3) {
                DP.getInstance().e(CLASS_TAG, "Got an edge without two particpating vertices. Skipping edge and trying to go on. "
                        + "Please report this bug to the developers.");
                continue;
            } 
            BigDecimal normalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(edge[2], chainLengths.get(edge[0]), chainLengths.get(edge[1]), 
                    weightType);

            edges.add(new Edge(edge[0], edge[1], edge[2], weightType, normalizedWeight));
        }
        
        this.chainLengths = new HashMap<>(chainLengths);  // Shallow copy to not change chain lengths outside of this class
        this.weightType = weightType;
        this.labelMap = labelMap;
    }
    
    /**
     * Used to create copy by newInstance
     * @param aggloClust 
     */
    public AgglomerativeClustering(AgglomerativeClustering aggloClust){
        this.edges = new ArrayList<>(aggloClust.edges);
        this.chainLengths = new HashMap<>(aggloClust.chainLengths);
        this.weightType = aggloClust.weightType;
        this.labelMap = new HashMap<>(aggloClust.labelMap);
        this.curStepNum = 0;
    }
            
    /**
     * Performs a hirachical clustering algorithm on the given complex graph. 
     * @return Clustered Tree
     */
    public ClusteringResult chainLengthClustering() {
        ClusteringResult clusteringResult = new ClusteringResult(true);
    
        while (edges.size() > 0) {
            chainLengthClusteringStep(clusteringResult);
            
            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1) {
                System.out.println("[DEBUG LV 2] Merges after step " + curStepNum);
                System.out.println(clusteringResult.toString());
            }
        }
   
        return clusteringResult;
    }
    
    private void chainLengthClusteringStep(ClusteringResult clusteringResult) {
        // 1) sort by normalized weight
        Collections.sort(edges, Collections.reverseOrder(new EdgeComparatorNormalizedWeight()));
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2) {
            System.out.println("[DEBUG LV 3] Sorted edge list from step " + curStepNum);
            System.out.println("[DEBUG LV 3] Assembly type: " + 
                                Settings.getInteger("PTGLgraphComputation_B_type_assembly_prediction"));
            edges.forEach(edge -> {
                System.out.println("  " + edge.toString());
            });
            
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1) {
            System.out.println("[DEBUG LV 2] Edge with highest normalized weight: " + edges.get(0).toString());
        }
        
        // choose edge to merge: greedy or interactive
        int mergeEdgeIndex;
        switch(Settings.getInteger("PTGLgraphComputation_B_type_assembly_prediction")){
            case(1): mergeEdgeIndex = getEdgeIndexFromUser(6, clusteringResult); break;
            case(2): mergeEdgeIndex = getEdgeIndexStochastic(); break;
            default: mergeEdgeIndex = 0; break;
        }
        
        // add vertices of edges to merges of cluster result
        Edge curEdge = edges.get(mergeEdgeIndex);
        clusteringResult.addMerge(
                IntStream.of(curEdge.getVertices()).boxed().toArray(Integer[]::new), 
                curEdge.absoluteWeight,
                chainLengths.get(curEdge.v1),
                chainLengths.get(curEdge.v2));
        
        int v1 = curEdge.v1;
        int v2 = curEdge.v2;
        
        // 2.1) remove edge between merged vertices
        edges.remove(mergeEdgeIndex);
        
        // 2.2) update chain lengths
        int combinedChainLength = chainLengths.get(v1) + chainLengths.get(v2);
        chainLengths.replace(v1, combinedChainLength);  // v1 is representative
        chainLengths.remove(v2);

        // 2.3 update edges
        HashMap<Edge, Edge> edgeMap = new HashMap<>();
        edges.forEach(edge -> {
            Edge updatedEdge = edge;  // Set it to edge and update below if necessary
            if (IntStream.of(edge.getVertices()).anyMatch(v -> v == v2 || v == v1)) {
                // possibly update indices
                int updatedV1 = (edge.getVertices()[0] == v2 ? v1 : edge.getVertices()[0]);
                int updatedV2 = (edge.getVertices()[1] == v2 ? v1 : edge.getVertices()[1]);
                int absoluteWeight = edge.getAbsoluteWeight();

                BigDecimal updatedNomalizedWeight = ComplexGraphEdgeWeightTypes.computeLengthNormalization(absoluteWeight, chainLengths.get(updatedV1), chainLengths.get(updatedV2), 
                    weightType);
                updatedEdge = new Edge(updatedV1, updatedV2, absoluteWeight, weightType, updatedNomalizedWeight);
            } 
            
            // try to add updated edge
            if (edgeMap.containsKey(updatedEdge)) {
                // duplicate edge -> merge the two edges
                edgeMap.get(updatedEdge).mergeEdgeIntoThis(updatedEdge, chainLengths.get(updatedEdge.v1), chainLengths.get(updatedEdge.v2));
            } else {
                edgeMap.put(updatedEdge, updatedEdge);
            }
        });
        
        edges = new ArrayList<>(edgeMap.values());    
        curStepNum++;
    }
    
    
    /**
     * Print a table of all edges and ask user for an input until it is a valid edge index.
     * @return valid edge index as int
     */
    private int getEdgeIndexFromUser(int numberIndentSpaces, ClusteringResult clusteringResult) {
        Integer edgeIndex = null;
        Scanner lineReader = new Scanner(System.in);
        System.out.println(" ".repeat(numberIndentSpaces) + "Step " + curStepNum + ": Please choose an edge for merging");
        for (String line : prettyFormatEdgeList(numberIndentSpaces + 2, clusteringResult).split("\n")) {
            System.out.println(line);
        }
        while (edgeIndex == null) {
            String input = lineReader.nextLine();
            try {
                edgeIndex = Integer.parseInt(input.trim());
                if (edgeIndex < 0 || edgeIndex > edges.size() - 1) {
                    System.out.println("Please choose a number between zero and " + (edges.size() - 1));
                    edgeIndex = null;
                }
            } catch (NumberFormatException e) {
                System.out.println("Please provide a valid number");
            }
        }
        return edgeIndex;
    }
    
    /**
     * Chooses an Edge by performing a roulett algorithm.
     * Probabilities are equal to edge weights. 
     * @return index of choosen edge in array. 
     */
    private int getEdgeIndexStochastic(){
        BigDecimal summed = new BigDecimal(0);
        for (Edge edge: edges){
            summed = summed.add(edge.normalizedWeight);
        }
        
        int edgeIndex = 0; 
        // generate a random BigDecimal in a range: max * random
        BigDecimal choice = summed.multiply(new BigDecimal(RANDOM.nextDouble()));
        
        summed = edges.get(0).normalizedWeight;
        // do binary search
        while(summed.compareTo(choice) < 0){
            summed = summed.add(edges.get(edgeIndex).normalizedWeight);
            edgeIndex++;
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 2) {
            System.out.println("[DEBUG LV 3] Summed Weight: " + summed.toString());
            System.out.println("[DEBUG LV 3] Chosen Weight: " + choice.toString());
            System.out.println("[DEBUG LV 3] Chosen Index: " + edgeIndex);
        }
        return edgeIndex; 
    }
    
    
    /**
     * Return a formatted string representation of the current edge list.
     * @return multi-line string
     */
    private String prettyFormatEdgeList(int numberIndentSpaces, ClusteringResult clusteringResult) {
        String formattedString = " ".repeat(numberIndentSpaces + 4) + "rounded normalized weight | #res-res contacts | sub assembly 1 - sub assembly 2\n";
        Integer lineNum = 0;
        for (Edge edge: edges) {
            formattedString += " ".repeat(numberIndentSpaces + 2) + lineNum.toString() + ": " + 
                    edge.normalizedWeight.setScale(4, RoundingMode.HALF_UP) + " | " + 
                    edge.absoluteWeight + " | " +
                    clusteringResult.getVerticesFromRepresentative(edge.v1).stream().map(ele -> labelMap.get(ele)).collect(toList()) + " - " +
                    clusteringResult.getVerticesFromRepresentative(edge.v2).stream().map(ele -> labelMap.get(ele)).collect(toList()) + "\n";
            lineNum++;
        }
        return formattedString;
    }
    
    /**
     * Tests agglomerative clustering implementation.
     */
    public static void main(String[] args) {
        System.out.println("CLUSTERING TEST");
        Settings.set("PTGLgraphComputation_I_debug_level", "3");
        Settings.set("PTGLgraphComputation_B_type_assembly_prediction", "2");
        
        // Example 1
        //   Additive length normalization merges
        //     2,3
        //     0,1
        //     0,2
        /*
        int[][] edgeList = {
            {0, 1, 10},
            {0, 2, 10},
            {1, 2, 5},
            {2, 3, 60}
        };
        
        Map<Integer, Integer> chainLengths = new HashMap<>();
        chainLengths.put(0, 5);
        chainLengths.put(1, 5);
        chainLengths.put(2, 15);
        chainLengths.put(3, 10);
        */
        
        // Example 2
        //   Additive length normalization merges
        //     5,8
        //     1,2
        //     5,7
        //     1,3
        //     5,6
        //     0,4
        //     0,5
        //     0,1
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

        AgglomerativeClustering clustering = new AgglomerativeClustering(edgeList, chainLengths, EdgeWeightType.ADDITIVE_LENGTH_NORMALIZATION, labelMap);
        
        ClusteringResult cr = clustering.chainLengthClustering();        
        
        System.out.println(cr.toNewickString());
        
        clustering = new AgglomerativeClustering(edgeList, chainLengths, EdgeWeightType.MULTIPLICATIVE_LENGTH_NORMALIZATION, labelMap);
        clustering.chainLengthClustering();
        
        clustering = new AgglomerativeClustering(edgeList, chainLengths, EdgeWeightType.ABSOLUTE_WEIGHT, labelMap);
        clustering.chainLengthClustering();
        
        System.out.println("FINISHED CLUSTERING TEST AND EXITING");
        System.exit(0);
    }
}
