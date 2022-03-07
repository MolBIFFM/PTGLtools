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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;
import proteingraphs.ComplexGraphEdgeWeightTypes.EdgeWeightType;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import tools.DP;
import settings.Settings;

/**
 *
 * @author jnw
 */
public class AgglomerativeClustering {
    private final static String CLASS_TAG = "CG";
    
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
     * Constructor of an agglomerative clustering. Class cannot be static so that it can hold instances of nested class Edge.
     * @param edgeList Array of [vertex ID 1, vertex ID 2, edge weight]
     * @param chainLengths
     * @param weightType
     */
    public AgglomerativeClustering(int[][] edgeList, Map<Integer, Integer> chainLengths, EdgeWeightType weightType, Map<Integer, String> labelMap) {
        curStepNum = 0;
        
        edges = new ArrayList<>();
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
     * 
     * @return 
     */
    public ClusteringResult chainLengthClustering() {
        ClusteringResult clusteringResult = new ClusteringResult(false);
    
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
            edges.forEach(edge -> {
                System.out.println("  " + edge.toString());
            });
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") > 1) {
            System.out.println("[DEBUG LV 2] Edge with highest normalized weight: " + edges.get(0).toString());
        }
        
        // choose edge to merge: greedy or interactive
        int mergeEdgeIndex;
        if (Settings.getBoolean("PTGLgraphComputation_B_interactive_assembly_prediction")) {
            // interactive
            int lineNum = 0;
            for (String line : prettyFormatEdgeList().split("\n")) {
                System.out.println(String.valueOf(lineNum) + ": " + line);
                lineNum++;
            }
        } else {
            // greedy
            mergeEdgeIndex = 0;
        }
        
        // add vertices of edges with highest normalized weight to merges of cluster result
        clusteringResult.addMerge(IntStream.of(edges.get(0).getVertices()).boxed().toArray(Integer[]::new));
        
        int v1 = edges.get(0).v1;
        int v2 = edges.get(0).v2;
        
        // 2.1) remove edge between merged vertices
        edges.remove(0);
        
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
     * Return a formatted string representation of the current edge list.
     * @return multi-line string
     */
    private String prettyFormatEdgeList() {
        String formattedString = "";
        for (Edge edge: edges) {
            formattedString += labelMap.get(edge.v1) + "-" + labelMap.get(edge.v2) + ", #res-res contacts: " + edge.absoluteWeight;
            if (weightType != EdgeWeightType.ABSOLUTE_WEIGHT) {
                formattedString += ", " + weightType.name + ": " + edge.normalizedWeight;
            }
            formattedString += "\n";
        }
        return formattedString;
    }
    
    /**
     * Tests agglomerative clustering implementation.
     */
    public static void main() {
        System.out.println("CLUSTERING TEST");
        Settings.set("PTGLgraphComputation_I_debug_level", "3");
        
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
