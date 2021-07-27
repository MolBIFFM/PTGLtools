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
import java.util.Collections;
import java.util.Comparator;
import java.util.Map;
import proteingraphs.ComplexGraphEdgeWeightTypes.EdgeWeightTypes;
import proteingraphs.ComplexGraphEdgeWeightTypes;
import tools.DP;

/**
 *
 * @author jnw
 */
public class AgglomerativeClustering {
    private final static String CLASS_TAG = "CG";
    
    private ArrayList<Edge> edges;
    
    class Edge implements Comparable<Edge> {
        int v1;
        int v2;
        int absoluteWeight;
        BigDecimal normalizedWeight;
        EdgeWeightTypes normalizationType;
        
        Edge (int v1, int v2, int absoluteWeight, EdgeWeightTypes normalizationType, BigDecimal normalizedWeight) {
            this.v1 = v1;
            this.v2 = v2;
            this.absoluteWeight = absoluteWeight;
            this.normalizationType = normalizationType;
            this.normalizedWeight = normalizedWeight;            
        }
        
        @Override
        public int compareTo(Edge otherEdge) {
            return normalizedWeight.compareTo(otherEdge.normalizedWeight);
        }
        
        @Override
        public String toString() {
            return "Edge {" + v1 + "," + v2 + "}, weight type: " + normalizationType.name + ", normalized weight: " + normalizedWeight.toString();
        }
    }
    
    class EdgeComparatorNormalizedWeight implements Comparator<Edge> {
        @Override
        public int compare(Edge e1, Edge e2) {
            return e1.compareTo(e2);
        }
    }
    
        /**
     * 
     * @param edgeList Array of [vertex ID 1, vertex ID 2, edge weight]
     * @param chainLengths
     * @param weightType
     * @return 
     */
    public AgglomerativeClustering(int[][] edgeList, Map<Integer, Integer> chainLengths, EdgeWeightTypes weightType) {
        edges = new ArrayList<>();
        for (int[] vertexIDs : edgeList) {
            if (vertexIDs.length < 3) {
                DP.getInstance().e(CLASS_TAG, "Got an edge without two particpating vertices. Skipping edge and trying to go on. "
                        + "Please report this bug to the developers.");
                continue;
            }
            BigDecimal normalizedWeight;
            
            switch (weightType) {
                case ADDITIVE_LENGTH_NORMALIZATION:
                    normalizedWeight = ComplexGraphEdgeWeightTypes.computeAdditiveLengthNormlization(BigDecimal.valueOf(vertexIDs[2]), chainLengths.get(vertexIDs[0]), chainLengths.get(vertexIDs[1]));
                    break;
                case MULTIPLICATIVE_LENGTH_NORMALIZATION:
                    normalizedWeight = ComplexGraphEdgeWeightTypes.computeMultiplicativeLengthNormlization(BigDecimal.valueOf(vertexIDs[2]), chainLengths.get(vertexIDs[0]), chainLengths.get(vertexIDs[1]));
                    break;
                case ABSOLUTE_WEIGHT:
                default:
                    normalizedWeight = BigDecimal.valueOf(vertexIDs[2]);
            }
            
            edges.add(new Edge(vertexIDs[0], vertexIDs[1], vertexIDs[2], weightType, normalizedWeight));
        }
    }
    
    
    public ClusteringResult chainLengthClustering() {
        ClusteringResult clustering = new ClusteringResult(false);

        // TODELETE
        System.out.println("Edges before sorting");
        edges.forEach(edge -> {
            System.out.println(edge.toString());
        });

        // TODO sort
        Collections.sort(edges, Collections.reverseOrder(new EdgeComparatorNormalizedWeight()));
        
        // TODELETE
        System.out.println("\nEdges after sorting");
        edges.forEach(edge -> {
            System.out.println(edge.toString());
        });
        
        
        return clustering;
    }
    
    private void chainLengthClusteringStep(ArrayList<Edge> edges) {
        
    }
}
