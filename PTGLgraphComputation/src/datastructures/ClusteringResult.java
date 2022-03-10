/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package datastructures;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import tools.DP;
import proteingraphs.ComplexGraphEdgeWeightTypes;

/**
 *
 * @author jnw
 */
public class ClusteringResult {
    private final static String CLASS_TAG = "CG";
    
    private final ArrayList<Integer[]> merges = new ArrayList<>();
    private final Map<Integer, ArrayList<Integer>> mapRepresentativeToVertices = new HashMap<>();
    private BigDecimal consecutiveLargeInterfaceScore = BigDecimal.ZERO;
    private final Boolean binary;
    private int maxNumberMergedAtOnce = 0;
    
    public ClusteringResult(Boolean binary) {
        this.binary = binary;
    }
    
    public Boolean addMerge(Integer[] merge, BigDecimal weight) {
        if (binary) {
            if (merge.length > 2) {
                DP.getInstance().e(CLASS_TAG, "Tried adding a merge of more than two vertices to a binary clustering result. Rejecting merge and trying to go on. "
                        + "Please report this bug to the developers.");
            }
        }
        merges.add(merge);
        
        // update representative map
        ArrayList<Integer> leafs = new ArrayList<>();
        for (Integer vertexID : merge) {
            if (mapRepresentativeToVertices.containsKey(vertexID)) {
                leafs.addAll(mapRepresentativeToVertices.get(vertexID));
            } else {
                leafs.add(vertexID);
            }
        }
        Collections.sort(leafs);
        mapRepresentativeToVertices.put(merge[0], leafs);
        
        // update score
        consecutiveLargeInterfaceScore = consecutiveLargeInterfaceScore.add(scoreFromThisMerge(weight, leafs.size()));
        
        maxNumberMergedAtOnce = Math.max(maxNumberMergedAtOnce, merge.length);
        return true;
    }
    
    
    private BigDecimal scoreFromThisMerge(BigDecimal weight, int numberLeafsUnderMerge) {
        return weight.divide(BigDecimal.valueOf(numberLeafsUnderMerge), ComplexGraphEdgeWeightTypes.PRECISION, RoundingMode.HALF_UP);
    }
    
    
    private ArrayList<String[]> renameVertexIdsInMerges(Map<Integer, String> mapVertexIdToName) {
        ArrayList<String[]> renamedMerges = new ArrayList<>();
        for (Integer[] merge : merges) {
            String[] renamedMerge = new String[merge.length];
            for (int vertexIndex = 0; vertexIndex < merge.length; vertexIndex++) {
                if (mapVertexIdToName.containsKey(merge[vertexIndex])) {
                    renamedMerge[vertexIndex] = mapVertexIdToName.get(merge[vertexIndex]);
                } else {
                    renamedMerge[vertexIndex] = Integer.toString(merge[vertexIndex]);
                }
            }
            renamedMerges.add(renamedMerge);
        }
        return renamedMerges;
    }
       
    /**
     * Returns a list of arrays as string, print of one array per line.
     * @param <T> inner array
     * @param list of arrays
     * @return string representation of each array per line
     */
    private <T> String listOfArrayToString(ArrayList<T[]> list) {
        String result = "";
        for (T[] innerArr : list) {
            result += Arrays.toString(innerArr);
            result += "\n";
        }
        return result.strip();  // strip to remove last newline
    }
    
    
    @Override
    public String toString() {
        return listOfArrayToString(merges);
    }
    
    /**
     * Same as toString but replacing vertex IDs by Strings of given map.
     * @param mapVertexIdToName maps integer vertex IDs to Strings, e.g., name
     * @return 
     */
    public String toString(Map<Integer, String> mapVertexIdToName) {
        return listOfArrayToString(renameVertexIdsInMerges(mapVertexIdToName));
    }
    
    
    /**
     * Returns the clustering result as Newick string of format 9.
     * @param <T>
     * @param merges as ArrayList of arrays of arbitrary type
     * @return clustering result as Newick string of format 9
     */
    private <T> String toNewickString(ArrayList<T[]> merges) {
        Map<T, String> subStrings = new HashMap<>();
        T lastRepresentative = null;
        for (T[] merge : merges) {
            lastRepresentative = merge[0];  //  first vertex is always representative
            ArrayList<String> verticeStrings = new ArrayList<>();
            for (T vertex : merge) {
                if (subStrings.containsKey(vertex)) {
                    verticeStrings.add(subStrings.get(vertex));
                    subStrings.remove(vertex);
                } else {
                    verticeStrings.add(vertex.toString());
                }
            }
            subStrings.put(merge[0], "(" + String.join(",", verticeStrings) + ")");
        }
        return subStrings.get(lastRepresentative) + ";";
    }
    
    public String toNewickString(Map<Integer, String> mapVertexIdToName) {
        return toNewickString(renameVertexIdsInMerges(mapVertexIdToName));
    }
    
    public String toNewickString() {
        return toNewickString(merges);
    }
    
    public ArrayList<Integer> getVerticesFromRepresentative(Integer vertexID) {
        ArrayList<Integer> vertices = mapRepresentativeToVertices.get(vertexID);
        return (vertices == null ? new ArrayList<>(Arrays.asList(vertexID)) : vertices);
    }
    
    public BigDecimal getScore() {
        return consecutiveLargeInterfaceScore;
    }
    
    /**
     * Tests clusteringResult implementation.
     */
    public static void main() {
        // test toNewick for non-binary dendrogram
        //   Result: ((((1,2),3),(4,5,6),7));
        ClusteringResult clusteringResult = new ClusteringResult(false);
        clusteringResult.addMerge(new Integer[]{0,1}, BigDecimal.ONE);
        clusteringResult.addMerge(new Integer[]{1,3}, BigDecimal.ONE);
        clusteringResult.addMerge(new Integer[]{4,5,6}, BigDecimal.ONE);
        clusteringResult.addMerge(new Integer[]{1,4,7}, BigDecimal.ONE);
        assert "((((1,2),3),(4,5,6),7));".equals(clusteringResult.toNewickString());
        
        clusteringResult = new ClusteringResult(false);
        clusteringResult.addMerge(new Integer[]{0,1}, BigDecimal.valueOf(100));
        clusteringResult.addMerge(new Integer[]{0,3}, BigDecimal.valueOf(99));
        clusteringResult.addMerge(new Integer[]{0,2}, BigDecimal.valueOf(196));
        clusteringResult.addMerge(new Integer[]{0,4}, BigDecimal.valueOf(194));
        BigDecimal expectedScore = BigDecimal.valueOf(170.8);
        System.out.println("Expected score: " + expectedScore.toString());
        System.out.println("Computed score: " + clusteringResult.getScore().toString());
        assert clusteringResult.getScore().equals(expectedScore);
        
        clusteringResult = new ClusteringResult(false);
        clusteringResult.addMerge(new Integer[]{1,3}, BigDecimal.valueOf(99));
        clusteringResult.addMerge(new Integer[]{1,2}, BigDecimal.valueOf(196));
        clusteringResult.addMerge(new Integer[]{1,4}, BigDecimal.valueOf(194));
        clusteringResult.addMerge(new Integer[]{0,1}, BigDecimal.valueOf(100));
        expectedScore = BigDecimal.valueOf(118).add(
                BigDecimal.valueOf(196).divide(BigDecimal.valueOf(3), ComplexGraphEdgeWeightTypes.PRECISION, RoundingMode.HALF_UP));
        System.out.println("Expected score: " + expectedScore.toString());
        System.out.println("Computed score: " + clusteringResult.getScore().toString());
        assert clusteringResult.getScore().equals(expectedScore);
    }
}
