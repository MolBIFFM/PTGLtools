/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import tools.DP;

/**
 *
 * @author jnw
 */
public class ClusteringResult {
    private final static String CLASS_TAG = "CG";
    
    private final ArrayList<Integer[]> merges = new ArrayList<>();
    private final Boolean binary;
    private int maxNumberMergedAtOnce = 0;
    
    public ClusteringResult(Boolean binary) {
        this.binary = binary;
    }
    
    public Boolean addMerge(Integer[] merge) {
        if (binary) {
            if (merge.length > 2) {
                DP.getInstance().e(CLASS_TAG, "Tried adding a merge of more than two vertices to a binary clustering result. Rejecting merge and trying to go on. "
                        + "Please report this bug to the developers.");
            }
        }
        merges.add(merge);
        maxNumberMergedAtOnce = Math.max(maxNumberMergedAtOnce, merge.length);
        return true;
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
        return listOfArrayToString(renamedMerges);
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
}
