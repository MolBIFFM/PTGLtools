/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package datastructures;

import java.util.ArrayList;
import tools.DP;

/**
 *
 * @author jnw
 */
public class ClusteringResult {
    private final static String CLASS_TAG = "CG";
    
    private ArrayList<int[]> merges = new ArrayList<>();
    private final Boolean binary;
    
    public ClusteringResult(Boolean binary) {
        this.binary = binary;
    }
    
    public Boolean addMerge(int[] merge) {
        if (binary) {
            if (merge.length > 2) {
                DP.getInstance().e(CLASS_TAG, "Tried adding a merge of more than two vertices to a binary clustering result. Rejecting merge and trying to go on. "
                        + "Please report this bug to the developers.");
            }
        }
        merges.add(merge);
        return true;
    }
}
