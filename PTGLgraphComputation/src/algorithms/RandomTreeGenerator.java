/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package algorithms;

import java.util.ArrayList;

/**
 *
 * @author jnw
 */
public class RandomTreeGenerator {
    
    class TreeVertex {
        int ID;
        Boolean leaf;
        
        TreeVertex(int ID, Boolean leaf) {
            this.ID = ID;
            this.leaf = leaf;
        }
    }
    
    class TreeEdge {
        TreeVertex v1;
        TreeVertex v2;
        
        TreeEdge (TreeVertex v1, TreeVertex v2) {
            this.v1 = v1;
            this.v2 = v2;
        }
    }
    
    public String randomRootedBinaryTree(int numberVertices) {
        // consider trivial cases
        if (numberVertices < 1) { return ""; }
        switch (numberVertices) {
            case 1:
                return "(1);";
            case 2:
                return "(1,2);";
        }

        // usual cases: more than two vertices
        //  the following algorithm is based on Furnas, 1984, Journal of Classification
        ArrayList<TreeEdge> edges = new ArrayList<>();
        
        // 1) "doublet" tree connected by a single edge
        edges.add(new TreeEdge(new TreeVertex(0, true), new TreeVertex(1, true)));
        
        return "Return value (WIP)";
    }
}
