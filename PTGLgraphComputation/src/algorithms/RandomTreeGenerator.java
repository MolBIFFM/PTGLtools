/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package algorithms;

import java.util.ArrayList;
import tools.DP;
import settings.Settings;

/**
 *
 * @author jnw
 */
public class RandomTreeGenerator {
    private final static String CLASS_TAG = "rndTree"; 
    
    class TreeVertex {
        int ID;
        Boolean leaf;
        Integer leafID;
        
        TreeVertex(int ID, Boolean leaf, int leafID) {
            this.ID = ID;
            this.leaf = leaf;
            this.leafID = leafID;
        }
        
        TreeVertex(int ID) {
            this.ID = ID;
            leaf = false;
            leafID = null;
        }
        
        @Override
        public String toString() {
            return "ID: " + ID + ", isLeaf: " + leaf + ", leafID: " + leafID;
        }
    }
    
    class TreeEdge {
        TreeVertex v1;
        TreeVertex v2;
        
        TreeEdge (TreeVertex v1, TreeVertex v2) {
            this.v1 = v1;
            this.v2 = v2;
        }
        
        @Override
        public String toString() {
            return v1.toString() + " -- " + v2.toString();
        }
    }
    
    class Tree {
        ArrayList<TreeEdge> edges;
        int maxID = 0;
        int maxLeafID = 1;
        
        Tree(Boolean doubletTree, int numExpectedEdges) {
            edges = new ArrayList<>(numExpectedEdges);
            if (doubletTree) {
                edges.add(new TreeEdge(new TreeVertex(0, true, 0), new TreeVertex(1, true, 1)));
                maxID = 1;
                maxLeafID = 1;
            }
        }
    }
    
    public String generateRandomRootedBinaryTree(int numberVertices) {
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
        
        // 1) "doublet" tree connected by a single edge
        Tree tree = new Tree(true, numberVertices * 2 - 1);  // initial capacity is 2n-1 since unrooted binary trees with n > 3 have that many edges
        
        for (int step = 2; step < numberVertices + 1; step++) {
            if (Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                System.out.println("[DEBUG LV 2] Step: " + step);
                System.out.println("[DEBUG LV 2] List of edges: ");
                for (TreeEdge edge : tree.edges) {
                    System.out.println("  " + edge.toString());
                }
            }
            
            // 2) choose a random edge
            int randomEdgeID = (int) ((Math.random() * (tree.edges.size() - 1 )));
            addNewVertexToBinaryTree(tree, randomEdgeID);
        }
        
        
        return "Return value (WIP)";
    }
    
    
    private void addNewVertexToBinaryTree(Tree tree, int edgeNumber) {
        ArrayList<TreeEdge> edges = tree.edges;

        if (edgeNumber < 0 || edgeNumber >= edges.size()) {
            DP.getInstance().w(CLASS_TAG, "Cannot add vertex at an edge position '" + edgeNumber  + "' for edge list of size '" + edges.size() + "'. "
                        + "Ignoring this step, but this might render the random tree useless. Please report this bug to the developers.");
            return;
        }
                
        // (v1,v2) -> (v1,newInnerVertex), (v2,newInnerVertex), (newLeaf,newInnerVertex)
        TreeVertex newInnerVertex = new TreeVertex(tree.maxID + 1);
        TreeVertex newLeaf = new TreeVertex(tree.maxID + 2, true, tree.maxLeafID + 1);
        TreeEdge newEdge1 = new TreeEdge(edges.get(edgeNumber).v2, newInnerVertex);  // (v2,newInnerVertex)
        TreeEdge newEdge2 = new TreeEdge(newLeaf, newInnerVertex);  // (newLeaf,newInnerVertex)
        
        edges.add(newEdge1);
        edges.add(newEdge2);
        edges.get(edgeNumber).v2 = newInnerVertex;  // (v1,v2) -> (v1,newInnerVertex)
        
        tree.maxID += 2;
        tree.maxLeafID += 1;
    }
}
