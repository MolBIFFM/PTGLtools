/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package algorithms;

import java.util.ArrayList;
import java.util.HashMap;
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
        
        public TreeVertex[] getVertices () {
            return new TreeVertex[] {v1, v2};
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
    
    /**
     * Generates a random binary tree from a uniform distribution. Algorithm based on Furnas, 1984, Journal of Classification.
     * @param rooted Whether tree should be rooted
     * @param numberLeaves How many leafs the tree should have
     * @return Tree as Newick string
     */
    public String generateRandomBinaryTree(Boolean rooted, int numberLeaves) {
        // consider trivial cases
        if (numberLeaves < 1) { return ""; }
        switch (numberLeaves) {
            case 1:
                return "(1);";
            case 2:
                return "(1,2);";
        }

        // usual cases: more than two vertices
        //  the following algorithm is based on Furnas, 1984, Journal of Classification
        
        // 1) "doublet" tree connected by a single edge
        Tree tree = new Tree(true, numberLeaves * 2 - 3);  // initial capacity is 2n-3 since unrooted binary trees with n > 3 have that many edges
        
        for (int step = 2; step < numberLeaves; step++) {
            // 2) choose a random edge
            int randomEdgeID = (int) ((Math.random() * (tree.edges.size() - 1 )));
            addNewVertexToBinaryTree(tree, randomEdgeID);
        }
        
        if (Settings.getInteger("PTGLgraphComputation_I_debug_level") >= 2) {
                System.out.println("[DEBUG LV 2] List of edges: ");
                for (TreeEdge edge : tree.edges) {
                    System.out.println("  " + edge.toString());
                }
            }
        
        return treeToNewick(tree, rooted);
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
    
    
    private String treeToNewick(Tree tree, Boolean rootTree) {
        // approach 2: build up a edge list as map and traverse from root
        HashMap<TreeVertex, TreeVertex[]> adjacencyMap = new HashMap<>();
        
        if (rootTree) {
            int rndEdgeID = (int) ((Math.random() * (tree.edges.size() - 1 )));
            TreeEdge edgeAsRoot = tree.edges.get(rndEdgeID);
            
            // TODELETE
            System.out.println("Root at: " + edgeAsRoot.toString());
            
            for (TreeEdge edge : tree.edges) {
                if (edge.equals(edgeAsRoot)) { continue; }  // skip root edge
                
                // TODELETE
                System.out.println("Adding to map: " + edge.toString());
                
                addEdgeToAdjMap(adjacencyMap, edge);               
            }
            
            // TODELETE
            System.out.println("final map");
            for (TreeVertex v : adjacencyMap.keySet()) {
                System.out.println(v.ID + ": ");
                for (TreeVertex child : adjacencyMap.get(v)) {
                    if (child != null) {
                        System.out.print(child.toString() + " | ");
                    } else {
                        System.out.print("null" + " | ");
                    }
                }
                System.out.println("");
            }
            
            return "(" + neighborsToSubNewick(edgeAsRoot.v1, adjacencyMap, edgeAsRoot.v1) + "," + neighborsToSubNewick(edgeAsRoot.v2, adjacencyMap, edgeAsRoot.v2) + ");";
        } else {
            DP.getInstance().e(CLASS_TAG, "Returning of Newick string for random unrooted tree not implemented, yet. Returning empty string.");
            return "";
        }
    }
    
    private void addEdgeToAdjMap(HashMap<TreeVertex, TreeVertex[]> map, TreeEdge edge) {
        addDirectedEdgeToAdjMap(map, edge.v1, edge.v2);
        addDirectedEdgeToAdjMap(map, edge.v2, edge.v1);
    }
    
    private void addDirectedEdgeToAdjMap(HashMap<TreeVertex, TreeVertex[]> map, TreeVertex startVertex, TreeVertex endVertex) {
        // TODELETE
        System.out.println("startVertex: " + startVertex.ID);
        System.out.println("endVertex: " + endVertex.ID);
        
        if (map.containsKey(startVertex)) {
            if (map.get(startVertex)[1] == null) {
                map.get(startVertex)[1] = endVertex;  // we do not need to test for 0th element, it gets set when key is added
            } else {
                map.get(startVertex)[2] = endVertex;  // can only be 2nd element then
            }
        } else {
            map.put(startVertex, new TreeVertex[]{endVertex, null, null});
        }
        
        // TODELETE
            System.out.println("Map after addDirectedEdgeToAdjMap");
            for (TreeVertex v : map.keySet()) {
                System.out.println(v.ID + ": ");
                for (TreeVertex child : map.get(v)) {
                    if (child != null) {
                        System.out.print(child.toString() + " | ");
                    } else {
                        System.out.print("null" + " | ");
                    }
                }
                System.out.println("");
            }
    }
    
    private String neighborsToSubNewick(TreeVertex v, HashMap<TreeVertex, TreeVertex[]> adjMap, TreeVertex parent) {
        TreeVertex[] neighbors = adjMap.get(v);
        if (v.leaf) {
            // v is leaf and we end here
            return v.leafID.toString();
        } else {
            TreeVertex[] neighborsWithoutParent = removeParentFromNeighbors(neighbors, parent);
            return "(" + neighborsToSubNewick(neighborsWithoutParent[0], adjMap, v) + "," + neighborsToSubNewick(neighborsWithoutParent[1], adjMap, v) + ")";
        }
    }
    
    
    private TreeVertex[] removeParentFromNeighbors(TreeVertex[] neighbors, TreeVertex parent) {
        TreeVertex[] neighborsWithoutParent = new TreeVertex[neighbors.length - 1];
        int arrayCounter = 0;
        for (TreeVertex v : neighbors) {
            if (!( v == null) && ! (v.equals(parent))) {
                neighborsWithoutParent[arrayCounter] = v;
                arrayCounter++;
            }
        }
        return neighborsWithoutParent;
    }

    
    private String edgeToNewick(TreeEdge edge) {
        return "(" + edge.v1.ID + "," + edge.v2.ID + ")";
    }
}
