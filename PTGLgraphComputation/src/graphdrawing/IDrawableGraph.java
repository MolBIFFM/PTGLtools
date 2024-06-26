/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2015. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package graphdrawing;

import datastructures.SparseGraph;
import java.util.List;

/**
 * A graph that can be drawn using the protein graph drawing method of the PTGL. The graph
 * needs to be edge typed and vertex typed, and the types must match the PTGL labels.
 * @author spirit
 */
public interface IDrawableGraph {
    
    public String getPropertyString(String name);
    public List<IDrawableEdge> getDrawableEdges();
    public List<IDrawableVertex> getDrawableVertices();
    public Boolean containsEdge(Integer i, Integer j);
    public String getSpatRelOfEdge(Integer i, Integer j);
    public String getFGNotationOfVertex(Integer i);
    public SparseGraph<String, String> toSparseGraph();
    
}
