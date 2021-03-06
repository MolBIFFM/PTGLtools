/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2015. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package linnottograph;

import graphdrawing.IDrawableEdge;
import graphdrawing.IDrawableGraph;
import graphdrawing.IDrawableVertex;
import java.util.List;

/**
 *
 * @author spirit
 */
public interface ILinnotToGraph {
    public List<IDrawableEdge> getEdges();
    public List<IDrawableVertex> getVertices();
    public IDrawableGraph getDrawableGraph();
    
}
