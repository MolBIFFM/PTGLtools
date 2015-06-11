/*
 * This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.
 *
 * Copyright Tim Schäfer 2015. VPLG is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package linnottograph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import plcc.SSEGraph;
import tools.DP;

/**
 * Parser for ADJ linnot strings. Note that this parses a linnot string,
 * but constructs a PG (not an FG) from it. You can use the PG function to compute FGs in order to get
 * an FG from it (there may be multiple FGs, but the other ones will be isolated vertices. If the linnot one is an isolated vertex, there will obviously be only one FG in the list.).
 * @author spirit
 */
public class LinnotParserADJ extends LinnotParserRED implements ILinnotParser {
    
    public LinnotParserADJ(String linnot, String graphType) {
        super(linnot, graphType);
    }
    
    @Override
    public List<String> getVertexTypesNtoC() {
        List<String> vtypes = new ArrayList<>();
        Integer maxShift = this.getMaxShiftLeft();
        List<Integer> visited = this.getVisitPath();
        int numVerts = this.getAllVisitedVertices().size();
        Map<Integer, String> m = new HashMap<>();
        
        //System.out.println("visited " + visited.size() + ". maxShift = " + maxShift + "." );
        
        Integer indexNtoC;
        Integer reltoStartPos;
        for(int i = 0; i < visited.size(); i++) {
            reltoStartPos = visited.get(i);
            String vtype = LinnotParserTools.getSSETypeFromToken(tokens[i], SSEGraph.GRAPHTYPE_ALBELIG);
            indexNtoC = reltoStartPos - maxShift;
            //System.out.println("indexNtoC=" + indexNtoC);
            m.put(indexNtoC, vtype);
        }
        
        for(int i = 0; i < numVerts; i++) {
            vtypes.add(m.get(i));
        }                
        
        
        
        return vtypes;
    }
    
    @Override
    public List<Integer[]> getOutGraphEdges() {
        List<Integer[]> shiftedEdges = this.getNonZEdges();
        List<Integer[]> finalEdges = new ArrayList<>();
        Integer maxShift = this.getMaxShiftLeft();
                       
        Integer [] oe;
        for(Integer[] e : shiftedEdges) {
            oe = new Integer[]{ e[0] - maxShift, e[1] - maxShift, e[2] };
            finalEdges.add(oe);
        }
        
        if(! this.distancesMakeSense()) {
            DP.getInstance().w("LinnotParserADJ", "getOutGraphEdges: Distances make no sense, linnot may not be a valid ADJ linnot string.");
        }
        return finalEdges;
    }
    
    
    
}
