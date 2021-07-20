/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package proteingraphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Contains information of the types of edge weights and normalizations of Complex Graphs.
 * @author jnw
 */
public final class ComplexGraphEdgeWeightTypes {
    
    // List here the weights / normalizations and assign short but descriptive spaced names
    public static final String EdgeWeightTypeNames[] = new String[] {
        "absolute weight",
        "additive length normalization",
        "multiplicative length normalization",
        "lucid multiplicative length normalization"
    };
    
    // List here the weights / normalizations withr reference to the array EdgeWeightTypeNames if you want an own constant
    public static final String ABSOLUTE_WEIGHT = EdgeWeightTypeNames[0];
    public static final String ADDITIVE_NORMALIZATION = EdgeWeightTypeNames[1];
    public static final String MULTIPLICATIVE_NORMALIZATION = EdgeWeightTypeNames[2];
    public static final String LUCID_MULTIPLICATIVE_NORMALIZATION = EdgeWeightTypeNames[3];
    
    // Assign here a description to each weight / normalization
    public static final HashMap<String, String> MAP_TYPE_DESCRIPTION = new HashMap<>();
    static {
        MAP_TYPE_DESCRIPTION.put(ABSOLUTE_WEIGHT, "number of residue-residue contacts");
        MAP_TYPE_DESCRIPTION.put(ADDITIVE_NORMALIZATION, ABSOLUTE_WEIGHT + " / (length chain 1 + chain 2)");
        MAP_TYPE_DESCRIPTION.put(MULTIPLICATIVE_NORMALIZATION, ABSOLUTE_WEIGHT + " / (length chain 1 * chain 2)");
        MAP_TYPE_DESCRIPTION.put(LUCID_MULTIPLICATIVE_NORMALIZATION, MULTIPLICATIVE_NORMALIZATION + " / smallest multiplicative weight");
    }
}
