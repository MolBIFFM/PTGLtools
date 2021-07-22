/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package proteingraphs;

import java.math.BigDecimal;
import java.util.HashMap;

/**
 * Contains information of the types of edge weights and normalizations of Complex Graphs.
 * @author jnw
 */
public final class ComplexGraphEdgeWeightTypes {
    
    // List here the weights / normalizations and assign short but descriptive spaced names
    public static enum EdgeWeightTypes {
        ABSOLUTE_WEIGHT,
        ADDITIVE_LENGTH_NORMALIZATION,
        MULTIPLICATIVE_LENGTH_NORMALIZATION,
        LUCID_MULTIPLICATIVE_NORMALIZATION;
        
        public String name;
        static {
            ABSOLUTE_WEIGHT.name = "absolute weight";
            ADDITIVE_LENGTH_NORMALIZATION.name = "additive length normalization";
            MULTIPLICATIVE_LENGTH_NORMALIZATION.name = "multiplicative length normalization";
            LUCID_MULTIPLICATIVE_NORMALIZATION.name = "lucid multiplicative length normalization";
        }
        
        public String description;
        static {
            ABSOLUTE_WEIGHT.description = "number of residue-residue contacts";
            ADDITIVE_LENGTH_NORMALIZATION.description = ABSOLUTE_WEIGHT.name + " / (length chain 1 + chain 2)";
            MULTIPLICATIVE_LENGTH_NORMALIZATION.description = ABSOLUTE_WEIGHT.name + " / (length chain 1 * chain 2)";
            LUCID_MULTIPLICATIVE_NORMALIZATION.description = MULTIPLICATIVE_LENGTH_NORMALIZATION.name + " / smallest multiplicative weight";
        }
    }
}
