/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2021. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */
package proteingraphs;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.HashMap;
import java.util.Map;

/**
 * Contains information of the types of edge weights and normalizations of Complex Graphs.
 * @author jnw
 */
public final class ComplexGraphEdgeWeightTypes {
    
    private final static int PRECISION = 35;  // used as precision for the BigDecimal normalized edge weight, i.e., number of digits left and right of decimal point
    
    // List here the weights / normalizations and assign short but descriptive spaced names
    public static enum EdgeWeightType {
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
    
    public static BigDecimal computeAdditiveLengthNormlization(BigDecimal numResContacts, int lengthChainA, int lengthChainB) {
        return numResContacts.divide(BigDecimal.valueOf(lengthChainA).add(BigDecimal.valueOf(lengthChainB)), PRECISION, RoundingMode.HALF_UP);
    }
    
    public static BigDecimal computeMultiplicativeLengthNormlization(BigDecimal numResContacts, int lengthChainA, int lengthChainB) {
        return numResContacts.divide(BigDecimal.valueOf(lengthChainA).multiply(BigDecimal.valueOf(lengthChainB)), PRECISION, RoundingMode.HALF_UP);
    }
    
    public static BigDecimal computeLucidMultiplicativeLengthNormlization(BigDecimal multNorm, BigDecimal smallestMultNorm) {
        return multNorm.divide(smallestMultNorm, PRECISION, RoundingMode.HALF_UP);
    }
    
    public static BigDecimal computeLengthNormalization(int numResContacts, int lengthChainA, int lengthChainB, EdgeWeightType weightType) {
        BigDecimal normalizedWeight;
        switch (weightType) {
                case ADDITIVE_LENGTH_NORMALIZATION:
                    normalizedWeight = computeAdditiveLengthNormlization(BigDecimal.valueOf(numResContacts), lengthChainA, lengthChainB);
                    break;
                case MULTIPLICATIVE_LENGTH_NORMALIZATION:
                    normalizedWeight = computeMultiplicativeLengthNormlization(BigDecimal.valueOf(numResContacts), lengthChainA, lengthChainB);
                    break;
                case ABSOLUTE_WEIGHT:
                default:
                    normalizedWeight = BigDecimal.valueOf(numResContacts);
            }
        return normalizedWeight;
    }
}
