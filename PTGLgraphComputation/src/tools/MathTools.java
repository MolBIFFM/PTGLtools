/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2015. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package tools;

import java.util.Arrays;
import java.util.List;

/**
 * Some math tools for PTGLgraphComputation.
 * @author spirit
 */
public class MathTools {
    
    
    /**
     * Floor division of integers, inspired by Java8.
     * @param divident the divident
     * @param divisor the divisor, must not be null
     * @return the quotient
     */
    public static Integer floorDiv(Integer divident, Integer divisor) {
        Double divRes = (divident.doubleValue()) / divisor.doubleValue();
        Integer quotient = ((Double)Math.floor(divRes)).intValue();
        return quotient;
    }
    
    
    /**
     * Sums up any number of integer arrays element wise.
     * @param arrays any number of Integer arrays
     * @return 
     */
    public static Integer[] elementWiseSum(Integer[]... arrays) {
        // check for trivial cases
        if (arrays.length == 0) {
            return null;
        } else if (arrays.length == 1) {
            return arrays[0];
        }  // from now on: arrays length > 1
        
        // find the max and min length
        Integer maxLength = arrays[0].length;
        Integer minLength = arrays[0].length;
        for (Integer[] arr : arrays) {
            maxLength = Math.max(maxLength, arr.length);
            minLength = Math.min(minLength, arr.length);
        }
        
        Integer[] sumArr = new Integer[maxLength];
        
        // now sum the values of existing
        for (Integer[] arr : arrays) {
            for (int i = 0; i < arr.length; i++) {
                sumArr[i] = (sumArr[i] == null ? arr[i] : sumArr[i] + arr[i]);  // assign if nothing set yet otherwise sum
            }
        }
        
        return sumArr;
    }
    
    /**
     * Sums up any number of Double arrays element wise.
     * @param arrays any number of Integer arrays
     * @return 
     */
    public static Double[] elementWiseSumDouble(Double[]... arrays) {
        // check for trivial cases
        if (arrays.length == 0) {
            return null;
        } else if (arrays.length == 1) {
            return arrays[0];
        }  // from now on: arrays length > 1
        
        // find the max and min length
        Integer maxLength = arrays[0].length;
        Integer minLength = arrays[0].length;
        for (Double[] arr : arrays) {
            maxLength = Math.max(maxLength, arr.length);
            minLength = Math.min(minLength, arr.length);
        }
        
        Double[] sumArr = new Double[maxLength];
        
        // now sum the values of existing
        for (Double[] arr : arrays) {
            for (int i = 0; i < arr.length; i++) {
                sumArr[i] = (sumArr[i] == null ? arr[i] : sumArr[i] + arr[i]);  // assign if nothing set yet otherwise sum
            }
        }
        
        return sumArr;
    }
    
    /**
     * Divides each element of an array by the same divisor.
     * @param arr Integer array
     * @param divisor
     * @param round
     * @return 
     */
    public static float[] elementWiseDivision(Integer[] arr, float divisor, Boolean round) {
        // check for trivial cases
        if (divisor == 0) {
            DP.getInstance().w("Tried to divide by zero in tools.MathTools.elementWiseDivision. Returning null.");
            return null;
        }
        
        float[] quotientArr = new float[arr.length];
        
        for (int i = 0; i < arr.length; i++) {
            quotientArr[i] = arr[i] / divisor;
            if (round) { quotientArr[i] = Math.round(quotientArr[i]); }
        }

        return quotientArr;
    }
}
