/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import settings.Settings;
import tools.DP;

/**
 *
 * @author Marvin
 */
public class SseParser {
    
    private static boolean silent = Settings.getBoolean("PTGLgraphComputation_B_silent");
    
    private static String sseFile;
    
    private static HashMap<String,Integer> colHeaderPosMap = new HashMap<>();
    private static String[] lineData;
    private static boolean inStructConf;
    
    protected static void setSseFile(String sf){
        sseFile = sf;
    }
    
    /**
     * largely copied from CifParser. Reads file and creates SseRanges from _struct_conf data-lines.
     */
    protected static void createSseRangesFromDsspCifFile(){
        String columnHeader;
        String line;
        
        try{
            BufferedReader in = new BufferedReader(new FileReader(sseFile));
            while ((line = in.readLine()) != null){
                System.out.println("============================================================================\n" + line);
                if (line.startsWith("#")){
                    inStructConf = false;  // '#' closes a category, so we can't possibly be in the struct_conf category
                    lineData = null;  // clear lineData. so if it is accessed after the struct_conf category was parsed, it doesn't contain any misleading data
                }
                else if (line.startsWith("_struct_conf.")){  // check for '_struct_conf.' (this category contains the Sse-Info) and save the column header (the value after the dot
                    inStructConf = true;
                    columnHeader = line.split("\\.")[1].trim();
                    colHeaderPosMap.put(columnHeader, colHeaderPosMap.size());
                    System.out.println("added " + columnHeader + " to colHeaderPosMap");
                    System.out.println(colHeaderPosMap);
                }
                else if (inStructConf && ! line.startsWith("_")){  // if we are in the struct_conf category and the line does not start with a '_', is must be a line containing data --> create SSE from that line
                    lineData = CifParser.lineToArray(line);
                    System.out.println("parsed the line to: " + Arrays.toString(lineData) + " now calling handle");
                    CifParser.handleStructConfLine(lineData, colHeaderPosMap);
                }
            }
        } catch (IOException e) {
            System.err.println("ERROR: Could not parse SSE file.");
            System.err.println("ERROR: Message: " + e.getMessage());
            System.exit(1);
        }
    }
}
