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
                    handleStructConfLine();
                }
            }
        } catch (IOException e) {
            System.err.println("ERROR: Could not parse SSE file.");
            System.err.println("ERROR: Message: " + e.getMessage());
            System.exit(1);
        }
    }
    
    /**
     * Handles lines starting with _struct_conf.
     * Creates a String[] representing the sse and adds that to chainwiseSseDict of the CifParser.
     */
    private static void handleStructConfLine(){
        // md
        // zum testen: die 7tim.cif file von PDB manuell umgeschrieben: dssp _struct_conf am Ende angehängt, ursprüngliche _struct_conf in struct_conf2 umbenannt
        // man muss eig. noch sicherstellen, dass es auch eine von dssp erstellte Kategorie ist. Denn _struct_conf kann auch so existieren, enthält dann aber ggf andere Daten (z.B. nur Helix-Daten in der 7tim.cif file)
        // todo: es kann passieren, dass die auth-cols auf PDB-cols gemapped werden. Dann sind wir hier aber inkonsistent, denn ich benutze immer auth cols
        // in dssp V4 schauen, welche conf_types es alles gibt und dann die HashMap bzw Funktion anpasen
        // wird die dssp-Reduction (alle Helices als H, manche Sachen zu E zusammenfassen) später schon von einer Funktion gemacht?
        // eine SSE kann nur in einer Chain sein richtig? dann muss ich nicht endChain übergeben
        // Das matchen während den ATOM-lines machen oder später mit s_molecules? Dann aber darauf achten, dass auch rna und lig in s_molecules sein können
        
        // extract important data from the _struct_conf-line (which SSE as well as it's start and end position)
        String SseStructure = lineData[colHeaderPosMap.get("conf_type_id")];
        String sseChain = lineData[colHeaderPosMap.get("beg_auth_asym_id")];
        String startIndexStr = lineData[colHeaderPosMap.get("beg_auth_seq_id")];
        String startICode = lineData[colHeaderPosMap.get("pdbx_beg_PDB_ins_code")];
        String endIndexStr = lineData[colHeaderPosMap.get("end_auth_seq_id")];
        String endICode = lineData[colHeaderPosMap.get("pdbx_end_PDB_ins_code")];
        
        // if there is no iCode present, set the variable to " " to be consistent and compatible with the existing code
        if (startICode.equals("?") || startICode.equals(".")){startICode = " ";}
        if (endICode.equals("?") || endICode.equals(".")){endICode = " ";}
        
        // add the current SSE to the ArrayList of the other SSE-Ranges
        String[] currRange = {SseStructure, sseChain, startIndexStr, startICode, endIndexStr, endICode};
        CifParser.addSseToSseDict(currRange);
    }
}
