/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */

package io;

import java.awt.Rectangle;
import java.io.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.image.PNGTranscoder;
import org.apache.commons.compress.archivers.ArchiveException;
import org.apache.commons.compress.archivers.ArchiveStreamFactory;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.utils.IOUtils;
import org.w3c.dom.Document;
import resources.Resources;
import tools.DP;

/**
 * This static class holds various methods related to input and output. It makes use of the Apache commons compress
 * library to extract .tar.gz. files.
 * 
 * @author ts
 */
public class IO {
    
    public static final String CLASS_TAG = "IO";
   
    public static Integer[] mapStringIntegerToArraySortedByMapKeys(Map<String, Integer> maxDiams) {
        List<String> sortedKeys = new ArrayList<>(maxDiams.keySet());
        Collections.sort(sortedKeys);
        
        Integer[] res = new Integer[maxDiams.keySet().size()];
        
        int i = 0;
        for(String key : sortedKeys) {
            res[i] = maxDiams.get(key);
            i++;
        }
        return res;
    }
    
    public static double[] mapStringDoubleToArraySortedByMapKeys(Map<String, Double> maxDiams) {
        List<String> sortedKeys = new ArrayList<>(maxDiams.keySet());
        Collections.sort(sortedKeys);
        
        double[] res = new double[maxDiams.keySet().size()];
        
        int i = 0;
        for(String key : sortedKeys) {
            res[i] = maxDiams.get(key);
            i++;
        }
        return res;
    }
    
    public static double[] integerArrayToDoubleArray(Integer[] in) {
        double[] out = new double[in.length];
        for(int i = 0; i < in.length; i++) {
            out[i] = in[i].doubleValue();
        }
        return out;
    } 
    
    /**
     * Transforms a path into a web server path. Currently, this replaces the Windows file separator '\' with the unix/web file separator '/' in the input string and returns the result as a new string.
     * Note that the original string is NOT altered.
     * @param path the input path
     * @return the modified web path version
     */
    public static String pathToWebPath(String path) {
        if(path == null) { return null; }
        String copy = path.replace("\\", "/");
        return copy; 
    }
    
    /**
     * Determines whether a string (repr. a path) ends with / or \
     * @param baseOutputDir the input path
     * @return whether it ends with / or \
     */
    public static Boolean pathEndsWithFsAlready(String baseOutputDir) {
        if(baseOutputDir.length() == 0) {
            return false;
        }
        else {
            if(baseOutputDir.endsWith("/") || baseOutputDir.endsWith("\\")) {
                return true;
            }
            else {
                return false;
            }
        }
    }
    
    /**
     * Generates PDB-style subdir tree for PDB ID and chain ('ic/8icd/A'), creates required dir if it does not exist, and returns relative path to it.
     * @param baseOutputDir the base output dir to work from, e.g., "." or "data"
     * @param pdbid the PDB ID, e.g., "8icd"
     * @param chain the chain name, e.g., "A"
     * @return the path of the subdir tree, without (=excluding) the baseOutputDir. So if you gave a baseOutputDir="data/", pdbid="8icd" and chain="A", if will return "ic/8icd/A/". If the directory could NOT be created, it will return null.
     */
    public static String createSubDirTreeDir(String baseOutputDir, String pdbid, String chain) {
        File targetDir = IO.generatePDBstyleSubdirTreeName(new File(baseOutputDir), pdbid, chain);
        String expDir = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, chain);  // something like 'ic/8icd/A'
        if(targetDir != null) {
            ArrayList<String> errors = IO.createDirIfItDoesntExist(targetDir);
            if( ! errors.isEmpty()) {
                for(String err : errors) {
                    DP.getInstance().e("IO", "dirOrDie: " + err);
                }
            } else {                
                String fsOrNot = File.separator;
                if(IO.pathEndsWithFsAlready(baseOutputDir)) { fsOrNot = ""; }
                String expDirRelToOutputDir = fsOrNot + expDir + File.separator;
                String expDirRelToWorkDir = baseOutputDir + fsOrNot + expDir + File.separator;
                File expDirFile = new File(expDirRelToWorkDir);
                if(expDirFile != null) {
                    if(expDirFile.canRead() && expDirFile.isDirectory()) {
                        return expDirRelToOutputDir;
                    }
                }
            }                    
        } else {
            DP.getInstance().e("IO", "dirOrDie: Could not determine PDB-style subdir path name.");                
        }
        return null;
    }
    
    
    /**
     * Parses a property file, assumes that all keys are integer values and that all values are strings. Stores them in a map and returns it.
     * @param mappingsFile the input file name
     * @return the result map
     * @throws IOException if IO stuff goes wrong
     */
    public static Map<Integer, String> parseMappingsFile(String mappingsFile) throws IOException {
        Map<Integer, String> m = new HashMap<>();
        
        Properties mappings = new Properties();

        BufferedInputStream stream = new BufferedInputStream(new FileInputStream(mappingsFile));
        mappings.load(stream);
        stream.close();

         for (Map.Entry<Object, Object> entry : mappings.entrySet()) {
            String key = (String)entry.getKey();
            String value = (String)entry.getValue();
            m.put(Integer.valueOf(key), value);
        }     
        
        return m;
    }
    
    /**
     * Converts an area of interest within an input SVG image to PNG format and writes it to the specified file system location.
     * @param inputFileSVG the input file, has to be in SVG format (e.g., "infile.svg")
     * @param outputFilePNG the output path for the PNG image, including file extension (e.g., "outfile.png")
     * @param aoiX x coordinate of the upper left corner of the area of interest in the SVG image
     * @param aoiY y coordinate of the upper left corner of the area of interest in the SVG image
     * @param aoiWidth width of the area of interest
     * @param aoiHeight height of the area of interest
     * @return true if it worked out, false if some error occurred
     */
    public static Boolean convertSVGtoPNG(String inputFileSVG, String outputFilePNG, Integer aoiX, Integer aoiY, Integer aoiWidth, Integer aoiHeight) {
        
        
        PNGTranscoder png = new PNGTranscoder();
        
        // define area of interest
        Rectangle aoi = new Rectangle();
        aoi.x = aoiX;
        aoi.y = aoiY;
        aoi.width = aoiWidth;
        aoi.height = aoiHeight;
        
        png.addTranscodingHint(PNGTranscoder.KEY_AOI, aoi);
        png.addTranscodingHint(PNGTranscoder.KEY_WIDTH, new Float(aoi.width) );
        png.addTranscodingHint(PNGTranscoder.KEY_HEIGHT, new Float(aoi.height) );

        BufferedReader istream;
        
        try {
            istream = new BufferedReader(new FileReader(inputFileSVG));
        } catch (FileNotFoundException ex) {
            System.err.println("ERROR: convertSVGtoPNG(): SVG input file  '" + inputFileSVG + "' not found.");
            return(false);
        }
        
        TranscoderInput input = new TranscoderInput(istream);
        FileOutputStream fout;
        TranscoderOutput fileOutput;
        
        try {
            fout = new FileOutputStream(outputFilePNG);
            fileOutput= new TranscoderOutput(fout);
            png.transcode(input, fileOutput);
        } catch (Exception ex) {
            System.err.println("ERROR: convertSVGtoPNG(): Can't write to output file '" + outputFilePNG + "'.");
            return(false);
        }
                                   
        try {
            fout.flush();
            fout.close();
        } catch (IOException ex) {
            System.err.println("ERROR: convertSVGtoPNG(): Writing of output file '" + outputFilePNG + "' failed, could not flush buffer.");
            return(false);
        }
        return(true);
    }
    
    /**
     * Checks whether the two lists a and b contain the same values (not Integer objects!) in the same order.
     * @param a the first list
     * @param b the second list
     * @return true if the two lists a and b contain the same values (not Integer objects!) in the same order, false otherwise
     */
    public static Boolean intListsContainsSameValuesInSameOrder(List<Integer> a, List<Integer> b) {
        if(a.size() != b.size()) {
            return false;
        }
        
        int u, v;
        for(int i = 0; i < a.size(); i++) {
            u = a.get(i).intValue();
            v = b.get(i).intValue();
            if( u != v) {
                return false;
            }
        }
        return true;
    }
    
    public static String listOfintegerArraysToString(List<Integer[]> list) {
        StringBuilder sb = new StringBuilder();
        for(Integer[] in : list) {
            sb.append("[");
            for(int i = 0; i < in.length; i++) {
                sb.append(in[i]);
                if(i < in.length - 1) {
                    sb.append(",");
                }
            }
            sb.append("]");
        }
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of an array.
     * @param in the array
     * @return the string
     */
    public static String integerArrayToString(Integer[] in) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < in.length; i++) {
            sb.append(in[i]);
            if(i < in.length - 1) {
                sb.append(",");
            }
        }
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of an array.
     * @param in the array
     * @return the string
     */
    public static String intArrayToString(int[] in) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < in.length; i++) {
            sb.append(in[i]);
            if(i < in.length - 1) {
                sb.append(",");
            }
        }
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of an array.
     * @param in the array
     * @return the string
     */
    public static String floatArrayToString(float[] in) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < in.length; i++) {
            //sb.append(in[i]);
            sb.append(String.format(Locale.ENGLISH, "%.2f", in[i]));
            if(i < in.length - 1) {
                sb.append(", ");
            }
        }
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of an array.
     * @param in the array
     * @return the string
     */
    public static String stringArrayToString(String[] in) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < in.length; i++) {
            sb.append(in[i]);
            if(i < in.length - 1) {
                sb.append(",");
            }
        }
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of an array.
     * @param in the array
     * @return the string
     */
    public static String doubleArrayToString(Double[] in) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < in.length; i++) {
            sb.append(in[i]);
            if(i < in.length - 1) {
                sb.append(",");
            }
        }
        return sb.toString();
    }
    
    
    /**
     * Debug function to get a string representation of a Map.
     * @param map theMap
     * @return the string
     */
    public static String mapIntegerToString(Map<Integer, Integer> map) {
        StringBuilder sb = new StringBuilder();
        for(Entry e : map.entrySet()) {
            sb.append(e.getKey()).append("=>").append(e.getValue()).append(" ");
        }
        return sb.toString();
    }
    
    
    /**
     * Creates a string of a map with one key-value pair per line.
     * @param <K> datatype of the keys
     * @param <V> datatype of the values
     * @param map LinkedHashMap so that the order is remained
     * @param prefix added before each key-value pair, e.g., "\t"
     * @param connector added between key and value, e.g. " " or "="
     * @return 
     */
    public static <K, V> String mapToKeyValueString(LinkedHashMap<K,V> map, String prefix, String connector) {
        String keyValueString = "";
        for (K key : map.keySet()) {
            keyValueString += prefix + key.toString() + connector + map.get(key) + "\n";
        }
        return keyValueString;
    }
    
    
    /**
     * Debug function to get a string representation of a Set.
     * @param set the Set
     * @return the string
     */
    public static String setIntegerToString(Set<Integer> set) {
        StringBuilder sb = new StringBuilder();
        for(Integer i : set) {
            sb.append(i).append(" ");
        }
        return sb.toString();
    }
    
    public static String hashMapValuesGreater0ToString(HashMap<Integer, Integer> map) {
        StringBuilder sb = new StringBuilder();
        for(Entry e : map.entrySet()) {
            if((Integer)e.getValue() > 0) {
                sb.append(e.getKey()).append("=>").append(e.getValue()).append(" ");
            }
        }
        return sb.toString();
    }
    
    public static String hashMapValuesGreater0OfKeysToString(HashMap<Integer, Integer> map, List<Integer> interestedInKeys) {
        StringBuilder sb = new StringBuilder();
        Integer value;
        for(Integer key : interestedInKeys) {
            value = map.get(key);
            if(value > 0) {
                sb.append(key).append("=>").append(value).append(" ");
            }
        }
        return sb.toString();
    }
    
    
    /**
     * Ugly hack to remove the prefix "./" or ".\" from a string (if it exists).
     * @param path the input string, which may or may not start with "./" or ".\"
     * @return if the input string starts with one of the mentioned patterns, returns the input string with the pattern removed. Otherwise, the output string is the input string.
     */
    public static String stripTrailingShitFromPathIfThere(String path) {
        String np = path;
        if(np.startsWith("./") || np.startsWith(".\\")) {
            np = np.substring(2);
        }
        return np;
    }
    
    public static String removeSlashAtEndIfThere(String path) {
        String np = path;
        if(np.length() == 0) {
            return np;
        }
        if(np.endsWith("/") || np.endsWith("\\")) {
            np = np.substring(0, np.length() - 1);
        }
        return np;
    }
    
    /** Untar an input file into an output file.

     * The output file is created in the output folder, having the same name
     * as the input file, minus the '.tar' extension. 
     * 
     * @param inputFile     the input .tar file
     * @param outputDir     the output directory file. 
     * @throws IOException 
     * @throws FileNotFoundException
     *  
     * @return  The {@link List} of {@link File}s with the untared content.
     * @throws ArchiveException 
     */
    public static List<File> unTar(final File inputFile, final File outputDir) throws FileNotFoundException, IOException, ArchiveException {
        
        final List<File> untaredFiles = new LinkedList<File>();
        final InputStream is = new FileInputStream(inputFile); 
        final TarArchiveInputStream debInputStream = (TarArchiveInputStream) new ArchiveStreamFactory().createArchiveInputStream("tar", is);
        TarArchiveEntry entry = null; 
        while ((entry = (TarArchiveEntry)debInputStream.getNextEntry()) != null) {
            final File outputFile = new File(outputDir, entry.getName());
            if (entry.isDirectory()) {
                if (!outputFile.exists()) {
                    if (!outputFile.mkdirs()) {
                        throw new IllegalStateException(String.format("Couldn't create directory %s.", outputFile.getAbsolutePath()));
                    }
                }
            } else {
                final OutputStream outputFileStream = new FileOutputStream(outputFile); 
                IOUtils.copy(debInputStream, outputFileStream);
                outputFileStream.close();
            }
            untaredFiles.add(outputFile);
        }
        debInputStream.close(); 

        return untaredFiles;
    }
    
    public static boolean fileExistsIsFileAndCanRead(File f) {
        if(f != null) {
            if(f.isFile() && f.canRead()) {
                return true;
            }            
        }
        return false;
    }
    
    public static boolean dirExistsIsDirectoryAndCanWrite(File f) {
        if(f != null) {
            if(f.isDirectory() && f.canWrite()) {
                return true;
            }            
        }
        return false;
    }

    /**
     * Ungzip an input file into an output file.
     * <p>
     * The output file is created in the output folder, having the same name
     * as the input file, minus the '.gz' extension. 
     * 
     * @param inputFile     the input .gz file
     * @param outputDir     the output directory file. 
     * @throws IOException 
     * @throws FileNotFoundException
     *  
     * @return  The file with the ungzipped content.
     */
    public static File unGzip(final File inputFile, final File outputDir) throws FileNotFoundException, IOException {

        final File outputFile = new File(outputDir, inputFile.getName().substring(0, inputFile.getName().length() - 3));

        final GZIPInputStream in = new GZIPInputStream(new FileInputStream(inputFile));
        final FileOutputStream out = new FileOutputStream(outputFile);

        for (int c = in.read(); c != -1; c = in.read()) {
            out.write(c);
        }

        in.close();
        out.close();

        return outputFile;
    }
    
    /**
     * Tries to delete all Files in the list.
     * @param files a list of files
     * @return the number of successfully deleted files
     */ 
    public static Integer deleteFiles(ArrayList<File> files) {
        
        Integer numSuccess = 0;
        
        for (File f : files) {
            if(f.delete()) {
                numSuccess++;
            }
        }
        
        return(numSuccess);
    }
    
    
    /**
     * Writes a SVG Doc to a PNG format file.
     * @param doc
     * @param outputFilename
     * @param aoi
     * @throws Exception 
     */
    public static void writeSVGDOC2PNG (Document doc, String outputFilename, Rectangle aoi) throws Exception {
        
        PNGTranscoder trans = new PNGTranscoder();
        
        trans.addTranscodingHint(PNGTranscoder.KEY_WIDTH, new Float(aoi.width));
        trans.addTranscodingHint(PNGTranscoder.KEY_HEIGHT, new Float(aoi.height));
        trans.addTranscodingHint(PNGTranscoder.KEY_AOI , aoi);

        TranscoderInput input = new TranscoderInput(doc);
        OutputStream ostream = new FileOutputStream(outputFilename);
        TranscoderOutput output = new TranscoderOutput(ostream);
        trans.transcode(input, output);

        ostream.flush();
        ostream.close();
    }
    
    
    /**
     * Writes the string 'text' to the text file 'targetFile'. Tries to create the file and overwrite stuff in it.
     * @return true if it worked out, false otherwise. Will spit warning to STDERR if things go wrong.
     */ 
    public static Boolean stringToTextFile(String targetFile, String text) {
        String file = targetFile;
        FileWriter fw = null;
        PrintWriter pw = null;

        try {
            fw = new FileWriter(file);
            pw = new PrintWriter(fw);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + file + "': " + e.getMessage() + ".");            
            return(false);
        }

        pw.print(text);      
        pw.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            fw.close();
        } catch(Exception ex) {
            DP.getInstance().w("Could not close FileWriter for file '" + file + "': " + ex.getMessage() + ".");
            return(false);
        }
        
        return(true);
    }
    
    
    
    public static String space(int s) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < s; i++) {
            sb.append("    ");
        }
        return sb.toString();
    }
    
    public static String treeSpace(int s) {
        if(s == 0) { return ""; }
        StringBuilder sb = new StringBuilder();        
        for(int i = 0; i < s - 1; i++) {
            sb.append("     ");
        }
        sb.append("+----");
        return sb.toString();
    }
    
    /**
     * Debug function to get a string representation of a list.
     * @param ar the list
     * @param fs the field separator (e.g., " ")
     * @return the string
     */
    public static String intListToString(List<Integer> ar, String fs) {
        StringBuilder s = new StringBuilder();
        for(int i = 0; i < ar.size(); i++) {
            s.append(ar.get(i));
            if(i < ar.size() - 1) {
                s.append(fs);
            }
        }
        return s.toString();
    }
    
    
    /**
     * Prints a map to a string.
     * @param m the Map
     * @param keyValueConnectionString the string to use to connect a key with its value (try ':' or '=' if in doubt)
     * @param pairSeparator the string to print between key-value pairs (try ' ' or ', ' if in doubt)
     * @return the string
     */
    public static String intMapToString(Map<Integer, Integer> m, String keyValueConnectionString, String pairSeparator) {
        StringBuilder s = new StringBuilder();
        for(Integer i : m.keySet()) {
            s.append(i).append(keyValueConnectionString).append(m.get(i)).append(pairSeparator);
        }
        return s.toString();
    }
    
    /**
     * Debug function to get a string representation of a list, using space as the separator.
     * @param ar the list
     * @return the string
     */
    public static String intListToString(List<Integer> ar) {        
        return IO.intListToString(ar, " ");
    }
    
    /**
     * Debug function to get a string representation of a list.
     * @param ar the list
     * @param fs the field separator (e.g., " ")
     * @return the resulting string
     */
    public static String stringListToString(List<String> ar, String fs) {
        StringBuilder s = new StringBuilder();
        for(int i = 0; i < ar.size(); i++) {
            s.append(ar.get(i));
            if(i < ar.size() - 1) {
                s.append(fs);
            }
        }
        return s.toString();
    }
    
    
    /**
     * Debug function to get a string representation of a list, using space as the separator.
     * @param ar the list
     * @return the string
     */
    public static String stringListToString(List<String> ar) {        
        return IO.stringListToString(ar, " ");
    }
    
    /**
     * Debug function to get a string representation of a list.
     * @param ar the list
     * @return the string
     */
    public static String intListToString(List<Integer> ar, String start, String end) {
        StringBuilder sb = new StringBuilder();
        sb.append(start);
        for(Integer i : ar) {
            sb.append(i).append(" ");
        }
        sb.append(end);
        return sb.toString();
    }
    
    
    /**
     * Tries to create targetDir if it does not yet exist. Uses the existing dir (without deleting anything in there) if it does
     * already exist.
     * @param targetDir the directory, may be a dir structure
     * @return a list of error messages. If this list is empty, everything worked and the dir is ready. If not, something went wrong.
     */
    public static ArrayList<String> createDirIfItDoesntExist(File targetDir) {
        
        ArrayList<String> problems = new ArrayList<String>();
        
        if(targetDir == null) {
            problems.add("Directory is null.");
            return problems;
        }
        
        if(targetDir.isDirectory()) {
            // dir already exsts
            if( ! targetDir.canWrite()) {
                problems.add("Cannot write to existing output directory '" + targetDir.getAbsolutePath() + "'.");
            } else {
                // all fine, we cqan write to it and will just use the existing one
                return new ArrayList<String>();
            }
        } else {
            if(targetDir.isFile()) {
                problems.add("Cannot create output directory '" + targetDir.getAbsolutePath() + "', file (not a directory!) with that name exists.");
            }
            else {
                try {
                    Boolean resMkdir = targetDir.mkdirs();
                    if(resMkdir) {
                        // all ok, we created it (and thus can write to it)
                        return new ArrayList<String>();                        
                    }
                }catch(Exception e) {
                    problems.add("Could not create required directory structure.");
                }
            }
        }
        
        return problems;
    }
    
    
    /**
     * Generates the PDB-style subdir path name for a given PDB ID and base directory
     * (like on the RCSB PDB FTP-server), but does NOT create the path in the file system.
     * 
     * @param pdbid the pdbid, a 4 character string
     * @param baseOutputDir the base dir under which to create the path
     * @return the path or null if no such path could be generated from the input data
     */
    public static File generatePDBstyleSubdirTreeName(String pdbid, File baseOutputDir) {
        
        if(baseOutputDir == null) {
            return null;
        }
        //System.out.println("generatePDBstyleSubdirTreeName: baseOutputDir=" + baseOutputDir.getAbsolutePath() + ".");
        
        File dirStructure;
        String fs = System.getProperty("file.separator");
        
        if(! (pdbid.length() == 4)) {
            //System.err.println("ERROR: PDB ID of length 4 required to output images in directory tree, using '" + baseOutputDir + "'.");
            dirStructure = null;
            //System.exit(1);
        } else {                    
            String mid2Chars = pdbid.substring(1, 3);                    
            dirStructure = new File(baseOutputDir.getAbsolutePath() + fs + mid2Chars + fs + pdbid);
        }
        
        return dirStructure;
    }
    
    /**
     * Generates the PDB-style subdir path name for a given PDB ID, chain and base directory
     * (similar to on the RCSB PDB FTP-server), but does NOT create the path in the file system.
     * 
     * @param pdbid the pdbid, a 4 character string
     * @param chain the pdb chain id, a 1 character string
     * @param baseOutputDir the base dir under which to create the path
     * @return the path or null if no such path could be generated from the input data
     */
    public static File generatePDBstyleSubdirTreeNameWithChain(File baseOutputDir, String pdbid, String chain) {
        
        if(baseOutputDir == null) {
            return null;
        }
        //System.out.println("generatePDBstyleSubdirTreeName: baseOutputDir=" + baseOutputDir.getAbsolutePath() + ".");
        
        File dirStructure;
        String fs = System.getProperty("file.separator");
        
        if(! (pdbid.length() == 4)) {
            //System.err.println("ERROR: PDB ID of length 4 required to output images in directory tree, using '" + baseOutputDir + "'.");
            dirStructure = null;
            //System.exit(1);
        } else {                    
            String relPathToBaseOutputDir = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, chain);
            dirStructure = new File(baseOutputDir.getAbsolutePath() + fs + relPathToBaseOutputDir);
        }
        
        return dirStructure;
    }
    
    
    /**
     * Generates the PDB-FTP server style-sub path from the pidb id and chain. This path consists of 3 sub directory
     * names: the mid 2 chars, the PDB id and the chain. Example: for 8icd, chain A, this is 'ic/8icd/A'.
     * @param pdbid the PDB id, e.g., 8icd
     * @param chain the PDB chain name, e.g., A
     * @return the PDB-style sub dir path, e.g., 'ic/8icd/A'. Note that there are no slashes at the start and end.
     */
    public static String getRelativeOutputPathtoBaseOutputDir(String pdbid, String chain) {
        String fs = System.getProperty("file.separator");
        String res;
        if(! (pdbid.length() == 4)) {
            DP.getInstance().w("IO.getRelativeOutputPathtoBaseOutputDir()", "The length of the PDB ID is not 4, returning empty path.");
            res = "";
        } else {                    
            String mid2Chars = pdbid.substring(1, 3);  
            res = mid2Chars + fs + pdbid + fs + chain;
        }
        
        return res;
    }
    
    /**
     * Generates the PDB-style subdir path name for a given PDB ID and base directory
     * (similar to on the RCSB PDB FTP-server), but does NOT create the path in the file system.
     * 
     * @param pdbid the pdbid, a 4 character string
     * @param baseOutputDir the base dir under which to create the path
     * @return the path or null if no such path could be generated from the input data
     */
    public static File generatePDBstyleSubdirTreeName(File baseOutputDir, String pdbid) {
        
        if(baseOutputDir == null) {
            return null;
        }
        //System.out.println("generatePDBstyleSubdirTreeName: baseOutputDir=" + baseOutputDir.getAbsolutePath() + ".");
        
        File dirStructure;
        String fs = System.getProperty("file.separator");
        
        if(! (pdbid.length() == 4)) {
            //System.err.println("ERROR: PDB ID of length 4 required to output images in directory tree, using '" + baseOutputDir + "'.");
            dirStructure = null;
            //System.exit(1);
        } else {                    
            String mid2Chars = pdbid.substring(1, 3);                    
            dirStructure = new File(baseOutputDir.getAbsolutePath() + fs + mid2Chars + fs + pdbid);
        }
        
        return dirStructure;
    }
    
    
   /**
    * Generates the PDB-style subdir path name for a given PDB ID, chain and base directory
    * (simlar to on the RCSB PDB FTP-server), but does NOT create the path in the file system.
    *
    * @param pdbid the pdbid, a 4 character string
    * @param chain the pdb chain id, a 1 character string
    * @param baseOutputDir the base dir under which to create the path
    * @return the path or null if no such path could be generated from the input data
    */
    public static File generatePDBstyleSubdirTreeName(File baseOutputDir, String pdbid, String chain) {
        if(baseOutputDir == null) {
            return null;
        }
        
        //System.out.println("generatePDBstyleSubdirTreeName: baseOutputDir=" + baseOutputDir.getAbsolutePath() + ".");
        File dirStructure;
        String fs = System.getProperty("file.separator");
        
        if(! (pdbid.length() == 4)) {
            //System.err.println("ERROR: PDB ID of length 4 required to output images in directory tree, using '" + baseOutputDir + "'.");
            dirStructure = null;
            //System.exit(1);
        } else {
            String mid2Chars = pdbid.substring(1, 3);
            dirStructure = new File(baseOutputDir.getAbsolutePath() + fs + mid2Chars + fs + pdbid + fs + chain);
        }
        return dirStructure;
    } 
    
    
    /**
     * Copies a file from this JAR's resources to a file system destination.
     * @param pathToResourceFile the path to the resource inside the JAR, e.g., "resources/vplg_logo.png".
     * @param targetFile the destination file
     * @return true it the file was copied, false otherwise
     * @throws Exception if the resource was not found
     */
    public static boolean copyResourceFileToFileSystemLocation(String pathToResourceFile, File targetFile) throws Exception {
        boolean found = false;
        InputStream is = Resources.class.getClassLoader().getResourceAsStream(pathToResourceFile);
        if(is != null) {
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            FileOutputStream fos = new FileOutputStream(targetFile);
            byte[] buffer = new byte[1024];
            int bytesRead;
            
            while ((bytesRead = is.read(buffer)) != -1) {
                fos.write(buffer, 0, bytesRead);
            }
            fos.flush();
            br.close();
            is.close();
            found = true;
        } else {
            throw new Exception ("Resource '" + pathToResourceFile + "' not found.");
        }
        return found;
    }
    
    
    /**
     * Converts a matrix to a string in Parek NET format. See http://gephi.github.io/users/supported-graph-formats/pajek-net-format/. Note though that vertices start with 1 (not 0).
     * @param matrix the input matrix (a graph). The 2 dimensions must be equal (length of outer and inner matrix). It must not be null.
     * @return a string in Parek NET format representing the graph
     */
    public static String intMatrixToPajekFormat(Integer[][] matrix) {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("*Vertices ").append(matrix.length).append("\n");
        
        // We do not print the vertex list since we do not have vertex labels in the int[][] matrix and the number of verts in enough to know.         
        //for(int i = 0; i < matrix.length; i++) {
        //    sb.append("").append(i+1).append("\n");
        //}
        
        // So go with the edges:        
        for(int i = 0; i < matrix.length; i++) {
            for(int j = i+1; j < matrix.length; j++) {
                if(matrix[i][j] > 0) {
                    sb.append("").append(i+1).append(" ").append(j+1).append("\n");
                }
            }
        }
        
        return sb.toString();
    }

    /**
     * Writes a string to a file.
     * @param content of the file
     * @param path of the file
     * @param overwrite whether to overwrite if file already exists
     * @return whether it worked
     */
    public static Boolean writeStringToFile(String content, String path, Boolean overwrite) {

        // check if file exists and we should not overwrite
        File tmpFile = new File(path);
        if (tmpFile.isFile() && !overwrite) {
            DP.getInstance().w(CLASS_TAG, "Tried creating file '" + path + "', but file already existed and "
                    + "overwrite flag was set to false. Returning failure and going on.");
            return false;
        }
        
        // try to write file
        try {
            FileWriter tmpWriter = new FileWriter(path);
            tmpWriter.write(content);
            tmpWriter.flush();
            tmpWriter.close();
        } catch (IOException e) {
            DP.getInstance().e(CLASS_TAG, "Could not create file '" + path + "'. Error: " + e.toString() + ". "
                    + "Returning failure and going on.");
            return false;
        }
        
        return true;  // if code reaches here, it worked
    }
    
    
    public static Properties readPropertiesFromFile(String path) {
        Properties tempProperties = new Properties();
        try {
            BufferedInputStream stream = new BufferedInputStream(new FileInputStream(path));
            tempProperties.load(stream);
            stream.close();
            return tempProperties;
        } catch (IOException e) {
            DP.getInstance().w(CLASS_TAG, "Could not read properties from file '" + path + "'. Error: " + e.toString() + ". " +
                    "Returning null and going on." );
            return null;
        }
    }
    
    
    public static String readFileToString(String path) {
        try {
            return Files.readString(Path.of(path));
        } catch (IOException e) {
            DP.getInstance().w(CLASS_TAG, "Could not read file '" + path + "' to String. Error: " + e.toString() + ". "
                    + "Returning empty string and going on.");
        }
        return "";
    }
    
    
    /**
     * Expects a path with filename and extends the name if already existing.
     * @param absPath absolute filepath
     * @return 
     */
    public static String getUniqueFilename(String absPath) {
        File uniqueFile = new File(absPath);
        Integer extension = 0;  // count which is applied to file name to create an unique one
        while (uniqueFile.exists()) {
            extension ++;
            uniqueFile = new File(absPath + "_" + extension.toString());
        }
        return uniqueFile.getAbsolutePath();
    }
}
