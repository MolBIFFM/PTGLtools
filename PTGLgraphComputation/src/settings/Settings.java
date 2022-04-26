/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2020. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */

package settings;

import graphdrawing.DrawTools.IMAGEFORMAT;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Properties;
import tools.DP;
import io.IO;

/**
 *
 * @author jnw
 */
public class Settings {
    
    // PTGLgraphComputation's version as MAJOR.MINOR.PATCH
    //   major: big (re-)implementation, new user interface, new overall architecture
    //   minor: new functions / features, git merges of branches
    //   patch: (hot-)fixes, small changes
    //   no version change: fix typos, changes to comments, debug prints, small changes to non-result output, changes within git branch
    // -> only increment with commit / push / merge not while programming
    // -> when incrementing, reset lower levels to zero, e.g., 3.2.1 -> 3.3.0
    static final private String VERSION = "3.5.0";
    
    static final private String PROGRAM_NAME = "PTGLgraphComputation";
    
    static ArrayList<Section> sections;  // holds the setting sections in their correct order
    // WARNING: Following data structures MUST be filled after sections have been filled
    static HashMap<String, String> mapSettingNameToSectionName;  // Maps settings names (=keys) to the section names they are in
    static HashMap<String, Integer> mapSectionNameToSectionIndex;  // Maps the section names (=key) to their index in the ArrayList
    
    static int numTotalSettings;  // number of all existing settings
    
    static Boolean createdDefaultFile = false;
    static int numLoadedSettings = 0;  // number of settings loaded from the default file
    
    static final String PACKAGE_TAG = "SETTINGS";
    static final String FS = System.getProperty("file.separator");
    static private final File DEFAULT_FILE = new File(System.getProperty("user.home") + FS  + "PTGLgraphComputation_settings.txt");
    
    // init all the static fields
    static {
        // first, init sections
        initSections();
        
        // WARNING: Following functions filling the data structures MUST be called after creation of sections
        // second, init names map
        initNamesMap();
        
        // third, init index map
        initIndexMap();
        
        numTotalSettings = mapSettingNameToSectionName.size();
    }

    /**
     * Creates all the settings sections.
     * @return settings sections.
     */
    static private void initSections() {
        sections = new ArrayList<>();
        
        // add new sections here
        // USER SECTIONS
        sections.add(new Section("General settings", "User"));
        sections.add(new Section("Amino Acid Graphs (AAG)", "User"));
        sections.add(new Section("Protein Graph (PG)", "User"));
        sections.add(new Section("Graph types", "User"));
        sections.add(new Section("Folding Graph (FG)", "User"));
        sections.add(new Section("Complex Graph (CG)", "User"));
        sections.add(new Section("Structure visualization", "User"));
        
        // ADVANCED SECTIONS
        sections.add(new Section("Prints / Error handling", "Advanced"));
        sections.add(new Section("Performance", "Advanced"));
        sections.add(new Section("Parser", "Advanced"));
        sections.add(new Section("Secondary Structure Element (SSE)", "Advanced"));
        sections.add(new Section("SSE orientation", "Advanced"));
        sections.add(new Section("Contact definition", "Advanced"));
        sections.add(new Section("AAG", "Advanced"));
        sections.add(new Section("Database (DB) connection", "Advanced"));
        sections.add(new Section("DB settings", "Advanced"));
        sections.add(new Section("DB structure search", "Advanced"));
        sections.add(new Section("FG", "Advanced"));
        sections.add(new Section("PG", "Advanced"));
        sections.add(new Section("CG", "Advanced"));
        sections.add(new Section("Motifs", "Advanced"));
        sections.add(new Section("Cluster mode", "Advanced"));
        sections.add(new Section("DB: representative chains", "Advanced"));
        sections.add(new Section("Image settings", "Advanced"));
        sections.add(new Section("Output settings", "Advanced"));
        sections.add(new Section("DB: similarity search", "Advanced"));        
        
        // DEVELOPER SECTIONS
        sections.add(new Section("Deprecated", "Developer"));
        sections.add(new Section("Debug", "Developer"));
        sections.add(new Section("Disabled", "Developer"));
        sections.add(new Section("Alternate AAG", "Developer"));
    }
    
    
    /**
     * Creates a map for settings names to section names. Always call after sections are created!
     */
    static private void initNamesMap() {
        mapSettingNameToSectionName = new HashMap<>();
        for (Section tmpSection : sections) {
            for (Setting tmpSetting : tmpSection.settings) {
                mapSettingNameToSectionName.put(tmpSetting.getName(), tmpSection.name);
            }
        }
    }
    
    
    /**
     * Creates a map for section names to their index. Always call after sections are created!
     */
    static private void initIndexMap() {
        mapSectionNameToSectionIndex = new HashMap<>();
        for (int i = 0; i < sections.size(); i++) {
            mapSectionNameToSectionIndex.put(sections.get(i).name, i);
        }
    }
    
    
    /**
     * Reads the settings from a file to set overwritten settings.
     * @param filepath use default file path if empty
     * @param rewrite whether file should be rewritten if not conform with file format
     * @return number of settings loaded
     */
    public static int loadFromFile(String filepath, Boolean rewrite) {
        filepath = (filepath.equals("") ? DEFAULT_FILE.getAbsolutePath() : filepath);  // use default filepath if not specified
        int numLoaded = 0;
        
        Properties readProperties = IO.readPropertiesFromFile(filepath);
        @SuppressWarnings("unchecked")  // skip warning associated to next line, according to https://www.boraji.com/how-to-iterate-properites-in-java
        Enumeration<String> enumer = (Enumeration<String>) readProperties.propertyNames();  // seems we have to create it outside the loop head
        while (enumer.hasMoreElements()) {
            String key = enumer.nextElement();
            String value = readProperties.getProperty(key);
            if (set(key, value)) { numLoaded++; } 
        }
        
        if (rewrite) {
            // check whether default file is correctly formatted and contains all existing settings
            String readFile = IO.readFileToString(filepath);
            if (! readFile.equals(asFormattedString())) {
                // file differs from how it should look: provide some info, copy old one and write new one
                // info
                ArrayList<String> unoverwrittenSettings = getUnoverwrittenSettings();
                if (unoverwrittenSettings.size() > 0) {
                    DP.getInstance().w(PACKAGE_TAG, "Settings file is missing following settings:");
                    for (String tmpSettingName : unoverwrittenSettings) {
                        DP.getInstance().w(tmpSettingName, 2);
                    }
                } else {
                    DP.getInstance().w(PACKAGE_TAG, "Settings file contains all settings, but is not formatted properly");
                }

                // copy old file
                String copiedFilepath = IO.getUniqueFilename(filepath + "_copied");
                DP.getInstance().w("Copying old settings file to '" + copiedFilepath + "' and writing correctly formatted settings file "
                        + "with all settings and your previous values to '" + filepath + "'. Nothing to thank ;-)", 2);
                IO.writeStringToFile(readFile, copiedFilepath, false);

                // rewrite
                IO.writeStringToFile(asFormattedString(), filepath, true);
            }
        }
        
        return numLoaded;
    }
    
    
    public static void init() {
        // read settings in if existing and create otherwise
        if (IO.fileExistsIsFileAndCanRead(DEFAULT_FILE)) {
            numLoadedSettings = loadFromFile("", true);
        } else {
            // create default settings file
            DP.getInstance().w(PACKAGE_TAG, "Could not load settings from properties file, trying to create it.");
            DP.getInstance().i(PACKAGE_TAG, "If you have used PTGLgraphComputation previously and just upgraded to a new version, you might want to use your "
                    + "previous settings file. It should be located at your home directory and be named '.PTGLgraphComputation_settings'. Simply rename it to "
                    + "'PTGLgraphComputation_settings.txt' and rerun PTGLgraphComputation once. It creates the newly formatted settings file for you with your previous settings.");
            if (IO.writeStringToFile(asFormattedString(), DEFAULT_FILE.getAbsolutePath(), false)) {
                createdDefaultFile = true;
            } else {
                DP.getInstance().w(PACKAGE_TAG, "Did not found default file and could not create it at '" + DEFAULT_FILE.getAbsolutePath() + "'. " +
                        "Trying to go on and using default values.");
            }
        }
    }
    
    
    static public String asFormattedString() {
        String formattedString = "";
        String lastSectionType = "";  // used to track in which type of section we are: User, Advanced, Developer

        // header
        formattedString += "##### PTGLgraphComputation SETTINGS #####\n\n";
        formattedString += "# This file contains the settings for PTGLtools' PTGLgraphComputation as key-value pairs per line.\n"
                + "# The character after 'PTGLgraphComputation' indicates which data type is expected: "
                + "B(oolean), S(tring), I(nteger) or F(loat)\n";
        formattedString += "# The file is structured in sections which either belong to user, advanced or developer settings.\n\n";
        
        // sections
        for (int i = 0; i < sections.size(); i++) {
            Section tmpSection = sections.get(i);
            
            // print type if encountering new one
            if (! lastSectionType.equals(tmpSection.type)) {
                formattedString += "### " + tmpSection.type + " settings ###\n\n";
                lastSectionType = tmpSection.type;
            }
            
            formattedString += tmpSection.asFormattedString() + "\n";
        }
        
        return formattedString;
    }
    
    
    static private Setting getSettingByName(String name) {
        if (mapSettingNameToSectionName.containsKey(name)) {
            return sections.get(mapSectionNameToSectionIndex.get(mapSettingNameToSectionName.get(name))).getSetting(name);
        } else {
            return null;
        }
    }
    
    
    /**
     * Retrieves the setting with key 'key' from the settings and returns it as a String. Note that it is considered a fatal error if no such key exists. Ask first using 'contains()' if you're not sure. :)
     * @param key the key to get
     * @return the value of the specified key
     */
    public static String get(String key) {
        if(mapSettingNameToSectionName.containsKey(key)) {
            return(getSettingByName(key).getValue());
        }
        else {
            DP.getInstance().e(PACKAGE_TAG, "No config file or default value for setting '" + key + "' exists, setting invalid. "
                    + "Please inform a developer.  Exiting now.");
            System.exit(1);                
            return("ERROR");    // for the IDE                    
        }
        
    }
    
    
    /**
     * Tries to cast the value of the property key 'key' to Integer and return it. If this fails it is considered a fatal error.
     * @param key the key of the properties hashmap
     * @return the value of the key as an Integer
     */
    public static Integer getInteger(String key) {
        Integer i = null;
        String s = get(key);

        try {
            i = Integer.valueOf(s);
        }
        catch (NumberFormatException e) {
            DP.getInstance().e(PACKAGE_TAG, "Could not load setting '" + key + "' from settings as an Integer, invalid format. Exiting now.");
            System.exit(1);
        }
        return(i);
    }
    
    
    /**
     * Tries to cast the value of the property key 'key' to Float and return it. If this fails it is considered a fatal error.
     * @param key the key of the properties hashmap
     * @return the value of the key as a Float
     */
    public static Float getFloat(String key) {
        Float f = null;
        String s = get(key);

        try {
            f = Float.valueOf(s);
        }
        catch (NumberFormatException e) {
            DP.getInstance().e(PACKAGE_TAG, "Could not load setting '" + key + "' from settings as an Float, invalid format. Exiting now.");
            System.exit(1);
        }
        return(f);
    }
    
    
    /**
     * Tries to extract the value of the property key 'key' as a Boolean and return it. If this fails it is considered a fatal error.
     * The only accepted string representations of Booleans are "true" and "false".
     * @param key the key of the properties hashmap
     * @return the value of the key as a Boolean
     */
    public static Boolean getBoolean(String key) {
        Boolean b = null;
        String s = null;

        s = get(key);

        if(s.toLowerCase().equals("true")) {
            return(true);
        }
        else if(s.toLowerCase().equals("false")) {
            return(false);
        }
        else {
            DP.getInstance().e(PACKAGE_TAG, "Could not load setting '" + key + "' from settings as an Boolean, invalid format. Exiting now.");
            System.exit(1);
            return(false);      // never reached
        }
    }
    
    
    static public Boolean set(String key, String value) {
        Setting targetSetting = getSettingByName(key);
        if (targetSetting != null) {
            return getSettingByName(key).setOverwrittenValue(value);
        } else {
            DP.getInstance().w(PACKAGE_TAG, "Could not find setting '" + key + "' to overwrite its value. Going on without changes to the setting.");
            return false;
        }
    }
    
    
    /**
     * Returns the version string. This is NOT guaranteed to be a number.
     * @return the PTGLgraphComputation version
     */
    public static String getVersion() {
        return(VERSION);
    }
    
    
    /**
     * Returns the application tag that is printed as a prefix for all output lines.
     * @return the apptag
     */
    public static String getApptag() {
        return("[PTGLgraphComputation] ");
    }
    
    
    /**
     * Returns the default config file location. Note that this may or may not be in use atm. Use getConfigFile() instead if you
     * need the file that is currently in use.
     * @return the path as a String
     */
    public static String getDefaultConfigFilePath() {
        return(DEFAULT_FILE.getAbsolutePath());
    }
    
    
    public static ArrayList<String> getUnoverwrittenSettings() {
        ArrayList<String> unoverwrittenSettings = new ArrayList<>();
        for (Section tmpSection : sections) {
            for (Setting tmpSetting : tmpSection.settings) {
                if (! tmpSetting.isOverwritten()) {
                    unoverwrittenSettings.add(tmpSetting.getName());
                }
            }
        }
        return unoverwrittenSettings;
    }
    
    
    /**
     * Creates an array of the output image formats for protein graphs which are set in the settings.
     * @return the output formats, collected from settings like 'PTGLgraphComputation_B_img_output_format_PNG'
     */
    public static IMAGEFORMAT[] getProteinGraphOutputImageFormats() {
        ArrayList<IMAGEFORMAT> formats = new ArrayList<>();
        
        if(getBoolean("PTGLgraphComputation_B_img_output_format_PNG")) {
            formats.add(IMAGEFORMAT.PNG);
        }
        if(getBoolean("PTGLgraphComputation_B_img_output_format_PDF")) {
            formats.add(IMAGEFORMAT.PDF);
        }
        // --- ignore SVG because it is always produced ---
        //if(Settings.getBoolean("PTGLgraphComputation_B_img_output_format_SVG")) {
        //    formats.add(IMAGEFORMAT.SVG);
        //}
        
        return (IMAGEFORMAT[])formats.toArray(new IMAGEFORMAT[formats.size()]);
    }
    
    /**
     * Creates an array of the output image formats for aa graphs which are set in the settings.
     * @return the output formats, collected from settings like 'PTGLgraphComputation_B_img_AAG_output_format_PNG'
     */
    public static IMAGEFORMAT[] getAminoAcidGraphOutputImageFormats() {
        ArrayList<IMAGEFORMAT> formats = new ArrayList<>();
        
        if(getBoolean("PTGLgraphComputation_B_img_AAG_output_format_PNG")) {
            formats.add(IMAGEFORMAT.PNG);
        }
        if(getBoolean("PTGLgraphComputation_B_img_AAG_output_format_PDF")) {
            formats.add(IMAGEFORMAT.PDF);
        }
       
        
        return (IMAGEFORMAT[])formats.toArray(new IMAGEFORMAT[formats.size()]);
    }
    
    /**
     * Creates an array of the output image formats for folding graphs which are set in the settings.
     * @return the output formats, collected from settings like 'PTGLgraphComputation_B_img_FG_output_format_PNG'
     */
    public static IMAGEFORMAT[] getFoldingGraphOutputImageFormats() {
        ArrayList<IMAGEFORMAT> formats = new ArrayList<>();
        
        if(getBoolean("PTGLgraphComputation_B_img_FG_output_format_PNG")) {
            formats.add(IMAGEFORMAT.PNG);
        }
        if(getBoolean("PTGLgraphComputation_B_img_FG_output_format_PDF")) {
            formats.add(IMAGEFORMAT.PDF);
        }
        // --- ignore SVG because it is always produced ---
        //if(Settings.getBoolean("PTGLgraphComputation_B_img_FG_output_format_SVG")) {
        //    formats.add(IMAGEFORMAT.SVG);
        //}
        
        return (IMAGEFORMAT[])formats.toArray(new IMAGEFORMAT[formats.size()]);
    }
    
    /**
     * Creates an array of the output image formats for complex graphs which are set in the settings.
     * @return the output formats, collected from settings like 'PTGLgraphComputation_B_img_CG_output_format_PNG'
     */
    public static IMAGEFORMAT[] getComplexGraphOutputImageFormats() {
        ArrayList<IMAGEFORMAT> formats = new ArrayList<IMAGEFORMAT>();
        
        if(getBoolean("PTGLgraphComputation_B_img_CG_output_format_PNG")) {
            formats.add(IMAGEFORMAT.PNG);
        }
        if(getBoolean("PTGLgraphComputation_B_img_CG_output_format_PDF")) {
            formats.add(IMAGEFORMAT.PDF);
        }
        // --- ignore SVG because it is always produced ---
        //if(Settings.getBoolean("PTGLgraphComputation_B_img_CG_output_format_SVG")) {
        //    formats.add(IMAGEFORMAT.SVG);
        //}
        
        return (IMAGEFORMAT[])formats.toArray(new IMAGEFORMAT[formats.size()]);
    }
    
    
    /**
     * Returns the line for images containing program name and version.
     * @return line containing program name and version
     */
    public static String getImagesCreatorVersionLine() {
        return "Created with " + getProgramName() + " version " + getVersion();
    }
    
    
    // ### simple getter and setter
    public static int getNumLoadedSettings() { return numLoadedSettings; }
    public static Boolean getCreatedDefaultFile() { return createdDefaultFile; }
    public static String getProgramName() { return PROGRAM_NAME; }
}
