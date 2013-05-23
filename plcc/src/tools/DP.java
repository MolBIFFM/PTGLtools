/*
 * This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.
 *
 * Copyright Tim Schäfer 2012. VPLG is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package tools;

import plcc.Settings;

/**
 * A debug printer. This is a singleton.
 * @author ts
 */
public class DP {
    private static DP instance = null;
    
    public static final String appTag = "[PLCC] ";
    
    public static final String errorTag = "[ERROR] ";
    public static final String warningTag = "[WARNING] ";
    public static final String infoTag = "[INFO] ";
    public static final String debugTag = "[DEBUG] ";
    public static final String messageTag = "[MSG] ";
    
    protected DP() {
        // prevent instantiation
    }
    
    public static DP getInstance() {
      if(instance == null) {
         instance = new DP();
         instance.init();
      }
      return instance;
   }
    
    private void init() {
        // init stuff here
    }
    
    public void w(String msg) {
        Boolean doWarn = true;
        try {
            doWarn = Settings.getBoolean("plcc_B_no_warn");
        } catch(Exception e) {
            // the settings have not been inited yet, so assume that we should warn
            //System.err.println("WARNING: No settings yet.");
        }
        if(doWarn) {
            System.err.println(DP.appTag + DP.warningTag + msg);
        }
    }
    
    /**
     * Prints an error. Newline is added at the end.
     * @param msg the message to print
     */
    public void e(String msg) {
        System.err.println(DP.appTag + DP.errorTag + msg);
    }
    
    
    /**
     * Prints a standard output message. Newline is added at the end.
     * @param msg the message to print
     */
    public void p(String msg) {
        System.out.println(DP.appTag + DP.messageTag + msg);
    }
    
    /**
     * Prints a standard output message to STDOUT without adding a newline. No other tags are added.
     * @param msg the message to print
     */
    public void pp(String msg) {
        System.out.print(msg);
    }
    
    /**
     * Prints an error message to STDERR without adding a newline. No other tags are added.
     * @param msg the message to print
     */
    public void ee(String msg) {
        System.err.print(msg);
    }
    
    
    /**
     * Prints an info message. Newline is added at the end.
     * @param msg the message to print
     */
    public void i(String msg) {
        System.out.println(DP.appTag + DP.infoTag + msg);
    }
    
    /**
     * Prints a debug message. Newline is added at the end.
     * @param msg the message to print
     */
    public void d(String msg) {
        System.out.println(DP.appTag + DP.debugTag + msg);
    }
    
    
}
