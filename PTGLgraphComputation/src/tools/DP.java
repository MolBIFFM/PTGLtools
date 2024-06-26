/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package tools;

import settings.Settings;

/**
 * A debug printer. This is a singleton.
 * @author ts
 */
public class DP {
    private static DP instance = null;
    
    public static final String appTag = Settings.getApptag();
    
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
    
    /**
     * Prints a warning. Newline is added at the end.
     * @param msg the message to print
     * @param optionalNumberLeadingSpaces vararg where first element is used to create that many leading spaces before warning
     */
    public void w(String msg, int... optionalNumberLeadingSpaces) {
        Boolean doWarn = true;
        try {
            doWarn = ( ! Settings.getBoolean("PTGLgraphComputation_B_no_warn"));
        } catch(Exception e) {
            // the settings have not been inited yet, so assume that we should warn
            //System.err.println("WARNING: No settings yet.");
        }
              
        // if optionalNumberLeadingSpaces is given, take first number for number of spaces. Zero spaces otherwise
        String leadingSpaces = new String(new char[optionalNumberLeadingSpaces.length > 0? optionalNumberLeadingSpaces[0] : 0]).replace("\0", " ");
        
        if(doWarn) {
            System.err.println(leadingSpaces + DP.appTag + DP.warningTag + msg);
        }
    }
    
    /**
     * Prints a warning. Newline is added at the end.
     * @param srcTag the source tag to print before the actual message
     * @param msg the message to print
     */
    public void w(String srcTag, String msg, int... optionalNumberLeadingSpaces) {
        Boolean doWarn = true;
        try {
            doWarn = ( ! Settings.getBoolean("PTGLgraphComputation_B_no_warn"));
        } catch(Exception e) {
            // the settings have not been inited yet, so assume that we should warn
            //System.err.println("WARNING: No settings yet.");
        }
        
        // if optionalNumberLeadingSpaces is given, take first number for number of spaces. Zero spaces otherwise
        String leadingSpaces = new String(new char[optionalNumberLeadingSpaces.length > 0? optionalNumberLeadingSpaces[0] : 0]).replace("\0", " ");
        
        if(doWarn) {
            System.err.println(leadingSpaces + DP.appTag + DP.warningTag + "[" + srcTag + "] " + msg);
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
     * Prints an error. Newline is added at the end.
     * @param srcTag the source tag to print before the actual message
     * @param msg the message to print
     */
    public void e(String srcTag, String msg) {
        System.err.println(DP.appTag + DP.errorTag + "[" + srcTag + "] " + msg);
    }
    
    /**
     * Prints a critical error and exits. Newline is added at the end.
     * @param msg the message to print
     */
    public void c(String srcTag, String msg) {
        System.err.println(DP.appTag + DP.errorTag + "[" + srcTag + "] " + msg);
        System.exit(1);
    }
    
    /**
     * Prints a critical error and exits. Newline is added at the end.
     * @param msg the message to print
     */
    public void c(String msg) {
        System.err.println(DP.appTag + DP.errorTag + msg);
        System.exit(1);
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
     * Prints an info message. Newline is added at the end.
     * @param msg the message to print
     */
    public void i(String srcTag, String msg) {
        System.out.println(DP.appTag + DP.infoTag + "[" + srcTag + "] " + msg);
    }
    
    /**
     * Prints a debug message. Newline is added at the end.
     * @param msg the message to print
     */
    public void d(String msg) {
        System.out.println(DP.appTag + DP.debugTag + msg);
    }
    
    /**
     * Prints a debug message. Newline is added at the end.
     * @param msg the message to print
     */
    public void d(String srcTag, String msg) {
        System.out.println(DP.appTag + DP.debugTag + "[" + srcTag + "] " + msg);
    }
    
    /**
     * Does nothing anymore in this implementation.
     */
    public void flush() {
        // placeholder function for backwards compatibility only, does nothing
    }
    
    
}
