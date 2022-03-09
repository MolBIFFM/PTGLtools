/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2020. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */

package settings;

import java.util.ArrayList;
import tools.DP;

/**
 *
 * @author jnw
 */
class Setting {
    final private String name;  // name of the setting, should have format: PTGLgraphComputation_D_NAME, where D is the data type
    final private char dataType;  // values can only be saved as string, so we need to know, which data type they should have, capital letter
    final private String defaultValue;  // hard coded default value
    final private String documentation;  // documentation string
    private String overwrittenValue = UNOVERWRITTEN_PLACEHOLDER;  // the value of the setting that can be changed by the user via a settings file or command line options
    
    // contains the allowed data types for the check in the constructor as capital letters
    static final ArrayList<Character> ALLOWED_DATA_TYPES = new ArrayList<Character>() {
        {
            add('B');
            add('I');
            add('F');
            add('S');
            add('L');
        }
    };
    
    static final String UNOVERWRITTEN_PLACEHOLDER = "$UNSET$";

    
    /**
     * Constructor for SingleSetting.
     * @param name name of the setting, should have format: PTGLgraphComputation_D_NAME, where D is the data type
     * @param dataType data type of the setting
     * @param defaultValue hard coded default value of the setting
     * @param documentation
     */
    public Setting(String name, char dataType, String defaultValue, String documentation) {
        this.name = name;
        
        // check that data type is allowed to avoid sloppy programming
        if (ALLOWED_DATA_TYPES.contains(dataType)) {
            this.dataType = dataType;
        } else {
            DP.getInstance().e(Settings.PACKAGE_TAG, "Creating the setting '" + name + "' failed, because the data type is not allowed or correctly formatted. "
                    + "This is an error in the code, please inform the developer of this software. Exiting now.");
            this.dataType = 'S';  // for the IDE
            System.exit(1);
        }
        
        // if developer do no mistakes, this check is unneccessary - that is why we do the check
        if (checkDataType(defaultValue)) { 
            this.defaultValue = defaultValue;
        } else {
            DP.getInstance().e(Settings.PACKAGE_TAG, "Creating the setting '" + name + "' failed, because the default value '" + defaultValue + "' "
                    + "seems not to be of the required data type '" + getDataTypeString() + "'. "
                    + "This is an error in the code, please inform the developer of this software. Exiting now.");
            System.exit(1);
            this.defaultValue = "";  // for the IDE
        }
        
        this.documentation = documentation;
    }
    
    
    private Boolean checkDataType(String value) {
        switch (dataType) {
            case 'S':
            case 'L':
                return true;  // String is always String and every string is a list
            case 'I':
                if (tools.TextTools.isInteger(value)) {
                    return true;
                }
                break;
            case 'F':
                if (tools.TextTools.isFloat(value)) {
                    return true;
                }
                break;
            case 'B':
                if (tools.TextTools.isTrueOrFalse(value)) {
                    return true;
                }
                break;
        }
        return false;
    }
    
    
    /**
     * Sets the overwritten value after a type check.
     * @param value
     * @return 
     */
    Boolean setOverwrittenValue(String value) {
        if (checkDataType(value)) {
            if (! value.equals(UNOVERWRITTEN_PLACEHOLDER)) {
                this.overwrittenValue = value;
                return true;
            } else {
                DP.getInstance().w(Settings.PACKAGE_TAG, "Could not overwrite setting '" + name + "', because value '" + value + "' is used as internal "
                        + " placeholder for unset values. Please use something else. ");
            }
        } else {
            DP.getInstance().w(Settings.PACKAGE_TAG, "Could not overwrite setting '" + name + "', because value '" + value + "' is of wrong format. "
                    + "Required for this setting: " + getDataTypeString());
        }
        return false;
    }
    
    
    /**
     * Returns the value set by the user and otherwise the default value.
     * @return value as String
     */
    String getValue() {
        if (isOverwritten()) {
            return overwrittenValue;
        } else {
            // default value not overwritten
            // TODO if PTGLgraphComputation_B_warn_cfg_fallback_to_default warn
            // System.out.println("INFO: Settings: Using internal default value '" + s + "' for setting '" + key + "'. Edit config file to override.");
            return defaultValue;
        }
    }
    
    private String getDataTypeString() {
        switch (dataType) {
            case 'S': { return "String"; }
            case 'I': { return "Integer"; }
            case 'F': { return "Float"; }
            case 'B': { return "Boolean"; }
            default: { return "UNSPECIFIED DATATYPE"; }
        }
    }
    
    Boolean isOverwritten() { 
        return ! overwrittenValue.equals(UNOVERWRITTEN_PLACEHOLDER); }
    
    @Override
    public String toString() {
        return "[SETTING] " + name + " | default: '" + defaultValue + "' | overwritten: '" + overwrittenValue + "'";
    }
    
    String asFormattedString() {
        String formattedString = "";
        formattedString += "# " + documentation + "\n";  // documentation
        formattedString += "# Default: '" + defaultValue + "'\n";  // default value
        formattedString += name + "=" + getValue() + "\n";  // key-value pair
        return formattedString;
    }
    
    // ### simple getter / setter ###
    String getDefaultValue() { return defaultValue; }
    String getDocumentation() { return documentation; }
    String getName() { return name; }
}
