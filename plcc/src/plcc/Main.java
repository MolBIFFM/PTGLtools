/*
 * This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.
 *
 * Copyright Tim Schäfer 2012. VPLG is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */


package plcc;

// imports
import algorithms.ConnectedComponents;
import proteingraphs.PTGLNotations;
import resultcontainers.ProteinResults;
import resultcontainers.PTGLNotationFoldResult;
import resultcontainers.ProteinFoldingGraphResults;
import resultcontainers.ProteinChainResults;
import resultcontainers.ComplexGraphResult;
import io.IO;
import io.FileParser;
import io.DBManager;
import proteinstructure.ProtMetaInfo;
import proteingraphs.FoldingGraphComparator;
import proteingraphs.MolContactInfo;
import proteingraphs.ComplexGraph;
import proteingraphs.ProtGraphs;
import proteingraphs.SSEComparator;
import proteingraphs.SSEGraph;
import proteingraphs.ContactMatrix;
import proteingraphs.ProtGraph;
import proteingraphs.FoldingGraph;
import graphdrawing.ProteinGraphDrawer;
import graphdrawing.DrawTools;
import graphformats.GraphFormats;
import proteinstructure.Model;
import proteinstructure.Residue;
import proteinstructure.Chain;
import proteinstructure.AminoAcid;
import proteinstructure.Atom;
import proteinstructure.Molecule;
import proteinstructure.RNA;
import proteinstructure.SSE;
import algorithms.GraphMetrics;
import algorithms.GraphPropResults;
import algorithms.GraphProperties;
import algorithms.GraphRandomizer;
import datastructures.AAGraph;
import datastructures.AAInteractionNetwork;
import datastructures.PPIGraph;
import datastructures.SimpleGraphInterface;
import datastructures.SparseGraph;
import htmlgen.CssGenerator;
import htmlgen.HtmlGenerator;
import htmlgen.JmolTools;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;
import java.awt.image.*;
import java.io.*;
import static java.lang.System.exit;
import java.nio.channels.*;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Locale;
import java.util.logging.Level;
import javax.imageio.ImageIO;
import javax.xml.parsers.ParserConfigurationException;
//import java.net.*;
//import org.jgrapht.*;
//import org.jgrapht.graph.*;
//import org.jgrapht.alg.ConnectivityInspector;
//import com.google.gson.*;
//import datastructures.UndirectedGraph;
import net.sourceforge.spargel.writers.GMLWriter;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.xml.sax.SAXException;
import graphdrawing.DrawTools.IMAGEFORMAT;
import graphdrawing.DrawableGraph;
import graphdrawing.IDrawableGraph;
import graphdrawing.SimpleGraphDrawer;
import java.util.concurrent.TimeUnit;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import parsers.GMLGraphParser;
import parsers.IGraphParser;
import similarity.CompareOneToDB;
import similarity.Similarity;
import similarity.SimilarityByGraphlets;
import tools.DP;
import tools.PiEffectCalculations;
import tools.PlccUtilities;
import tools.XMLContentHandlerPDBRepresentatives;
import tools.XMLErrorHandlerJAX;
import tools.XMLParserJAX;

/**
 * This is the Main class of plcc.
 *
 * The program 'plcc' (Protein ligand contact calculation) does a bit more than its name implies. It reads a PDB file
 * and calculates the contacts between all protein residues, treating their atoms as solid spheres with a radius of 2
 * Angstroem. It also calculates the contacts of HET groups (potential ligands) with the protein, but ignores stuff
 * like water (HOH, DOD), artificial atoms (Q) and solvent like methanol. Hydrogen atoms are also ignored if present
 * in the PDB file.
 *
 * If told so via command line parameters, the program also calculates the secondary structure elements (SSEs) based
 * on the results of the DSSP algorithm (the program 'dsspcmbi', see Kabsch & Sander 1983, 'Dictionary of protein secondary
 * structure: pattern recognition of hydrogen-bonded and geometrical features') and the contacts between them.
 *
 * It then creates a graph of the SSE connections, similar to the description in Koch 1997, 'Ein graphentheoretischer Ansatz
 * zum paarweisen und multiplen Vergleich von Proteinstrukturen' and draws it.
 *
 */

public class Main {


    // declare class vars
    //Logger logger = LogManager.getLogger(Main.class.getName());
    
    // the array used to store statistics on contacts between the different AA types
    public static final Integer NUM_AAs = 20 + 1 + 1;    // These are 20 real AAs, 1 extra for ligands and 1 for the total number count at index 0 (1st AA starts at index 1).
    public static final Integer MAX_ATOMS_PER_AA = 15;   // As everybody knows, TRP is the AA with most atoms and
                                            //  has 27 of them, but we ignore 'H' atoms (they are not
                                            //  included in most PDB files and if so, we filter them out.
                                            //  So only 14 of those atoms remain (but index starts at 1).
    
    /** The number of different contact types which are stored for a pair of residues. See calculateAtomContactsBetweenResidues() and
     the MolContactInfo class for details and usage. */
    public static final Integer NUM_RESIDUE_PAIR_CONTACT_TYPES = 12;
    
    /**
     * The number of different contacts types according to the alternative contact model which are stored for a pair of residues.
     * See calculateAtomContactsBetweenResiduesAlternativeModel() and the MolContactInfo class for details and usage.
     */
    public static final Integer NUM_RESIDUE_PAIR_CONTACT_TYPES_ALTERNATIVE_MODEL = 41;

    /**
     * The contacts of a chainName. The 4 fields are: AA 1 index, AA 2 index, atom index in AA 1, atom index in chainName 2.
     */
    static Integer[][][][] contact;

    static Integer globalMaxCenterSphereRadius;


    /** This is required for the additional speedup during the calculation of residue contacts that
      allows us to skip the next few residues if the distance between sequential neighbor residues
      is large (I. Koch). It is set by setGlobalMaxSeqNeighborResDist() and used by the function
      Main.calculateAllContacts(). */
    static Integer globalMaxSeqNeighborResDist;
    
    static final String version = Settings.getVersion();
    
    /**
     * Files added to this array during the run get deletes on program exit.
     */
    static ArrayList<File> deleteFilesOnExit;
    
    /** Whether the PDB file name given on the command line is used. This is not the case for command lines which only operate on the database or which need no input file (e.g., --recreate-tables). */
    static Boolean useFileFromCommandline = true;
    
    // Lists of residues and RNAs. They are initilized as null and created once from molecules the first time the function
    //   resFromMolecules or rnaFromMolecules is called
    static ArrayList<Residue> residues = null;
    static ArrayList<RNA> rnas = null;
    

    public static void checkArgsUsage(String[] args, Boolean[] argsUsed) {
        for(int i = 0; i < argsUsed.length; i++) {
            Boolean b = argsUsed[i];
            if( ! b) {
                System.err.println("WARNING: Command line argument #" + i +": '" + args[i] + "' not recognized or used.");
            }
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
                
        Date computationStartTime = new Date();
        StringBuilder outputToBePrintedUnlessSilent = new StringBuilder();
        
        outputToBePrintedUnlessSilent.append("[======================== plcc -- Protein-Ligand Contact Calculation ========================]\n");
        outputToBePrintedUnlessSilent.append("Init... (Version ").append(version).append(")\n");
        
        // *************************************************** load default settings from config file *************************************

        // The settings are defined in the Settings class. They are loaded from the config file below and can then be overwritten
        //  by command line arguments.
        Settings.init();
        
        // Check library dir, warn if not there        
        // NOTE: This is not necessary anymore, because we can store the libs inside the JAR using the jar-for-store target added in the build.xml file. Run 'ant ci' to generate the file, it will then be in store/plcc.jar.
        /*
        File libDir = new File("lib");
        if( ! libDir.exists()) {
            DP.getInstance().w("Library directory '" + libDir.toString() + "' not found. Libraries missing, will crash if functionality of them is required during this run.");
            DP.getInstance().w("Was this program executed from its installation directory? If not you need to copy the lib/ directory to the working directory.");
            DP.getInstance().w("INFO: The Java library path is set to: '" + System.getProperty("java.library.path") + "'.");
            DP.getInstance().w("INFO: The Java classloader path is set to: '" + System.getProperty("java.class.path") + "'.");
        }
        */

        
        int numSettingsLoaded = Settings.load("");
        if(numSettingsLoaded > 0) {             // Empty string means that the default file of the Settings class is used
            outputToBePrintedUnlessSilent.append("  Loaded ").append(numSettingsLoaded).append(" settings from properties file.\n");
            
        }
        else {
            DP.getInstance().w("Could not load settings from properties file, trying to create it.");
            if(Settings.createDefaultConfigFile()) {
                System.out.println("  Default config file created, will use it from now on.");
            } else {
                DP.getInstance().w("Could not create default config file, check permissions. Using internal default settings.");
            }
            Settings.resetAll();        // init settings with internal defaults for this run
        }

        //Settings.printAll();


        String fs = System.getProperty("file.separator");

        

        // ****************************************************    declarations    **********************************************************

        ArrayList<Model> models = new ArrayList<Model>();
        ArrayList<Chain> chains = new ArrayList<Chain>();
        ArrayList<Molecule> molecules = new ArrayList<Molecule>();
        ArrayList<Residue> residuesWithoutLigands;
        HashMap<Character, ArrayList<Integer>> sulfurBridges = new HashMap<Character, ArrayList<Integer>>();
        ArrayList<Atom> atoms = new ArrayList<Atom>();
        ArrayList<SSE> dsspSSEs = new ArrayList<SSE>();
        ArrayList<SSE> ptglSSEs = new ArrayList<SSE>();
        ArrayList<String> allModelsIDsOfWholePDBFile = new ArrayList<String>();
        deleteFilesOnExit = new ArrayList<File>();
        
        String pdbid = "";
        String dsspFile = "";
        String dsspLigFile = "";
        String pdbFile = "";
        String outputDir = ".";
        String pdbIdDotGeoFile = "";
        String pdbIdDotGeoLigFile = "";
        String conDotSetFile = "";
        String chainsFile = "";
        String ligandsFile = "";
        String modelsFile = "";
        String sseMappingsFile = "";
        String resMapFile = "";
        String convertModelsToChainsInputFile = "";
        String convertModelsToChainsOutputFile = "";
        
        Boolean compareResContacts = false;
        String compareResContactsFile = "";
        Boolean compareSSEContacts = false;
        String compareSSEContactsFile = "";                


        ArrayList<MolContactInfo> cInfo;
        
        // init contact statistics array
        contact = new Integer[NUM_AAs][NUM_AAs][MAX_ATOMS_PER_AA][MAX_ATOMS_PER_AA];
        for(Integer i = 0; i < NUM_AAs; i++) {
            for(Integer j = 0; j < NUM_AAs; j++) {
                for(Integer k = 0; k < MAX_ATOMS_PER_AA; k++) {
                    for(Integer l = 0; l < MAX_ATOMS_PER_AA; l++) {
                        contact[i][j][k][l] = 0;
                    }
                }
            }
        }

        
        // ****************************************************    parse args    **********************************************************
        if(args.length > 0) {

            Boolean[] argsUsed = new Boolean[args.length];
            Arrays.fill(argsUsed, Boolean.FALSE);
            
            if(args[0].equals("-h") || args[0].equals("--help") || args[0].equals("--helpdev")) {
                usage();
                if(args[0].equals("--helpdev")) {
                    usagedev();
                }
                argsUsed[0] = true;
                checkArgsUsage(args, argsUsed);
                System.exit(0);
            }

            // get pdbid from first arg
            pdbid = args[0];
            argsUsed[0] = true;
            
            if(! pdbid.equals("NONE")) {
                DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
                Calendar cal = Calendar.getInstance();
                outputToBePrintedUnlessSilent.append("  Starting computation for PDB ID '" + pdbid + "' at " + dateFormat.format(cal.getTime()) + ".\n");
            }

            final Integer expectedLengthPDBID = 4;
            if(pdbid.length() != expectedLengthPDBID) {
                DP.getInstance().w("The given PDB identifier '" + pdbid + "' has an unusual length of " + pdbid.length() + " characters, expected " + expectedLengthPDBID + ".");
            }
            
            // set default file names from the pdb id (these may be overwritten by args later)
            pdbFile = pdbid + ".pdb";
            dsspFile = pdbid + ".dssp";

            // parse the rest of the arguments, if any
            if(args.length > 1) {

                for (Integer i = 1; i < args.length; i++) {

                    String s = args[i];

                    if(s.equals("-h") || s.equals("--help")) {
                        usage();
                        argsUsed[i] = true;
                        checkArgsUsage(args, argsUsed);
                        System.exit(0);
                    }
                    
                    
                    if(s.equals("--graphlet-sim-method")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            String gm;
                            gm = args[i+1];
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            
                            if(gm.equals(SimilarityByGraphlets.GRAPHLET_SIM_METHOD_RGF)) {
                                Settings.set("plcc_S_search_similar_graphlet_scoretype", SimilarityByGraphlets.GRAPHLET_SIM_METHOD_RGF);
                            }
                            else if(gm.equals(SimilarityByGraphlets.GRAPHLET_SIM_METHOD_CUSTOM)) {
                                Settings.set("plcc_S_search_similar_graphlet_scoretype", SimilarityByGraphlets.GRAPHLET_SIM_METHOD_CUSTOM);
                            }
                            else {
                                System.err.println("ERROR: Invalid graphlet similarity method given. Use 'rgf' or 'cus'.");
                                syntaxError();
                            }                                
                        }
                    }

                    if(s.equals("-d") || s.equals("--dsspfile")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            dsspFile = args[i+1];
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    if(s.equals("-Y") || s.equals("--skip-vast")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            
                            try {
                                int tmp = Integer.parseInt(args[i+1]);
                            } catch(Exception e) {
                                syntaxError();
                            }
                            Settings.set("plcc_B_skip_too_large", "true");
                            Settings.set("plcc_I_skip_num_atoms_threshold", args[i+1]);
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    
                    
                    if(s.equals("--gz-dsspfile")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            String dsspFilenameGZ = args[i+1];
                            File dsspFileGZ = new File(dsspFilenameGZ);                            
                            File dsspFileUnpacked = null;
                            try {
                                // gunzip the file
                                //String tmpDir = dsspFileGZ.getAbsoluteFile().getParent();
                                String tmpDir = Settings.get("plcc_S_temp_dir");                                
                                dsspFileUnpacked = IO.unGzip(dsspFileGZ, new File(tmpDir));
                                deleteFilesOnExit.add(dsspFileUnpacked);
                                dsspFile = dsspFileUnpacked.toString();
                            }
                            catch(Exception e) {
                                System.err.println("ERROR: Could not extract input DSSP file in gz format at '" + dsspFilenameGZ + "'.");
                                //System.err.println("ERROR: The message was '" + e.getMessage() + "'.");                                
                                System.exit(1);
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    if(s.equals("--gz-pdbfile")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            String pdbFilenameGZ = args[i+1];
                            File pdbFileGZ = new File(pdbFilenameGZ);
                            File pdbFileUnpacked = null;                                                        
                            try {
                                //String tmpDir = pdbFileGZ.getAbsoluteFile().getParent();
                                String tmpDir = Settings.get("plcc_S_temp_dir");
                                //System.out.println("sourcedir='" + sourceDir + "'");
                                // gunzip the file
                                pdbFileUnpacked = IO.unGzip(pdbFileGZ, new File(tmpDir));
                                deleteFilesOnExit.add(pdbFileUnpacked);
                                pdbFile = pdbFileUnpacked.toString();
                            }
                            catch(Exception e) {
                                System.err.println("ERROR: Could not extract input PDB file in gz format at '" + pdbFilenameGZ + "'.");
                                //System.err.println("ERROR: The message was '" + e.getMessage() + "'.");
                                System.exit(1);
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }

                    
                    if(s.equals("--report-db-proteins")) {
                        useFileFromCommandline = false;
                        Settings.set("plcc_B_report_db_proteins", "true");
                        argsUsed[i] = true;
                    }
                    
                    
                                        
                    if(s.equals("--convert-models-to-chains")) {
                        convertModelsToChainsInputFile = null;
                        convertModelsToChainsOutputFile = null;
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            useFileFromCommandline = false;
                            Settings.set("plcc_B_convert_models_to_chains", "true");                        
                            convertModelsToChainsInputFile = pdbFile;
                            argsUsed[i] = true;
                            convertModelsToChainsOutputFile = args[i+1];
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    
                    if(s.equals("--alt-aa-contacts")) {
                        useFileFromCommandline = true;
                        Settings.set("plcc_B_alternate_aminoacid_contact_model", "true");
                        Settings.set("plcc_B_handle_hydrogen_atoms_from_reduce", "true");
                        //adjust other settings here
                        Settings.set("plcc_B_AAgraph_allchainscombined", "true");
                        Settings.set("plcc_B_aminoacidgraphs_include_ligands", "true");
                        Settings.set("plcc_B_quit_after_aag", "true");
                        outputToBePrintedUnlessSilent.append("  Using alternate residue contact model. Only computing PPI interface graphs.\n");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--alt-aa-contacts-ligands")) {
                        useFileFromCommandline = true;
                        Settings.set("plcc_B_alternate_aminoacid_contact_model_with_ligands", "true");
                        Settings.set("plcc_B_handle_hydrogen_atoms_from_reduce", "true");
                        //adjust other settings here
                        Settings.set("plcc_B_AAgraph_allchainscombined", "true");
                        Settings.set("plcc_B_aminoacidgraphs_include_ligands", "true");
                        Settings.set("plcc_B_quit_after_aag", "true");
                        outputToBePrintedUnlessSilent.append("  Using alternate residue contact model. Only computing PPI interface graphs.\n");
                        argsUsed[i] = true;
                    }
                    
                    
                    
                    if(s.equals("-p") || s.equals("--pdbfile")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            pdbFile = args[i+1];
                        }
                        argsUsed[i] = true;
                        argsUsed[i+1] = true;
                    }
                    
                    
                    if(s.equals("-C") || s.equals("--create-config")) {
                        useFileFromCommandline = false;
                        // The config file has already been created before parsing the command line if it did not exist, so we just do nothing here.
                        // We intentionally keep this option so we do not need to run other commands to get a config.
                        System.out.println("Tried to create PLCC config file at '" + Settings.getDefaultConfigFilePath() + "' (see above). Exiting.");
                        argsUsed[i] = true;
                        checkArgsUsage(args, argsUsed);
                        System.exit(0);
                    }
                    
                    if(s.equals("-N") || s.equals("--no-warn")) {
                        Settings.set("plcc_B_no_warn", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--cluster")) {                        
                        Settings.set("plcc_B_clustermode", "true");
                        Settings.set("plcc_B_silent", "true");
                        //Settings.set("plcc_B_write_chains_file", "true");
                        Settings.set("plcc_B_compute_motifs", "true");
                        Settings.set("plcc_B_complex_graphs", "true");
                        Settings.set("plcc_B_separate_contacts_by_chain", "false");
                        Settings.set("plcc_B_folding_graphs", "true");
                        Settings.set("plcc_B_draw_folding_graphs", "true");
                        Settings.set("plcc_B_output_images_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree", "true");
                        Settings.set("plcc_B_useDB", "true");
                        //Settings.set("plcc_B_AAgraph_allchainscombined", "true");
                        //Settings.set("plcc_B_AAgraph_perchain", "true");
                        Settings.set("plcc_B_no_warn", "true");
                        argsUsed[i] = true;
                        
                        //jnw_2019
                        System.out.println("NOTE: Following settings have been excluded from ");
                        System.out.println("  --cluster option:");
                        System.out.println("    plcc_B_write_chains_file");
                        System.out.println("    plcc_B_AAgraph_allchainscombined");
                        System.out.println("    plcc_B_AAgraph_perchain");
                        System.out.println("  Set them to true manually if you want them");
                    }
                    
                    if(s.equals("--write-chains-file")) {
                        Settings.set("plcc_B_write_chains_file", "true");
                        argsUsed[i] = true;
                    }                                                
                            
                    if(s.equals("-Z") || s.equals("--silent")) {
                        Settings.set("plcc_B_silent", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--less-output")) {
                        Settings.set("plcc_B_only_essential_output", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--verbose")) {
                        Settings.set("plcc_B_only_essential_output", "false");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--motifs")) {
                        Settings.set("plcc_B_compute_motifs", "true");
                        argsUsed[i] = true;
                    }
                    if(s.equals("--no-motifs")) {
                        Settings.set("plcc_B_compute_motifs", "false");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--db-batch")) {
                        Settings.set("plcc_B_db_use_batch_inserts", "true");
                        argsUsed[i] = true;
                    }
                    if(s.equals("--no-db-batch")) {
                        Settings.set("plcc_B_db_use_batch_inserts", "false");
                        argsUsed[i] = true;
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--compute-whole-db-graphlet-similarities")) {
                        useFileFromCommandline = false;
                        Settings.set("plcc_B_useDB", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_pg", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_cg", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_aag", "true");
                        argsUsed[i] = true;
                    }
                    if(s.equals("--compute-whole-db-graphlet-similarities-pg")) {
                        useFileFromCommandline = false;
                        Settings.set("plcc_B_useDB", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_pg", "true");
                        argsUsed[i] = true;
                    }
                    if(s.equals("--compute-whole-db-graphlet-similarities-cg")) {
                        useFileFromCommandline = false;
                        Settings.set("plcc_B_useDB", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_cg", "true"); 
                        argsUsed[i] = true;
                    }
                    if(s.equals("--compute-whole-db-graphlet-similarities-aag")) {
                        useFileFromCommandline = false;
                        Settings.set("plcc_B_useDB", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities", "true");
                        Settings.set("plcc_B_compute_graphlet_similarities_aag", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--compute-graph-metrics")) {
                        Settings.set("plcc_B_compute_graph_metrics", "true");
                        argsUsed[i] = true;
                    }
                    
                    
                    if(s.equals("--set-pdb-representative-chains-pre")) {
                        useFileFromCommandline = false;
                        if(args.length <= i+2 ) {
                            syntaxError("The --set-pdb-representative-chains-pre option requires an XML file to read the data from AND the info whether to remove old labels. For the latter, valid values are 'keep' or 'remove'.");
                        }
                        else {
                            useFileFromCommandline = false;                                                        
                            Settings.set("plcc_B_useDB", "true");
                            Settings.set("plcc_B_set_pdb_representative_chains_pre", "true");
                            Settings.set("plcc_S_representative_chains_xml_file", args[i+1]);
                            if(args[i+2].equals("keep")) {
                                Settings.set("plcc_B_set_pdb_representative_chains_remove_old_labels_pre", "false");
                            }
                            else if(args[i+2].equals("remove")) {
                                Settings.set("plcc_B_set_pdb_representative_chains_remove_old_labels_pre", "true");
                            }
                            else {
                                syntaxError("The --set-pdb-representative-chains-pre option requires an XML file to read the data from AND the info whether to remove old labels. For the latter, valid values are 'keep' or 'remove'.");
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            argsUsed[i+2] = true;
                        }
                    }
                    
                    if(s.equals("--set-pdb-representative-chains-post")) {
                        useFileFromCommandline = false;
                        if(args.length <= i+2 ) {
                            syntaxError("The --set-pdb-representative-chains-post option requires an XML file to read the data from AND the info whether to remove old labels. For the latter, valid values are 'keep' or 'remove'.");
                        }
                        else {
                            useFileFromCommandline = false;                                                        
                            Settings.set("plcc_B_useDB", "true");
                            Settings.set("plcc_B_set_pdb_representative_chains_post", "true");
                            Settings.set("plcc_S_representative_chains_xml_file", args[i+1]);
                            if(args[i+2].equals("keep")) {
                                Settings.set("plcc_B_set_pdb_representative_chains_remove_old_labels_post", "false");
                            }
                            else if(args[i+2].equals("remove")) {
                                Settings.set("plcc_B_set_pdb_representative_chains_remove_old_labels_post", "true");
                            }
                            else {
                                syntaxError("The --set-pdb-representative-chains-post option requires an XML file to read the data from AND the info whether to remove old labels. For the latter, valid values are 'keep' or 'remove'.");
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            argsUsed[i+2] = true;
                        }
                    }
                    
                    
                    
                    if(s.equals("-L") || s.equals("--lig-filter")) {
                        if(args.length <= i+2 ) {
                            syntaxError("The --lig-filter option requires two arguments, you can set one of them to zero if you want only one border.");
                        }
                        else {
                            Integer min = 0;
                            Integer max = 0;
                            try {
                                min = Integer.parseInt(args[i+1]);
                                max = Integer.parseInt(args[i+2]);                                
                            } catch (Exception e) {
                                syntaxError("The minimum and maximum atom count for ligands have to be numeric / integers.");
                            }
                            
                            if(min > max && max != 0) {
                                syntaxError("Minimum ligand atom number must not be greater than maximum ligand atom number (unless max = 0).");
                            } else {
                                Settings.set("plcc_I_lig_min_atoms", min.toString());
                                Settings.set("plcc_I_lig_max_atoms", max.toString());
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            argsUsed[i+2] = true;
                        }
                    }
                    
                    
                    if(s.equals("-S") || s.equals("--sim-measure")) {
                        System.out.println("Setting similarity measure to " + args[i+1] + ".");
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            Settings.set("plcc_S_search_similar_method", args[i+1]);
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;                            
                        }
                    }
                    
                    
                    
                    
                    if(s.equals("-M") || s.equals("--similar")) {
                        useFileFromCommandline = false;
                        if(args.length <= i+3 ) {
                            syntaxError();
                        }
                        else {                            
                            Settings.set("plcc_B_search_similar", "true");
                            Settings.set("plcc_B_search_similar_PDBID", args[i+1]);
                            Settings.set("plcc_B_search_similar_chainID", args[i+2]);
                            Settings.set("plcc_S_search_similar_graphtype", args[i+3]);                            
                        }
                        argsUsed[i] = true;
                        argsUsed[i+1] = true;
                        argsUsed[i+2] = true;
                        argsUsed[i+3] = true;
                    }
                    
                                                                               
                    if(s.equals("-v") || s.equals("--del-db-protein")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            String pdbidToDelete = args[i+1];
                            Integer numRows = 0;
                            System.out.println("Deleting protein with PDB identifier '" + pdbidToDelete + "' from database...");
                            
                            if(DBManager.init(Settings.get("plcc_S_db_name"), Settings.get("plcc_S_db_host"), Settings.getInteger("plcc_I_db_port"), Settings.get("plcc_S_db_username"), Settings.get("plcc_S_db_password"), Settings.getBoolean("plcc_B_db_use_autocommit"))) {
                                try {
                                    
                                    numRows = DBManager.deletePdbidFromDB(pdbidToDelete);
                                } catch(Exception e) {
                                    System.err.println("ERROR: Deleting protein failed: '" + e.getMessage() + "'. Exiting.");
                                    //System.exit(1);
                                }
                            }
                                                        
                            if(numRows > 0) {
                                System.out.println("Protein deleted from database, " + numRows + " rows affected. Exiting.");
                            }
                            else {
                                System.out.println("Protein was not in the database, " + numRows + " rows affected. Exiting.");
                            }
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            checkArgsUsage(args, argsUsed);
                            System.exit(0);
                        }
                    }
                    
                    if(s.equals("--check-whether-in-db")) {
                        argsUsed[i] = true;                        
                        if(DBManager.init(Settings.get("plcc_S_db_name"), Settings.get("plcc_S_db_host"), Settings.getInteger("plcc_I_db_port"), Settings.get("plcc_S_db_username"), Settings.get("plcc_S_db_password"), Settings.getBoolean("plcc_B_db_use_autocommit"))) {
                            Boolean found;
                            try {                                    
                                found = DBManager.proteinExistsInDB(pdbid);
                                if(found) {
                                    if( ! Settings.getBoolean("plcc_B_silent")) {
                                        System.out.println("Protein '" + pdbid + "' found in database (exiting with return code 0).");
                                    }
                                    checkArgsUsage(args, argsUsed);
                                    System.exit(0);
                                }
                                else {
                                    if( ! Settings.getBoolean("plcc_B_silent")) {
                                        System.out.println("Protein '" + pdbid + "' NOT found in database (exiting with return code 1).");
                                    }
                                    checkArgsUsage(args, argsUsed);
                                    System.exit(1);
                                }
                            } catch(Exception e) {
                                System.err.println("ERROR: Checking for protein failed: '" + e.getMessage() + "' (exiting with return code 2)");
                                checkArgsUsage(args, argsUsed);
                                System.exit(2);
                            }
                        }
                        else {
                            System.err.println("ERROR: Could not check whether protein exists in database, connection failed (exiting with return code 2).");
                            checkArgsUsage(args, argsUsed);
                            System.exit(2);
                        }
                    }
                    
                    
                    if(s.equals("-j") || s.equals("--ddb")) {
                        if(args.length <= i+4 ) {
                            syntaxError();
                        }
                        else {
                            Settings.set("plcc_B_useDB", "true");
                            String a_pdbid = args[i+1];
                            String a_chain = args[i+2];
                            String a_gt = args[i+3];
                            String a_outFile = args[i+4];
                                                        
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            argsUsed[i+2] = true;
                            argsUsed[i+3] = true;
                            argsUsed[i+4] = true;
                            System.out.println("Retrieving " + a_gt + " graph of PDB entry " + a_pdbid + ", chain " + a_chain + " from database.");
                            
                            if(DBManager.init(Settings.get("plcc_S_db_name"), Settings.get("plcc_S_db_host"), Settings.getInteger("plcc_I_db_port"), Settings.get("plcc_S_db_username"), Settings.get("plcc_S_db_password"), Settings.getBoolean("plcc_B_db_use_autocommit"))) {
                                drawPlccGraphFromDB(a_pdbid, a_chain, a_gt, a_outFile + Settings.get("plcc_S_img_output_fileext"), false);
                                System.out.println("Handled " + a_gt + " graph of PDB entry " + a_pdbid + ", chain " + a_chain + ", exiting.");
                                checkArgsUsage(args, argsUsed);
                                System.exit(0);
                            }
                            else {
                                checkArgsUsage(args, argsUsed);
                                System.exit(1);
                            }                                                                                    
                        }
                    }
                    
                    
                    if(s.equals("-x") || s.equals("--check-rescts")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            compareResContacts = true;
                            compareResContactsFile = args[i+1];
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    if(s.equals("-X") || s.equals("--check-ssects")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            Settings.set("plcc_B_debug_compareSSEContacts", "true");
                            Settings.set("plcc_S_debug_compareSSEContactsFile", args[i+1]);
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }                                        
                    
                    
                    if(s.equals("-y") || s.equals("--write-geodat")) {
                        Settings.set("plcc_B_ptgl_geodat_output", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--force")) {
                        Settings.set("plcc_F_abort_if_pdb_resolution_worse_than", "-1.0");
                        Settings.set("plcc_I_abort_if_num_molecules_below", "-1");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-E") || s.equals("--separate-contacts")) {
                        Settings.set("plcc_B_separate_contacts_by_chain", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-G") || s.equals("--complex-graphs")) {
                        Settings.set("plcc_B_complex_graphs", "true");
                        Settings.set("plcc_B_separate_contacts_by_chain", "false");
                        argsUsed[i] = true;
                    }
                    
                    
                    if(s.equals("-z") || s.equals("--ramaplot")) {
                        Settings.set("plcc_B_ramachandran_plot", "true");
                        argsUsed[i] = true;
                    }
                    

                    if(s.equals("-o") || s.equals("--outputdir")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            outputDir = args[i+1];
                            Settings.set("plcc_S_output_dir", outputDir);
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }

                    //if(s.equals("-i") || s.equals("--ignore_ligands_geolig")) {
                    //    System.out.println("WARNING: Ignoring of ligands not fully implemented yet. Contacts still include them.");
                    //    System.setProperty("plcc.useLigands", "false");
                    //}

                    if(s.equals("-r") || s.equals("--recreate-tables") || s.equals("--recreate-tables-empty")) {
                        Boolean fillTypes = ( ! s.equals("--recreate-tables-empty"));
                        if(fillTypes) {
                            System.out.println("Recreating DB tables and adding base type data...");
                        } else {
                            System.out.println("Recreating DB tables, not adding any base type data...");
                        }
                        if(DBManager.init(Settings.get("plcc_S_db_name"), Settings.get("plcc_S_db_host"), Settings.getInteger("plcc_I_db_port"), Settings.get("plcc_S_db_username"), Settings.get("plcc_S_db_password"), true)) {
                            
                            if(DBManager.dropTables()) {
                                System.out.println("  DB: Tried to drop statistics tables (no error messages => OK).");
                            }
                            if(DBManager.createTables(fillTypes)) {
                                System.out.println("  DB: Tried to create statistics tables (no error messages => OK).");
                            }
                            
                            boolean printTables = true;
                            if(printTables) {
                                StringBuilder sb = new StringBuilder();
                                ArrayList<String> t = DBManager.getPlccTablesCurrentlyInDatabase();
                                sb.append("There are now " + t.size() + " tables in the database " + DBManager.getConnectionInfoDatabaseName() + ": ");
                                for (int j = 0; j < t.size(); j++) {
                                    sb.append(t.get(j));
                                    if(j < t.size() - 1) {
                                        sb.append(", ");
                                    }
                                }
                                sb.append("\n");
                                sb.append("NOTE: You may want to set representative chains (PRE) now if needed.");
                                System.out.println(sb.toString());
                            }
                        }
                        else {
                            System.err.println("ERROR: Could not modify tables, DB connection failed.");
                        }
                        // exit
                        System.out.println("Done recreating DB tables, exiting.");
                        argsUsed[i] = true;
                        Main.checkArgsUsage(args, argsUsed);
                        Main.doExit(0);
                    }
                    
                    if(s.equals("-s") || s.equals("--draw-linnots")) {
                        Settings.set("plcc_B_folding_graphs", "true");
                        Settings.set("plcc_B_draw_folding_graphs", "true");
                        argsUsed[i] = true;
                    }

                    if(s.equals("-a") || s.equals("--include-coils")) {
                        Settings.set("plcc_B_include_coils", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-B") || s.equals("--force-backbone")) {
                        Settings.set("plcc_B_forceBackboneContacts", "true");
                        argsUsed[i] = true;
                    }
                   

                    if(s.equals("-w") || s.equals("--dont-write-images")) {
                        Settings.set("plcc_B_draw_graphs", "false");
                        argsUsed[i] = true;
                    }

                    if(s.equals("-c") || s.equals("--dont-calc-graphs")) {
                        Settings.set("plcc_B_calc_draw_graphs", "false");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--no-ptgl")) {
                        Settings.set("plcc_B_strict_ptgl_behaviour", "false");
                        argsUsed[i] = true;
                    }
                    

                    if(s.equals("-u") || s.equals("--use-database")) {
                        Settings.set("plcc_B_useDB", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-f") || s.equals("--folding-graphs")) {
                        Settings.set("plcc_B_folding_graphs", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-k") || s.equals("--img-dir-tree")) {
                        Settings.set("plcc_B_output_images_dir_tree", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-K") || s.equals("--graph-dir-tree")) {
                        Settings.set("plcc_B_output_textfiles_dir_tree", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--output-subdir-tree")) {
                        Settings.set("plcc_B_output_images_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-W") || s.equals("--output-www")) {
                        Settings.set("plcc_B_output_images_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree_html", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-H") || s.equals("--output-www-with-core")) {
                        Settings.set("plcc_B_output_images_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree_html", "true");
                        Settings.set("plcc_B_output_textfiles_dir_tree_core_html", "true");
                        argsUsed[i] = true;
                    }
                    
                    
                    
                    if(s.equals("-m") || s.equals("--image-format")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            String format = args[i+1].toUpperCase();
                            
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            
                            if((format.equals("PNG"))) {
                                Settings.set("plcc_S_img_output_format", "PNG");
                                Settings.set("plcc_S_img_output_fileext", ".png");
                            }
                            else if((format.equals("JPG"))) {
                                Settings.set("plcc_S_img_output_format", "JPG");
                                Settings.set("plcc_S_img_output_fileext", ".jpg");
                            }
                            else if((format.equals("TIF"))) {
                                Settings.set("plcc_S_img_output_format", "TIF");
                                Settings.set("plcc_S_img_output_fileext", ".tif");
                            }
                            else {
                                System.err.println("ERROR: Requested image output format '" + format + "' invalid. Use 'PNG' or 'SVG' for bitmap or vector output.");
                                checkArgsUsage(args, argsUsed);
                                syntaxError();
                            }                            
                        }
                    }
                    
                    
                    if(s.equals("--contact-level-debugging")) {                                                                        
                        Settings.set("plcc_B_contact_debug_dysfunct", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--draw-aag")) {                                                                        
                        Settings.set("plcc_B_draw_aag", "true");
                        argsUsed[i] = true;
                    }
                    
                    

                    if(s.equals("-n") || s.equals("--textfiles")) {
                        Settings.set("plcc_B_ptgl_text_output", "true");
                        argsUsed[i] = true;
                    }
                    
                    
                    if(s.equals("-q") || s.equals("--fg-notations")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            // We want to calculate folding graphs!
                            Settings.set("plcc_B_folding_graphs", "true");
                            
                            // If the user specifies the folding graph types manually, all those
                            // which are NOT listed default to 'off':                                                        
                            Settings.set("plcc_B_foldgraphtype_KEY", "false");
                            Settings.set("plcc_B_foldgraphtype_ADJ", "false");
                            Settings.set("plcc_B_foldgraphtype_RED", "false");
                            Settings.set("plcc_B_foldgraphtype_SEQ", "false");
                            
                            // Now add the listed ones back:
                            String types = args[i+1].toLowerCase();
                            
                            Integer nv = 0; // number of valid folding graph identifiers
                            
                            if(types.contains("k")) { Settings.set("plcc_B_foldgraphtype_KEY", "true"); nv++; }
                            if(types.contains("a")) { Settings.set("plcc_B_foldgraphtype_ADJ", "true"); nv++; }
                            if(types.contains("r")) { Settings.set("plcc_B_foldgraphtype_RED", "true"); nv++; }
                            if(types.contains("s")) { Settings.set("plcc_B_foldgraphtype_SEQ", "true"); nv++; }
                            
                            // sanity check
                            if(nv != types.length()) {
                                DP.getInstance().w("List of folding graph notations given on command line '" + types + "' contains invalid chars (" + types.length() + " given, " + nv + " valid).");
                                DP.getInstance().w("Valid chars: 'k' => KEY, 'a' => ADJ, 'r' => RED, 's' => SEQ. Example: '-q kr'");
                                
                                if(nv <= 0) {
                                    syntaxError();
                                }
                            }
                        }
                        
                    }
                    

                    if(s.equals("-g") || s.equals("--sse-graphtypes")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            // We want to calculate SSE graphs!
                            Settings.set("plcc_B_calc_draw_graphs", "true");
                            
                            // If the user specifies the folding graph types manually, all those
                            // which are NOT listed default to 'off':                                                        
                            Settings.set("plcc_B_graphtype_albe", "false");
                            Settings.set("plcc_B_graphtype_albelig", "false");
                            Settings.set("plcc_B_graphtype_alpha", "false");
                            Settings.set("plcc_B_graphtype_alphalig", "false");
                            Settings.set("plcc_B_graphtype_beta", "false");
                            Settings.set("plcc_B_graphtype_betalig", "false");
                            
                            // Now add the listed ones back:
                            String types = args[i+1].toLowerCase();
                            
                            Integer nv = 0; // number of valid folding graph identifiers
                            
                            if(types.contains("a")) { Settings.set("plcc_B_graphtype_alpha", "true"); nv++; }
                            if(types.contains("b")) { Settings.set("plcc_B_graphtype_beta", "true"); nv++; }
                            if(types.contains("c")) { Settings.set("plcc_B_graphtype_albe", "true"); nv++; }
                            if(types.contains("d")) { Settings.set("plcc_B_graphtype_alphalig", "true"); nv++; }
                            if(types.contains("e")) { Settings.set("plcc_B_graphtype_betalig", "true"); nv++; }
                            if(types.contains("f")) { Settings.set("plcc_B_graphtype_albelig", "true"); nv++; }
                            
                            // sanity check
                            if(nv != types.length()) {
                                DP.getInstance().w("List of folding graph notations given on command line '" + types + "' contains invalid chars (" + types.length() + " given, " + nv + " valid).");
                                DP.getInstance().w("Valid chars: 'a' => alpha, 'b' => beta, 'c' => albe, 'd' => alphalig, 'e' => betalig, 'f' => albelig. Example: '-g ace'");
                                
                                if(nv <= 0) {
                                    syntaxError();
                                }
                            }
                        }
                        
                    }
                    
                    
                    if(s.equals("-O") || s.equals("--outputformats")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            // If the user specifies the graph output formats manually, all those
                            // which are NOT listed default to 'off':                                                        
                            Settings.set("plcc_B_output_GML", "false");
                            Settings.set("plcc_B_output_TGF", "false");
                            Settings.set("plcc_B_output_DOT", "false");
                            Settings.set("plcc_B_output_kavosh", "false");
                            Settings.set("plcc_B_output_eld", "false");
                            Settings.set("plcc_B_output_plcc", "false");
                            Settings.set("plcc_B_output_perlfg", "false");
                            Settings.set("plcc_B_output_json", "false");
                            Settings.set("plcc_B_output_xml", "false");
                            Settings.set("plcc_B_output_gexf", "false");
                            Settings.set("plcc_B_output_cytoscapejs", "false");
                            
                            
                            
                            // Now add the listed ones back:
                            String types = args[i+1].toLowerCase();
                            
                            if(types.equals("x")) {
                                // Do nothing more, x means the user wants none of the formats
                            } else {
                            
                                Integer nv = 0; // number of valid folding graph identifiers

                                if(types.contains("g")) { Settings.set("plcc_B_output_GML", "true"); nv++; }
                                if(types.contains("t")) { Settings.set("plcc_B_output_TGF", "true"); nv++; }
                                if(types.contains("d")) { Settings.set("plcc_B_output_DOT", "true"); nv++; }
                                if(types.contains("k")) { Settings.set("plcc_B_output_kavosh", "true"); nv++; }
                                if(types.contains("e")) { Settings.set("plcc_B_output_eld", "true"); nv++; }
                                if(types.contains("p")) { Settings.set("plcc_B_output_plcc", "true"); nv++; }
                                if(types.contains("l")) { Settings.set("plcc_B_output_perlfg", "true"); nv++; }
                                if(types.contains("j")) { Settings.set("plcc_B_output_json", "true"); nv++; }
                                if(types.contains("m")) { Settings.set("plcc_B_output_xml", "true"); nv++; }
                                if(types.contains("f")) { Settings.set("plcc_B_output_gexf", "true"); nv++; }
                                if(types.contains("c")) { Settings.set("plcc_B_output_cytoscapejs", "true"); nv++; }
                                                                                                

                                // sanity check
                                if(nv != types.length()) {
                                    DP.getInstance().w("List of output formats given on command line '" + types + "' contains invalid chars (" + types.length() + " given, " + nv + " valid).");
                                    DP.getInstance().w("Valid chars: 'g' => GML, 't' => TGF, 'd' => DOT lang, 'k' => kavosh, 'e' => edge list, 'p' => PLCC, 'l' => Perf FG, 'j' => JSON, 'm' => XML, 'f' => GEXF, 'c' => CytoscapeJS. Example: '-O tgp'");

                                    if(nv <= 0) {
                                        syntaxError();
                                    }
                                }
                            }
                        }
                        
                    }
                    
                    
                    

                    if(s.equals("--debug") || s.equals("-D")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            Settings.set("plcc_I_debug_level", args[i+1]);
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    if(s.equals("--force-chain") || s.equals("-e")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            Settings.set("plcc_B_force_chain", "true");
                            Settings.set("plcc_S_forced_chain_id", args[i+1]);    
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                        }
                    }
                    
                    
                    

                    if(s.equals("-t") || s.equals("--draw-tgf-graph")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            Settings.set("plcc_B_graphimg_header", "false");
                            System.out.println("Drawing custom graph in TGF format from file '" + args[i+1] + "'.");
                            IMAGEFORMAT[] formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG };
                            drawTGFGraph(args[i+1], args[i+1], formats, new HashMap<Integer, String>());
                            System.out.println("Done drawing TGF graph, exiting.");
                            checkArgsUsage(args, argsUsed);
                            System.exit(1);
                        }
                    }
                    
                    if(s.equals("--draw-gml-graph") || s.equals("--props-gml-graph")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            Settings.set("plcc_B_graphimg_header", "false");
                            
                            Boolean computeSparseGraphProps = true;
                            Boolean alsoComputePropsForLargestCC = true;    // will only be done if the graph consists of more than 1 CC
                            Boolean skipDrawing = true;
                            
                            if(s.equals("--draw-gml-graph")){
                                DP.getInstance().i("Drawing custom graph in GML format from file '" + args[i+1] + "'.");
                                computeSparseGraphProps = false;
                                skipDrawing = false;
                            }
                            else {
                                DP.getInstance().i("Computing custom graph properties for graph in GML format from file '" + args[i+1] + "'.");
                                computeSparseGraphProps = true;
                                skipDrawing = true;
                            }
                            
                            
                            IMAGEFORMAT[] formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG };
                            String gml = null;
                            
                            try {
                                gml = FileParser.slurpFileToString(args[i+1]);  
                                //System.out.println("Received GML string: " + gml);
                            } catch (Exception e) {
                                DP.getInstance().e("Could not read GML graph file: " + e.getMessage());
                                System.exit(1);
                            }
                            
                            Map<Integer, String> vertexMappings = new HashMap<>();
                            String outFilePathNoExt = args[i+1]; // Use the name of the input file as output. Note that the image extension will be appended to this.
                            
                            // read vertex mappings from file if a file has been specified
                            if(args.length > i+2) {
                                argsUsed[i+2] = true;
                                DP.getInstance().i("Parsing vertex mappings file...");
                                try {
                                    vertexMappings = IO.parseMappingsFile(args[i+2]);
                                } catch (Exception e) {
                                    DP.getInstance().e("Could not read vertex mappings file: " + e.getMessage());
                                    System.exit(1);
                                }
                                DP.getInstance().i("Received " + vertexMappings.keySet().size() + " vertex mappings.");
                            }                            
                            
                            if(args.length > i+3) {
                                argsUsed[i+3] = true;
                                outFilePathNoExt = args[i+3];
                                DP.getInstance().i("Using custom image output base file path '" + outFilePathNoExt + "'.");
                            }
                            
                            GMLGraphParser p = new GMLGraphParser(gml);
                            IDrawableGraph g = p.getDrawableGraph();
                            DP.getInstance().i("Received graph with " + g.getDrawableVertices().size() + " vertices and " + g.getDrawableEdges().size() + " edges.");
                            
                            
                            
                            
                            if(computeSparseGraphProps) {
                                Date propsComputationStartTime, propsComputationEndTime; long timeDiff;
                                String unique_name;
                                String description = null;
                                Boolean isForLargestConnectedComponent = false;
                                Integer num_verts;
                                Integer num_edges;
                                Integer min_degree;
                                Integer max_degree;
                                Integer num_connected_components;
                                Integer diameter;
                                Integer radius;
                                Double avg_cluster_coeff;
                                Double avg_shortest_path_length;
                                Integer[] degreedist;
                                Double avg_degree;
                                Double density;
                                Integer[] cumul_degreedist;
                                Long runtime_secs;
                                
                                propsComputationStartTime = new Date();
                                DP.getInstance().i("Main", "Constructing sparse graph...");
                                SparseGraph<String, String> sg = g.toSparseGraph();
                                DP.getInstance().i("Main", "Computing sparse graph properties for SG with " + sg.getNumVertices() +" verts and " + sg.getNumEdges() + " edges ...");
                                GraphProperties gp = new GraphProperties(sg);
                                //DP.getInstance().i("Main", "Props for SG constructed.");
                                GraphPropResults gpr = gp.getGraphPropResults();
                                Boolean alsoPrintGraphProps = false;
                                if(alsoPrintGraphProps) {                                    
                                    System.out.println(GraphProperties.getOverviewPropsString(true, gpr));
                                    //DP.getInstance().i("Main", "Props for SG computed and printed.");                                                                    
                                }
                                
                                unique_name = "g" + System.currentTimeMillis(); // has to be kept for the CC as well
                                
                                num_verts = gpr.numVertices;
                                num_edges = gpr.numEdges;
                                min_degree = gpr.maxDegree;
                                max_degree = gpr.minDegree;
                                num_connected_components = gpr.numConnectedComponents;
                                diameter = gpr.graphDiameter;
                                radius = gpr.graphRadius;
                                avg_cluster_coeff = gpr.averageClusterCoefficient;
                                avg_shortest_path_length = gpr.averageShortestPathLength;
                                degreedist = gp.getDegreeDistributionUpTo(50);
                                avg_degree = gpr.averageDegree;
                                density = gpr.density;
                                cumul_degreedist = gp.getCumulativeDegreeDistributionUpToAsArray(50);

                                
                                propsComputationEndTime = new Date();
                                timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);
                                System.out.println("Computed graph properties for custom graph with " + gp.getNumVertices() + " verts and " + gp.getNumEdges() + " edges in " + runtime_secs + " seconds.");
                                

                                if(Settings.getBoolean("plcc_B_useDB")) {
                                    try {
                                        DBManager.deleteCustomGraphStatsFromDBByUniqueName(unique_name);
                                        DBManager.writeCustomgraphStatsToDB(unique_name, description, isForLargestConnectedComponent, num_verts, num_edges, min_degree, max_degree, num_connected_components, diameter, radius, avg_cluster_coeff, avg_shortest_path_length, degreedist, avg_degree, density, cumul_degreedist, runtime_secs);
                                        System.out.println("Graph properties for custom graph written to database, graph has been assigned unique name '" + unique_name + "' in db.");
                                    }catch(SQLException e) {
                                        DP.getInstance().e("Main", "Could not write custom graph stats to db: '" + e.getMessage() + "'.");
                                    }
                                }
                                
                                if(alsoComputePropsForLargestCC) {
                                    if( ! gp.getGraphIsConnected()) {                                        
                                        propsComputationStartTime = new Date();
                                        SimpleGraphInterface lcc = gp.getLargestConnectedComponent();
                                        DP.getInstance().i("Main", "Graph is not connected. Also computing sparse graph properties for its largest CC with " + lcc.getSize() +" verts ...");
                                        
                                        GraphProperties lccp = new GraphProperties(lcc);
                                        
                                        isForLargestConnectedComponent = true;
                                        num_verts = lccp.getNumVertices();
                                        num_edges = lccp.getNumEdges();
                                        min_degree = lccp.getMinDegree();
                                        max_degree = lccp.getMaxDegree();
                                        num_connected_components = lccp.getConnectedComponents().size();
                                        diameter = lccp.getGraphDiameter();
                                        radius = lccp.getGraphRadius();
                                        avg_cluster_coeff = lccp.getAverageClusterCoefficient();
                                        avg_shortest_path_length = lccp.getAverageShortestPathLength();
                                        degreedist = lccp.getDegreeDistributionUpTo(50);
                                        avg_degree = lccp.getAverageDegree();
                                        density = lccp.getDensity();
                                        cumul_degreedist = lccp.getCumulativeDegreeDistributionUpToAsArray(50);
                                        
                                        propsComputationEndTime = new Date();
                                        timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                        runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);
                                        System.out.println("Computed graph properties for largest CC of custom graph with " + lccp.getNumVertices() + " verts and " + lccp.getNumEdges() + " edges in " + runtime_secs + " seconds.");
                                        
                                        if(alsoPrintGraphProps) {
                                            GraphPropResults gpr_lcc = lccp.getGraphPropResults();
                                            System.out.println(GraphProperties.getOverviewPropsString(true, gpr_lcc));
                                        }
                                        
                                        if(Settings.getBoolean("plcc_B_useDB")) {
                                            try {
                                                DBManager.writeCustomgraphStatsToDB(unique_name, description, isForLargestConnectedComponent, num_verts, num_edges, min_degree, max_degree, num_connected_components, diameter, radius, avg_cluster_coeff, avg_shortest_path_length, degreedist, avg_degree, density, cumul_degreedist, runtime_secs);
                                            }catch(SQLException e) {
                                                DP.getInstance().e("Main", "Could not write custom graph stats for largest CC to db: '" + e.getMessage() + "'.");
                                            }
                                        }
                                    }
                                }
                            }
                                                        
                            
                            if( ! skipDrawing) {
                                ProteinGraphDrawer.drawDrawableGraph(outFilePathNoExt, formats, g, vertexMappings);
                                DP.getInstance().i("Main", "Done drawing GML graph to base file '" + outFilePathNoExt + "', exiting.");
                            }                            
                            checkArgsUsage(args, argsUsed);
                            System.exit(0);
                        }
                    }
                    
                    if(s.equals("-l") || s.equals("--draw-plcc-graph")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            System.out.println("Drawing protein graph in plcc format from file '" + args[i+1] + "'.");
                            drawPlccGraphFromFile(args[i+1], args[i+1] + Settings.get("plcc_S_img_output_fileext"), false);
                            System.out.println("Handled plcc graph file '" + args[i+1] + "', exiting.");
                            checkArgsUsage(args, argsUsed);
                            System.exit(0);
                        }
                    }
                    
                    if(s.equals("-b") || s.equals("--draw-plcc-fgs")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            Settings.set("plcc_B_folding_graphs", "true");
                            Settings.set("plcc_B_draw_graphs", "true");                            
                            System.out.println("Drawing protein graph and folding graphs in plcc format from file '" + args[i+1] + "'.");
                            drawPlccGraphFromFile(args[i+1], args[i+1] + Settings.get("plcc_S_img_output_fileext"), true);
                            System.out.println("Handled plcc graph file '" + args[i+1] + " and folding graphs', exiting.");
                            checkArgsUsage(args, argsUsed);
                            System.exit(0);
                        }
                    }
                    
                    if(s.equals("-i") || s.equals("--aa-graphs")) {
                        Settings.set("plcc_B_AAgraph_allchainscombined", "true");
                        Settings.set("plcc_B_AAgraph_perchain", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--aa-graphs-pdb")) {
                        Settings.set("plcc_B_AAgraph_allchainscombined", "true");
                        Settings.set("plcc_B_AAgraph_perchain", "false");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--aa-graphs-chain")) {
                        Settings.set("plcc_B_AAgraph_allchainscombined", "false");
                        Settings.set("plcc_B_AAgraph_perchain", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("-I") || s.equals("--mmCIF-parser")) {
                        Settings.set("plcc_B_use_mmCIF_parser", "true");
                        argsUsed[i] = true;
                    }
                    
                    if(s.equals("--cg-threshold")) {
                        if(args.length <= i+1 ) {
                            syntaxError();
                        }
                        else {
                            argsUsed[i] = true;
                            argsUsed[i+1] = true;
                            Settings.set("plcc_I_cg_contact_threshold", args[i+1]);
                        }
                    }
                    
                    if(s.equals("--chain-spheres-speedup")) {
                        argsUsed[i] = true;
                        Settings.set("plcc_B_chain_spheres_speedup", "true");
                    }
                   
                    if(s.equals("--include-rna"))  {
                        argsUsed[i] = true;
                        Settings.set("plcc_include_rna", "true");
                    }
                    
                    
                    
                    if(s.equals("--matrix-structure-search")) {
                        if(args.length <= i+3 ) {
                            syntaxError();
                        }
                        Settings.set("plcc_S_linear_notation_type", args[i+1]);                        
                        Settings.set("plcc_S_linear_notation", args[i+2]);
                        Settings.set("plcc_S_linear_notation_graph_type", args[i+3]);
                        Settings.set("plcc_B_matrix_structure_search", "true");
                        argsUsed[i] = argsUsed[i+1] = argsUsed[i+2] = argsUsed[i+3] = true;
                                                
                    }
                    
                    if (s.equals("--matrix-structure-search-db")){
                        if (args.length <= i+3){
                            syntaxError();
                        }
                        Settings.set("plcc_S_linear_notation_type", args[i+1]);                        
                        Settings.set("plcc_S_linear_notation", args[i+2]);
                        Settings.set("plcc_S_linear_notation_graph_type", args[i+3]);
                        Settings.set("plcc_B_matrix_structure_search_db", "true");
                        argsUsed[i] = argsUsed[i+1] = argsUsed[i+2] = argsUsed[i+3] = true;
                    }
                    
                    if (s.equals("--settingsfile")) {
                        // as this may overwrite comman line arguments, this sould always preceed command lines
                        if (i > 0) {
                            DP.getInstance().w("Command line argument '--settingsfile' should always be first argument, as it may overwrite settings "
                                    + "from previous command line options.", 2);
                        }
                        
                        if (args.length <= i + 1) {
                            syntaxError();
                        }
                        numSettingsLoaded = Settings.load(args[i + 1]);
                        if (numSettingsLoaded > 0) {
                            outputToBePrintedUnlessSilent.append("  Loaded ").append(numSettingsLoaded).append(" settings as requested by command line argument from properties file ").append(args[i + 1]).append(".\n");
                        }
                        argsUsed[i] = argsUsed[i+1] = true;
                    }
                    
                    

                } //end for loop
                checkArgsUsage(args, argsUsed); // warn if there were extra command line args we do not know
            }

            
        } else {
            usage_short();      // the first argument (pdbid) is required!
            System.exit(1);
        }
        

        
        Boolean silent = false;
        if(Settings.getBoolean("plcc_B_silent")) {
            if(Settings.getBoolean("plcc_B_print_silent_notice")) {
                String startTime = new SimpleDateFormat("HH:mm:ss").format(new Date());
                if(Settings.getBoolean("plcc_B_no_warn")) {
                    System.out.println("[PLCC] [" + pdbid + "] [" + startTime + "] Silent mode and no-warn active, only errors and a final completion notice will be printed from now on. Bye.");
                } else {
                    System.out.println("[PLCC] [" + pdbid + "] [" + startTime + "] Silent mode active, only errors, warnings and a final completion notice will be printed from now on. Bye.");
                }
            }
            silent = true;
        } else {
            System.out.print(outputToBePrintedUnlessSilent.toString());
            if(Settings.getBoolean("plcc_B_no_warn")) {
                System.out.println("  No-warn active, no warnings will be printed.");
            }
        }
        
        if(Settings.getBoolean("plcc_B_use_mmCIF_parser")) {
            if (pdbFile.endsWith(".pdb")) {
                pdbFile = pdbFile.replace(".pdb", ".cif");
            }
            if (! silent) {
                    System.out.println("Using mmCIF parser and therefore looking for .cif file.");
                    System.out.println("Filename now is: " + pdbFile);
            }
        } else {
            // using old parser: everything that does not work with it here
            if (Settings.getBoolean("plcc_include_rna")) {
                DP.getInstance().w("Legacy PDB file parser and inclusion of RNA switched on, but the old parser does not suppoert this setting." +
                        " Use a mmCIF file and command line option '-I' to include RNA. Switching off RNA inclusion now to go on.");
                Settings.set("plcc_include_rna", "false");
            }
        }
        
        if(Settings.getBoolean("plcc_B_clustermode")) {
            if(! silent) {
                System.out.println("Cluster mode active.");
            }
        }
        
        if (Settings.getBoolean("plcc_B_centroid_method")) {
            if (! silent) {
                System.out.println("HINT: Using centroid instead of C_alpha for maxSeqNeighborDist calculation as requested by settings.");
            }
            if (! Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                DP.getInstance().w("Centroid method switched on with chain sphere speedup turned off. The centroid method is optimized for the chain " +
                    "sphere speedup and may produce wrong results without it!");
            }
        }
        
        // This check is rather useless and it will break PDB files that were split into multiple files (one for each
        //  model) and renames, e.g. "2kos_1.pdb" for model 1 of protein 2kos. It is therefore disabaled atm.
        if(pdbid.length() != 4) {
            if(! silent) {
                //System.out.println("WARNING: The PDB identifier '" + pdbid + "' should be 4 characters long (but is " + pdbid.length() + ").");            
            }
        }

        // ****************************************************    test for required files    **********************************************************

        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            print_debug_malfunction_warning();            
        }
        
        if(Settings.getInteger("plcc_I_debug_level") > 0) {
            if(! silent) {
                System.out.println("  Debug level set to " + Settings.getInteger("plcc_I_debug_level") + ".");
            }
        }
        
        if(useFileFromCommandline) {
            if(! silent) {
                if( ! pdbid.toUpperCase().equals("NONE")) {
                    System.out.println("  Using PDB file '" + pdbFile + "', dssp file '" + dsspFile + "', output directory '" + outputDir + "'.");
                }
            }
        }
        
        if(Settings.getInteger("plcc_I_debug_level") > 0) {
            Settings.writeDocumentedDefaultFile(System.getProperty("user.home") + System.getProperty("file.separator") + ".plcc_example_settings");
        }
        
        
        if(Settings.getBoolean("plcc_B_search_similar")) {
            if(! silent) {
                System.out.println("Searching for proteins similar to PDB ID '" + Settings.get("plcc_B_search_similar_PDBID") + "' chain '" + Settings.get("plcc_B_search_similar_chainID") + "' graph type '" + Settings.get("plcc_S_search_similar_graphtype") + "'.");
            }
            
            if(DBManager.init(Settings.get("plcc_S_db_name"), Settings.get("plcc_S_db_host"), Settings.getInteger("plcc_I_db_port"), Settings.get("plcc_S_db_username"), Settings.get("plcc_S_db_password"), Settings.getBoolean("plcc_B_db_use_autocommit"))) {
                
                if(Settings.get("plcc_S_search_similar_method").equals(Similarity.SIMILARITYMETHOD_STRINGSSE)) {
                    if(! silent) {
                        System.out.println("Using similarity method '" + Settings.get("plcc_S_search_similar_method") + "'.");                        
                    }
                
                    String patternSSEString = null;
                    try {
                        patternSSEString = DBManager.getSSEStringOfProteinGraph(Settings.get("plcc_B_search_similar_PDBID"), Settings.get("plcc_B_search_similar_chainID"), Settings.get("plcc_S_search_similar_graphtype"));
                    } catch (Exception e) {
                        System.err.println("ERROR: DB: Could not retrieve SSE string for requested graph from database, exiting.");
                        System.exit(1);
                    }

                    if(patternSSEString == null) {
                        System.err.println("ERROR: DB: SSE string for requested graph is not in the database, exiting.");
                        System.exit(1);
                    } else {
                        if(! silent) {
                            System.out.println("Using pattern SSEstring '" + patternSSEString + "'.");
                        }
                    }

                    CompareOneToDB.performSSEStringComparison(patternSSEString);
                } else if (Settings.get("plcc_S_search_similar_method").equals(Similarity.SIMILARITYMETHOD_GRAPHSET)) {
                    if(! silent) {
                        System.out.println("Using similarity method '" + Settings.get("plcc_S_search_similar_method") + "'.");
                    }
                    CompareOneToDB.performGraphSetComparison();                    
                }
                else if (Settings.get("plcc_S_search_similar_method").equals(Similarity.SIMILARITYMETHOD_GRAPHCOMPAT)) {
                    if(! silent) {
                        System.out.println("Using similarity method '" + Settings.get("plcc_S_search_similar_method") + "'.");
                    }
                    CompareOneToDB.performGraphCompatGraphComparison();                    
                }
                else {
                    System.err.println("ERROR: Invalid similarity method: '" + Settings.get("plcc_S_search_similar_method") + "'. Use --help for info on valid settings.");
                    System.exit(1);
                }
                
                System.exit(0);
            } else {
                System.err.println("ERROR: Could not connect to DB, exiting.");
                System.exit(1);
            }            
            
        }
        
        if(Settings.getBoolean("plcc_B_compute_graphlet_similarities")) {
                        
            int graphletStartIndex = Settings.getInteger("plcc_I_compute_all_graphlet_similarities_start_graphlet_index");
            int graphletEndIndex = Settings.getInteger("plcc_I_compute_all_graphlet_similarities_end_graphlet_index");
            int numToConsider = (graphletEndIndex - graphletStartIndex) + 1;
            
            if(numToConsider < 29) {
                System.err.println("WARNING: Considering only " + numToConsider + " graphlets. Ignore this warning if this is what you want.");
                if(numToConsider <= 0) {
                    System.err.println("ERROR: Configured to consider " + numToConsider + " graphlets in config file. This makes no sense, exiting.");
                    System.exit(1);
                }
            }
            
            if(! silent) {
                System.out.println("Computing pairwise graphlet similarities. Considering " + numToConsider + " graphlets from the array in database (index " + Settings.getInteger("plcc_I_compute_all_graphlet_similarities_start_graphlet_index") + " - " + Settings.getInteger("plcc_I_compute_all_graphlet_similarities_end_graphlet_index") + ").");
            }
            
            if(DBManager.initUsingDefaults()) {
                // protein graphs
                if(Settings.getBoolean("plcc_B_compute_graphlet_similarities_pg")) {
                    if(! silent) {
                        System.out.println("*Computing pairwise graphlet similarities for all protein graphs in the database. This will take a lot of time and memory for large databases...");
                    }
                    Long[] res = DBManager.computeGraphletSimilarityScoresForPGsWholeDatabaseAndStoreBest(ProtGraph.GRAPHTYPE_ALBE, Settings.getInteger("plcc_I_compute_all_graphlet_similarities_num_to_save_in_db"));
                    // numChainsFound, numGraphletsFound, numScoresComputed, numScoresSaved
                    if(! silent) {
                        System.out.println(" PG graphlet similarity done. Found " + res[0] + " PGs (used ALBE graph type) and " + res[1] + " graphlet counts for them in the DB. Computed " + res[2] + " similarity scores and saved " + res[3] + " of them to the DB.");
                    }
                    if(res[0] <= 1) {
                        System.err.println(" WARNING: Database needs to contain at least 2 protein graphs to compute PG graphlet similarity. Run PLCC to add protein graphs to the database (and then GraphletAnalyzer to compute graphlet counts for them).");
                    }
                    if(res[0] >= 2 && res[1] == 0) {
                        System.err.println(" WARNING: Found protein graphs but no graphlet counts for them in the database. (Did you run GraphletAnalyzer on the GML files to compute and save graphlet counts to the database?)");
                    }
                    if( ! DBManager.getAutoCommit()) {
                        DBManager.commit();
                    }
                }
                
                // complex graphs
                if(Settings.getBoolean("plcc_B_compute_graphlet_similarities_cg")) {
                    if(! silent) {
                        System.out.println("*Computing pairwise graphlet similarities for all complex graphs in the database. This will take a lot of time and memory for large databases...");
                    }
                    Long[] res = DBManager.computeGraphletSimilarityScoresForCGsWholeDatabaseAndStoreBest(Settings.getInteger("plcc_I_compute_all_graphlet_similarities_num_to_save_in_db"));
                    // numChainsFound, numGraphletsFound, numScoresComputed, numScoresSaved
                    if(! silent) {
                        System.out.println(" CG graphlet similarity done. Found " + res[0] + " CGs and " + res[1] + " graphlet counts for them in the DB. Computed " + res[2] + " similarity scores and saved " + res[3] + " of them to the DB.");
                    }
                    if(res[0] <= 1) {
                        System.err.println(" WARNING: Database needs to contain at least 2 complex graphs to compute CG graphlet similarity. Run PLCC to add complex graphs to the database (and then GraphletAnalyzer to compute graphlet counts for them).");
                    }
                    if(res[0] >= 2 && res[1] == 0) {
                        System.err.println(" WARNING: Found complex graphs but no graphlet counts for them in the database. (Did you run GraphletAnalyzer on the GML files to compute and save graphlet counts to the database?)");
                    }
                    if( ! DBManager.getAutoCommit()) {
                        DBManager.commit();
                    }
                }
                
                // amino acid graphs
                if(Settings.getBoolean("plcc_B_compute_graphlet_similarities_aag")) {
                    if(! silent) {
                        System.out.println("*Computing pairwise graphlet similarities for all amino acid graphs in the database. This will take a lot of time and memory for large databases...");
                    }
                    Integer numToSave = Settings.getInteger("plcc_I_compute_all_graphlet_similarities_num_to_save_in_db");
                    if(numToSave < 0) { numToSave = null; }
                    Long[] res = DBManager.computeGraphletSimilarityScoresForAAGsWholeDatabaseAndStoreBest(numToSave);
                    // numChainsFound, numGraphletsFound, numScoresComputed, numScoresSaved
                    if(! silent) {
                        System.out.println(" AAG graphlet similarity done. Found " + res[0] + " AAGs and " + res[1] + " graphlet counts for them in the DB. Computed " + res[2] + " similarity scores and saved " + res[3] + " of them to the DB.");
                    }
                    if(res[0] <= 1) {
                        System.err.println(" WARNING: Database needs to contain at least 2 amino acid graphs to compute AAG graphlet similarity. Run PLCC to add amino acid graphs to the database (and then GraphletAnalyzer to compute graphlet counts for them).");
                    }
                    if(res[0] >= 2 && res[1] == 0) {
                        System.err.println(" WARNING: Found amino acid graphs but no graphlet counts for them in the database. (Did you run GraphletAnalyzer on the GML files to compute and save graphlet counts to the database?)");
                    }
                    if( ! DBManager.getAutoCommit()) {
                        DBManager.commit();
                    }
                }
                if(! silent) {
                    System.out.println("All done. Exiting.");
                }
                System.exit(0);
            } else {
                System.err.println("ERROR: Could not connect to DB, exiting.");
                System.exit(1);
            }
            
            
        }
        
        // mark the rep chains in the DB, then exit
        if(Settings.getBoolean("plcc_B_set_pdb_representative_chains_pre") || Settings.getBoolean("plcc_B_set_pdb_representative_chains_post")) {
            
            File xmlFile = new File(Settings.get("plcc_S_representative_chains_xml_file"));
            
            if(! silent) {
                System.out.println("Marking all representative PDB chains in the database from data in XML file '" + Settings.get("plcc_S_representative_chains_xml_file") + "'...");
                System.out.println("  Parsing XML file...");
            }
            
            if( ! (xmlFile.isFile() && xmlFile.canRead())) {
                System.err.println("ERROR: Cannot read XML file '" + xmlFile.getAbsolutePath() + "' or not a normal file.");
                System.exit(1);
            }
            
            // get list by parsing XML
            List<String[]> repChains = new ArrayList<>();
            XMLParserJAX p;
            String[] sep;
            String xml = null;
            try {
                xml = FileParser.slurpFileToString(xmlFile.getAbsolutePath());
            }
            catch(IOException e) {
                System.err.println("ERROR: Failed to read XML file '" + xmlFile.getAbsolutePath() + "', exiting.");
            }
            
            try {
                p = new XMLParserJAX();
                p.setErrorHandler(new XMLErrorHandlerJAX(System.err));
                XMLContentHandlerPDBRepresentatives handler = new XMLContentHandlerPDBRepresentatives();            
                p.handleXML(xml, handler);
                List<String> pdbChains = handler.getPdbChainList();
                System.out.println("Received a list of " + pdbChains.size() + " chains from XML parser. Writing data to DB...");
                for(String ic : pdbChains) {
                    sep = PlccUtilities.parsePdbidAndChain(ic);
                    if(sep != null) {
                        repChains.add(sep);
                        //System.out.println("PDB ID: " + sep[0] + ", chainName " + sep[1] + "");
                    } else {
                        System.err.println("WARNING: Result from XML parsing could not be split into PDB ID and chain, skipping.");
                    }
                }

            } catch(ParserConfigurationException | SAXException | IOException e) {
                System.err.println("ERROR: '" + e.getMessage() + "'. Could not parse XML file, aborting.");
                System.exit(1);
            }
            
            
            // let the DB manager handle the list
            if(DBManager.initUsingDefaults()) {
                Integer[] res;
                Integer numOldLabelsRemoved = 0;

                Integer numChainsInDB = DBManager.countChainsInDB();
                
                // ----- update the existing chains table -----
                if(Settings.getBoolean("plcc_B_set_pdb_representative_chains_post")) {
                    if(numChainsInDB.equals(0)) {
                        System.err.println("WARNING: No chains in the database. (Did you intend to run the PRE version of this command?)");
                    }
                    System.out.println("This should be run after an update, when the database is full. It marks proteins as part of the non-redundant set in the chain table, and entries in this table only exists AFTER you have added data to the database.");
                
                    

                    // remove the old marking if required
                    if(Settings.getBoolean("plcc_B_set_pdb_representative_chains_remove_old_labels_post")) {
                        try {
                            numOldLabelsRemoved = DBManager.markAllChainsAsNonRepresentativeInChainsTable();
                            if( ! DBManager.getAutoCommit()) {
                                DBManager.commit();
                            }
                            // numChainsInList, numChainsUpdatedInDB
                            if(! silent) {
                                System.out.println("  Removed " + numOldLabelsRemoved + " old labels from the " + numChainsInDB + " chains in the DB.");
                            }
                        }
                        catch(SQLException e) {
                            System.err.println("ERROR: Removing old representative labels from chains in DB failed: '" + e.getMessage() + "'.");
                            System.exit(1);
                        }
                    }


                    // mark the chains as representative in the chains table (for already existing chains)
                    try {
                        res = DBManager.markAllRepresentativeExistingChainsInChainsTableFromList(repChains);
                        // numChainsInList, numChainsUpdatedInDB
                        if(! silent) {
                            System.out.println("Done. Found and updated " + res[1] + " of the " + res[0] + " labels from the XML file in the chains table of the DB. (The DB contains " + numChainsInDB + " chains atm.)");
                        }
                    }
                    catch(SQLException e) {
                        System.err.println("ERROR: Updating representatives in chains table of DB failed: '" + e.getMessage() + "'.");
                        System.exit(1);
                    }
                }

                // ----- update the info table -----
                if(Settings.getBoolean("plcc_B_set_pdb_representative_chains_pre")) {
                    if(numChainsInDB > 0) {
                        System.err.println("WARNING: Database already contains " + numChainsInDB + " chains. (Did you intend to run the POST version of this command?)");
                    }
                    System.out.println("This command should be run BEFORE adding data to the database. It is only required if you want PLCC to compute graph statistics for representative chains during the update. You still need to explicitely tell PLCC to do this via command line options.");
                    
                    // remove the old marking in list table if required
                    if(Settings.getBoolean("plcc_B_set_pdb_representative_chains_remove_old_labels_pre")) {
                        try {
                            numOldLabelsRemoved = DBManager.markAllChainsAsNonRepresentativeInInfoTable();
                            if( ! DBManager.getAutoCommit()) {
                                DBManager.commit();
                            }
                            // numChainsInList, numChainsUpdatedInDB
                            if(! silent) {
                                System.out.println("  Removed " + numOldLabelsRemoved + " old labels from the info table in the DB.");
                            }
                        }
                        catch(SQLException e) {
                            System.err.println("ERROR: Removing old representative labels from the info table in DB failed: '" + e.getMessage() + "'.");
                            System.exit(1);
                        }
                    }
                    
                    
                    // mark the chains as representative in the chains table (for already existing chains)
                    try {
                        res = DBManager.markAllRepresentativeExistingChainsInInfoTableFromList(repChains);
                        if(! silent) {
                            System.out.println("Done. Found and updated " + res[1] + " of the " + res[0] + " labels from the XML file in the DB in info table. " + res[2] + " were already in the DB (or duplicated in the XML file).");
                        }
                    }
                    catch(SQLException e) {
                        System.err.println("ERROR: Updating representatives in list table of DB failed: '" + e.getMessage() + "'.");
                        System.exit(1);
                    }
                   
                }
                
                
                if( ! DBManager.getAutoCommit()) {
                    DBManager.commit();
                }
                
                // fill the table of representative chains, only used to compute statistics for these chains WHEN THEY ARE INSERTED LATER
                
                
                System.exit(0);
            } else {
                System.err.println("ERROR: Could not connect to DB, exiting.");
                System.exit(1);
            }   
            
        }
                
        
        if(Settings.getBoolean("plcc_B_report_db_proteins")) {
            String reportFileName = "db_contents_proteins.txt";
            if(! silent) {                
                System.out.println("Reporting list of proteins which are currently in the database to file '" + reportFileName + "'.");
            }
            
            if(DBManager.initUsingDefaults()) {
                List<String> pdbids = DBManager.getAllPDBIDsInTheDB();
                if(IO.stringToTextFile(reportFileName, IO.stringListToString(pdbids, " "))) {
                    if(! silent) {
                        System.out.println("Done. Reported " + pdbids.size() + " proteins in file.");
                    }
                }
                else {
                    System.err.println("ERROR: Writing output file failed.");
                }
                // numChainsFound, numGraphletsFound, numScoresComputed, numScoresSaved
                
                System.exit(0);
            } else {
                System.err.println("ERROR: Could not connect to DB, exiting.");
                System.exit(1);
            }               
        }
        
        // convert pdb file with multiple models to pdb file with multiple chains only
        // the models will be converted to separated chains
        if(Settings.getBoolean("plcc_B_convert_models_to_chains")) {
                
                System.out.println("  Converting models of input PDB file '" + convertModelsToChainsInputFile + "' to chains and storing new PDB file at '" + convertModelsToChainsOutputFile + "'...");
                Boolean status = FileParser.convertPdbModelsToChains(convertModelsToChainsInputFile, convertModelsToChainsOutputFile);
                if(status) {
                    System.out.println("  Conversion done. The resulting new PDB file is at '" + convertModelsToChainsOutputFile + "'. Note that only structural data (atom positions) have been converted, HEADER information has not been adapted.");
                    System.exit(0);
                }
                else {
                    DP.getInstance().w("Main", "Conversion of PDB file '" + convertModelsToChainsInputFile + "' failed.");
                    System.exit(1);
                }
       }
        
        
        // **************************************    Protein structure search in the database    ******************************************
        
        if(Settings.getBoolean("plcc_B_matrix_structure_search_db")) {
            String linnotGraphType = Settings.get("plcc_S_linear_notation_graph_type");

            // check input parameters first
            String viableLinnotGraphTypes[] = new String[] {"alpha", "beta", "albe"};
            if (! Arrays.asList(viableLinnotGraphTypes).contains(linnotGraphType)) {
                DP.getInstance().e("Main", "Unrecognized linear notation graph type '" + linnotGraphType +
                        "'. Allowed values: " + Arrays.toString(viableLinnotGraphTypes) + " Exiting now.");
                System.exit(1);
            }
            // NOTE linnot type not checked as currently only red used 
            
            if(DBManager.initUsingDefaults()) {
                System.out.println("Start searching the linear notation " + Settings.get("plcc_S_linear_notation") + " in the PTGL database.");

                ArrayList<ArrayList<String>> results = new ArrayList<ArrayList<String>>();

                results = DBManager.matrixSearchDB(Settings.get("plcc_S_linear_notation"), linnotGraphType);

                int count_results = results.size();
                System.out.println("  The structure was found in " + count_results + " protein chains.");
                if (count_results > 0){

                    //writing results into a text file
                    try {
                        System.out.println("Saving all proteins (pdbid and chain) in the file 'matrix_search_db_results.lst'. ");

                        File file_results = new File("matrix_search_db_results.lst");
                        FileWriter writer = new FileWriter(file_results);

                        for (ArrayList<String> r : results){ //print the results = pdbid and chain of all proteins, that contain the linear notation from the input
                            writer.write(r.get(0) + r.get(1)); //write results to file
                            writer.write(System.getProperty("line.separator")); //add a new line
                        }
                        writer.flush();
                        writer.close();

                    } catch (IOException ex) {
                        java.util.logging.Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                System.out.println("Exiting now");

                System.exit(0);
            } else {
                System.err.println("ERROR: Could not connect to DB, exiting.");
                System.exit(1);
            }
        }
        

        // ************************************** check database connection ******************************************        

        if(Settings.getBoolean("plcc_B_useDB")) {
            String plcc_db_name = Settings.get("plcc_S_db_name");
            String plcc_db_host = Settings.get("plcc_S_db_host");
            Integer plcc_db_port = Settings.getInteger("plcc_I_db_port");
            String plcc_db_username = Settings.get("plcc_S_db_username");
            String plcc_db_password = Settings.get("plcc_S_db_password");
            if(! silent) {
                System.out.println("  Checking database connection to host '" + plcc_db_host + "' on port '" + plcc_db_port + "'...");
            }
            if(DBManager.initUsingDefaults()) {
                if(! silent) {
                    System.out.println("  -> Database connection OK.");
                }
            }
            else {
                System.out.println("  -> Database connection FAILED.");
                DP.getInstance().w("Could not establish database connection, not writing anything to the DB.");
                Settings.set("plcc_B_useDB", "false");
            }
        }
        else {
            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.println("  Not using the database as requested by options.");
            }
        }
        

        // **************************************    file checks    ******************************************
        File input_file;
        File output_dir;
        
        // pdb file
        input_file = new File(pdbFile);

        if(! (input_file.exists() && input_file.isFile())) {
            System.err.println("ERROR: pdbfile '" + pdbFile + "' not found. Exiting.");
            System.exit(1);
        }
                

        // dssp file
        input_file = new File(dsspFile);
        if(! (input_file.exists() && input_file.isFile())) {
            System.err.println("ERROR: dsspfile '" + dsspFile + "' not found. Exiting.");
            System.exit(1);
        }

        output_dir = new File(outputDir);
        if(! (output_dir.exists() && output_dir.isDirectory())) {
            System.err.println("ERROR: output directory '" + outputDir + "' not found. Exiting.");
            System.exit(1);
        }


        if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
            System.out.println("  Checked required files and directories, looks good.");
        }

        pdbIdDotGeoFile = output_dir + fs + pdbid.toLowerCase() + ".geo";               // holds info on contacts between residues of the PDB file
        pdbIdDotGeoLigFile = output_dir + fs + pdbid.toLowerCase() + ".geolig";         // holds info on contacts between residues + ligands of the PDB file
        conDotSetFile = output_dir + fs + pdbid.toLowerCase() + ".contactstats";        // holds statistics on atom contacts by residue type
        dsspLigFile = output_dir + fs + pdbid.toLowerCase() + ".dssplig";
        chainsFile = output_dir + fs + pdbid.toLowerCase() + ".chains";
        ligandsFile = output_dir + fs + pdbid.toLowerCase() + ".ligands";
        modelsFile = output_dir + fs + pdbid.toLowerCase() + ".models";
        resMapFile = output_dir + fs + pdbid.toLowerCase();         // chainName and file extension is added later 

        if(Settings.getBoolean("plcc_B_ptgl_text_output")) {
            if(! silent) {
                System.out.println("  Using output files:\n    * " + pdbIdDotGeoFile + " for contact data\n    * " + pdbIdDotGeoLigFile + " for lig contact data");
                System.out.println("    * " + conDotSetFile + " for contact statistics\n    * " + dsspLigFile + " for DSSP ligand file.");
            }
        }


        // **************************************    here we go: parse files and get data    ******************************************
        if(! silent && ! Settings.getBoolean("plcc_B_only_essential_output")) {
            System.out.println("Getting data...");
        }
        
        FileParser.initData(pdbFile, dsspFile);
         
        if (Settings.getBoolean("plcc_B_debug_only_parse")) {
            System.out.println("Exiting now as requested by settings.");
            System.exit(0);
        }
        
        allModelsIDsOfWholePDBFile = FileParser.getAllModelIDsFromWholePdbFile();

        if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
            // print all model IDs from the PDB file (not just the handled model)
            System.out.println("    PDB: Found the following NMR models in the whole PDB file:");
            System.out.print("    PDB:");

            if(allModelsIDsOfWholePDBFile.size() > 0) {

                for(Integer i = 0; i < allModelsIDsOfWholePDBFile.size(); i++) {
                    System.out.print("  " + allModelsIDsOfWholePDBFile.get(i));
                }
                System.out.print("\n");
            }
            else {
                System.out.print(" <None> (No NMR data, an artificial default model was added.)\n");
            }
        }
        
        models = FileParser.getModels();    // doesn't do much anymore since only the PDB lines of model 1 are currently in there
        chains = FileParser.getChains();       
        molecules = FileParser.getMolecule();
        atoms = FileParser.getAtoms();
        sulfurBridges = FileParser.getSulfurBridges();
        
        if(Settings.getBoolean("plcc_B_skip_too_large")) {
            if(atoms.size() > Settings.getInteger("plcc_I_skip_num_atoms_threshold")) {
                System.err.println("===== Terminating: 'plcc_B_skip_too_large' is enabled in settings and this protein has too many atoms =====");
                if(! silent) {
                    System.out.println("Maximal number of atoms allowed in PDB file is set to " + Settings.getInteger("plcc_I_skip_num_atoms_threshold") + " and PDB " + pdbid + " contains " + atoms.size() + ".");
                    System.out.println("This is a batch/cluster feature to ignore very large PDB files which take ages to compute.");
                    System.out.println("A max_atoms setting of around 80.000 makes most sense (see PDB statistics).");
                    System.out.println("Set 'plcc_B_skip_too_large' to false in settings to avoid this. Exiting.");
                }
                System.exit(0);
            }
        }
        
        // check whether we need to abort processing because of bad resolution or too few residues
        Float badResolution = Settings.getFloat("plcc_F_abort_if_pdb_resolution_worse_than");
        if(badResolution >= 0.0) {
            HashMap<String, String> md;
            md = FileParser.getMetaData();
            
            Double resolution = -1.0;
            try {
                resolution = Double.valueOf(md.get("resolution"));
            } catch (Exception e) {
                resolution = -1.0;
                DP.getInstance().w("Could not determine resolution of PDB file for protein '" + pdbid + "', assuming NMR with resolution '" + resolution + "'.");            
            }
            
            if(resolution > badResolution) {
                DP.getInstance().e("Main", "Aborting further processing of PDB '" + pdbid + "': resolution '" + resolution + "' too bad, must be at '" + badResolution + "' or better. (Set 'plcc_F_abort_if_pdb_resolution_worse_than' to a negative float value in the config file to prevent this behaviour or use --force.) Exiting now.");
                System.exit(0);
            }
        }
        
        Integer minNumberOfResidues = Settings.getInteger("plcc_I_abort_if_num_molecules_below");
        if(molecules.size() < minNumberOfResidues) {
            DP.getInstance().e("Main", "Aborting further processing of PDB '" + pdbid + "': molecule count '" + molecules.size() + "' too low, must be at least '" + minNumberOfResidues + "'. (Set 'plcc_I_abort_if_num_molecules_below' to a negative int value in the config file to prevent this behaviour or use --force.) Exiting now.");
            System.exit(0);
        }
        
        
        
        String sBondString;
        if(sulfurBridges.size() > 0) {
            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.print("    DSSP: Protein contains " + sulfurBridges.size() + " disulfide bridges: ");
            
                for(Character key : sulfurBridges.keySet()) {
                    sBondString = "(";
                    ArrayList<Integer> dsspResidueIDs = sulfurBridges.get(key);

                    if(dsspResidueIDs.size() != 2) {
                        DP.getInstance().w("Disulfide bridge has " + dsspResidueIDs.size() + " member residues.");
                    }

                    for(Integer dsspID : dsspResidueIDs) {
                        sBondString += (dsspID + " ");
                    }
                    sBondString += ")";
                    System.out.print(sBondString);
                }
                
                System.out.print("\n"); 
            }
        } else {
            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.println("    DSSP: Protein contains no disulfide bridges.");
            }
        }

        // This is now done, separate for each chainName, by a function in Main.java below
        // dsspSSEs = FileParser.getDsspSSEs(); 
        // ptglSSEs = FileParser.getPtglSSEs();
        

        if(chains.size() < 1) { System.out.println("WARNING: Input files contain no chains."); }
        if(molecules.size() < 1) { System.out.println("WARNING: Input files contain no molecules."); }
        if(atoms.size() < 1) { System.out.println("WARNING: Input files contain no atoms."); }



        // DEBUG: print all SSEs
        //System.out.println("Printing all " + SSEs.size() + " SSEs according to DSSP definition:");
        //for(Integer i = 0; i < SSEs.size(); i++) {
        //    System.out.println("  " + SSEs.get(i));
        //}
        
        // If only one chain is present disable chain spheres speedup to avoid
        // costly pre processing
        //  --> really? Since seq neigh skip it may also be faster for single chains (see speedtest)
        //      --> yeah, still seems so
        if (Settings.getBoolean("plcc_B_chain_spheres_speedup") && chains.size() == 1) {
            Settings.set("plcc_B_chain_spheres_speedup", "false");
            Settings.set("plcc_B_centroid_method", "false");  // not optimized for old contact computation
            if (! silent) {
                System.out.println("  Note: Chain spheres speedup was turned on in settings, but only one chain was detected. " +
                    "To save time, setting was turned off for this structure.");
            }
        }

        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            // Use fake values to disable residue skipping and prevent contact calculations before the
            //  residue level contact function is called for all residues (less cluttered debug output).
            globalMaxCenterSphereRadius = 1000;
            globalMaxSeqNeighborResDist = 1000;
        }
        else {
            
            if (Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                // does currently not use Res skipping, so no preprocessing needed
                // set the values to be sure that everything works properly
                globalMaxCenterSphereRadius = Integer.MAX_VALUE;
                globalMaxSeqNeighborResDist = Integer.MAX_VALUE; 
            } else {
                // All residues exist, we can now calculate their maximal center sphere radius
                globalMaxCenterSphereRadius = getGlobalMaxCenterSphereRadius(molecules);
                if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                    System.out.println("  Maximal center sphere radius for all residues is " + globalMaxCenterSphereRadius + ".");
                }

                // ... and the maximal distance between neighbors in the AA sequence.
                // Note that this is a lot less useful with ligands enabled since they are always listed at the 
                //  end of the chainName and may be far (in 3D) from their predecessor in the sequence.
                globalMaxSeqNeighborResDist = getGlobalMaxSeqNeighborResDist(molecules);     
                if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                    System.out.println("  Maximal distance between residues that are sequence neighbors is " + globalMaxSeqNeighborResDist + ".");
                }
            }        
        }
        // ... and fill in the frequencies of all AAs in this protein.
        getAADistribution(resFromMolecules(molecules));

        if(! silent) {
            if (! Settings.getBoolean("plcc_include_rna")) {
                // RNA off
                System.out.println("Received all data (" + models.size() + " Models, " + chains.size() + " Chains, " + molecules.size() + 
                        " Residues, " + atoms.size() + " Atoms).");
            } else {
                // RNA on
                System.out.println("Received all data (" + models.size() + " Models, " + chains.size() + " Chains, " + molecules.size() + 
                        " Molecules (" + resFromMolecules(molecules).size() + " Residues and " + rnaFromMolecules(molecules).size() +
                        " RNAs), " + atoms.size() + " Atoms).");
            }
        }


        // **************************************    calculate the contacts    ******************************************

        // JSON test 
        /*
        Boolean jsonTest = false;
        if(jsonTest) {
            Gson gson = new Gson();
            gson.toJson(1);
            gson.toJson("abcd");
            gson.toJson(new Long(10));
            int[] values = { 1 };
            gson.toJson(values);

            //(Deserialization)
            int one = gson.fromJson("1", int.class);
            Integer one1 = gson.fromJson("1", Integer.class);
            Long one2 = gson.fromJson("1", Long.class);
            Boolean bfalse = gson.fromJson("false", Boolean.class);
            String str = gson.fromJson("\"abc\"", String.class);
            String anotherStr = gson.fromJson("[\"abc\"]", String.class);
        }
        */

        if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
            System.out.println("Calculating residue contacts...");
        }

        
        boolean separateContactsByChain = Settings.getBoolean("plcc_B_separate_contacts_by_chain");
        
        if(separateContactsByChain) {
            if(! silent) {
                System.out.println("Separating atom contact computation by chains. Will not compute inter-chain contacts.");
            }
        }         
        else {
            if( ! Settings.getBoolean("plcc_B_complex_graphs")) {
                if((chains.size() > 1 && atoms.size() > 50000) || (chains.size() > 2 && atoms.size() > 15000)) {
                    if(! silent) {
                        System.out.println("INFO: This multi-chain protein is large (" + atoms.size() + " atoms, " + chains.size() + " chains).");
                        System.out.println("INFO:  Using chain separation (the '-E' command line switch) will speed up the computation a lot.");
                        if (! Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                            System.out.println("INFO:  Turning on chain_spheres_speedup may shorten the runtime significantly.");
                        }
                    }
                }
            }
        }
        
        ArrayList<MolContactInfo> cInfoThisChain;
        ProteinResults.getInstance().setPdbid(pdbid);
        
        if(separateContactsByChain) {
            cInfoThisChain = new ArrayList<MolContactInfo>();   // will be computed separately for each chainName later
            cInfo = null;                                       // will not be used in this case (separateContactsByChain=on)
        } else {        
            if(Settings.getBoolean("plcc_B_alternate_aminoacid_contact_model") || Settings.getBoolean("plcc_B_alternate_aminoacid_contact_model_with_ligands")) {
                cInfo = calculateAllContactsAlternativeModel(resFromMolecules(molecules));
            }
            else {
                if (Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                    cInfo = calculateAllContactsChainSphereSpeedup(chains);
                } else {
                    cInfo = calculateAllContacts(molecules);
                }
            }
            if(! silent) {
                System.out.println("Received data on " + cInfo.size() + " residue contacts that have been confirmed on atom level.");
            }
            cInfoThisChain = null;
        }
               
        if (Settings.getBoolean("plcc_B_debug_only_contact_comp")) {
            System.out.println("Exiting now as requested by settings.");
            System.exit(0);
        }
        
        // DEBUG: compare computed contacts with those from a geom_neo file
        if(compareResContacts) {
            if(separateContactsByChain) {
                DP.getInstance().w("Cannot compare residue contacts for whole PDB file in compute-contacts-per-chain mode.");
            }
            else {
                if(! silent) {
                    System.out.println("DEBUG: Comparing calculated residue contacts with those in the file '" + compareResContactsFile + "'.");
                }
                FileParser.compareResContactsWithPdbidDotGeoFile(compareResContactsFile, false, cInfo);
            }
        }               
        
        
        
        
        // **************************************    output of the results of atom/residue level contact computation   ******************************************
        
        // print overview to STDOUT
        if(Settings.getBoolean("plcc_B_print_contacts")) {
            if(separateContactsByChain) {
                DP.getInstance().w("Cannot show contact overview for whole PDB file in compute-contacts-per-chain mode.");
            } 
            else {
                if(! silent) {
                    System.out.println("Showing contact overview...");
                }
                showContactOverview(cInfo);
            }
        }

        if(Settings.getBoolean("plcc_B_print_contacts")) {
            if(separateContactsByChain) {
                DP.getInstance().w("Cannot show contact statistics for whole PDB file in compute-contacts-per-chain mode.");
            }
            else {
                if(! silent) {
                    System.out.println("Showing statistics overview of residues...");
                }
                showContactStatistics(resFromMolecules(molecules).size());
            }
            
        }
        
        // write the detailed and formated results to the output file
        if(Settings.getBoolean("plcc_B_ptgl_text_output")) {
            
            if(separateContactsByChain) {
                DP.getInstance().w("Cannot write residue contact info file for whole PDB file in compute-contacts-per-chain mode.");
            }
            else {
                if(! silent) {
                    System.out.println("Writing residue contact info file...");
                }
                writeContacts(cInfo, pdbIdDotGeoFile, false);
            }

            if(Settings.getBoolean("plcc_B_write_lig_geolig")) {
                if(separateContactsByChain) {
                    DP.getInstance().w("Cannot show contact statistics for whole PDB file in compute-contacts-per-chain mode.");
                }   
                else {
                    if(! silent) {
                        System.out.println("Writing full contact info file including ligands...");
                    }
                    writeContacts(cInfo, pdbIdDotGeoLigFile, true);
                }
            }

            if(separateContactsByChain) {
                DP.getInstance().w("Cannot write statistics file for whole PDB file in compute-contacts-per-chain mode.");
            }
            else {
                if(! silent) {
                    System.out.println("Writing statistics file for residues...");
                }
                writeStatistics(conDotSetFile, pdbid, resFromMolecules(molecules).size());
            }

            // write the dssplig file
            if(! silent) {
                System.out.println("Writing DSSP ligand file for residues...");
            }
            //writeDsspLigFile(dsspFile, dsspLigFile, cInfo, residues);
            writeOrderedDsspLigFile(dsspFile, dsspLigFile, resFromMolecules(molecules));

            // write chains file
            if(! silent) {
                System.out.println("Writing chain file...");
            }
            writeChains(chainsFile, pdbid, chains);

            // write ligand file
            if(! silent) {
                System.out.println("Writing ligand file for residues...");
            }
            writeLigands(ligandsFile, pdbid, resFromMolecules(molecules));
            
            // write residue mapping files
            if(! silent) {
                System.out.println("Writing residue mapping files for all chains...");
            }
            for(Chain c : chains) {
                writeResMappings(resMapFile + "_" + c.getPdbChainID() + ".resmap", c);
            }
                                    

            // write models file
            if(! silent) {
                System.out.println("Writing models file...");
            }
            writeModels(modelsFile, pdbid, allModelsIDsOfWholePDBFile);

            // generate the PyMol selection script
            if(separateContactsByChain) {
                DP.getInstance().w("Cannot write PyMol script for whole PDB file in compute-contacts-per-chain mode.");
            }
            else {
                if(! silent) {
                    System.out.println("Generating Pymol script to highlight protein-ligand contacts...");
                }
                // String pms = getPymolSelectionScript(cInfo);             // old version: all protein residues on 1 selection
                String pms = getPymolSelectionScriptByLigand(cInfo);        // new version: a separate selection of protein residues for each ligand


                //System.out.println("***** Pymol script follows: *****");
                //System.out.print(pms);
                //System.out.println("***** End of Pymol script. *****");

                String pmsFile = outputDir + fs + pdbid + ".pymol";

                if(writeStringToFile(pmsFile, pms)) {
                    if(! silent) {
                        System.out.println("  PyMol script written to file '" + pmsFile + "'.");
                    }
                }
                else {
                    System.err.println("ERROR: Could not write PyMol script to file '" + pmsFile + "'.");
                }
            }

        }
        else {
            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.println("  Not writing any interim results to text files as requested (geom_neo compatibility mode off).");
            }
        }
        
        if(Settings.getBoolean("plcc_B_write_chains_file")) {
            // write the chains file if it has not yet been written. Used in cluster mode for the GraphletAnlayzer software to know the names of all graphs of all chains.
            if( ! Settings.getBoolean("plcc_B_ptgl_text_output")) {
                if(! silent) {
                    System.out.println("Writing chain file...");
                }
                writeChains(chainsFile, pdbid, chains);
            }
        }
        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            print_debug_malfunction_warning();
            System.out.println("WARNING: ABORTING execution here due to DEBUG settings, results incomplete.");
            System.exit(1);
        }
        
        ArrayList<Chain> handleChains = new ArrayList<Chain>();
        if(Settings.getBoolean("plcc_B_force_chain")) {
        
            if(! silent) {
                System.out.println(" Forced handling of the chain with chain ID '" + Settings.get("plcc_S_forced_chain_id") + "' only.");
            }
            
            for(Chain c : chains) {
                if(c.getPdbChainID().equals(Settings.get("plcc_S_forced_chain_id"))) {
                    handleChains.add(c);
                }
            }
            
        }
        else {
            handleChains = chains;
        }
        
        if(handleChains.size() < 1) {
            DP.getInstance().w("No chains to handle found in input data (" + chains.size() + " chains total).");
            Main.doExit(1);
        }
                
        

        // *************************************** ramachandran plots for chains  *********************************//
        
        Boolean drawRPlots = Settings.getBoolean("plcc_B_ramachandran_plot");
        String plotPath, label;
        if(drawRPlots) {
            for(Chain c : handleChains) {
                plotPath = outputDir + fs + pdbid + "_" + c.getPdbChainID() + "_plot";
                label = "Ramachandran plot of PDB entry " + pdbid + ", chain " + c.getPdbChainID() + "";
                drawRamachandranPlot(plotPath, c.getResidues(), label);
            }
        }
                        
        

        // ********************************************** SSE stuff and graph computation ******************************************//

        // sanity check for settings
        if(Settings.getInteger("plcc_I_aag_min_residue_seq_distance_for_contact") > 0 && Settings.getInteger("plcc_I_aag_max_residue_seq_distance_for_contact") > 0) {
            DP.getInstance().w("Main", "Settings 'plcc_I_aag_min_residue_seq_distance_for_contact' and 'plcc_I_aag_max_residue_seq_distance_for_contact' should NOT be used together: inter-chain contacts cannot pass both checks. Set one of them to 0 to disable it.");
        }
        
        if(Settings.getBoolean("plcc_B_calc_draw_graphs")) {
            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output") )) {
                System.out.println("Calculating SSE graphs.");
                System.out.println("Calculating SSEs for all chains of protein " + pdbid + "...");
            }
            
            if(Settings.getBoolean("plcc_B_useDB")) {
                writeProteinDataToDatabase(pdbid, residues.size());
                if( ! DBManager.getAutoCommit()) {
                    DBManager.commit();
                }
            }
            
            if(Settings.getBoolean("plcc_B_AAgraph_allchainscombined")) {
                if(separateContactsByChain || cInfo == null) {
                    System.err.println("Cannot compute amino acid level contact graph for all chains combined, contact separation is on.");
                } else {
                    if(! silent) {
                        System.out.println("Computing amino acid level contact graph for all chains combined.");
                    }
                    
                    
                    String subDirTree = "";
                    if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                        subDirTree = IO.createSubDirTreeDir(outputDir, pdbid, "ALL");
                        if(subDirTree == null) { 
                            DP.getInstance().e("Main", "Could not create subdir tree (outputDir='" + outputDir + "', pdbid='" + pdbid + "'). Missing file system level access rights?"); 
                            System.exit(1); 
                        }
                    }
                    
                    if (Settings.getBoolean("plcc_B_alternate_aminoacid_contact_model") || Settings.getBoolean("plcc_B_alternate_aminoacid_contact_model_with_ligands")) {

                        PPIGraph ppig;
                        ppig = new PPIGraph(residues, cInfo);
                        ppig.setPdbid(pdbid);
                        ppig.setChainid(AAGraph.CHAINID_ALL_CHAINS);
                        // write the PPI graph to disc
                        String ppigFile = outputDir + fs + subDirTree + pdbid + "_aagraph.gml";
                        if (writeStringToFile(ppigFile, ppig.toGraphModellingLanguageFormat())) {
                            if (!silent) {
                                System.out.println("  PPIGraph for all chains written to file '" + ppigFile + "'.");
                            }
                        } else {
                            System.err.println("ERROR: Could not write PPIGraph for all chains to file '" + ppigFile + "'.");
                        }

                        // write the simple PPI graph to disc
                        String simplePpigFile = outputDir + fs + subDirTree + pdbid + "_aagraph_simple.fanmod";
                        if (writeStringToFile(simplePpigFile, ppig.toFanMod().get(0))) {
                            if (!silent) {
                                System.out.println("  Simple PPIGraph for all chains written to file '" + simplePpigFile + "'.");
                            }
                        } else {
                            System.err.println("ERROR: Could not write simple PPIGraph for all chains to file '" + simplePpigFile + "'.");
                        }
                        // Write corresponding index file to disc
                        String simplePpigIndexFile = outputDir + fs + subDirTree + pdbid + "_aagraph_simple.id";
                        if (writeStringToFile(simplePpigIndexFile, ppig.toFanMod().get(1))) {
                            if (!silent) {
                                System.out.println("  Simple PPIGraph index for all chains written to file '" + simplePpigIndexFile + "'.");
                            }
                        } else {
                            System.err.println("ERROR: Could not write simple PPIGraph index for all chains to file '" + simplePpigIndexFile + "'.");
                        }

                        writePPIstatistics(cInfo, pdbid);

                        // Writes and saves a python script that can be used to visualize bonds with PyMol
                        if (getPymolSelectionScriptPPI(cInfo, pdbid)) {
                            System.out.println("[PYMOL] Python script successfully written.");
                        } else {
                            System.out.println("[PYMOL] Error: Python script could not be written.");
                        }

                    } else {
                    
                    //DEBUG
                    //for(Residue r : residues) {
                    //    if(r.isOtherRes()) {System.out.println("##########OTHER! type="+ r.getType() ); }
                    //    else { System.out.println("###########SSE: " + r.getNonEmptySSEString() + " type=" + r.getType()); }
                    //}
                    residuesWithoutLigands = new ArrayList<>();
                    AAGraph aag;
                    
                    Boolean skipLigandsForAAGraphs = ( ! Settings.getBoolean("plcc_B_aminoacidgraphs_include_ligands"));
                    if(skipLigandsForAAGraphs) {
                        for(Residue r : residues) {
                            if(r.isAA()) { residuesWithoutLigands.add(r); }
                        }
                        aag = new AAGraph(residuesWithoutLigands, cInfo);
                    }
                    else {
                        aag = new AAGraph(residues, cInfo);
                    }
                    aag.setPdbid(pdbid);
                    aag.setChainid(AAGraph.CHAINID_ALL_CHAINS);
                    //String gml = aag.toGraphModellingLanguageFormat();
                    //IO.stringToTextFile("graph.gml", gml);
                    
                    

                    if(Settings.getBoolean("plcc_B_useDB")) {
                        try {
                            if(Settings.getBoolean("plcc_B_write_graphstrings_to_database_aag")) {
                                DBManager.writeAminoAcidGraphToDB(pdbid, AAGraph.CHAINID_ALL_CHAINS, aag.toGraphModellingLanguageFormat(), aag.getNumVertices(), aag.getNumEdges());
                            }
                            else {
                                DBManager.writeAminoAcidGraphToDB(pdbid, AAGraph.CHAINID_ALL_CHAINS, null, aag.getNumVertices(), aag.getNumEdges());
                            }
                            
                            if(! silent) {
                                System.out.println("Wrote all chains AA graph of " + pdbid + " to DB.");
                            }
                        } catch(SQLException e) {
                            DP.getInstance().w("Main", "Could not write all chains AA graph to DB: '" + e.getMessage() + "'.");
                        }
                    } 

                    
                    if(Settings.getBoolean("plcc_B_compute_graph_metrics")) {    
                        //aag.selfCheck();
                        GraphProperties gp = new GraphProperties(aag);
                        SparseGraph lcc = (SparseGraph)gp.getLargestConnectedComponent();
                        GraphProperties gp_lcc = new GraphProperties(lcc);
                        Date propsComputationEndTime; Date propsComputationStartTime; Long timeDiff; Long runtime_secs;

                        if(Settings.getBoolean("plcc_B_useDB")) {

                            if( ! DBManager.getAutoCommit()) {
                                DBManager.commit();
                            }

                            try {
                                Long graph_db_id = DBManager.getDBAminoacidgraphID(pdbid);
                                if(graph_db_id > 0L) {
                                    //System.out.println("Found aa graph " + pdbid + ".");
                                    // write graph properties
                                    propsComputationStartTime = new Date();                                    
                                    GraphPropResults gpr = gp.getGraphPropResults();                                    
                                    propsComputationEndTime = new Date();
                                    timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                    runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);                                    
                                    DBManager.writeAminoacidgraphStatsToDB(graph_db_id, Boolean.FALSE, gpr.numVertices, gpr.numEdges, gpr.minDegree, gpr.maxDegree, gpr.numConnectedComponents, gpr.graphDiameter, gpr.graphRadius, gpr.averageClusterCoefficient, gpr.averageShortestPathLength, gp.getDegreeDistributionUpTo(50), gpr.averageDegree, gpr.density, gp.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);
                                    
                                    Boolean rewireGraphsToCompareWithRandom = true;
                                    Boolean writeRandomGraphPropsToDatabase = true;
                                    Double edgeRewireProbability = 1.0d;
                                    String rand_graph_name = pdbid + "_AA_random";  // same for graph and lcc
                                    
                                    if(rewireGraphsToCompareWithRandom) {   // this was only needed for graph analyses carried out for the Ph.D. thesis of TS, ignore                                                                                
                                        //System.out.println("###TEST-AAG-GP-BEFORE-REWIRING: \n" + GraphProperties.getOverviewPropsString(true, gpr) + "###");
                                        //aag.selfCheck();
                                        
                                        GraphRandomizer gr = new GraphRandomizer(aag, edgeRewireProbability);
                                        GraphProperties gp_rand = new GraphProperties(aag); // now changed
                                        propsComputationStartTime = new Date();
                                        GraphPropResults gpr_rand = gp_rand.getGraphPropResults();
                                        propsComputationEndTime = new Date();
                                        timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                        runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);
                                        //System.out.println("###TEST-AAG-GP-AFTER-RANDOMIZATION with p="+edgeRewireProbability+": \n" + GraphProperties.getOverviewPropsString(false, gpr_rand) + "###");
                                        //aag.selfCheck();
                                        
                                        
                                        
                                        if(writeRandomGraphPropsToDatabase) {                                            
                                            String rand_graph_description = "full random graph";
                                            DBManager.deleteCustomGraphStatsFromDBByUniqueName(rand_graph_name);
                                            DBManager.writeCustomgraphStatsToDB(rand_graph_name, rand_graph_description, Boolean.FALSE, gpr_rand.numVertices, gpr_rand.numEdges, gpr_rand.minDegree, gpr_rand.maxDegree, gpr_rand.numConnectedComponents, gpr_rand.graphDiameter, gpr_rand.graphRadius, gpr_rand.averageClusterCoefficient, gpr_rand.averageShortestPathLength, gp_rand.getDegreeDistributionUpTo(50), gpr_rand.averageDegree, gpr_rand.density, gp_rand.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);                                            
                                        }
                                    }
                                    
                                    Boolean selectResiduesWithPropertiesForVisualization = false;
                                    if(selectResiduesWithPropertiesForVisualization) { //determine and print max ecc vertex set, for Ph.D. thesis only (output can be used to mark the residues in PyMol)
                                        
                                        //Set<Integer> s = gp.determineMaxEccVertexSet(); System.out.println("Max ECC vertices: " + IO.setIntegerToString(s) + ".");
                                        //Set<Integer> s = gp.determineVertexSetWithEccAtLeast(12);
                                        //Set<Integer> s = gp.determineVertexSetWithClCAtLeast(0.8);
                                        //Set<Integer> s = gp.determineVertexSetWithClCAtMost(0.25);
                                        //Set<Integer> s = gp.determineVertexSetWithDegreeAtLeast(14);
                                        
                                        Set<Integer> s = aag.findResidues("A", 153, null);
                                        
                                        List<Residue> l = aag.getResiduesFromSetByIndex(s);
                                        System.out.println("Pymol select script for " + l.size() + " of the " + aag.getSize() + " residues: '" + Main.getPymolSelectionScriptForResidues(l) + "'.");
                                    }
                                    
                                    
                                    Boolean showResProps = false;
                                    if(showResProps) {
                                        Map<String, Integer> maxDegrees = aag.getMaxDegreeByAminoacidType();
                                        Map<String, Integer> maxDiams = aag.getMaxFerretDiameterByAminoacidType();
                                        Map<String, Double> avgDiams = aag.getAverageFerretDiameterByAminoacidType();
                                        Map<String, Double> avgDegrees = aag.getAverageDegreeByAminoacidType();

                                        List<String> sortedKeys = new ArrayList<>(maxDegrees.keySet());
                                        Collections.sort(sortedKeys);
                                        StringBuilder sb = new StringBuilder();
                                        sb.append("Residue types: \n");
                                        for(String key : sortedKeys) {
                                            sb.append(key + "\t");
                                        }
                                        sb.append("\n");
                                        sb.append("Atom count by residue type: \n");
                                        Integer[] numAtoms = new Integer[20];
                                        int tmp = 0;
                                        for(String key : sortedKeys) {
                                            sb.append(AminoAcid.atomCountOfName1(key) + "\t");
                                            numAtoms[tmp] = AminoAcid.atomCountOfName1(key);
                                            tmp++;
                                        }
                                        sb.append("\n");
                                        sb.append("maxDegree by residue type: \n");
                                        for(String key : sortedKeys) {
                                            sb.append(maxDegrees.get(key) + "\t");
                                        }
                                        sb.append("\n");
                                        sb.append("maxDiam by residue type: \n");
                                        for(String key : sortedKeys) {
                                            sb.append(maxDiams.get(key) + "\t");
                                        }
                                        sb.append("\n");
                                        System.out.println(sb.toString());
                                        PearsonsCorrelation pc = new PearsonsCorrelation();
                                        
                                        double[] maxDiamsAr = IO.integerArrayToDoubleArray(IO.mapStringIntegerToArraySortedByMapKeys(maxDiams));
                                        double[] numAtomsAr = IO.integerArrayToDoubleArray(numAtoms);
                                        double[] maxDegreesAr = IO.integerArrayToDoubleArray(IO.mapStringIntegerToArraySortedByMapKeys(maxDegrees));
                                        double[] avgDiamsAr = IO.mapStringDoubleToArraySortedByMapKeys(avgDiams);
                                        double[] avgDegreesAr = IO.mapStringDoubleToArraySortedByMapKeys(avgDegrees);
                                        
                                        Double corrMaxDiam2NumAtoms = pc.correlation(maxDiamsAr, numAtomsAr);
                                        Double corrMaxDiam2Degrees = pc.correlation(maxDiamsAr, maxDegreesAr);
                                        Double corrAvgDiam2NumAtoms = pc.correlation(avgDiamsAr, numAtomsAr);
                                        Double corrAvgDiam2MaxDegrees = pc.correlation(avgDiamsAr, maxDegreesAr);
                                        Double corrAvgDiam2AvgDegrees = pc.correlation(avgDiamsAr, avgDegreesAr);
                                        Double corrNumAtoms2MaxDegrees = pc.correlation(numAtomsAr, maxDegreesAr);
                                        Double corrNumAtoms2AvgDegrees = pc.correlation(numAtomsAr, avgDegreesAr);
                                        
                                        System.out.println("Pearson correlation, MaxDiameter to NumAtoms: " + corrMaxDiam2NumAtoms);
                                        System.out.println("Pearson correlation, MaxDiameter to MaxDegrees: " + corrMaxDiam2Degrees);
                                        System.out.println("Pearson correlation, AvgDiameter to NumAtoms: " + corrAvgDiam2NumAtoms);
                                        System.out.println("Pearson correlation, AvgDiameter to MaxDegrees: " + corrAvgDiam2MaxDegrees);
                                        System.out.println("Pearson correlation, AvgDiameter to AvgDegrees: " + corrAvgDiam2AvgDegrees);
                                        System.out.println("Pearson correlation, NumAtoms to MaxDegrees: " + corrNumAtoms2MaxDegrees);
                                        System.out.println("Pearson correlation, NumAtoms to AvgDegrees: " + corrNumAtoms2AvgDegrees);
                                    }
                                    
                                    // write properties of largest CC of graph, if the graph consists of more than one CC
                                    if(gpr.numConnectedComponents > 1) {
                                        propsComputationStartTime = new Date();                                    
                                        GraphPropResults gpr_lcc = gp_lcc.getGraphPropResults();                                    
                                        propsComputationEndTime = new Date();
                                        timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                        runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);                                    

                                        DBManager.writeAminoacidgraphStatsToDB(graph_db_id, Boolean.TRUE, gpr_lcc.numVertices, gpr_lcc.numEdges, gpr_lcc.minDegree, gpr_lcc.maxDegree, gpr_lcc.numConnectedComponents, gpr_lcc.graphDiameter, gpr_lcc.graphRadius, gpr_lcc.averageClusterCoefficient, gpr_lcc.averageShortestPathLength, gp_lcc.getDegreeDistributionUpTo(50), gpr_lcc.averageDegree, gpr_lcc.density, gp_lcc.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);

                                        if(rewireGraphsToCompareWithRandom) {
                                            // write properties of random equiv of largest CC of graph
                                            GraphRandomizer gr_lcc = new GraphRandomizer(lcc, edgeRewireProbability);
                                            GraphProperties gp_lcc_rand = new GraphProperties(lcc);
                                            propsComputationStartTime = new Date();
                                            GraphPropResults gpr_lcc_rand = gp_lcc_rand.getGraphPropResults();
                                            propsComputationEndTime = new Date();
                                            timeDiff = propsComputationEndTime.getTime() - propsComputationStartTime.getTime();//as given
                                            runtime_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiff);                                    
                                            if(writeRandomGraphPropsToDatabase) {
                                                String rand_graph_description = "largest CC random";
                                                //DBManager.deleteCustomGraphStatsFromDBByUniqueName(rand_graph_name);
                                                DBManager.writeCustomgraphStatsToDB(rand_graph_name, rand_graph_description, Boolean.TRUE, gpr_lcc_rand.numVertices, gpr_lcc_rand.numEdges, gpr_lcc_rand.minDegree, gpr_lcc_rand.maxDegree, gpr_lcc_rand.numConnectedComponents, gpr_lcc_rand.graphDiameter, gpr_lcc_rand.graphRadius, gpr_lcc_rand.averageClusterCoefficient, gpr_lcc_rand.averageShortestPathLength, gp_lcc_rand.getDegreeDistributionUpTo(50), gpr_lcc_rand.averageDegree, gpr_lcc_rand.density, gp_lcc_rand.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);                                            
                                            }
                                        }
                                    }
                                }
                                else {
                                    DP.getInstance().e("Main", "Could not write aa graph properties to DB, graph not found in database.");
                                }
                            } catch(SQLException e) {
                                DP.getInstance().e("SQL error while trying to store aa graph stats: '" + e.getMessage()+ "'.");
                            }
                        }

                    }

                    // write the AA graph to disc
                    String aagFile = outputDir + fs + subDirTree + pdbid + "_aagraph.gml";
                    if(writeStringToFile(aagFile, aag.toGraphModellingLanguageFormat())) {
                        if(! silent) {
                            System.out.println("  AAGraph for all chains written to file '" + aagFile + "'.");
                        }
                    }
                    else {
                        System.err.println("ERROR: Could not write AAGraph for all chains to file '" + aagFile + "'.");
                    }
                    
                    
                    
                    if(Settings.getBoolean("plcc_B_draw_aag")) {
                        Map<Integer, Color> cmap = new HashMap<>();
                        // fill color map
                        Color c;
                        Boolean use3Colors = false;
                        for(int i = 0; i < aag.getNumVertices(); i++) {
                            Residue r = aag.getVertex(i);
                            if(r != null) {
                                if(use3Colors) {
                                    c = ProteinGraphDrawer.getChemProp3Color(r.getChemicalProperty3OneLetterString());
                                } else {
                                    c = ProteinGraphDrawer.getChemProp5Color(r.getChemicalProperty5OneLetterString());
                                }
                                if(c != null) {
                                    cmap.put(i, c);
                                }
                            }
                        }
                        
                        // fill label map with 1 letter AA code
                        Map<Integer, String> lmap = new HashMap<>();
                        String l;
                        for(int i = 0; i < aag.getNumVertices(); i++) {
                            Residue r = aag.getVertex(i);
                            if(r != null) {
                                l = r.getAAName1();
                                if(l != null) {
                                    lmap.put(i, l);
                                }
                            }
                        }
                        
                                                
                        String aagDrawFileNoExt = outputDir + fs + subDirTree + pdbid + "_aagraph_vis";
                        //System.out.println("Drawing aag to base file '" + aagDrawFileNoExt + "'.");
                        SimpleGraphDrawer.drawSimpleGraphGrid(aagDrawFileNoExt, Settings.getAminoAcidGraphOutputImageFormats(), aag, cmap, lmap);                        
                    }
                    
                    // write the AA contact statistics matrix (by AA type, not single AA)
                    String aaMatrixFile = outputDir + fs + subDirTree + pdbid + "_aatypematrix.gml";
                    if(writeStringToFile(aaMatrixFile, aag.getAminoAcidTypeInteractionMatrixGML())) {
                        if(! silent) {
                            System.out.println("  AA type contact stats matrix for all chains written to file '" + aaMatrixFile + "'.");
                        }
                    }
                    else {
                        System.err.println("ERROR: Could not write AA type contact stats matrix for all chains to file '" + aaMatrixFile + "'.");
                    }
                    
                    //List<String> aaNames = new ArrayList<>();
                    //aaNames.addAll(Arrays.asList(AminoAcid.names3));
                    //AAInteractionNetwork aai = new AAInteractionNetwork(aaNames, aag.getAminoAcidTypeInteractionMatrix(false));
                    if(Settings.getBoolean("plcc_B_useDB")) {

                        if( ! DBManager.getAutoCommit()) {
                            DBManager.commit();
                        }
                        
                        if(Settings.getBoolean("plcc_B_compute_graph_metrics")) {
                            
                            int[] contactCountsByAATypeAbsolute = aag.getTotalContactCountByAAType();
                            int[] aacounts = aag.getAATypeCountsNoSum();
                            if(contactCountsByAATypeAbsolute.length != 20 || aacounts.length != 20) {
                                DP.getInstance().e("Main", "Wrong lenghts of amino acid stats arrays detected (" + contactCountsByAATypeAbsolute.length + "/" +  aacounts.length + "). This is a bug.");
                            }
                            
                            double[] contactCountsByAATypeNormalized = new double[20];
                            int sumContactCounts = 0;
                            for(int i = 0; i < contactCountsByAATypeAbsolute.length; i++) {
                                sumContactCounts += contactCountsByAATypeAbsolute[i];
                            }                                                        
                            
                            try {
                                DBManager.writeAbsoluteAATypeInteractionCountsToDB(pdbid, contactCountsByAATypeAbsolute);
                            } catch(SQLException e) {
                                DP.getInstance().e("Main", "Could not write absolute AA type interaction counts to database, skipping: '" + e.getMessage() + "'.");
                            }
                            
                            int sumAACounts = 0;    // the total number of AA residues
                            Boolean containsZeroCounts = false;
                            for(int i = 0; i < aacounts.length; i++) {
                                if(aacounts[i] == 0) {
                                    containsZeroCounts = true;
                                }
                                sumAACounts += aacounts[i];
                            }
                            
                            if( ! containsZeroCounts) { // avoid division by zero
                                for(int i = 0; i < contactCountsByAATypeNormalized.length; i++) {
                                    contactCountsByAATypeNormalized[i] = contactCountsByAATypeAbsolute[i] / aacounts[i];
                                }
                                
                                try {
                                    DBManager.writeNormalizedAATypeInteractionCountsToDB(pdbid, contactCountsByAATypeNormalized);
                                } catch(SQLException e) {
                                    DP.getInstance().e("Main", "Could not write normalized AA type interaction counts to database, skipping: '" + e.getMessage() + "'.");
                                }
                            }
                            
                            

                        }
                    }
                                       

                }
                }


            }
                        
            
            if(separateContactsByChain || Settings.getBoolean("plcc_B_AAgraph_perchain")) {
                String chainID;
                ArrayList<Chain> theChain;
                int numChainsHandled = 0;
                for(Chain c : handleChains) {
                    // add current chainName
                    theChain = new ArrayList<Chain>();   // it is a list, but only contains this single chainName
                    theChain.add(c);
                    
                    // compute chainName contacts
                    cInfoThisChain = calculateAllContactsLimitedByChain(residues, c.getPdbChainID());
                    
                    if(Settings.getBoolean("plcc_B_AAgraph_perchain")) {
                        AAGraph aag = new AAGraph(c.getResidues(), cInfoThisChain);
                        aag.setPdbid(pdbid);
                        aag.setChainid(c.getPdbChainID());
                        
                        if(Settings.getBoolean("plcc_B_useDB")) {
                            try {
                                if(Settings.getBoolean("plcc_B_write_graphstrings_to_database_aag")) {
                                    DBManager.writeAminoAcidGraphToDB(pdbid, c.getPdbChainID(), aag.toGraphModellingLanguageFormat(), aag.getNumVertices(), aag.getNumEdges());
                                }
                                else {
                                    DBManager.writeAminoAcidGraphToDB(pdbid, c.getPdbChainID(), null, aag.getNumVertices(), aag.getNumEdges());
                                }
                                
                                if(! silent) {
                                    System.out.println("Wrote chain '" + c.getPdbChainID() + "' AA graph of " + pdbid + " to DB.");
                                }
                            } catch(SQLException e) {
                                DP.getInstance().w("Main", "Could not write chain '" + c.getPdbChainID() + "' AA graph to DB: '" + e.getMessage() + "'.");
                            }
                        } 
                        
                        // write AA graph
                        String aagFile = outputDir + fs + pdbid + "_aagraph_chain_" + c.getPdbChainID() + ".gml";
                        if(writeStringToFile(aagFile, aag.toGraphModellingLanguageFormat())) {
                            if(! silent) {
                                System.out.println("  AAGraph for chain " + c.getPdbChainID() + " written to file '" + aagFile + "'.");
                            }
                        }
                        else {
                            System.err.println("ERROR: Could not write AAGraph for chain " + c.getPdbChainID() + " to file '" + aagFile + "'.");
                        }
                        
                        // write the AA contact statistics matrix (by AA type, not single AA)
                        String aaMatrixFile = outputDir + fs + pdbid + "_aatypematrix_chain_" + c.getPdbChainID() + ".gml";
                        if(writeStringToFile(aaMatrixFile, aag.getAminoAcidTypeInteractionMatrixGML())) {
                            if(! silent) {
                                System.out.println("  AA type contact stats matrix for chain " + c.getPdbChainID() + " written to file '" + aaMatrixFile + "'.");
                            }
                        }
                        else {
                            System.err.println("ERROR: Could not write AA type contact stats matrix for chain " + c.getPdbChainID() + " to file '" + aaMatrixFile + "'.");
                        }
                    }
                    
                    if(Settings.getBoolean("plcc_B_quit_after_aag")) {
                        System.out.println("Quitting after AAG computation as requested by settings.");
                        Main.doExit(0);
                    }
                    
                    if(separateContactsByChain) {
                        calculateSSEGraphsForChains(theChain, residues, cInfoThisChain, pdbid, outputDir);
                    }
                    
                    if(Settings.getBoolean("plcc_B_useDB")) {
                        if( ! DBManager.getAutoCommit()) {
                            DBManager.commit();
                        }
                    }
                    
                    numChainsHandled++;
                }
            }
            
            if(Settings.getBoolean("plcc_B_quit_after_aag")) {
                System.out.println("Quitting after AAG computation as requested by settings.");
                Main.doExit(0);
            }
            
            if( ! separateContactsByChain){  // no chainName separation active                
                calculateSSEGraphsForChains(handleChains, residues, cInfo, pdbid, outputDir);
                //calculateComplexGraph(handleChains, residues, cInfo, pdbid, outputDir);
                if(Settings.getBoolean("plcc_B_useDB")) {
                    if( ! DBManager.getAutoCommit()) {
                        DBManager.commit();
                    }
                }
            }
            
            if(! silent) {
                System.out.println("All " + handleChains.size() + " chains done.");
            }
            
            // discard contact info to free memory
            cInfo = null;
            cInfoThisChain = null;
        }
        else {
            if(! silent) {
                System.out.println("  Not calculating SSEs and not drawing graphs as requested.");
            }
        }        
        
        // ******************************************** HTML output ******************************************************* //
        
        if(Settings.getBoolean("plcc_B_output_textfiles_dir_tree_html")) {
            if(! silent) {
                System.out.println(" Producing web pages for PDB " + pdbid + ".");
            }
            File outputBaseDir = new File(outputDir);   // the path where the subdir tree will be created (this contains 'ic/8icd/A'-like directories in the end)
            File outputDirProtein = IO.generatePDBstyleSubdirTreeName(new File(outputDir), pdbid);  // looks like '$OUTPUTBASEDIR/ic/8icd/'
            if( ! IO.dirExistsIsDirectoryAndCanWrite(outputDirProtein)) {
                IO.createDirIfItDoesntExist(outputDirProtein);
            }                        
            
            HtmlGenerator htmlGen = new HtmlGenerator(outputDirProtein);
            CssGenerator cssGen = new CssGenerator();
            String cssFilePath = outputBaseDir.getAbsolutePath() + fs + "vplgweb.css";            
            String cssFilePathRed = outputBaseDir.getAbsolutePath() + fs + "vplgweb_red.css";
            String cssFilePathBlue = outputBaseDir.getAbsolutePath() + fs + "vplgweb_blue.css";
            String cssFilePathGreen = outputBaseDir.getAbsolutePath() + fs + "vplgweb_green.css";
            String fsWeb = "/"; // the internet is UNIX                        
            htmlGen.setRelativeCssFilePathsFromBasedir(new String[] { fsWeb + "vplgweb.css", fsWeb + "vplgweb_red.css", fsWeb + "vplgweb_blue.css", fsWeb + "vplgweb_green.css" });        // no, this ain't beautiful            
            htmlGen.setCssTitles(new String[] { "default", "red", "blue", "green" });
            
            if(Settings.getBoolean("plcc_B_output_textfiles_dir_tree_core_html")) {  
                if(! silent) {
                    System.out.println("  Writing core webpages. The base output directory is '" + outputBaseDir.getAbsolutePath() + "'.");
                }
                if(cssGen.writeDefaultCssFileTo(new File(cssFilePath), CssGenerator.COLORSCHEME_RED)) {
                    if(! silent) {
                    System.out.println("   Wrote default CSS file for web pages to '" + cssFilePath + "'.");
                    }
                }
                if(cssGen.writeDefaultCssFileTo(new File(cssFilePathRed), CssGenerator.COLORSCHEME_RED)) {
                    if(! silent) {                    
                        System.out.println("   Wrote red CSS file for web pages to '" + cssFilePathRed + "'.");
                    }
                }
                if(cssGen.writeDefaultCssFileTo(new File(cssFilePathBlue), CssGenerator.COLORSCHEME_BLUE)) {
                    if(! silent) {
                        System.out.println("   Wrote blue CSS file for web pages to '" + cssFilePathBlue + "'.");
                    }
                }
                if(cssGen.writeDefaultCssFileTo(new File(cssFilePathGreen), CssGenerator.COLORSCHEME_GREEN)) {
                    if(! silent) {
                        System.out.println("   Wrote green CSS file for web pages to '" + cssFilePathGreen + "'.");
                    }
                }
                htmlGen.generateCoreWebpages(outputBaseDir);
            }
                        
            if(! silent) {
                System.out.println("  Writing protein-specific webpages. The protein-specific output directory is '" + outputDirProtein.getAbsolutePath() + "'.");
            }
            htmlGen.generateAllWebpagesForResult(ProteinResults.getInstance());
            if(! silent) {
                System.out.println("  Web pages done for PDB " + pdbid + ".");
                System.out.println(" Produced all web pages.");
            }
        }

        //drawTGFGraph("graph.tgf", "graph.tgf.png");       //DEBUG
        //();        //DEBUG
        
        // writing contacts to the database
        /*
        if(Settings.getBoolean("plcc_B_useDB")) {
            String somePDBID = "7tim";
            String chain1 = "A";      
            String chain2 = "B";            
            System.out.println("[DEBUG] =====Writing complex contacts.=====");
            try {
                if(DBManager.writeInterchainContactsToDB(somePDBID, chain1, chain2, 5, 10, 15, 20, 17, 13, 0)) {
                    System.out.println("DEBUG: Wrote interchain contacts between chains " + chain1 + " and " + chain2 + " of PDB " + somePDBID + " to DB.");
                } else {
                    System.err.println("ERROR: Could not write interchain contacts between chains " + chain1 + " and " + chain2 + " of PDB " + somePDBID + " to DB. Skipping.");
                }
            } catch(Exception e) {
                System.err.println("ERROR: Could not write interchain contacts between chains " + chain1 + " and " + chain2 + " of PDB " + somePDBID + " to DB. Skipping. Error was: '" + e.getMessage() + "'.");
                //System.out.println("ERROR: Could not write interchain contacts between chains " + chain1 + " and " + chain2 + " of PDB " + somePDBID + " to DB. Skipping. Error was: '" + e.getMessage() + "'.");
            }
            System.out.println("[DEBUG] =====Done writing complex contacts.=====");
        }
        */
        

        // ****************************************************    all done    ********************************************************** //
        
        
        
        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            print_debug_malfunction_warning();            
        }
        
        if(deleteFilesOnExit.size() > 0) {
            if(! silent) {
                System.out.print("Deleting " + deleteFilesOnExit.size() + " temporary files... ");
            }
            Integer numDel = IO.deleteFiles(deleteFilesOnExit);
            if(! silent) {
                System.out.print(numDel + " ok.\n");
            }
        }
        
        Date totalComputationEndTime = new Date();
        long timeDiffTotal = totalComputationEndTime.getTime() - computationStartTime.getTime();//as given
        long runtimeTotal_secs = TimeUnit.MILLISECONDS.toSeconds(timeDiffTotal);
        int[] compTimes = splitDurationToComponentTimes(runtimeTotal_secs);
        
        if(Settings.getBoolean("plcc_B_useDB")) {
            try {
                DBManager.updateProteinTotalRuntimeInDB(pdbid, runtimeTotal_secs);
            }catch(SQLException e) {
                DP.getInstance().w("Main", "Could not set total runtime for PDB file '" + pdbid + "' in database.");
            }
            
            try {
                DBManager.updateProteinRunCompletedInDB(pdbid, Boolean.TRUE);
            }catch(SQLException e) {
                DP.getInstance().w("Main", "Could not set total protein inserted state for PDB file '" + pdbid + "' in database.");
            }
        }
        
        if(Settings.getBoolean("plcc_B_useDB")) {
            if( ! DBManager.getAutoCommit()) {
                DBManager.commit();            
            }
            DBManager.closeConnection();
        }
        
        if(! silent) {
            System.out.println("(Too much clutter? Try the '--silent' command line option.)");
            System.out.println("All done, exiting. Total runtime was " + runtimeTotal_secs + " seconds ("+compTimes[0]+":" + String.format("%02d:%02d", compTimes[1], compTimes[2])+" hms) for " + residues.size() + " residues.");
        } else {
            if(Settings.getBoolean("plcc_B_print_silent_notice")) {
                String endTime = new SimpleDateFormat("HH:mm:ss").format(new Date());
                System.out.println("[PLCC] [" + pdbid + "] [" + endTime + "] " + "All done silently, exiting. Total runtime was " + runtimeTotal_secs + " seconds ("+compTimes[0]+":" + String.format("%02d:%02d", compTimes[1], compTimes[2])+" hms). "+residues.size()+ "residues.");
            }
        }
        System.exit(0);

    }

    
    /**
     * Takes a duration given in seconds, and gives it in the hours + minutes + seconds format.
     * @param secsIn the duration in number of seconds
     * @return the duration as hours+mins+secs, in an int array (in the order just given)
     */
    public static int[] splitDurationToComponentTimes(long secsIn) {
        long longVal = secsIn;
        int hours = (int) longVal / 3600;
        int remainder = (int) longVal - hours * 3600;
        int mins = remainder / 60;
        remainder = remainder - mins * 60;
        int secs = remainder;

        int[] ints = {hours , mins , secs};
        return ints;
    }

    /**
     * Draws the image of a graph from the file 'tpgFile' (which is expected to contain a graph in the Trivial Graph Format) and writes it to the PNG file img.
     * @param tgfFile the path to the TGF file which contains the graph to draw
     * @param imgFilePathNoExt the output path where the resulting image should be written
     * @param formats a list of image formats
     */
    public static void drawTGFGraph(String tgfFile, String imgFilePathNoExt, IMAGEFORMAT[] formats, HashMap<Integer, String> vertexMarkings) {

        //System.out.println("Testing tgf implementation using file '" + tgfFile + "'.");

        ProtGraph pg = ProtGraphs.fromTrivialGraphFormatFile(tgfFile);
        System.out.println("  Loaded graph with " + pg.numVertices() + " vertices and " + pg.numEdges() + " edges, drawing to base file '" + imgFilePathNoExt + "'.");
        ProteinGraphDrawer.drawProteinGraph(imgFilePathNoExt, true, formats, pg, vertexMarkings, new ArrayList<String>());
        //pg.print();
        System.out.println("  Graph image written to base file '" + imgFilePathNoExt + "'.");
    }
    
    
    
    /**
     * Draws the image of a graph from the file 'plccGraphFile' (which is expected to contain a graph in PLCC format) and writes it to the PNG file 'img'.
     * @param plccGraphFile the input file in plcc format
     * @param imgNoExt the path where to write the output image (without the file extension)
     */
    public static void drawPlccGraphFromFile(String plccGraphFile, String imgNoExt, Boolean drawFoldingGraphsAsWell) {

        //System.out.println("Testing plcc graph reading implementation using file '" + plccGraphFile + "'.");

        String graphString = FileParser.slurpFileSingString(plccGraphFile);
        
        if(Settings.getInteger("plcc_I_debug_level") > 0) {
            ProtGraphs.printPlccMetaData(graphString);
        }
                
        IMAGEFORMAT[] formats = new IMAGEFORMAT[] { DrawTools.IMAGEFORMAT.PNG };
        ProtGraph pg = ProtGraphs.fromPlccGraphFormatString(graphString);
        System.out.println("  Loaded graph with " + pg.numVertices() + " vertices and " + pg.numEdges() + " edges, drawing to base file '" + imgNoExt + "'.");
        ProteinGraphDrawer.drawProteinGraph(imgNoExt, false, formats, pg, new HashMap<Integer, String>(), new ArrayList<String>()); 
        System.out.println("  Protein graph image written to base file '" + imgNoExt + "'.");

        //if(drawFoldingGraphsAsWell) {
        //    calculateFoldingGraphsForSSEGraph(pg, Settings.get("plcc_S_output_dir"));
        //}        
        //else {
        //    System.err.println("ERROR: Drawing of graph failed.");
        //}       
    }
    
    
    /**
     * Draws the image of a graph from the file 'plccGraphFile' (which is expected to contain a graph in PLCC format) and writes it to the PNG file 'img'.
     * @param g_pdbid the PDB id
     * @param g_chainid the chainName ID
     * @param drawFoldingGraphsAsWell whether to draw FGs
     * @param g_graphtype the graph type
     * @param outputImgNoExt the path where to write the output image (without the file extension)
     */
    public static void drawPlccGraphFromDB(String g_pdbid, String g_chainid, String g_graphtype, String outputImgNoExt, Boolean drawFoldingGraphsAsWell) {

        System.out.println("Retrieving " + g_graphtype + " graph for PDB entry " + g_pdbid + " chain " + g_chainid + " from DB.");
        
        String graphString = null;
        try { graphString = DBManager.getGraphStringGML(g_pdbid, g_chainid, g_graphtype); }
        catch (SQLException e) { System.err.println("ERROR: SQL: Drawing of graph from DB failed: '" + e.getMessage() + "'."); return; }
        
        if(graphString == null) {
            System.err.println("ERROR: No such graph in the database to draw.");
            return;
        }
        
        ProtGraph pg = ProtGraphs.fromPlccGraphFormatString(graphString);
        IMAGEFORMAT[] formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG };
        System.out.println("  Loaded graph with " + pg.numVertices() + " vertices and " + pg.numEdges() + " edges, drawing to base file '" + outputImgNoExt + "'.");
        ProteinGraphDrawer.drawProteinGraph(outputImgNoExt, false, formats, pg, new HashMap<Integer, String>(), new ArrayList<String>());
        System.out.println("  Protein graph image written to base file '" + outputImgNoExt + "'.");

        //if(drawFoldingGraphsAsWell) {
        //    calculateFoldingGraphsForSSEGraph(pg, Settings.get("plcc_S_output_dir"));
        //}
    }


    /**
     * Just a test function for the DBManager class, does not do anything useful. Creates some tables, inserts stuff and deletes it when it is done.
     */
    public static void dbTesting() {


        System.out.println("Starting DB tests.");
        ArrayList<ArrayList<String>> tableData;
        Integer numRows = null;

        if(Settings.getBoolean("plcc_B_useDB")) {
            System.out.println("DB OK");

            System.out.println("Starting test: INSERT");
            numRows = DBManager.doInsertQuery("INSERT INTO test (name, age) VALUES ('john', 29);");
            System.out.println("Test INSERT done, " + numRows + " affected rows.");
            numRows = DBManager.doInsertQuery("INSERT INTO test (name, age) VALUES ('john', 28);");
            System.out.println("Test INSERT done, " + numRows + " affected rows.");
            numRows = DBManager.doInsertQuery("INSERT INTO test (name, age) VALUES ('john', 27);");
            System.out.println("Test INSERT done, " + numRows + " affected rows.");


            System.out.println("Starting test: UPDATE");
            numRows = DBManager.doUpdateQuery("UPDATE test SET age = 30 WHERE name = 'john' AND age = 29;");
            System.out.println("Test UPDATE done, " + numRows + " affected rows.");


            System.out.println("Starting test: SELECT");            
            try {
                tableData = DBManager.doSelectQuery("SELECT name, age FROM test WHERE name = 'john';");
                
                System.out.print("----------\n");
                
                for(ArrayList<String> row : tableData) {
                    for(String value : row) {
                        System.out.print(value + " ");
                    }
                    System.out.print("\n");
                }
                
                System.out.print("----------\n");
            } catch(Exception e) {
                System.err.println("ERROR: SELECT query failed.");
                e.printStackTrace();
            }


            System.out.println("Starting test: DELETE");
            numRows = DBManager.doDeleteQuery("DELETE FROM test WHERE name = 'john' AND age = 30;");
            System.out.println("Test DELETE done, " + numRows + " affected rows.");


            try {
                tableData = DBManager.doSelectQuery("SELECT name, age FROM test WHERE name = 'john';");

                System.out.print("----------\n");

                for(ArrayList<String> row : tableData) {
                    for(String value : row) {
                        System.out.print(value + " ");
                    }
                    System.out.print("\n");
                }

                System.out.print("----------\n");
            } catch(Exception e) {
                System.err.println("ERROR: SELECT query failed.");
                e.printStackTrace();
            }

            
            numRows = DBManager.doDeleteQuery("DELETE FROM test WHERE name = 'john';");
            System.out.println("Test DELETE done, " + numRows + " affected rows. All tests done.");
        }
        else {
            DP.getInstance().w("dbTesting(): Not using DB or DB connection failed. Skipping DB tests.");
            System.exit(1);
        }
    }


    /**
     * Writes data on the current protein to the database and deletes any old data in there which has the same PDB identifier.
     * @param pdbid the PDB ID to use
     */
    public static void writeProteinDataToDatabase(String pdbid, int numResidues) {
        
        Boolean silent = Settings.getBoolean("plcc_B_silent");
        
        // meta Data in CIF files is created while parsing the file once
        //     -> no need to call a function to create it, just get it!
        HashMap<String, String> md = FileParser.getMetaData();

        Double resolution = -1.0;
        try {
            resolution = Double.valueOf(md.get("resolution"));
        } catch (Exception e) {
            resolution = -1.0;
            DP.getInstance().w("Could not determine resolution of PDB file for protein '" + pdbid + "', assuming NMR with resolution '" + resolution + "'.");
            
        }
                
        //pdb_id, title, header, keywords, experiment, resolution
        if(Settings.getBoolean("plcc_B_useDB")) {
            // Try to delete the protein from the DB in case it is already in there. This won't hurt if it is not.
            
            int numDel = 0;            
            numDel = DBManager.deletePdbidFromDB(pdbid);
            if(! silent) {
                System.out.println("  Deleted " + numDel + " protein(s) with PDB ID '" + pdbid + "' from DB.");
            }
                        
            try {
                DBManager.writeProteinToDB(pdbid, md.get("title"), md.get("header"), md.get("keywords"), md.get("experiment"), resolution, numResidues, Boolean.parseBoolean(md.get("isLarge")));
                if(! silent) {
                    System.out.println("  Info on protein '" + pdbid + "' written to DB.");
                }
            }catch (Exception e) {
                DP.getInstance().w("Could not write info on protein '" + pdbid + "' to DB: '" + e.getMessage() + "'.");
            }
        }        
    }
    
    /**
     * Calculates all SSE graph types which are configured in the config file for all given chains.
     * @param allChains a list of chains, each chainName will be handled separately
     * @param resList a list of residues
     * @param resContacts a list of residue contacts
     * @param pdbid the PDBID of the protein, required to name files properly etc.
     * @param outputDir where to write the output files. the filenames are deduced from graph type and pdbid.
     */
    public static void calculateSSEGraphsForChains(List<Chain> allChains, List<Residue> resList, ArrayList<MolContactInfo> resContacts, String pdbid, String outputDir) {
        Boolean silent = Settings.getBoolean("plcc_B_silent");
               
        //System.out.println("calculateSSEGraphsForChains: outputDir='" + outputDir + "'.");
        Chain c;
        List<SSE> chainDsspSSEs = new ArrayList<>();
        List<SSE> chainLigSSEs = new ArrayList<>();
        List<SSE> chainPtglSSEs = new ArrayList<>();
        List<SSE> allChainSSEs = new ArrayList<>();
        String fs = System.getProperty("file.separator");
               
        // meta Data in CIF files is created while parsing the file once
        //     -> no need to call a function to create it, just get it!
        HashMap<String, String> md = FileParser.getMetaData();        
        
        // check which chains belong to the same macro molecule
        Map<String, List<String>> macroMoleculesOfPDBfileToChains = new HashMap<>();    // key of the map is the MOL_ID
        Map<String, Map<String, String>> macroMolecules = new HashMap<>();   // each inner hashmap contains the properties of a macromolecule, "name" => the_name, "id" => MOL_ID, .... The outer string is the mol_ID
       
        
        // handle all chains
        ProteinChainResults pcr;
        for(Integer i = 0; i < allChains.size(); i++) {
            c = allChains.get(i);
            String chain = c.getPdbChainID();
            if(! silent) {
                System.out.println("  +++++ Handling chain '" + chain + "'. +++++");
            }
                        
            // CIF parser does not parse all protein chain meta information, return in these cases default values
            ProtMetaInfo pmi = FileParser.getMetaInfo(pdbid, chain);
                       
            //pmi.print();
            
            
            c.setMacromolID(pmi.getMacromolID());
            c.setMacromolName(pmi.getMolName());
            
            // collect macromol data, will be used to write MM to database after this chainName loop
            Map<String, String> tmpMacroMol = new HashMap<>();
            tmpMacroMol.put("pdb_mol_id", pmi.getMacromolID()); // not strictly needed, it is also put as the key for this MM later
            tmpMacroMol.put("pdb_mol_name", pmi.getMolName());
            tmpMacroMol.put("pdb_org_sci", pmi.getOrgScientific());
            tmpMacroMol.put("pdb_org_common", pmi.getOrgCommon());
            tmpMacroMol.put("pdb_all_chains", pmi.getAllMolChains());
            tmpMacroMol.put("pdb_ec_number", pmi.getECNumber());
            
            macroMolecules.put(pmi.getMacromolID(), tmpMacroMol);
            
            md.put("pdb_mol_id", pmi.getMacromolID());
            md.put("pdb_mol_name", pmi.getMolName());
            md.put("pdb_org_sci", pmi.getOrgScientific());
            md.put("pdb_org_common", pmi.getOrgCommon());
            md.put("pdb_all_chains", pmi.getAllMolChains());
            md.put("pdb_ec_number", pmi.getECNumber());
            
            if(! silent) {
                System.out.println("    Chain '" + chain + "' MOL_ID is '" + pmi.getMacromolID() + "', MOL_NAME is '" + pmi.getMolName() + "'.");
            }
            
            // keep track of mm
            if( ! macroMoleculesOfPDBfileToChains.containsKey(pmi.getMacromolID())) {
                macroMoleculesOfPDBfileToChains.put(pmi.getMacromolID(), new ArrayList<String>());
            }
            macroMoleculesOfPDBfileToChains.get(pmi.getMacromolID()).add(allChains.get(i).getPdbChainID());
            
            pcr = new ProteinChainResults(c.getPdbChainID());
            // register results for chainName
            ProteinResults.getInstance().addProteinChainResults(pcr, chain);
            pcr.setChainMetaData(pmi);

            
            
            
            if(Settings.getBoolean("plcc_B_useDB")) {
                String ligName3Trimmed;
                try {
                    if(DBManager.writeChainToDB(chain, pdbid, pmi.getMacromolID(), pmi.getMolName(), pmi.getOrgScientific(), pmi.getOrgCommon())) {
                        if(! silent) {
                            System.out.println("    Info on chain '" + chain + "' of protein '" + pdbid + "' written to DB.");
                        }
                        
                        Long chainDbId = DBManager.getDBChainID(pdbid, chain);
                        
                        if(chainDbId >= 1) {
                            for(Residue ligand : c.getAllLigandResidues()) {
                                ligName3Trimmed = ligand.getTrimmedName3();
                                
                                // TODO: these 3 lines are not thread-safe, running several plcc instances in parallel may lead to race conditions
                                DBManager.writeLigandToDBUnlessAlreadyThere(ligName3Trimmed, ligand.getLigName(), ligand.getLigFormula(), ligand.getLigSynonyms());
                                if( ! DBManager.assignmentLigandToProteinChainExistsInDB(chainDbId, ligName3Trimmed)) {
                                    DBManager.assignLigandToProteinChain(chainDbId, ligName3Trimmed);
                                }
                            }
                        }
                        else {
                            DP.getInstance().w("Could not retrieve chain info from DB, writing it failed it seems.");
                        }
                    }
                    else {
                        DP.getInstance().w("Could not write info on chain '" + chain + "' of protein '" + pdbid + "' to DB.");
                    }
                }
                catch(Exception e) {
                    DP.getInstance().w("DB: Could not reset DB connection: '" + e.getMessage() + "'.");
                }
            }

            // determine SSEs for this chainName
            if(! silent) {
                System.out.println("    Creating all SSEs for chain '" + chain + "' consisting of " + c.getResidues().size() + " residues.");
            }
            chainDsspSSEs = createAllDsspSSEsFromResidueList(c.getResidues());
            
            if(chainDsspSSEs.isEmpty()) {
                if(Settings.getBoolean("plcc_B_skip_empty_chains")) {
                    if(! silent) {
                        System.out.println("  +++++ Skipping chain " + chain + " due to empty residue list. +++++");
                        
                    }
                    continue;
                }
            }
            
            if(Settings.getInteger("plcc_I_debug_level") > 0) {
                printSSEList(chainDsspSSEs, "DSSP");
            }
            
            if(Settings.getBoolean("plcc_B_ptgl_text_output")) {
                String sseMappingsFile = Settings.get("plcc_S_output_dir") + fs + pdbid.toLowerCase()  + "_" + chain + ".ssemap";
                writeSSEMappings(sseMappingsFile, c, pdbid);
            }
            
            chainPtglSSEs = createAllPtglSSEsFromDsspSSEList(chainDsspSSEs);
            
            if(Settings.getInteger("plcc_I_debug_level") > 0) {
                printSSEList(chainPtglSSEs, "PTGL");
            }
            
            chainLigSSEs = createAllLigandSSEsFromResidueList(c.getResidues(), chainDsspSSEs);
            allChainSSEs = mergeSSEs(chainPtglSSEs, chainLigSSEs);
            if(! silent) {
                System.out.println("    Added " + chainLigSSEs.size() + " ligand SSEs to the SSE list, now at " + allChainSSEs.size() + " SSEs.");
                System.out.print("    SSEs: ");

                for(Integer j = 0; j < allChainSSEs.size(); j++) {
                    System.out.print(allChainSSEs.get(j).getSseType());
                }
                System.out.print("\n");
            }
                        
            // SSEs have been calculated, now assign the PTGL labels and sequential numbers on the chainName
            for(Integer j = 0; j < allChainSSEs.size(); j++) {
                allChainSSEs.get(j).setSeqSseChainNum(j + 1);   // This is the correct value, determined from the list of all valid SSEs of this chainName
                allChainSSEs.get(j).setSseIDPtgl(getPtglSseIDForNum(j));

                if(Settings.getBoolean("plcc_B_useDB")) {
                    
                    //if( ! Settings.getBoolean("plcc_B_db_use_batch_inserts")) {
                                                               
                        try {
                           SSE ssej = allChainSSEs.get(j);
                           Integer ssePositionInChain = j + 1;
                           Long insertID = DBManager.writeSSEToDB(pdbid, chain, ssej.getStartDsspNum(), ssej.getEndDsspNum(), ssej.getStartPdbResID(), ssej.getEndPdbResID(), ssej.getAASequence(), ssej.getSSETypeInt(), ssej.getTrimmedLigandName3(), ssePositionInChain); 
                           //System.out.println("  Info on SSE #" + (j + 1) + " of chainName '" + c.getPdbChainID() + "' of protein '" + pdbid + "' written to DB.");
                           if(insertID > 0) {
                               DBManager.writeEmptySecondatEntryForSSE(insertID);
                           }
                           else {
                               DP.getInstance().w("Main", "Insert ID of SSE is < 0, insert failed. Cannot write the secondat entry for the SSE to the DB.");
                           }
                        }
                        catch(Exception e) {
                            DP.getInstance().w("Could not write info on SSE # " + j + " of chain '" + chain + "' of protein '" + pdbid + "' to DB.");
                        }
                    //}
                }
            }
            
            // batch insert all SSEs at once of appropriate
            /*
            if(Settings.getBoolean("plcc_B_useDB") &&  Settings.getBoolean("plcc_B_db_use_batch_inserts")) {
                try {
                    int insertCount = DBManager.writeAllSSEsOfChainToDB(pdbid, chainName, allChainSSEs);
                    if(insertCount != allChainSSEs.size()) {
                        DP.getInstance().e("Main", "Only " + insertCount + " of the " + allChainSSEs.size() + " SSEs were written to the DB. Exiting.");
                        Main.doExit(1);
                    }
                } catch(SQLException e) {
                    DP.getInstance().e("Main", "Writing all chainName SSE list to DB failed: '" + e.getMessage() + "'.");
                }
            }
            */


            //printSSEList(chainDsspSSEs, "DSSP SSEs of chainName '" + c.getPdbChainID() + "'");
            //printSSEList(chainPtglSSEs, "PTGL SSEs of chainName '" + c.getPdbChainID() + "'");
            //printSSEList(chainLigSSEs, "Ligand SSEs of chainName '" + c.getPdbChainID() + "'");
            //printSSEList(allChainSSEs, "All SSEs of chainName '" + c.getPdbChainID() + "'");


            // ************* Calculate the different graph types *************** //
            //List<String> graphTypes = Arrays.asList("albe", "albelig", "beta", "betalig", "alpha", "alphalig");       // old hardcoded stuff
            //List<String> graphTypes = Arrays.asList("albelig");                                                       // old hardcoded stuff
            
            // read the list of requested graph types from the settings
            List<String> graphTypes = new ArrayList<String>();

            if(Settings.getBoolean("plcc_B_graphtype_albe")) { graphTypes.add(SSEGraph.GRAPHTYPE_ALBE); }
            if(Settings.getBoolean("plcc_B_graphtype_albelig")) { graphTypes.add(SSEGraph.GRAPHTYPE_ALBELIG); }
            if(Settings.getBoolean("plcc_B_graphtype_alpha")) { graphTypes.add(SSEGraph.GRAPHTYPE_ALPHA); }
            if(Settings.getBoolean("plcc_B_graphtype_alphalig")) { graphTypes.add(SSEGraph.GRAPHTYPE_ALPHALIG); }
            if(Settings.getBoolean("plcc_B_graphtype_beta")) { graphTypes.add(SSEGraph.GRAPHTYPE_BETA); }
            if(Settings.getBoolean("plcc_B_graphtype_betalig")) { graphTypes.add(SSEGraph.GRAPHTYPE_BETALIG); }
            
            
            String fileNameWithExtension = null;
            String fileNameWithoutExtension = null;
            String filePathImg = null;
            String filePathGraphs = null;
            String filePathHTML = null;
            String imgFile = null;


            for(String gt : graphTypes) {
                // create the protein graph for this graph type
                //System.out.println("SSEs: " + allChainSSEs);                

                ProtGraph pg = calcGraphType(gt, allChainSSEs, c, resContacts, pdbid);
                pg.setInfo(pdbid, chain, c.getMacromolID(), gt);
                pg.addMetadata(md);
                
                pcr.addProteinGraph(pg, gt);
                
                
                if(Settings.getBoolean("plcc_B_debug_compareSSEContacts")) {
                    if(gt.equals(SSEGraph.GRAPHTYPE_ALBE)) {
                        if(! silent) {
                            System.out.println("Comparing calculated SSE contacts with those in the file '" + Settings.get("plcc_S_debug_compareSSEContactsFile") + "'...");
                        }
                        FileParser.compareSSEContactsWithGeoDatFile(Settings.get("plcc_S_debug_compareSSEContactsFile"), pg);
                    }        
                    else {
                        if(! silent) {
                            System.out.println("INFO: SSE contact comparison request ignored since this is not an albe graph.");
                        }
                    }
                }
                
                Integer isoLig = pg.numIsolatedLigands();
                String coilsUsed = "";
                if(Settings.getBoolean("plcc_B_include_coils")) {
                    coilsUsed = " including coils";
                }
                if(isoLig > 0) {
                    if(! silent) {
                        System.out.println("      The " + gt + " graph of " + pdbid + " chain " + chain + coilsUsed + " contains " + isoLig + " isolated ligands.");
                    }
                }
                
                

                // draw the protein graph image

                filePathImg = outputDir;
                filePathGraphs = outputDir;
                filePathHTML = outputDir;
                String coils = "";
                if(Settings.getBoolean("plcc_B_include_coils")) {
                    //System.out.println("  Considering coils, this may fragment SSEs.");
                    coils = "_coils";
                }
                fileNameWithoutExtension = pdbid + "_" + chain + "_" + gt + coils + "_PG";
                fileNameWithExtension = fileNameWithoutExtension + Settings.get("plcc_S_img_output_fileext");
                
                //pg.toFile(file + ".ptg");
                //pg.print();                
                // Create the file in a subdir tree based on the protein meta data if requested
                if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                   
                    File targetDir = IO.generatePDBstyleSubdirTreeNameWithChain(new File(outputDir), pdbid, chain);
                    if(targetDir != null) {
                        ArrayList<String> errors = IO.createDirIfItDoesntExist(targetDir);
                        if( ! errors.isEmpty()) {
                            for(String err : errors) {
                                System.err.println("ERROR: " + err);
                            }
                        } else {
                            filePathImg = targetDir.getAbsolutePath();
                            filePathGraphs = targetDir.getAbsolutePath();
                            filePathHTML = targetDir.getAbsolutePath();
                        }                    
                    } else {
                        System.err.println("ERROR: Could not determine PDB-style subdir path name.");
                    }
                }
                
                String gmlFileNoPath = "";    // for DB path reconstruction later
                String jsonFileNoPath = "";
                String xmlFileNoPath = "";
                String dotlanguageFileNoPath = "";
                String kavoshFileNoPath = "";
                String plccFileNoPath = "";
                String gexfFileNoPath = "";
                String cytoscapejsFileNoPath = "";
                String msvgFileNoPath = "";
                
                HashMap<String, String> writtenFormatsDBFilesNoPath = new HashMap<>();
                
                String graphFormatsWritten = "";
                Integer numFormatsWritten = 0;
                if(Settings.getBoolean("plcc_B_output_GML")) {
                    String gmlFile = filePathGraphs + fs + fileNameWithoutExtension + ".gml";
                    gmlFileNoPath = fileNameWithoutExtension + ".gml";
                    if(IO.stringToTextFile(gmlFile, pg.toGraphModellingLanguageFormat())) {
                        graphFormatsWritten += "gml "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_GML, gmlFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_GML, new File(gmlFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_TGF")) {
                    String tgfFile = filePathGraphs + fs + fileNameWithoutExtension + ".tgf";
                    if(IO.stringToTextFile(tgfFile, pg.toTrivialGraphFormat())) {
                        graphFormatsWritten += "tgf "; numFormatsWritten++;
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_TGF, new File(tgfFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_DOT")) {
                    String dotLangFile = filePathGraphs + fs + fileNameWithoutExtension + ".gv";
                    dotlanguageFileNoPath = fileNameWithoutExtension + ".gv";
                    if(IO.stringToTextFile(dotLangFile, pg.toDOTLanguageFormat())) {
                        graphFormatsWritten += "gv "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_DOTLANGUAGE, dotlanguageFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_DOTLANGUAGE, new File(dotLangFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_kavosh")) {
                    String kavoshFile = filePathGraphs + fs + fileNameWithoutExtension + ".kavosh";
                    kavoshFileNoPath = fileNameWithoutExtension + ".kavosh";
                    if(IO.stringToTextFile(kavoshFile, pg.toKavoshFormat())) {
                        graphFormatsWritten += "kavosh "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_KAVOSH, kavoshFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_KAVOSH, new File(kavoshFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_eld")) {
                    String elFile = filePathGraphs + fs + fileNameWithoutExtension + ".el_edges";
                    String nodeTypeListFile = filePathGraphs + fs + fileNameWithoutExtension + ".el_ntl";
                    if(IO.stringToTextFile(elFile, pg.toEdgeList()) && IO.stringToTextFile(nodeTypeListFile, pg.getNodeTypeList())) {
                        graphFormatsWritten += "el "; numFormatsWritten++;
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_EDGELIST, new File(elFile));
                    }
                }
                // write the SSE info text file for the image (plcc graph format file)
                if(Settings.getBoolean("plcc_B_output_plcc")) {
                    String plccGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".plg";
                    plccFileNoPath = fileNameWithoutExtension + ".plg";
                    if(IO.stringToTextFile(plccGraphFile, pg.toVPLGGraphFormat())) {
                        graphFormatsWritten += "plg "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_VPLG, plccFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_VPLG, new File(plccGraphFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_perlfg")) {
                    String perlGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".graph";
                    if(IO.stringToTextFile(perlGraphFile, pg.toPTGLGraphFormatPerl())) {
                        graphFormatsWritten += "perlfg "; numFormatsWritten++;
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_PERLFOLDINGGRAPHSCRIPT, new File(perlGraphFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_json")) {
                    String jsonGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".json";
                    jsonFileNoPath = fileNameWithoutExtension + ".json";
                    if(IO.stringToTextFile(jsonGraphFile, pg.toJSONFormat())) {
                        graphFormatsWritten += "json "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_JSON, jsonFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_JSON, new File(jsonGraphFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_msvg")) {
                    String msvgGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".man.svg";
                    msvgFileNoPath = fileNameWithoutExtension + ".man.svg";
                    if(IO.stringToTextFile(msvgGraphFile, pg.toManualSVGFormat())) {
                        graphFormatsWritten += "msvg "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_MANUALSVG, msvgFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_MANUALSVG, new File(msvgGraphFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_gexf")) {
                    String gexfGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".gexf";
                    gexfFileNoPath = fileNameWithoutExtension + ".gexf";
                    if(IO.stringToTextFile(gexfGraphFile, pg.toGEXFFormat())) {
                        graphFormatsWritten += "gexf "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_GEXF, gexfFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_GEXF, new File(gexfGraphFile));
                    }
                }
                if(Settings.getBoolean("plcc_B_output_xml")) {
                    String xmlGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".xml";
                    xmlFileNoPath = fileNameWithoutExtension + ".xml";
                    if(IO.stringToTextFile(xmlGraphFile, pg.toXMLFormat())) {
                        graphFormatsWritten += "xml "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_XML, xmlFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_XML, new File(xmlGraphFile));
                    }
                    else {
                        DP.getInstance().w("Main", "Failed to write PG file in XML format.");
                    }
                }
                if(Settings.getBoolean("plcc_B_output_cytoscapejs")) {
                    String cytoscapejsGraphFile = filePathGraphs + fs + fileNameWithoutExtension + ".cyjs";
                    cytoscapejsFileNoPath = fileNameWithoutExtension + ".cyjs";
                    if(IO.stringToTextFile(cytoscapejsGraphFile, pg.toCytoscapeJSFormat())) {
                        graphFormatsWritten += "cyjs "; numFormatsWritten++; writtenFormatsDBFilesNoPath.put(GraphFormats.GRAPHFORMAT_CYTOSCAPEJS, cytoscapejsFileNoPath);
                        pcr.addProteinGraphOutputFile(gt, GraphFormats.GRAPHFORMAT_CYTOSCAPEJS, new File(cytoscapejsGraphFile));
                    }
                }
                
                
                /*
                Boolean jsonTest = false;
                if(jsonTest) {
                    String jsonFile = "graph_test.json";
                    if(IO.stringToTextFile(jsonFile, testJSONFormat())) {
                        System.out.println("Wrote json test file.");
                    }
                }
                */
                
                
                if(numFormatsWritten > 0) {
                    if(! (Settings.getBoolean("plcc_B_silent") || Settings.getBoolean("plcc_B_only_essential_output"))) {
                        System.out.println("      Exported protein ligand graph in " + numFormatsWritten + " formats (" + graphFormatsWritten + ") to '" + new File(filePathGraphs).getAbsolutePath() + fs + "'.");
                    }
                }
                
                imgFile = filePathImg + fs + fileNameWithExtension;
                String imgFileNoExt = filePathImg + fs + fileNameWithoutExtension;
                                
                
                // But we may need to write the graph to the database
                if(Settings.getBoolean("plcc_B_useDB")) {
                                        
                    
                    try { 
                        Boolean res;
                        if(Settings.getBoolean("plcc_B_write_graphstrings_to_database_pg")) {
                            res = DBManager.writeProteinGraphToDB(pdbid, chain, ProtGraphs.getGraphTypeCode(gt), pg.toGraphModellingLanguageFormat(), pg.toVPLGGraphFormat(), pg.toKavoshFormat(), pg.toDOTLanguageFormat(), pg.toJSONFormat(), pg.toXMLFormat(), pg.getSSEStringSequential(), pg.containsBetaBarrel());
                        }
                        else {                            
                            res = DBManager.writeProteinGraphToDB(pdbid, chain, ProtGraphs.getGraphTypeCode(gt), null, null, null, null, null, null, pg.getSSEStringSequential(), pg.containsBetaBarrel()); 
                        }
                        
                        if((! silent) && res) {
                            System.out.println("      Inserted '" + gt + "' graph of PDB ID '" + pdbid + "' chain '" + chain + "' into DB.");
                        }
                    }
                    catch(SQLException e) { 
                        DP.getInstance().e("Main", "Failed to insert '" + gt + "' graph of PDB ID '" + pdbid + "' chain '" + chain + "' into DB: '" + e.getMessage() + "'."); 
                    }
                    
                    // update GML path
                    if(writtenFormatsDBFilesNoPath.size() > 0) {
                        Long graphDBID = -1L;
                        try {
                            graphDBID = DBManager.getDBProteinGraphID(pdbid, chain, gt);
                        } catch(SQLException ex) {
                            DP.getInstance().e("Main", "Could not find graph in database to update GML file path: '" + ex.getMessage() + "'.");
                        }
                        if(graphDBID > 0) {                                                                                                               
                            
                            for(String format : writtenFormatsDBFilesNoPath.keySet()) {
                                String fileDBPath = writtenFormatsDBFilesNoPath.get(format);
                                
                                if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                                    fileDBPath = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, chain) + fs + fileDBPath;
                                }
                                
                                try {
                                    DBManager.updateProteinGraphTextformatPathInDB(graphDBID, format, IO.pathToWebPath(IO.stripTrailingShitFromPathIfThere(fileDBPath)));
                                } catch(SQLException ex) {
                                    DP.getInstance().e("Main", "Could not update format '" + format + "' file path of graph in database: '" + ex.getMessage() + "'.");
                                }
                            }                                                        
                        }
                    }
                    
                    // assign SSEs in database
                    try {
                        int numAssigned = DBManager.assignSSEsToProteinGraphInOrder(pg.getVertices(), pdbid, chain, ProtGraphs.getGraphTypeCode(gt));
                        if(! (Settings.getBoolean("plcc_B_silent") || Settings.getBoolean("plcc_B_only_essential_output"))) {
                            System.out.println("      Assigned " + numAssigned + " SSEs to " + gt + " graph of PDB ID '" + pdbid + "' chain '" + chain + "' in the DB.");
                        }
                    } catch(SQLException ex) {
                       DP.getInstance().e("Main", "Could not assign SSEs to graph in the database: '" + ex.getMessage() + "'.");
                    }
                }

                if(Settings.getBoolean("plcc_B_draw_graphs")) {
                    
                    IMAGEFORMAT[] formats;
                    // formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG, DrawTools.IMAGEFORMAT.PDF };                    
                    formats = Settings.getProteinGraphOutputImageFormats();

                    HashMap<IMAGEFORMAT, String> filesByFormatCurNotation = ProteinGraphDrawer.drawProteinGraph(imgFileNoExt, false, formats, pg, new HashMap<Integer, String>(), new ArrayList<String>());
                    //if(! silent) {
                    //    System.out.println("      Image of graph written to file '" + imgFile + "'.");
                    //}
                    
                    //pcr.addProteinGraphImageBitmap(gt, new File(imgFile));
                    for(IMAGEFORMAT f : filesByFormatCurNotation.keySet()) {
                        pcr.addProteinGraphOutputImage(gt, f.toString(), new File(filesByFormatCurNotation.get(f)));
                    }

                    // set image location in database if required
                    if(Settings.getBoolean("plcc_B_useDB")) {
                        Long graphDBID = -1L;
                        try {
                            graphDBID = DBManager.getDBProteinGraphID(pdbid, chain, gt);
                        } catch(SQLException ex) {
                            DP.getInstance().e("Main", "Could not find graph in database: '" + ex.getMessage() + "'.");
                        }
                        if(graphDBID > 0) {

                            String dbImagePath;
                            for(IMAGEFORMAT format : filesByFormatCurNotation.keySet()) {
                                dbImagePath = fileNameWithoutExtension + DrawTools.getFileExtensionForImageFormat(format);


                                if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                                    dbImagePath = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, chain) + fs + dbImagePath;
                                }
                                //DP.getInstance().d("dbImagePath is '" + dbImagePath + "'.");


                                try {
                                    DBManager.updateProteinGraphImagePathInDB(graphDBID, format, IO.pathToWebPath(dbImagePath));
                                } catch(SQLException e) {
                                    DP.getInstance().e("Main", "Could not update graph image path in database: '" + e.getMessage() + "'.");
                                }
                            }
                        } else {
                            DP.getInstance().e("Main", "Could not find " + gt + " graph for PDB " + pdbid + " chain " + chain + " in database to set image path.");
                        }
                        
                    }                    
                                     
                }
                else {
                    if(! silent) {
                        System.out.println("      Image and graph output disabled, not drawing and writing protein graph files.");
                    }
                }
                
                
                
                if(Settings.getInteger("plcc_I_debug_level") > 0) {
                    if(! silent) {
                        System.out.println("      Graph plus string is '" + pg.getGraphPlusString() + "'.");
                    }
                }
                
                // commands to draw the graph in 3D into the PDB coords in JMOL
                if(Settings.getBoolean("plcc_B_Jmol_graph_vis_commands")) {                                                            
                    String graphVisualizationFileJmolCommands = filePathGraphs + fs + fileNameWithoutExtension + ".jmol";
                    if(IO.stringToTextFile(graphVisualizationFileJmolCommands, JmolTools.visualizeGraphCommands(pg, true, true))) {
                        if(Settings.getBoolean("plcc_B_output_textfiles_dir_tree_html")) {
                            pcr.addProteinGraphVisJmolCommandFile(gt, new File(graphVisualizationFileJmolCommands));
                        }
                    }
                }
                
                // commands to color the SSEs of the graph blue in 3D in JMOL
                if(Settings.getBoolean("plcc_B_Jmol_graph_vis_resblue_commands")) {                                        
                    
                    String graphVisualizationResBlueFileJmolCommands = filePathGraphs + fs + fileNameWithoutExtension + "_resblue" + ".jmol";
                    if(IO.stringToTextFile(graphVisualizationResBlueFileJmolCommands, JmolTools.visualizeGraphSubsetSSEsInBlue(pg, pg.getVertices(), true, true))) {
                        if(Settings.getBoolean("plcc_B_output_textfiles_dir_tree_html")) {
                            pcr.addProteinGraphVisResBlueJmolCommandFile(gt, new File(graphVisualizationResBlueFileJmolCommands));
                        }
                    }
                }
                
                // ###TEST-PG-METRICS
                if(Settings.getBoolean("plcc_B_compute_graph_metrics")) {
                    
                    if(pg.getSize() > 0) {
                    
                        pg.computeConnectedComponents();
                        FoldingGraph fg = pg.getLargestConnectedComponent();                        

                        GraphProperties gp = new GraphProperties(pg);
                        GraphProperties sgp = new GraphProperties(fg);

                        if(Settings.getBoolean("plcc_B_useDB")) {

                            if( ! DBManager.getAutoCommit()) {
                                DBManager.commit();
                            }

                            try {
                                Long graph_db_id = DBManager.getDBProteinGraphID(pdbid, chain, gt);
                                if(graph_db_id > 0L) {
                                    //System.out.println("Found graph " + pdbid + " " + chain + " " + gt + " with ID " + graph_db_id + ".");
                                    // write graph properties
                                    Long runtime_secs = null;
                                    DBManager.writeProteingraphStatsToDB(graph_db_id, Boolean.FALSE, gp.getNumVertices(), gp.getNumEdges(), gp.getMinDegree(), gp.getMaxDegree(), gp.getConnectedComponents().size(), gp.getGraphDiameter(), gp.getGraphRadius(), gp.getAverageClusterCoefficient(), gp.getAverageShortestPathLength(), gp.getDegreeDistributionUpTo(50), gp.getAverageDegree(), gp.getDensity(), gp.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);
                                    // write properties of largest CC of graph
                                    DBManager.writeProteingraphStatsToDB(graph_db_id, Boolean.TRUE, sgp.getNumVertices(), sgp.getNumEdges(), sgp.getMinDegree(), sgp.getMaxDegree(), sgp.getConnectedComponents().size(), sgp.getGraphDiameter(), sgp.getGraphRadius(), sgp.getAverageClusterCoefficient(), sgp.getAverageShortestPathLength(), sgp.getDegreeDistributionUpTo(50), sgp.getAverageDegree(), sgp.getDensity(), sgp.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);
                                }
                                else {
                                    DP.getInstance().e("Main", "Could not write graph properties to DB, graph not found in database.");
                                }
                            } catch(SQLException e) {
                                DP.getInstance().e("SQL error while trying to store graph stats: '" + e.getMessage()+ "'.");
                            }
                        }
                    
                    }
                
                }
                
                
                /* ----------------------------------------------- Folding graphs ---------------------------------------------- */

                if(Settings.getBoolean("plcc_B_folding_graphs")) {
                    //if(gt.equals(ProtGraphs.GRAPHTYPE_STRING_ALPHA) || gt.equals(ProtGraphs.GRAPHTYPE_STRING_BETA) || gt.equals(ProtGraphs.GRAPHTYPE_STRING_ALBE)) {
                        
                        if( ! (Settings.getBoolean("plcc_B_silent") || Settings.getBoolean("plcc_B_only_essential_output"))) {
                            System.out.println("      Computing " + gt + " folding graphs.");
                        }
                        
                        //System.out.println("!!!!!!calling for path '" + filePathImg + "'.");
                        ProteinFoldingGraphResults fgRes = calculateFoldingGraphsForSSEGraph(pg, filePathImg);                                    
                        pcr.addProteinFoldingGraphResults(gt, fgRes);
                        
                        if (Settings.getBoolean("plcc_B_matrix_structure_search") && Settings.get("plcc_S_linear_notation_graph_type").equals(pg.getGraphType())){
                            
                            // turn the linear notation into an adjacencymatrix and search it in the adjacencymatrix of the protein
                            
                            PTGLNotations p = new PTGLNotations(pg);
                            p.stfu();
                            //p.adjverbose = true;
                            List<PTGLNotationFoldResult> resultsPTGLNotations = p.getResults(); //generate linear notation for proteingraph
                            
                            if (!silent){
                                System.out.println("          Linear notation of the folding graph: " + resultsPTGLNotations.get(0).adjNotation + " ---");
                            }
                            
                            String gt_new = gt;
                            //Changing graphtype, because the function parseRedOrAdjToMatrix doesn't need any information about ligands
                            switch(gt){
                                case "alphalig":
                                    gt_new = "alpha";
                                case "betalig":
                                    gt_new = "beta";
                                case "albelig":
                                    gt_new = "albe";
                            }

                            ArrayList<ArrayList<Character>> matrix = new ArrayList<>(); //the adjacencymatrix for the linear notation
                            matrix = DBManager.parseRedOrAdjToMatrix(resultsPTGLNotations.get(0).adjNotation, gt_new);
                            
                            //save the linear notation from the input in the adjacencymatrix "pattern"
                            ArrayList<ArrayList<Character>> pattern = new ArrayList<>(); 
                            
                            pattern = DBManager.parseRedOrAdjToMatrix(Settings.get("plcc_S_linear_notation"), Settings.get("plcc_S_linear_notation_graph_type"));
                            
                            if (pattern.size() <= matrix.size()){
                                //start searching
                                if (!silent){
                                    System.out.println("      --- Start searching the linear notation " + Settings.get("plcc_S_linear_notation") +" in the folding graph. ---");
                                }
                                int[] output_array = new int [2]; //saves the indexes in matrix, where the pattern was found
                                output_array = DBManager.matrixSearch(pattern, matrix);
                                
                                if (!silent && output_array[0] != -1){ //if the pattern wasn't found, output_array[0] = -1
                                    System.out.println("     **** Linear notation found at indexes (" + output_array[0] + ", " + output_array[1] + ") of the adjacency matrix from the folding graph. ****");
                                } 
                            }
                        }
                        
                    //} else {
                    //    if( ! silent) {
                    //        System.out.println("      Handling folding graphs, but skipping graph type '" + gt + "'.");
                    //    }
                    //}
                }
                else {
                    if( ! (Settings.getBoolean("plcc_B_silent") || Settings.getBoolean("plcc_B_only_essential_output"))) {
                        System.out.println("      Not handling folding graphs.");
                    }
                }
              
            }
            
            
            if(! silent) {
                System.out.println("  +++++ All " + graphTypes.size() + " protein graphs of chain " + c.getPdbChainID() + " handled. +++++");
            }
            
            if(Settings.getBoolean("plcc_B_useDB") && Settings.getBoolean("plcc_B_folding_graphs") && Settings.getBoolean("plcc_B_compute_motifs")) {
                Integer numAssigned = 0;
                try {
                    numAssigned = DBManager.checkAndAssignChainToAllMotifsInDatabase(pdbid, chain);
                    if(! silent) {
                        System.out.println("      Computed SSE motifs for chain " + chain + ", found " + numAssigned + " motifs in all folding graph linear notations.");
                    }
                } catch(Exception e) {
                    DP.getInstance().e("Main", "Computing SSE motifs failed for PDB " + pdbid + " chain " + chain + ": '" + e.getMessage() + "'.");
                    //e.printStackTrace();
                }
            }
            else {
                if(! silent) {
                    System.out.println("      Not computing any motifs for chain " + chain + " (disabled, requires FGs and database).");
                }
            }
            
            // Testing only
            //System.out.println("Chain allRes chemProps: " + c.getChainChemProps5StringAllResidues());
            //String[] str = c.getChainChemPropsStringSSEResiduesOnly(" ");
            //System.out.println("Chain SSE chemProps: " + str[0]);
            //System.out.println("Chain SSE sourceSSE: " + str[1]);
            
        }
        
        DBManager.commit();
        
        if( ! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
            System.out.print("The following macromolecules exist in the PDB file (format: MOL_ID(list of chains) ...):");
            for(String key : macroMoleculesOfPDBfileToChains.keySet()) {
                System.out.print(" " + key + "(");
                System.out.print(IO.stringListToString(macroMoleculesOfPDBfileToChains.get(key)) + ")");
            }
            System.out.print("\n");
        }
        
        // we may need to write the macromolecules to the DB
        if(Settings.getBoolean("plcc_B_useDB")) {
            for(String molID : macroMolecules.keySet()) {
                Map<String, String> mol = macroMolecules.get(molID);
                //String pdb_id, String molIDPDBfile, String molName, String molECNumber, String orgScientific, String orgCommon
                //tmpMacroMol.put("pdb_mol_id", pmi.getMacromolID()); // not strictly needed, it is also put as the key for this MM later
                //tmpMacroMol.put("pdb_mol_name", pmi.getMolName());
                //tmpMacroMol.put("pdb_org_sci", pmi.getOrgScientific());
                //tmpMacroMol.put("pdb_org_common", pmi.getOrgCommon());
                //tmpMacroMol.put("pdb_all_chains", pmi.getAllMolChains());
                //tmpMacroMol.put("pdb_ec_number", pmi.getECNumber());
                try {
                    DBManager.writeMacromoleculeToDB(pdbid, mol.get("pdb_mol_id"), mol.get("pdb_mol_name"), mol.get("pdb_ec_number"), mol.get("pdb_org_sci"), mol.get("pdb_org_common"), mol.get("pdb_all_chains"));
                    if( ! silent) {
                        System.out.println("  Writing macromolecule '" + mol.get("pdb_mol_name") + "' with PDB MOL_ID '" + mol.get("pdb_mol_id") + "' to database. Consists of chains '" + mol.get("pdb_all_chains") + "'.");
                    }
                } catch(Exception e) {
                    DP.getInstance().e("Main", "Failed to write macromolecule to database: '" + e.getMessage() + "'.");
                }
            }
            
            DBManager.commit();
            
            // now assign the chains to their respective macromolecule
            int numAssigned = 0;
            for(Integer i = 0; i < allChains.size(); i++) {
                c = allChains.get(i);
                String chainName = c.getPdbChainID();
                Long chain_db_id = DBManager.getDBChainID(pdbid, chainName);
                Long mm_db_id = DBManager.getDBMacromoleculeID(pdbid, c.getMacromolID());
                if(chain_db_id >= 1L && mm_db_id >= 1L) {
                    try {
                        DBManager.assignChainToMacromolecule(chain_db_id, mm_db_id);
                        numAssigned++;
                    }catch(SQLException e) {
                        DP.getInstance().e("Main", "SQL error while assigning a chain to a macromolecule: '" + e.getMessage() + "'.");
                    }
                }
                else {
                    DP.getInstance().e("Main", "Assigning chain '" + chainName + "' to macromolecule with PDB file MOL_ID '" + c.getMacromolID() + "' failed. Chain or MM not in database.");
                }
            }
            
            if( ! silent) {
                System.out.println("  Assigned the " + numAssigned + " chains in the PDB data to the " + macroMolecules.keySet().size() + " different macromolecules.");
            }
            
            DBManager.commit();
            
        }
        
        // Calculate Complex Graph
        if(Settings.getBoolean("plcc_B_complex_graphs")) {
            // calculate ALBELIG CG           
            calculateComplexGraph(allChains, resList, resContacts, pdbid, outputDir, SSEGraph.GRAPHTYPE_ALBELIG);
            
            // test: also calculate ALBE CG
            //calculateComplexGraph(allChains, resList, resContacts, pdbid, outputDir, SSEGraph.GRAPHTYPE_ALBE);
        }
        
        
        
        
        //System.out.println("All " + allChains.size() + " chains done.");                
    }
    
    
    
    /**
     * Calculates all requested folding graphs for the protein graph (or 'SSE graph') pg. Which graphs are drawn is determined by the 
     * setting on the command line / configuration file.
     * @param pg the protein graphs
     * @param outputDir the file system path where to write the image files. Has to exist and be writable.
     * @return the protein folding graph results, which gives access to the graphs and output files
     */
    public static ProteinFoldingGraphResults calculateFoldingGraphsForSSEGraph(ProtGraph pg, String outputDir) {
        //System.out.println("Searching connected components in " + graphType + " graph of chainName " + c.getPdbChainID() + ".");
        boolean silent = Settings.getBoolean("plcc_B_silent");
        boolean essentialOutOnly = Settings.getBoolean("plcc_B_only_essential_output");
        //ArrayList<FoldingGraphComputationResult> fgcs = pg.getFoldingGraphComputationResults();
        //ArrayList<FoldingGraph> ccs = pg.getConnectedComponents();
        //ArrayList<FoldingGraph> foldingGraphs = FoldingGraphComputationResult.getFoldingGraphsREDandKEYFromFGCR(fgcs);
        ArrayList<FoldingGraph> foldingGraphs = pg.getConnectedComponents();
        Collections.sort(foldingGraphs, new FoldingGraphComparator());
        
        PTGLNotations p = new PTGLNotations(pg);
        p.stfu();
        //p.adjverbose = true;
        List<PTGLNotationFoldResult> resultsPTGLNotations = p.getResults();
        
                
        HashMap<Integer, FoldingGraph> ccsList = new HashMap<Integer, FoldingGraph>();
        int fgMinSizeDraw = Settings.getInteger("plcc_I_min_fgraph_size_draw");
        int numFGsWithMinSize = 0;
        FoldingGraph fg = null;           // A connected component of a protein graph is a folding graph
        for (int i = 0; i < foldingGraphs.size(); i++) {
            fg = foldingGraphs.get(i);
            ccsList.put(i, fg);
            if(fg.getSize() >= fgMinSizeDraw) {
                numFGsWithMinSize++;
            }
        }        
        ProteinFoldingGraphResults fgRes = new ProteinFoldingGraphResults(ccsList);
        fg = null;
        
        String fgFile = null;
        String fs = System.getProperty("file.separator");

        //System.out.println("Found " + ccs.size() + " connected components in " + graphType + " graph of chainName " + c.getPdbChainID() + ".");
        if(! silent) {
            System.out.println("      --- Handling all " + foldingGraphs.size() + " " + pg.getGraphType() + " folding graphs of the " + pg.getGraphType() + " protein graph (" + numFGsWithMinSize + " with >= " + fgMinSizeDraw + " verts) ---");
        }
        for(Integer j = 0; j < foldingGraphs.size(); j++) {            
            fg = foldingGraphs.get(j);
            Integer fg_number = fg.getFoldingGraphNumber();
            String pdbid = fg.getPdbid();
            String chain = fg.getChainid();
            String gt = fg.getGraphType();

            PTGLNotationFoldResult pnfr = resultsPTGLNotations.get(j);
            if(!Objects.equals(pnfr.getFoldNumber(), fg_number)) {
                DP.getInstance().e("Main", "calculateFoldingGraphsForSSEGraph(): fg_number of PTGLNotationFoldResult does not match current fg_number.");
                Main.doExit(1);
            }
            
            
            
            
            // graph strings in GML format and others, does NOT include PTGL linear notations
            if(fg.numAlphaBetaVertices() >= Settings.getInteger("plcc_I_min_fgraph_size_write_to_file")) {
                writeFGGraphStrings(fg, outputDir, fg.getFoldingGraphNumber());
                if(Settings.getBoolean("plcc_B_output_fg_linear_notations_to_file")) {
                    writeFGLinearNotationStrings(pnfr.getFoldingGraph(), outputDir, pnfr.getFoldNumber(), pnfr);
                }                     
            }
                
            // Draw all folding graphs in all notations
            //List<String> notations = Arrays.asList("KEY", "ADJ", "RED", "SEQ");
            ArrayList<String> notations = new ArrayList<String>();

            if(Settings.getBoolean("plcc_B_foldgraphtype_KEY")) { notations.add(FoldingGraph.FG_NOTATION_KEY); }
            if(Settings.getBoolean("plcc_B_foldgraphtype_ADJ")) { notations.add(FoldingGraph.FG_NOTATION_ADJ); }
            if(Settings.getBoolean("plcc_B_foldgraphtype_RED")) { notations.add(FoldingGraph.FG_NOTATION_RED); }
            if(Settings.getBoolean("plcc_B_foldgraphtype_SEQ")) { notations.add(FoldingGraph.FG_NOTATION_SEQ); }                                                                        
            if(Settings.getBoolean("plcc_B_foldgraphtype_SEQ")) { notations.add(FoldingGraph.FG_NOTATION_DEF); }
            
            // We may need to write the folding graph to the database
            Long fgDbId = -1L;
            if(Settings.getBoolean("plcc_B_useDB")) { 
                
                if(fg.numAlphaBetaVertices() < Settings.getInteger("plcc_I_min_fgraph_size_write_to_db")) {
                    if(! (silent || essentialOutOnly)) {
                        System.out.println("       *Not writing " + gt + " folding graph #" + j + " of size " + fg.numVertices() + " to DB (minimum size is " + Settings.getInteger("plcc_I_min_fgraph_size_write_to_db") + ").");
                    }                   
                }
                else {
                                
                    try { 
                        if(Settings.getBoolean("plcc_B_write_graphstrings_to_database_fg")) {
                            fgDbId = DBManager.writeFoldingGraphToDB(pdbid, chain, ProtGraphs.getGraphTypeCode(gt), fg_number, FoldingGraph.getFoldNameOfFoldNumber(fg_number), fg.getMinimalVertexIndexInParentGraph(), fg.toGraphModellingLanguageFormat(), fg.toVPLGGraphFormat(), fg.toKavoshFormat(), fg.toDOTLanguageFormat(), fg.toJSONFormat(), fg.toXMLFormat(), fg.getSSEStringSequential(), fg.containsBetaBarrel()); 
                        }
                        else {
                            fgDbId = DBManager.writeFoldingGraphToDB(pdbid, chain, ProtGraphs.getGraphTypeCode(gt), fg_number, FoldingGraph.getFoldNameOfFoldNumber(fg_number), fg.getMinimalVertexIndexInParentGraph(), null, null, null, null, null, null, fg.getSSEStringSequential(), fg.containsBetaBarrel()); 
                        }
                        
                        if(! (silent || essentialOutOnly)) {
                            System.out.println("        Inserted '" + gt + "' folding graph # " + fg_number + " of PDB ID '" + pdbid + "' chain '" + chain + "' into DB.");
                        }
                    }
                    catch(SQLException e) { 
                        DP.getInstance().e("Main", "Failed to insert '" + gt + "' folding graph # " + fg_number + " of PDB ID '" + pdbid + "' chain '" + chain + "' into DB: '" + e.getMessage() + "'."); 
                    }

                    // assign SSEs in database
                    try {
                        fgDbId = DBManager.getDBFoldingGraphID(pdbid, chain, gt, fg_number);
                    } catch(SQLException sqlex) {
                        DP.getInstance().e("Main", "Folding graph #" + fg_number + " not found in DB: '" + sqlex.getMessage() + "'.");
                        fgDbId = -1L;
                    }
                    if(fgDbId >= 1) {
                        try {
                            Integer[] numAssigned = DBManager.assignSSEsToFoldingGraphInOrderWithSecondat(fg.getVertices(), fgDbId, gt, fg_number, FoldingGraph.getFoldNameOfFoldNumber(fg_number));
                            if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                                if(Objects.equals(numAssigned[0], numAssigned[1])) {
                                    System.out.println("        Assigned " + numAssigned[0] + " SSEs to " + gt + " folding graph # " + fg_number + " of PDB ID '" + pdbid + "' chain '" + chain + "' in the DB.");
                                } else {
                                    DP.getInstance().w("Main", "SSE to folding graph assignment: " + numAssigned[0] + " SSEs assigned to fg # " + fg_number + " in sse2fg table, but " + numAssigned[1] + " in secondat table -- expected same value.");
                                }
                            }


                        } catch(SQLException ex) {
                           DP.getInstance().e("Main", "Could not assign SSEs to graph in the database: '" + ex.getMessage() + "'.");
                        }
                    } else {
                        DP.getInstance().e("Main", "Cannot assign SSEs to folding graph #" + fg_number + ", FG not found in DB.");
                    }


                    Long insertIDSingleTable = -1L;
                    try {
                        //insertID = DBManager.writeFoldLinearNotationsToDatabaseMultiTable(fg.getPdbid(), fg.getChainid(), fg.getGraphType(), fg.getFoldingGraphNumber(), pnfr);
                        insertIDSingleTable = DBManager.writeFoldLinearNotationsToDatabaseSingleTable(fg.getPdbid(), fg.getChainid(), fg.getGraphType(), fg.getFoldingGraphNumber(), pnfr);
                    } catch (SQLException ex) {
                        DP.getInstance().e("Main", "Could not write " + fg.getGraphType() + " FG # " + fg.getFoldingGraphNumber() + " linear notations to DB: '" + ex.getMessage() + "'.");
                    }


                    if(insertIDSingleTable < 0) {
                        DP.getInstance().e("Main", "Could not write " + fg.getGraphType() + " FG # " + fg.getFoldingGraphNumber() + " linear notations to single linnot table in DB, empty insert ID.");                        
                    }
                    else {
                        if(! (silent || essentialOutOnly)) {
                            System.out.println("        Wrote PTGL linear notations for " + fg.getGraphType() + " FG # " + fg.getFoldingGraphNumber() + " to database.");
                        }
                    } 
                
                }
                
            }
            
            // draw folding graphs                                   
            if(Settings.getBoolean("plcc_B_draw_folding_graphs")) {
                
                if(fg.numAlphaBetaVertices() >= Settings.getInteger("plcc_I_min_fgraph_size_draw")) {
                
                    if(! (silent || essentialOutOnly)) {
                        System.out.println("        Drawing all " + notations.size() + " notations of the " + pg.getGraphType() + " FG #" + j + " of size " + fg.getSize() + ".");
                    }
                
                    //System.out.println("          At FG #" + j + "(fg_number=" + fg_number + ").");
                    for(String notation : notations) {                                                

                        String fileNameWithoutExtension = pg.getPdbid() + "_" + pg.getChainid() + "_" + pg.getGraphType() + "_FG_" + j + "_" + notation;
                        //String fileNameWithExtension = fileNameWithoutExtension + ".png";
                        
                        
                        fgFile = outputDir + System.getProperty("file.separator") + fileNameWithoutExtension; //Settings.get("plcc_S_img_output_fileext");

                        Boolean drawingSucceeded = false;
                        IMAGEFORMAT[] formats = Settings.getFoldingGraphOutputImageFormats();
                        HashMap<IMAGEFORMAT, String> filesByFormatCurNotation = new HashMap<>();

                        if(notation.equals(FoldingGraph.FG_NOTATION_ADJ)) {     
                            filesByFormatCurNotation = ProteinGraphDrawer.drawFoldingGraphADJ(fgFile, false, formats, pnfr);                        
                        }
                        else if(notation.equals(FoldingGraph.FG_NOTATION_RED)) {
                            filesByFormatCurNotation = ProteinGraphDrawer.drawFoldingGraphRED(fgFile, false, formats, pnfr);                                                
                        }
                        else if(notation.equals(FoldingGraph.FG_NOTATION_SEQ)) {
                            filesByFormatCurNotation = ProteinGraphDrawer.drawFoldingGraphSEQ(fgFile, false, formats, pnfr);                                                
                        }
                        else if(notation.equals(FoldingGraph.FG_NOTATION_KEY)) {                            
                            filesByFormatCurNotation = ProteinGraphDrawer.drawFoldingGraphKEY(fgFile, false, formats, pnfr);                                                                            
                        }
                        else if(notation.equals(FoldingGraph.FG_NOTATION_DEF)) {                            
                            filesByFormatCurNotation = ProteinGraphDrawer.drawFoldingGraphDEF(fgFile, false, formats, pnfr);                                                                            
                        }

                        drawingSucceeded = ( ! filesByFormatCurNotation.isEmpty());                                                

                        //if(fg.drawFoldingGraph(notation, fgFile)) {
                        if(drawingSucceeded) {
                            if(! silent) {
                                //System.out.println("         -Folding graph #" + j + " of the " + pg.getGraphType() + " graph of chainName " + pg.getChainid() + " written to file '" + fgFile + "' in " + notation + " notation.");
                            }

                            // save image path to database if required
                            if(Settings.getBoolean("plcc_B_useDB")) {
                                
                                //DP.getInstance().d("dbImagePath is '" + dbImagePath + "'.");                            

                                                              
                                for(IMAGEFORMAT format : filesByFormatCurNotation.keySet()) {
                                    String dbImagePath = fileNameWithoutExtension;
                                    if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                                        dbImagePath = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, chain) + fs + fileNameWithoutExtension;
                                    }
                                    
                                    dbImagePath += DrawTools.getFileExtensionForImageFormat(format);
                                    
                                    int numAff = 0;
                                    try {
                                        numAff = DBManager.updateFoldingGraphImagePathInDB(fgDbId, format, notation, IO.pathToWebPath(dbImagePath));
                                    } catch(SQLException e) {
                                        DP.getInstance().e("Main", "Could not update format " + format + " folding graph image path in database: '" + e.getMessage() + "'.");
                                    }
                                    
                                    if(numAff == 0) {
                                        DP.getInstance().e("Main", "Could not update format " + format + " folding graph image path in database, 0 rows affected.");
                                    }
                                    else {
                                        if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
                                            System.out.println("          Updated FG " + notation + " notation " + format + " format image path in database.");
                                        }
                                    }
                                }
                                
                                /*
                                try {
                                    DBManager.updateFoldingGraphImagePathInDB(fgDbId, DrawTools.IMAGEFORMAT.PNG, notation, dbImagePath);
                                } catch(SQLException e) {
                                    DP.getInstance().e("Main", "Could not update " + notation + " notation folding graph image path in database: '" + e.getMessage() + "'.");
                                }
                                */

                            }                                                
                        }
                        else {
                            if(Settings.getInteger("plcc_I_debug_level") > 0) {
                                System.err.println("NOTE: Could not draw notation " + notation + " of folding graph #" + j + " of the " + pg.getGraphType() + " graph of chain " + pg.getChainid() + ". (Tried to write to file '" + fgFile + "'.)");
                            }
                        }
                    } // notations
                    
                } else {
                    // FG too small
                    if(! (silent  || essentialOutOnly)) {
                        System.out.println("        Not drawing any of the " + notations.size() + " notations of the " + pg.getGraphType() + " FG #" + j + " (fg_number=" + fg_number + "): albe size is " + fg.numAlphaBetaVertices() +", min is " + Settings.getInteger("plcc_I_min_fgraph_size_draw") + ".");
                    }
                }
                
            }
            else {
                if(! silent) {
                    System.out.println("        Not drawing folding graphs, disabled.");
                }
            }
        }    
        return fgRes;
    }
    

    /**
     * Writes the folding graph strings to files on the HDD. This does NOT include the PTGL notations, only GML etc.
     * @param fg the folding graph
     * @param outputDir the directory where to write the files
     * @param fgNumber the CC number of this FG (within the parent PG)
     */
    public static void writeFGGraphStrings(FoldingGraph fg, String outputDir, int fgNumber) {
        
        boolean silent = Settings.getBoolean("plcc_B_silent");
        
        if(! (silent || Settings.getBoolean("plcc_B_only_essential_output"))) {
            System.out.println("       *Handling " + fg.getGraphType() + " folding Graph #" + fgNumber + " containing " + fg.numVertices() + " vertices and " + fg.numEdges() + " edges (" + fg.numSSEContacts() + " SSE contacts).");
        }
        
        String fs = File.separator;
        String fileNameWithoutExtension = fg.getPdbid() + "_" + fg.getChainid() + "_" + fg.getGraphType() + "_FG_" + fgNumber;
        String graphFormatsWritten = "";
        Integer numFormatsWritten = 0;
        if(Settings.getBoolean("plcc_B_output_GML")) {
            String gmlfFile = outputDir + fs + fileNameWithoutExtension + ".gml";
            if(IO.stringToTextFile(gmlfFile, fg.toGraphModellingLanguageFormat())) {
                graphFormatsWritten += "gml "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_TGF")) {
            String tgfFile = outputDir + fs + fileNameWithoutExtension + ".tgf";
            if(IO.stringToTextFile(tgfFile, fg.toTrivialGraphFormat())) {
                graphFormatsWritten += "tgf "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_DOT")) {
            String dotLangFile = outputDir + fs + fileNameWithoutExtension + ".gv";
            if(IO.stringToTextFile(dotLangFile, fg.toDOTLanguageFormat())) {
                graphFormatsWritten += "gv "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_kavosh")) {
            String kavoshFile = outputDir + fs + fileNameWithoutExtension + ".kavosh";
            if(IO.stringToTextFile(kavoshFile, fg.toKavoshFormat())) {
                graphFormatsWritten += "kavosh "; numFormatsWritten++;
            }
        }
        // write the SSE info text file for the image (plcc graph format file)
        if(Settings.getBoolean("plcc_B_output_plcc")) {
            String plccGraphFile = outputDir + fs + fileNameWithoutExtension + ".plg";
            if(IO.stringToTextFile(plccGraphFile, fg.toVPLGGraphFormat())) {
                graphFormatsWritten += "plg "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_json")) {
            String jsonGraphFile = outputDir + fs + fileNameWithoutExtension + ".json";
            if(IO.stringToTextFile(jsonGraphFile, fg.toJSONFormat())) {
                graphFormatsWritten += "json "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_gexf")) {
            String gexfGraphFile = outputDir + fs + fileNameWithoutExtension + ".gexf";
            if(IO.stringToTextFile(gexfGraphFile, fg.toGEXFFormat())) {
                graphFormatsWritten += "gexf "; numFormatsWritten++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_xml")) {
            String xmlGraphFile = outputDir + fs + fileNameWithoutExtension + ".xml";
            if(IO.stringToTextFile(xmlGraphFile, fg.toXMLFormat())) {
                graphFormatsWritten += "xml "; numFormatsWritten++;
            }
        }


        if(numFormatsWritten > 0) {
            if(! (silent  || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.println("        Exported folding graph #" + fgNumber + " in " + numFormatsWritten + " formats (" + graphFormatsWritten + ") to '" + new File(outputDir).getAbsolutePath() + fs + "'.");
            }
        }
    }
    
    
    /**
     * Writes the PTGL linear notations strings.
     * @param fg the folding graph (a PG connected component)
     * @param outputDir the output base directory
     * @param fgNumber the number of the FG (in the parent graph connected component ordering, see FoldingGraphComparator class for details). You can use fg.getFoldingGraphNumber() to determine this.
     * @param pnfr the PTGL linear notation result for this fold/connected component/FG. The PTGLNotations class computes this.
     */
    public static void writeFGLinearNotationStrings(FoldingGraph fg, String outputDir, int fgNumber, PTGLNotationFoldResult pnfr) {
        
        boolean silent = Settings.getBoolean("plcc_B_silent");
        
        String fs = File.separator;
        String fileNameWithoutExtension = fg.getPdbid() + "_" + fg.getChainid() + "_" + fg.getGraphType() + "_FG_" + fgNumber;
        Integer numFormatsWritten = 0;
        StringBuilder sb = new StringBuilder();
        
        sb.append("# line format is '<notation_type>:<start_vertex_index_in_parent_protein_graph>:<number_of_sses_in_connected_component>:<notation_string>'\n");
        sb.append("ADJ:").append(pnfr.adjStart).append(":").append(pnfr.adjSize).append(":").append(pnfr.adjNotation).append("\n");
        sb.append("RED:").append(pnfr.redStart).append(":").append(pnfr.redSize).append(":").append(pnfr.redNotation).append("\n");
        sb.append("KEY:").append(pnfr.keyStartFG).append(":").append(pnfr.keySize).append(":").append(pnfr.keyNotation).append("\n");
        sb.append("SEQ:").append(pnfr.seqStart).append(":").append(pnfr.seqSize).append(":").append(pnfr.seqNotation).append("\n");
        
        String linearNotationsFile = outputDir + fs + fileNameWithoutExtension + ".ptgllinnot";
        if(IO.stringToTextFile(linearNotationsFile, sb.toString())) {
            numFormatsWritten++;
        }
        
        if(numFormatsWritten > 0) {
            if(! (silent  || Settings.getBoolean("plcc_B_only_essential_output"))) {
                System.out.println("        Exported linear notations of " + fg.getGraphType() + " folding graph #" + fgNumber + " in PTGL format to dir '" + new File(outputDir).getAbsolutePath() + fs + "'.");
            }
        }
        else {
            if(! silent) {
                System.err.println("        Could not export linear notations of " + fg.getGraphType() + " folding graph #" + fgNumber + " in PTGL format to dir '" + new File(outputDir).getAbsolutePath() + fs + "'.");
            }
        }
    }
    
    /**
     * Determines the PTGL SSE ID for the SSE at sequential position seqNum in the primary sequence.
     *
     * @param seqNum the sequential SSE number
     * @return The PTGL SSE ID (e.g. "A" for the SSE# 1, "D" for SSE #4).
     */
    public static String getPtglSseIDForNum(Integer seqNum) {
        String alph = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        if(seqNum < alph.length()) {
            return("" + alph.charAt(seqNum));
        }
        else {
            return("" + alph.charAt(seqNum % alph.length()));
        }
    }
    
    /*
    public static String testJSONFormat() {
        Gson gson = new Gson();
        String json = gson.toJson(new UndirectedGraph());  
        return json;
    }
    */
    

    
    /**
     * Writes the String 'data' into the file  at 'filePath', overwriting the file if it exists
     * and creating it if not.
     * @param filePath the target file system location for the file
     * @param data the String to write to the file
     * @return true if it worked out, false otherwise
     */
    public static Boolean writeStringToFile(String filePath, String data) {

        FileWriter fw = null;
        PrintWriter pw = null;
        Boolean res = false;

        try {
            fw = new FileWriter(filePath);
            pw = new PrintWriter(fw);
            pw.print(data);
            pw.close();
            res = true;
        }
        catch (Exception e) {
            DP.getInstance().w("Could not write to file '" + filePath + "'.");
            res = false;
        }
       
        try {
            fw.close();
        } catch(Exception ex) {
            DP.getInstance().w("Could not close FileWriter for file '" + filePath + "'.");
        }
        return(res);
    }

    
    /**
     * Calculates the requested type of protein graph for the given list of SSEs (usually all SSEss of a chainName) and residue contacts.
     * @param graphType the requested graph type
     * @param allChainSSEs the SSEs of a chainName
     * @param c the chainName
     * @param resContacts a list of residue contacts (between residues of c)
     * @param pdbid the PDBID of the protein the chainName c belongs to
     * @return the resulting protein graph
     */
    public static ProtGraph calcGraphType(String graphType, List<SSE> allChainSSEs, Chain c, List<MolContactInfo> resContacts, String pdbid) {

        ContactMatrix chainCM;
        Boolean silent = Settings.getBoolean("plcc_B_silent");

        if(! silent) {
            System.out.println("    ----- Calculating " + graphType + " protein graph of chain " + c.getPdbChainID() + ". -----");
        }

        List<String> keepSSEs = new ArrayList<String>();
        List<SSE> filteredChainSSEs;

        // Check whether coils should be kept
        if(Settings.getBoolean("plcc_B_include_coils")) {
            keepSSEs.add(Settings.get("plcc_S_coilSSECode"));
        }

        // Filter SSEs depending on the requested graph type
        if(graphType.equals("albe")) {
            keepSSEs.add("E");
            keepSSEs.add("H");
        }
        else if(graphType.equals("alpha")) {
            keepSSEs.add("H");
        }
        else if(graphType.equals("beta")) {
            keepSSEs.add("E");            
        }
        else if(graphType.equals("albelig")) {
            keepSSEs.add("E");
            keepSSEs.add("H");
            keepSSEs.add("L");
        }
        else if(graphType.equals("alphalig")) {
            keepSSEs.add("H");
            keepSSEs.add("L");
        }
        else if(graphType.equals("betalig")) {
            keepSSEs.add("E");
            keepSSEs.add("L");
        }
        else {
            System.err.println("ERROR: calcGraphType(): Graph type '" + graphType + "' invalid.");            
            Main.doExit(1);
        }

        // Filters have been configured, now do the actual filtering.
        filteredChainSSEs = filterAllSSEsButList(allChainSSEs, keepSSEs);

        // The SSEs should already be ordered anyway
        Collections.sort(filteredChainSSEs, new SSEComparator());
        
        // SSE list has been filtered, let's go

        // Calculate SSE level contacts
        chainCM = new ContactMatrix(filteredChainSSEs, pdbid);
        chainCM.restrictToChain(c.getPdbChainID());
        chainCM.fillFromContactList(resContacts, keepSSEs);
        //chainCM.printTotalContactMatrix("TT");
        chainCM.calculateSSEContactMatrix();
        //chainCM.printSSEMatrix();
        
        Boolean computeAll = false;
        
        chainCM.calculateSSESpatialRelationMatrix(resContacts, computeAll);
        //chainCM.printSpatialRelationMatrix();
        //chainCM.printResContMatrix();
        //chainCM.printTotalContactMatrix("TT");
        
        if(Settings.getBoolean("plcc_B_ptgl_geodat_output")) {
            String gdf = Settings.get("plcc_S_output_dir") + System.getProperty("file.separator") + pdbid + "_" + c.getPdbChainID() + "_" + graphType + ".geodat";            
            if(writeStringToFile(gdf, chainCM.toGeodatFormat(false, true))) {
                if(! silent) {
                    System.out.println("  Wrote SSE level contacts for chain " + chainCM.getChain() + " in geo.dat format to file '" + gdf + "'.");
                }
            }
            else {
                DP.getInstance().w("Failed to write SSE level contacts in geo.dat format to file '" + gdf + "'. Check permissions.");
            }
        }
        

        //UndirectedGraph<String, LabeledEdge> g = chainCM.toJGraph();

        // We only write the SSE contacts for the albelig graph because it contains all SSEs we are interested in.
        //  Writing them for all makes them appear multiple times.
        if(Settings.getBoolean("plcc_B_useDB") && graphType.equals("albelig")) {
            
            if(Settings.getBoolean("plcc_B_db_use_batch_inserts")) {
                chainCM.batchWriteContactStatisticsToDB();
            } else {
                chainCM.writeContactStatisticsToDB();
            }
        }


        ProtGraph pg = chainCM.toProtGraph();
        pg.declareProteinGraph();

        
        if(Settings.getBoolean("plcc_B_forceBackboneContacts")) {
            if(! silent) {
                System.out.println("      Adding backbone contacts to consecutive SSEs of the " + graphType + " graph.");
            }
            pg.addFullBackboneContacts();            
        }

        //System.out.println("    ----- Done with " + graphType + " graph of chainName " + c.getPdbChainID() + ". -----");
        return(pg);

    }


    
    /**
     * Filters a list of SSEs by type, returning a list containing all SSEs that have of the SSE types defined
     * in keepSSEs.
     * @param sses the input list
     * @param keepSSEs the list of types to keep
     * @return A filtered list of SSEs.
     */
    public static List<SSE> filterAllSSEsButList(List<SSE> sses, List<String> keepSSEs) {
        
        ArrayList<SSE> kept = new ArrayList<SSE>();

        if(keepSSEs.size() < 1) {
            DP.getInstance().w("filterAllSSEsButList(): list keepSSEs is empty, removing *all* SSEs.");
        }

        for(Integer i = 0; i < sses.size(); i++) {
            if(keepSSEs.contains(sses.get(i).getSseType())) {
                kept.add(sses.get(i));
            }
        }

        return(kept);
    }


    /**
     * Prints short info on how to get help on command line options and exits the program.
     */
    public static void syntaxError() {
        System.err.println("ERROR: Invalid command line. Use '-h' or --help' for info on how to run this program.");
        Main.doExit(1);
    }
    
    /**
     * Prints short info on how to get help on command line options and exits the program.
     */
    public static void syntaxError(String hint) {
        System.err.println("ERROR: Invalid command line. Use '-h' or --help' for info on how to run this program.");
        System.err.println("ERROR: Hint: '" + hint + "'");
        Main.doExit(1);
    }

    /**
     * Writes the a priori frequencies of the amino acids into the contact[][][][] array.
     * @param res the residue list which is used to determine the AA frequencies
     */
    public static void getAADistribution(List<Residue> res) {

        Residue r = null;
        Integer rID = -1;

        for(Integer i = 0; i < res.size(); i++) {
            r = res.get(i);
            rID = r.getInternalAAID();

            contact[rID][0][0][0]++;

        }
    }

    /**
     * Computes residue contacts following the definition by Andreas for interfaces between chains of multi-chain proteins.
     * @param res a list of residues
     * @return A list of MolContactInfo objects, each representing a pair of residues that are in contact.
     */
    public static ArrayList<MolContactInfo> calculateAllContactsAlternativeModel(List<Residue> res) {

        System.out.println("\n *** Calculation of interchain contacts. Still WIP! ***");
        FileParser.silent = false;
        Residue a, b;
        Integer rs = res.size();
        
        // jnw_2019: switch rna off for alternative model
        if (Settings.getBoolean("plcc_B_include_rna")) {
            DP.getInstance().w("Inclusion of RNA not implemented for alternative contacts model. Ignoring RNA.");
        }
        
        Integer numResContactsChecked, numResContactsPossible, numResContactsImpossible;
        numResContactsChecked = numResContactsPossible = numResContactsImpossible = 0;

        MolContactInfo rci;
        ArrayList<MolContactInfo> contactInfo = new ArrayList<MolContactInfo>();
        

        for(Integer i = 0; i < rs; i++) {
            a = res.get(i);
            
            //DEBUG
            //System.out.println("calculatePiEffects: Residue " + a.getFancyName() + " has " + a.getHydrogenAtoms().size() + " H atoms.");

            for(Integer j = i + 1; j < rs; j++) {
                b = res.get(j);
                numResContactsChecked++;

                
                // We only need to check on atom level if the center spheres overlap
                //if(a.contactPossibleWithResidue(b)) {                                        
                if(a.interchainContactPossibleWithResidue(b)) {                                        
                    numResContactsPossible++;
                    
                    rci = calculateAtomContactsBetweenResiduesAlternativeModel(a, b);
                    
                    if( rci != null) {
                        // There were atoms contacts!
                        contactInfo.add(rci);
                        /* System.out.println("#!# BBHB " + rci.getNumContactsBBHB());
                        System.out.println("#!# BBBH " + rci.getNumContactsBBBH());
                        System.out.println("#!# IVDW " + rci.getNumContactsIVDW());
                        System.out.println("#!# ISS " + rci.getNumContactsISS());
                        System.out.println("#!# BCHB " + rci.getNumContactsBCHB());
                        System.out.println("#!# BCBH " + rci.getNumContactsBCBH());
                        System.out.println("#!# CBHB " + rci.getNumContactsCBHB());
                        System.out.println("#!# CBBH " + rci.getNumContactsCBBH());
                        System.out.println("#!# CCHB " + rci.getNumContactsCCHB());
                        System.out.println("#!# CCBH " + rci.getNumContactsCCBH());
                        
                        System.out.println("#!# BB " + rci.getNumContactsBB());
                        System.out.println("#!# CB " + rci.getNumContactsCB());
                        System.out.println("#!# BC " + rci.getNumContactsBC());
                        System.out.println("#!# CC " + rci.getNumContactsCC());
                        System.out.println("#!# BL " + rci.getNumContactsBL());
                        System.out.println("#!# LB " + rci.getNumContactsLB());
                        System.out.println("#!# CL " + rci.getNumContactsCL());
                        System.out.println("#!# LC " + rci.getNumContactsLC());
                        System.out.println("#!# LL " + rci.getNumContactsLL());
                        
                        System.out.println("#!# NHPI " + rci.getNumContactsNHPI());
                        System.out.println("#!# PINH " + rci.getNumContactsPINH());
                        System.out.println("#!# CAHPI " + rci.getNumContactsCAHPI());
                        System.out.println("#!# PICAH " + rci.getNumContactsPICAH());
                        System.out.println("#!# CNHPI " + rci.getNumContactsCNHPI());
                        System.out.println("#!# PICNH " + rci.getNumContactsPICNH());
                        System.out.println("#!# SHPI " + rci.getNumContactsSHPI());
                        System.out.println("#!# PISH " + rci.getNumContactsPISH());
                        System.out.println("#!# XOHPI " + rci.getNumContactsXOHPI());
                        System.out.println("#!# PIXOH " + rci.getNumContactsPIXOH());
                        System.out.println("#!# PROCDHPI " + rci.getNumContactsPROCDHPI());
                        System.out.println("#!# PIPROCDH " + rci.getNumContactsPIPROCDH());
                        System.out.println("#!# CCACOH " + rci.getNumContactsCCACOH());
                        System.out.println("#!# CCOCAH " + rci.getNumContactsCCOCAH());
                        System.out.println("#!# BCACOH " + rci.getNumContactsBCACOH());
                        System.out.println("#!# BCOCAH " + rci.getNumContactsBCOCAH());
                        */
                    }
                }
                else {
                    numResContactsImpossible++;
                    //System.out.println("    DSSP res# " + a.getDsspNum() + "/" + b.getDsspNum() + ": No atom contact possible, skipping atom level checks.");
                }
            }
        }
        
        if(! FileParser.silent) {
            System.out.println("  Checked " + numResContactsChecked + " contacts for " + rs + " residues: " + numResContactsPossible + " possible, " + contactInfo.size() + " found, " + numResContactsImpossible + " impossible (collision spheres check).");
        }
        
        
        
        return(contactInfo);
    }
    
    
    private static long calculateSkipNeighborNum(Residue res1, Residue res2, int maxSequenceNeighborDist, int currentSeqPos, int SeqLength) {
        // jnw_2019: following taken from old contact computation and adopted such that maxSeqNeighDist without ligands and within each chain
        //   See there for comments how sequence neigbhor skip works in general (removed here for brevity)
        
        if (maxSequenceNeighborDist > 0) {
            int justToBeSure = 2;  // Account for rounding errors up to justToBeSure Angström (eg. 2 means 0.2)
            int combinedAtomRadius = 0;
            long spaceBetweenResidues, numResToSkip;

            combinedAtomRadius += (res1.isLigand()) ? Settings.getInteger("plcc_I_lig_atom_radius") : Settings.getInteger("plcc_I_atom_radius");
            combinedAtomRadius += (res2.isLigand()) ? Settings.getInteger("plcc_I_lig_atom_radius") : Settings.getInteger("plcc_I_atom_radius");

            spaceBetweenResidues = res1.distTo(res2) - (combinedAtomRadius + res1.getSphereRadius() + res2.getSphereRadius() + justToBeSure);

            //DEBUG
            /*
            System.out.println("ResCenterDist: " + res1.distTo(res2));
            System.out.println("Res1 CenterSphereRadius: " + res1.getSphereRadius());
            System.out.println("Res2 CenterSphereRadius: " + res2.getSphereRadius());
            System.out.println("combinedAtomRadius: " + combinedAtomRadius);
            System.out.println("spaceBetweenResidues: " + spaceBetweenResidues);
            System.out.println("chainMaxSeqNeighborAAResDist " + maxSequenceNeighborDist);
            */

            if(spaceBetweenResidues > maxSequenceNeighborDist) {           
                numResToSkip = spaceBetweenResidues / maxSequenceNeighborDist;

                if(Settings.getInteger("plcc_I_debug_level") >= 2) {
                    System.out.println("  [DEBUG LV 2] Residue skipping kicked in for DSSP res " + res1.getDsspNum() + ", skipped " + numResToSkip + " residues after " + res2.getDsspNum() + " in distance " + res1.distTo(res2));
                }

                // preserve correct statistics if skip
                // +1 b/c j is incremented after loop body
                if (currentSeqPos + numResToSkip + 1 > SeqLength) {
                    // only add which are really skipped
                    return SeqLength - currentSeqPos - 1;
                } else {
                    return numResToSkip;
                }
            } else {
                return 0;
            }
        } else {
            // Means that there are no residues in the chain, ergo no skip to be calculated -> REALLY?
            return 0;
        }
    }
    
    
     /**
     * Calculates all atom contacts between all chains which are in contact. Speed up for atom
     * contact computation, especially for structures of many chains.
     * @param chains list of chains
     * @return ArrayList of ResContactInfo holding all the contact information
     */
    public static ArrayList<MolContactInfo> calculateAllContactsChainSphereSpeedup(List<Chain> chains) {
        Boolean silent = Settings.getBoolean("plcc_B_silent");
        Chain chainA, chainB;
        Integer chainCount = chains.size();
        int numberResTotal = 0;
        MolContactInfo rci;
        ArrayList<MolContactInfo> contactInfo = new ArrayList<>();
        Residue res1, res2;
        //Molecule mol1, mol2;

        // variables for statistics
        long numResContactsChecked, numResContactsPossible, numResContactsImpossible, chainSkippedRes, seqNeighSkippedResIntraChain, seqNeighSkippedResInterChain;
        int chainChainSkipped = 0;
        numResContactsChecked = numResContactsPossible = numResContactsImpossible = chainSkippedRes = seqNeighSkippedResIntraChain = seqNeighSkippedResInterChain = 0;
        Integer numIgnoredLigandContacts = 0;
        long numResToSkip;  // also for skipping
        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            chainCount = 2;
            System.out.println("DEBUG: Warning: Limiting residue contact computation to the first " + chainCount + " chains and their residues.");            
        } 
        
        // loop over chains and residues
        for (int k = 0; k < chainCount; k++) {
            chainA = chains.get(k);
            
            // skip chain if no (protein) atoms in it
            if (chainA.getRadiusFromCentroid() == -1) {
                continue;
            }
            
            int chainANumberResidues = chainA.getResidues().size();
            ArrayList<Residue> AAResiduesA = new ArrayList<>();
            ArrayList<Residue> ligResiduesA = new ArrayList<>();
            AAResiduesA.addAll(chainA.getAllAAResidues());
            ligResiduesA.addAll(chainA.getAllLigandResidues());
            
            numberResTotal += chainANumberResidues;
            
            Integer chainAMaxSeqNeighborAADist = chainA.getMaxSeqNeighborAADist();
            
            // - - - contacts within chains (incl sequence neighbor skip)  - - -
            //
            // multiple loops: 
            // 1) res1 = AA
            //   1.1) res2 = AA -> seq neigh skip
            //   HINT: 1.2) res2 = lig NOT needed (same as Lig-AA, but no skip possible)
            // 2) res1 = lig
            //   2.1) res2 = AA -> seq neigh skip
            //   2.2) res2 = lig -> no skip possible b/c we do not calculate MaxSeqNeighborDist for ligands
            
            // 1)
            for (int i = 0; i < AAResiduesA.size(); i++) {
                res1 = AAResiduesA.get(i);
                
                // 1.1)
                for (int j = i + 1; j < AAResiduesA.size(); j++) {
                    res2 = AAResiduesA.get(j);

                    if(Settings.getInteger("plcc_I_debug_level") >= 1) {
                        System.out.println("  [DEBUG LV 1] Checking DSSP pair (loop 1.1) " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                    }
                    
                    numResContactsChecked++;

                    // We only need to check on atom level if the center spheres overlap
                    if (res1.contactPossibleWithResidue(res2)) {                                        
                        numResContactsPossible++;

                        rci = calculateAtomContactsBetweenResidues(res1, res2);
                        if( rci != null) {
                            // There were atoms contacts!
                            // there cannot be a lig in this contact -> always add without checking for plcc_B_write_lig_geolig
                            contactInfo.add(rci);
                        }
                    }
                    else {
                        numResContactsImpossible++;
                        numResToSkip = calculateSkipNeighborNum(res1, res2, chainAMaxSeqNeighborAADist, j, AAResiduesA.size());
                        j += numResToSkip;
                        seqNeighSkippedResIntraChain += numResToSkip;
                    }
                }
            }
            
            // 2)
            // can be skipped if plcc_B_write_lig_geolig = false
            if (Settings.getBoolean("plcc_B_write_lig_geolig")) {
                for (int i = 0; i < ligResiduesA.size(); i++) {
                    res1 = ligResiduesA.get(i);
                                       
                    // 2.1)
                    // here we have to start at j = 0 b/c it is another list!
                    for (int j = 0; j < AAResiduesA.size(); j++) {

                        res2 = AAResiduesA.get(j);

                        if (Settings.getInteger("plcc_I_debug_level") >= 1) {
                            System.out.println("  [DEBUG LV 1] Checking DSSP pair (loop 2.1) " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                        }

                        numResContactsChecked++;

                        // We only need to check on atom level if the center spheres overlap
                        if (res1.contactPossibleWithResidue(res2)) {                                        
                            numResContactsPossible++;

                            rci = calculateAtomContactsBetweenResidues(res1, res2);
                            if( rci != null) {
                                // There were atoms contacts!
                                contactInfo.add(rci);
                            }
                        }
                        else {
                            numResContactsImpossible++;
                            numResToSkip = calculateSkipNeighborNum(res1, res2, chainAMaxSeqNeighborAADist, j, AAResiduesA.size());
                            j += numResToSkip;
                            seqNeighSkippedResIntraChain += numResToSkip;
                        }
                    }

                    // 2.2)
                    for(int j = i + 1; j < ligResiduesA.size(); j++) {

                        res2 = ligResiduesA.get(j);

                        if(Settings.getInteger("plcc_I_debug_level") >= 1) {
                            System.out.println("  [DEBUG LV 1] Checking DSSP (loop 2.2) pair " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                        }

                        numResContactsChecked++;

                        // We only need to check on atom level if the center spheres overlap
                        if(res1.contactPossibleWithResidue(res2)) {                                        
                            numResContactsPossible++;

                            rci = calculateAtomContactsBetweenResidues(res1, res2);
                            if( rci != null) {
                                // There were atoms contacts!
                                contactInfo.add(rci);
                            }
                        }
                        else {
                            numResContactsImpossible++;

                            // lig-lig no seq neigh skip possible
                        }
                    }
                }
            }
            
            // - - - contacts between chains (only if chains overlap!)(incl seq neighbor skip) - - -
            //
            // 1 & 2 may not correspond to chainA & B, instead reorder them depending on where we expect most skips to happen:
            //   -> longer chain as 2 (b/c the longer the chain the more possible skips)
            //   (if (nearly) equally long, chain as 2 which has smaller maxSeqNeighborDist (b/c this allows more skips in theory))
            //
            // multiple loops: 
            // 1) res1 = AA
            //   1.1) res2 = AA -> seq neigh skip
            //   HINT: 1.2) res2 = lig NOT needed (same as Lig-AA, but no skip possible)
            // 2) res1 = lig
            //   2.1) res2 = AA -> seq neigh skip
            //   2.2) res2 = lig -> no skip possible b/c we do not calculate MaxSeqNeighborDist for ligands
            // 3) res1 = lig (from other chain)
            //   3.1) res2 = AA (from THIS chain) -> seq neigh skip
            
            for (int l = k + 1; l < chainCount; l++) {
                
                chainB = chains.get(l);
                int chainBNumberResidues = chainB.getResidues().size();
                
                // skip chain if no (protein) atoms in it
                if (chainB.getRadiusFromCentroid() == -1) {
                    continue;
                }
                
                // vars for possible swapping of inner and outer loop to maximize skips -> only for loop 1.1) b/c we have no maxSeqNeighborDist for ligands!
                ArrayList<Residue> innerLoopChainAAs = new ArrayList<>();
                ArrayList<Residue> outerLoopChainAAs = new ArrayList<>();
                int innerChainMaxSeqNeighborAADist;
                innerChainMaxSeqNeighborAADist = 0;  // value needs to be initialized for Netbeans, just take something small

                // check chain overlap
                if (chainA.contactPossibleWithChain(chainB)) {
                    ArrayList<Residue> AAResiduesB = new ArrayList<>();
                    ArrayList<Residue> ligResiduesB = new ArrayList<>();
                    AAResiduesB.addAll(chainB.getAllAAResidues());
                    ligResiduesB.addAll(chainB.getAllLigandResidues());

                    Integer chainBMaxSeqNeighborAADist = chainB.getMaxSeqNeighborAADist();
                                       
                    // decide which chain as 2 (used for skip in loop 1.1), atm just take longer chain as 2
                    if (chainA.getAllAAResidues().size() > chainB.getAllAAResidues().size()) {
                        outerLoopChainAAs.addAll(AAResiduesB);
                        innerLoopChainAAs.addAll(AAResiduesA);
                        innerChainMaxSeqNeighborAADist = chainBMaxSeqNeighborAADist;
                    } else {
                        outerLoopChainAAs.addAll(AAResiduesA);
                        innerLoopChainAAs.addAll(AAResiduesB);
                        innerChainMaxSeqNeighborAADist = chainAMaxSeqNeighborAADist;
                    }
                                       
                    // 1)
                    for (int i = 0; i < outerLoopChainAAs.size(); i++) {

                        res1 = outerLoopChainAAs.get(i);
                        
                        // 1.1)
                        // NOTE: we cant just go from j = i + 1 on now or we would miss some contacts!
                        for (int j = 0; j < innerLoopChainAAs.size(); j++) {
                            res2 = innerLoopChainAAs.get(j);

                            if (Settings.getInteger("plcc_I_debug_level") >= 1) {
                                if(! silent) {
                                    System.out.println("  Checking DSSP pair " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                                }
                            }                

                            numResContactsChecked++;

                            // We only need to check on atom level if the center spheres overlap
                            if (res1.contactPossibleWithResidue(res2)) {                                        
                                numResContactsPossible++;

                                rci = calculateAtomContactsBetweenResidues(res1, res2);
                                if( rci != null) {
                                    // There were atoms contacts!
                                    // there cannot be a lig in this contact -> always add without checking for plcc_B_write_lig_geolig
                                    contactInfo.add(rci);
                                }
                            }
                            else {
                                numResContactsImpossible++;
                                numResToSkip = calculateSkipNeighborNum(res1, res2, innerChainMaxSeqNeighborAADist, j, innerLoopChainAAs.size());
                                j += numResToSkip;
                                seqNeighSkippedResInterChain += numResToSkip;                            
                            }
                        }
                    }
                    
                    // can be skipped if plcc_B_write_lig_geolig = false
                    if (Settings.getBoolean("plcc_B_write_lig_geolig")) {
                        
                        // 2)
                        for (int i = 0; i < ligResiduesA.size(); i++) {

                            res1 = ligResiduesA.get(i);
                            
                            // 2.1)
                            // NOTE: we cant just go from j = i + 1 on now or we would miss some contacts!
                            for (int j = 0; j < AAResiduesB.size(); j++) {

                                res2 = AAResiduesB.get(j);

                                if (Settings.getInteger("plcc_I_debug_level") >= 1) {
                                    if(! silent) {
                                        System.out.println("  Checking DSSP pair " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                                    }
                                }                

                                numResContactsChecked++;

                                // We only need to check on atom level if the center spheres overlap
                                if (res1.contactPossibleWithResidue(res2)) {                                        
                                    numResContactsPossible++;

                                    rci = calculateAtomContactsBetweenResidues(res1, res2);
                                    if( rci != null) {
                                        // There were atoms contacts!
                                        // there cannot be a lig in this contact -> always add without checking for plcc_B_write_lig_geolig
                                        contactInfo.add(rci);
                                    }
                                }
                                else {
                                    numResContactsImpossible++;
                                    numResToSkip = calculateSkipNeighborNum(res1, res2, chainBMaxSeqNeighborAADist, j, AAResiduesB.size());
                                    j += numResToSkip;
                                    seqNeighSkippedResInterChain += numResToSkip;                            
                                }
                            }
                            
                            // 2.2)
                            // NOTE: we cant just go from j = i + 1 on now or we would miss some contacts!
                            for (int j = 0; j < ligResiduesB.size(); j++) {

                                res2 = ligResiduesB.get(j);

                                if (Settings.getInteger("plcc_I_debug_level") >= 1) {
                                    if(! silent) {
                                        System.out.println("  Checking DSSP pair " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                                    }
                                }                

                                numResContactsChecked++;

                                // We only need to check on atom level if the center spheres overlap
                                if (res1.contactPossibleWithResidue(res2)) {                                        
                                    numResContactsPossible++;

                                    rci = calculateAtomContactsBetweenResidues(res1, res2);
                                    if( rci != null) {
                                        // There were atoms contacts!
                                        // there cannot be a lig in this contact -> always add without checking for plcc_B_write_lig_geolig
                                        contactInfo.add(rci);
                                    }
                                }
                                else {
                                    numResContactsImpossible++;
                                    
                                    // no skip possible for lig-lig
                                }
                            }
                        }
                        
                        // 3)
                        // HINT: for this nested loop outer and inner are swapped: this way we can profit from the lig-AA skipping
                        for (int i = 0; i < ligResiduesB.size(); i++) {

                            res1 = ligResiduesB.get(i);

                            // 3.1)
                            // NOTE: we cant just go from j = i + 1 on now or we would miss some contacts!
                            for (int j = 0; j < AAResiduesA.size(); j++) {

                                res2 = AAResiduesA.get(j);

                                if (Settings.getInteger("plcc_I_debug_level") >= 1) {
                                    if(! silent) {
                                        System.out.println("  Checking DSSP pair " + res1.getDsspNum() + "/" + res2.getDsspNum() + "...");
                                    }
                                }                

                                numResContactsChecked++;

                                // We only need to check on atom level if the center spheres overlap
                                if (res1.contactPossibleWithResidue(res2)) {                                        
                                    numResContactsPossible++;

                                    rci = calculateAtomContactsBetweenResidues(res1, res2);
                                    if( rci != null) {
                                        // There were atoms contacts!
                                        // there cannot be a lig in this contact -> always add without checking for plcc_B_write_lig_geolig
                                        contactInfo.add(rci);
                                    }
                                }
                                else {
                                    numResContactsImpossible++;
                                    numResToSkip = calculateSkipNeighborNum(res1, res2, chainAMaxSeqNeighborAADist, j, AAResiduesA.size());
                                    j += numResToSkip;
                                    seqNeighSkippedResInterChain += numResToSkip;                            
                                }
                            }
                        }
                    }
                } else {
                    chainChainSkipped++;
                    chainSkippedRes += (chainANumberResidues * chainBNumberResidues);
                }
            }
        }
        
        // - - - statistics - - -
        //
        int maxChainChainContactsPossible;
        if (chains.size() > 1) {
            maxChainChainContactsPossible = (chains.size() * (chains.size() - 1)) / 2;
        } else {
            maxChainChainContactsPossible = 0;
        }

        if(! FileParser.silent) {
            System.out.println("  Skipped " + chainChainSkipped + " chain-chain contacts (and " + chainSkippedRes + " otherwise checked residue contacts) of " + maxChainChainContactsPossible + " maximal contacts due to chain sphere check.");
            System.out.println("  Skipped " + seqNeighSkippedResIntraChain + " intra chain and " + seqNeighSkippedResInterChain + " inter chain residue contacts due to sequence neighbor skip.");
            System.out.println("  Checked " + numResContactsChecked + " contacts for " + numberResTotal + " residues: " + numResContactsPossible + " possible, " + contactInfo.size() + " found, " + numResContactsImpossible + " impossible (collison spheres check).");
        }

        if( ! Settings.getBoolean("plcc_B_write_lig_geolig")) {
            if(! FileParser.silent) {
                if (Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                    System.out.println("  Configured to ignore ligands. Because of chain sphere speedup we do not even know how much were ignored.");
                } else {
                    System.out.println("  Configured to ignore ligands, ignored " + numIgnoredLigandContacts + " ligand contacts.");
                }
            }
        }
        
        return contactInfo;
    }
    
    
    /**
     * Calculates all contacts between the residues in mols.
     * @param mols A list of Residue objects.
     * @return A list of MolContactInfo objects, each representing a pair of residues that are in contact.
     */
    public static ArrayList<MolContactInfo> calculateAllContacts(ArrayList<Molecule> mols) {
        
        
        Boolean silent = Settings.getBoolean("plcc_B_silent");
        
        Molecule a, b;
        Integer rs = mols.size();
        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            rs = 2;
            System.out.println("DEBUG: Warning: Limiting residue contact computation to the first " + rs + " residues.");            
        }        

        long numResContactsChecked, numResContactsPossible, numResContactsImpossible, numCmpSkipped;
        numResContactsChecked = numResContactsPossible = numResContactsImpossible = numCmpSkipped = 0;

        long numResToSkip, spaceBetweenResidues;
        MolContactInfo rci;
        ArrayList<MolContactInfo> contactInfo = new ArrayList<MolContactInfo>();

        Integer atomRadius = Settings.getInteger("plcc_I_atom_radius");
        Integer atomRadiusLig = Settings.getInteger("plcc_I_lig_atom_radius");

        
        //System.out.println("  Atom radius set to " + atomRadius + " for protein atoms, " + atomRadiusLig + " for ligand atoms (unit is 1/10th Angstroem).");

        Integer globalMaxCollisionRadius = globalMaxCenterSphereRadius + atomRadius;
        Integer globalMaxCenterSphereDiameter = globalMaxCollisionRadius * 2;
        Integer numIgnoredLigandContacts = 0;

        for(int i = 0; i < rs; i++) {

            a = mols.get(i);
            numResToSkip = 0L;

            
            
            for(int j = i + 1; j < rs; j++) {

                b = mols.get(j);
                
                // DEBUG
                if(Settings.getInteger("plcc_I_debug_level") >= 1) {
                    if(! silent) {
                        System.out.println("  Checking DSSP pair " + a.getDsspNum() + " (Chain " + 
                                a.getChainID() + " Residue " + a.getPdbNum() + ") and " + 
                                b.getDsspNum() + " (Chain " +  b.getChainID() + " Residue " +
                                b.getPdbNum() + ") ...");
                    }
                    //System.out.println("    " + a.getAtomsString());
                    //System.out.println(a.atomInfo());
                    //System.out.println("    " + b.getAtomsString());
                    //System.out.println(b.atomInfo());
                }                

                numResContactsChecked++;

                // We only need to check on atom level if the center spheres overlap
                if(a.contactPossibleWithResidue(b)) {                                        
                    numResContactsPossible++;

                    //System.out.println("    DSSP mols# " + a.getDsspNum() + "/" + b.getDsspNum() + ": Collision spheres overlap, checking on atom level.");

                    rci = calculateAtomContactsBetweenResidues(a, b);
                    if( rci != null) {
                        // There were atoms contacts!

                        if(Settings.getBoolean("plcc_B_write_lig_geolig")) {
                            // Just add it without asking questions about the residue types
                            contactInfo.add(rci);
                        }
                        else {
                            // We should ignore ligand contacts
                            if(a.getType().equals(1) || b.getType().equals(1)) {
                                // This IS a ligand contact so ignore it
                                numIgnoredLigandContacts++;
                                // System.out.println("  Ignored ligand contact between DSSP residues " + a.getDsspNum() + " and " + b.getDsspNum() + ".");
                            }
                            else {
                                // This is NOT a ligand contact so add it
                                contactInfo.add(rci);
                            }
                        }
                    }
                }
                else {
                    numResContactsImpossible++;
                    //System.out.println("    DSSP mols# " + a.getDsspNum() + "/" + b.getDsspNum() + ": No atom contact possible, skipping atom level checks.");
                       
                    // Further speedup: If the distance of a residue to another Residue is very large we
                    //  may be able to skip some of the next residues (I. Koch):
                    //  If the distance between them is

                    spaceBetweenResidues = a.distTo(b) - (2 * atomRadius + a.getSphereRadius() + b.getSphereRadius());                  
                    if(spaceBetweenResidues > globalMaxSeqNeighborResDist) {
                        // In this case we can skip at least one residue.

                        // How often does the max dist between to sequential neighbor residues fit into the space between these two residues 'a' and 'b' ? If it fits in there n times we can
                        //  skip the next n residues: even if they were arranged in a straight line from 'a' to 'b' they could not reach it!
                        // numResToSkip = spaceBetweenResidues / globalMaxCenterSphereDiameter;
                        numResToSkip = spaceBetweenResidues / globalMaxSeqNeighborResDist;

                        if(Settings.getInteger("plcc_I_debug_level") >= 2) {
                            System.out.println("  [DEBUG LV 2] Residue skipping kicked in for DSSP res " + a.getDsspNum() + ", skipped " + numResToSkip + " residues after " + b.getDsspNum() + " in distance " + a.distTo(b));
                        }
                        
                        j += numResToSkip;
                        numCmpSkipped += numResToSkip;

                    }
                }

            }

        }


        if(! FileParser.silent) {
            System.out.println("  Checked " + numResContactsChecked + " contacts for " + rs + " residues: " + numResContactsPossible + " possible, " + contactInfo.size() + " found, " + numResContactsImpossible + " impossible (collison spheres check).");
            System.out.println("  Did not check " + numCmpSkipped + " contacts (skipped by seq neighbors check), would have been " + (numResContactsChecked + numCmpSkipped)  + ".");
        }

        if( ! Settings.getBoolean("plcc_B_write_lig_geolig")) {
            if(! FileParser.silent) {
                System.out.println("  Configured to ignore ligands, ignored " + numIgnoredLigandContacts + " ligand contacts.");
            }
        }

        return(contactInfo);

    }
    
    /**
     * Calculates all contacts between the residues in res.
     * @param res A list of Residue objects.
     * @return A list of MolContactInfo objects, each representing a pair of residues that are in contact.
     */
    public static ArrayList<MolContactInfo> calculateAllContactsLimitedByChain(List<Residue> res, String handledChain) {
        
        Boolean silent = Settings.getBoolean("plcc_B_silent");
                
        Residue a, b;
        Integer rs = res.size();
        String chainTag = "Chain " + handledChain + ": ";
        
        if(Settings.getBoolean("plcc_B_contact_debug_dysfunct")) {
            rs = 2;
            System.out.println("DEBUG: " + chainTag + "Limiting residue contact computation to the first " + rs + " residues.");            
        }        

        Integer numResContactsChecked, numResContactsPossible, numResContactsImpossible, numCmpSkipped;
        numResContactsChecked = numResContactsPossible = numResContactsImpossible = numCmpSkipped = 0;

        Integer numResToSkip, spaceBetweenResidues;
        MolContactInfo rci;
        ArrayList<MolContactInfo> contactInfo = new ArrayList<MolContactInfo>();

        Integer atomRadius = Settings.getInteger("plcc_I_atom_radius");
        Integer atomRadiusLig = Settings.getInteger("plcc_I_lig_atom_radius");

        if( ! silent) {
            System.out.println("  " + chainTag + "Atom radius set to " + atomRadius + " for protein atoms, " + atomRadiusLig + " for ligand atoms (unit is 1/10th Angstroem).");
        }

        Integer globalMaxCollisionRadius = globalMaxCenterSphereRadius + atomRadius;
        Integer globalMaxCenterSphereDiameter = globalMaxCollisionRadius * 2;
        Integer numIgnoredLigandContacts = 0;
        Integer numResPairsSkippedWrongChain = 0;
        
        Boolean includeLigandsFromOtherChains = Settings.getBoolean("plcc_B_consider_all_ligands_for_each_chain");
        /*
        if(includeLigandsFromOtherChains) {
            if(!silent) {
                DP.getInstance().i("Main", "calculateAllContactsLimitedByChain: Also considering valid 3D ligand contacts from ligands associated to other chains in PDB file.");
            }
        }
        */
        
        int numLigandsFromOtherChainsKept = 0;
        
        
        for(Integer i = 0; i < rs; i++) {

            a = res.get(i);
            numResToSkip = 0;
            
            if( ! a.getChainID().equals(handledChain)) {
                if( ! (includeLigandsFromOtherChains && a.isLigand())) {
                    numResPairsSkippedWrongChain += (rs - i);   // skipped all this_atom -- other pairs
                    continue;
                } else {
                    numLigandsFromOtherChainsKept++;
                }
            }

            
            
            for(Integer j = i + 1; j < rs; j++) {

                b = res.get(j);
                if( ! b.getChainID().equals(handledChain)) {
                    if( ! (includeLigandsFromOtherChains && b.isLigand())) {
                        numResPairsSkippedWrongChain++;
                        continue;
                    } 
                }
                
                // DEBUG
                if(Settings.getInteger("plcc_I_debug_level") >= 1) {
                    System.out.println("  " + chainTag + "Checking DSSP pair " + a.getDsspNum() + "/" + b.getDsspNum() + "...");
                    //System.out.println("    " + a.getAtomsString());
                    //System.out.println(a.atomInfo());
                    //System.out.println("    " + b.getAtomsString());
                    //System.out.println(b.atomInfo());
                }                

                numResContactsChecked++;

                // We only need to check on atom level if the center spheres overlap
                if(a.contactPossibleWithResidue(b)) {                                        
                    numResContactsPossible++;

                    //System.out.println("    DSSP res# " + a.getDsspNum() + "/" + b.getDsspNum() + ": Collision spheres overlap, checking on atom level.");

                    rci = calculateAtomContactsBetweenResidues(a, b);
                    if( rci != null) {
                        // There were atoms contacts!

                        if(Settings.getBoolean("plcc_B_write_lig_geolig")) {
                            // Just add it without asking questions about the residue types
                            contactInfo.add(rci);
                        }
                        else {
                            // We should ignore ligand contacts
                            if(a.getType().equals(Residue.RESIDUE_TYPE_LIGAND) || b.getType().equals(Residue.RESIDUE_TYPE_LIGAND)) {
                                // This IS a ligand contact so ignore it
                                numIgnoredLigandContacts++;
                                // System.out.println("  Ignored ligand contact between DSSP residues " + a.getDsspNum() + " and " + b.getDsspNum() + ".");
                            }
                            else {
                                // This is NOT a ligand contact so add it
                                contactInfo.add(rci);
                            }
                        }
                    }
                }
                else {
                    numResContactsImpossible++;
                    //System.out.println("    DSSP res# " + a.getDsspNum() + "/" + b.getDsspNum() + ": No atom contact possible, skipping atom level checks.");

                    // Further speedup: If the distance of a residue to another Residue is very large we
                    //  may be able to skip some of the next residues (I. Koch):
                    //  If the distance between them is

                    spaceBetweenResidues = a.distTo(b) - (2 * atomRadius + a.getSphereRadius() + b.getSphereRadius());
                    if(spaceBetweenResidues > globalMaxSeqNeighborResDist) {
                        // In this case we can skip at least one residue.

                        // How often does the max dist between to sequential neighbor residues fit into the space between these two residues 'a' and 'b' ? If it fits in there n times we can
                        //  skip the next n residues: even if they were arranged in a straight line from 'a' to 'b' they could not reach it!
                        // numResToSkip = spaceBetweenResidues / globalMaxCenterSphereDiameter;
                        numResToSkip = spaceBetweenResidues / globalMaxSeqNeighborResDist;

                        // System.out.println("  Residue skipping kicked in for DSSP res " + a.getDsspNum() + ", skipped " + numResToSkip + " residues after " + b.getDsspNum() + " in distance " + a.resCenterDistTo(b) + ".");
                        j += numResToSkip;
                        numCmpSkipped += numResToSkip;

                    }
                }

            }

        }


        if( ! silent) {
            System.out.println("  " + chainTag + "Checked " + numResContactsChecked + " contacts for " + rs + " residues: " + numResContactsPossible + " possible, " + contactInfo.size() + " found, " + numResContactsImpossible + " impossible (collison spheres check).");
            System.out.println("  " + chainTag + "Did not check " + numResPairsSkippedWrongChain + " residue pairs because they were part of different chains.");
            System.out.println("  " + chainTag + "Did not check " + numCmpSkipped + " contacts (skipped by seq neighbors check), would have been " + (numResContactsChecked + numCmpSkipped)  + ".");
            if(includeLigandsFromOtherChains) {
                System.out.println("  " + chainTag + "Checked residues of this chain for contacts with " + numLigandsFromOtherChainsKept + " ligands from other chains.");
            } else {
                System.out.println("  " + chainTag + "Ignored ligands assigned to other chains in the PDB file when computing contacts of residues of this chain.");
            }
            
            if( ! Settings.getBoolean("plcc_B_write_lig_geolig")) {
                System.out.println("  " + chainTag + "Configured to ignore ligands, ignored " + numIgnoredLigandContacts + " ligand contacts.");
            }
        }
        
        

        return(contactInfo);

    }

    
    /**
     * Calculates the atom contacts between the residues 'a' and 'b'.
     * @param a one of the residues of the residue pair
     * @param b one of the residues of the residue pair
     * @return A MolContactInfo object with information on the atom contacts between 'a' and 'b'.
     */
    public static MolContactInfo calculateAtomContactsBetweenResidues(Molecule a, Molecule b) {

        
        ArrayList<Atom> atoms_a = a.getAtoms();
        ArrayList<Atom> atoms_b = b.getAtoms();

        Atom x, y;
        Integer dist = null;
        Integer CAdist = a.distTo(b);
        MolContactInfo result = null;


        Integer[] numPairContacts = new Integer[Main.NUM_RESIDUE_PAIR_CONTACT_TYPES];
        // The positions in the numPairContacts array hold the number of contacts of each type for a pair of residues:
        // Some cheap vars to make things easier to understand (a replacement for #define):
        /*
        Integer TT = 0;         //  0 = total number of contacts            (all residue type combinations)
        Integer BB = 1;         //  1 = # of backbone-backbone contacts     (protein - protein only)
        Integer CB = 2;         //  2 = # of sidechain-backbone contacts    (protein - protein only)
        Integer BC = 3;         //  3 = # of backbone-sidechain contacts    (protein - protein only)
        Integer CC = 4;         //  4 = # of sidechain-sidechain contacts   (protein - protein only)
        Integer HB = 5;         //  5 = # of H-bridge contacts 1, N=>0      (protein - protein only)
        Integer BH = 6;         //  6 = # of H-bridge contacts 2, 0=>N      (protein - protein only)
        Integer BL = 7;         //  7 = # of backbone-ligand contacts       (protein - ligand only)
        Integer LB = 8;         //  8 = # of ligand-backbone contacts       (protein - ligand only)
        Integer CL = 9;         //  9 = # of sidechain-ligand contacts      (protein - ligand only)
        Integer LC = 10;        // 10 = # of ligand-sidechain contacts      (protein - ligand only)
        Integer LL = 11;        // 11 = # of ligand-ligand contacts         (ligand - ligand only)
        */


        Integer numTotalLigContactsPair = 0;



        Integer[] minContactDistances = new Integer[numPairContacts.length];
        // Holds the minimal distances of contacts of the appropriate type (see numPairContacts, index 0 is unused)

        Integer[] contactAtomNumInResidueA = new Integer[numPairContacts.length];
        // Holds the number Atom x has in its residue a for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)

        Integer[] contactAtomNumInResidueB = new Integer[numPairContacts.length];
        // Holds the number Atom y has in its residue b for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)
        

        for(Integer k = 0; k < numPairContacts.length; k++) {      // init all arrays
            numPairContacts[k] = 0;
            minContactDistances[k] = -1;    // We are looking for the smallest distance >= 0 in this function so do NOT set '0' as the initial value or everything will get fucked up!
            contactAtomNumInResidueA[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            contactAtomNumInResidueB[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            // The initial value of '-1' is required for the atom index arrays because our index
            // starts at '0', but geom_neo treads a '0' in a line of <pdbid>.geo as 'no contact' because
            // it starts its index at '1'.
            // This problem is solved by the functions in MolContactInfo: they return (our_index + 1). This
            // means that:
            //   1) If no contact was detected, our_index is -1 and they return 0, which means 'no contact' to geom_neo.
            //   2) If a contact was detected, our_index is converted to the geom_neo index. :)
        }

        // We assume that the first 5 atoms (index 0..4) in a residue that is an AA are backbone atoms,
        //  while all other (6..END) are assumed to be side chainName atoms.
        //  The backbone atoms should have atom names ' N  ', ' CA ', ' C  ' and ' O  ' but we don't check
        //  this atm because geom_neo doesn't do that and we want to stay compatible.
        //  Of course, all of this only makes sense for resides that are AAs, not for ligands. We care for that.
        //  jnw_2019: changed to 3 b/c there are only 4 bb atoms from 0..3: N, CA, C, O
        Integer numOfLastBackboneAtomInResidue = 3;
        Integer atomIndexOfBackboneN = 0;       // backbone nitrogen atom index
        Integer atomIndexOfBackboneO = 3;       // backbone oxygen atom index

        Integer aIntID = a.getInternalAAID();     // Internal AA ID (ALA=1, ARG=2, ...)
        Integer bIntID = b.getInternalAAID();
        Integer statAtomIDi, statAtomIDj;



        // Iterate through all atoms of the two residues and check contacts for all pairs
        outerloop:
        for(Integer i = 0; i < atoms_a.size(); i++) {
            
            
            if(i >= MAX_ATOMS_PER_AA && a.isAA()) {
                DP.getInstance().w("calculateAtomContactsBetweenResidues(): The AA residue " + a.getUniquePDBName() + " of type " + a.getName3() + " has more atoms than allowed, skipping atom #" + i + ".");
                break;
            }
            
            x = atoms_a.get(i);

            innerloop:
            for(Integer j = 0; j < atoms_b.size(); j++) {
                                
                if(j >= MAX_ATOMS_PER_AA && b.isAA()) {
                    DP.getInstance().w("calculateAtomContactsBetweenResidues(): The AA residue " + b.getUniquePDBName() + " of type " + b.getName3() + " has more atoms than allowed, skipping atom #" + j + ".");
                    //continue;
                    break outerloop;
                }
                
                y = atoms_b.get(j);
                                

                // Check whether a contact exist. If so, classify it. Note that the code of geom_neo works based on the
                //  position of an atom in the atom list of its residue (e.g., it assumes that the 2nd atom of an AA is
                //  the C alpha atom. While this seems to hold for many PDB files it will produce wrong results if atoms
                //  are missing from the PDB file or other weird stuff is going on in there.

                //System.out.println("      Checking atom pair " + x.getChemSym() + " and " + y.getChemSym() + " with residue index " + i + " and " + j + " of residues " + a.getFancyName() + " and " + b.getFancyName() + ".");
                //System.out.println("        " + x);
                //System.out.println("        " + y);

                dist = x.distToAtom(y);

                if(x.atomContactTo(y)) {             // If a contact is detected, Atom.atomContactTo() returns true
                    

                    // The van der Waals radii spheres overlap, contact found.
                    numPairContacts[MolContactInfo.TT]++;   // update total number of contacts for this residue pair
                    
                    // DEBUG
                    //System.out.println("DEBUG: Atom contact in distance " + dist + " between atom " + x + " and " + y + ".");


                    // Update contact statistics.
                    statAtomIDi = i + 1;    // The field '0' is used for all contacs and we need to follow geom_neo conventions so we start the index at 1 instead of 0.
                    statAtomIDj = j + 1;
                    if(x.isLigandAtom()) { statAtomIDi = 1; }       // Different ligands can have different numbers of atoms and separating them just makes no sense. We assign all contacts to the first atom.
                    if(y.isLigandAtom()) { statAtomIDj = 1; }
                    
                    try {

                        contact[0][0][0][0]++;                 // update global total number of contacts
                        contact[aIntID][bIntID][0][0]++;       // contacts AA type a <-> AA type b
                        contact[bIntID][aIntID][0][0]++;       // contacts AA type b <-> AA type a


                        //System.out.println("DEBUG: a=" + aIntID + ",b=" + bIntID + ",i=" + statAtomIDi + ",j=" + statAtomIDj + ". Residues=" + a.getFancyName() + "," + b.getFancyName() + ".");
                        contact[aIntID][bIntID][statAtomIDi][statAtomIDj]++;       // contacts of atoms of AAs
                        contact[bIntID][aIntID][statAtomIDj][statAtomIDi]++;

                        contact[aIntID][bIntID][statAtomIDi][0]++;                 // total number of contacts for atom x of this AA
                        contact[bIntID][aIntID][0][statAtomIDi]++;

                        contact[aIntID][bIntID][statAtomIDj][0]++;                 // total number of contacts for atom x of this AA
                        contact[bIntID][aIntID][0][statAtomIDj]++;
                    } catch(java.lang.ArrayIndexOutOfBoundsException e) {
                        //DP.getInstance().w("calculateAtomContactsBetweenResidues():Contact statistics array out of bounds. Residues with excessive number of atoms detected: " + e.getMessage() + ".");
                        DP.getInstance().w("calculateAtomContactsBetweenResidues(): Atom count for residues too high (" + e.getMessage() + "), ignoring contacts for these atoms (aIntID=" + aIntID + ", bIntID=" + bIntID + ", statAtomIDi=" + statAtomIDi + ", statAtomIDj=" + statAtomIDj + ").");
                        continue;
                    }
                    
                    // Determine the contact type.                    
                    if(x.isProteinAtom() && y.isProteinAtom()) {
                        // *************************** protein - protein contact *************************


                        // Check the exact contact type
                        if(i <= numOfLastBackboneAtomInResidue && j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - backbone contact
                            numPairContacts[MolContactInfo.BB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BB] < 0) || dist < minContactDistances[MolContactInfo.BB]) {
                                minContactDistances[MolContactInfo.BB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BB] = i;
                                contactAtomNumInResidueB[MolContactInfo.BB] = j;
                            }

                        }
                        else if(i > numOfLastBackboneAtomInResidue && j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a chainName - backbone contact
                            numPairContacts[MolContactInfo.CB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CB] < 0) || dist < minContactDistances[MolContactInfo.CB]) {
                                minContactDistances[MolContactInfo.CB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CB] = i;
                                contactAtomNumInResidueB[MolContactInfo.CB] = j;
                            }

                        }
                        else if(i <= numOfLastBackboneAtomInResidue && j > numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - chainName contact
                            numPairContacts[MolContactInfo.BC]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BC] < 0) || dist < minContactDistances[MolContactInfo.BC]) {
                                minContactDistances[MolContactInfo.BC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BC] = i;
                                contactAtomNumInResidueB[MolContactInfo.BC] = j;
                            }
                        }
                        else if(i > numOfLastBackboneAtomInResidue && j > numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a chainName - chainName contact
                            numPairContacts[MolContactInfo.CC]++;          // 'C' instead of 'S' for side chainName pays off

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CC] < 0) || dist < minContactDistances[MolContactInfo.CC]) {
                                minContactDistances[MolContactInfo.CC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CC] = i;
                                contactAtomNumInResidueB[MolContactInfo.CC] = j;
                            }
                        }
                        else {
                            System.err.println("ERROR: Congrats, you found a bug in the atom contact type determination code (res " + a.getPdbNum() + " atom " + i + " / res " + b.getPdbNum() + " atom " + j + ").");
                            System.err.println("ERROR: Atom types are: i (PDB atom #" + x.getPdbAtomNum() + ") => " + x.getAtomType() + ", j (PDB atom #" + y.getPdbAtomNum() + ") => " + y.getAtomType() + ".");
                            Main.doExit(1);
                        }

                        // Check for H bridges separately
                        if(i.equals(atomIndexOfBackboneN) && j.equals(atomIndexOfBackboneO)) {
                            // H bridge from backbone atom 'N' of residue a to backbone atom 'O' of residue b.
                            numPairContacts[MolContactInfo.HB]++;
                            // There can only be one of these so if we found it, simply update the distance.
                            minContactDistances[MolContactInfo.HB] = dist;
                        }

                        if(i.equals(atomIndexOfBackboneO) && j.equals(atomIndexOfBackboneN)) {
                            // H bridge from backbone atom 'O' of residue a to backbone atom 'N' of residue b.
                            numPairContacts[MolContactInfo.BH]++;
                            // There can only be one of these so if we found it, simply update the distance.
                            minContactDistances[MolContactInfo.BH] = dist;
                        }

                    }
                    
                    else if(x.isProteinAtom() && y.isLigandAtom()) {
                        // *************************** protein - ligand contact *************************
                        numTotalLigContactsPair++;

                        // Check the exact contact type
                        if(i <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - ligand contact
                            numPairContacts[MolContactInfo.BL]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BL] < 0) || dist < minContactDistances[MolContactInfo.BL]) {
                                minContactDistances[MolContactInfo.BL] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BL] = i;
                                contactAtomNumInResidueB[MolContactInfo.BL] = j;
                            }

                        }
                        else {
                            // to be precise, this is a side chainName - ligand contact
                            numPairContacts[MolContactInfo.CL]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CL] < 0) || dist < minContactDistances[MolContactInfo.CL]) {
                                minContactDistances[MolContactInfo.CL] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CL] = i;
                                contactAtomNumInResidueB[MolContactInfo.CL] = j;
                            }
                        }

                    }
                    else if(x.isLigandAtom() && y.isProteinAtom()) {
                        // *************************** ligand - protein contact *************************
                        numTotalLigContactsPair++;

                        // Check the exact contact type
                        if(j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a ligand - backbone contact
                            numPairContacts[MolContactInfo.LB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.LB] < 0) || dist < minContactDistances[MolContactInfo.LB]) {
                                minContactDistances[MolContactInfo.LB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.LB] = i;
                                contactAtomNumInResidueB[MolContactInfo.LB] = j;
                            }

                        }
                        else {
                            // to be precise, this is a ligand - side chainName contact
                            numPairContacts[MolContactInfo.LC]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.LC] < 0) || dist < minContactDistances[MolContactInfo.LC]) {
                                minContactDistances[MolContactInfo.LC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.LC] = i;
                                contactAtomNumInResidueB[MolContactInfo.LC] = j;
                            }
                        }
                            
                    }
                    else if(x.isLigandAtom() && y.isLigandAtom()) {
                        // *************************** ligand - ligand contact *************************
                        numTotalLigContactsPair++;

                        // no choices here, ligands have no sub type
                        numPairContacts[MolContactInfo.LL]++;
                        
                        // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                        if((minContactDistances[MolContactInfo.LL] < 0) || dist < minContactDistances[MolContactInfo.LL]) {
                            minContactDistances[MolContactInfo.LL] = dist;
                            contactAtomNumInResidueA[MolContactInfo.LL] = i;
                            contactAtomNumInResidueB[MolContactInfo.LL] = j;
                        }
                        

                    }
                    else {
                        // *************************** unknown contact, wtf? *************************
                        // This branch should never be hit because atoms of type OTHER are ignored while creating the list of Atom objects
                        System.out.println("WARNING: One of the atoms " + x.getPdbAtomNum() + " and " + y.getPdbAtomNum() + " is of type UNKNOWN. Bug?");
                    }                                                            
                }                
                else {
                    // No atom contact for these 2 atoms, but there could be contacts between others
                }
            }
        }


        // Iteration through all atoms of the two residues is done
        if(numPairContacts[MolContactInfo.TT] > 0) {
            result = new MolContactInfo(numPairContacts, minContactDistances, contactAtomNumInResidueA, contactAtomNumInResidueB, a, b, CAdist, numTotalLigContactsPair);
        }
        else {
            result = null;
        }
        
        return(result);         // This is null if no contact was detected
        
    }
 
    /**
     * Checks if a pi effect between the given X-H and aromatic ring occurs.
     * @param x Atom which is has a H.
     * @param h Atom which is bound to X.
     * @param acceptor ArrayList of Atom containing the atoms of the aromatic ring. Must be at least a five-ring.
     * @param type optional boolean(s). Index 0:true then X=N gets treated as sidechain N of Arg/Lys and X=C gets treated as Pro-CD. Index 1:
     *  true then the normal gets flipped. MAY BE DELETED.
     * @return double distance as 10th of Angstrom from X to ring midpoint (=> parameter 1). Returns -1 if no Pi effect occurs or not enough atoms
     * in ring
     */
    public static double calculateDistancePiEffect(Atom x, Atom h, ArrayList<Atom> acceptor, boolean... type) {
        double[] ringMidKoord = new double[3];
        double[] spanningVectorA = new double[3];
        double[] spanningVectorB = new double[3];
        double[] HXVector = new double[3];
        double[] HMidpointVector = new double[3];
        double[] normal = new double[3];
        double[] XMidpointVector = new double[3];
        double[] parameter = new double[4];
        
        if (x == null || h == null) {
            return -1;
        }
        
        if (acceptor.size() >= 5) {
            //calculate mid-point of ring
            ringMidKoord = PiEffectCalculations.calculateMidpointOfAtoms(acceptor);
            
            //calculate parameter 1: X-midpoint distance
            parameter[0] = PiEffectCalculations.calculateDistanceAtomToPoint(x, ringMidKoord);

            //calculate parameter 2: H-midpoint distance
            parameter[1] = PiEffectCalculations.calculateDistanceAtomToPoint(h, ringMidKoord);
            
            //calculate parameter 3: N-H-midpoint angle
            HXVector[0] = x.getCoordX() - h.getCoordX();
            HXVector[1] = x.getCoordY() - h.getCoordY();
            HXVector[2] = x.getCoordZ() - h.getCoordZ();
            
            HMidpointVector[0] = ringMidKoord[0] - h.getCoordX();
            HMidpointVector[1] = ringMidKoord[1] - h.getCoordY();
            HMidpointVector[2] = ringMidKoord[2] - h.getCoordZ();
            
            parameter[2] = PiEffectCalculations.calculateAngleBetw3DVecs(HXVector, HMidpointVector);
            
            //converte to degree
            parameter[2] = PiEffectCalculations.converteRadianToDegree(parameter[2]);

            //calculate normal vector
            spanningVectorA[0] = acceptor.get(1).getCoordX() - acceptor.get(0).getCoordX();
            spanningVectorA[1] = acceptor.get(1).getCoordY() - acceptor.get(0).getCoordY();
            spanningVectorA[2] = acceptor.get(1).getCoordZ() - acceptor.get(0).getCoordZ();

            spanningVectorB[0] = acceptor.get(2).getCoordX() - acceptor.get(0).getCoordX();
            spanningVectorB[1] = acceptor.get(2).getCoordY() - acceptor.get(0).getCoordY();
            spanningVectorB[2] = acceptor.get(2).getCoordZ() - acceptor.get(0).getCoordZ();
            
            normal = PiEffectCalculations.calculateNormalOfPlane(spanningVectorA, spanningVectorB);
            
            if (type.length == 2) {
                if (type[1]) {
                    normal[0] *= -1;
                    normal[1] *= -1;
                    normal[2] *= -1;
                }
            }
            
            //not in use due to change above
            //check if normal vector points to X-H group. Flip if not.
            //May be replaced with better calculation of normal vector
            //normal = PiEffectCalculations.checkDirectionNormal(normal, ringMidKoord, h);

            //calculate vector from N to ring mid-point
            XMidpointVector[0] = ringMidKoord[0] - x.getCoordX();
            XMidpointVector[1] = ringMidKoord[1] - x.getCoordY();
            XMidpointVector[2] = ringMidKoord[2] - x.getCoordZ();

            //calculate parameter 4: N-midpoint-normal angle
            parameter[3] = PiEffectCalculations.calculateAngleBetw3DVecs(XMidpointVector, normal);

            //converte to degree
            parameter[3] = PiEffectCalculations.converteRadianToDegree(parameter[3]);
           
            //determine if pi effect occurs
            //choose threshold by given X (altered if type[0]=true, see javadocs for more information)
            if ("N".equals(x.getAtomShortName())) {
                if (type.length > 0) {
                    if (type[0]) { // Lys sidechain N-H...Pi
                        if ((parameter[0] / 10) <= 4.0 && (parameter[1] / 10) <= 3.8 && parameter[2] >= 10 && parameter[3] <= 30) {
                            return parameter[0];
                        }
                    } else { // backbone N-H...Pi
                        if ((parameter[0] / 10) <= 4.3 && (parameter[1] / 10) <= 3.5 && parameter[2] >= 120 && parameter[3] <= 30) {
                            return parameter[0];
                        }
                    } 
                } else { // backbone N-H...Pi
                    if ((parameter[0] / 10) <= 4.3 && (parameter[1] / 10) <= 3.5 && parameter[2] >= 120 && parameter[3] <= 30) {
                        return parameter[0];
                    }
                }
            }
            
            //CA-HA...Pi and Pro-CD-Hd...Pi
            if("CA".equals(x.getAtomShortName()) || "CD".equals(x.getAtomShortName())) { 
                if ((parameter[0] / 10) <= 4.3 && (parameter[1] / 10) <= 3.8 && parameter[2] >= 120 && parameter[3] <= 30) {
                    return parameter[0];
                }
            }
            
            //S-H...Pi
            if ("SG".equals(x.getAtomShortName())) {
                if ((parameter[0] / 10) <= 4.0 && (parameter[1] / 10) <= 3.5 && parameter[2] >= 120 && parameter[3] <= 30) {
                    return parameter[0];
                } 
            }
            
            //O-H...Pi
            if (x.getAtomShortName().contains("OG") || "OH".equals(x.getAtomShortName())) {
                if ((parameter[0] / 10) <= 3.8 && (parameter[1] / 10) <= 3.5 && parameter[2] >= 120 && parameter[3] <= 30) {
                    return parameter[0];
                } 
            }
            
            if (type.length > 0) {
                if (type[0]) { // Arg sidechain N-H...Pi
                    if ("NZ".equals(x.getAtomShortName()) || "NE".equals(x.getAtomShortName()) || "NH1".equals(x.getAtomShortName()) || "NH2".equals(x.getAtomShortName())) {
                        if ((parameter[0] / 10) <= 4.0 && (parameter[1] / 10) <= 3.8 && parameter[2] >= 10 && parameter[3] <= 30) {
                                    return parameter[0];
                        }
                    }
                }
            }
        }
        // No pi effect occurs or not at least 5 atoms in acceptor (Warning has been printed previously).
        return -1;  
    }
    
    /**
     * Checks if a pi effect occurs for the specific group of CA-H...O=C effects (overloaded method).
     * @param ca Atom backbone C alpha
     * @param h Atom Hydrogen bound to ca
     * @param c Atom C of another residue
     * @param o Atom oxygen with double bond to c
     * @return double distance as 10th of Angstrom from CA to O (=> parameter 1). Returns -1 if no Pi effect occurs
     */
    public static double calculateDistancePiEffect(Atom ca, Atom h, Atom c, Atom o) {
        double[] parameter = new double[4];
        double[] hCAVector = new double[3];
        double[] hOVector = new double[3];
        double[] oHVector = new double[3];
        double[] oCVector = new double[3];
        
        //parameter 1
        parameter[0] = ca.distToAtom(o);
        
        //parameter 2
        parameter[1] = h.distToAtom(o);
        
        //parameter 3
        hCAVector[0] = ca.getCoordX() - h.getCoordX();
        hCAVector[1] = ca.getCoordY() - h.getCoordY();
        hCAVector[2] = ca.getCoordZ() - h.getCoordZ();
        
        hOVector[0] = o.getCoordX() - h.getCoordX();
        hOVector[1] = o.getCoordY() - h.getCoordY();
        hOVector[2] = o.getCoordZ() - h.getCoordZ();
        
        parameter[2] = PiEffectCalculations.calculateAngleBetw3DVecs(hCAVector, hOVector);
        
        //converte to degree
        parameter[2] = PiEffectCalculations.converteRadianToDegree(parameter[2]);
        
        //parameter 4
        oHVector[0] = h.getCoordX() - o.getCoordX();
        oHVector[1] = h.getCoordY() - o.getCoordY();
        oHVector[2] = h.getCoordZ() - o.getCoordZ();
        
        oCVector[0] = c.getCoordX() - o.getCoordX();
        oCVector[1] = c.getCoordY() - o.getCoordY();
        oCVector[2] = c.getCoordZ() - o.getCoordZ();
        
        parameter[3] = PiEffectCalculations.calculateAngleBetw3DVecs(oHVector, oCVector);
        
        //converte to degree
        parameter[3] = PiEffectCalculations.converteRadianToDegree(parameter[3]);
        
        //determine if pi effect occurs
        if ((parameter[0] / 10) <= 3.8 && (parameter[1] / 10) <= 3.3 && parameter[2] >= 120 && parameter[3] >= 90) {
            return parameter[0];
        }
        
        //No pi effect occurs.
        return -1;
    }
    
    /**
     * Calculates the planarity of an aromatic ring.
     * First a least-squares plane is calculated for the atoms of the aromatic ring
     * and then the RMSD from all those atoms to this plane is calculated. The lower
     * the RMSD the more planar the aromatic ring.
     * @param a residue for whom the aromatic ring planarity will be calculated
     * @return String containing information about the residue, the RMSD value, and the
     * atoms forming the aromatic ring.
     */
    public static String checkAromaticRingPlanarity(Residue a) {
        
        ArrayList<Atom> atoms_a = a.getAtoms();
        
        ArrayList<String> sidechainPiRings = new ArrayList<String>();
        sidechainPiRings.add("TRP");
        sidechainPiRings.add("TYR");
        sidechainPiRings.add("PHE");
        
        ArrayList<Atom> six_ring = new ArrayList<Atom>();
        ArrayList<Atom> five_ring = new ArrayList<Atom>(); //in case of TRP
        
        StringBuilder sb = new StringBuilder();
        
        
        // atm this is used to estimate how planar the aromatic rings actually are
        // (in a mathematical sense) to see if we have to worry about the calculation
        // of the normal vector for those aromatic rings
        if (sidechainPiRings.contains(a.getName3())) {           
            //get atoms of six-ring (and five-ring in case of Trp)
            six_ring.clear();
            five_ring.clear();
            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {                
                if (atoms_a.size() >= 11) {
                    for (Integer k = 5; k < 11; k++) {
                        six_ring.add(atoms_a.get(k));
                    }
                    double rmsd = PiEffectCalculations.calculateAromaticRingPlanarity(six_ring);
                    sb.append("[RMSD] ").append(a.getChainID()).append(" ").append(a.getPdbNum()).append(" ").append(a.getName3()).append(" ").append(String.valueOf(rmsd));
                    sb.append("[RMSD|ATOMS] ").append(six_ring.toString());
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() + ") contains not enough atoms.");
                }

            } else if ("TRP".equals(a.getName3())) {
                           
                //five_ring includes CG, CD1, NE1, CE2, CD2
                if (atoms_a.size() >= 9) {
                    for (Integer k = 5; k < 10; k++) {
                        five_ring.add(atoms_a.get(k));    
                    }
                    double rmsd = PiEffectCalculations.calculateAromaticRingPlanarity(five_ring);
                    sb.append("[RMSD] ").append(a.getChainID()).append(" ").append(a.getPdbNum()).append(" ").append(a.getName3()).append(" ").append(String.valueOf(rmsd));
                    sb.append("[RMSD|ATOMS] ").append(five_ring.toString());
                } else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms for 5-ring.");
                }
                
                //six_ring includes CD2, CE2, CE3, CZ2, CZ3, CH2
                if (atoms_a.size() >= 13) {
                    six_ring.add(atoms_a.get(7)); //CD2
                    for (Integer k = 9; k < 14; k++) {
                        six_ring.add(atoms_a.get(k));
                    }
                    double rmsd = PiEffectCalculations.calculateAromaticRingPlanarity(six_ring);
                    sb.append("[RMSD] ").append(a.getChainID()).append(" ").append(a.getPdbNum()).append(" ").append(a.getName3()).append(" ").append(String.valueOf(rmsd));
                    sb.append("[RMSD|ATOMS] ").append(six_ring.toString());
                } else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms for 6-ring.");
                }
            }
        }
        return sb.toString();
    }
    
    
    /** Calculate pi-effects between residue 'a' and 'b'.
     * Calculates pi-effects occurring between  residues of different chains.
     * The results can be integrated into the PPI detection.
     * @param a one of the residues of the residue pair
     * @param b one of the residues of the residue pair
     * @return A MolContactInfo object with information on the pi-effects between 'a' and 'b'.
     */
    public static MolContactInfo calculatePiEffects(Residue a, Residue b) {
        
        ArrayList<Atom> atoms_a = a.getAtoms();
        ArrayList<Atom> atoms_b = b.getAtoms();
        
        Integer CAdist = a.distTo(b);
        MolContactInfo result = null;
        ArrayList<Atom[]> atomAtomContacts = new ArrayList<Atom[]>();
        Atom[] donorAcceptor = new Atom[2];
        ArrayList<String> atomAtomContactType = new ArrayList<String>();
        ArrayList<String> sidechainOAAs = new ArrayList<String>();
        sidechainOAAs.add("ASP");
        sidechainOAAs.add("GLU");
        sidechainOAAs.add("ASN");
        sidechainOAAs.add("GLN");
        ArrayList<String> sidechainOHAAs = new ArrayList<String>();
        sidechainOHAAs.add("SER");
        sidechainOHAAs.add("THR");
        sidechainOHAAs.add("TYR");
        ArrayList<String> sidechainNHAAs = new ArrayList<String>();
        sidechainNHAAs.add("ARG");
        sidechainNHAAs.add("HIS");
        sidechainNHAAs.add("LYS");
        sidechainNHAAs.add("ASN");
        sidechainNHAAs.add("GLN");
        sidechainNHAAs.add("TRP");
        ArrayList<String> sidechainPiRings = new ArrayList<String>();
        sidechainPiRings.add("TRP");
        sidechainPiRings.add("TYR");
        sidechainPiRings.add("PHE");

        
        ArrayList<Atom> six_ring = new ArrayList<Atom>();
        ArrayList<Atom> five_ring = new ArrayList<Atom>(); //in case of TRP
        
        
        Integer[] numPairContacts = new Integer[Main.NUM_RESIDUE_PAIR_CONTACT_TYPES_ALTERNATIVE_MODEL];
        Integer numTotalLigContactsPair = 0;

        Integer[] minContactDistances = new Integer[numPairContacts.length];
        // Holds the minimal distances of contacts of the appropriate type (see numPairContacts, index 0 is unused)

        Integer[] contactAtomNumInResidueA = new Integer[numPairContacts.length];
        // Holds the number Atom x has in its residue a for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)

        Integer[] contactAtomNumInResidueB = new Integer[numPairContacts.length];
        // Holds the number Atom y has in its residue b for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)
        

        for(Integer k = 0; k < numPairContacts.length; k++) {      // init all arrays
            numPairContacts[k] = 0;
            minContactDistances[k] = -1;    // We are looking for the smallest distance >= 0 in this function so do NOT set '0' as the initial value or everything will get fucked up!
            contactAtomNumInResidueA[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            contactAtomNumInResidueB[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            // The initial value of '-1' is required for the atom index arrays because our index
            // starts at '0', but geom_neo treads a '0' in a line of <pdbid>.geo as 'no contact' because
            // it starts its index at '1'.
            // This problem is solved by the functions in MolContactInfo: they return (our_index + 1). This
            // means that:
            //   1) If no contact was detected, our_index is -1 and they return 0, which means 'no contact' to geom_neo.
            //   2) If a contact was detected, our_index is converted to the geom_neo index. :)
        }

        //non-canonical interactions (new implementation)
        int piDist = -1;
        if (! (a.isLigand() || b.isLigand())) {
        //residue b includes aromatic ring (acceptor) and res a is donor
        if (sidechainPiRings.contains(b.getName3())) {           
            //get atoms of six-ring (and five-ring in case of Trp)
            six_ring.clear();
            five_ring.clear();
            if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {                
                if (atoms_b.size() >= 11) {
                    for (Integer k = 5; k < 11; k++) {
                        six_ring.add(atoms_b.get(k));
                    }
                }
                else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() + ") contains not enough atoms.");
                }

            } else if ("TRP".equals(b.getName3())) {
                           
                //five_ring includes CG, CD1, NE1, CE2, CD2
                if (atoms_b.size() >= 10) {
                    for (Integer k = 5; k < 10; k++) {
                        five_ring.add(atoms_b.get(k));    
                    }
                } else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms for 5-ring.");
                }
                
                //six_ring includes CD2, CE2, CE3, CZ2, CZ3, CH2
                if (atoms_b.size() >= 14) {
                    six_ring.add(atoms_b.get(7)); //CD2
                    for (Integer k = 9; k < 14; k++) {
                        six_ring.add(atoms_b.get(k));
                    }
                } else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms for 6-ring.");
                }
            }
            
            //NHPI
            if (atoms_a.size() > 0 && a.getHydrogenAtoms().size() > 0) {
                if (! "PRO".equals(a.getName3())) { //Pro contains no backbone N-H
                    for (Atom h : a.getHydrogenAtoms()) {
                        if (h.getAtomName().replaceAll("\\s+","").equals("H")) {
                            piDist = (int)(calculateDistancePiEffect(atoms_a.get(0), h, six_ring) / 10);
                      
                            if ( piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.NHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.NHPI] < 0 || piDist > minContactDistances[MolContactInfo.NHPI]) {
                                    minContactDistances[MolContactInfo.NHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.NHPI] = 0; //backbone N
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New NHPI: " + atoms_a.get(0).toString() + "/" + h.toString());
                                atomAtomContactType.add("NHPI");
                                donorAcceptor[0] = atoms_a.get(0);
                                atomAtomContacts.add(donorAcceptor);
                            }
                     

                            piDist = (int)(calculateDistancePiEffect(atoms_a.get(0), h, six_ring, false, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.NHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.NHPI] < 0 || piDist > minContactDistances[MolContactInfo.NHPI]) {
                                    minContactDistances[MolContactInfo.NHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.NHPI] = 0; //backbone N
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 5; //CG of six_ring
                                        donorAcceptor[0] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[0] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New NHPI: " + atoms_a.get(0).toString() + "/" + h.toString());
                                atomAtomContactType.add("NHPI");
                                donorAcceptor[0] = atoms_a.get(0);
                                atomAtomContacts.add(donorAcceptor);
                            }

                            if ("TRP".equals(b.getName3())) {
                                piDist = (int)(calculateDistancePiEffect(a.getAtoms().get(0), h, five_ring) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.NHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.NHPI] < 0 || piDist > minContactDistances[MolContactInfo.NHPI]) {
                                        minContactDistances[MolContactInfo.NHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.NHPI] = 0; //backbone N
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New NHPI: " + atoms_a.get(0).toString() + "/" + h.toString());
                                    atomAtomContactType.add("NHPI");
                                    donorAcceptor[0] = atoms_a.get(0);
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }

                                piDist = (int)(calculateDistancePiEffect(a.getAtoms().get(0), h, five_ring, false, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.NHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.NHPI] < 0 || piDist > minContactDistances[MolContactInfo.NHPI]) {
                                        minContactDistances[MolContactInfo.NHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.NHPI] = 0; //backbone N
                                        contactAtomNumInResidueB[MolContactInfo.NHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New NHPI: " + atoms_a.get(0).toString() + "/" + h.toString());
                                    atomAtomContactType.add("NHPI");
                                    donorAcceptor[0] = atoms_a.get(0);
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }

                            }
                        }
                    }
                }
            }
            else {
                if (! (atoms_a.size() > 0)) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains no atoms."
                        + " Continue search in next residues.");
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains no hydrogen."
                        + " Continue search in next residues.");
                }
            }

            //CNHPI
            if ("LYS".equals(a.getName3())) {
                if (atoms_a.size() >= 9 && a.getHydrogenAtoms().size() > 0) {
                    for (Atom hz : a.getHydrogenAtoms()) {
                        if (hz.getAtomName().contains("HZ")) {
                            piDist = (int)(calculateDistancePiEffect(atoms_a.get(8), hz, six_ring, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CNHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.CNHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CNHPI] = 8; //sidechain N
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }                             
//                                System.out.println("New CNHPI: " + atoms_a.get(8).toString() + "/" + hz.toString());
                                atomAtomContactType.add("CNHPI");
                                donorAcceptor[0] = atoms_a.get(8);
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                            piDist = (int)(calculateDistancePiEffect(atoms_a.get(8), hz, six_ring, true, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CNHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.CNHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CNHPI] = 8; //sidechain N
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CNHPI: " + atoms_a.get(8).toString() + "/" + hz.toString());
                                atomAtomContactType.add("CNHPI");
                                donorAcceptor[0] = atoms_a.get(8);
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                            if ("TRP".equals(b.getName3())) {
                                piDist = (int)(calculateDistancePiEffect(atoms_a.get(8), hz, five_ring, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.CNHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.CNHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.CNHPI] = 8; //sidechain N
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New CNHPI: " + atoms_a.get(8).toString() + "/" + hz.toString());
                                    atomAtomContactType.add("CNHPI");
                                    donorAcceptor[0] = atoms_a.get(8);
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                                
                                piDist = (int)(calculateDistancePiEffect(atoms_a.get(8), hz, five_ring, true, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.CNHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.CNHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.CNHPI] = 8; //sidechain N
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New CNHPI: " + atoms_a.get(8).toString() + "/" + hz.toString());
                                    atomAtomContactType.add("CNHPI");
                                    donorAcceptor[0] = atoms_a.get(8);
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                            }
                        }
                    }
                }
                else {
                    if (! (atoms_a.size() >= 9)) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms."
                        + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains no hydrogen."
                            + " Continue search in next residues.");
                    }
                }
            }
            
            if ("ARG".equals(a.getName3())) {
                Atom argN;
                
                if (atoms_a.size() >= 11 && a.getHydrogenAtoms().size() > 0) {
                    for (Atom argH : a.getHydrogenAtoms()) {
                        //check if H is bond to a sidechain N
                        if (argH.getAtomName().contains("HE") || argH.getAtomName().contains("HH")) {
                            argN = atoms_a.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                            
                            piDist = (int)(calculateDistancePiEffect(argN, argH, six_ring, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CNHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.CNHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CNHPI] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CNHPI: " + argN.toString() + "/" + argH.toString());
                                atomAtomContactType.add("CNHPI");
                                donorAcceptor[0] = atoms_a.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                            piDist = (int)(calculateDistancePiEffect(argN, argH, six_ring, true, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CNHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.CNHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CNHPI] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CNHPI: " + argN.toString() + "/" + argH.toString());
                                atomAtomContactType.add("CNHPI");
                                donorAcceptor[0] = atoms_a.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                            if ("TRP".equals(b.getName3())) {
                                piDist = (int)(calculateDistancePiEffect(argN, argH, five_ring, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.CNHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.CNHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.CNHPI] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New CNHPI: " + argN.toString() + "/" + argH.toString());
                                    atomAtomContactType.add("CNHPI");
                                    donorAcceptor[0] = atoms_a.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                                
                                piDist = (int)(calculateDistancePiEffect(argN, argH, five_ring, true, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.CNHPI]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.CNHPI] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.CNHPI] = piDist;
                                        contactAtomNumInResidueA[MolContactInfo.CNHPI] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                        contactAtomNumInResidueB[MolContactInfo.CNHPI] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New CNHPI: " + argN.toString() + "/" + argH.toString());
                                    atomAtomContactType.add("CNHPI");
                                    donorAcceptor[0] = atoms_a.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                    donorAcceptor[1] = atoms_b.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                            }
                        }
                    }
                }
                else {
                    if (! (atoms_a.size() >= 10)) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms."
                        + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains no hydrogen."
                            + " Continue search in next residues.");
                    }
                }
            }
            
            //CAHPI            
            if ((atoms_a.size() > 0 && a.getHydrogenAtoms().size() > 0)) {
                    
                Atom ca = null;
                for (Atom c : atoms_a) {
                    if (c.getAtomName().contains("CA")) {
                        ca = c;
                    }
                }
                if (ca == null) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                        + " Continue search in next residues.");
                }
                else {
                
                for (Atom ha : a.getHydrogenAtoms()) {
                    if (ha.getAtomName().contains("HA")) {

                        piDist = (int)(calculateDistancePiEffect(ca, ha, six_ring) / 10);

                        if ( piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.CAHPI]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.CAHPI] < 0 || piDist > minContactDistances[MolContactInfo.CAHPI]) {
                                minContactDistances[MolContactInfo.CAHPI] = piDist;
                                contactAtomNumInResidueA[MolContactInfo.CAHPI] = 1; //CA
                                if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                    contactAtomNumInResidueB[MolContactInfo.CAHPI] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_b.get(5);
                                }
                                else {
                                    contactAtomNumInResidueB[MolContactInfo.CAHPI] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_b.get(7);
                                }
                            }
//                            System.out.println("New CAHPI: " + ca.toString() + "/" + ha.toString());
                            atomAtomContactType.add("CAHPI");
                            donorAcceptor[0] = atoms_a.get(1);
                            atomAtomContacts.add(donorAcceptor);
                        }
                    
                            piDist = (int)(calculateDistancePiEffect(ca, ha, six_ring, false, true) / 10);
                    
                            if ( piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CAHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CAHPI] < 0 || piDist > minContactDistances[MolContactInfo.CAHPI]) {
                                    minContactDistances[MolContactInfo.CAHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CAHPI] = 1; //CA
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CAHPI: " + ca.toString() + "/" + ha.toString());
                                atomAtomContactType.add("CAHPI");
                                donorAcceptor[0] = atoms_a.get(1);
                                atomAtomContacts.add(donorAcceptor);
                            }

                    if ("TRP".equals(b.getName3())) {
                            piDist = (int)(calculateDistancePiEffect(ca, ha, five_ring) / 10);
                            if ( piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CAHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CAHPI] < 0 || piDist > minContactDistances[MolContactInfo.CAHPI]) {
                                    minContactDistances[MolContactInfo.CAHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CAHPI] = 1; //CA
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CAHPI: " + ca.toString() + "/" + ha.toString());
                                atomAtomContactType.add("CAHPI");
                                donorAcceptor[0] = atoms_a.get(1);
                                atomAtomContacts.add(donorAcceptor);
                            }

                            piDist = (int)(calculateDistancePiEffect(ca, ha, five_ring, false, true) / 10);
                            if ( piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CAHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CAHPI] < 0 || piDist > minContactDistances[MolContactInfo.CAHPI]) {
                                    minContactDistances[MolContactInfo.CAHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CAHPI] = 1; //CA
                                    if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_b.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CAHPI] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                }
//                                System.out.println("New CAHPI: " + ca.toString() + "/" + ha.toString());
                                atomAtomContactType.add("CAHPI");
                                donorAcceptor[0] = atoms_a.get(1);
                                atomAtomContacts.add(donorAcceptor);
                            }
                        }
                    }
                }
            }
            }
            else {
                DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (No CA)."
                        + " Continue search in next residues.");
            }
            //PROCDHPI
            if (a.getName3().equals("PRO")) {
                if (atoms_a.size() >= 7 && a.getHydrogenAtoms().size() > 0) {
                    for (Atom hd : a.getHydrogenAtoms()) {
                        if (hd.getAtomName().contains("HD")) {
                            piDist = (int)(calculateDistancePiEffect(atoms_a.get(6), hd, six_ring) / 10);
                        if ( piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PROCDHPI]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PROCDHPI] < 0 || piDist > minContactDistances[MolContactInfo.PROCDHPI]) {
                                minContactDistances[MolContactInfo.PROCDHPI] = piDist;
                                contactAtomNumInResidueA[MolContactInfo.PROCDHPI] = 6; //CD
                                if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_b.get(5);
                                }
                                else {
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_b.get(7);
                                }
                            }
//                            System.out.println("New PROCDHPI: " + atoms_a.get(6).toString() + "/" + hd.toString());
                            atomAtomContactType.add("PROCDHPI");
                            donorAcceptor[0] = atoms_a.get(6);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(atoms_a.get(6), hd, six_ring, false, true) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PROCDHPI]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PROCDHPI] < 0 || piDist > minContactDistances[MolContactInfo.PROCDHPI]) {
                                minContactDistances[MolContactInfo.PROCDHPI] = piDist;
                                contactAtomNumInResidueA[MolContactInfo.PROCDHPI] = 6; //CD
                                if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_b.get(5);
                                }
                                else {
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_b.get(7);
                                }
                            }
//                            System.out.println("New PROCDHPI: " + atoms_a.get(6).toString() + "/" + hd.toString());
                            atomAtomContactType.add("PROCDHPI");
                            donorAcceptor[0] = atoms_a.get(6);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        if ("TRP".equals(b.getName3())) {
                            piDist = (int)(calculateDistancePiEffect(a.getAtoms().get(6), hd, five_ring) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PROCDHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PROCDHPI] < 0 || piDist > minContactDistances[MolContactInfo.PROCDHPI]) {
                                    minContactDistances[MolContactInfo.PROCDHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.PROCDHPI] = 6; //CD
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 5; //CG of five_ring
                                }
//                                System.out.println("New PROCDHPI: " + atoms_a.get(6).toString() + "/" + hd.toString());
                                atomAtomContactType.add("PROCDHPI");
                                donorAcceptor[0] = atoms_a.get(6);
                                donorAcceptor[1] = atoms_b.get(5);
                                atomAtomContacts.add(donorAcceptor);
                            }

                            piDist = (int)(calculateDistancePiEffect(a.getAtoms().get(6), hd, five_ring, false, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PROCDHPI]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PROCDHPI] < 0 || piDist > minContactDistances[MolContactInfo.PROCDHPI]) {
                                    minContactDistances[MolContactInfo.PROCDHPI] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.PROCDHPI] = 6; //CD
                                    contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = 5; //CG of five_ring
                                }
//                                System.out.println("New PROCDHPI: " + atoms_a.get(6).toString() + "/" + hd.toString());
                                atomAtomContactType.add("PROCDHPI");
                                donorAcceptor[0] = atoms_a.get(6);
                                donorAcceptor[1] = atoms_b.get(5);
                                atomAtomContacts.add(donorAcceptor);
                            }

                        }
                    }
                }
            }
                else {
                    if (! (atoms_a.size() >= 7)) {
                        DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms."
                            + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogens."
                            + " Continue search in next residues.");
                    }
                }
            }
            
            //SHPI
            if (a.getName3().equals("CYS")) { 
                if (atoms_a.size() >= 6 && a.getHydrogenAtoms().size() > 0) {
                    for(Atom hg : a.getHydrogenAtoms()) {
                        if (hg.getAtomName().contains("HG")) {
                    piDist = (int)(calculateDistancePiEffect(atoms_a.get(5), hg, six_ring) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.SHPI]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.SHPI] < 0 || piDist > minContactDistances[MolContactInfo.SHPI]) {
                            minContactDistances[MolContactInfo.SHPI] = piDist;
                            contactAtomNumInResidueA[MolContactInfo.SHPI] = 5;
                            if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_b.get(5);
                            }
                            else {
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_b.get(7);
                            }
                        }
//                        System.out.println("New SHPI: " + atoms_a.get(5).toString() + "/" + hg.toString());
                        atomAtomContactType.add("SHPI");
                        donorAcceptor[0] = atoms_a.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    piDist = (int)(calculateDistancePiEffect(atoms_a.get(5), hg, six_ring, false, true) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.SHPI]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.SHPI] < 0 || piDist > minContactDistances[MolContactInfo.SHPI]) {
                            minContactDistances[MolContactInfo.SHPI] = piDist;
                            contactAtomNumInResidueA[MolContactInfo.SHPI] = 5;
                            if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_b.get(5);
                            }
                            else {
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_b.get(7);
                            }
                        }
//                        System.out.println("New SHPI: " + atoms_a.get(5).toString() + "/" + hg.toString());
                        atomAtomContactType.add("SHPI");
                        donorAcceptor[0] = atoms_a.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    if ("TRP".equals(b.getName3())) {
                        piDist = (int)(calculateDistancePiEffect(atoms_a.get(5), hg, five_ring) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.SHPI]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.SHPI] < 0 || piDist > minContactDistances[MolContactInfo.SHPI]) {
                                minContactDistances[MolContactInfo.SHPI] = piDist;
                                contactAtomNumInResidueA[MolContactInfo.SHPI] = 5;
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 5; //CG of five_ring
                            }
//                            System.out.println("New SHPI: " + atoms_a.get(5).toString() + "/" + hg.toString());
                            atomAtomContactType.add("SHPI");
                            donorAcceptor[0] = atoms_a.get(5);
                            donorAcceptor[1] = atoms_b.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(atoms_a.get(5), hg, five_ring, false, true) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.SHPI]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.SHPI] < 0 || piDist > minContactDistances[MolContactInfo.SHPI]) {
                                minContactDistances[MolContactInfo.SHPI] = piDist;
                                contactAtomNumInResidueA[MolContactInfo.SHPI] = 5;
                                contactAtomNumInResidueB[MolContactInfo.SHPI] = 5; //CG of five_ring
                            }
//                            System.out.println("New SHPI: " + atoms_a.get(5).toString() + "/" + hg.toString());
                            atomAtomContactType.add("SHPI");
                            donorAcceptor[0] = atoms_a.get(5);
                            donorAcceptor[1] = atoms_b.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }
                    }   
                }
                    }
                }
                else {
                    if (! (atoms_a.size() >= 5)) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no SG)."
                        + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogen atoms."
                            + " Continue search in next residues.");
                    }
                }         
            }
            
            //XOHPI
            if (sidechainOHAAs.contains(a.getName3())) {
                Atom OHAA_X = null;
                Atom OHAA_H = null;
                
                if ("SER".equals(a.getName3())) {
                    if (atoms_a.size() >= 6 && a.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_a.get(5);
                        for (Atom h : a.getHydrogenAtoms()){
                            if (h.getAtomName().contains("HG")) {
                                OHAA_H = h;
                            }
                        }
                    }
                    else {
                        if (! (atoms_a.size() >= 6)) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no OG).");
                        }
                        else {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogens (no HG).");
                        }
                    }
                }
                else if ("THR".equals(a.getName3())) {
                    if (atoms_a.size() >= 6 && a.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_a.get(5);
                        for (Atom h : atoms_a){
                            if (h.getAtomName().contains("HG1")) {
                                OHAA_H = h;
                            }
                        }
                    }
                    else {
                        if (! (atoms_a.size() >= 6)) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no OG1).");
                        }
                        else {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogens (no HG1).");
                        }
                    }
                }
                else {
                    if (atoms_a.size() >= 12 && a.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_a.get(11);
                        for (Atom h : a.getHydrogenAtoms()) {
                            if (h.getAtomName().contains("HH")) {
                                OHAA_H = h;
                            }
                        }
                    }
                    else {
                        if (! (atoms_a.size() >= 12)) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no OH).");
                        }
                        else {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogens (no HH).");
                        }
                    }
                }
                
                piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, six_ring) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.XOHPI]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.XOHPI] < 0 || piDist > minContactDistances[MolContactInfo.XOHPI]) {
                        minContactDistances[MolContactInfo.XOHPI] = piDist;
                        if ("TYR".equals(a.getName3())) {
                            contactAtomNumInResidueA[MolContactInfo.XOHPI] = 11;
                            donorAcceptor[0] = atoms_a.get(11);
                        }
                        else {
                            contactAtomNumInResidueA[MolContactInfo.XOHPI] = 5;
                            donorAcceptor[0] = atoms_a.get(5);
                        }
                        if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 5; //CG of six_ring
                            donorAcceptor[1] = atoms_b.get(5);
                        }
                        else {
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 7; //CD2 of six_ring
                            donorAcceptor[1] = atoms_b.get(7);
                        }
                    }
//                    System.out.println("New XOHPI: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                    atomAtomContactType.add("XOHPI");
                    atomAtomContacts.add(donorAcceptor);
                }

                piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, six_ring, false, true) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.XOHPI]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.XOHPI] < 0 || piDist > minContactDistances[MolContactInfo.XOHPI]) {
                        minContactDistances[MolContactInfo.XOHPI] = piDist;
                        if ("TYR".equals(a.getName3())) {
                            contactAtomNumInResidueA[MolContactInfo.XOHPI] = 11;
                            donorAcceptor[0] = atoms_a.get(11);
                        }
                        else {
                            contactAtomNumInResidueA[MolContactInfo.XOHPI] = 5;
                            donorAcceptor[0] = atoms_a.get(5);
                        }
                        if ("TYR".equals(b.getName3()) || "PHE".equals(b.getName3())) {
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 5; //CG of six_ring
                            donorAcceptor[1] = atoms_b.get(5);
                        }
                        else {
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 7; //CD2 of six_ring
                            donorAcceptor[1] = atoms_b.get(7);
                        }
                    }
//                    System.out.println("New XOHPI: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                    atomAtomContactType.add("XOHPI");
                    atomAtomContacts.add(donorAcceptor);
                }

                if ("TRP".equals(b.getName3())) {
                    piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, five_ring) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.XOHPI]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.XOHPI] < 0 || piDist > minContactDistances[MolContactInfo.XOHPI]) {
                            minContactDistances[MolContactInfo.XOHPI] = piDist;
                            if ("TYR".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.XOHPI] = 11;
                                donorAcceptor[0] = atoms_a.get(11);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.XOHPI] = 5;
                                donorAcceptor[0] = atoms_a.get(5);
                            }
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 5; //CG of five_ring
                        }
//                        System.out.println("New XOHPI: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                        atomAtomContactType.add("");
                        donorAcceptor[1] = atoms_b.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, five_ring, false, true) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.XOHPI]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.XOHPI] < 0 || piDist > minContactDistances[MolContactInfo.XOHPI]) {
                            minContactDistances[MolContactInfo.XOHPI] = piDist;
                            if ("TYR".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.XOHPI] = 11;
                                donorAcceptor[0] = atoms_a.get(11);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.XOHPI] = 5;
                                donorAcceptor[0] = atoms_a.get(5);
                            }
                            contactAtomNumInResidueB[MolContactInfo.XOHPI] = 5; //CG of five_ring
                        }
//                        System.out.println("New XOHPI: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                        atomAtomContactType.add("");
                        donorAcceptor[1] = atoms_b.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }
                }   
            }
        }
        
        //C=O in sidechain of aminoacid b acts as acceptor
        //CCAHCO
        if (sidechainOAAs.contains(b.getName3())) {
            //check for a valid res b
            if ( (("GLU".equals(b.getName3()) || "GLN".equals(b.getName3())) && atoms_b.size() >= 8) ||
                    (("ASP".equals(b.getName3()) || "ASN".equals(b.getName3())) && atoms_b.size() >= 7) ) {
                if ( ("GLU".equals(b.getName3()) || "GLN".equals(b.getName3()) && atoms_b.get(6).getAtomShortName().equals("CD") && 
                        atoms_b.get(7).getAtomShortName().equals("OE1")) ||
                       (("ASP".equals(b.getName3()) || "ASN".equals(b.getName3())) && atoms_b.get(5).getAtomShortName().equals("CG") &&
                        atoms_b.get(6).getAtomShortName().equals("OD1")) ) {  
            
                    if ((atoms_a.size() >= 2 && a.getHydrogenAtoms().size() > 0) || ("PRO".equals(a.getName3()) && atoms_a.size() >= 2 && a.getHydrogenAtoms().size() > 0)) {
                        
                        Atom ca = null;
                        for (Atom c : atoms_a) {
                            if (c.getAtomName().contains("CA")) {
                                ca = c;
                            }
                        }
                        if (ca == null) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                                + " Continue search in next residues.");
                        }
                        else {
                        for (Atom ha : a.getHydrogenAtoms()) {
                            if (ha.getAtomName().contains("HA")) {                  
                            
                                if ("GLU".equals(b.getName3()) || "GLN".equals(b.getName3())) {
                                    piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_b.get(6), atoms_b.get(7)) / 10);
                                }
                                else {
                                    piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_b.get(5), atoms_b.get(6)) / 10);
                                }
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CCAHCO]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CCAHCO] < 0 || piDist > minContactDistances[MolContactInfo.CCAHCO]) {
                                    minContactDistances[MolContactInfo.CCAHCO] = piDist;
                                    contactAtomNumInResidueA[MolContactInfo.CCAHCO] = 1; //CA
                                    if ("GLU".equals(b.getName3()) || "GLN".equals(b.getName3())) {
                                        contactAtomNumInResidueB[MolContactInfo.CCAHCO] = 7; //OE1
                                        donorAcceptor[1] = atoms_b.get(7);
                                    }
                                    else {
                                        contactAtomNumInResidueB[MolContactInfo.CCAHCO] = 6; //OD1
                                        donorAcceptor[1] = atoms_b.get(6);
                                    }
                                }
//                                System.out.println("New CCAHCO: " + ca.toString() + "/" + ha.toString());
                                atomAtomContactType.add("CCAHCO");
                                donorAcceptor[0] = atoms_a.get(1);
                                atomAtomContacts.add(donorAcceptor);
                            }
                        }
                        }
                    }
                    }
                    else {
                        if ((! (atoms_a.size() > 2)) || "PHE".equals(a.getName3()) && !(atoms_a.size() > 2)) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                                + " Continue search in next residues.");
                        }
                        else {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough hydrogen (no HA)."
                                + " Continue search in next residues.");
                        }
                    }
            
                }
                else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not the expected atom names for C=O."
                        + " Continue seach in next residues.");
                }
            }
            else {
                DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no C=O)."
                    + " Continue search in next residues.");
            }
        }
        
        //C=O in backbone of aminoacid b acts as acceptor
        //BCAHCO
        if (atoms_b.size() >= 4 && atoms_a.size() > 0 && a.getHydrogenAtoms().size() > 0) {
            
            if (atoms_b.get(2).getAtomShortName().equals("C") && atoms_b.get(3).getAtomShortName().equals("O")) {
                
                Atom ca = null;
                for (Atom c : atoms_a) {
                    if (c.getAtomName().contains("CA")) {
                        ca = c;
                    }
                }
                
                if (ca == null) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                        + " Continue search in next residues.");
                }
                else {
                for (Atom ha : a.getHydrogenAtoms()) {
                    if (ha.getAtomName().contains("HA")) {
                
                piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_b.get(2), atoms_b.get(3)) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.BCAHCO]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.BCAHCO] < 0 || piDist > minContactDistances[MolContactInfo.BCAHCO]) {
                        minContactDistances[MolContactInfo.BCAHCO] = piDist;
                        contactAtomNumInResidueA[MolContactInfo.BCAHCO] = 1; //CA
                        contactAtomNumInResidueB[MolContactInfo.BCAHCO] = 3; //O         
                    }
//                    System.out.println("New BCAHCO: " + ca.toString() + "/" + ha.toString());
                    atomAtomContactType.add("BCAHCO");
                    donorAcceptor[0] = atoms_a.get(1);
                    donorAcceptor[1] = atoms_b.get(3);
                    atomAtomContacts.add(donorAcceptor);
                }
                
                    }
                }
            }
            }
            else {
                if (! (atoms_b.get(2).getAtomShortName().equals("C") && atoms_b.get(3).getAtomShortName().equals("O"))) {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") "
                            + "contains not the expected atom name (C or O)"
                            + " for calculation of BCAHCO. Continue search in next residue.");
                }
                else { // no HA at expected position
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") "
                            + "contains not the expected atom name (HA or HA2 in case of GLY)"
                            + " for calculation of BCAHCO. Continue search in next residue.");
                }
            }
        }
        else {
            if (! (atoms_a.size() >= 2)) {
                DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") "
                        + "contains not enough atoms (no CA for BCAHCO)."
                        + " Continue search in next residue.");
            }
            else if (! (a.getHydrogenAtoms().size() >= 2)) {
                DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") "
                        + "contains not enough hydrogens (no HA for BCAHCO)."
                        + " Continue search in next residue.");
            }
            else { // ! atoms_b.size() >= 4
                DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") "
                        + "ontains not enough atoms (no C=O for BCAHCO)."
                        + " Continue search in next residue.");
            }
        }    
        
        //residue a includes aromatic ring and res b is donor
        if (sidechainPiRings.contains(a.getName3())) {
            six_ring.clear();
            five_ring.clear();
            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                if (atoms_a.size() >= 11) {
                    //CG, CD1, CD2, CE1, CE2, CZ
                    for (Integer k = 5; k < 11; k++) {
                        six_ring.add(atoms_a.get(k));
                    }
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (six ring).");
                }

            } else if ("TRP".equals(a.getName3())) {
                
                //five_ring includes CG, CD1, NE1, CE2, CD2
                if (atoms_a.size() >= 10) {
                    for (Integer k = 5; k < 10; k++) {
                        five_ring.add(atoms_a.get(k));      
                    }
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (five ring).");
                }
                
                //six_ring includes CD2, CE2, CE3, CZ2, CZ3, CH2
                if (atoms_a.size() >= 14) { 
                    six_ring.add(atoms_a.get(7)); //CD2
                    for (Integer k = 9; k < 14; k++) {
                        six_ring.add(atoms_a.get(k));
                    }
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (six ring).");
                }
            }
            
            //PINH
            if (atoms_b.size() > 0 && b.getHydrogenAtoms().size() > 0) {
                if (! "PRO".equals(b.getName3())) { //Pro contains no backbone N-H
                    for (Atom h : b.getHydrogenAtoms()) {
                        if (h.getAtomName().replaceAll("\\s+","").equals("H")) {
                        piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(0), h, six_ring) / 10);
                      
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PINH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PINH] < 0 || piDist > minContactDistances[MolContactInfo.PINH]) {
                            minContactDistances[MolContactInfo.PINH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PINH] = 0; //backbone N
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PINH: " + b.getAtoms().get(0).toString() + "/" + h.toString());
                        atomAtomContactType.add("PINH");
                        donorAcceptor[0] = atoms_b.get(0);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(0), h, six_ring, false, true) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PINH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PINH] < 0 || piDist > minContactDistances[MolContactInfo.PINH]) {
                            minContactDistances[MolContactInfo.PINH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PINH] = 0; //backbone N
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PINH: " + b.getAtoms().get(0).toString() + "/" + h.toString());
                        atomAtomContactType.add("PINH");
                        donorAcceptor[0] = atoms_b.get(0);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    if ("TRP".equals(b.getName3())) {
                        piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(0), h, five_ring) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PINH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PINH] < 0 || piDist > minContactDistances[MolContactInfo.PINH]) {
                                minContactDistances[MolContactInfo.PINH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PINH] = 0; //backbone N
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 5; //CG of five_ring
                            }
//                            System.out.println("New PINH: " + b.getAtoms().get(0).toString() + "/" + h.toString());
                            atomAtomContactType.add("PINH");
                            donorAcceptor[0] = atoms_b.get(0);
                            donorAcceptor[1] = atoms_a.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(0), h, five_ring, false, true) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PINH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PINH] < 0 || piDist > minContactDistances[MolContactInfo.PINH]) {
                                minContactDistances[MolContactInfo.PINH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PINH] = 0; //backbone N
                                contactAtomNumInResidueA[MolContactInfo.PINH] = 5; //CG of five_ring
                            }
//                            System.out.println("New PINH: " + b.getAtoms().get(0).toString() + "/" + h.toString());
                            atomAtomContactType.add("PINH");
                            donorAcceptor[0] = atoms_b.get(0);
                            donorAcceptor[1] = atoms_a.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }
                    }
                      }
                    }
                }
            }
            else {
                if (! (atoms_b.size() > 0)) {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains no atoms."
                        + " Continue search in next residues.");
                    }
                else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains no hydrogens."
                        + " Continue search in next residues.");
                }
            }
            
            //PICNH
            if ("LYS".equals(b.getName3())) {
                if (atoms_b.size() >= 9 && b.getHydrogenAtoms().size() > 0) {
                   for (Atom hz : b.getHydrogenAtoms()) {
                       if (hz.getAtomName().contains("HZ")) {
                           piDist = (int)(calculateDistancePiEffect(atoms_b.get(8), hz, six_ring, true) / 10);
                           if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PICNH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.PICNH]) {
                                    minContactDistances[MolContactInfo.PICNH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PICNH] = 8; //sidechain N
                                    if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_a.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_a.get(7);
                                    }
                                }
//                                System.out.println("New PICNH: " + b.getAtoms().get(8).toString() + "/" + hz.toString());
                                atomAtomContactType.add("PICNH");
                                donorAcceptor[0] = atoms_b.get(8);
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                            piDist = (int)(calculateDistancePiEffect(atoms_b.get(8), hz, six_ring, true, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PICNH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.PICNH]) {
                                    minContactDistances[MolContactInfo.PICNH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PICNH] = 8; //sidechain N
                                    if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_a.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_a.get(7);
                                    }
                                }
//                                System.out.println("New PICNH: " + b.getAtoms().get(8).toString() + "/" + hz.toString());
                                atomAtomContactType.add("PICNH");
                                donorAcceptor[0] = atoms_b.get(8);
                                atomAtomContacts.add(donorAcceptor);
                            }

                            if ("TRP".equals(a.getName3())) {
                                piDist = (int)(calculateDistancePiEffect(atoms_b.get(8), hz, five_ring, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.PICNH]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.PICNH]) {
                                        minContactDistances[MolContactInfo.PICNH] = piDist;
                                        contactAtomNumInResidueB[MolContactInfo.PICNH] = 8; //sidechain N
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New PICNH: " + b.getAtoms().get(8).toString() + "/" + hz.toString());
                                    atomAtomContactType.add("PICNH");
                                    donorAcceptor[0] = atoms_b.get(8);
                                    donorAcceptor[1] = atoms_a.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                                
                                piDist = (int)(calculateDistancePiEffect(atoms_b.get(8), hz, five_ring, true, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.PICNH]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.PICNH]) {
                                        minContactDistances[MolContactInfo.PICNH] = piDist;
                                        contactAtomNumInResidueB[MolContactInfo.PICNH] = 8; //sidechain N
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New PICNH: " + b.getAtoms().get(8).toString() + "/" + hz.toString());
                                    atomAtomContactType.add("PICNH");
                                    donorAcceptor[0] = atoms_b.get(8);
                                    donorAcceptor[1] = atoms_a.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                            }
                        }        
                   }
                }
                else {
                    if (! (atoms_b.size() >= 9)) {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no NZ)."
                            + " Continue search in next residues.");
                        }
                    else {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains no hydrogens."
                            + " Continue search in next residues.");
                    }
                }
            }
                
            if ("ARG".equals(b.getName3())) {
                Atom argN;

                if (atoms_b.size() >= 11 && b.getHydrogenAtoms().size() > 0) {
                    for (Atom argH : b.getHydrogenAtoms()) {

                        //check if H is bond to a sidechain N
                        if (argH.getAtomName().contains("HE") || argH.getAtomName().contains("HH")) {
                            if (argH.getAtomName().contains("HE")) {
                                argN = atoms_b.get(7);
                            }
                            else if (argH.getAtomName().contains("HH1")) {
                                argN = atoms_b.get(9);
                            }
                            else {
                                argN = atoms_b.get(10);
                            }

                            piDist = (int)(calculateDistancePiEffect(argN, argH, six_ring, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PICNH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.PICNH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PICNH] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                    if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_a.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_a.get(7);
                                    }
                                }
//                                System.out.println("New PICNH: " + argN.toString() + "/" + argH.toString());
                                atomAtomContactType.add("PICNH");
                                donorAcceptor[0] = atoms_b.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                atomAtomContacts.add(donorAcceptor);
                            }

                            piDist = (int)(calculateDistancePiEffect(argN, argH, six_ring, true, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PICNH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                    minContactDistances[MolContactInfo.PICNH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PICNH] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                    if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of six_ring
                                        donorAcceptor[1] = atoms_a.get(5);
                                    }
                                    else {
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 7; //CD2 of six_ring
                                        donorAcceptor[1] = atoms_a.get(7);
                                    }
                                }
//                                System.out.println("New PICNH: " + argN.toString() + "/" + argH.toString());
                                atomAtomContactType.add("PICNH");
                                donorAcceptor[0] = atoms_b.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                atomAtomContacts.add(donorAcceptor);
                            }

                            if ("TRP".equals(b.getName3())) {
                                piDist = (int)(calculateDistancePiEffect(argN, argH, five_ring, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.PICNH]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.PICNH] = piDist;
                                        contactAtomNumInResidueB[MolContactInfo.PICNH] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New PICNH: " + argN.toString() + "/" + argH.toString());
                                    atomAtomContactType.add("PICNH");
                                    donorAcceptor[0] = atoms_b.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                    donorAcceptor[1] = atoms_a.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }

                                piDist = (int)(calculateDistancePiEffect(argN, argH, five_ring, true, true) / 10);
                                if (piDist > 0) {
                                    numPairContacts[MolContactInfo.TT]++;
                                    numPairContacts[MolContactInfo.PICNH]++;
                                    donorAcceptor = new Atom[2];
                                    if (minContactDistances[MolContactInfo.PICNH] < 0 || piDist > minContactDistances[MolContactInfo.CNHPI]) {
                                        minContactDistances[MolContactInfo.PICNH] = piDist;
                                        contactAtomNumInResidueB[MolContactInfo.PICNH] = PiEffectCalculations.giveAtomNumOfNBondToArgH(argH); //sidechain N bond to argH
                                        contactAtomNumInResidueA[MolContactInfo.PICNH] = 5; //CG of five_ring
                                    }
//                                    System.out.println("New PICNH: " + argN.toString() + "/" + argH.toString());
                                    atomAtomContactType.add("PICNH");
                                    donorAcceptor[0] = atoms_b.get(PiEffectCalculations.giveAtomNumOfNBondToArgH(argH));
                                    donorAcceptor[1] = atoms_a.get(5);
                                    atomAtomContacts.add(donorAcceptor);
                                }
                            }
                        }
                    }
                }
                else {
                    if (! (atoms_b.size() >= 10)) {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms."
                            + " Continue search in next residues.");
                        }
                    else {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains no hydrogens."
                            + " Continue search in next residues.");
                    }
                }  
            }
            
            //PICAH            
            if ((atoms_b.size() > 0) && b.getHydrogenAtoms().size() >  0) {
                
                Atom ca = null;
                for (Atom c : atoms_b) {
                    if (c.getAtomName().contains("CA")) {
                        ca = c;
                    }
                }
                if (ca == null) {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                        + " Continue search in next residues.");
                }
                else {
                for (Atom ha : b.getHydrogenAtoms()) {
                    if (ha.getAtomName().contains("HA")) {                 
                    
                    
                    piDist = (int)(calculateDistancePiEffect(ca, ha, six_ring) / 10);
                    if ( piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PICAH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PICAH] < 0 || piDist > minContactDistances[MolContactInfo.PICAH]) {
                            minContactDistances[MolContactInfo.PICAH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PICAH] = 1; //CA
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PICAH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PICAH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PICAH: " + ca.toString() + "/" + ha.toString());
                        atomAtomContactType.add("PICAH");
                        donorAcceptor[0] = atoms_b.get(1);
                        atomAtomContacts.add(donorAcceptor);
                    }
                    
                        piDist = (int)(calculateDistancePiEffect(ca, ha, six_ring, false, true) / 10);
                    
                    if ( piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PICAH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PICAH] < 0 || piDist > minContactDistances[MolContactInfo.CAHPI]) {
                            minContactDistances[MolContactInfo.PICAH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PICAH] = 1; //CA
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PICAH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PICAH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PICAH: " + ca.toString() + "/" + ha.toString());
                        atomAtomContactType.add("PICAH");
                        donorAcceptor[0] = atoms_b.get(1);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    if ("TRP".equals(a.getName3())) {
                        piDist = (int)(calculateDistancePiEffect(ca, ha, five_ring) / 10);
                        if ( piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PICAH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PICAH] < 0 || piDist > minContactDistances[MolContactInfo.PICAH]) {
                                minContactDistances[MolContactInfo.PICAH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PICAH] = 1; //CA
                                if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                    contactAtomNumInResidueA[MolContactInfo.PICAH] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_a.get(5);
                                }
                                else {
                                    contactAtomNumInResidueA[MolContactInfo.PICAH] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_a.get(7);
                                }
                            }
//                            System.out.println("New PICAH: " + ca.toString() + "/" + ha.toString());
                            atomAtomContactType.add("PICAH");
                            donorAcceptor[0] = atoms_b.get(1);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(ca, ha, five_ring, false, true) / 10);
                        if ( piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PICAH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PICAH] < 0 || piDist > minContactDistances[MolContactInfo.PICAH]) {
                                minContactDistances[MolContactInfo.PICAH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PICAH] = 1; //CA
                                if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                    contactAtomNumInResidueA[MolContactInfo.PICAH] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_a.get(5);
                                }
                                else {
                                    contactAtomNumInResidueA[MolContactInfo.PICAH] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_a.get(7);
                                }
                            }
//                            System.out.println("New PICAH: " + ca.toString() + "/" + ha.toString());
                            atomAtomContactType.add("PICAH");
                            donorAcceptor[0] = atoms_b.get(1);
                            atomAtomContacts.add(donorAcceptor);
                        }

                    }
                }
                }
            }
            }
            else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no CA)."
                        + " Continue search in next residues.");
            }
            
            //PIPROCDH
            if (b.getName3().equals("PRO")) {       
                if (atoms_b.size() >= 7 && b.getHydrogenAtoms().size() > 0) {
                    for (Atom hd : b.getHydrogenAtoms()) {
                        if (hd.getAtomName().contains("HD")) {
                        piDist = (int)(calculateDistancePiEffect(atoms_b.get(6), hd, six_ring) / 10);
                        if ( piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PIPROCDH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PIPROCDH] < 0 || piDist > minContactDistances[MolContactInfo.PIPROCDH]) {
                                minContactDistances[MolContactInfo.PIPROCDH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PIPROCDH] = 6; //CD
                                if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_a.get(5);
                                }
                                else {
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_a.get(7);
                                }
                            }
//                            System.out.println("New PIPROCDH: " + atoms_b.get(6).toString() + "/" + hd.toString());
                            atomAtomContactType.add("PIPROCDH");
                            donorAcceptor[0] = atoms_b.get(6);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(atoms_b.get(6), hd, six_ring, false, true) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PIPROCDH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PIPROCDH] < 0 || piDist > minContactDistances[MolContactInfo.PIPROCDH]) {
                                minContactDistances[MolContactInfo.PIPROCDH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PIPROCDH] = 6; //CD
                                if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 5; //CG of six_ring
                                    donorAcceptor[1] = atoms_a.get(5);
                                }
                                else {
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 7; //CD2 of six_ring
                                    donorAcceptor[1] = atoms_a.get(7);
                                }
                            }
//                            System.out.println("New PIPROCDH: " + atoms_b.get(6).toString() + "/" + hd.toString());
                            atomAtomContactType.add("PIPROCDH");
                            donorAcceptor[0] = atoms_b.get(6);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        if ("TRP".equals(a.getName3())) {
                            piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(6), hd, five_ring) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PIPROCDH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PIPROCDH] < 0 || piDist > minContactDistances[MolContactInfo.PIPROCDH]) {
                                    minContactDistances[MolContactInfo.PIPROCDH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PIPROCDH] = 6; //CD
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 5; //CG of five_ring
                                }
//                                System.out.println("New PIPROCDH: " + atoms_b.get(6).toString() + "/" + hd.toString());
                                atomAtomContactType.add("PIPROCDH");
                                donorAcceptor[0] = atoms_b.get(6);
                                donorAcceptor[1] = atoms_a.get(5);
                                atomAtomContacts.add(donorAcceptor);
                            }

                            piDist = (int)(calculateDistancePiEffect(b.getAtoms().get(6), hd, five_ring, false, true) / 10);
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.PIPROCDH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.PIPROCDH] < 0 || piDist > minContactDistances[MolContactInfo.PIPROCDH]) {
                                    minContactDistances[MolContactInfo.PIPROCDH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.PIPROCDH] = 6; //CD
                                    contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = 5; //CG of five_ring
                                }
//                                System.out.println("New PIPROCDH: " + atoms_b.get(6).toString() + "/" + hd.toString());
                                atomAtomContactType.add("PIPROCDH");
                                donorAcceptor[0] = atoms_b.get(6);
                                donorAcceptor[1] = atoms_a.get(5);
                                atomAtomContacts.add(donorAcceptor);
                            }

                        }
                    }
                }
                }
                else {
                    if (! (atoms_a.size() >= 7)) {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms."
                            + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogens."
                            + " Continue search in next residues.");
                    }
                }
            }
            
            //PISH
            if (b.getName3().equals("CYS")) {          
                if (atoms_b.size() >= 6 && b.getHydrogenAtoms().size() > 0) {
                    for (Atom hg : b.getHydrogenAtoms()) {
                        if (hg.getAtomName().contains("HG")) {
                            piDist = (int)(calculateDistancePiEffect(atoms_b.get(5), hg, six_ring) / 10);
                            
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PISH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PISH] < 0 || piDist > minContactDistances[MolContactInfo.PISH]) {
                            minContactDistances[MolContactInfo.PISH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PISH] = 5;
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PISH: " + atoms_b.get(5).toString() + "/" + hg.toString());
                        atomAtomContactType.add("PISH");
                        donorAcceptor[0] = atoms_b.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    piDist = (int)(calculateDistancePiEffect(atoms_b.get(5), hg, six_ring, false, true) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PISH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PISH] < 0 || piDist > minContactDistances[MolContactInfo.PISH]) {
                            minContactDistances[MolContactInfo.PISH] = piDist;
                            contactAtomNumInResidueB[MolContactInfo.PISH] = 5;
                            if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 5; //CG of six_ring
                                donorAcceptor[1] = atoms_a.get(5);
                            }
                            else {
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 7; //CD2 of six_ring
                                donorAcceptor[1] = atoms_a.get(7);
                            }
                        }
//                        System.out.println("New PISH: " + atoms_b.get(5).toString() + "/" + hg.toString());
                        atomAtomContactType.add("PISH");
                        donorAcceptor[0] = atoms_b.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }
                    
                    if ("TRP".equals(a.getName3())) {
                        piDist = (int)(calculateDistancePiEffect(atoms_b.get(5), hg, five_ring) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PISH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PISH] < 0 || piDist > minContactDistances[MolContactInfo.PISH]) {
                                minContactDistances[MolContactInfo.PISH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PISH] = 5;
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 5; //CG of five_ring
                            }
//                            System.out.println("New PISH: " + atoms_b.get(5).toString() + "/" + hg.toString());
                            atomAtomContactType.add("PISH");
                            donorAcceptor[0] = atoms_b.get(5);
                            donorAcceptor[1] = atoms_a.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }

                        piDist = (int)(calculateDistancePiEffect(atoms_b.get(5), hg, five_ring, false, true) / 10);
                        if (piDist > 0) {
                            numPairContacts[MolContactInfo.TT]++;
                            numPairContacts[MolContactInfo.PISH]++;
                            donorAcceptor = new Atom[2];
                            if (minContactDistances[MolContactInfo.PISH] < 0 || piDist > minContactDistances[MolContactInfo.PISH]) {
                                minContactDistances[MolContactInfo.PISH] = piDist;
                                contactAtomNumInResidueB[MolContactInfo.PISH] = 5;
                                contactAtomNumInResidueA[MolContactInfo.PISH] = 5; //CG of five_ring
                            }
//                            System.out.println("New PISH: " + atoms_b.get(5).toString() + "/" + hg.toString());
                            atomAtomContactType.add("PISH");
                            donorAcceptor[0] = atoms_b.get(5);
                            donorAcceptor[1] = atoms_a.get(5);
                            atomAtomContacts.add(donorAcceptor);
                        }
                    }
                        
                      }
                    }
                }
                else {
                    if (! (atoms_a.size() >= 5)) {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no SG)."
                        + " Continue search in next residues.");
                    }
                    else {
                        DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogen (no HG)."
                            + " Continue search in next residues.");
                    }
                }
            }
            
            //PIXOH
            if (sidechainOHAAs.contains(b.getName3())) {
                Atom OHAA_X = null;
                Atom OHAA_H = null;
                
                if ("SER".equals(b.getName3())) {
                    if (atoms_b.size() >= 6 && b.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_b.get(5);
                        for (Atom hg : b.getHydrogenAtoms()) {
                            if (hg.getAtomName().contains("HG")) {
                                OHAA_H = hg;
                            }
                        }
                    }
                    else {
                        if (! (atoms_b.size() >= 6)) {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no OG).");
                        }
                        else {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogens (no HG).");
                        }
                    }
                }
                else if ("THR".equals(b.getName3())) {
                    if (atoms_b.size() >= 6 && b.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_b.get(5);
                        for (Atom hg : b.getHydrogenAtoms()) {
                            if (hg.getAtomName().contains("HG")) {
                                OHAA_H = hg;
                            }
                        }
                    }
                    else {
                        if (! (atoms_b.size() >= 6)) {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no OG1).");
                        }
                        else {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogens (no HG1).");
                        }
                    }
                }
                else {
                    if (atoms_b.size() >= 12 && b.getHydrogenAtoms().size() > 0) {
                        OHAA_X = atoms_b.get(11);
                        for (Atom hh : b.getHydrogenAtoms()) {
                            if (hh.getAtomName().contains("HH")) {
                                OHAA_H = hh;
                            }
                        }
                    }
                    else {
                        if (! (atoms_b.size() >= 12)) {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no OH).");
                        }
                        else {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogens (no HH).");
                        }
                    }
                }
                
                piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, six_ring) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.PIXOH]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.PIXOH] < 0 || piDist > minContactDistances[MolContactInfo.PIXOH]) {
                        minContactDistances[MolContactInfo.PIXOH] = piDist;
                        if ("TYR".equals(b.getName3())) {
                            contactAtomNumInResidueB[MolContactInfo.PIXOH] = 11;
                            donorAcceptor[0] = atoms_b.get(11);
                        }
                        else {
                            contactAtomNumInResidueB[MolContactInfo.PIXOH] = 5;
                            donorAcceptor[0] = atoms_b.get(5);
                        }
                        if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 5; //CG of six_ring
                            donorAcceptor[1] = atoms_a.get(5);
                        }
                        else {
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 7; //CD2 of six_ring
                            donorAcceptor[1] = atoms_a.get(7);
                        }
                    }
//                    System.out.println("New PIXOH: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                    atomAtomContactType.add("PIXOH");
                    atomAtomContacts.add(donorAcceptor);
                }

                piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, six_ring, false, true) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.PIXOH]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.PIXOH] < 0 || piDist > minContactDistances[MolContactInfo.PIXOH]) {
                        minContactDistances[MolContactInfo.PIXOH] = piDist;
                        if ("TYR".equals(b.getName3())) {
                            contactAtomNumInResidueB[MolContactInfo.PIXOH] = 11;
                            donorAcceptor[0] = atoms_b.get(11);
                        }
                        else {
                            contactAtomNumInResidueB[MolContactInfo.PIXOH] = 5;
                            donorAcceptor[0] = atoms_b.get(5);
                        }
                        if ("TYR".equals(a.getName3()) || "PHE".equals(a.getName3())) {
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 5; //CG of six_ring
                            donorAcceptor[1] = atoms_a.get(5);
                        }
                        else {
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 7; //CD2 of six_ring
                            donorAcceptor[1] = atoms_a.get(7);
                        }
                    }
//                    System.out.println("New PIXOH: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                    atomAtomContactType.add("PIXOH");
                    atomAtomContacts.add(donorAcceptor);
                }

                if ("TRP".equals(a.getName3())) {
                    piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, five_ring) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PIXOH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PIXOH] < 0 || piDist > minContactDistances[MolContactInfo.PIXOH]) {
                            minContactDistances[MolContactInfo.PIXOH] = piDist;
                            if ("TYR".equals(b.getName3())) {
                                contactAtomNumInResidueB[MolContactInfo.PIXOH] = 11;
                                donorAcceptor[0] = atoms_b.get(11);
                            }
                            else {
                                contactAtomNumInResidueB[MolContactInfo.PIXOH] = 5;
                                donorAcceptor[0] = atoms_b.get(5);
                            }
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 5; //CG of five_ring
                        }
//                        System.out.println("New PIXOH: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                        atomAtomContactType.add("PIXOH");
                        donorAcceptor[1] = atoms_a.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }

                    piDist = (int)(calculateDistancePiEffect(OHAA_X, OHAA_H, five_ring, false, true) / 10);
                    if (piDist > 0) {
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.PIXOH]++;
                        donorAcceptor = new Atom[2];
                        if (minContactDistances[MolContactInfo.PIXOH] < 0 || piDist > minContactDistances[MolContactInfo.PIXOH]) {
                            minContactDistances[MolContactInfo.PIXOH] = piDist;
                            if ("TYR".equals(b.getName3())) {
                                contactAtomNumInResidueB[MolContactInfo.PIXOH] = 11;
                                donorAcceptor[0] = atoms_b.get(11);
                            }
                            else {
                                contactAtomNumInResidueB[MolContactInfo.PIXOH] = 5;
                                donorAcceptor[0] = atoms_b.get(5);
                            }
                            contactAtomNumInResidueA[MolContactInfo.PIXOH] = 5; //CG of five_ring
                        }
//                        System.out.println("New PIXOH: " + OHAA_X.toString() + "/" + OHAA_H.toString());
                        atomAtomContactType.add("PIXOH");
                        donorAcceptor[1] = atoms_a.get(5);
                        atomAtomContacts.add(donorAcceptor);
                    }
                }   
            }
        }
        
        //C=O in sidechain of aminoacid a acts as acceptor
        //CCOCAH
        if (sidechainOAAs.contains(a.getName3())) {
            //check for a valid res a
            if ( (("GLU".equals(a.getName3()) || "GLN".equals(a.getName3())) && atoms_a.size() >= 8) ||
                    (("ASP".equals(a.getName3()) || "ASN".equals(a.getName3())) && atoms_a.size() >= 7) ) {
                if ( ("GLU".equals(a.getName3()) || "GLN".equals(a.getName3()) && atoms_a.get(6).getAtomShortName().equals("CD") && 
                        atoms_a.get(7).getAtomShortName().equals("OE1")) ||
                       (("ASP".equals(a.getName3()) || "ASN".equals(a.getName3())) && atoms_a.get(5).getAtomShortName().equals("CG") &&
                        atoms_a.get(6).getAtomShortName().equals("OD1")) ) {  
            
                    if ((atoms_b.size() >= 2) && b.getHydrogenAtoms().size() > 0) {
                        Atom ca = null;
                        for (Atom c : atoms_b) {
                            if (c.getAtomName().contains("CA")) {
                                ca = c;
                            }
                        }
                        if (ca == null) {
                            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                                + " Continue search in next residues.");
                        }
                        else {

                        for (Atom ha : b.getHydrogenAtoms()) {
                            if (ha.getAtomName().contains("HA")) {                
                            
                                
                                if ("GLU".equals(a.getName3()) || "GLN".equals(a.getName3())) {
                                    piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_a.get(6), atoms_a.get(7)) / 10);
                                }
                                else {
                                    piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_a.get(5), atoms_a.get(6)) / 10);
                                }
                            if (piDist > 0) {
                                numPairContacts[MolContactInfo.TT]++;
                                numPairContacts[MolContactInfo.CCOCAH]++;
                                donorAcceptor = new Atom[2];
                                if (minContactDistances[MolContactInfo.CCOCAH] < 0 || piDist > minContactDistances[MolContactInfo.CCOCAH]) {
                                    minContactDistances[MolContactInfo.CCOCAH] = piDist;
                                    contactAtomNumInResidueB[MolContactInfo.CCOCAH] = 1; //CA
                                    if ("GLU".equals(a.getName3()) || "GLN".equals(a.getName3())) {
                                        contactAtomNumInResidueA[MolContactInfo.CCOCAH] = 7; //OE1
                                        donorAcceptor[1] = atoms_a.get(7);
                                    }
                                    else {
                                        contactAtomNumInResidueA[MolContactInfo.CCOCAH] = 6; //OD1
                                        donorAcceptor[1] = atoms_a.get(6);
                                    }
                                }
//                                System.out.println("New CCOCAH: " + ca.toString() + "/" + ha.toString());
                                atomAtomContactType.add("CCOCAH");
                                donorAcceptor[0] = atoms_b.get(1);
                                atomAtomContacts.add(donorAcceptor);
                            }
                            
                        }
                        }
                    }
                    }
                    else {
                            DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no CA)."
                                + " Continue search in next residues.");
                    }
            
                }
                else {
                    DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not the expected atom names for C=O."
                        + " Continue seach in next residues.");
                }
            }
            else {
                DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no C=O)."
                    + " Continue search in next residues.");
            }
        }
        
        //C=O in backbone of aminoacid a acts as acceptor
        //BCOCAH
        if (atoms_a.size() >= 4 && atoms_b.size() >  0 && b.getHydrogenAtoms().size() > 0) {
            if (atoms_a.get(2).getAtomShortName().equals("C") && atoms_a.get(3).getAtomShortName().equals("O")) {
            
            Atom ca = null;
            for (Atom c : atoms_b) {
                if (c.getAtomName().contains("CA")) {
                    ca = c;
                }
            }
            
            if (ca == null) {
                DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no CA)."
                    + " Continue search in next residues.");
                }
            else {
            for (Atom ha : b.getHydrogenAtoms()) {
                if (ha.getAtomName().contains("HA")) {
                
                    piDist = (int)(calculateDistancePiEffect(ca, ha, atoms_a.get(2), atoms_a.get(3)) / 10);
                if (piDist > 0) {
                    numPairContacts[MolContactInfo.TT]++;
                    numPairContacts[MolContactInfo.BCOCAH]++;
                    donorAcceptor = new Atom[2];
                    if (minContactDistances[MolContactInfo.BCOCAH] < 0 || piDist > minContactDistances[MolContactInfo.BCOCAH]) {
                        minContactDistances[MolContactInfo.BCOCAH] = piDist;
                        contactAtomNumInResidueB[MolContactInfo.BCOCAH] = 1; //CA
                        contactAtomNumInResidueA[MolContactInfo.BCOCAH] = 3; //O         
                    }
//                    System.out.println("New BCOCAH: " + ca.toString() + "/" + ha.toString());
                    atomAtomContactType.add("BCICAH");
                    donorAcceptor[0] = atoms_b.get(1);
                    donorAcceptor[1] = atoms_a.get(3);
                    atomAtomContacts.add(donorAcceptor);
                }
                
            }
        }
            }
        }
        else {
            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no C=O)."
                + " Continue search in next residues.");
        }
        }
        else {
            if (! (atoms_a.size() >= 4)) {
            DP.getInstance().w("main", a.getName3() + " (" + a.getFancyName() + " chain: " + a.getChainID() +  ") contains not enough atoms (no C=O)."
                + " Continue search in next residues.");
            }
            else {
                if (!(atoms_b.size() > 0)) {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough atoms (no CA)."
                + " Continue search in next residues.");
                }
                else {
                    DP.getInstance().w("main", b.getName3() + " (" + b.getFancyName() + " chain: " + b.getChainID() +  ") contains not enough hydrogen (no HA)."
                + " Continue search in next residues.");
                }
            }
        }
        }
        
         // Iteration through all atoms of the two residues is done
        if(numPairContacts[MolContactInfo.TT] > 0) {
            result = new MolContactInfo(numPairContacts, minContactDistances, contactAtomNumInResidueA, contactAtomNumInResidueB, a, b, CAdist, numTotalLigContactsPair);
        }
        else {
            result = null;
        }
        return(result);         // This is null if no contact was detected
    }
    
    
    /**
     * Alternative model to calculate the atom contacts between residue 'a' and 'b'.
     * This alternative model is used to calculate interchain contacts between
     * atoms that are used to detect protein-protein interactions.
     * @param a one of the residues of the residue pair
     * @param b one of the residues of the residue pair
     * @return A MolContactInfo object with information on the atom contacts between 'a' and 'b'.
     */
    public static MolContactInfo calculateAtomContactsBetweenResiduesAlternativeModel(Residue a, Residue b) {

        ArrayList<Atom> atoms_a = a.getAtoms();
        ArrayList<Atom> atoms_b = b.getAtoms();
        
        Atom x, y;
        Integer dist = null;
        Integer CAdist = a.distTo(b);
        MolContactInfo result = null;
        
        ArrayList<Atom[]> atomAtomContacts = new ArrayList<Atom[]>();
        Atom[] donorAcceptor = new Atom[2];
        ArrayList<String> atomAtomContactType = new ArrayList<String>();
        ArrayList<String> sidechainOAAs = new ArrayList<String>();
        sidechainOAAs.add("ASP");
        sidechainOAAs.add("GLU");
        sidechainOAAs.add("ASN");
        sidechainOAAs.add("GLN");
        ArrayList<String> sidechainOHAAs = new ArrayList<String>();
        sidechainOHAAs.add("SER");
        sidechainOHAAs.add("THR");
        sidechainOHAAs.add("TYR");
        ArrayList<String> sidechainNHAAs = new ArrayList<String>();
        sidechainNHAAs.add("ARG");
        sidechainNHAAs.add("HIS");
        sidechainNHAAs.add("LYS");
        sidechainNHAAs.add("ASN");
        sidechainNHAAs.add("GLN");
        sidechainNHAAs.add("TRP");
        ArrayList<String> sidechainPiRings = new ArrayList<String>();
        sidechainPiRings.add("TRP");
        sidechainPiRings.add("TYR");
        sidechainPiRings.add("PHE");
        
        // Only used for testing purposes right now.
        // System.out.println(checkAromaticRingPlanarity(a)); 
        // System.out.println(checkAromaticRingPlanarity(b)); 
       
        Integer[] numPairContacts = new Integer[Main.NUM_RESIDUE_PAIR_CONTACT_TYPES_ALTERNATIVE_MODEL];
        // The positions in the numPairContacts array hold the number of contacts of each type for a pair of residues:
        // Some cheap vars to make things easier to understand (a replacement for #define):
        /*
        Integer TT = 0;         //  0 = total number of contacts            (all residue type combinations)
        Integer BB = 1;         //  1 = # of backbone-backbone contacts     (protein - protein only)
        Integer CB = 2;         //  2 = # of sidechain-backbone contacts    (protein - protein only)
        Integer BC = 3;         //  3 = # of backbone-sidechain contacts    (protein - protein only)
        Integer CC = 4;         //  4 = # of sidechain-sidechain contacts   (protein - protein only)
        Integer HB = 5;         //  5 = # of H-bridge contacts 1, N=>0      (protein - protein only)
        Integer BH = 6;         //  6 = # of H-bridge contacts 2, 0=>N      (protein - protein only)
        Integer BL = 7;         //  7 = # of backbone-ligand contacts       (protein - ligand only)
        Integer LB = 8;         //  8 = # of ligand-backbone contacts       (protein - ligand only)
        Integer CL = 9;         //  9 = # of sidechain-ligand contacts      (protein - ligand only)
        Integer LC = 10;        // 10 = # of ligand-sidechain contacts      (protein - ligand only)
        Integer LL = 11;        // 11 = # of ligand-ligand contacts         (ligand - ligand only)
        Integer DISULFIDE = 12  // 12 = # of disulfide bridges
        ---------------- alternative model contact types ----------------------
        Integer IHB = 13        // 13 = # of interchain H-bridge contacts 1, N=>O
        Integer IBH = 14        // 14 = # of interchain H-bridge contacts 2, O=>N
        Integer IVDW = 15       // 15 = # of interchain van der Waals interactions
        Integer ISS = 16        // 16 = # of interchain disulfide bridges
        Integer IPI = 17        // 17 = # of interchain pi-effects
        Integer ISB = 18        // 18 = # of interchain salt bridges
        */

        
        Integer numTotalLigContactsPair = 0;



        Integer[] minContactDistances = new Integer[numPairContacts.length];
        // Holds the minimal distances of contacts of the appropriate type (see numPairContacts, index 0 is unused)

        Integer[] contactAtomNumInResidueA = new Integer[numPairContacts.length];
        // Holds the number Atom x has in its residue a for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)

        Integer[] contactAtomNumInResidueB = new Integer[numPairContacts.length];
        // Holds the number Atom y has in its residue b for the contact with minimal distance of that type.
        // See minContactDistances and numPairContacts; index 0 is unused; index 5 + 6 are also unused (atom is obvious + always the same)
        

        for(Integer k = 0; k < numPairContacts.length; k++) {      // init all arrays
            numPairContacts[k] = 0;
            minContactDistances[k] = -1;    // We are looking for the smallest distance >= 0 in this function so do NOT set '0' as the initial value or everything will get fucked up!
            contactAtomNumInResidueA[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            contactAtomNumInResidueB[k] = -1;   // We HAVE to assign '-1', NOT any other value here! See comments below for details.
            // The initial value of '-1' is required for the atom index arrays because our index
            // starts at '0', but geom_neo treads a '0' in a line of <pdbid>.geo as 'no contact' because
            // it starts its index at '1'.
            // This problem is solved by the functions in MolContactInfo: they return (our_index + 1). This
            // means that:
            //   1) If no contact was detected, our_index is -1 and they return 0, which means 'no contact' to geom_neo.
            //   2) If a contact was detected, our_index is converted to the geom_neo index. :)
        }

        // We assume that the first 5 atoms (index 0..4) in a residue that is an AA are backbone atoms,
        //  while all other (6..END) are assumed to be side chainName atoms.
        //  The backbone atoms should have atom names ' N  ', ' CA ', ' C  ' and ' O  ' but we don't check
        //  this atm because geom_neo doesn't do that and we want to stay compatible.
        //  Of course, all of this only makes sense for resides that are AAs, not for ligands. We care for that.
        Integer numOfLastBackboneAtomInResidue = 4;
        Integer atomIndexOfBackboneN = 0;       // backbone nitrogen atom index
        Integer atomIndexOfBackboneO = 3;       // backbone oxygen atom index
        Integer atomIndexOfBackboneC = 2;       // backbone carbon atom index
        Integer atomIndexOfBackboneCa = 1;    // backbone C-alpha atom index

        Integer aIntID = a.getInternalAAID();     // Internal AA ID (ALA=1, ARG=2, ...)
        Integer bIntID = b.getInternalAAID();
        Integer statAtomIDi, statAtomIDj;

        
        
        // Add contact type interchain disulfide bridge.
        // Interchain disulfide bridges are parsed from the the DSSP file itself earlier.
        // This adds the detected disulfide bridge contacts to the statistics.
        HashMap<Character, ArrayList<Integer>> interchainSB = FileParser.getInterchainSulfurBridges();
        
        // Check for the current residue pair if there is a sulfur bridge between them.
        // If that's the case, then add this contact to the statistics.
        for(ArrayList<Integer> valuePair : interchainSB.values()) {
                if((a.getDsspNum().equals(valuePair.get(0)) && b.getDsspNum().equals(valuePair.get(1))) || (a.getDsspNum().equals(valuePair.get(1)) && b.getDsspNum().equals(valuePair.get(0)))) {
                    numPairContacts[MolContactInfo.TT]++;   // update total number of contacts for this residue pair
                    numPairContacts[MolContactInfo.ISS]++;   // update disulfide bridge number of contacts for this residue pair
                }
        }
        

        // Iterate through all atoms of the two residues and check contacts for all pairs
        outerloop:
        for(Integer i = 0; i < atoms_a.size(); i++) {
            
            
            if(i >= MAX_ATOMS_PER_AA && a.isAA()) {
                DP.getInstance().w("calculateAtomContactsBetweenResidues(): The AA residue " + a.getUniquePDBName() + " of type " + a.getName3() + " has more atoms than allowed, skipping atom #" + i + ".");
                break;
            }
            
            x = atoms_a.get(i);

            innerloop:
            for(Integer j = 0; j < atoms_b.size(); j++) {
                                
                if(j >= MAX_ATOMS_PER_AA && b.isAA()) {
                    DP.getInstance().w("calculateAtomContactsBetweenResidues(): The AA residue " + b.getUniquePDBName() + " of type " + b.getName3() + " has more atoms than allowed, skipping atom #" + j + ".");
                    //continue;
                    break outerloop;
                }
                
                y = atoms_b.get(j);
                                
                
                // Check whether a contact exist. If so, classify it. Note that the code of geom_neo works based on the
                //  position of an atom in the atom list of its residue (e.g., it assumes that the 2nd atom of an AA is
                //  the C alpha atom. While this seems to hold for many PDB files it will produce wrong results if atoms
                //  are missing from the PDB file or other weird stuff is going on in there.

                //System.out.println("      Checking atom pair " + x.getChemSym() + " and " + y.getChemSym() + " with residue index " + i + " and " + j + " of residues " + a.getFancyName() + " and " + b.getFancyName() + ".");
                //System.out.println("        " + x);
                //System.out.println("        " + y);

                dist = x.distToAtom(y);
                
                // H-bonds
                if (dist < 39) {
                    // Backbone - Backbone H-bonds
                    // NH -> 0
                    if (i.equals(atomIndexOfBackboneN) && j.equals(atomIndexOfBackboneO)) {
                        if (x.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                            for (Atom h : a.getHydrogenAtoms()) {
                                if (h.getAtomName().equals(" H  ")) {
                                    if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && h.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.BBHB]++;
                                        if ((minContactDistances[MolContactInfo.BBHB] < 0) || dist < minContactDistances[MolContactInfo.BBHB]) {
                                            minContactDistances[MolContactInfo.BBHB] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.BBHB] = i;
                                            contactAtomNumInResidueB[MolContactInfo.BBHB] = j;
                                        }
//                                        System.out.println("New BB NHO: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("BBNHO");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = x;
                                        donorAcceptor[1] = y;
                                        atomAtomContacts.add(donorAcceptor);
                                          
                                    }
                                }
                            }
                        }
                    }
                    if (i.equals(atomIndexOfBackboneO) && j.equals(atomIndexOfBackboneN)) {
                        if (y.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                            for (Atom h : b.getHydrogenAtoms()) {
                                if (h.getAtomName().equals(" H  ")) {
                                    if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && h.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.BBBH]++;
                                        if ((minContactDistances[MolContactInfo.BBBH] < 0) || dist < minContactDistances[MolContactInfo.BBBH]) {
                                            minContactDistances[MolContactInfo.BBBH] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.BBBH] = i;
                                            contactAtomNumInResidueB[MolContactInfo.BBBH] = j;
                                        }
//                                        System.out.println("New BB ONH: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("BBONH");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = y;
                                        donorAcceptor[1] = x;
                                        atomAtomContacts.add(donorAcceptor);
                                    }
                                }
                            }
                        }
                    }
                    
                                
                    // Backbone - Sidechain H-bonds
                    // O -> NH
                    if (i.equals(atomIndexOfBackboneO) && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" N"))) {
                        if (sidechainNHAAs.contains(b.getName3()) && y.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                            for (Atom h : b.getHydrogenAtoms()) {
                                if (b.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                        || b.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                        || b.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                        || b.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                        || b.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                        || b.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                    if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && h.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.BCBH]++;
                                        if (minContactDistances[MolContactInfo.BCBH] < 0 || dist < minContactDistances[MolContactInfo.BCBH]) {
                                            minContactDistances[MolContactInfo.BCBH] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.BCBH] = i;
                                            contactAtomNumInResidueB[MolContactInfo.BCBH] = j;
                                        }
//                                        System.out.println("New BS ONH: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("BSONH");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = y;
                                        donorAcceptor[1] = x;
                                        atomAtomContacts.add(donorAcceptor);
                                    }
                                }
                            }
                        }
                    }
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" N")) && j.equals(atomIndexOfBackboneO)) {
                        if (sidechainNHAAs.contains(a.getName3()) && x.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                            for (Atom h : a.getHydrogenAtoms()) {
                                if (a.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                        || a.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                        || a.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                        || a.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                        || a.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                        || a.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                    if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && h.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.CBHB]++;
                                        if ((minContactDistances[MolContactInfo.CBHB] < 0) || dist < minContactDistances[MolContactInfo.CBHB]) {
                                            minContactDistances[MolContactInfo.CBHB] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.CBHB] = i;
                                            contactAtomNumInResidueB[MolContactInfo.CBHB] = j;
                                        }
//                                        System.out.println("New SB ONH: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("SBONH");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = x;
                                        donorAcceptor[1] = y;
                                        atomAtomContacts.add(donorAcceptor);
                                    }
                                }
                            }
                        }
                    }
                    
                    // O -> OH
                    if (i.equals(atomIndexOfBackboneO) && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {
                        // Check if the amino acid is Serine, Threonine, or Tyrosine, as only those have hydroxy groups in their side chain.
                        // All other amino acids only have oxygen in the side chain as part of carbonyl groups.
                        if (sidechainOHAAs.contains(b.getName3()) && y.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                            for (Atom h : b.getHydrogenAtoms()) {
                                if (b.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                        || b.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                        || b.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                    if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && h.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.BCBH]++;
                                        if ((minContactDistances[MolContactInfo.BCBH] < 0) || dist < minContactDistances[MolContactInfo.BCBH]) {
                                            minContactDistances[MolContactInfo.BCBH] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.BCBH] = i;
                                            contactAtomNumInResidueB[MolContactInfo.BCBH] = j;
                                        }
//                                        System.out.println("New BS OOH: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("BSOOH");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = y;
                                        donorAcceptor[1] = x;
                                        atomAtomContacts.add(donorAcceptor);
                                    }
                                }
                            }
                        }
                    }
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O")) && j.equals(atomIndexOfBackboneO)) {
                        // Check if the amino acid is serine, threonine or tyrosine, as only those have hydroxy groups in their side chain.
                        // All other amino acids only have oxygen in the side chain as part of carbonyl groups.
                        if (sidechainOHAAs.contains(a.getName3()) && x.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                            for (Atom h : a.getHydrogenAtoms()) {
                                if (a.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                        || a.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                        || a.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                    if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && h.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneC))) {
                                        numPairContacts[MolContactInfo.TT]++;
                                        numPairContacts[MolContactInfo.CBHB]++;
                                        if ((minContactDistances[MolContactInfo.CBHB] < 0) || dist < minContactDistances[MolContactInfo.CBHB]) {
                                            minContactDistances[MolContactInfo.CBHB] = dist;
                                            contactAtomNumInResidueA[MolContactInfo.CBHB] = i;
                                            contactAtomNumInResidueB[MolContactInfo.CBHB] = j;
                                        }
//                                        System.out.println("New SB ONH: " + x.toString() + "/" + y.toString());
                                        atomAtomContactType.add("SBONH");
                                        donorAcceptor = new Atom[2];
                                        donorAcceptor[0] = x;
                                        donorAcceptor[1] = y;
                                        atomAtomContacts.add(donorAcceptor);
                                    }
                                }
                            }
                        }
                    }
                    
                    // NH -> O
                    if (i.equals(atomIndexOfBackboneN) && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {
                        // Check if amino acid is Asparagine, Glutamine, Glutamic Acid, or Aspartic Acid, as only those have carbonyl groups in their side chain.
                        try {
                            if ((b.getName3().equals("ASP") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                    || (b.getName3().equals("GLU") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                    || (b.getName3().equals("ASN") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                    || (b.getName3().equals("GLN") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))) {

                                for (Atom h : a.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(y, atoms_b.get(6))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.BCHB]++;
                                            if ((minContactDistances[MolContactInfo.BCHB] < 0) || dist < minContactDistances[MolContactInfo.BCHB]) {
                                                minContactDistances[MolContactInfo.BCHB] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.BCHB] = i;
                                                contactAtomNumInResidueB[MolContactInfo.BCHB] = j;
                                            }
//                                            System.out.println("New BS NHO: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("BSNHO");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = x;
                                            donorAcceptor[1] = y;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }

                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> O backbone-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }                    
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O")) && j.equals(atomIndexOfBackboneN)) {
                        // Check if amino acid is Asparagine, Glutamine, Glutamic Acid, or Aspartic Acid, as only those have carbonyl groups in their side chain.
                        try {
                            if ((a.getName3().equals("ASP") && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                    || (a.getName3().equals("GLU") && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                    || (a.getName3().equals("ASN") && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                    || (a.getName3().equals("GLN") && y.hbondAtomAngleBetween(x, atoms_a.get(6)))) {

                                for (Atom h : b.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(x, atoms_a.get(6))) {

                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.CBBH]++;
                                            if ((minContactDistances[MolContactInfo.CBBH] < 0) || dist < minContactDistances[MolContactInfo.CBBH]) {
                                                minContactDistances[MolContactInfo.CBBH] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.CBBH] = i;
                                                contactAtomNumInResidueB[MolContactInfo.CBBH] = j;
                                            }
//                                            System.out.println("New SB NHO: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("SBNHO");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = y;
                                            donorAcceptor[1] = x;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> O sidechain-backbone calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }
                    
                    // NH -> OH
                    if (i.equals(atomIndexOfBackboneN) && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {
                        // Check if amino acid is Serine, Threonine, or Tyrosine, as only those have hydroxy groups in their side chain.

                        // NH as donor
                        try {
                            if ((b.getName3().equals("SER") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue)))
                                    || (b.getName3().equals("THR") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue)))
                                    || (b.getName3().equals("TYR") && x.hbondAtomAngleBetween(y, atoms_b.get(10)))) {

                                for (Atom h : a.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TYR") && h.hbondAtomAngleBetween(y, atoms_b.get(10))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.BCHB]++;
                                            if ((minContactDistances[MolContactInfo.BCHB] < 0) || dist < minContactDistances[MolContactInfo.BCHB]) {
                                                minContactDistances[MolContactInfo.BCHB] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.BCHB] = i;
                                                contactAtomNumInResidueB[MolContactInfo.BCHB] = j;
                                            }
//                                            System.out.println("New BS NHOH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("BSNHOH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = x;
                                            donorAcceptor[1] = y;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            } // NH as acceptor
                            else if (sidechainOHAAs.contains(b.getName3()) && y.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneCa))) {
                                for (Atom h : b.getHydrogenAtoms()) {
                                    if (b.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                            || b.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                            || b.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                        if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && h.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneCa))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.BCBH]++;
                                            if ((minContactDistances[MolContactInfo.BCBH] < 0) || dist < minContactDistances[MolContactInfo.BCBH]) {
                                                minContactDistances[MolContactInfo.BCBH] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.BCBH] = i;
                                                contactAtomNumInResidueB[MolContactInfo.BCBH] = j;
                                            }
//                                            System.out.println("New BS NHOH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("BSNHOH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = y;
                                            donorAcceptor[1] = x;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> OH backbone-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }                    
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O")) && j.equals(atomIndexOfBackboneN)) {
                        // Check if amino acid is Serine, Threonine, or Tyrosine, as only those have hydroxy groups in their side chain.

                        // NH as donor
                        try {
                            if ((a.getName3().equals("SER") && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                    || (a.getName3().equals("THR") && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                    || (a.getName3().equals("TYR") && y.hbondAtomAngleBetween(x, atoms_a.get(10)))) {

                                for (Atom h : b.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TYR") && h.hbondAtomAngleBetween(x, atoms_a.get(10))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.CBBH]++;
                                            if ((minContactDistances[MolContactInfo.CBBH] < 0) || dist < minContactDistances[MolContactInfo.CBBH]) {
                                                minContactDistances[MolContactInfo.CBBH] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.CBBH] = i;
                                                contactAtomNumInResidueB[MolContactInfo.CBBH] = j;
                                            }
//                                            System.out.println("New SB NHOH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("SBNHOH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = y;
                                            donorAcceptor[1] = x;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            } // NH as acceptor
                            else if (sidechainOHAAs.contains(a.getName3()) && x.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneCa))) {
                                for (Atom h : a.getHydrogenAtoms()) {
                                    if (a.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                            || a.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                            || a.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                        if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && h.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneCa))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.CBHB]++;
                                            if ((minContactDistances[MolContactInfo.CBHB] < 0) || dist < minContactDistances[MolContactInfo.CBHB]) {
                                                minContactDistances[MolContactInfo.CBHB] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.CBHB] = i;
                                                contactAtomNumInResidueB[MolContactInfo.CBHB] = j;
                                            }
//                                            System.out.println("New SB NHOH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("SBNHOH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = x;
                                            donorAcceptor[1] = y;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }

                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> OH sidechain-backbone calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }
                    
                    // NH -> NH
                    if (i.equals(atomIndexOfBackboneN) && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" N"))) {
                        // Check if amino acid is Arginine, Histidine, Lysine, Asparagine, Glutamine, or Tryptophan, as only those have amine groups in their side chain.

                        // Backbone NH as donor
                        try {
                            if ((b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                    || (b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                    || (b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                    || (b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                    || (b.getName3().equals("LYS") && x.hbondAtomAngleBetween(y, atoms_b.get(7)))
                                    || (b.getName3().equals("ASN") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                    || (b.getName3().equals("GLN") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                    || (b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                    || (b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(9)))) {

                                for (Atom h : a.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(6))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("LYS") && h.hbondAtomAngleBetween(y, atoms_b.get(7))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(y, atoms_b.get(8))
                                                || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(9))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.BCHB]++;
                                            if ((minContactDistances[MolContactInfo.BCHB] < 0) || dist < minContactDistances[MolContactInfo.BCHB]) {
                                                minContactDistances[MolContactInfo.BCHB] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.BCHB] = i;
                                                contactAtomNumInResidueB[MolContactInfo.BCHB] = j;
                                            }
//                                            System.out.println("New BS NHNH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("BSNHNH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = x;
                                            donorAcceptor[1] = y;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        // Backbone NH as acceptor
                            else if (sidechainNHAAs.contains(b.getName3()) && y.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneCa))) {
                                for (Atom h : b.getHydrogenAtoms()) {
                                    if (b.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                            || b.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                            || b.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                            || b.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                            || b.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                            || b.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                        if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && h.hbondAtomAngleBetween(x, atoms_a.get(atomIndexOfBackboneCa))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.BCBH]++;
                                            if ((minContactDistances[MolContactInfo.BCBH] < 0) || dist < minContactDistances[MolContactInfo.BCBH]) {
                                                minContactDistances[MolContactInfo.BCBH] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.BCBH] = i;
                                                contactAtomNumInResidueB[MolContactInfo.BCBH] = j;
                                            }
//                                            System.out.println("New BS NHNH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("BSNHNH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = y;
                                            donorAcceptor[1] = x;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> NH backbone-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }                    
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" N")) && j.equals(atomIndexOfBackboneN)) {
                        // Check if amino acid is Arginine, Histidine, Lysine, Asparagine, Glutamine, or Tryptophan, as only those have amine groups in their side chain.

                        // Backbone NH as donor
                        try {
                            if ((a.getName3().equals("ARG") && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                    || (a.getName3().equals("ARG") && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                    || (a.getName3().equals("HIS") && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                    || (a.getName3().equals("HIS") && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                    || (a.getName3().equals("LYS") && y.hbondAtomAngleBetween(x, atoms_a.get(7)))
                                    || (a.getName3().equals("ASN") && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                    || (a.getName3().equals("GLN") && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                    || (a.getName3().equals("TRP") && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                    || (a.getName3().equals("TRP") && y.hbondAtomAngleBetween(x, atoms_a.get(9)))) {

                                for (Atom h : b.getHydrogenAtoms()) {
                                    if (h.getAtomName().equals(" H  ")) {
                                        if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(6))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("LYS") && h.hbondAtomAngleBetween(x, atoms_a.get(7))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(x, atoms_a.get(8))
                                                || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(9))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.CBBH]++;
                                            if ((minContactDistances[MolContactInfo.CBBH] < 0) || dist < minContactDistances[MolContactInfo.CBBH]) {
                                                minContactDistances[MolContactInfo.CBBH] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.CBBH] = i;
                                                contactAtomNumInResidueB[MolContactInfo.CBBH] = j;
                                            }
//                                            System.out.println("New SB NHNH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("SBNHNH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = y;
                                            donorAcceptor[1] = x;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        // Backbone NH as acceptor
                            else if (sidechainNHAAs.contains(a.getName3()) && x.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneCa))) {
                                for (Atom h : a.getHydrogenAtoms()) {
                                    if (a.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                            || a.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                            || a.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                            || a.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                            || a.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                            || a.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                        if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && h.hbondAtomAngleBetween(y, atoms_b.get(atomIndexOfBackboneCa))) {
                                            numPairContacts[MolContactInfo.TT]++;
                                            numPairContacts[MolContactInfo.CBHB]++;
                                            if ((minContactDistances[MolContactInfo.CBHB] < 0) || dist < minContactDistances[MolContactInfo.CBHB]) {
                                                minContactDistances[MolContactInfo.CBHB] = dist;
                                                contactAtomNumInResidueA[MolContactInfo.CBHB] = i;
                                                contactAtomNumInResidueB[MolContactInfo.CBHB] = j;
                                            }
//                                            System.out.println("New SB NHNH: " + x.toString() + "/" + y.toString());
                                            atomAtomContactType.add("SBNHNH");
                                            donorAcceptor = new Atom[2];
                                            donorAcceptor[0] = x;
                                            donorAcceptor[1] = y;
                                            atomAtomContacts.add(donorAcceptor);
                                        }
                                    }
                                }
                            }
                        } catch (java.lang.IndexOutOfBoundsException e) {
                            DP.getInstance().w("main", "NH -> NH sidechain-backbone calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                        }
                    }
                    
                    // Sidechain-Sidechain H-bond
                    
                    // O -> NH
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" N"))) {

                        // O is always acceptor
                        for (String sidechainNHAA : sidechainNHAAs) {
                            try {
                                if ((a.getName3().equals("ASP") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                        || (a.getName3().equals("GLU") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                        || (a.getName3().equals("ASN") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                        || (a.getName3().equals("GLN") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))) {

                                    for (Atom h : b.getHydrogenAtoms()) {
                                        if (b.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                || b.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                || b.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                || b.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                || b.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                || b.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                            
                                            if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                    || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(x, atoms_a.get(6))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCBH]++;
                                                if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                    minContactDistances[MolContactInfo.CCBH] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                }
//                                                System.out.println("New SS ONH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSONH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = y;
                                                donorAcceptor[1] = x;
                                                atomAtomContacts.add(donorAcceptor);
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "O -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                    }                   
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" N"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {

                        // O is always acceptor
                        for (String sidechainNHAA : sidechainNHAAs) {
                            try {
                                if ((b.getName3().equals("ASP") && a.getName3().equals(sidechainNHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (b.getName3().equals("GLU") && a.getName3().equals(sidechainNHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (b.getName3().equals("ASN") && a.getName3().equals(sidechainNHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (b.getName3().equals("GLN") && a.getName3().equals(sidechainNHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(6)))) {

                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                || a.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                || a.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                || a.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                || a.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                || a.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {
                                            
                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(y, atoms_b.get(6))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS ONH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSONH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "O -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                    }
                            
                    // O -> OH
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {

                        // O is always acceptor
                        // O in residue a, OH in residue B
                        for (String sidechainOHAA : sidechainOHAAs) {
                            try {

                                if ((a.getName3().equals("ASP") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                        || (a.getName3().equals("GLU") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                        || (a.getName3().equals("ASN") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                        || (a.getName3().equals("GLN") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))) {

                                    for (Atom h : b.getHydrogenAtoms()) {
                                        if (b.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                || b.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                || b.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                            if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                    || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(x, atoms_a.get(6))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCBH]++;
                                                if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                    minContactDistances[MolContactInfo.CCBH] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                }
//                                                System.out.println("New SS OOH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSOOH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = y;
                                                donorAcceptor[1] = x;
                                                atomAtomContacts.add(donorAcceptor);
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "O -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                        // O in residue b, OH in residue A
                        for (String sidechainOHAA : sidechainOHAAs) {
                            try {
                                if ((b.getName3().equals("ASP") && a.getName3().equals(sidechainOHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (b.getName3().equals("GLU") && a.getName3().equals(sidechainOHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (b.getName3().equals("ASN") && a.getName3().equals(sidechainOHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (b.getName3().equals("GLN") && a.getName3().equals(sidechainOHAA) && x.hbondAtomAngleBetween(y, atoms_b.get(6)))) {

                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                || a.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                || a.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ASP|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("GLU|GLN") && h.hbondAtomAngleBetween(y, atoms_b.get(6))) {

                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS OOH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSOOH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "O -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                    }
                    
                    // OH -> NH
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" N"))) {

                        Boolean contactFound = false;

                        // OH as donor
                        for (String sidechainOHAA : sidechainOHAAs) {
                            try {
                                if ((a.getName3().equals(sidechainOHAA) && b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("LYS") && x.hbondAtomAngleBetween(y, atoms_b.get(7)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("ASN") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("GLN") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainOHAA) && b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(9)))) {

                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                || a.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                || a.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(6))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("LYS") && h.hbondAtomAngleBetween(y, atoms_b.get(7))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(y, atoms_b.get(8))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(9))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS OHNH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSOHNH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                                contactFound = true;
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "OH -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                        
                        if (!contactFound) {
                            // OH as acceptor
                            for (String sidechainNHAA : sidechainNHAAs) {
                                try {
                                    if ((a.getName3().equals("SER") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                            || (a.getName3().equals("THR") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                            || (a.getName3().equals("TYR") && b.getName3().equals(sidechainNHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(10)))) {

                                        for (Atom h : b.getHydrogenAtoms()) {
                                            if (b.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                    || b.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                    || b.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                    || b.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                    || b.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                    || b.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {

                                                if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TYR") && h.hbondAtomAngleBetween(x, atoms_a.get(10))) {
                                                    numPairContacts[MolContactInfo.TT]++;
                                                    numPairContacts[MolContactInfo.CCBH]++;
                                                    if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                        minContactDistances[MolContactInfo.CCBH] = dist;
                                                        contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                        contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                    }
//                                                    System.out.println("New SS OHNH: " + x.toString() + "/" + y.toString());
                                                    atomAtomContactType.add("SSOHNH");
                                                    donorAcceptor = new Atom[2];
                                                    donorAcceptor[0] = y;
                                                    donorAcceptor[1] = x;
                                                    atomAtomContacts.add(donorAcceptor);
                                                }
                                            }
                                        }
                                    }
                                } catch (java.lang.IndexOutOfBoundsException e) {
                                    DP.getInstance().w("main", "OH -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                                }
                            }
                        }
                    }
                    
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" N"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {

                        Boolean contactFound = false;

                        // NH as donor
                        for (String sidechainNHAA : sidechainNHAAs) {
                            try {
                                if ((a.getName3().equals(sidechainNHAA) && b.getName3().equals("SER") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("THR") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("TYR") && x.hbondAtomAngleBetween(y, atoms_b.get(10)))) {

                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                || a.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                || a.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                || a.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                || a.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                || a.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {

                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TYR") && h.hbondAtomAngleBetween(y, atoms_b.get(10))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS OHNH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSOHNH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                                contactFound = true;
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "NH -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                        
                        if (!contactFound) {
                            // NH as acceptor
                            for (String sidechainOHAA : sidechainOHAAs) {
                                try {
                                    if ((a.getName3().equals("ARG") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("ARG") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                            || (a.getName3().equals("HIS") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                            || (a.getName3().equals("HIS") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                            || (a.getName3().equals("LYS") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(7)))
                                            || (a.getName3().equals("ASN") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                            || (a.getName3().equals("GLN") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("TRP") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("TRP") && b.getName3().equals(sidechainOHAA) && y.hbondAtomAngleBetween(x, atoms_a.get(9)))) {

                                        for (Atom h : b.getHydrogenAtoms()) {
                                            if (b.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                    || b.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                    || b.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {
                                                if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(6))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("LYS") && h.hbondAtomAngleBetween(x, atoms_a.get(7))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(x, atoms_a.get(8))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(9))) {
                                                    numPairContacts[MolContactInfo.TT]++;
                                                    numPairContacts[MolContactInfo.CCBH]++;
                                                    if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                        minContactDistances[MolContactInfo.CCBH] = dist;
                                                        contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                        contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                    }
//                                                    System.out.println("New SS OHNH: " + x.toString() + "/" + y.toString());
                                                    atomAtomContactType.add("SSOHNH");
                                                    donorAcceptor = new Atom[2];
                                                    donorAcceptor[0] = y;
                                                    donorAcceptor[1] = x;
                                                    atomAtomContacts.add(donorAcceptor);
                                                }
                                            }
                                        }
                                    }
                                } catch (java.lang.IndexOutOfBoundsException e) {
                                    DP.getInstance().w("main", "NH -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                                }
                            }
                        }
                    }
                    
                    
                    // OH -> OH
                    if(((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" O")) 
                        && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" O"))) {
                        
                        Boolean contactFound = false;
                        
                        // first OH as donor
                        for(String sidechainOHAA : sidechainOHAAs) {
                        try {
                        if((a.getName3().equals(sidechainOHAA) && b.getName3().equals("SER") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue))) ||
                           (a.getName3().equals(sidechainOHAA) && b.getName3().equals("THR") && x.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue))) ||
                           (a.getName3().equals(sidechainOHAA) && b.getName3().equals("TYR") && x.hbondAtomAngleBetween(y, atoms_b.get(10)))    ) {
                            
                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                || a.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                || a.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {

                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(y, atoms_b.get(numOfLastBackboneAtomInResidue))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TYR") && h.hbondAtomAngleBetween(y, atoms_b.get(10))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS OHOH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSOHOH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                                contactFound = true;
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "OH -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                        
                        if (!contactFound) {
                            for (String sidechainOHAA2 : sidechainOHAAs) {
                                try {
                                    // first OH as acceptor
                                    if ((a.getName3().equals("SER") && b.getName3().equals(sidechainOHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                            || (a.getName3().equals("THR") && b.getName3().equals(sidechainOHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue)))
                                            || (a.getName3().equals("TYR") && b.getName3().equals(sidechainOHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(10)))) {

                                        for (Atom h : b.getHydrogenAtoms()) {
                                            if (b.getName3().equals("SER") && h.getAtomName().equals(" HG ")
                                                    || b.getName3().equals("THR") && h.getAtomName().equals(" HG1")
                                                    || b.getName3().equals("TYR") && h.getAtomName().matches(" HH ")) {

                                                if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("SER|THR") && h.hbondAtomAngleBetween(x, atoms_a.get(numOfLastBackboneAtomInResidue))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TYR") && h.hbondAtomAngleBetween(x, atoms_a.get(10))) {
                                                    numPairContacts[MolContactInfo.TT]++;
                                                    numPairContacts[MolContactInfo.CCBH]++;
                                                    if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                        minContactDistances[MolContactInfo.CCBH] = dist;
                                                        contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                        contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                    }
//                                                    System.out.println("New SS OHOH: " + x.toString() + "/" + y.toString());
                                                    atomAtomContactType.add("SSOHOH");
                                                    donorAcceptor = new Atom[2];
                                                    donorAcceptor[0] = y;
                                                    donorAcceptor[1] = x;
                                                    atomAtomContacts.add(donorAcceptor);
                                                }
                                            }
                                        }
                                    }
                                } catch (java.lang.IndexOutOfBoundsException e) {
                                    DP.getInstance().w("main", "OH -> OH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                                }
                            }
                        }
                    }
                    
                    
                    // NH -> NH
                    if (((i > numOfLastBackboneAtomInResidue) && x.getChemSym().equals(" N"))
                            && ((j > numOfLastBackboneAtomInResidue) && y.getChemSym().equals(" N"))) {

                        Boolean contactFound = false;

                        // first NH as donor
                        for (String sidechainNHAA : sidechainNHAAs) {
                            try {
                                if ((a.getName3().equals(sidechainNHAA) && b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("ARG") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("HIS") && x.hbondAtomAngleBetween(y, atoms_b.get(8)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("LYS") && x.hbondAtomAngleBetween(y, atoms_b.get(7)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("ASN") && x.hbondAtomAngleBetween(y, atoms_b.get(5)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("GLN") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(6)))
                                        || (a.getName3().equals(sidechainNHAA) && b.getName3().equals("TRP") && x.hbondAtomAngleBetween(y, atoms_b.get(9)))) {

                                    for (Atom h : a.getHydrogenAtoms()) {
                                        if (a.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                || a.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                || a.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                || a.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                || a.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                || a.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {

                                            if (h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(y, atoms_b.get(5))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(6))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("LYS") && h.hbondAtomAngleBetween(y, atoms_b.get(7))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(y, atoms_b.get(8))
                                                    || h.distToAtom(y) < 25 && x.hbondAtomAngleBetween(h, y) && b.getName3().equals("TRP") && h.hbondAtomAngleBetween(y, atoms_b.get(9))) {
                                                numPairContacts[MolContactInfo.TT]++;
                                                numPairContacts[MolContactInfo.CCHB]++;
                                                if ((minContactDistances[MolContactInfo.CCHB] < 0) || dist < minContactDistances[MolContactInfo.CCHB]) {
                                                    minContactDistances[MolContactInfo.CCHB] = dist;
                                                    contactAtomNumInResidueA[MolContactInfo.CCHB] = i;
                                                    contactAtomNumInResidueB[MolContactInfo.CCHB] = j;
                                                }
//                                                System.out.println("New SS NHNH: " + x.toString() + "/" + y.toString());
                                                atomAtomContactType.add("SSNHNH");
                                                donorAcceptor = new Atom[2];
                                                donorAcceptor[0] = x;
                                                donorAcceptor[1] = y;
                                                atomAtomContacts.add(donorAcceptor);
                                                contactFound = true;
                                            }
                                        }
                                    }
                                }
                            } catch (java.lang.IndexOutOfBoundsException e) {
                                DP.getInstance().w("main", "NH -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                            }
                        }
                         
                        if (!contactFound) {
                            for (String sidechainNHAA2 : sidechainNHAAs) {
                                try {
                                    // first NH as acceptor
                                    if ((a.getName3().equals("ARG") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("ARG") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                            || (a.getName3().equals("HIS") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                            || (a.getName3().equals("HIS") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(8)))
                                            || (a.getName3().equals("LYS") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(7)))
                                            || (a.getName3().equals("ASN") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(5)))
                                            || (a.getName3().equals("GLN") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("TRP") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(6)))
                                            || (a.getName3().equals("TRP") && b.getName3().equals(sidechainNHAA2) && y.hbondAtomAngleBetween(x, atoms_a.get(9)))) {

                                        for (Atom h : b.getHydrogenAtoms()) {
                                            if (b.getName3().equals("ARG") && h.getAtomName().equals(" HE ")
                                                    || b.getName3().equals("HIS") && h.getAtomName().equals(" HE2")
                                                    || b.getName3().equals("LYS") && h.getAtomName().matches(" HZ1| HZ2| HZ3")
                                                    || b.getName3().equals("ASN") && h.getAtomName().matches("HD21|HD22")
                                                    || b.getName3().equals("GLN") && h.getAtomName().matches("HE21|HE22")
                                                    || b.getName3().equals("TRP") && h.getAtomName().equals(" HE1")) {

                                                if (h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("HIS|ASN") && h.hbondAtomAngleBetween(x, atoms_a.get(5))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|GLN|TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(6))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("LYS") && h.hbondAtomAngleBetween(x, atoms_a.get(7))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().matches("ARG|HIS") && h.hbondAtomAngleBetween(x, atoms_a.get(8))
                                                        || h.distToAtom(x) < 25 && y.hbondAtomAngleBetween(h, x) && a.getName3().equals("TRP") && h.hbondAtomAngleBetween(x, atoms_a.get(9))) {
                                                    numPairContacts[MolContactInfo.TT]++;
                                                    numPairContacts[MolContactInfo.CCBH]++;
                                                    if ((minContactDistances[MolContactInfo.CCBH] < 0) || dist < minContactDistances[MolContactInfo.CCBH]) {
                                                        minContactDistances[MolContactInfo.CCBH] = dist;
                                                        contactAtomNumInResidueA[MolContactInfo.CCBH] = i;
                                                        contactAtomNumInResidueB[MolContactInfo.CCBH] = j;
                                                    }
//                                                    System.out.println("New SS NHNH: " + x.toString() + "/" + y.toString());
                                                    atomAtomContactType.add("SSNHNH");
                                                    donorAcceptor = new Atom[2];
                                                    donorAcceptor[0] = y;
                                                    donorAcceptor[1] = x;
                                                    atomAtomContacts.add(donorAcceptor);
                                                }
                                            }
                                        }
                                    }
                                } catch (java.lang.IndexOutOfBoundsException e) {
                                    DP.getInstance().w("main", "NH -> NH sidechain-sidechain calculation failed for residues " + x.getPdbResNum() + " and " + y.getPdbResNum() + ". Possibly because of missing atoms of the sidechain in the pdb file.");
                                }
                            }
                        }
                    }
                } // end of if(dist < 39)
                    
                    
                    
                    
                
                
                
                if(x.vdwAtomContactTo(y)) {             // If a contact is detected, Atom.vdwAtomContactTo() returns true

                    // The van der Waals radii spheres overlap, contact found.
                    // NOT updating total number of contacts for this residue pair at this point, since we need to check later if we want to include ligand contacts
                    
                    // DEBUG
                    //System.out.println("DEBUG: Atom contact in distance " + dist + " between atom " + x + " and " + y + ".");


                    // Update contact statistics.
                    statAtomIDi = i + 1;    // The field '0' is used for all contacts and we need to follow geom_neo conventions so we start the index at 1 instead of 0.
                    statAtomIDj = j + 1;
                    if(x.isLigandAtom()) { statAtomIDi = 1; }       // Different ligands can have different numbers of atoms and separating them just makes no sense. We assign all contacts to the first atom.
                    if(y.isLigandAtom()) { statAtomIDj = 1; }
                    
                    try {

                        contact[0][0][0][0]++;                 // update global total number of contacts
                        contact[aIntID][bIntID][0][0]++;       // contacts AA type a <-> AA type b
                        contact[bIntID][aIntID][0][0]++;       // contacts AA type b <-> AA type a


                        //System.out.println("DEBUG: a=" + aIntID + ",b=" + bIntID + ",i=" + statAtomIDi + ",j=" + statAtomIDj + ". Residues=" + a.getFancyName() + "," + b.getFancyName() + ".");
                        contact[aIntID][bIntID][statAtomIDi][statAtomIDj]++;       // contacts of atoms of AAs
                        contact[bIntID][aIntID][statAtomIDj][statAtomIDi]++;

                        contact[aIntID][bIntID][statAtomIDi][0]++;                 // total number of contacts for atom x of this AA
                        contact[bIntID][aIntID][0][statAtomIDi]++;

                        contact[aIntID][bIntID][statAtomIDj][0]++;                 // total number of contacts for atom x of this AA
                        contact[bIntID][aIntID][0][statAtomIDj]++;
                    } catch(java.lang.ArrayIndexOutOfBoundsException e) {
                        //DP.getInstance().w("calculateAtomContactsBetweenResidues():Contact statistics array out of bounds. Residues with excessive number of atoms detected: " + e.getMessage() + ".");
                        DP.getInstance().w("calculateAtomContactsBetweenResidues(): Atom count for residues too high (" + e.getMessage() + "), ignoring contacts for these atoms (aIntID=" + aIntID + ", bIntID=" + bIntID + ", statAtomIDi=" + statAtomIDi + ", statAtomIDj=" + statAtomIDj + ").");
                        continue;
                    }
                    
                    // Determine the contact type.                    
                    if(x.isProteinAtom() && y.isProteinAtom()) {
                        // *************************** protein - protein contact *************************
                        numPairContacts[MolContactInfo.TT]++;
                        numPairContacts[MolContactInfo.IVDW]++;
                        if(minContactDistances[MolContactInfo.IVDW] < 0 || dist < minContactDistances[MolContactInfo.IVDW]) {
                            minContactDistances[MolContactInfo.IVDW] = dist;
                            contactAtomNumInResidueA[MolContactInfo.IVDW] = i;
                            contactAtomNumInResidueB[MolContactInfo.IVDW] = j;
                        }
//                        System.out.println("New IVDW: " + x.toString() + "/" + y.toString());
                        atomAtomContactType.add("IVDW");
                        donorAcceptor = new Atom[2];
                        donorAcceptor[0] = x;
                        donorAcceptor[1] = y;
                        atomAtomContacts.add(donorAcceptor);


                        // Check the exact contact type
                        if(i <= numOfLastBackboneAtomInResidue && j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - backbone contact
                            numPairContacts[MolContactInfo.BB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BB] < 0) || dist < minContactDistances[MolContactInfo.BB]) {
                                minContactDistances[MolContactInfo.BB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BB] = i;
                                contactAtomNumInResidueB[MolContactInfo.BB] = j;
                            }

                        }
                        else if(i > numOfLastBackboneAtomInResidue && j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a chainName - backbone contact
                            numPairContacts[MolContactInfo.CB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CB] < 0) || dist < minContactDistances[MolContactInfo.CB]) {
                                minContactDistances[MolContactInfo.CB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CB] = i;
                                contactAtomNumInResidueB[MolContactInfo.CB] = j;
                            }

                        }
                        else if(i <= numOfLastBackboneAtomInResidue && j > numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - chainName contact
                            numPairContacts[MolContactInfo.BC]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BC] < 0) || dist < minContactDistances[MolContactInfo.BC]) {
                                minContactDistances[MolContactInfo.BC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BC] = i;
                                contactAtomNumInResidueB[MolContactInfo.BC] = j;
                            }
                        }
                        else if(i > numOfLastBackboneAtomInResidue && j > numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a chainName - chainName contact
                            numPairContacts[MolContactInfo.CC]++;          // 'C' instead of 'S' for side chainName pays off

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CC] < 0) || dist < minContactDistances[MolContactInfo.CC]) {
                                minContactDistances[MolContactInfo.CC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CC] = i;
                                contactAtomNumInResidueB[MolContactInfo.CC] = j;
                            }
                        }
                        else {
                            System.err.println("ERROR: Congrats, you found a bug in the atom contact type determination code (res " + a.getPdbNum() + " atom " + i + " / res " + b.getPdbNum() + " atom " + j + ").");
                            System.err.println("ERROR: Atom types are: i (PDB atom #" + x.getPdbAtomNum() + ") => " + x.getAtomType() + ", j (PDB atom #" + y.getPdbAtomNum() + ") => " + y.getAtomType() + ".");
                            Main.doExit(1);
                        }
                        
                        /*
                        if(i.equals(atomIndexOfBackboneN) && j.equals(atomIndexOfBackboneO)) {
                            // H bridge from backbone atom 'N' of residue a to backbone atom 'O' of residue b.
                            numPairContacts[MolContactInfo.HB]++;
                            // There can only be one of these so if we found it, simply update the distance.
                            minContactDistances[MolContactInfo.HB] = dist;
                        }

                        if(i.equals(atomIndexOfBackboneO) && j.equals(atomIndexOfBackboneN)) {
                            // H bridge from backbone atom 'O' of residue a to backbone atom 'N' of residue b.
                            numPairContacts[MolContactInfo.BH]++;
                            // There can only be one of these so if we found it, simply update the distance.
                            minContactDistances[MolContactInfo.BH] = dist;
                        }*/
                    }
                    
                    else if(Settings.getBoolean("plcc_B_alternate_aminoacid_contact_model_with_ligands")) {
                    
                    if(x.isProteinAtom() && y.isLigandAtom()) {
                        // *************************** protein - ligand contact *************************
                        numPairContacts[MolContactInfo.TT]++;
                        numTotalLigContactsPair++;
//                        System.out.println("New ProtLig: " + x.toString() + "/" + y.toString());
                        atomAtomContactType.add("PROTLIG");
                        donorAcceptor = new Atom[2];
                        donorAcceptor[0] = x;
                        donorAcceptor[1] = y;
                        atomAtomContacts.add(donorAcceptor);

                        // Check the exact contact type
                        if(i <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a backbone - ligand contact
                            numPairContacts[MolContactInfo.BL]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.BL] < 0) || dist < minContactDistances[MolContactInfo.BL]) {
                                minContactDistances[MolContactInfo.BL] = dist;
                                contactAtomNumInResidueA[MolContactInfo.BL] = i;
                                contactAtomNumInResidueB[MolContactInfo.BL] = j;
                            }

                        }
                        else {
                            // to be precise, this is a side chainName - ligand contact
                            numPairContacts[MolContactInfo.CL]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.CL] < 0) || dist < minContactDistances[MolContactInfo.CL]) {
                                minContactDistances[MolContactInfo.CL] = dist;
                                contactAtomNumInResidueA[MolContactInfo.CL] = i;
                                contactAtomNumInResidueB[MolContactInfo.CL] = j;
                            }
                        }

                    }
                    else if(x.isLigandAtom() && y.isProteinAtom()) {
                        // *************************** ligand - protein contact *************************
                        numPairContacts[MolContactInfo.TT]++;
                        numTotalLigContactsPair++;
//                        System.out.println("New LigProt: " + x.toString() + "/" + y.toString());
                        atomAtomContactType.add("LIGPROT");
                        donorAcceptor = new Atom[2];
                        donorAcceptor[0] = x;
                        donorAcceptor[1] = y;
                        atomAtomContacts.add(donorAcceptor);
                        
                        // Check the exact contact type
                        if(j <= numOfLastBackboneAtomInResidue) {
                            // to be precise, this is a ligand - backbone contact
                            numPairContacts[MolContactInfo.LB]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.LB] < 0) || dist < minContactDistances[MolContactInfo.LB]) {
                                minContactDistances[MolContactInfo.LB] = dist;
                                contactAtomNumInResidueA[MolContactInfo.LB] = i;
                                contactAtomNumInResidueB[MolContactInfo.LB] = j;
                            }

                        }
                        else {
                            // to be precise, this is a ligand - side chainName contact
                            numPairContacts[MolContactInfo.LC]++;

                            // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                            if((minContactDistances[MolContactInfo.LC] < 0) || dist < minContactDistances[MolContactInfo.LC]) {
                                minContactDistances[MolContactInfo.LC] = dist;
                                contactAtomNumInResidueA[MolContactInfo.LC] = i;
                                contactAtomNumInResidueB[MolContactInfo.LC] = j;
                            }
                        }
                            
                    }
                    else if(x.isLigandAtom() && y.isLigandAtom()) {
                        // *************************** ligand - ligand contact *************************
                        numPairContacts[MolContactInfo.TT]++;
                        numTotalLigContactsPair++;
//                        System.out.println("New LigLig: " + x.toString() + "/" + y.toString());
                        atomAtomContactType.add("LIGLIG");
                        donorAcceptor = new Atom[2];
                        donorAcceptor[0] = x;
                        donorAcceptor[1] = y;
                        atomAtomContacts.add(donorAcceptor);

                        // no choices here, ligands have no sub type
                        numPairContacts[MolContactInfo.LL]++;
                        
                        // update data if this is the first contact of this type or if it is better (smaller distance) than the old contact
                        if((minContactDistances[MolContactInfo.LL] < 0) || dist < minContactDistances[MolContactInfo.LL]) {
                            minContactDistances[MolContactInfo.LL] = dist;
                            contactAtomNumInResidueA[MolContactInfo.LL] = i;
                            contactAtomNumInResidueB[MolContactInfo.LL] = j;
                        }
                        

                        }
                    }
                    else if (x.isLigandAtom() || y.isLigandAtom()) {
                        // Nothing to do as we computing ligand contacts is disabled by the config.
                        // We still need to catch ligand molecules here, otherwise "detect" unknown contacts
                    }
                    else {
                        // *************************** unknown contact, wtf? *************************
                        // This branch should never be hit because atoms of type OTHER are ignored while creating the list of Atom objects
                        System.out.println("WARNING: One of the atoms " + x.getPdbAtomNum() + " and " + y.getPdbAtomNum() + " is of type UNKNOWN. Bug?");
                    }                                                            
                }                
                else {
                    // No atom contact for these 2 atoms, but there could be contacts between others
                }
            }
        }
        
        MolContactInfo piResults = calculatePiEffects(a, b);
        if(piResults != null) {
            
             numPairContacts[MolContactInfo.TT] += piResults.getNumContactsTotal();
             numPairContacts[MolContactInfo.NHPI] = piResults.getNumContactsNHPI();
             numPairContacts[MolContactInfo.PINH] = piResults.getNumContactsPINH();
             numPairContacts[MolContactInfo.CAHPI] = piResults.getNumContactsCAHPI();
             numPairContacts[MolContactInfo.PICAH] = piResults.getNumContactsPICAH();
             numPairContacts[MolContactInfo.CNHPI] = piResults.getNumContactsCNHPI();
             numPairContacts[MolContactInfo.PICNH] = piResults.getNumContactsPICNH();
             numPairContacts[MolContactInfo.SHPI] = piResults.getNumContactsSHPI();
             numPairContacts[MolContactInfo.PISH] = piResults.getNumContactsPISH();
             numPairContacts[MolContactInfo.XOHPI] = piResults.getNumContactsXOHPI();
             numPairContacts[MolContactInfo.PIXOH] = piResults.getNumContactsPIXOH();
             numPairContacts[MolContactInfo.PROCDHPI] = piResults.getNumContactsPROCDHPI();
             numPairContacts[MolContactInfo.PIPROCDH] = piResults.getNumContactsPIPROCDH();
             numPairContacts[MolContactInfo.CCAHCO] = piResults.getNumContactsCCACOH();
             numPairContacts[MolContactInfo.CCOCAH] = piResults.getNumContactsCCOCAH();
             numPairContacts[MolContactInfo.BCAHCO] = piResults.getNumContactsBCACOH();
             numPairContacts[MolContactInfo.BCOCAH] =  piResults.getNumContactsBCOCAH();
             
             minContactDistances[MolContactInfo.NHPI] = piResults.getNHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.NHPI] = piResults.getNHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.NHPI] = piResults.getNHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PINH] = piResults.getPINHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PINH] = piResults.getPINHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PINH] = piResults.getPINHContactAtomNumB();
             minContactDistances[MolContactInfo.CAHPI] = piResults.getCAHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.CAHPI] = piResults.getCAHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.CAHPI] = piResults.getCAHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PICAH] = piResults.getPICAHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PICAH] = piResults.getPICAHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PICAH] = piResults.getPICAHContactAtomNumB();
             minContactDistances[MolContactInfo.CNHPI] = piResults.getCNHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.CNHPI] = piResults.getCNHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.CNHPI] = piResults.getCNHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PICNH] = piResults.getPICNHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PICNH] = piResults.getPICNHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PICNH] = piResults.getPICNHContactAtomNumB();
             minContactDistances[MolContactInfo.SHPI] = piResults.getSHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.SHPI] = piResults.getSHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.SHPI] = piResults.getSHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PISH] = piResults.getPISHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PISH] = piResults.getPISHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PISH] = piResults.getPISHContactAtomNumB();
             minContactDistances[MolContactInfo.XOHPI] = piResults.getXOHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.XOHPI] = piResults.getXOHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.XOHPI] = piResults.getXOHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PIXOH] = piResults.getPIXOHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PIXOH] = piResults.getPIXOHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PIXOH] = piResults.getPIXOHContactAtomNumB();
             minContactDistances[MolContactInfo.PROCDHPI] = piResults.getPROCDHPIContactDist();
             contactAtomNumInResidueA[MolContactInfo.PROCDHPI] = piResults.getPROCDHPIContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PROCDHPI] = piResults.getPROCDHPIContactAtomNumB();
             minContactDistances[MolContactInfo.PIPROCDH] = piResults.getPIPROCDHContactDist();
             contactAtomNumInResidueA[MolContactInfo.PIPROCDH] = piResults.getPIPROCDHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.PIPROCDH] = piResults.getPIPROCDHContactAtomNumB();
             minContactDistances[MolContactInfo.CCAHCO] = piResults.getCCAHCOContactDist();
             contactAtomNumInResidueA[MolContactInfo.CCAHCO] = piResults.getCCAHCOContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.CCAHCO] = piResults.getCCAHCOContactAtomNumB();
             minContactDistances[MolContactInfo.CCOCAH] = piResults.getCCOCAHContactDist();
             contactAtomNumInResidueA[MolContactInfo.CCOCAH] = piResults.getCCOCAHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.CCOCAH] = piResults.getCCOCAHContactAtomNumB();
             minContactDistances[MolContactInfo.BCAHCO] = piResults.getBCAHCOContactDist();
             contactAtomNumInResidueA[MolContactInfo.BCAHCO] = piResults.getBCAHCOContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.BCAHCO] = piResults.getBCAHCOContactAtomNumB();
             minContactDistances[MolContactInfo.BCOCAH] = piResults.getBCOCAHContactDist();
             contactAtomNumInResidueA[MolContactInfo.BCOCAH] = piResults.getBCOCAHContactAtomNumA();
             contactAtomNumInResidueB[MolContactInfo.BCOCAH] = piResults.getBCOCAHContactAtomNumB();
        }
                               
        // Iteration through all atoms of the two residues is done
        if(numPairContacts[MolContactInfo.TT] > 0) {
            result = new MolContactInfo(numPairContacts, minContactDistances, contactAtomNumInResidueA, contactAtomNumInResidueB, a, b, CAdist, numTotalLigContactsPair, atomAtomContactType, atomAtomContacts);
        }
        else {
            result = null;
        }
        return(result);         // This is null if no contact was detected
        
       
     }


    /**
     * Prints an overview of all contacts to STDOUT.
     * @param rciList the residue contact information to consider
     */
    public static void showContactOverview(ArrayList<MolContactInfo> rciList) {

        ArrayList<MolContactInfo> contacts = rciList;
        MolContactInfo rci;
        Integer contactNum;

        System.out.println("  Numbr TypeA TypeB ResA# ResB# ResA ResB Dist #Cont");

        for(Integer i = 0; i < rciList.size(); i++) {

            rci = rciList.get(i);
            contactNum = i + 1;

            
            System.out.printf("   %4d   %3s   %3s   %3d   %3d  %3s  %3s  %3d    %2d\n", contactNum, rci.getResTypeStringA(), rci.getResTypeStringB(), rci.getDsspNumA(), rci.getDsspNumB(), rci.getName3A(), rci.getName3B(), rci.getMolPairDist(), rci.getNumContactsTotal());

        }

    }



    /**
     * Prints an overview of the contact statistics to STDOUT.
     * @param tnr the total number of residues, only required to be printed.
     */
    public static void showContactStatistics(final Integer tnr) {

        final Integer numTotalRes = tnr;

        System.out.printf("%10d - total number of contacts\n%10d - total number of residues\n", contact[0][0][0][0], numTotalRes);

        // ***** Print a priori frequencies row *****
        System.out.println("Amino acid a priori frequencies:");

        // header
        for(Integer i = 1; i < NUM_AAs; i++) {
            // header
            System.out.printf("%3s ", AminoAcid.intIDToName3(i));
        }
        System.out.print("\n");
        
        // frequencies
        for(Integer i = 1; i < NUM_AAs; i++) {
            System.out.printf("%3d ", contact[i][0][0][0]);
        }
        System.out.print("\n\n");


        // ***** Print contacts statistics matrix *****
        System.out.println("Amino acid atom level contact matrix:");


        // header
        System.out.print("    ");
        for(Integer i = 1; i < NUM_AAs; i++) {
            // header
            System.out.printf("%3s ", AminoAcid.intIDToName3(i));
        }
        System.out.print("\n");

        // Iterate through all amino acid pairs...
        for(Integer i = 1; i < NUM_AAs; i++) {

            System.out.printf("%3s ", AminoAcid.intIDToName3(i));

            for(Integer j = 1; j < NUM_AAs; j++ ) {

                System.out.printf("%3d ", contact[i][j][0][0]);

            }

            System.out.print("\n");
        }
    }


    /**
     * Writes formated contact data to file '<pdbid>.geo'.
     * @param rciList the list of resContactInfo objects representing the contacts to consider
     * @param gf the path to the geoFile
     * @param useLigands whether to use .geolig format (instead of .geo format) and consider ligand contacts (they will be 
     * ignored if this is 'false'). If this is true, each line will have additional fields at the
     * end which contain the ligand info.
     */
    public static void writeContacts(ArrayList<MolContactInfo> rciList, String gf, Boolean useLigands) {

        ArrayList<MolContactInfo> contacts = rciList;
        MolContactInfo rci = null;
        Integer contactNum = 0;
        String geoFile = gf;
        FileWriter geoFW = null;
        PrintWriter geoFH = null;

        Integer fHB1Dist, fHB2Dist, fCenterSphereRadiusResA, fCenterSphereRadiusResB, fResPairDist;
        Integer fBBContactDist, fBCContactDist, fCBContactDist, fCCContactDist;
        Integer fBLContactDist, fLBContactDist, fCLContactDist, fLCContactDist, fLLContactDist;

        System.out.println("  Contact output file set to '" + geoFile + "'.");

        try {
            geoFW = new FileWriter(geoFile);
            geoFH = new PrintWriter(geoFW);

            // Remove this header, it is only there for debugging purposes
            //if(Settings.getInteger("plcc_I_debug_level") >= 1) {
            //    geoFH.print("CNUM RNA RNB ID RA ID RA DST HB1 HB2 BBD IA IB CBD IA IB BCD IA IB CCD IA IB TT\n");
            //}

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + geoFile + "'.");
            e.printStackTrace();
            //System.exit(1);
        }
            
        for(Integer i = 0; i < contacts.size(); i++) {

            rci = contacts.get(i);
            contactNum = i + 1;

            //System.out.println("    Writing contact " + contactNum + ".");

            // The format of the resulting file is explained in detail in the file 'doc/help_pdbid.geo'.

            try {
                // prepare the values that need to be adapted to the geom_neo output format
                fHB1Dist = rci.getHB1Dist(); if(fHB1Dist > 999) { fHB1Dist = 999; } if(fHB1Dist < 0) { fHB1Dist = 0; }
                fHB2Dist = rci.getHB2Dist(); if(fHB2Dist > 999) { fHB2Dist = 999; } if(fHB2Dist < 0) { fHB2Dist = 0; }
                fCenterSphereRadiusResA = rci.getCenterSphereRadiusResA(); if(fCenterSphereRadiusResA > 99) { fCenterSphereRadiusResA = 99; }
                fCenterSphereRadiusResB = rci.getCenterSphereRadiusResB(); if(fCenterSphereRadiusResB > 99) { fCenterSphereRadiusResB = 99; }
                fResPairDist = rci.getMolPairDist(); if(fResPairDist > 999) { fResPairDist = 999; }

                fBBContactDist = rci.getBBContactDist(); if(fBBContactDist > 999) { fBBContactDist = 999; } if(fBBContactDist < 0) { fBBContactDist = 0; }
                fBCContactDist = rci.getBCContactDist(); if(fBCContactDist > 999) { fBCContactDist = 999; } if(fBCContactDist < 0) { fBCContactDist = 0; }
                fCBContactDist = rci.getCBContactDist(); if(fCBContactDist > 999) { fCBContactDist = 999; } if(fCBContactDist < 0) { fCBContactDist = 0; }
                fCCContactDist = rci.getCCContactDist(); if(fCCContactDist > 999) { fCCContactDist = 999; } if(fCCContactDist < 0) { fCCContactDist = 0; }


                // print first part of the line
                geoFH.printf("%-4d %3d %-3d %2d %-2d %2d %-2d %3d %3d %-3d ",
                            contactNum,
                            rci.getDsspNumA(),
                            rci.getDsspNumB(),
                            rci.getAAIDResA(),
                            fCenterSphereRadiusResA,
                            rci.getAAIDResB(),
                            fCenterSphereRadiusResB,
                            fResPairDist,
                            fHB1Dist,
                            fHB2Dist
                );
                // print 2nd part of the line (split into 2 function calls only to increase readability a bit)

                geoFH.printf("%3d %2d %-2d %3d %2d %-2d %3d %2d %-2d %3d %2d %-2d %2d",
                            fBBContactDist,
                            rci.getBBContactAtomNumA(),
                            rci.getBBContactAtomNumB(),

                            fCBContactDist,
                            rci.getCBContactAtomNumA(),
                            rci.getCBContactAtomNumB(),

                            fBCContactDist,
                            rci.getBCContactAtomNumA(),
                            rci.getBCContactAtomNumB(),

                            fCCContactDist,
                            rci.getCCContactAtomNumA(),
                            rci.getCCContactAtomNumB(),

                            rci.getNumContactsTotal()
                );

                // print the ligand contact part of the line
                if(useLigands) {

                    geoFH.printf(" ");

                    // prepare ligand values
                    fBLContactDist = rci.getBLContactDist(); if(fBLContactDist > 999) { fBLContactDist = 999; } if(fBLContactDist < 0) { fBLContactDist = 0; }
                    fLBContactDist = rci.getLBContactDist(); if(fLBContactDist > 999) { fLBContactDist = 999; } if(fLBContactDist < 0) { fLBContactDist = 0; }
                    fCLContactDist = rci.getCLContactDist(); if(fCLContactDist > 999) { fCLContactDist = 999; } if(fCLContactDist < 0) { fCLContactDist = 0; }
                    fLCContactDist = rci.getLCContactDist(); if(fLCContactDist > 999) { fLCContactDist = 999; } if(fLCContactDist < 0) { fLCContactDist = 0; }
                    fLLContactDist = rci.getLLContactDist(); if(fLLContactDist > 999) { fLLContactDist = 999; } if(fLLContactDist < 0) { fLLContactDist = 0; }

                    // write ligand data
                    geoFH.printf("%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d",
                            fBLContactDist,
                            rci.getBLContactAtomNumA(),
                            rci.getBLContactAtomNumB(),

                            fLBContactDist,
                            rci.getLBContactAtomNumA(),
                            rci.getLBContactAtomNumB(),

                            fCLContactDist,
                            rci.getCLContactAtomNumA(),
                            rci.getCLContactAtomNumB(),

                            fLCContactDist,
                            rci.getLCContactAtomNumA(),
                            rci.getLCContactAtomNumB(),

                            fLLContactDist,
                            rci.getLLContactAtomNumA(),
                            rci.getLLContactAtomNumB(),

                            rci.getNumLigContactsTotal()
                    );
                }
                geoFH.print("\n");


            } catch (Exception ew) {
                System.err.println("ERROR: Could not write info on contact " + contactNum + " to file '" + geoFile + "'.");
                ew.printStackTrace();
                //System.exit(1);
            }


        }

        // Remove this footer, is it there for DEBUG only
        //if(Settings.getInteger("plcc_I_debug_level") >= 1) {
        //    geoFH.print("EOF\n");
        //}
        

        geoFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!
        try {
            geoFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + geoFile + "'.");
            ex.printStackTrace();
            //System.exit(1);
        }

        System.out.println("  Wrote contact info on " + contactNum + " residue pairs to file '" + geoFile + "'.");

    }


    /**
     * Writes formated statistics to file 'cf', which usually should be 'con.set' (even though Windows may cry at this).
     * @param cf the output path where to write the conset file
     * @param pdb the pdbid of the protein, required for output only
     * @param tnr the total number of residues, required for output only
     */
    public static void writeStatistics(String cf, final String pdb, final Integer tnr) {



        String consetFile = cf;
        FileWriter conFW = null;
        PrintWriter conFH = null;
        Integer numTotalRes = tnr;

        System.out.println("  Statistics output file set to '" + consetFile + "'.");

        try {
            conFW = new FileWriter(consetFile);
            conFH = new PrintWriter(conFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + consetFile + "': '" + e.getMessage() + "'.");
            //e.printStackTrace();
            //System.exit(1);
        }

        conFH.println("AA distribution and contact statistics for protein " + pdb + ".\n");
        conFH.printf("%10d - total number of atom level contacts\n%10d - total number of residues\n\n", contact[0][0][0][0], numTotalRes);
  
        for(Integer i = 1; i < NUM_AAs; i++) {
            conFH.printf("%10d - a priori frequency of amino acid %2d (%s)\n", contact[i][0][0][0], i, AminoAcid.intIDToName3(i));
        }
        conFH.print("\n");

        // Iterate through all amino acid pairs...
        for(Integer i = 1; i < NUM_AAs; i++) {

            for(Integer j = 1; j < NUM_AAs; j++ ) {

                conFH.printf("%6d - contacts between amino acids %2d (%s) and %2d (%s)\n", contact[i][j][0][0], i, AminoAcid.intIDToName3(i), j, AminoAcid.intIDToName3(j));

                // ... and through all their atoms
                for(Integer k = 0; k <= AminoAcid.atomCountOfID(i); k++ ) {

                    conFH.printf("\n");

                    if( k == 1 ) { conFH.printf("\n"); }

                    for(Integer l = 0; l <= AminoAcid.atomCountOfID(i); l++) {

                        if( l == 1 ) { conFH.printf("      "); }

                        // prints contact statistics for the atoms
                        conFH.printf("%6d ", contact[i][j][k][l]);
                    }
                }
                conFH.printf("\n\n");
            }
        }

      
        conFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            conFW.close();
        } catch(Exception ex) {
            System.err.println("WARNING: Could not close FileWriter for file '" + consetFile + "'.");
            //ex.printStackTrace();
            //System.exit(1);
        }

        System.out.println("  Wrote contact statistics to file '" + consetFile + "'.");
        
    }

    /** Write formated statistics to file out_file.
     * Takes a MolContactInfo object as input that stores all the necessary information to generate output files.
 Two output files are generated: One with information about the different contact types, and another one 
 with information about the detected atom-atom contacts.
     * @param results The calculated results that should be saved.
     * @param out_file Where the files are saved.
     */
    public static void writePPIstatistics(ArrayList<MolContactInfo> results, String out_file) {
        
        // Initialize all different contact types
        int BBHB, BBBH, IVDW, ISS, BCHB, BCBH, CBHB, CBBH, CCHB, CCBH, BB, CB, 
                BC, CC, BL, LB, CL, LC, LL, NHPI, PINH, CAHPI, PICAH, CNHPI, 
                PICNH, SHPI, PISH, XOHPI, PIXOH, PROCDHPI, PIPROCDH, CCACOH, 
                CCOCAH, BCACOH, BCOCAH;
        BBHB = BBBH = IVDW = ISS = BCHB = BCBH = CBHB = CBBH = CCHB = CCBH = BB = CB = 
                BC = CC = BL = LB = CL = LC = LL = NHPI = PINH = CAHPI = PICAH = CNHPI = 
                PICNH = SHPI = PISH = XOHPI = PIXOH = PROCDHPI = PIPROCDH = CCACOH = 
                CCOCAH = BCACOH = BCOCAH = 0;
        
        ArrayList<String> atomAtomContactTypes = new ArrayList<String>();
        ArrayList<Atom[]> atomAtomContacts = new ArrayList<Atom[]>();
        
        // Go through the calculated results and get the important information.
        for(MolContactInfo result : results) {
            BBHB += result.getNumContactsBBHB();
            BBBH += result.getNumContactsBBBH();
            IVDW += result.getNumContactsIVDW();
            ISS += result.getNumContactsISS();
            BCHB += result.getNumContactsBCHB();
            BCBH += result.getNumContactsBCBH();
            CBHB += result.getNumContactsCBHB();
            CBBH += result.getNumContactsCBBH();
            CCHB += result.getNumContactsCCHB();
            CCBH += result.getNumContactsCCBH();
            BB += result.getNumContactsBB();
            CB += result.getNumContactsCB();
            BC += result.getNumContactsBC();
            CC += result.getNumContactsCC();
            BL += result.getNumContactsBL();
            LB += result.getNumContactsLB();
            CL += result.getNumContactsCL();
            LC += result.getNumContactsLC();
            LL += result.getNumContactsLL();
            NHPI += result.getNumContactsNHPI();
            PINH += result.getNumContactsPINH();
            CAHPI += result.getNumContactsCAHPI();
            PICAH += result.getNumContactsPICAH();
            CNHPI += result.getNumContactsCNHPI();
            PICNH += result.getNumContactsPICNH();
            SHPI += result.getNumContactsSHPI();
            PISH += result.getNumContactsPISH();
            XOHPI += result.getNumContactsXOHPI();
            PIXOH += result.getNumContactsPIXOH();
            PROCDHPI += result.getNumContactsPROCDHPI();
            PIPROCDH += result.getNumContactsPIPROCDH();
            CCACOH += result.getNumContactsCCACOH();
            CCOCAH += result.getNumContactsCCOCAH();
            BCACOH += result.getNumContactsBCACOH();
            BCOCAH += result.getNumContactsBCOCAH();
            
            atomAtomContactTypes.addAll(result.getAtomAtomContactTypes());
            atomAtomContacts.addAll(result.getAtomAtomContacts());
            
            
        }
        
        StringBuilder sb = new StringBuilder();
        String lineSep = System.lineSeparator();
        
        // Write the atom-atom contact information.
        for(int x = 0; x < atomAtomContactTypes.size(); x++) {
            sb.append(atomAtomContactTypes.get(x));
            sb.append(";");
            sb.append(atomAtomContacts.get(x)[0]);
            sb.append(";");
            sb.append(atomAtomContacts.get(x)[1]);
            sb.append(lineSep);
        }
        
        File aacs = new File(out_file + "_atom_atom_contacts.csv");
        try {
            FileWriter fw = new FileWriter(aacs.getAbsoluteFile());
            try (BufferedWriter bw = new BufferedWriter(fw)) {
                bw.append(sb);
            }
        } catch (IOException ex) {
            DP.getInstance().e("main", "Unable to write atom_atom_contacts.csv file.");
        }
        
        
        // Write the contact type information.
        sb = new StringBuilder();
        sb.append("BBHB: ");
        sb.append(BBHB);
        sb.append(lineSep);
        sb.append("BBBH: ");
        sb.append(BBBH);
        sb.append(lineSep);
        sb.append("IVDW: ");
        sb.append(IVDW);
        sb.append(lineSep);
        sb.append("ISS: ");
        sb.append(ISS);
        sb.append(lineSep);
        sb.append("BCHB: ");
        sb.append(BCHB);
        sb.append(lineSep);
        sb.append("BCBH: ");
        sb.append(BCBH);
        sb.append(lineSep);
        sb.append("CBHB: ");
        sb.append(CBHB);
        sb.append(lineSep);
        sb.append("CBBH: ");
        sb.append(CBBH);
        sb.append(lineSep);
        sb.append("CCHB: ");
        sb.append(CCHB);
        sb.append(lineSep);
        sb.append("CCBH: ");
        sb.append(CCBH);
        sb.append(lineSep);
        sb.append("BB: ");
        sb.append(BB);
        sb.append(lineSep);
        sb.append("CB: ");
        sb.append(CB);
        sb.append(lineSep);
        sb.append("BC: ");
        sb.append(BC);
        sb.append(lineSep);
        sb.append("CC: ");
        sb.append(CC);
        sb.append(lineSep);
        sb.append("BL: ");
        sb.append(BL);
        sb.append(lineSep);
        sb.append("LB: ");
        sb.append(LB);
        sb.append(lineSep);
        sb.append("CL: ");
        sb.append(CL);
        sb.append(lineSep);
        sb.append("LC: ");
        sb.append(LC);
        sb.append(lineSep);
        sb.append("LL: ");
        sb.append(LL);
        sb.append(lineSep);
        sb.append("NHPI: "); 
        sb.append(NHPI); 
        sb.append(lineSep);
        sb.append("PINH: "); 
        sb.append(PINH); 
        sb.append(lineSep);
        sb.append("CAHPI: "); 
        sb.append(CAHPI); 
        sb.append(lineSep);
        sb.append("PICAH: ");
        sb.append(PICAH);
        sb.append(lineSep);
        sb.append("CNHPI: "); 
        sb.append(CNHPI); 
        sb.append(lineSep);
        sb.append("PICNH: "); 
        sb.append(PICNH); 
        sb.append(lineSep);
        sb.append("SHPI: "); 
        sb.append(SHPI); 
        sb.append(lineSep);
        sb.append("PISH: "); 
        sb.append(PISH); 
        sb.append(lineSep);
        sb.append("XOHPI: "); 
        sb.append(XOHPI); 
        sb.append(lineSep);
        sb.append("PIXOH: "); 
        sb.append(PIXOH); 
        sb.append(lineSep);
        sb.append("PROCDHPI: "); 
        sb.append(PROCDHPI); 
        sb.append(lineSep);
        sb.append("PIPROCDH: "); 
        sb.append(PIPROCDH); 
        sb.append(lineSep);
        sb.append("CCACOH: "); 
        sb.append(CCACOH); 
        sb.append(lineSep);
        sb.append("CCOCAH: "); 
        sb.append(CCOCAH); 
        sb.append(lineSep);
        sb.append("BCACOH: "); 
        sb.append(BCACOH); 
        sb.append(lineSep);
        sb.append("BCOCAH: "); 
        sb.append(BCOCAH); 
        sb.append(lineSep);
        sb.append(lineSep);
        
        File stats = new File(out_file + ".stats");
        try {
            FileWriter fw = new FileWriter(stats.getAbsoluteFile());
            try (BufferedWriter bw = new BufferedWriter(fw)) {
                bw.append(sb);
            }
        } catch (IOException ex) {
            DP.getInstance().e("main", "Unable to write .stats file.");
        }
        
        
        // Information on residue level.
        /*sb = new StringBuilder();
        sb.append("Checked residue contacts: ");
        sb.append(numResContactsChecked);
        sb.append(lineSep);
        sb.append("Resdues: ");
        sb.append(rs);
        sb.append(lineSep);
        sb.append("Posible: ");
        sb.append(numResContactsPossible);
        sb.append(lineSep);
        sb.append("Found: ");
        sb.append(contactInfo.size());
        sb.append(lineSep);
        sb.append("Impossible (collision spheres check): ");
        sb.append(numResContactsImpossible);
        
        File res = new File(out_file + ""_res.stats");
        try {
            FileWriter fw = new FileWriter(res.getAbsoluteFile());
            try (BufferedWriter bw = new BufferedWriter(fw)) {
                bw.append(sb);
            }
        } catch (IOException ex) {
            DP.getInstance().e("main", "Unable to write _res.stats file.");
        }*/
    }

    /**
     * Determines the maximum center sphere radius of all residues. This is the distance from the (center of the) central atom to the (center of the) atom which
     * is farthest from that atom.
     * @param mols a residue list
     * @return the radius of the largest center sphere
     */
    public static Integer getGlobalMaxCenterSphereRadius(List<Molecule> mols) {

        Molecule m = null;
        Integer maxRad, curRad;
        maxRad = curRad = 0;

        for(Integer i = 0; i < mols.size(); i++) {
            m = mols.get(i);
            curRad = m.getSphereRadius();

            if(curRad > maxRad) {
                maxRad = curRad;
            }
        }
     
        return(maxRad);
    }

    
    /**
     * Determines the maximum 3D distance between residues which are sequential neighbors. Note that the distance can be pretty large because
     * of ligands which are always listed at the end even though they may not be in the vicinity of the last protein residue, i.e., the one
     * directly preceeding them.
     * @param mols a list of molecules which to consider
     * @param centroidMethod optional argument (default: false) for using centroids instead of C_alpha as center
     * @return the maximum distance between two consecutive protein residues in the primary sequence
     */
    public static Integer getGlobalMaxSeqNeighborResDist(List<? extends Molecule> mols) {
              
        Molecule r, s;
        r = s = null;
        Integer maxDist, curDist, rID, sID, rT, sT;
        maxDist = curDist  = rID = sID = rT = sT = 0;


        // Iterate through residues in sequential order (DSSP numbering) and determine
        //  the maximal distance (center to center) of all pairs of residues that are
        //  neighbors in the amino acid sequence.
        for(Integer i = 0; i < mols.size() - 1; i++) {

            r = mols.get(i);         // Is this really DSSP ordering? how can we check quickly?
            s = mols.get(i + 1);     // NOTE: Yes, it is DSSP ordering since residues are parsed from the
                                    //        DSSP file and residues in there are in DSSP ordering.
                                    
            // despite the note from above, simply check that this holds true when chain sphere speedup is on (since lists are passed to the function differently)
            if (Settings.getBoolean("plcc_B_chain_spheres_speedup")) {
                if (s.getDsspNum() != r.getDsspNum() + 1) {
                    if (! Settings.getBoolean("plcc_B_no_warn")) {
                        DP.getInstance().w("Function getGlobalMaxSeqNeighborResDist: " + r.getFancyName() + " (chain " + r.getChainID() + ") and " + 
                                s.getFancyName() + " (chain " + s.getChainID() + ") " +
                                "are not consecutive by DSSP order! This may interfere with sequence neighbor speedup. Results should be correct, though.");
                    }
                }
            }
            
            curDist = r.distTo(s);
            if(curDist > maxDist) {
                maxDist = curDist;
                rID = r.getDsspNum();
                rT = r.getType();
                sID = s.getDsspNum();
                sT = s.getType();
            }
        }

        if (Settings.getInteger("plcc_I_debug_level") >= 3) {
            System.out.println("  Neighbor residues " + rID + " (type " + rT + ") and " + sID + " (type " + sT + ") found in distance " + maxDist + ".");
        }
        
        return(maxDist);
    }


    /**
     * Writes the chains file, a file containing a list of all chainName identifiers, separated by spaces.
     * @param chainsFile the path to the output file
     * @param pdbid the PDB ID of the protein, required for output only
     * @param chains a list of Chain objects, should be all chains of this protein
     */
    public static void writeChains(String chainsFile, String pdbid, ArrayList<Chain> chains) {

        FileWriter chainFW = null;
        PrintWriter chainFH = null;

        // open files
        try {
            chainFW = new FileWriter(chainsFile);
            chainFH = new PrintWriter(chainFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + chainsFile + "'.");
            e.printStackTrace();
            Main.doExit(1);
        }


        // finally: write stuff
        //chainFH.printf("%s: ", pdbid);

        for(Integer i = 0; i < chains.size(); i++) {
            chainFH.printf("%s ", chains.get(i).getPdbChainID());
        }

        //chainFH.printf("%d\n", chains.size());
        //chainFH.print("# Info is always parsed from the first non-comment line of this file so don't mess that up.\n");
        //chainFH.print("# Line format is as follows: '<pdbid>: [<chainID> ...] <total_number_of_chains>'\n");

        chainFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            chainFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + chainsFile + "'.");
            //ex.printStackTrace();
            //Main.doExit(1);
        }

        if(! Settings.getBoolean("plcc_B_silent")) {
            System.out.println("  Wrote chain info to file '" + chainsFile + "'.");       
        }

    }
    
    
    /**
     * Writes the residue info file that maps PDB residue IDs to DSSP residue IDs for all residues of the given chainName. 
     * @param mapFile the path to the output file
     * @param c the chainName to consider (all residues of this chainName will be used)
     */
    public static void writeSSEMappings(String mapFile, Chain c, String pdbid) {
        String s = "# SSE mappings for protein " + pdbid + " chain " + c.getPdbChainID() + " follow in format <PDB res number> <DSSP res number> <DSSP assignment> <PLCC assignment>";
        ArrayList<Residue> res = c.getResidues();
                
        
        for (Residue r : res) {
            s += "" + r.getPdbNum() + "|" + r.getDsspNum() + "|" + r.getSSEStringDssp() + "|" + r.getSSETypePlcc() + "\n";
        }
        
        IO.stringToTextFile(mapFile, s);
    }
    
    
    /**
     * Writes the residue info file that maps PDB residue IDs to DSSP residue IDs for all residues of the given chainName. 
     * @param mapFile the path to the output file
     * @param c the chainName to consider (all residues of this chainName will be used)
     */
    public static void writeResMappings(String mapFile, Chain c) {

        FileWriter mapFW = null;
        PrintWriter mapFH = null;
        ArrayList<Residue> res = c.getResidues();

        // open files
        try {
            mapFW = new FileWriter(mapFile);
            mapFH = new PrintWriter(mapFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + mapFile + "'.");
            //e.printStackTrace();            
            Main.doExit(1);
        }


        for (Residue r : res) {
            mapFH.print("PDB|" + r.getPdbNum() + "|DSSP|" + r.getDsspNum() + "\n");
        }
        
        //chainFH.printf("%d\n", chains.size());
        //chainFH.print("# Info is always parsed from the first non-comment line of this file so don't mess that up.\n");
        //chainFH.print("# Line format is as follows: '<pdbid>: [<chainID> ...] <total_number_of_chains>'\n");

        mapFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            mapFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + mapFile + "'.");
            //ex.printStackTrace();
            //System.exit(1);
        }

        System.out.println("  Wrote PDB/DSSP residue mapping info to file '" + mapFile + "'.");       

    }
    

     /**
     * Writes the ligand file.
     * @param ligFile the path to the output ligand file
     * @param pdbid the PDB ID of the current protein
     * @param ligands a list of residues which are expected to be ligands
     */
    public static void writeLigands(String ligFile, String pdbid, List<Residue> ligands) {

        FileWriter ligFW = null;
        PrintWriter ligFH = null;
        Residue r = null;

        // open files
        try {
            ligFW = new FileWriter(ligFile);
            ligFH = new PrintWriter(ligFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + ligFile + "'.");
            e.printStackTrace();
            Main.doExit(1);
        }


        // finally: write stuff
        for(Integer i = 0; i < ligands.size(); i++) {

            r = ligands.get(i);

            if(r.isLigand()) {
                ligFH.printf("%s %s %s%d %d %s %s %s %s\n", r.getModelID(), r.getChainID(), r.getChainID(), r.getPdbNum(), r.getDsspNum(), r.getName3(), r.getLigName(), r.getLigFormula(), r.getLigSynonyms());
            }
        }

        ligFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            ligFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + ligFile + "'.");
            ex.printStackTrace();
            Main.doExit(1);
        }

        System.out.println("  Wrote ligand info to file '" + ligFile + "'.");

    }



    /**
     * Writes the model file.
     * @param modelsFile the output file path
     * @param pdbid the PDB ID of the current protein
     * @param modelIDs a list of String which represent the IDs of all models contained in the current PDB file
     */
    public static void writeModels(String modelsFile, String pdbid, ArrayList<String> modelIDs) {

        FileWriter modelFW = null;
        PrintWriter modelFH = null;

        // open files
        try {
            modelFW = new FileWriter(modelsFile);
            modelFH = new PrintWriter(modelFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + modelsFile + "'.");
            e.printStackTrace();
            System.exit(1);
        }

        // finally: write stuff
        for(Integer i = 0; i < modelIDs.size(); i++) {
            modelFH.printf("%s ", modelIDs.get(i));
        }

        modelFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            modelFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + modelsFile + "'.");
            //ex.printStackTrace();
            //System.exit(1);
        }

        System.out.println("  Wrote model info to file '" + modelsFile + "'.");
    }


    /**
     * Writes a modified version of the original DSSP file. This version includes info on the ligand residues that is appended after the last line of the regular DSSP file.
     * @param dsspFile the input DSSP file that will be parsed for info
     * @param dsspLigFile the output path of the DSSPLig file, which is generated by this function by adding ligand info lines to the DSSP data
     * @param contacts the contacts to consider
     * @param res the residues to consider
     */
    public static Boolean writeDsspLigFile(String dsspFile, String dsspLigFile, ArrayList<MolContactInfo> contacts, ArrayList<Residue> res) {
        
        DP.getInstance().w("writeDsspLigFile(): This function is deprecated, use writeOrderedDsspLigFile() instead.\n");
        
        FileWriter dsspLigFW = null;
        PrintWriter dsspLigFH = null;

        // copy the DSSP file first, we will then append the ligand-specific lines to the copied file
        try {
            copyFile(new File(dsspFile), new File(dsspLigFile));
        } catch (Exception ef) {
            System.err.println("ERROR: Could not copy file '" + dsspFile + "' to '" + dsspLigFile + "'.");
            ef.printStackTrace();
            Main.doExit(1);
        }

        System.out.println("  DSSP ligand output file set to '" + dsspLigFile + "', file created.");

        // open files
        try {
            dsspLigFW = new FileWriter(dsspLigFile, true);  // The 2nd (boolean) parameter sets APPEND mode
            dsspLigFH = new PrintWriter(dsspLigFW);

        }
        catch (Exception e) {
            System.err.println("ERROR: Could not write to file '" + dsspLigFile + "'.");
            e.printStackTrace();
            Main.doExit(1);
        }

        // now get and write the ligand lines
        Residue r;
        Locale loc = Locale.ENGLISH;
        for(Integer i = 0; i < res.size(); i++) {

            r = res.get(i);

            if(r.isLigand()) {

                // Now write the line for this ligand in DSSP format (see http://swift.cmbi.ru.nl/gv/dssp/).
                // We do this in several printf() commands to make this code easier to read even though it may be slower (I guess the FileWrite has buffering so it shouldn't be that bad).
                // example line follows after header line (between the ticks):
                //         #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
                //      '   47   47 A E  E    S-c   29   0A  71    -19,-1.8   -17,-2.8     2,-0.0     2,-0.3  -0.981  70.4-156.7-136.2 150.2   43.9  -12.9   14.8'

                // Print DSSP residue number, PDB residue number, chainName, AA name in 1 letter code and SSE summary letter for ligand
                //      '   47   47 A E  E'
                dsspLigFH.printf(loc, "  %3d  %3d %1s %1s  %1s", r.getDsspNum(), r.getPdbNum(), r.getChainID(), r.getAAName1(), Settings.get("plcc_S_ligSSECode"));

                // Print structure detail block (empty for ligand), beta bridge 1 partner residue number (always 0 for ligands), beta bridge 2 partner residue number (always 0 for ligands),
                //  bet sheet label (empty (" ") for ligands) and solvent accessible surface (SAS) of this residue (not required by PTGL, just set to some value)
                //      '   47   47 A E  E    S-c   29   0A  71'
                dsspLigFH.printf(loc, " %7s %3d %3d%1s %3d", "       ", 0, 0, " ", Settings.getInteger("plcc_I_ligSAS"));

                // Print the 4 possible H-bridges: N-H-->O, O-->H-N, N-H-->O, O-->H-N. First value (%4d) is the relative residue to which the H bond leads (e.g.: '-15' := 15 reasidues before this one, '4':= 4 residues after this one).
                //  The float value (%3.1f) is the strength of the H bond in kcal/mol. We dont use these and leave them empty for now.
                //      '   47   47 A E  E    S-c   29   0A  71    -19,-1.8   -17,-2.8     2,-0.0     2,-0.3'
                dsspLigFH.printf(loc, "   %4d,%4.1f  %4d,%4.1f  %4d,%4.1f  %4d,%4.1f", 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0);

                // Print the TCO, KAPPA, ALPHA, PHI and PSI angles. Then the center atom coordinates. That's it.
                dsspLigFH.printf(loc, "  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f\n", -0.5, 55.5, 55.5, 55.5, 55.5, (r.getCenterAtom().getCoordX() / 10.0f), (r.getCenterAtom().getCoordY() / 10.0f), (r.getCenterAtom().getCoordZ() / 10.0f));
            }

        }

        dsspLigFH.close();  // If it isn't closed it won't flush the buffer and large parts of the file will be missing!

        try {
            dsspLigFW.close();
        } catch(Exception ex) {
            System.err.println("ERROR: Could not close FileWriter for file '" + dsspLigFile + "'.");
            ex.printStackTrace();
            Main.doExit(1);
        }

        System.out.println("  Wrote DSSP ligand info to file '" + dsspLigFile + "'.");


        return(true);
    }

    
    /**
     * Writes an ordered dssplig file. Ordered means that the ligand residues are not simply appended to the end of the file, but each ligand is inserted as the
 last residue of the chainName it is associated with. This is required for PTGL compatibility (don't ask).
     * @param dsspFile the input DSSP file that will be parsed for info
     * @param dsspLigFile the output path of the DSSPLig file, which is generated by this function by adding ligand info lines to the DSSP data
     * @param res the residues to consider
     * @return true if it worked out. Note though that this is considered critical.
     * 
     */
    public static Boolean writeOrderedDsspLigFile(String dsspFile, String dsspLigFile, List<Residue> res) {

        File dFile = new File(dsspFile);
        File dligFile = new File(dsspLigFile);

        // copy the DSSP file first, we will then append the ligand-specific lines to the copied file
        try {
            copyFile(dFile, dligFile);
        } catch (Exception ef) {
            System.err.println("ERROR: Could not copy file '" + dsspFile + "' to '" + dsspLigFile + "': '" + ef.getMessage() + "'.");
            ef.printStackTrace();
            System.exit(1);
        }

        System.out.println("  DSSP ligand output file set to '" + dsspLigFile + "', file created.");

        Residue r;

        for(Integer i = 0; i < res.size(); i++) {

            r = res.get(i);

            if(r.isLigand()) {
                
                dligFile = new File(dsspLigFile);

                try {
                    insertLigandLineIntoDsspligFile(dligFile, getLastLineOfChain(dligFile, r.getChainID()), r);
                } catch(Exception cf) {
                    System.err.println("ERROR: Failed to insert line for ligand residue '" + r.getFancyName() + "' into dssplig file: '" + cf.getMessage() + "'.");
                    cf.printStackTrace();
                    System.exit(1);
                }

            }
        }

        return(true);

    }


    /**
     * Determines the last line number in a DSSP file that contains information on a residue that is part of chainName 'chainID' (i.e., on the last residue of that chainName).
     * @param dsspligFile the DSSPLig file to parse
     * @param chainID the chainName ID to consider
     * @return the line number of the last line belonging to chainName chainID
     */
    public static Integer getLastLineOfChain(File dsspligFile, String chainID) {

        Integer ln = -1;
        String line = "";
        String curChain = "";
        Boolean dataFound = false;

        ArrayList<String> dl = FileParser.slurpFile(dsspligFile.getAbsolutePath());

        for(Integer i = 0; i < dl.size(); i++) {
            line = dl.get(i);

            if(dataFound) {
                curChain = line.substring(11, 12);
                if(curChain.equals(chainID)) {
                    ln = i;
                }
            }
            else {
                if(line.charAt(2) == '#') {
                    dataFound = true;
                }
            }
            

        }

        if(ln == - 1) {
            System.err.println("ERROR: getLastLineofChain(): DSSPLIG file '" + dsspligFile.getName() + "' contains no information on requested chain '" + chainID + "'.");
            System.exit(1);
        }
        else {
            //System.out.println("DEBUG: Last line of chainName '" + chainID + "' in DSSPLIG file '" + dsspligFile.getName() + "' is " + ln + ".");
        }
        return(ln);
    }



    /**
     * Inserts information on the residue 'r' (which is expected to be a ligand residue) into the DSSP file 'inFile' at line number 'lineno'. Obviously, the line numbers of all 
     * following lines will be increased by 1.
     * @param inFile the input file that will be changed
     * @param lineno the line number
     * @param r the residue that should be added at the specified line number
     */
    public static void insertLigandLineIntoDsspligFile(File inFile, int lineno, Residue r)
       throws Exception {
     // temp file
     File outFile = new File("dssplig.tmp");

     // input
     FileInputStream fis  = new FileInputStream(inFile);
     BufferedReader in = new BufferedReader
         (new InputStreamReader(fis));

     // output
     FileOutputStream fos = new FileOutputStream(outFile);
     PrintWriter out = new PrintWriter(fos);

     String thisLine = "";
     Locale loc = Locale.ENGLISH;       // required for proper format!
     int i = 0;
     while ((thisLine = in.readLine()) != null) {
         
       out.println(thisLine);

       if(i == lineno) {
           // Now write the line for this ligand in DSSP format (see http://swift.cmbi.ru.nl/gv/dssp/).
           
            // We do this in several printf() commands to make this code easier to read even though it may be slower (I guess the FileWrite has buffering so it shouldn't be that bad).
            // example line follows after header line (between the ticks):
            //         #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
            //      '   47   47 A E  E    S-c   29   0A  71    -19,-1.8   -17,-2.8     2,-0.0     2,-0.3  -0.981  70.4-156.7-136.2 150.2   43.9  -12.9   14.8'

            // Print DSSP residue number, PDB residue number, chainName, AA name in 1 letter code and SSE summary letter for ligand
            //      '   47   47 A E  E'
            out.printf(loc, "  %3d  %3d %1s %1s  %1s", r.getDsspNum(), r.getPdbNum(), r.getChainID(), r.getAAName1(), Settings.get("plcc_S_ligSSECode"));

            // Print structure detail block (empty for ligand), beta bridge 1 partner residue number (always 0 for ligands), beta bridge 2 partner residue number (always 0 for ligands),
            //  bet sheet label (empty (" ") for ligands) and solvent accessible surface (SAS) of this residue (not required by PTGL, just set to some value)
            //      '   47   47 A E  E    S-c   29   0A  71'
            out.printf(loc, " %7s %3d %3d%1s %3d", "       ", 0, 0, " ", Settings.getInteger("plcc_I_ligSAS"));

            // Print the 4 possible H-bridges: N-H-->O, O-->H-N, N-H-->O, O-->H-N. First value (%4d) is the relative residue to which the H bond leads (e.g.: '-15' := 15 reasidues before this one, '4':= 4 residues after this one).
            //  The float value (%3.1f) is the strength of the H bond in kcal/mol. We dont use these and leave them empty for now.
            //      '   47   47 A E  E    S-c   29   0A  71    -19,-1.8   -17,-2.8     2,-0.0     2,-0.3'
            out.printf(loc, "   %4d,%4.1f  %4d,%4.1f  %4d,%4.1f  %4d,%4.1f", 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0);

            // Print the TCO, KAPPA, ALPHA, PHI and PSI angles. Then the center atom coordinates. That's it.
            out.printf(loc, "  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f\n", -0.5, 55.5, 55.5, 55.5, 55.5, (r.getCenterAtom().getCoordX() / 10.0f), (r.getCenterAtom().getCoordY() / 10.0f), (r.getCenterAtom().getCoordZ() / 10.0f));
       }
       
       i++;
    }

    out.flush();
    out.close();
    in.close();

    inFile.delete();
    outFile.renameTo(inFile);
  }


    /**
     * Copies file 'in' to file 'out' on file system level.
     * @param in the input file to be copied
     * @param out the output file
     */
    public static void copyFile(File in, File out)
        throws IOException
    {
        FileChannel inChannel = new
            FileInputStream(in).getChannel();
        FileChannel outChannel = new
            FileOutputStream(out).getChannel();
        try {
            inChannel.transferTo(0, inChannel.size(),
                    outChannel);
        }
        catch (IOException e) {
            throw e;
        }
        finally {
            if (inChannel != null) inChannel.close();
            if (outChannel != null) outChannel.close();
        }
    }


    /**
     * Prints a short note on how to use this program to STDOUT. Called if the user runs the program without command line arguments (most likely because he does not know about them).
     */
    public static void usage_short() {
        System.out.println("This program is part of VPLG, http://vplg.sourceforge.net. Copyright MolBI Group 2012 - 2014.");
        System.out.println("VPLG is free software and comes without any warranty. See LICENSE for details.");        
        System.out.println("USAGE:         java -jar plcc.jar <pdbid> [OPTIONS]");
        System.out.println("DETAILED HELP: java -jar plcc.jar --help");
    }
    
    
    /**
     * Prints detailed development usage info to STDOUT. The options printed here are internal experimental options, for development only.     
     */
    public static void usagedev() {
        System.out.println("DEV OPTIONS:");
        System.out.println("     Unsupported development-only options follow. For internal usage only.");
        System.out.println("       --alt-aa-contacts          : use alternate AA contact model by AS. Exits after AAGs.");
        System.out.println("       --alt-aa-contacts-ligands  : use alternate AA contact model including ligands by AS. Exits after AAGs.");
    }
    
    /**
     * Prints detailed usage info to STDOUT.
     */
    public static void usage() {
        System.out.println("USAGE: java -jar plcc.jar <pdbid> [OPTIONS]");
        System.out.println("       java -jar plcc.jar --help");
        System.out.println("valid OPTIONS are: ");
        System.out.println("-a | --include-coils       : convert the SSE type of all ignored residues to C (coil) and include coils in the graphs (may split other SSEs)");
        System.out.println("-b | --draw-plcc-fgs <f>   : read graph in plcc format from file <f> and draw it and all its folding graphs, then exit (pdbid will be ignored)*");
        System.out.println("-B | --force-backbone      : add contacts of special type 'backbone' between all SSEs of a graph in sequential order (N to C terminus)");
        System.out.println("-c | --dont-calc-graphs    : do not calculate SSEs contact graphs, stop after residue level contact computation");
        System.out.println("-C | --create-config       : creates a default config file if none exists yet, then exits.*");
        System.out.println("-D | --debug <level>       : set debug level (0: off, 1: normal debug output. >=2: detailed debug output, up to 4 currently)  [DEBUG]");
        System.out.println("-d | --dsspfile <dsspfile> : use input DSSP file <dsspfile> (instead of assuming '<pdbid>.dssp')");
        System.out.println("     --gz-dsspfile <f>     : use gzipped input DSSP file <f>.");
        System.out.println("-e | --force-chain <c>     : only handle the chain with chain ID <c>.");
        System.out.println("-E | --separate-contacts   : separate contact computation by chain (way faster but disables all functions which require inter-chain contacts (stats, complex graphs)");
        System.out.println("-f | --folding-graphs      : also handle foldings graphs and compute their linear notations. See -s if you also want them drawn.");
        System.out.println("-g | --sse-graphtypes <l>  : compute only the SSE graphs in list <l>, e.g. 'abcdef' = alpha, beta, alhpabeta, alphalig, betalig and alphabetalig.");
        System.out.println("-G | --complex-graphs      : compute and output complex graphs. Disables contact separation (see -E) if used.");
        System.out.println("-h | --help                : show this help message and exit");
        System.out.println("-H | --output-www-with-core: add HTML navigation files and core files to the subdir tree (implies -k, -W). ");
        //System.out.println("-i | --ignoreligands       : ignore ligand contacts in geom_neo format output files [DEBUG]");
        System.out.println("-i | --aa-graphs           : compute and output amino acid-based graphs as well. In these graphs, each vertes is an amino acid instead of a SSE. Computes per-chain and per-PDB file AAGs.");        
        System.out.println("     --aa-graphs-pdb       : computes amino acid graphs (see above), but only the graphs for the whole PDB file (all chains combined into one graph).");        
        System.out.println("     --aa-graphs-chain     : computes amino acid graphs (see above), but only the graphs which model a single chain each (n graphs for a PDB file with n chains).");                
        System.out.println("-I | --mmCIF-parser        : uses mmCIF parser for provided file (looks for .cif file)");
        System.out.println("-j | --ddb <p> <c> <gt> <f>: get the graph type <gt> of chain <c> of pdbid <p> from the DB and draw it to file <f> (omit the file extension)*");        
        System.out.println("-k | --output-subdir-tree  : write all output files to a PDB-style sudbir tree of the output dir (e.g., <OUTDIR>/ic/8icd/<outfile>). ");                
        System.out.println("-l | --draw-plcc-graph <f> : read graph in plcc format from file <f> and draw it to <f>.png, then exit (<pdbid> will be ignored)*");                
        System.out.println("-L | --lig-filter <i> <a>  : only consider ligands which have at least <i> and at most <a> atoms. A setting of zero means no limit.");                
        System.out.println("-m | --image-format <f>    : write output images in format <f>, which can be 'PNG' or 'JPG' (SVG vector format is always written).");
        System.out.println("-M | --similar <p> <c> <g> : find the proteins which are most similar to pdbid <p> chain <c> graph type <g> in the database.");
        System.out.println("-n | --textfiles           : write meta data, debug info and interim results like residue contacts to text files (slower)");
        System.out.println("-N | --no-warn             : do not print any warnings (intended for cluster use to reduce job output in logs).");
        System.out.println("-o | --outputdir <dir>     : write output files to directory <dir> (instead of '.', the current directory)");
        System.out.println("-O | --outputformats <list>: write only graph output formats in <list>, where g=GML, t=TGF, d=DOT language, e=Kavosh edge list, p=PLCC, j=json. Specify 'x' for none.");
        System.out.println("-p | --pdbfile <pdbfile>   : use input PDB file <pdbfile> (instead of assuming '<pdbid>.pdb')");
        System.out.println("-P | --write-chains-file   : write an info file containing all chain names of the handled PDB file (for other software to know output file names)");
        System.out.println("     --gz-pdbfile <f>      : use gzipped input PDB file <f>.");
        System.out.println("-q | --fg-notations <list> : draw only the folding graph notations in <list>, e.g. 'kars' = KEY, ADJ, RED and SEQ.");
        System.out.println("-r | --recreate-tables     : drop and recreate database tables and add base type data, then exit (see -u). Creates a reday-to-use database.*");
        System.out.println("     --recreate-tables-empty : drop and recreate database tables without adding base type data, then exit (see -u). Creates a completely empty database, suitable for restoring a dump exported from another server.*");        
        System.out.println("-s | --draw-linnots        : not only compute the folding graph linear notations, but draw all 4 of them to image files");
        System.out.println("-S | --sim-measure <m>     : use similarity measure <m>. Valid settings include 'string_sse', 'graph_set' and 'graph_compat'.");
        System.out.println("-t | --draw-tgf-graph <f>  : read graph in TGF format from file <f> and draw it to <f>.png, then exit (pdbid will be ignored)*");
        System.out.println("     --draw-gml-graph <f>  : read graph in GML format from file <f> and draw it to <f>.png, then exit (pdbid will be ignored)*");
        System.out.println("     --props-gml-graph <f> : read graph in GML format from file <f> and compute graph properties, then exit (pdbid will be ignored)*");
        System.out.println("-u | --use-database        : write SSE contact data to database [requires DB credentials in cfg file]");                       
        System.out.println("-v | --del-db-protein <p>  : delete the protein chain with PDBID <p> from the database [requires DB credentials in cfg file]");
        System.out.println("-w | --dont-write-images   : do not draw the SSE graphs and write them to image files [DEBUG]");                             
        System.out.println("-W | --output-www          : add HTML navigation files to the subdir tree (implies -k). ");
        //System.out.println("-x | --check-rescts <f>    : compare the computed residue level contacts to those in geom_neo format file <f> and print differences");
        //System.out.println("-X | --check-ssects <f>    : compare the computed SSE level contacts to those in bet_neo format file <f> and print differences");
        //System.out.println("-y | --write-geodat        : write the computed SSE level contacts in geo.dat format to a file (file name: <pdbid>_<chain>.geodat)");        
        System.out.println("-Y | --skip-vast <atoms>   : abort program execution for PDB files with more than <atoms> atoms before contact computation (for cluster mode, try 80000).");        
        System.out.println("-z | --ramaplot            : draw a ramachandran plot of each chain to the file '<pdbid>_<chain>_plot.svg'");        
        System.out.println("-Z | --silent              : silent mode. do not write output to STDOUT."); 
        System.out.println("     --verbose             : verbose mode. more detailed output."); 
        System.out.println("   --compute-graph-metrics : compute graph metrics like cluster coefficient for PGs. Slower!");
        System.out.println("   --check-whether-in-db   : check whether the PDB file exists in the database and exit. Returns 0 if it does, 1 if not, value >1 on error.");
        System.out.println("   --draw-aag              : visualize amino acid graphs, test only");
        System.out.println("   --force                 : process a PDB file even if the resolution is too bad or the number of residues is too low according to settings.");
        System.out.println("   --cluster               : Set all options for cluster mode. Equals '-f -u -k -s -G -i -Z -P'.");
        System.out.println("   --cg-threshold <Int>    : Overwrites setting for contact thresholds for edges in complex graphs.");
        System.out.println("   --chain-spheres-speedup : speedup for contact computation based on comparison of chain spheres");
        System.out.println("   --include-rna           : Parse RNA and include in graph formalism and visualization");
        System.out.println("   --matrix-structure-search <nt> <ln> <gt>: search a structure <ln> in linear notation in a Proteingraph; <nt> = type of linnot; <gt> = graphtype of linnot");
        System.out.println("   --matrix-structure-search-db <nt> <ln> <gt>: search a structure <ln> in linear notation in the whole database; <nt> = type of linnot; <gt> = graphtype of linnot");
        System.out.println("   --settingsfile <f>      : load settings from file <f>.");
        System.out.println("");
        System.out.println("The following options only make sense for database maintenance:");
        System.out.println("--set-pdb-representative-chains-pre <file> <k> : Set non-redundant chain status for all chains in DB from XML file <file>. <k> determines what to do with existing flags, valid options are 'keep' or 'remove'. Get the file from PDB REST API. Run this pre-update, BEFORE new data will be added.");
        System.out.println("--set-pdb-representative-chains-post <file> <k> : Set non-redundant chain status for all chains already existing in the chains table of the DB from XML file <file>. <k> determines what to do with existing flags, valid options are 'keep' or 'remove'. Get the file from PDB REST API. Run this post-update, after new data has been added.");
        System.out.println("");
        System.out.println("The following options are tools integrated into this software:");
        System.out.println("--convert-models-to-chains <infile> <outfile>: Rewrite the input PDB file, transforming models into chains. Useful when treating PDB files which contain biological assemblies and you want all parts of a protein complex at once, e.g., to analyze the interactions between the chains.");
        System.out.println("");
        
        System.out.println("EXAMPLES: java -jar plcc.jar 8icd");
        System.out.println("          java -jar plcc.jar 8icd -D 2 -d /tmp/dssp/8icd.dssp -p /tmp/pdb/8icd.pdb");
        System.out.println("          java -jar plcc.jar 8icd -o /tmp");
        System.out.println("          java -jar plcc.jar 1o1d -E");
        System.out.println("          java -jar plcc.jar none -l prot_graph_3kmf_A.plg");
        System.out.println("          java -jar plcc.jar none -m PNG -ddb 8icd A albelig ~/img/protein_graph");
        System.out.println("          java -jar plcc.jar 6cbe -I --cg-threshold 2");
        System.out.println("");
        System.out.println("REQUIRED INPUT FILES: This program requires the PDB file and the DSSP file of a protein.");
        System.out.println("                      This does not apply to options that don't use it (marked with * above), of course.");
        System.out.println("                      A PDBID still has to be given as first argument, it will be ignored though (use 'NONE').");
        System.out.println("");
        System.out.println("NOTES: ");
        System.out.println("       -The DSSP program assumes that the input PDB file only has a single model.");
        System.out.println("        You have to split PDB files with multiple models up BEFORE running DSSP (use the 'splitpdb' tool).");
        System.out.println("        If you don't do this, the broken DSSP file will get this program into trouble.");
        System.out.println("       -See the config file '" + Settings.getConfigFile() + "' in your userhome to set advanced options.");
        System.out.println("       -Try 'java -jar plcc.jar <PDBID>' for a start.");
    }

    /**
     * Generates a script for the PyMol protein visualization software that selects all ligands and the residues they are in contact with
     * according to the residue contact calculations of this program.
     * @param contacts the contacts to consider for the script.
     * @return the script as a single string. note that the string may consist of multiple lines.
     */
    public static String getPymolSelectionScript(ArrayList<MolContactInfo> contacts) {

        ArrayList<Residue> protRes = new ArrayList<Residue>();
        ArrayList<Residue> ligRes = new ArrayList<Residue>();

        String script = "";
        String scriptProt = "";
        String scriptLig = "";

        MolContactInfo c = null;
        // Select all residues of the protein that have ligand contacts
        
        for (Integer i = 0; i < contacts.size(); i++) {
            c = contacts.get(i);
            if(c.getNumLigContactsTotal() > 0) {
                // This is a ligand contact, add the residues to one of the lists depending on their type (ligand or protein residue)
                // c getMolA und c getMolB -> typecast
                // Handle resA
                
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                    
                if(c.getMolA().isAA()) {
                    if( ! protRes.contains(c.getMolA())) {
                        protRes.add((Residue)c.getMolA());
                    }
                }
                if(c.getMolA().isLigand()) {
                    if( ! ligRes.contains(c.getMolA())) {
                        ligRes.add((Residue)c.getMolA());
                    }
                }

                // Handle resB
                if(c.getMolB().isAA()) {
                    if( ! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                }
                if(c.getMolB().isLigand()) {
                    if( ! ligRes.contains(c.getMolB())) {
                        ligRes.add((Residue)c.getMolB());
                    }
                }

            }
        }
        }

        // Now create the scripts...

        // Handle protein residues
        if(protRes.isEmpty()) {
            scriptProt = "select none";
        } else {

            scriptProt = "select contact_res, resi ";
            for(Integer j = 0; j < protRes.size(); j++) {
                scriptProt += protRes.get(j).getPdbNum();

                if(j < (protRes.size() - 1)) {
                    scriptProt += "+";
                }
            }
        }

        // Handle ligands
        if(ligRes.isEmpty()) {
            scriptLig = "select none";
        } else {

            scriptLig = "";
            for(Integer k = 0; k < ligRes.size(); k++) {

                scriptLig += "select lig_" + ligRes.get(k).getName3().trim() + ", resi " + ligRes.get(k).getPdbNum() + "\n";

            }
        }


        script = scriptProt + "\n" + scriptLig;
        return(script);
    }
    
    /**
     * Generates a script for the PyMol protein visualization software that selects all ligands and the residues they are in contact with
     * according to the residue contact calculations of this program.
     * @param protRes the protein residues
     * @return the script as a single string. note that the string may consist of multiple lines.
     */
    public static String getPymolSelectionScriptForResidues(List<Residue> protRes) {


        StringBuilder scriptProt = new StringBuilder();

        
        // Handle protein residues
        if(protRes.isEmpty()) {
            scriptProt.append("select none");
        } else {

            scriptProt.append("select prot_res, resi ");
            for(Integer j = 0; j < protRes.size(); j++) {
                scriptProt.append(protRes.get(j).getPdbNum());

                if(j < (protRes.size() - 1)) {
                    scriptProt.append("+");
                }
            }
        }       

        return(scriptProt.toString());
    }


    /**
     * Creates the PyMOL script, but with a separate selection of contact residues for each ligand.
     * @param contacts the contacts to consider for the script.
     * @return the PyMol script as a string, which may consist of more than one line.
     */
    public static String getPymolSelectionScriptByLigand(ArrayList<MolContactInfo> contacts) {
        
        //the program assumes that we are working with object of class Molecule. 
        //However, there are methods that just need a residue as input and/or output parameter, 
        //so we always have to query which variables belong to which instance and do a typecast.

        ArrayList<Residue> protRes = new ArrayList<Residue>();
        ArrayList<Residue> ligRes = new ArrayList<Residue>();
        ArrayList<Residue> ligCont = new ArrayList<Residue>();
     
        String scriptLig = "";
        String scriptThisLigCont = "";

        MolContactInfo c = null;
        Residue r = null;
        // Select all residues of the protein that have ligand contacts

        for (Integer i = 0; i < contacts.size(); i++) {
            c = contacts.get(i);
            if(c.getNumLigContactsTotal() > 0) {
                // This is a ligand contact, add the residues to one of the lists depending on their type (ligand or protein residue)
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                    // Handle resA
                if(c.getMolA().isAA()) {
                    if( ! protRes.contains(c.getMolA())) {
                        protRes.add((Residue)c.getMolA());
                    }
                }
                if(c.getMolA().isLigand()) {
                    if( ! ligRes.contains(c.getMolA())) {
                        ligRes.add((Residue)c.getMolA());
                    }
                 
                }
                
                // Handle resB
                if(c.getMolB().isAA()) {
                    if( ! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                }
                if(c.getMolB().isLigand()) {
                    if( ! ligRes.contains(c.getMolB())) {
                        ligRes.add((Residue)c.getMolB());
                    }
                }

               }

                

            }
        }

        // Now create the scripts...

        
        // Handle ligands
        if(ligRes.isEmpty()) {
            // None of the residue contacts include ligands.
            return("select none");
        } else {

            scriptLig = "";
            for(Integer k = 0; k < ligRes.size(); k++) {

                r = ligRes.get(k);
                scriptLig += "select lig_" + r.getName3().trim() + r.getPdbNum() + ", chain " + r.getChainID() + " and resi " + r.getPdbNum() + "\n";

                // create the list of contact residues for this ligand
                ligCont = new ArrayList<Residue>();
                for(Integer j = 0; j < contacts.size(); j++) {

                    c = contacts.get(j);
                    
                    if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                    }
                    else{
                        if(c.getDsspNumA().equals(r.getDsspNum())) {
                        
                        
                        // first residue A is this ligand, so the other one is the contact residue
                        ligCont.add((Residue)c.getMolB());
                    }
                    else if(c.getDsspNumB().equals(r.getDsspNum())) {
                        // second residue B is this ligand, so the other one is the contact residue
                        ligCont.add((Residue)c.getMolA());
                    }
                    else {
                        // The current ligand is not involved in this contact
                    } 
                    }

                   
                }

                
                // Now create the PyMOL script of the contact residues of this ligand
                if(ligCont.isEmpty()) {
                    // The list of contact residues should not be empty at this time.
                    scriptThisLigCont = "";
                    DP.getInstance().w("getPymolSelectionScriptByLigand(): Residue without contacts in contact list. Bug?");
                } else {

                    scriptThisLigCont = "select lig_" + r.getName3().trim() + r.getPdbNum() + "_contacts,";

                    for(Integer j = 0; j < ligCont.size(); j++) {
                        scriptThisLigCont += " (resi " + ligCont.get(j).getPdbNum() + " and chain " + ligCont.get(j).getChainID() + ")";

                        if(j < (ligCont.size() - 1)) {
                            scriptThisLigCont += " or";
                        }
                    }

                    scriptLig += scriptThisLigCont + "\n";
                }

                // That's it for this ligand.
                // System.out.println("Ligand " + r.getFancyName() + " (chainName " + r.getChainID() + ") has " + ligCont.size() + " contacts on residue level.");

            }
        }

        return(scriptLig);
    }

    /**
     * Creates a Python script for PyMol to visualize bonds calculated with the alternative contact model.
     * @param contacts the contacts to consider for this script.
     * @return true if python file could be written, otherwise false.
     */
    public static Boolean getPymolSelectionScriptPPI (ArrayList<MolContactInfo> contacts, String pdbid) {
        ArrayList<Residue> protRes = new ArrayList<Residue>();  // all residues of interchain protein contacts
        ArrayList<Residue> ligRes = new ArrayList<Residue>();   // all residues of ligand contacts
        ArrayList<Residue> ivdwRes = new ArrayList<Residue>();  // all residues of interchain van der Waals contacts
        ArrayList<Residue> issRes = new ArrayList<Residue>();   // all residues of interchain sulfur bridge contacts
        ArrayList<Residue> bbRes = new ArrayList<Residue>();    // all residues of interchain backbone-backbone h-bridge contacts
        ArrayList<Residue> bcRes = new ArrayList<Residue>();    // all residues of interchain backbone-sidechain h-bridge contacts
        ArrayList<Residue> cbRes = new ArrayList<Residue>();      // all residues of interchain sidechain-backbone h-bridge contacts
        ArrayList<Residue> ccRes = new ArrayList<Residue>();    // all residues of interchain sidechain-sidechain h-bridge contacts
        
        String script = "";
        String scriptProt = "";
        String scriptLig = "";
        String scriptIvdw = "";
        String scriptIss = "";
        String scriptBB = "";
        String scriptBC = "";
        String scriptCB = "";
        String scriptCC = "";
        
        ArrayList<ArrayList<Integer>> bondsIvdw = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainIvdw = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsIvdwAtoms = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsAtomNamesIvdw = new ArrayList<ArrayList<String>>();
        
        ArrayList<ArrayList<Integer>> bondsIss = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainIss = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsIssAtoms = new ArrayList<ArrayList<Integer>>();
        
        ArrayList<ArrayList<Integer>> bondsBB = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainBB = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsBBAtoms = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsAtomNamesBB = new ArrayList<ArrayList<String>>();
        
        ArrayList<ArrayList<Integer>> bondsBC = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainBC = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsBCAtoms = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsAtomNamesBC = new ArrayList<ArrayList<String>>();
        
        ArrayList<ArrayList<Integer>> bondsCB = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainCB = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsCBAtoms = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsAtomNamesCB = new ArrayList<ArrayList<String>>();
        
        ArrayList<ArrayList<Integer>> bondsCC = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsChainCC = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> bondsCCAtoms = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<String>> bondsAtomNamesCC = new ArrayList<ArrayList<String>>();
        
        Boolean fileWriteSuccess = false;
        
        MolContactInfo c = null;
        // Select all residues of the protein that have ligand contacts
        for (Integer i = 0; i < contacts.size(); i++) {
            c = contacts.get(i);
            
            if(c.getNumLigContactsTotal() > 0) {
                // This is a ligand contact, add the residues to one of the lists depending on their type (ligand or protein residue)
                
                // Handle resA
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                    
                    
                    if(c.getMolA().isAA()) {
                    if( ! protRes.contains(c.getMolA())) {
                        protRes.add((Residue)c.getMolA());
                    }
                }
                if(c.getMolA().isLigand()) {
                    if( ! ligRes.contains(c.getMolA())) {
                        ligRes.add((Residue)c.getMolA());
                    }
                }
                
                 // Handle resB
                if(c.getMolB().isAA()) {
                    if( ! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                }
                if(c.getMolB().isLigand()) {
                    if( ! ligRes.contains(c.getMolB())) {
                        ligRes.add((Residue)c.getMolB());
                    }
                }
                }
                
            }
            
            // Select all residues of the protein that are protein-protein contacts
            
            if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                if((c.getNumContactsBB() > 0) || (c.getNumContactsBC() > 0) || (c.getNumContactsCB() > 0)|| (c.getNumContactsCC() > 0)) {
                if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                }
                if(! protRes.contains(c.getMolB())) {
                    protRes.add((Residue)c.getMolB());
                }
            }
            }
            
            
            // Select all residues of the protein that have interchain vdW contacts
            if(c.getNumContactsIVDW() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                  if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                }   
                if(! protRes.contains(c.getMolB())) {
                    protRes.add((Residue)c.getMolB());
                }
                if(! ivdwRes.contains(c.getMolA())) {
                    ivdwRes.add((Residue)c.getMolA());
                }
                if(! ivdwRes.contains(c.getMolB())) {
                    ivdwRes.add((Residue)c.getMolB());
                }  
                
                ArrayList<Integer> tmp = new ArrayList<Integer>();
                tmp.add(c.getPdbNumA());
                tmp.add(c.getPdbNumB());
                bondsIvdw.add(tmp);
                
                ArrayList<String> chains = new ArrayList<String>();
                chains.add(c.getMolA().getChainID());
                chains.add(c.getMolB().getChainID());
                bondsChainIvdw.add(chains);
                
                ArrayList<Integer> atoms = new ArrayList<Integer>();
                atoms.add(c.getIVDWContactAtomNumA() - 1);
                atoms.add(c.getIVDWContactAtomNumB() - 1);
                bondsIvdwAtoms.add(atoms);
                
                ArrayList<String> atomNames = new ArrayList<String>();
                atomNames.add(c.getMolA().getAtoms().get(c.getIVDWContactAtomNumA() - 1).getAtomName());
                atomNames.add(c.getMolB().getAtoms().get(c.getIVDWContactAtomNumB() - 1).getAtomName());
                bondsAtomNamesIvdw.add(atomNames);
                }
                 
                
            }
            
            // Select all residues of the protein that have interchain backbone-backbone h-bridge contacts
            if(c.getNumContactsBBHB()> 0 || c.getNumContactsBBBH() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                  if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                    }
                    if(! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                    if(! bbRes.contains(c.getMolA())) {
                        bbRes.add((Residue)c.getMolA());
                    }
                    if(! bbRes.contains(c.getMolB())) {
                        bbRes.add((Residue)c.getMolB());
                    } 
                    
                    ArrayList<Integer> tmp = new ArrayList<Integer>();
                tmp.add(c.getPdbNumA());
                tmp.add(c.getPdbNumB());
                bondsBB.add(tmp);
                
                ArrayList<String> chains = new ArrayList<String>();
                chains.add(c.getMolA().getChainID());
                chains.add(c.getMolB().getChainID());
                bondsChainBB.add(chains);
                
                if(c.getNumContactsBBHB() > 0) {
                    ArrayList<Integer> atoms = new ArrayList<Integer>();
                    atoms.add(c.getBBHBContactAtomNumA() - 1);
                    atoms.add(c.getBBHBContactAtomNumB() - 1);
                    bondsBBAtoms.add(atoms);
                
                ArrayList<String> atomNames = new ArrayList<String>();
                atomNames.add(c.getMolA().getAtoms().get(c.getBBHBContactAtomNumA() - 1).getAtomName());
                atomNames.add(c.getMolB().getAtoms().get(c.getBBHBContactAtomNumB() - 1).getAtomName());
                bondsAtomNamesBB.add(atomNames);
                }
                
                if(c.getNumContactsBBBH() > 0) {
                    ArrayList<Integer> atoms = new ArrayList<Integer>();
                    atoms.add(c.getBBBHContactAtomNumA() - 1);
                    atoms.add(c.getBBBHContactAtomNumB() - 1);
                    bondsBBAtoms.add(atoms);
                
                ArrayList<String> atomNames = new ArrayList<String>();
                atomNames.add(c.getMolA().getAtoms().get(c.getBBBHContactAtomNumA() - 1).getAtomName());
                atomNames.add(c.getMolB().getAtoms().get(c.getBBBHContactAtomNumB() - 1).getAtomName());
                bondsAtomNamesBB.add(atomNames);
                }
                }
                
                
                
                
                
            }
            
            // Select all residues of the protein that have interchain backbone-sidechain h-bridge contacts
            if(c.getNumContactsBCHB()> 0 || c.getNumContactsBCBH() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                    if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                    }
                    if(! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                    if(! bcRes.contains(c.getMolA())) {
                        bcRes.add((Residue)c.getMolA());
                    }
                    if(! bcRes.contains(c.getMolB())) {
                        bcRes.add((Residue)c.getMolB());
                    }

                    ArrayList<Integer> tmp = new ArrayList<Integer>();
                    tmp.add(c.getPdbNumA());
                    tmp.add(c.getPdbNumB());
                    bondsBC.add(tmp);

                    ArrayList<String> chains = new ArrayList<String>();
                    chains.add(c.getMolA().getChainID());
                    chains.add(c.getMolB().getChainID());
                    bondsChainBC.add(chains);

                    if(c.getNumContactsBCHB() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getBCHBContactAtomNumA() - 1);
                        atoms.add(c.getBCHBContactAtomNumB() - 1);
                        bondsBCAtoms.add(atoms);
                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getBCHBContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getBCHBContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesBC.add(atomNames);

                    }

                    if(c.getNumContactsBCBH() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getBCBHContactAtomNumA() - 1);
                        atoms.add(c.getBCBHContactAtomNumB() - 1);
                        bondsBCAtoms.add(atoms);
                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getBCBHContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getBCBHContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesBC.add(atomNames);
                    }
                }   
            }
            
            // Select all residues of the protein that have interchain sidechain-backbone h-bridge contacts
            if(c.getNumContactsCBHB()> 0 || c.getNumContactsCBBH() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                   if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                    }
                    if(! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                    if(! cbRes.contains(c.getMolA())) {
                        cbRes.add((Residue)c.getMolA());
                    }
                    if(! cbRes.contains(c.getMolB())) {
                        cbRes.add((Residue)c.getMolB());
                    }

                    ArrayList<Integer> tmp = new ArrayList<Integer>();
                    tmp.add(c.getPdbNumA());
                    tmp.add(c.getPdbNumB());
                    bondsCB.add(tmp);

                    ArrayList<String> chains = new ArrayList<String>();
                    chains.add(c.getMolA().getChainID());
                    chains.add(c.getMolB().getChainID());
                    bondsChainCB.add(chains);

                    if(c.getNumContactsCBHB() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getCBHBContactAtomNumA() - 1);
                        atoms.add(c.getCBHBContactAtomNumB() - 1);
                        bondsCBAtoms.add(atoms);

                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getCBHBContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getCBHBContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesCB.add(atomNames);
                    }

                    if(c.getNumContactsCBBH() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getCBBHContactAtomNumA() - 1);
                        atoms.add(c.getCBBHContactAtomNumB() - 1);
                        bondsCBAtoms.add(atoms);

                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getCBBHContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getCBBHContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesCB.add(atomNames);
                    } 
                }
            }
            
            // Select all residues of the protein that have interchain sidechain-sidechain h-bridge contacts
            if(c.getNumContactsCCHB()> 0 || c.getNumContactsCCBH() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                    if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                    }
                    if(! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                    if(! ccRes.contains(c.getMolA())) {
                        ccRes.add((Residue)c.getMolA());
                    }
                    if(! ccRes.contains(c.getMolB())) {
                        ccRes.add((Residue)c.getMolB());
                    }

                    ArrayList<Integer> tmp = new ArrayList<Integer>();
                    tmp.add(c.getPdbNumA());
                    tmp.add(c.getPdbNumB());
                    bondsCC.add(tmp);

                    ArrayList<String> chains = new ArrayList<String>();
                    chains.add(c.getMolA().getChainID());
                    chains.add(c.getMolB().getChainID());
                    bondsChainCC.add(chains);

                    if(c.getNumContactsCCHB() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getCCHBContactAtomNumA() - 1);
                        atoms.add(c.getCCHBContactAtomNumB() - 1);
                        bondsCCAtoms.add(atoms);

                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getCCHBContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getCCHBContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesCC.add(atomNames);
                    }

                    if(c.getNumContactsCCBH() > 0) {
                        ArrayList<Integer> atoms = new ArrayList<Integer>();
                        atoms.add(c.getCCBHContactAtomNumA() - 1);
                        atoms.add(c.getCCBHContactAtomNumB() - 1);
                        bondsCCAtoms.add(atoms);

                        ArrayList<String> atomNames = new ArrayList<String>();
                        atomNames.add(c.getMolA().getAtoms().get(c.getCCBHContactAtomNumA() - 1).getAtomName());
                        atomNames.add(c.getMolB().getAtoms().get(c.getCCBHContactAtomNumB() - 1).getAtomName());
                        bondsAtomNamesCC.add(atomNames);
                    }
                }    
            }
            
            // Select all residues of the protein that have interchain sulfur contacts
            if(c.getNumContactsISS() > 0) {
                if( !(c.getMolA() instanceof Residue) || !(c.getMolB() instanceof Residue)){
                    continue;
                }
                else{
                   if(! protRes.contains(c.getMolA())) {
                    protRes.add((Residue)c.getMolA());
                    }
                    if(! protRes.contains(c.getMolB())) {
                        protRes.add((Residue)c.getMolB());
                    }
                    if(! issRes.contains(c.getMolA())) {
                        issRes.add((Residue)c.getMolA());
                    }
                    if(! issRes.contains(c.getMolB())) {
                        issRes.add((Residue)c.getMolB());
                    }

                    ArrayList<Integer> tmp = new ArrayList<Integer>();
                    tmp.add(c.getPdbNumA());
                    tmp.add(c.getPdbNumB());
                    bondsIss.add(tmp);

                    ArrayList<String> chains = new ArrayList<String>();
                    chains.add(c.getMolA().getChainID());
                    chains.add(c.getMolB().getChainID());
                    bondsChainIss.add(chains); 
                }
            }
        }

        // Now create the scripts...

        // Handle protein residues
        if(protRes.isEmpty()) {
            scriptProt = "select none";
        } else {

            scriptProt = "select protein_contact_res, ";
            for(Integer j = 0; j < protRes.size(); j++) {
                scriptProt += "chain " + protRes.get(j).getChainID() + " and resi " + protRes.get(j).getPdbNum();

                if(j < (protRes.size() - 1)) {
                    scriptProt += " + ";
                }
            }
        }

        // Handle ligands
        if(ligRes.isEmpty()) {
            scriptLig = "select none";
        } else {

            scriptLig = "select lig_contact_res, ";
            for(Integer k = 0; k < ligRes.size(); k++) {
                scriptLig += "chain " + ligRes.get(k).getChainID() + " and resi " + ligRes.get(k).getPdbNum();
                
                if(k < (ligRes.size() - 1)) {
                    scriptLig += " + ";
                }

            }
        }

        if(ivdwRes.isEmpty()) {
            scriptIvdw = "select none";
        } else {
            scriptIvdw = "select ivdw_contact_res, ";
            for(Integer i = 0; i < ivdwRes.size(); i++) {
                scriptIvdw += "chain " + ivdwRes.get(i).getChainID() + " and resi " + ivdwRes.get(i).getPdbNum();
                
                if(i < (ivdwRes.size() - 1)) {
                    scriptIvdw += " + ";
                }
            }
        }
        
        if(bbRes.isEmpty()) {
            scriptBB = "select none";
        } else {
            scriptBB = "select bb_h_bridge_contact_res, ";
            for(Integer i = 0; i < bbRes.size(); i++) {
                scriptBB += "chain " + bbRes.get(i).getChainID() + " and resi " + bbRes.get(i).getPdbNum();
                
                if(i < (bbRes.size() - 1)) {
                    scriptBB += " + ";
                }
            }
        }
        
        if(bcRes.isEmpty()) {
            scriptBC = "select none";
        } else {
            scriptBC = "select bc_h_bridge_contact_res, ";
            for(Integer i = 0; i < bcRes.size(); i++) {
                scriptBC += "chain " + bcRes.get(i).getChainID() + " and resi " + bcRes.get(i).getPdbNum();
                
                if(i < (bcRes.size() - 1)) {
                    scriptBC += " + ";
                }
            }
        }
        
        if(cbRes.isEmpty()) {
            scriptCB = "select none";
        } else {
            scriptCB = "select cb_h_bridge_contact_res, ";
            for(Integer i = 0; i < cbRes.size(); i++) {
                scriptCB += "chain " + cbRes.get(i).getChainID() + " and resi " + cbRes.get(i).getPdbNum();
                
                if(i < (cbRes.size() - 1)) {
                    scriptCB += " + ";
                }
            }
        }
        
        if(ccRes.isEmpty()) {
            scriptCC = "select none";
        } else {
            scriptCC = "select cc_h_bridge_contact_res, ";
            for(Integer i = 0; i < ccRes.size(); i++) {
                scriptCC += "chain " + ccRes.get(i).getChainID() + " and resi " + ccRes.get(i).getPdbNum();
                
                if(i < (ccRes.size() - 1)) {
                    scriptCC += " + ";
                }
            }
        }
        
        if(issRes.isEmpty()) {
            scriptIss = "select none";
        } else {
            scriptIss = "select iss_contact_res, ";
            
            for(Integer x = 0; x < issRes.size(); x++) {
                scriptIss += "chain " + issRes.get(x).getChainID() + " and resi " + issRes.get(x).getPdbNum();
                
                if(x < (issRes.size()) - 1) {
                    scriptIss += " + ";
                }
            }
        }
        
        // Create python file that can be executed from pymol
        StringBuilder sb = new StringBuilder();
        String lineSep = System.lineSeparator();
        String blankLine = lineSep + lineSep;
        
        sb.append("#!/usr/bin/env python");
        sb.append(blankLine);
        
        sb.append("# This file is part of the Visualization of Protein Ligand Graphs (VPLG) software package.");
        sb.append(lineSep);
        sb.append("# *");
        sb.append(lineSep);
        sb.append("# * Copyright Tim Schaefer 2012. VPLG is free software, see the LICENSE and README files for details");
        sb.append(lineSep);
        sb.append("# * author Andreas Scheck");
        sb.append(lineSep);
        sb.append("# *");
        sb.append(blankLine);
        sb.append("\"\"\"");
        sb.append(lineSep);
        sb.append("    This is a small script for visualizing different types of chemical bonds, calculated with the VPLG software, in PyMol.");
        sb.append(lineSep);
        sb.append("    Fetch the protein you want to inspect in PyMol and then run this python script, created for that specific protein,");
        sb.append(lineSep);
        sb.append("    within PyMol (File -> Run... -> /path/to/this/file).");
        sb.append(lineSep);
        sb.append("\"\"\"");
        sb.append(blankLine);
        
        sb.append("import pymol");
        sb.append(blankLine);
        sb.append("select_ligands = \'" + scriptLig + "\'");
        sb.append(blankLine);
        sb.append("select_prot = \'" + scriptProt + "\'");
        sb.append(blankLine);
        sb.append("select_ivdw = \'" + scriptIvdw + "\'");
        sb.append(blankLine);
        sb.append("select_bb_h_bridge = \'" + scriptBB + "\'");
        sb.append(blankLine);
        sb.append("select_bc_h_bridge = \'" + scriptBC + "\'");
        sb.append(blankLine);
        sb.append("select_cb_h_bridge = \'" + scriptCB + "\'");
        sb.append(blankLine);
        sb.append("select_cc_h_bridge = \'" + scriptCC + "\'");
        sb.append(blankLine);
        sb.append("select_iss = \'" + scriptIss + "\'");
        sb.append(blankLine);
        sb.append("bonds_ivdw = " + bondsIvdw.toString());
        sb.append(blankLine);
        sb.append("bonds_bb = " + bondsBB.toString());
        sb.append(blankLine);
        sb.append("bonds_bc = " + bondsBC.toString());
        sb.append(blankLine);
        sb.append("bonds_cb = " + bondsCB.toString());
        sb.append(blankLine);
        sb.append("bonds_cc = " + bondsCC.toString());
        sb.append(blankLine);
        
        sb.append("bonds_chain_ivdw = [");
        for(ArrayList<String> valuePair : bondsChainIvdw) {
                sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_chain_bb = [");
        for(ArrayList<String> valuePair : bondsChainBB) {
                sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_chain_bc = [");
        for(ArrayList<String> valuePair : bondsChainBC) {
                sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_chain_cb = [");
        for(ArrayList<String> valuePair : bondsChainCB) {
                sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_chain_cc = [");
        for(ArrayList<String> valuePair : bondsChainCC) {
                sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_iss = " + bondsIss.toString());
        sb.append(blankLine);
        
        sb.append("bonds_chain_iss = [");
        for(ArrayList<String> valuePair : bondsChainIss) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_atom_names_ivdw = [");
        for(ArrayList<String> valuePair : bondsAtomNamesIvdw) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_atom_names_bb = [");
        for(ArrayList<String> valuePair : bondsAtomNamesBB) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_atom_names_bc = [");
        for(ArrayList<String> valuePair : bondsAtomNamesBC) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_atom_names_cb = [");
        for(ArrayList<String> valuePair : bondsAtomNamesCB) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("bonds_atom_names_cc = [");
        for(ArrayList<String> valuePair : bondsAtomNamesCC) {
            sb.append("[\'" + valuePair.get(0) + "\',\'" + valuePair.get(1) + "\'],");
        }
        sb.append("]");
        sb.append(blankLine);
        
        sb.append("pymol.cmd.do(\'{}\'.format(select_ivdw))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_iss))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_bb_h_bridge))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_bc_h_bridge))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_cb_h_bridge))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_cc_h_bridge))");
        sb.append(blankLine);
        sb.append("pymol.cmd.do(\'{}\'.format(select_ligands))");
        sb.append(blankLine);
        
        sb.append("for x, y in enumerate(bonds_ivdw):");
        sb.append(lineSep);
        sb.append("    pymol.cmd.distance(\'ivdw_distance\', \'chain {} and resi {} and name {}\'.format(bonds_chain_ivdw[x][0], y[0], bonds_atom_names_ivdw[x][0]), 'chain {} and resi {} and name {}'.format(bonds_chain_ivdw[x][1], y[1], bonds_atom_names_ivdw[x][1]))");
        sb.append(lineSep);
        sb.append("    pymol.cmd.color(\'red\', \'ivdw_distance\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_width, 7\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_length, 0.5\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_gap, 0.2\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set label_distance_digits, 2\')");
        sb.append(blankLine);
        
        sb.append("for x, y in enumerate(bonds_bb):");
        sb.append(lineSep);
        sb.append("    pymol.cmd.distance(\'bb_distance\', \'chain {} and resi {} and name {}\'.format(bonds_chain_bb[x][0], y[0], bonds_atom_names_bb[x][0]), 'chain {} and resi {} and name {}'.format(bonds_chain_bb[x][1], y[1], bonds_atom_names_bb[x][1]))");
        sb.append(lineSep);
        sb.append("    pymol.cmd.color(\'blue\', \'bb_distance\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_width, 7\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_length, 0.5\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_gap, 0.2\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set label_distance_digits, 2\')");
        sb.append(blankLine);
        
        sb.append("for x, y in enumerate(bonds_bc):");
        sb.append(lineSep);
        sb.append("    pymol.cmd.distance(\'bc_distance\', \'chain {} and resi {} and name {}\'.format(bonds_chain_bc[x][0], y[0], bonds_atom_names_bc[x][0]), 'chain {} and resi {} and name {}'.format(bonds_chain_bc[x][1], y[1], bonds_atom_names_bc[x][1]))");
        sb.append(lineSep);
        sb.append("    pymol.cmd.color(\'magenta\', \'bc_distance\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_width, 7\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_length, 0.5\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_gap, 0.2\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set label_distance_digits, 2\')");
        sb.append(blankLine);
        
        sb.append("for x, y in enumerate(bonds_cb):");
        sb.append(lineSep);
        sb.append("    pymol.cmd.distance(\'cb_distance\', \'chain {} and resi {} and name {}\'.format(bonds_chain_cb[x][0], y[0], bonds_atom_names_cb[x][0]), 'chain {} and resi {} and name {}'.format(bonds_chain_cb[x][1], y[1], bonds_atom_names_cb[x][1]))");
        sb.append(lineSep);
        sb.append("    pymol.cmd.color(\'violet\', \'cb_distance\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_width, 7\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_length, 0.5\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_gap, 0.2\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set label_distance_digits, 2\')");
        sb.append(blankLine);
        
        sb.append("for x, y in enumerate(bonds_cc):");
        sb.append(lineSep);
        sb.append("    pymol.cmd.distance(\'cc_distance\', \'chain {} and resi {} and name {}\'.format(bonds_chain_cc[x][0], y[0], bonds_atom_names_cc[x][0]), 'chain {} and resi {} and name {}'.format(bonds_chain_cc[x][1], y[1], bonds_atom_names_cc[x][1]))");
        sb.append(lineSep);
        sb.append("    pymol.cmd.color(\'cyan\', \'cc_distance\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_width, 7\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_length, 0.5\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set dash_gap, 0.2\')");
        sb.append(lineSep);
        sb.append("pymol.cmd.do(\'set label_distance_digits, 2\')");
        sb.append(blankLine);
        
        File pythonScript = new File("./" + pdbid + "_visualize_bonds_pymol.py");
        try {
            FileWriter fw = new FileWriter(pythonScript.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            bw.append(sb);
            bw.close();
            fileWriteSuccess = true;
        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        

        DP.getInstance().i("main", "Writing python script to visualize bonds within PyMol.");
        
        //System.out.println("Pymol selection scripts: protein, ligand, ivdw, iss\n");
        //script = scriptProt + "\n" + scriptLig + "\n" + scriptIvdw + "\n" + scriptIss + "\n" + bondsIvdw.toString() + "\n" + bondsChainIvdw.toString() + "\n" + bondsAtomNamesIvdw.toString() + "\n" + bondsIss.toString() + "\n" + bondsChainIss.toString();
        return(fileWriteSuccess);
    }
    
    
    /**
     * Creates all SSEs according to DSSP definition from a list of residues. The list is expected to be ordered like in the DSSP file. (Note that residues which are no alpha helix or beta strand are not assigned to any SSE default. You can change this in the settings by settings, see the 'plcc_B_include_coils' option.)
     * @param resList the residue list to consider
     * @return a list of secondary structure elements (SSEs)
     */
    private static ArrayList<SSE> createAllDsspSSEsFromResidueList(ArrayList<Residue> resList) {

        ArrayList<SSE> dsspSSElist = new ArrayList<SSE>();
        String lastResString = "some initial value that is not a valid residue string";
        String curResString = "";           // Doesn't matter, will be overwritten before 1st comparison
        SSE curSSE, lastSSE;
        curSSE = lastSSE = null;

        if(resList.size() < 1) {
            DP.getInstance().w("Main", "createAllDsspSSEsFromResidueList(): Creating empty list of SSEs: residue list is empty.");
            //System.out.println("      Found " + dsspSSElist.size() + " SSEs according to DSSP definition.");
            return(dsspSSElist);
        }

        // Note that the list of residues is ORDERED in dssp order because the residues were read from the DSSP file.
        Residue curResidue, lastResidue;
        curResidue = lastResidue = null;
        String coil = Settings.get("plcc_S_coilSSECode");

        for(Integer i = 0; i < resList.size(); i++) {

            curResidue = resList.get(i);
            curResString = curResidue.getSSEString();
            curResidue.setSSEStringDssp(curResString);

            // If coiled regions should be included, this line keeps them from being ignored below.
            // In this case, we also assign the SSE type "coil" to all SSEs which would otherwise be
            // ignored later by the getImportantSSEs() filter function. This way, each residue of the
            // protein is assigned some SSE type (H, G, L or C) and all residues appear in the graph.
            if(Settings.getBoolean("plcc_B_include_coils")) {

                if(curResString.equals(" ") || curResString.equals("B") || curResString.equals("S") || curResString.equals("T") || curResString.equals("I") || curResString.equals("G")) {
                    curResString = coil;
                    curResidue.setSSEString(coil);
                }
                
            }

            //System.out.println("   *At DSSP residue " + curResidue.getDsspNum() + ", PDB name is " + curResidue.getFancyName() + ", SSE string is '" + curResString + "'.");

            if( ! curResString.equals(lastResString)) {

                // The SSE string is different so the old SSE has ended.
                //System.out.println("    New SSE starts at residue " + curResidue.getFancyName() + " (may be invalid SSE though).");

                

                // Create a new SSE for the residue if this residue is part of any SSE (DSSP does NOT assign an SSE to all residues)
                if( ! curResidue.getSSEString().equals(" ")) {

                    //System.out.println("      Residue " + curResidue.getFancyName() + " is part of a valid SSE of type '" + curResidue.getSSEString() + "'.");
                    // This is a valid SSE according to DSSP so create it
                    curSSE = new SSE(curResidue.getSSEString());

                    //System.out.println("      Created new SSE.");

                    // set SSE properties
                    curSSE.addResidue(curResidue);
                    curSSE.setSeqSseNumDssp(dsspSSElist.size() + 1);
                    curSSE.setSseType(curResidue.getSSEString());

                    // set Residue properties
                    curResidue.setSSE(curSSE);
                    curResidue.setDsspSseState(true);


                    // Add the new SSE to the list of SSEs. Ignore ligands though, they are handled later because they are special.
                    if(! curSSE.getSseType().equals(SSE.SSECLASS_STRING_LIGAND)) {
                        dsspSSElist.add(curSSE);
                        //System.out.println("      Added new SSE of type '" + curSSE.getSseType() + "', " + s_allDsspSSEs.size() + " SSEs found so far.");
                    }

                }
                else {

                    // This residue is not part of any SSE according to the DSSP definition
                    curResidue.setDsspSseState(false);
                    //System.out.println("      Residue " + curResidue.getFancyName() + " is NOT part of any valid SSE.");
                }

            }
            else {
                // We're still in the old SSE.
                //System.out.println("    Residue " + curResidue.getFancyName() + " is NOT the start of a new SSE.");

                // Set the proper SSE for the residue if this residue is part of any SSE (DSSP does NOT assign an SSE to all residues)
                if( ! curResidue.getSSEString().equals(" ")) {

                    // This is a valid SSE according to DSSP
                    curResidue.setSSE(curSSE);
                    curResidue.setDsspSseState(true);
                    curSSE.addResidue(curResidue);
                    //System.out.println("      Residue " + curResidue.getFancyName() + " added to an existing SSE.");
                }
                else {

                    // This residue is not part of any SSE according to the DSSP definition
                    curResidue.setDsspSseState(false);
                    //System.out.println("      Residue " + curResidue.getFancyName() + " is NOT part of any valid SSE (neither was the last residue).");
                }
            }


            // update for next iteration of loop
            lastResString = curResString;
        }

        //System.out.println("      Found " + dsspSSElist.size() + " SSEs according to DSSP definition.");
        return(dsspSSElist);

    }

    
    /**
     * Returns a list containing all residues from resList that are part of chainName cID.
     * @param resList the residue list
     * @param cID the chainName ID
     * @return the filtered list of residues
     */
    public static ArrayList<Residue> filterResidueListByChain(ArrayList<Residue> resList, String cID) {
        ArrayList<Residue> filteredList = new ArrayList<Residue>();
        Residue r;

        for(Integer i = 0; i < resList.size(); i++) {
            r = resList.get(i);
            if(r.getChainID().equals(cID)) {
                filteredList.add(r);
            }
        }
        
        return(filteredList);
    }


    /**
     * Creates all the ligand SSEs and adds them to the output list of ligand SSEs (the input list is not changed, it is only required to determine the SSE numbers for the new ligand SSEs). Assumes that the ligands
     * are already in the list of residues (i.e., createAllLigandResiduesFromPdbData() has already been called).
     * @param resList the residue list
     * @param dsspSSElist the list of all SSEs according to DSSP definition
     * @return a list of all ligand SSEs that are created from the ligand residues in the residue list.
     */
    private static List<SSE> createAllLigandSSEsFromResidueList(List<Residue> resList, List<SSE> dsspSSElist) {

        List<SSE> ligSSElist = new ArrayList<SSE>();
        Residue r;
        SSE s;
        Integer ligSSECount = 1;
        
        // Use the min and max atoms setting for ligands
        Integer ligMinAtoms = Settings.getInteger("plcc_I_lig_min_atoms");
        Integer ligMaxAtoms = Settings.getInteger("plcc_I_lig_max_atoms");
        Boolean noMax = false;
        if(ligMinAtoms <= 0) {
            ligMinAtoms = 1;    // 1 means no filtering because every ligand has at least one atom
        }
        if(ligMaxAtoms <= 0) {
            noMax = true;
        }
        if( (ligMinAtoms > ligMaxAtoms) && !noMax) {
            System.out.println("WARNING: Setting for minimum number of ligand atoms is > maximum setting, won't ever match.");
        }
        
        Integer numAtoms;
        Integer numLigIgnoredAtomChecks = 0;
        Integer numLigIgnoredAtomChecksTooFewAtoms = 0;
        Integer numLigIgnoredAtomChecksTooManyAtoms = 0;

        for(Integer i = 0; i < resList.size(); i++) {

            r = resList.get(i);

            if(r.isLigand()) {
                
                numAtoms = r.getNumAtoms();
                if( (numAtoms < ligMinAtoms) || ((numAtoms > ligMaxAtoms) && !noMax) ) {
                    if( (numAtoms < ligMinAtoms)) {
                        numLigIgnoredAtomChecksTooFewAtoms++;
                    }
                    if((numAtoms > ligMaxAtoms) && !noMax) {
                        numLigIgnoredAtomChecksTooManyAtoms++;
                    }
                    // Ligand did NOT pass the atom check, ignore it
                    numLigIgnoredAtomChecks++;
                    continue;
                }

                ligSSECount++;
                s = new SSE(SSE.SSECLASS_STRING_LIGAND);

                // set SSE properties
                s.addResidue(r);
                s.setSeqSseNumDssp(dsspSSElist.size() + ligSSECount);
                s.setSseType(SSE.SSECLASS_STRING_LIGAND);
                s.setChain(r.getChain());

                // set Residue properties
                r.setSSE(s);
                r.setDsspSseState(true);


                // Add the new SSE to the list of SSEs!
                //System.out.println("      Adding ligand SSE '" + s + "'.");
                ligSSElist.add(s);

            }
        }

        //System.out.println("    Found " + ligSSElist.size() + " ligand SSEs.");
        if(numLigIgnoredAtomChecks > 0) {
            if(! Settings.getBoolean("plcc_B_silent")) {
                System.out.println("    Ignored " + numLigIgnoredAtomChecks + " ligands total due to atom number constraints: ligMinAtoms=" + ligMinAtoms + ", ligMaxAtoms" + ligMaxAtoms + ". " + numLigIgnoredAtomChecksTooFewAtoms + " had too few atoms, " + numLigIgnoredAtomChecksTooManyAtoms + " had too many atoms.");
            }
        }
        
        return(ligSSElist);
    }
    
    


    /**
     * This function creates a modified list of SSEs for the PTGL. It does this by filtering all SSEs
     * which are too short and all SSEs which are not of interest for the PTGL (those which are neither helices nor beta-strands).
     * @param inputSSEs the list of input SSEs, i.e., all SSEs of this chainName which are 
     */
    private static ArrayList<SSE> createAllPtglSSEsFromDsspSSEList(List<SSE> inputSSEs) {

        ArrayList<SSE> outputSSEs = new ArrayList<SSE>();

        // If the list is emtpy we're already done. :)
        if(inputSSEs.size() < 1) {
            DP.getInstance().w("Main", "createAllPtglSSEsFromDsspSSEList(): Creating empty list of PTGL SSEs, list of DSSP SSEs is empty.");
            return(outputSSEs);
        }

        // Let's do some pre-filtering to get all the SSEs we don't want out of there first
        List<SSE> impSSEs = getImportantSSETypes(inputSSEs);
        List<SSE> consideredSSEs = removeShortSSEs(impSSEs, Settings.getInteger("plcc_I_min_SSE_length"));
        
        if(Settings.getInteger("plcc_I_debug_level") > 0) {
            printSSEList(consideredSSEs, "Considered");
        }
        
        
        //System.out.println("    Pre-filtered list of SSEs, " + consideredSSEs.size() + " left out of " + inputSSEs.size() + " (ignored SSEs of type 'B', 'T' and 'S').");

        SSE curSSE, nextSSE, afterNextSSE;
        curSSE = nextSSE = afterNextSSE = null;
        Integer numSSEsLeft;
        String cst;         // SSE type of the current SSE

        for(Integer i = 0; i < consideredSSEs.size(); i++) {

            // Fill the 3 SSE variables if enough SSEs are left
            curSSE = consideredSSEs.get(i);
            cst = curSSE.getSseType();
            numSSEsLeft = consideredSSEs.size() - (i + 1);
            //System.out.println("    At considered SSE #" + i + " of type " + curSSE.getSseType() + ", " + numSSEsLeft + " SSEs left after this one.");

            // Decide what to do based on the type of the current SSE:
            if(cst.equals("L")) {
                // We don't need to do anything special if this SSE is a ligand, just add it.
                outputSSEs.add(curSSE);
                continue;
            }
            else if(cst.equals("E")) {
                // Same goes for non-isolated beta strands
                outputSSEs.add(curSSE);
                continue;
            }
            else if(cst.equals("H") || cst.equals("I") || cst.equals("G")) {
                // We hit one of the 3 helix types. These may require merging:
                //  If we hit a helix of type I or G preceded or followed by a helix of type H, we merge them into
                //  a single helix of type H.
                //  Also we don't distinguish between the 3 types of helices, just claim they all were of type H.
                   
                // Decide what to do depending on the number of SSEs left (because we may need to merge with those).
                if(numSSEsLeft < 1) {
                    // This is the last SSE so no more merging will follow. Just add this one.
                    curSSE.setSseType("H");
                    outputSSEs.add(curSSE);
                    continue;
                }
                else {
                    // There is at least one more SSE.
                    nextSSE = consideredSSEs.get(i + 1);
                    Boolean currentSSEmergedIntoNext = false;

                    

                    // If it is a helix of type 'H', merge these two SSEs.
                    if(nextSSE.getSseType().equals("H") || nextSSE.getSseType().equals("I") || nextSSE.getSseType().equals("G")) {
                        
                        if(Settings.getBoolean("plcc_B_merge_helices")) {
                            
                            // we only merge SSEs if they are roughly adjacent, i.e., if their distance in the AA sequence is smaller than a certain threshold
                            if(curSSE.getPrimarySeqDistanceInAminoAcidsTo(nextSSE) <= Settings.getInteger("plcc_I_merge_helices_max_dist")) {
                                // To merge, we add all residues of the next SSE to this one and skip the next one.
                                // System.out.println("    Merging SSEs #" + i + " of type " + cst +  " and #" + (i + 1) + " of type " + nextSSE.getSseType()  + ".");
                                //curSSE.addResidues(nextSSE.getResidues());
                                nextSSE.addResidues(curSSE.getResidues());
                                currentSSEmergedIntoNext = true;
                                //i++;    // ignore the next SSE, we assigned its residues to this one already                                                            
                            }                                                                                    
                        }                        
                    }                    

                    // Add current SSE only if not merged with next one
                    if( ! currentSSEmergedIntoNext) {
                        curSSE.setSseType("H");
                        outputSSEs.add(curSSE);          
                    }
                }



            }
            else if(cst.equals("C")){       // coils are only considered if explicitely requested
                if(Settings.getBoolean("plcc_B_include_coils")) {
                    outputSSEs.add(curSSE);
                }
            }
            else {
                System.out.println("      WARNING: Unhandled SSE of type '" + cst + "' encountered. Ignored, this should not happen.");
            }


            if(consideredSSEs.size() > (i + 1)) {
                nextSSE = consideredSSEs.get(i + 1);
            }


        }

        //System.out.println("    PTGL adaptations of SSE list (merging etc) done, " + outputSSEs.size() +  " left out of the " + consideredSSEs.size() + " pre-filtered SSEs.");

        return(outputSSEs);

    }

    
    /**
     * Filters all SSEs that are not considered at all by the PTGL out and returns another
     * ArrayList that doesn't include those SSEs. It does NOT modify the original list.
     *
     * By default, it filters SSEs of type "B" (residue in isolated beta-bridge), "S" (bend) and "T" (hydrogen bonded turn).
     *
     */
    private static List<SSE> getImportantSSETypes(List<SSE> list) {

        List<SSE> filteredList = new ArrayList<SSE>();

        SSE s = null;
        String sT = null;
        for(Integer i = 0; i < list.size(); i++) {

            s = list.get(i);
            sT = s.getSseType();

            // Add the SSE to the list if it is some kind of helix (H,G,I), a non-isolated beta strand (E) or a ligand (L)
            if(sT.equals("H") || sT.equals("G") || sT.equals("I") || sT.equals("E") || sT.equals("L")) {
                filteredList.add(s);

                //System.out.println("    Considering SSE '" + s + "'.");
            }

            // We may want to consider coils.
            if(Settings.getBoolean("plcc_B_include_coils")) {
                if(sT.equals("C")) {
                    filteredList.add(s);
                }
            }
        }

        return(filteredList);
    }
    
    /**
     * Filters all SSEs which are shorter than minLength (i.e., consists of less than minLength residues) from the list. Returns a new list and does NOT modify the old one.
     * @param list the input list, it is not modified
     * @param minLength the minimal length of the SSEs to keep, in residues
     * @return a new list which is a filtered version of the old one
     */
    private static List<SSE> removeShortSSEs(List<SSE> list, Integer minLength) {

        List<SSE> filteredList = new ArrayList<SSE>();

        for(SSE s : list) {
            if(s.getLength() >= minLength) {
                filteredList.add(s);
            }            
        }

        return(filteredList);
    }

    
    /**
     * Prints a formatted list of SSEs to STDOUT.
     * @param sl the list of SSEs
     * @param title the title to use for the output 
     */
    public static void printSSEList(List<SSE> sl, String title) {

        // just a counting aid
        System.out.println("---[ " + title + " ] (" + sl.size() + " SSEs) ---");
        System.out.println("           1         2         3         4         5         6         7         8         9");
        System.out.println("  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890");

        if(sl.size() < 1) {
            System.out.println("  <no_SSEs>");
        }
        else {

            // SSE list, only the type
            System.out.print("  ");
            for(Integer j = 0; j < sl.size(); j++) {
                System.out.print(sl.get(j).getSseType());
            }
            System.out.print("\n");

            // details for each SSE
            System.out.print(" ");
            for(Integer j = 0; j < sl.size(); j++) {
                System.out.print(" " + sl.get(j).shortStringRep());
            }
            System.out.print("\n");
        }            
    }

    
    /**
     * Merges 2 lists of SSEs and returns a new list containing all SSEs from both lists.
     * @param listA the first SSE list
     * @param listB the second SSE list
     * @return the merged list
     */
    public static List<SSE> mergeSSEs(List<SSE> listA, List<SSE> listB) {
        List<SSE> listMerged = new ArrayList<>();

        for(Integer i = 0; i < listA.size(); i++) {
            listMerged.add(listA.get(i));
        }

        for(Integer j = 0; j < listB.size(); j++) {
            listMerged.add(listB.get(j));
        }

        return(listMerged);
    }
    
    
    /**
     * Draws a Ramachandran plot (phi and psi angels of the residues) to an image file. Note that the angles are parsed from the DSSP
     * file but consider our SSE modifications. Still, the plot shows more about DSSP than about plcc.
     * 
     * @param filePath the output path of the image without extension
     * @param residues the Residues, in an ArrayList
     * @param label the label
     * @return true if the file was written, false otherwise
     */
    public static Boolean drawRamachandranPlot(String filePath, ArrayList<Residue> residues, String label) {
        
        
        // Apache Batik SVG library, using W3C DOM tree implementation
        // Get a DOMImplementation.
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);
        // Create an instance of the SVG Generator.
        //SVGGraphics2D ig2 = new SVGGraphics2D(document);
        
        
        

        // All these values are in pixels
        // page setup
        Integer marginLeft = 40;
        Integer marginRight = 40;
        Integer marginTop = 40;
        Integer marginBottom = 40;

        // The header that contains the text describing of the graph.
        Integer headerHeight = 40;

        // The footer contains the vertex numbering.
        Integer footerHeight = 40;
        
        // Where to start drawing
        Integer headerStartX = marginLeft;
        Integer headerStartY = marginTop;
        Integer imgStartX = marginLeft;
        Integer imgStartY = marginTop + headerHeight;
        Integer footerStartX = marginLeft;
        
        
        Integer plotWidth = 800;
        Integer plotHeight = 800;
    
        // the image area: the part where the vertices and arcs are drawn
        Integer imgWidth = plotWidth;
        Integer imgHeight = plotHeight;              

        // where to start drawing the vertices
        Integer footerStartY = imgStartY + plotHeight + 40;
        
        // putting it all together
        Integer pageWidth = marginLeft + imgWidth + marginRight;
        Integer pageHeight = marginTop + headerHeight + imgHeight + footerHeight + marginBottom;


        
        try {

            // ------------------------- Prepare stuff -------------------------
            // TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed into integer pixels
            BufferedImage bi = new BufferedImage(pageWidth, pageHeight, BufferedImage.TYPE_INT_ARGB);
            Graphics2D ig2 = bi.createGraphics();
            ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);

            // make background 
            ig2.setPaint(Color.LIGHT_GRAY);
            ig2.fillRect(0, 0, pageWidth, pageHeight);            
            ig2.setPaint(Color.BLACK);


            // prepare font
            Font font = new Font("TimesRoman", Font.PLAIN, 20);
            ig2.setFont(font);
            FontMetrics fontMetrics = ig2.getFontMetrics();
            
            java.awt.Shape shape;

            // ------------------------- Draw header -------------------------

            // check width of header string
            String plotHeader = label;
            Integer stringWidth = fontMetrics.stringWidth(plotHeader);       // Should be around 300px for the text above
            Integer stringHeight = fontMetrics.getAscent();
            String sseNumber;
            
            ig2.drawString(plotHeader, headerStartX, headerStartY);
            // draw axis labels
            
            font = new Font("TimesRoman", Font.PLAIN, 16);
            ig2.setFont(font);
            
            Integer xAxisLabelStartY = headerStartY + 30;
            
            // horizontal labels
            ig2.drawString("-180°", headerStartX, xAxisLabelStartY);
            ig2.drawString("-phi", headerStartX + (plotWidth / 4), xAxisLabelStartY);
            ig2.drawString("0°", headerStartX + (plotWidth / 2), xAxisLabelStartY);
            ig2.drawString("+phi", headerStartX + (plotWidth * 3/4), xAxisLabelStartY);
            ig2.drawString("+180°", headerStartX + plotWidth, xAxisLabelStartY);
            
            // vertical labels
            xAxisLabelStartY += 20;
            Integer yAxisLabelStartX = 5;
            ig2.drawString("+180°", yAxisLabelStartX, xAxisLabelStartY);
            ig2.drawString("+psi", yAxisLabelStartX, xAxisLabelStartY  + (plotHeight / 4));
            ig2.drawString("0°", yAxisLabelStartX, xAxisLabelStartY + (plotHeight / 2));
            ig2.drawString("-psi", yAxisLabelStartX, headerStartX + (plotHeight * 3/4));
            ig2.drawString("-180°", yAxisLabelStartX, headerStartX + plotHeight);
            
            // footer
            ig2.drawString("Color codes indicate SSE type of residues: red=helix, black=beta strand, gray=coil/other.", footerStartX, footerStartY);
            
            // ------------------------- Draw area and the borders around the drawing area ----------------------------
            ig2.setPaint(Color.WHITE);
            ig2.fillRect(imgStartX, imgStartY, plotWidth, plotHeight);
            
            // draw it                        
            ig2.setStroke(new BasicStroke(1));
            Rectangle2D border = new Rectangle2D.Double(imgStartX, imgStartY, plotWidth, plotHeight);
            shape = ig2.getStroke().createStrokedShape(border);
            ig2.draw(shape);
            
            // draw x and y axis
            ig2.setPaint(Color.LIGHT_GRAY);
            ig2.drawLine(imgStartX, imgStartY + (plotHeight / 2), imgStartX + plotWidth, imgStartY +  + (plotWidth / 2)); // y
            ig2.drawLine(imgStartX + (plotWidth / 2), imgStartY, imgStartX + (plotWidth / 2), imgStartY + plotHeight); // x
            

            // ------------------------- Draw the plot -------------------------
            
            // Draw the edges as arcs
            
            Rectangle2D.Double dot;
            Residue res;
            Integer edgeType, arcCenterX, arcCenterY, leftVert, rightVert, leftVertPosX, rightVertPosX, arcWidth, arcHeight, arcTopLeftX, arcTopLeftY, spacerX, spacerY;
            for(Integer i = 0; i < residues.size(); i++) {
                res = residues.get(i);

                // Choose color
                if(res.getSSEString().equals("H")) { ig2.setPaint(Color.RED); }
                else if(res.getSSEString().equals("E")) { ig2.setPaint(Color.BLACK); }
                else if(res.getSSEString().equals("L")) { ig2.setPaint(Color.MAGENTA); }
                else if(res.getSSEString().equals("C")) { ig2.setPaint(Color.GRAY); }                
                else { ig2.setPaint(Color.LIGHT_GRAY); }
                
                // skip ligands and other non-AAs
                if(! res.isAA()) {
                    continue;
                }

                // determine the position of the dot within the drawing area
                Double phiNorm = (res.getPhi() + 180.0) / 360.0;
                Double psiNorm = (res.getPsi() + 180.0) / 360.0;
                Double posDrawX = phiNorm * imgWidth;
                Double posDrawY = psiNorm * imgHeight;
                
                // now determine the absolute position in the canvas
                Double posCanvasX = imgStartX + posDrawX;
                Double posCanvasY = imgStartY + posDrawY;
                
                Double dotWidth = 1.0;
                Double dotHeight = 1.0;
                

                // draw it                        
                ig2.setStroke(new BasicStroke(1));
                dot = new Rectangle2D.Double(posCanvasX, posCanvasY, dotWidth, dotHeight);
                shape = ig2.getStroke().createStrokedShape(dot);
                ig2.draw(shape);

            }
            

            
            // all done, write the image to disk
            ImageIO.write(bi, "PNG", new File(filePath + ".png"));
            //ig2.stream(new FileWriter(filePath + ".svg"), false);
            ig2.dispose();

        } catch (Exception e) {
            DP.getInstance().w("Could not write image file for ramachandran plot to file '" + filePath + "'. Check permissions.");
            return(false);
        }

        return(true);        
    }
    
    
    /**
     * Ends plcc execution via calling System.exit(), but does maintenance work like closing open DB connections before that.
     * @param exitCode the exit code to give to the System.exit() function
     */
    public static void doExit(int exitCode) {
        if(Settings.getBoolean("plcc_B_useDB")) {
            DBManager.closeConnection();
        }
        System.exit(exitCode);
    }
    
    /**
     * Informs the user that the currently selected DEBUG settings brake the output of this program, e.g. because
     * important computations are aborted after the first few residues.
     */
    public static void print_debug_malfunction_warning() {
        System.out.println("###################################################################################");
        System.out.println("# WARNING: Contact debug mode ACTIVE. Will abort after the first few residues.    #");
        System.out.println("# WARNING: This mode will produce no or WRONG SSE level results! Do NOT use them. #");
        System.out.println("###################################################################################");
    }
    
    /**
     * Calculates complex graph types which are configured in the config file for all given chains.
     * @param allChains a list of chains, each chainName will be handled separately
     * @param resList a list of residues
     * @param resContacts a list of residue contacts
     * @param pdbid the PDBID of the protein, required to name files properly etc.
     * @param outputDir where to write the output files. the filenames are deduced from graph type and pdbid.
     * @param graphType the graph type, one of the constants like SSEGraph.GRAPHTYPE_ALBE 
     */
    public static void calculateComplexGraph(List<Chain> allChains, List<Residue> resList, List<MolContactInfo> resContacts, String pdbid, String outputDir, String graphType) {
        
        Boolean silent = Settings.getBoolean("plcc_B_silent");
        
        ComplexGraphResult cgr = new ComplexGraphResult();
        
        if(! silent) {
            System.out.println("Calculating complex graph (CG) of type " + graphType + ".");
        }
        
        Chain c;
        List<SSE> chainDsspSSEs = new ArrayList<SSE>();
        List<SSE> chainPtglSSEs = new ArrayList<SSE>();
        List<SSE> allChainSSEs = new ArrayList<SSE>();
        List<SSE> oneChainSSEs = new ArrayList<SSE>();
        List<Integer> chainEnd = new ArrayList<Integer>();
        HashMap<String, List<SSE>> chainSSEMap = new HashMap<>();
        String fs = System.getProperty("file.separator");
        String fileNameSSELevelWithExtension = null;
        String fileNameSSELevelWithoutExtension = null;
        String fileNameChainLevelWithExtension = null;
        String fileNameChainLevelWithoutExtension = null;        
        String filePathImg = null;
        String filePathGraphs = null;
        String filePathHTML = null;
        String imgFile = null;
        String plccGraphFileSSELevel = null;

        if(! silent) {
            System.out.println("  Calculating CG SSEs for all chains of protein " + pdbid + "...");
        }

        HashMap<String, String> md = FileParser.getMetaData();
        
        Double res;
        try {
            res = Double.valueOf(md.get("resolution"));
        } catch (Exception e) {
            res = -1.0;
            System.err.println("WARNING: Could not determine resolution of PDB file for protein '" + pdbid + "', assuming NMR with resolution '" + res + "'.");
            
        }
        
        List<String> keepSSEs = new ArrayList<String>();
        List<SSE> filteredChainSSEs;

        // Check whether coils should be kept
        if(Settings.getBoolean("plcc_B_include_coils")) {
            keepSSEs.add(Settings.get("plcc_S_coilSSECode"));
        }

        // Filter SSEs depending on the requested graph type
        switch (graphType) {
            case SSEGraph.GRAPHTYPE_ALBE:
                keepSSEs.add(SSE.SSECLASS_STRING_STRAND);
                keepSSEs.add(SSE.SSECLASS_STRING_HELIX);
                break;
            case SSEGraph.GRAPHTYPE_ALPHA:
                keepSSEs.add(SSE.SSECLASS_STRING_HELIX);
                break;
            case SSEGraph.GRAPHTYPE_BETA:
                keepSSEs.add(SSE.SSECLASS_STRING_STRAND);
                break;
            case SSEGraph.GRAPHTYPE_ALBELIG:
                keepSSEs.add(SSE.SSECLASS_STRING_STRAND);
                keepSSEs.add(SSE.SSECLASS_STRING_HELIX);
                keepSSEs.add(SSE.SSECLASS_STRING_LIGAND);
                break;
            case SSEGraph.GRAPHTYPE_ALPHALIG:
                keepSSEs.add(SSE.SSECLASS_STRING_HELIX);
                keepSSEs.add(SSE.SSECLASS_STRING_LIGAND);
                break;
            case SSEGraph.GRAPHTYPE_BETALIG:
                keepSSEs.add(SSE.SSECLASS_STRING_STRAND);
                keepSSEs.add(SSE.SSECLASS_STRING_LIGAND);
                break;
            default:
                System.err.println("ERROR: calculateComplexGraph(): Graph type '" + graphType + "' invalid. Skipping.");
                return;
        }
        
        // print warning here so that it only appears once (not once per chain)
        if (Settings.getBoolean("plcc_B_use_mmCIF_parser")) {
            DP.getInstance().w("CIF setting enabled: Parsing of protein meta info not fully " + 
                        "implemented yet.");
        }

        // Get SSEs for all chains
        for(Integer i = 0; i < allChains.size(); i++) {
            c = allChains.get(i);
            
            //if(! silent) {
            //    System.out.println("   *Handling chainName " + c.getPdbChainID() + ".");
            //}

            // CIF parser currently does not parse ProtMetaInfo in all means
            //   -> just ignore it here, so it does not land in md and causes no trouble
            if (! Settings.getBoolean("plcc_B_use_mmCIF_parser")) {
                ProtMetaInfo pmi = FileParser.getMetaInfo(pdbid, c.getPdbChainID());
                //pmi.print();
                md.put("pdb_mol_name", pmi.getMolName());
                md.put("pdb_org_sci", pmi.getOrgScientific());
                md.put("pdb_org_common", pmi.getOrgCommon());
            }

            // determine SSEs for this chainName
            //System.out.println("    Creating all SSEs for chainName '" + c.getPdbChainID() + "' consisting of " + c.getResidues().size() + " residues.");
            chainDsspSSEs = createAllDsspSSEsFromResidueList(c.getResidues());
            
            oneChainSSEs = createAllPtglSSEsFromDsspSSEList(chainDsspSSEs);
            
            List<SSE> chainLigSSEs =  createAllLigandSSEsFromResidueList(c.getResidues(), chainDsspSSEs);
            oneChainSSEs = mergeSSEs(oneChainSSEs, chainLigSSEs);

            /*
            if(! silent) {
                StringBuilder sb = new StringBuilder();
                sb.append("    SSEs: ");
                for(Integer j = 0; j < oneChainSSEs.size(); j++) {
                    sb.append(oneChainSSEs.get(j).getSseType());
                }
            }
            */

            // SSEs have been calculated, now assign the PTGL labels and sequential numbers on the chainName
            for(Integer j = 0; j < oneChainSSEs.size(); j++) {
                oneChainSSEs.get(j).setSeqSseChainNum(j + 1);   // This is the correct value, determined from the list of all valid SSEs of this chainName
                oneChainSSEs.get(j).setSseIDPtgl(getPtglSseIDForNum(j));
            }


            // Filter SSEs.
            filteredChainSSEs = filterAllSSEsButList(oneChainSSEs, keepSSEs);
            if(chainEnd.size()>0) {
                chainEnd.add(chainEnd.get(chainEnd.size() - 1) + filteredChainSSEs.size());
            }
            else {
                chainEnd.add(filteredChainSSEs.size());
            }
            allChainSSEs.addAll(filteredChainSSEs);
            chainSSEMap.put(c.getPdbChainID(), oneChainSSEs);
        }
                
        ComplexGraph compGraph = new ComplexGraph(pdbid, allChains, resContacts, Settings.getBoolean("plcc_B_writeComplexContactCSV"));    
        
        if(Settings.getBoolean("plcc_B_useDB")) {
            if(! silent) {System.out.print("    Writing chain complex contact info to DB...");}
        
            if(compGraph.writeChainComplexContactInfoToDB()){
                if(! silent) { System.out.println(" successfull!");}
            } else {
                DP.getInstance().e("Main", "Writing chain complex contact info to database FAILED!"); 
            }
            
            if(! silent) {System.out.println("    Writing SSE complex contact info to DB...");}
        
            compGraph.writeSSEComplexContactInfoToDB(pdbid);
        }     

        cgr.setCompGraph(compGraph);     
        
        if( ! silent) {
            System.out.println("  Preparing to write complex graph files (CG verts / edges : " + compGraph.getVertices().size() + " / " + compGraph.getEdges().size() +  ").");
        }
        
        //write Residue contact info to csv. 
        if(Settings.getBoolean("plcc_B_writeComplexContactCSV")) {
            try {
                FileWriter writer = new FileWriter(pdbid+"_contact_info.csv");  // TODO: This does not respect a possible subdir tree ('/ti/7tim/...') setting yet
                for(String x : compGraph.getContactInfo()){
                    writer.append(x);
                    writer.append("\n");
                }
                writer.flush();
                writer.close();   
                if( ! silent) {
                    System.out.println("    Contact info CSV written.");
                }
            } catch(IOException e){
                DP.getInstance().e("Main", "Failed to write complex graph contact info CSV file: '" + e.getMessage() + "'.");
            }
        }
        
        
        
        // Calculate SSE level contacts
        ContactMatrix chainCM = new ContactMatrix(allChainSSEs, pdbid);
        chainCM.restrictToChain("ALL");
        chainCM.fillFromContactList(resContacts, keepSSEs);
        chainCM.calculateSSEContactMatrix();                       
        
        if( ! silent) {
            System.out.println("    SSE contact matrix calculated.");
        }
        
        chainCM.calculateSSESpatialRelationMatrix(resContacts, false);                

        if( ! silent) {
            System.out.println("    Spatial relation matrix calculated.");
        }
        
        // This graph is still required because it is used for drawing the VPLG-style picture
        ProtGraph SseCg = chainCM.toProtGraph();                        
        SseCg.declareProteinGraph();      // this is declared a CG after setting info, see below.  
        
        if( ! silent) {
            System.out.println("    Graph created and declared a CG.");
        }

        
        if(Settings.getBoolean("plcc_B_forceBackboneContacts")) {
            if( ! silent) {
                System.out.println("      Adding backbone contacts to consecutive SSEs of the " + graphType + " graph.");
            }
            SseCg.addFullBackboneContacts();            
        }

        //System.out.println("    Done creating SSE-level " + graphType + " CG.");
        SseCg.setInfo(pdbid, "ALL", "ALL", "complex_" + graphType);
        SseCg.addMetadata(md);
        SseCg.setComplexData(chainEnd, allChains);
        SseCg.declareComplexGraph(true);
            
        /*
        if(Settings.getBoolean("plcc_B_compute_graph_metrics")) {
            // ###TEST-CG-METRICS
            System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%% Computing CCs of CG %%%%%%%%%%%%%%%%%%%%%%%%%");
            ArrayList<FoldingGraph> conCompsOfCG = SseCg.getConnectedComponents();
            FoldingGraph cgSubgraph = SseCg.getLargestConnectedComponent();
            Map<Integer, Integer> degreeDistrGraph = GraphMetrics.degreeDistribution(SseCg, false);
            Map<Integer, Integer> degreeDistrLargestCC = GraphMetrics.degreeDistribution(cgSubgraph, false);
            System.out.println("Largest CC has size " + cgSubgraph.getSize());
        }
                */
        
        filePathImg = outputDir;
        filePathGraphs = outputDir;
        filePathHTML = outputDir;
        String coils = "";
        if(Settings.getBoolean("plcc_B_include_coils")) {
            //System.out.println("  Considering coils, this may fragment SSEs.");
            coils = "_coils";
        }        
        
        fileNameSSELevelWithoutExtension = pdbid + "_complex_sses_" + graphType + coils + "_CG";
        fileNameChainLevelWithoutExtension = pdbid + "_complex_chains"  + coils + "_CG";    // the chainName-level graph always uses the full contact data, no need to differentiate by graphType (could be implemented of course)
        fileNameSSELevelWithExtension = fileNameSSELevelWithoutExtension + Settings.get("plcc_S_img_output_fileext");
        fileNameChainLevelWithExtension = fileNameChainLevelWithoutExtension + Settings.get("plcc_S_img_output_fileext");

        // Create the file in a subdir tree based on the protein meta data if requested
        File targetDir = null;
        if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {

            targetDir = IO.generatePDBstyleSubdirTreeName(new File(outputDir), pdbid, "ALL");
            if(targetDir != null) {
                ArrayList<String> errors = IO.createDirIfItDoesntExist(targetDir);
                if( ! errors.isEmpty()) {
                    for(String err : errors) {
                        System.err.println("ERROR: " + err);
                    }
                } else {
                    filePathImg = targetDir.getAbsolutePath();
                    filePathGraphs = targetDir.getAbsolutePath();
                    filePathHTML = targetDir.getAbsolutePath();
                }                    
            } else {
                System.err.println("ERROR: Could not determine PDB-style subdir path name.");
            }
        }
        
        //System.out.println("CG TARGET DIR IS " + targetDir + "." );
        String dbImagePathCGNoExt = fileNameSSELevelWithoutExtension;
        if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
            dbImagePathCGNoExt = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, "ALL") + fs + fileNameSSELevelWithoutExtension;
        }
        
        // ------------- chainName level complex graphs -----------
        // the simple complex graph (one vertex is one chainName):
        File gmlFileChainLevel = null;
        String gmlFileNameChainLevel = null;
        try {
            gmlFileChainLevel = new File(filePathGraphs + fs + pdbid + "_complex_chains_" + graphType + "_CG.gml");
            gmlFileChainLevel.createNewFile();
            if(compGraph.writeToFileGML(gmlFileChainLevel)) {
                if(! silent) {
                    System.out.println("    Complex graph chain-level notation written to file '" + gmlFileChainLevel.getName() + "' in GML format.");
                }
                cgr.setComGraphFileGML(gmlFileChainLevel);
                gmlFileNameChainLevel = fileNameChainLevelWithoutExtension + ".gml";
            } else {
                System.err.println("ERROR: Could not write complex graph to file '" + gmlFileChainLevel.getAbsolutePath() + "'.");
            }
        } catch(IOException e){
            System.err.println("ERROR: Could not write complex graph to file '" + gmlFileChainLevel.getAbsolutePath() + ": '" + e.getMessage()+ "'.");
        }
        
        // ------------------ write the SSE level graphs -----------
        if( ! silent) {
            System.out.println("    Writing SSE level graph text files...");
        }
        
        // the detailed complex graph (each vertex is one SSE, vertices ordered by chainName):
        String graphFormatsWrittenSSELevel = "";        
        Integer numFormatsWrittenSSELevel = 0;
        if(Settings.getBoolean("plcc_B_output_compgraph_GML")) {
            String gmlfFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".gml";                       
            if(IO.stringToTextFile(gmlfFileSSELevel, SseCg.toGraphModellingLanguageFormat())) {
                graphFormatsWrittenSSELevel += "gml "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_TGF")) {
            String tgfFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".tgf";
            if(IO.stringToTextFile(tgfFileSSELevel, SseCg.toTrivialGraphFormat())) {
                graphFormatsWrittenSSELevel += "tgf "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_DOT")) {
            String dotLangFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".gv";
            if(IO.stringToTextFile(dotLangFileSSELevel, SseCg.toDOTLanguageFormat())) {
                graphFormatsWrittenSSELevel += "gv "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_kavosh")) {
            String kavoshFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".kavosh";
            if(IO.stringToTextFile(kavoshFileSSELevel, SseCg.toKavoshFormat())) {
                graphFormatsWrittenSSELevel += "kavosh "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_XML")) {
            String xmlFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".xml";
            if(IO.stringToTextFile(xmlFileSSELevel, SseCg.toXMLFormat())) {
                graphFormatsWrittenSSELevel += "xml "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_JSON")) {
            String jsonFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".json";
            if(IO.stringToTextFile(jsonFileSSELevel, SseCg.toJSONFormat())) {
                graphFormatsWrittenSSELevel += "json "; numFormatsWrittenSSELevel++;
            }
        }
        if(Settings.getBoolean("plcc_B_output_compgraph_eld")) {
            String elFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".el_edges";
            String nodeTypeListFile = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".el_ntl";
            if(IO.stringToTextFile(elFileSSELevel, SseCg.toEdgeList()) && IO.stringToTextFile(nodeTypeListFile, SseCg.getNodeTypeList())) {
                graphFormatsWrittenSSELevel += "el "; numFormatsWrittenSSELevel++;
            }
        }
        // write the SSE info text file for the image (plcc graph format file)
        if(Settings.getBoolean("plcc_B_output_compgraph_plcc")) {
            plccGraphFileSSELevel = filePathGraphs + fs + fileNameSSELevelWithoutExtension + ".plg";
            if(IO.stringToTextFile(plccGraphFileSSELevel, SseCg.toVPLGGraphFormat())) {
                graphFormatsWrittenSSELevel += "plg "; numFormatsWrittenSSELevel++;
            }
        }
        
        if(numFormatsWrittenSSELevel > 0) {
            if(! silent) {
                System.out.println("    Exported complex graph in " + numFormatsWrittenSSELevel + " formats (" + graphFormatsWrittenSSELevel + ") to '" + new File(filePathGraphs).getAbsolutePath() + fs + "'.");
            }
        }


        String imgFileNoExt = filePathImg + fs + fileNameSSELevelWithoutExtension;
        String imgFileChainComplexNoExt = filePathImg + fs + fileNameChainLevelWithoutExtension;
        //imgFile = filePathImg + fs + fileNameWithExtension;
                
        if(Settings.getBoolean("plcc_B_draw_graphs")) {
            IMAGEFORMAT[] formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG, DrawTools.IMAGEFORMAT.PDF };
            ProteinGraphDrawer.drawProteinGraph(imgFileNoExt, false, formats, SseCg, new HashMap<Integer, String>(), new ArrayList<String>());
            if(! silent) {
                System.out.println("    Image of complex graph written to base file '" + imgFileNoExt + "'.");
            }
            
            Map<String, String> molInfoForChains = new HashMap<>();
            for(Chain tc : SseCg.getAllChains()) {               
                molInfoForChains.put(tc.getPdbChainID(), tc.getMacromolID());
            }
            
            HashMap<DrawTools.IMAGEFORMAT, String> drawnFormats = ComplexGraph.drawComplexGraph(imgFileChainComplexNoExt, false, formats, compGraph, molInfoForChains);
            if(! silent) {
                for(IMAGEFORMAT f : drawnFormats.keySet()) {
                    System.out.println("    Complex graph drawn in format " + f + " to file '" + drawnFormats.get(f) + "'.") ;
                }
            }
        }
        else {
            if(! silent) {
                System.out.println("    Image output disabled, not drawing complex graphs.");
            }
        } 
        
        // database
        if(Settings.getBoolean("plcc_B_useDB")) {
            String dbImagePathCG = fileNameSSELevelWithoutExtension;
            String dbImagePathChainCG = fileNameChainLevelWithoutExtension;
            if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                dbImagePathCG = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, "ALL") + fs + fileNameSSELevelWithoutExtension;
                dbImagePathChainCG = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, "ALL") + fs + fileNameChainLevelWithoutExtension;
            }
            //System.out.println("dbImagePathCG = '" + dbImagePathCG + "'");
            //System.out.println("dbImagePathChainCG = '" + dbImagePathChainCG + "'");
            //dbImagePath += DrawTools.getFileExtensionForImageFormat(format);
            try {
                if(Settings.getBoolean("plcc_B_write_graphstrings_to_database_cg")) {
                    DBManager.writeComplexGraphToDB(pdbid, SseCg.toGraphModellingLanguageFormat(), null, SseCg.toXMLFormat(), null, SseCg.toKavoshFormat(), null, dbImagePathCG + ".svg", dbImagePathChainCG + ".svg", dbImagePathCG + ".png", dbImagePathChainCG + ".png", dbImagePathCG + ".pdf", dbImagePathChainCG + ".pdf");
                } else {
                    DBManager.writeComplexGraphToDB(pdbid, null, null, null, null, null, null, dbImagePathCG + ".svg", dbImagePathChainCG + ".svg", dbImagePathCG + ".png", dbImagePathChainCG + ".png", dbImagePathCG + ".pdf", dbImagePathChainCG + ".pdf");
                }
                
                if(! silent) {
                    System.out.println("Wrote complex graph of " + pdbid + " to DB.");
                }
                
            } catch(SQLException e) {
                DP.getInstance().w("Main", "Could not write complex graph to DB: '" + e.getMessage() + "'.");
            }
        } 
        
        if(Settings.getBoolean("plcc_B_compute_graph_metrics") && graphType.equals(SSEGraph.GRAPHTYPE_ALBELIG)) {
            GraphProperties gp = new GraphProperties(SseCg);
            GraphProperties sgp = new GraphProperties(gp.getLargestConnectedComponent());
            
            
            
            // DEBUG ---------------
            /*
            System.out.println("??????????????????????????????");
            System.out.println("CG size " + SseCg.getSize());
            SimpleGraphInterface testCG = (SimpleGraphInterface)SseCg;
            System.out.println("CG as SGI size:" + testCG.getSize());
            
            List<FoldingGraph> lt = SseCg.getConnectedComponents();
            System.out.print("FGs:");
            for(FoldingGraph t : lt) {
                System.out.print(t.getSize() + " ");
            }
            System.out.println("");
                  
            
            
            SimpleGraphInterface test = gp.getLargestConnectedComponent();
            System.out.println("largest CC size:" + test.getSize());
            List<SimpleGraphInterface> tests = gp.getConnectedComponents();
            System.out.print("All " + tests.size() + " CC sizes:");
            for(SimpleGraphInterface t : tests) {
                System.out.print(t.getSize() + " ");
            }
            System.out.println("");
            // DEBUG ----------------
              */
            
            if(Settings.getBoolean("plcc_B_useDB")) {

                if( ! DBManager.getAutoCommit()) {
                    DBManager.commit();
                }

                try {
                    Long graph_db_id = DBManager.getDBComplexgraphID(pdbid);
                    if(graph_db_id > 0L) {
                        //System.out.println("Found complex graph " + pdbid + " with ID " + graph_db_id + ".");
                        // write graph properties
                        Long runtime_secs = null;
                        DBManager.writeComplexgraphStatsToDB(graph_db_id, Boolean.FALSE, gp.getNumVertices(), gp.getNumEdges(), gp.getMinDegree(), gp.getMaxDegree(), gp.getConnectedComponents().size(), gp.getGraphDiameter(), gp.getGraphRadius(), gp.getAverageClusterCoefficient(), gp.getAverageShortestPathLength(), gp.getDegreeDistributionUpTo(50), gp.getAverageDegree(), gp.getDensity(), gp.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);
                        // write properties of largest CC of graph
                        DBManager.writeComplexgraphStatsToDB(graph_db_id, Boolean.TRUE, sgp.getNumVertices(), sgp.getNumEdges(), sgp.getMinDegree(), sgp.getMaxDegree(), sgp.getConnectedComponents().size(), sgp.getGraphDiameter(), sgp.getGraphRadius(), sgp.getAverageClusterCoefficient(), sgp.getAverageShortestPathLength(), sgp.getDegreeDistributionUpTo(50), sgp.getAverageDegree(), sgp.getDensity(), sgp.getCumulativeDegreeDistributionUpToAsArray(50), runtime_secs);
                    }
                    else {
                        DP.getInstance().e("Main", "Could not write complex graph properties to DB, graph not found in database.");
                    }
                } catch(SQLException e) {
                    DP.getInstance().e("SQL error while trying to store complex graph stats: '" + e.getMessage()+ "'.");
                }
            }
        }
        
        if(Settings.getBoolean("plcc_B_draw_graphs") && Settings.getBoolean("plcc_B_draw_ligandcomplexgraphs") && graphType.equals(SSEGraph.GRAPHTYPE_ALBELIG)) {
            if(! silent) {
                System.out.println("    Drawing ligand-centered complex graphs...");
            }
            
            
            
            // Determine all ligands of all chains
            List<SSE> ligandsAllChains = new ArrayList<>();            
            for(SSE s : allChainSSEs) {
                if( s.isLigandSSE()) {
                    ligandsAllChains.add(s);
                }
            }
            
            if(! silent) {
                System.out.println("     Found " + ligandsAllChains.size() + " ligand SSEs total in the " + SseCg.getAllChains().size() + " chains.");
            }
            
            String ligimgFileNoExt, ligName;
            for(SSE ligandSSE : ligandsAllChains) {
                // OK, now we handle the ligands
                Chain lc = ligandSSE.getChain();
                String ligChainName = ligandSSE.getChain().getPdbChainID();
                Integer ligRes = ligandSSE.getStartResidue().getPdbNum();
                String lign3 = ligandSSE.getTrimmedLigandName3();
                ligName = ligChainName + "-" + ligRes + "-" + lign3;  // something like "A-234-ICT", meaning isocitric acid, PDB residue 234 of chainName A
                
                // determine all chains the ligand has contacts with (based on the SSEs it has contacts with):
                Integer ligIndex = SseCg.getSSEIndex(ligandSSE);
                List<String> ligContactChains = new ArrayList<>();
                
                List<String> ignoreChains = new ArrayList<>();
                for(Chain ic : SseCg.getAllChains()) {
                    ignoreChains.add(ic.getPdbChainID()); // we will delete the ones we are interested in later, see below
                }
                
                if(ligIndex < 0) {
                    DP.getInstance().e("Main", "Could not get index of ligand SSE '" + ligName + "' from complex graph for ligand-centered graph computation, skipping CLG.");
                    continue;
                }
                List<Integer> contactSSEIndices = SseCg.neighborsOf(ligIndex);
                
                for(Integer sseIndex : contactSSEIndices) {
                    String contactChainName = SseCg.getChainNameOfSSE(sseIndex);                    
                    if( ! ligContactChains.contains(contactChainName)) { ligContactChains.add(contactChainName); }
                }
                
                for(String relChainName : ligContactChains) {
                    ignoreChains.remove(relChainName);  // remove the chains which are relevant for this ligand from the ignore list (which contained all chains before)
                }
                
                // TODO: now we need to restrict the following graphs so that they only consider the chains we determined
                if(! silent) {
                    List<String> contactSSENames = new ArrayList<>();
                    for(Integer sseIndex : contactSSEIndices) {
                        contactSSENames.add(SseCg.getChainNameOfSSE(sseIndex) + "-" + SseCg.getVertex(sseIndex).getSSESeqChainNum() + "-" +  SseCg.getVertex(sseIndex).getSSEClass());    // something like "A-1-H", meaning the first SSE of chainName A, a helix
                    }
                    if(! silent) {
                        System.out.println("     *Ligand '" + ligName + "' is in contact with the following " + contactSSENames.size() + " SSEs: '" + IO.stringListToString(contactSSENames) + "'.");
                        System.out.println("      Ligand '" + ligName + "' is in contact with the following " + ligContactChains.size() + " chains: '" + IO.stringListToString(ligContactChains) + "'.");
                        System.out.println("      Ligand '" + ligName + "' LCG ignores the following " + ignoreChains.size() + " chains: '" + IO.stringListToString(ignoreChains) + "'.");
                    }                                        
                }    
                
                if(ligContactChains.size() < 2) {
                    if(! silent) {
                        System.out.println("      Ligand '" + ligName + "' only has contacts to a single chain, skipping its ligand-centered complex graph (just use the normal albelig graph).");
                    }
                    continue;
                }
                
                String ligfileNameSSELevelWithoutExtension = pdbid + "_ligand_complex_sses_" + ligName + "_" + graphType + coils + "_LCG";
                ligimgFileNoExt = filePathImg + fs + ligfileNameSSELevelWithoutExtension;

                IMAGEFORMAT[] formats = new IMAGEFORMAT[]{ DrawTools.IMAGEFORMAT.PNG, DrawTools.IMAGEFORMAT.PDF };
                Map<Integer, String> sseDrawLabels = new HashMap<Integer, String>();
                sseDrawLabels.put(ligIndex, lign3 + "-" + ligRes);
                
                // change graph info, this is so that the label on the image gets set properly
                SseCg.setInfo(pdbid, "ALL", "ALL", "ligand_complex_" + lign3 + "-" + ligRes);
                
                HashMap<DrawTools.IMAGEFORMAT, String> outCLGImages = ProteinGraphDrawer.drawProteinGraph(ligimgFileNoExt, false, formats, SseCg, sseDrawLabels, ignoreChains); // draw ligand-based complex graph
                if(! silent) {
                    System.out.println("      Image of ligand-centered complex graph for ligand '" + ligName + "' written to base file '" + imgFileNoExt + "'.");
                }      
                
                // database
                if(Settings.getBoolean("plcc_B_useDB")) {
                    
                    // write info on the graph itself to the DB
                    String dbImagePathLCG = ligimgFileNoExt;
                    if(Settings.getBoolean("plcc_B_output_images_dir_tree") || Settings.getBoolean("plcc_B_output_textfiles_dir_tree")) {
                        dbImagePathLCG = IO.getRelativeOutputPathtoBaseOutputDir(pdbid, "ALL") + fs + ligfileNameSSELevelWithoutExtension;
                    }

                    Long chain_db_id = DBManager.getDBChainID(pdbid, ligChainName);
                    if(chain_db_id <= 0L) { 
                        DP.getInstance().e("Main", "Chain " + ligChainName + " of PDB " + pdbid + " not found in DB, cannot write LCG info to DB.");
                        continue;
                    }
                    Long sse_db_id = DBManager.getSSEDBID(chain_db_id, ligandSSE.getStartDsspNum());
                    if(sse_db_id <= 0L) { 
                        DP.getInstance().e("Main", "SSE with DSSP start " + ligandSSE.getStartDsspNum() + " of PDB " + pdbid + " chain " + ligChainName + " not found in DB, cannot write LCG info to DB.");
                        continue;
                    }
                    
                    Boolean graphOK = false;
                    try {
                        DBManager.writeLigandCenteredComplexGraphToDB(pdbid, sse_db_id, dbImagePathLCG + ".svg", dbImagePathLCG + ".png", dbImagePathLCG + ".pdf");
                        graphOK = true;
                        if(! silent) {
                            System.out.println("      Wrote ligand-centered complex graph of ligand " + ligName + " (PDB " + pdbid + ") to DB.");
                        }
                    } catch(SQLException e) {
                        DP.getInstance().w("Main", "Could not write ligand-centered omplex graph to DB: '" + e.getMessage() + "'.");
                    }
                    
                    // assign all chains to the graph in the DB (only makes sense if graph was entered)
                    if(graphOK) {
                        Long lcg_db_id = DBManager.getDBLCGID(sse_db_id);
                        if(lcg_db_id <= 0L) { 
                            DP.getInstance().e("Main", "Ligand-centered complex graph of ligand " + ligName + ", identified by database SSE ID " + sse_db_id + ", not found in DB, cannot assign LCG to chain in DB.");
                            continue;
                        }
                        
                        for(String contactChainName : ligContactChains) {
                            Long contact_chain_db_id = DBManager.getDBChainID(pdbid, contactChainName);
                            if(contact_chain_db_id <= 0L) { 
                                DP.getInstance().e("Main", "Chain " + contactChainName + " of PDB " + pdbid + " not found in DB, cannot assign LCG to chain " + contactChainName + " in DB.");
                                continue;
                            }
                            
                            try {
                                DBManager.assignLigandCenteredComplexGraphToChain(lcg_db_id, contact_chain_db_id);
                                if(! silent) {
                                    System.out.println("        Assigned ligand-centered complex graph of ligand " + ligName + " to chain " + contactChainName + " in DB.");
                                }
                            } catch(SQLException e) {
                                DP.getInstance().e("Main", "Could not assign ligand-centered complex graph of ligand " + ligName + " to chain " + ligChainName + " in DB: '" + e.getMessage() + "'.");
                            }
                        }
                    }
                } 
            }
        }            
        else {
            if(! silent) {
                System.out.println("    Not drawing ligand-centered complex graphs (" + graphType + ").");
            }
        }
                
        
        if(! silent) {
            System.out.println("Complex graph computation done.");
        }
                
        ProteinResults.getInstance().setCompGraphRes(cgr);
    }
    
    
    /**
     * Retrieves a list of all residues from a list of molecules.
     * @param molecules ArrayList of Molecule
     * @return ArrayList of Residue
     */
    private static ArrayList<Residue> resFromMolecules(ArrayList<Molecule> molecules) {
        if (residues == null) {
            residues = new ArrayList<>();
            for (Molecule m : molecules) {
                if (m instanceof Residue) {
                    residues.add((Residue) m);
                }
            }
            // check if empty
            if (residues.isEmpty()) {
                System.out.println("Detected no residues within molecules when first called resFromMolecules. Program will rely on that from now on.");
            }
        }
        return residues;
    }
    
    
    /**
     * Retrieves a list of all RNAs from a list of molecules.
     * @param molecules ArrayList of Molecule
     * @return ArrayList of Molecule
     */
    private static ArrayList<RNA> rnaFromMolecules(ArrayList<Molecule> molecules) {
        if (rnas == null) {
            rnas = new ArrayList<>();
            for (Molecule m : molecules) {
                if (m instanceof RNA) {
                    rnas.add((RNA) m);
                }
            }
            // check if empty
            if (rnas.isEmpty()) {
                System.out.println("Detected no RNA within molecules when first called rnaFromMolecules. Program will rely on that from now on.");
            }
        }
        return rnas;
    }
    
}
