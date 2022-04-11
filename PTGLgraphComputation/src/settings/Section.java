/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Jan Niclas Wolf 2020. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author jnw
 */

package settings;

import java.util.ArrayList;
import java.util.HashMap;
import tools.DP;
import proteingraphs.ComplexGraphEdgeWeightTypes;

/**
 *
 * @author jnw
 */
class Section {
    
    String name;  // name of the section
    String type;  // either user, advanced or developer
    ArrayList<Setting> settings;  // the settings that are part of this section
    // WARNING: Following data structure MUST be filled after settings have been filled
    HashMap<String, Integer> mapSettingNameToSettingIndex;  // Maps the setting names (=key) to their index in the ArrayList
    
    static final ArrayList<String> ALLOWED_TYPES = new ArrayList<String>() {
        {
            add("User");
            add("Advanced");
            add("Developer");
        }
    };
    
    /**
     * Creates a new settings section-
     * @param name
     * @param type either user, advanced or developer
     */
    public Section(String name, String type) {
        this.name = name;
        
        // check that data type is allowed to avoid sloppy programming
        if (ALLOWED_TYPES.contains(type)) {
            this.type = type;
        } else {
            DP.getInstance().e(Settings.PACKAGE_TAG, "Creating the setting '" + name + " failed, because the data type is not allowed or correctly formatted. "
                    + "This is an error in the code, please inform the developer of this software. Exiting now.");
            this.type = "user";  // for the IDE
            System.exit(1);
        }
        
        this.type = type;
        
        initSettings();
        initIndexMap();
    }
    
    
    /**
     * Initializes a settings section, i.e., creates the corresponding settings.
     * Should only be called when the name of the section is set.
     */
    private void initSettings() {
        settings = new ArrayList<>();
        switch(name) {
            case "General settings":
                settings.add(new Setting("PTGLgraphComputation_B_use_mmCIF_parser", 'B', "true", "Whether the mmCIF parser should be used on the provided PDB coordinates file."));
                settings.add(new Setting("PTGLgraphComputation_S_output_dir", 'S', ".", "Output directory for all created files."));
                settings.add(new Setting("PTGLgraphComputation_B_calc_draw_graphs", 'B', "true", "Whether graphs are computed and output."));
                settings.add(new Setting("PTGLgraphComputation_B_draw_graphs", 'B', "true", "Whether graph visualizations are output."));
                settings.add(new Setting("PTGLgraphComputation_B_force_chain", 'B', "false", "Whether to force parsing and processing only a certain PDB chain."));
                settings.add(new Setting("PTGLgraphComputation_S_forced_chain_id", 'S', "A", "The forced chain ID, only used when PTGLgraphComputation_B_force_chain is true."));
                break;
                
            case "Amino Acid Graphs (AAG)":
                settings.add(new Setting("PTGLgraphComputation_B_draw_aag", 'B', "false", "Whether Amino Acid Graph visualizations are output. Uses grid visualization."));
                settings.add(new Setting("PTGLgraphComputation_B_aminoacidgraphs_include_ligands", 'B', "false", "Whether Amino Acid Graphs should include ligands."));
                settings.add(new Setting("PTGLgraphComputation_B_quit_after_aag", 'B', "false", "Whether to quit the program after computation of Amino Acid Graphs."));
                settings.add(new Setting("PTGLgraphComputation_B_img_AAG_output_format_PNG", 'B', "false", "Whether to write Amino Acid Graph visualizations in PNG format."));
                settings.add(new Setting("PTGLgraphComputation_B_img_AAG_output_format_PDF", 'B', "false", "Whether to write Amino Acid graph visualizations in PDF format."));
                // Removed: Currently not used as SVG is base for conversion to other file formats and therefore always created.
                //settings.add(new Setting("PTGLgraphComputation_B_img_AAG_output_format_SVG", 'B', "true", "Whether to write amino acid graph visualizations in PDF format. Note that this setting currently has no effect, SVG is always generated. The other formats get converted from the SVG."));
                break;
            
            case "Protein Graph (PG)":
                settings.add(new Setting("PTGLgraphComputation_B_skip_empty_chains", 'B', "true", "Whether to skip chains without any DSSP SSEs (i.e., contain only ligands) in Protein Graphs."));
                settings.add(new Setting("PTGLgraphComputation_B_img_output_format_PNG", 'B', "true", "Whether to write Protein Graph visualizations in PNG format."));
                settings.add(new Setting("PTGLgraphComputation_B_img_output_format_PDF", 'B', "false", "Whether to write Protein Graph visualizations in PDF format."));
                // Removed: Currently not used as SVG is base for conversion to other file formats and therefore always created.
                //settings.add(new Setting("PTGLgraphComputation_B_img_output_format_SVG", 'B', "true", "Whether to write protein graph visualizations in SVG format. Note that this setting currently has no effect, SVG is always generated. The other formats get converted from the SVG."));
                break;
            
            case "Graph types":
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_albe", 'B', "true", "Whether alpha-beta graphs are output."));
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_albelig", 'B', "true", "Whether alpha-beta graphs with ligands are output."));
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_alpha", 'B', "true", "Whether alpha graphs are output."));
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_alphalig", 'B', "true", "Whether alpha graphs with ligands are output."));
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_beta", 'B', "true", "Whether beta graphs are output."));
                settings.add(new Setting("PTGLgraphComputation_B_graphtype_betalig", 'B', "true", "Whether beta graphs with ligands are output."));
                break;
                
            case "Folding Graph (FG)":
                settings.add(new Setting("PTGLgraphComputation_B_folding_graphs", 'B', "false", "Whether Folding Graphs (connected components of a Protein Graph) are computed. This does NOT mean they are drawn."));
                settings.add(new Setting("PTGLgraphComputation_B_draw_folding_graphs", 'B', "true", "Whether Folding Graphs visualizations are output."));
                settings.add(new Setting("PTGLgraphComputation_B_img_FG_output_format_PNG", 'B', "true", "Whether to write Folding Graph visualizations in PNG format."));
                settings.add(new Setting("PTGLgraphComputation_B_img_FG_output_format_PDF", 'B', "false", "Whether to write Folding Graph visualizations in PDF format."));
                // Removed: Currently not used as SVG is base for conversion to other file formats and therefore always created.
                //settings.add(new Setting("PTGLgraphComputation_B_img_FG_output_format_SVG", 'B', "true", "Whether to write Folding Graph visualizations in PNG format. Note that this setting currently has no effect, SVG is always generated. The other formats get converted from the SVG."));
                settings.add(new Setting("PTGLgraphComputation_B_foldgraphtype_KEY", 'B', "true", "Whether KEY notation of Folding Graphs is computed and output (only applies if 'PTGLgraphComputation_B_folding_graphs' is 'true')."));
                settings.add(new Setting("PTGLgraphComputation_B_foldgraphtype_ADJ", 'B', "true", "Whether ADJ notation of Folding Graphs is computed and output (only applies if 'PTGLgraphComputation_B_folding_graphs' is 'true')."));
                settings.add(new Setting("PTGLgraphComputation_B_foldgraphtype_RED", 'B', "true", "Whether RED notation of Folding Graphs is computed and output (only applies if 'PTGLgraphComputation_B_folding_graphs' is 'true')."));
                settings.add(new Setting("PTGLgraphComputation_B_foldgraphtype_SEQ", 'B', "true", "Whether SEQ notation of Folding Graphs is computed and output (only applies if 'PTGLgraphComputation_B_folding_graphs' is 'true')."));
                settings.add(new Setting("PTGLgraphComputation_B_foldgraphtype_DEF", 'B', "false", "Whether DEF notation of Folding Graphs is computed and output (only applies if 'PTGLgraphComputation_B_folding_graphs' is 'true')."));
                break;
                
            case "Complex Graph (CG)":
                settings.add(new Setting("PTGLgraphComputation_B_complex_graphs", 'B', "false", "Whether Complex Graphs are computed and drawn."));
                settings.add(new Setting("PTGLgraphComputation_I_CG_contact_threshold", 'I', "1", "The lowest number of interchain residue contacts where an edge in Complex Graphs is drawn."));
                settings.add(new Setting("PTGLgraphComputation_B_CG_ignore_ligands", 'B', "true", "Whether ligands should be ignored. If not, they may be part of the chain they are "
                        + "added to by the authors of the PDB file. This influences the (number of) contacts. May be superseded in the future by placing ligands in "
                        + "own chain."));
                settings.add(new Setting("PTGLgraphComputation_B_img_CG_output_format_PNG", 'B', "true", "Whether to write Complex Graph visualizations in PNG format."));
                settings.add(new Setting("PTGLgraphComputation_B_img_CG_output_format_PDF", 'B', "false", "Whether to write Complex Graph visualizations in PDF format."));
                // Removed: Currently not used as SVG is base for conversion to other file formats and therefore always created.
                //settings.add(new Setting("PTGLgraphComputation_B_img_CG_output_format_SVG", 'B', "true", "Whether to write Complex Graph visualizations in SVG format.  Note that this setting currently has no effect, SVG is always generated. The other formats get converted from the SVG."));
                break;
                
            case "Structure visualization":
                settings.add(new Setting("PTGLgraphComputation_B_Jmol_graph_vis_commands", 'B', "true", "Whether to compute and print Jmol commands to visualize Protein Graphs in 3D."));
                settings.add(new Setting("PTGLgraphComputation_B_Jmol_graph_vis_resblue_commands", 'B', "false", "Whether to compute and print Jmol commands to color the residues of Protein Graphs blue in 3D."));
                break;
                
            case "Assembly prediction (AP)":
                settings.add(new Setting("PTGLgraphComputation_B_interactive_assembly_prediction", 'B', "false", "Whether to interactively create assembly prediction instead of greedily choosing best edge at each step."));
                settings.add(new Setting("PTGLgraphComputation_S_assembly_prediction_leaf_labels", 'S', "", "Comma-seperated labels for the leaves in the interactive prediction and output newick format. "
                        + "The labels are applied in the same order as in the PDB file. Will use chain IDs if not specified."));
                settings.add(new Setting("PTGLgraphComputation_L_AP_weight_types", 'L', String.join(",", ComplexGraphEdgeWeightTypes.getShortNames()), "Comma-seperated list of the weight types to use for the assembly prediction. "
                        + "For each weight type, a prediction is done. Valid choices are: " + String.join(",", ComplexGraphEdgeWeightTypes.getShortNames()) + "."));
                break;
                
            case "Prints / Error handling":
                settings.add(new Setting("PTGLgraphComputation_B_silent", 'B', "false", "Whether to suppress all output. Cluster mode, not recommended for normal usage."));
                settings.add(new Setting("PTGLgraphComputation_B_only_essential_output", 'B', "true", "Whether to print only high-level status information."));
                settings.add(new Setting("PTGLgraphComputation_B_print_contacts", 'B', "false", "Whether to print residue contacts to stdout (slower)"));
                settings.add(new Setting("PTGLgraphComputation_B_no_warn", 'B', "false", "Whether to suppress all warnings."));
                settings.add(new Setting("PTGLgraphComputation_B_no_parse_warn", 'B', "true", "Whether to suppress all warnings related to parsing of atoms and other data from the PDB and DSSP input files."));
                settings.add(new Setting("PTGLgraphComputation_B_no_not_impl_warn", 'B', "true", "Whether to suppress all warnings related to not implemented function."));
                settings.add(new Setting("PTGLgraphComputation_B_no_chain_break_info", 'B', "false", "Whether to suppress chain break info while parsing DSSP file (handy for some DSSP files for CIF data)."));
                settings.add(new Setting("PTGLgraphComputation_B_print_silent_notice", 'B', "true", "Whether to print a single line informing the user that silent mode is set in silent mode (includes PDB ID of current file)."));
                settings.add(new Setting("PTGLgraphComputation_B_warn_cfg_fallback_to_default", 'B', "true", "Whether to print warnings when a setting is not defined in the config file and internal defaults are used."));
                settings.add(new Setting("PTGLgraphComputation_B_split_dsspfile_warning", 'B', "false", "Whether to show a warning about splitting the DSSP file when multiple models are detected in a PDB file."));
                break;
                
            case "Performance":
                settings.add(new Setting("PTGLgraphComputation_B_separate_contacts_by_chain", 'B', "false", "Whether to compute atom contacts separated by chain. Faster but does not detect contacts between different chains and thus cannot be used for Complex Graph computation."));
                settings.add(new Setting("PTGLgraphComputation_F_abort_if_pdb_resolution_worse_than", 'F', "10.0", "Abort all processing of PDB files with worse resolution than provided "
                        + "with this setting in Angstroem as a float. If set to 20.0, the program will terminate and NOT further process a PDB file which is detected to have a resolution of more than 20 A. "
                        + "Set to a negative value like -1.0 to disable this, and thus parse all PDB files, no matter the resolution."));
                settings.add(new Setting("PTGLgraphComputation_I_abort_if_num_molecules_below", 'I', "30", "Abort all processing of PDB files with too few molecules as defined by the given value of the setting as integer. "
                        + "If set to 30, the program will terminate and NOT further process a PDB file which is detected to have less than 30 molecules. "
                        + "Set to a negative value like -1 to disable this, and thus parse all PDB files, no matter the molecule count."));
                settings.add(new Setting("PTGLgraphComputation_B_skip_too_large", 'B', "false", "Whether to abort if the protein has more than 'PTGLgraphComputation_I_skip_num_atoms_threshold' atoms."));
                settings.add(new Setting("PTGLgraphComputation_I_skip_num_atoms_threshold", 'I', "80000", "The maximal number of atoms per PDB file if 'PTGLgraphComputation_B_skip_too_large' is true. In that case, PTGLgraphComputation will abort for PDB files with more atoms."));
                settings.add(new Setting("PTGLgraphComputation_B_chain_spheres_speedup", 'B', "true", "Whether to use contact computation speedup based on comparison of chain spheres."));
                settings.add(new Setting("PTGLgraphComputation_B_centroid_method", 'B', "true", "Whether to use centroid of atoms instead of C_alpha for contact computation. Recommended use only with PTGLgraphComputation_B_chain_spheres_speedup."));
                settings.add(new Setting("PTGLgraphComputation_B_round_coordinates", 'B', "true", "Whether 3D atom coordinates should be rounded or truncated one decimal place."));
                settings.add(new Setting("PTGLgraphComputation_S_temp_dir", 'S', ".", "The directory where temporary files can be created. You need write access to it, of course."));
                break;
                
            case "Parser":
                settings.add(new Setting("PTGLgraphComputation_B_include_rna", 'B', "true", "Whether RNA should be parsed and included in graph formalism and visualisation (WIP)."));
                settings.add(new Setting("PTGLgraphComputation_B_convert_models_to_chains", 'B', "false", "Whether the PDB file should be checked for multiple models and if so convert those models to chains."));
                settings.add(new Setting("PTGLgraphComputation_I_defaultModel", 'I', "1", "The model to use if multiple models exist in the PDB file."));
                settings.add(new Setting("PTGLgraphComputation_B_always_skip_models", 'B', "true", "Whether models should be skipped for all methods."));
                settings.add(new Setting("PTGLgraphComputation_S_ligAACode", 'S', "J", "The amino acid code used to mark a ligand residue."));
                settings.add(new Setting("PTGLgraphComputation_I_lig_min_atoms", 'I', "1", "The minimum number of atoms a ligand needs to consist of to count as an SSE."));
                settings.add(new Setting("PTGLgraphComputation_I_lig_max_atoms", 'I', "-1", "The maximum number of atoms a ligand has to consist of to count as an SSE. Set to <0 for unlimited."));
                break;
                
            case "Secondary Structure Element (SSE)":
                settings.add(new Setting("PTGLgraphComputation_S_ligSSECode", 'S', "L", "The SSE code used to mark an SSE of type ligand."));
                settings.add(new Setting("PTGLgraphComputation_S_rnaSseCode", 'S', "R", "The SSE code used to mark an SSE of type RNA"));
                settings.add(new Setting("PTGLgraphComputation_I_min_SSE_length", 'I', "3", "The minimal length in AAs a non-ligand SSE must have to be considered (PTGL-style filtering of very short SSEs)"));
                settings.add(new Setting("PTGLgraphComputation_B_change_dssp_sse_b_to_e", 'B', "true", "Whether all isolated beta bridges (DSSP SSE type: B) should be changed to betastrands (DSSP SSE type: E)"));
                settings.add(new Setting("PTGLgraphComputation_B_fill_gaps", 'B', "true", "Whether two strands should be fused to one strand, when there is only one AA (with DSSP SSE ' ') between them."));
                settings.add(new Setting("PTGLgraphComputation_B_include_coils", 'B', "false", "Whether coils (DSSP SSE type ' ') should be considered as own vertices."));
                settings.add(new Setting("PTGLgraphComputation_S_coilSSECode", 'S', "C", "The amino acid code used to mark a coiled region residue."));
                settings.add(new Setting("PTGLgraphComputation_B_merge_helices", 'B', "true", "whether to merge different helix types if they are adjacent in the primary structure."));
                settings.add(new Setting("PTGLgraphComputation_I_merge_helices_max_dist", 'I', "0", "The maximal distance in amino acids of the primary sequence where helices are merged when 'PTGLgraphComputation_B_merge_helices' is true. "
                        + "The (default) value 0 means only directly adjacent SSEs are merged."));
                break;
                
            case "SSE orientation":
                settings.add(new Setting("PTGLgraphComputation_B_spatrel_use_dd", 'B', "false", "Whether to use double difference mode for computation of orientation of SSEs. False invokes vector-mode instead."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_vector_num_res_centroids", 'I', "4", "Vector mode: How many residues to use for centroid computation of start and end point of vector for SSE."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_max_deg_parallel", 'I', "65", "Vector mode: Degrees to which SSEs are classified to be parallel."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_min_deg_antip", 'I', "115", "Vector mode: Degrees from which SSEs are classified to be antiparallel in vector mode."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_largest_antip_ee", 'I', "0", "Double difference mode: All values <= the one given here are considered antiparallel. This is for E/E (strand/strand) interactions."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_smallest_parallel_ee", 'I', "1", "Double difference mode: All values >= the one given here are considered parallel. This is for E/E (strand/strand) interactions. Note that the range in between these 2 values is considered mixed (none in the case of E/E)."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_largest_antip_hh", 'I', "-8", "Double difference mode: Same as above, but for H/H (helix/helix) interactions."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_smallest_parallel_hh", 'I', "8", "Double difference mode: Same as above, but for H/H (helix/helix) interactions."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_largest_antip_he", 'I', "-6", "Double difference mode: Same as above, but for H/E or E/H interactions."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_smallest_parallel_he", 'I', "6", "Double difference mode: Same as above, but for H/E or E/H interactions."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_largest_antip_def", 'I', "-7", "Double difference mode: Same as above, this is the default for other interactions (e.g., coil/helix)."));
                settings.add(new Setting("PTGLgraphComputation_I_spatrel_dd_smallest_parallel_def", 'I', "7", "Double difference mode: Same as above, this is the default for other interactions."));
                break;
                
            case "Contact definition":
                settings.add(new Setting("PTGLgraphComputation_I_aa_atom_radius", 'I', "20", "The atom radius of protein atoms in 10th part Angstroem (setting 20 here means 2 A)"));
                settings.add(new Setting("PTGLgraphComputation_I_rna_atom_radius", 'I', "20", "The atom radius of RNA atoms in 10th part Angstroem (setting 20 here means 2 A)"));
                settings.add(new Setting("PTGLgraphComputation_I_lig_atom_radius", 'I', "30", "The atom radius of ligand atoms in 10th part Angstroem (setting 30 here means 3 A)"));
                settings.add(new Setting("PTGLgraphComputation_I_max_contacts_per_type", 'I', "1", "The maximum number of contacts of a certain type that is counted for a residue pair. Set it to something very large if you don't want any limit (Integer.MAX_VALUE comes to mind). "
                        + "The PTGL uses a setting of 1 (so if a pair has 3 B/B contacts and 2 C/B contacts, it is is counted as 1 B/B and 1 C/B.)"));
                settings.add(new Setting("PTGLgraphComputation_transitive_contacts", 'B', "true", "Whether contacts transitive by ligands should be considered for the contact calculation, and subsequently in the complex graphs."));
                settings.add(new Setting("PTGLgraphComputation_B_forceBackboneContacts", 'B', "false", "Whether all amino acids of a protein graph should be connected sequentially, from N to C terminus, with contacts of type backbone."));
                break;
                
            case "AAG":
                settings.add(new Setting("PTGLgraphComputation_I_aag_min_residue_seq_distance_for_contact", 'I', "0", "If set to a value greater than zero, only contacts from residues more than or exactly this value apart in the sequence can have a contact."
                        + "Can be used to force network of long-range contacts, ignoring backbone connections. Contacts between residues from different chains always are accepted. This will split the graph into many connected components."));
                settings.add(new Setting("PTGLgraphComputation_I_aag_max_residue_seq_distance_for_contact", 'I', "0", "If set to a value greater than zero, only contacts from residues less than or exactly this value apart in the sequence can for a contact. "
                        + "Can be used to force network of short-range contacts, ignoring long-range connections. Also ignores all contacts to residues from other chains."));
                settings.add(new Setting("PTGLgraphComputation_B_AAgraph_allchainscombined", 'B', "false", "Whether to compute and output combined amino acid graphs (one for all chains) as well."));
                settings.add(new Setting("PTGLgraphComputation_B_AAgraph_perchain", 'B', "false", "Whether to compute and output amino acid graphs per chain as well."));
                break;
                
            case "Database (DB) connection":
                settings.add(new Setting("PTGLgraphComputation_B_useDB", 'B', "false", "Whether to write any data to the PostgreSQL database."));
                settings.add(new Setting("PTGLgraphComputation_S_db_name", 'S', "vplg", "Database name"));
                settings.add(new Setting("PTGLgraphComputation_S_db_host", 'S', "127.0.0.1", "Hostname or IP of the DB server"));
                settings.add(new Setting("PTGLgraphComputation_I_db_port", 'I', "5432", "DB server port"));
                settings.add(new Setting("PTGLgraphComputation_S_db_username", 'S', "vplg", "DB username"));
                settings.add(new Setting("PTGLgraphComputation_S_db_password", 'S', "", "DB password (empty if local is TRUST for this user)"));
                settings.add(new Setting("PTGLgraphComputation_B_db_use_autocommit", 'B', "false", "Whether the DB connection gets set to autocommit when created."));
                break;
                
            case "DB settings":
                settings.add(new Setting("PTGLgraphComputation_S_graph_image_base_path", 'S', "/srv/www/htdocs/vplgweb/data/", "[DEPRECATED] The base path for graph visualizations for the database. The relative path to the path given here is used to locate the graph image on disk."));
                settings.add(new Setting("PTGLgraphComputation_B_db_use_batch_inserts", 'B', "false", "Whether inserts into the database should use batch mode instead of many single queries whenever possible. Only implemented for a few queries for which it may make sense."));
                settings.add(new Setting("PTGLgraphComputation_B_write_graphstrings_to_database_pg", 'B', "false", "[DEPRECATED] Whether to write the protein graph strings to the database in the different formats like XML, GML, TGF, etc. "
                        + "They are currently written to disk only by default, and NULL values are inserted into the DB, because we do not use the DB fields for anything atm."));
                settings.add(new Setting("PTGLgraphComputation_B_write_graphstrings_to_database_fg", 'B', "false", "[DEPRECATED] Whether to write the folding graph strings to the database in the different formats like XML, GML, TGF, etc. "
                        + "They are currently written to disk only by default, and NULL values are inserted into the DB, because we do not use the DB fields for anything atm."));
                settings.add(new Setting("PTGLgraphComputation_B_write_graphstrings_to_database_cg", 'B', "false", "[DEPRECATED] Whether to write the complex graph strings to the database in the different formats like XML, GML, TGF, etc. "
                        + "They are currently written to disk only by default, and NULL values are inserted into the DB, because we do not use the DB fields for anything atm."));
                settings.add(new Setting("PTGLgraphComputation_B_write_graphstrings_to_database_aag", 'B', "false", "[DEPRECATED] Whether to write the amino acid graph strings to the database in the different formats like XML, GML, TGF, etc. "
                        + "They are currently written to disk only by default, and NULL values are inserted into the DB, because we do not use the DB fields for anything atm."));
                settings.add(new Setting("PTGLgraphComputation_B_report_db_proteins", 'B', "false", "Write a list of PDB IDs which are currently in the database to the file 'db_content_proteins.txt', then exit."));
                settings.add(new Setting("PTGLgraphComputation_B_uglySQLhacks", 'B', "false", "Whether to rewrite string data like ligand atom formulas before inserting into DB."));
                break;
                
            case "DB structure search":
                settings.add(new Setting("PTGLgraphComputation_B_matrix_structure_search", 'B', "false", "Search a linear notation (input) in a proteinstructure (input)"));
                settings.add(new Setting("PTGLgraphComputation_B_matrix_structure_search_db", 'B', "false", "Search a linear notation (input) in the database, ignores any given proteinstructure as input"));
                settings.add(new Setting("PTGLgraphComputation_S_linear_notation", 'S', "", "The linear notation of a PG for the matrix structure comparison."));
                settings.add(new Setting("PTGLgraphComputation_S_linear_notation_type", 'S', "", "The type of linear notation (adj or red) for matrix structure comparison."));
                settings.add(new Setting("PTGLgraphComputation_S_linear_notation_graph_type", 'S', "", "alpha, beta or albe = The graph type of the PG graph of the linear notation."));
                break;
                
            case "FG":
                settings.add(new Setting("PTGLgraphComputation_I_min_fgraph_size_draw", 'I', "3", "The minimum size of a Folding Graph to be drawn. Setting this to 1 or 0 draws all of them, including isolated vertices."));
                settings.add(new Setting("PTGLgraphComputation_I_min_fgraph_size_write_to_db", 'I', "1", "The minimum size of a Folding Graph to be written to the database. Setting this to 1 or 0 saves all of them, including isolated vertices."));
                settings.add(new Setting("PTGLgraphComputation_I_min_fgraph_size_write_to_file", 'I', "3", "The minimum size of a Folding Graph to be written to a file in formats like GML. Setting this to 1 or 0 exports all of them to files, including isolated vertices."));
                settings.add(new Setting("PTGLgraphComputation_B_key_use_alternate_arcs", 'B', "true", "Whether to use alternative crossover arcs in KEY notation. The alternative arcs cut through other SSEs, the default ones use a vertical central line and shift the center to avoid this."));
                settings.add(new Setting("PTGLgraphComputation_B_print_notations_on_fg_images", 'B', "false", "Whether to add the notation string to the FG images."));
                settings.add(new Setting("PTGLgraphComputation_B_output_fg_linear_notations_to_file", 'B', "true", "Whether to save linear notations of Folding Graphs to a text file."));
                settings.add(new Setting("PTGLgraphComputation_B_special_linnot_rules_for_bifurcated_adj_and_red", 'B', "true", "Whether special rules should be used for computing the ADJ and RED notations of FGs. These special rules were not used in the latest PTGL version by PM, but it seems in older versions."));
                break;
                
            case "PG":
                settings.add(new Setting("PTGLgraphComputation_B_output_GML", 'B', "true", "Whether to save Protein Graphs to text files in Graph Modelling Language (GML) format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_TGF", 'B', "false", "Whether to save Protein Graphs to text files in Trivial Graph Format (TGF)."));
                settings.add(new Setting("PTGLgraphComputation_B_output_DOT", 'B', "false", "Whether to save Protein Graphs to text files in DOT language format)."));
                settings.add(new Setting("PTGLgraphComputation_B_output_kavosh", 'B', "false", "Whether to save Protein Graphs to text files in Kavosh format."));
                settings.add(new Setting("PTGLgraphComputation_B_kavosh_format_directed", 'B', "true", "Whether to treat the graphs as directed for the Kavosh output. If set to true, each edge (a, b) will appear twice in the output file: once as (a, b) and again as (b, a)."));
                settings.add(new Setting("PTGLgraphComputation_B_output_eld", 'B', "false", "Whether to save Protein Graphs to text files in edge list format with a vertex type list file."));
                settings.add(new Setting("PTGLgraphComputation_B_output_plcc", 'B', "false", "Whether to save Protein Graphs to text files in PTGLgraphComputation format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_perlfg", 'B', "false", "Whether to save Protein Graphs to text files in the PTGL format used by the Perl script to compute Folding Graph notations."));
                settings.add(new Setting("PTGLgraphComputation_B_output_json", 'B', "true", "Whether to save Protein Graphs to text files in JSON format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_xml", 'B', "true", "Whether to save Protein Graphs to text files in XML format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_gexf", 'B', "false", "Whether to save Protein Graphs to text files in GEXF format."));
                settings.add(new Setting("PTGLgraphComputation_S_gexf_format_version", 'S', "1.1", "The version of the GEFX file format to use. Supported are '1.1' and '1.2'."));
                settings.add(new Setting("PTGLgraphComputation_B_output_msvg", 'B', "true", "Whether to save Protein Graphs to text files in SVG format for interactive mode."));
                settings.add(new Setting("PTGLgraphComputation_B_output_cytoscapejs", 'B', "false", "Whether to save Protein Graphs to text files in CytoscapeJS format."));
                settings.add(new Setting("PTGLgraphComputation_S_graph_metadata_splitstring", 'S', "|", "The field separator used when writing meta data to exported graphs."));
                settings.add(new Setting("PTGLgraphComputation_B_add_metadata_comments_GML", 'B', "false", "Whether to add meta data to exported GML format graphs in comments. Note: Not all programs parse comments correctly."));
                settings.add(new Setting("PTGLgraphComputation_B_add_metadata_comments_DOT", 'B', "false", "Whether to add meta data to exported DOT format graphs in comments. Note: Not all programs parse comments correctly."));
                settings.add(new Setting("PTGLgraphComputation_B_compute_graph_metrics", 'B', "false", "Whether to compute graph metrics such as cluster coefficient for PGs. Slower!"));
                break;
                
            case "CG":
                settings.add(new Setting("PTGLgraphComputation_B_draw_ligandcomplexgraphs", 'B', "false", "[EXPERIMENTAL] Whether to draw ligand-centered Complex Graphs. This means one graph for each ligand in a PDB file, each showing all SSEs (regardless of chain) the ligand is in contact with."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_GML", 'B', "true", "Whether to save Complex Graphs (including SSE info) to a text file in Graph Modelling Language (GML) format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_TGF", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in Trivial Graph Format (TGF) format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_DOT", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in DOT format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_kavosh", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in Kavosh format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_eld", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in edge list format with a vertex type list file."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_plcc", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in PTGLgraphComputation format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_JSON", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in JSON format."));
                settings.add(new Setting("PTGLgraphComputation_B_output_compgraph_XML", 'B', "false", "Whether to save Complex Graphs (including SSE info) to a text file in XGMML format."));
                settings.add(new Setting("PTGLgraphComputation_B_writeComplexContactCSV", 'B', "false", "Whether to write a CSV file containing all contacts used for Complex Graph computation."));
                break;
                
            case "Motifs":
                settings.add(new Setting("PTGLgraphComputation_B_compute_motifs", 'B', "true", "Whether to search for motifs in the FG linear notations (after computing and writing them to the DB)."));
                settings.add(new Setting("PTGLgraphComputation_B_compute_tim", 'B', "false", "Whether to search for tim barrel motif in FGs if PTGLgraphComputation_B_compute_motifs is true."));
                break;
                
            case "Cluster mode":
                settings.add(new Setting("PTGLgraphComputation_B_clustermode", 'B', "false", "Whether to write extra output files used only in cluster mode, like GML albe graph file list."));
                settings.add(new Setting("PTGLgraphComputation_B_output_images_dir_tree", 'B', "false", "Whether to write output images to a PDB-style sub directory structure under the output directory instead of writing them in there directly. "
                        + "This is useful if you want to process the whole PDB because most filesystems will get into trouble with tens of thousands of files in a single directory. "
                        + "The directory structure will be chosen from the meta data, i.e., PDB ID, chain, graph type, etc."));
                settings.add(new Setting("PTGLgraphComputation_B_output_textfiles_dir_tree", 'B', "false", "Whether to write output graph text files to a PDB-style sub directory structure under the output directory instead of writing them in there directly. "
                        + "This is useful if you want to process the whole PDB because most filesystems will get into trouble with tens of thousands of files in a single directory. "
                        + "The directory structure will be chosen from the meta data, i.e., PDB ID, chain, graph type, etc."));
                settings.add(new Setting("PTGLgraphComputation_B_output_textfiles_dir_tree_html", 'B', "false", "Whether to write HTML navigation files to the output directory tree. Only used if PTGLgraphComputation_B_output_textfiles_dir_tree is true as well."));
                settings.add(new Setting("PTGLgraphComputation_B_output_textfiles_dir_tree_core_html", 'B', "false", "Whether to write the core PTGLweb HTML files to the output directory tree. Only used if PTGLgraphComputation_B_output_textfiles_dir_tree is true as well. "
                        + "These are the main page, search form and other stuff which is only needed once for the whole website."));
                settings.add(new Setting("PTGLgraphComputation_B_html_add_complex_graph_data", 'B', "true", "Whether to write data on the complex graph to the protein result HTML webpage (if available)."));
                break;
                
            case "DB: representative chains":
                settings.add(new Setting("PTGLgraphComputation_B_set_pdb_representative_chains_pre", 'B', "false", "Whether this PTGLgraphComputation run should assign the representative PDB chains from the XML file in the info table of the database and then exit. Requires path to XML file. "
                        + "Should be run BEFORE filling on a therefore empty database."));
                settings.add(new Setting("PTGLgraphComputation_B_set_pdb_representative_chains_remove_old_labels_pre", 'B', "true", "Whether the old labels should be removed from all chains in the info table before the new ones are applied. Removed means all chains are considered NOT part of the representative set."));
                settings.add(new Setting("PTGLgraphComputation_B_set_pdb_representative_chains_post", 'B', "false", "Whether this PTGLgraphComputation run should assign the representative PDB chains from the XML file in the info table of the database and then exit. Requires path to XML file. "
                        + "Should be run AFTER filling on a therefore non-empty database."));
                settings.add(new Setting("PTGLgraphComputation_B_set_pdb_representative_chains_remove_old_labels_post", 'B', "true", "Whether the old labels should be removed from all chains in the chains table before the new ones are applied. Removed means all chains are considered NOT part of the representative set."));
                settings.add(new Setting("PTGLgraphComputation_S_representative_chains_xml_file", 'S', "representatives.xml", "The path to the XML file containing the representative PDB chains from the PDB. You can get the file from the RCSB PDB REST web service."));
                break;
                
            case "Output settings":
                settings.add(new Setting("PTGLgraphComputation_B_csv_number_residues_chains", 'B', "false", "Whether a csv file with the number of residues in each chain is created."));
                settings.add(new Setting("PTGLgraphComputation_B_csv_contacts_inter", 'B', "false", "Whether a csv file with the interchain residue-residue contacts is created."));
                settings.add(new Setting("PTGLgraphComputation_B_csv_contacts_intra_inter", 'B', "false", "Whether a csv file with the inter- and intrachain residue-residue contacts is created."));
                settings.add(new Setting("PTGLgraphComputation_B_ptgl_text_output", 'B', "false", "Whether the PTGL text files (e.g., those required by the bet_neo) are written. Not writing them is faster but this program cannot replace the PTGL tool 'geom_neo' anymore if this is deactivated."));
                settings.add(new Setting("PTGLgraphComputation_B_gml_snake_case", 'B', "false", "Whether keys in all GML files should be written in snake case, i.e., with under scores instead of camel case. HINT: Originally, GML does not support snake case and so may some parsers for GML files."));
                settings.add(new Setting("PTGLgraphComputation_B_ptgl_geodat_output", 'B', "false", "Whether the PTGL text files geo.dat for SSE level contacts are written to a text file."));
                settings.add(new Setting("PTGLgraphComputation_B_ramachandran_plot", 'B', "false", "Whether a Ramachandran plot is drawn to a file for each chain (slower)."));                
                break;
                
            case "Image settings":
                settings.add(new Setting("PTGLgraphComputation_B_graphimg_header", 'B', "true", "Whether the graph images are created with header line with info on the graph type, PDB and chain ID."));
                settings.add(new Setting("PTGLgraphComputation_B_graphimg_footer", 'B', "true", "Whether the graph images are created with footer line with info on the vertices."));
                settings.add(new Setting("PTGLgraphComputation_B_graphimg_legend", 'B', "true", "Whether the graph images are created with a legend that explains the color codes and SSE symbols. This is part of the footer."));
                settings.add(new Setting("PTGLgraphComputation_B_graphimg_legend_always_all", 'B', "false", "Whether the legend should contain all possible edge and vertex types, i.e., even those not occurring in the current image."));
                settings.add(new Setting("PTGLgraphComputation_B_graphimg_add_linnot_start_vertex", 'B', "false", "Whether the start vertex index of the linear notation in the parent graph is written to the image."));
                settings.add(new Setting("PTGLgraphComputation_I_img_margin_left", 'I', "80", "Size of the left image margin in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_margin_top", 'I', "40", "Size of the top image margin in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_margin_right", 'I', "40", "Size of the right image margin in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_margin_bottom", 'I', "40", "Size of the bottom image margin in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_vert_dist", 'I', "50", "The distance between two consecutive vertices in the output image, in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_vert_radius", 'I', "10", "The radius of a vertex in the output image, in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_header_height", 'I', "40", "The height of the header area in the output image, in pixels."));
                settings.add(new Setting("PTGLgraphComputation_I_img_footer_height", 'I', "180", "The height of the footer area in the output image, in pixels. The footer is used to print the legend."));
                settings.add(new Setting("PTGLgraphComputation_S_img_default_font", 'S', "TimesRoman", "The default font used in output image labels. This has to be a valid font name, of course."));
                settings.add(new Setting("PTGLgraphComputation_I_img_default_font_size", 'I', "16", "The default font size used in output images."));
                settings.add(new Setting("PTGLgraphComputation_I_img_legend_font_size", 'I', "16", "The legend font size used in output images."));
                settings.add(new Setting("PTGLgraphComputation_I_img_text_line_height", 'I', "40", "The vertical distance between two lines of text in the image, e.g., in the footer."));
                settings.add(new Setting("PTGLgraphComputation_I_img_min_img_height", 'I', "160", "The minimum size of the image area where the graph is drawn."));
                settings.add(new Setting("PTGLgraphComputation_I_img_min_arc_height", 'I', "100", "The minimum size of the arc area within the image area."));
                settings.add(new Setting("PTGLgraphComputation_I_img_minPageWidth", 'I', "800", "The minimum image width in pixels, used in output images."));
                settings.add(new Setting("PTGLgraphComputation_I_img_minPageHeight", 'I', "600", "The minimum image height in pixels, used in output images."));
                settings.add(new Setting("PTGLgraphComputation_S_img_output_fileext", 'S', ".png", "[DEPRECATED] File extension of images. Not used for graph images anymore; applies to Ramaplot etc only."));
                break;
                
            case "DB: similarity search":
                settings.add(new Setting("PTGLgraphComputation_B_search_similar", 'B', "false", "[EXPERIMENTAL] Whether to search for the most similar protein."));
                settings.add(new Setting("PTGLgraphComputation_S_search_similar_PDBID", 'S', "8icd", "Used only when PTGLgraphComputation_B_search_similar is true. The protein PDB ID to use as a pattern during the similarity search."));
                settings.add(new Setting("PTGLgraphComputation_S_search_similar_chainID", 'S', "A", "Used only when PTGLgraphComputation_B_search_similar is true. The protein chain ID to use as a pattern during the similarity search."));
                settings.add(new Setting("PTGLgraphComputation_S_search_similar_graphtype", 'S', "albelig", "Used only when PTGLgraphComputation_B_search_similar is true. The graph type to use as a pattern during the similarity search."));
                settings.add(new Setting("PTGLgraphComputation_I_search_similar_num_results", 'I', "5", "Used only when PTGLgraphComputation_B_search_similar is true. The number of results to print (e.g., 3 for the 3 most similar proteins in the DB)."));
                settings.add(new Setting("PTGLgraphComputation_S_search_similar_method", 'S', "string_sse", "Used only when PTGLgraphComputation_B_search_similar is true. The similarity measure to use, valid settings: string_sse, graph_set, graph_compat"));
                settings.add(new Setting("PTGLgraphComputation_B_compute_graphlet_similarities", 'B', "false", "Whether to compute graphlet similarities in the DB and then exit. It depends on other settings which graphlet sims are actually computed."));
                settings.add(new Setting("PTGLgraphComputation_S_search_similar_graphlet_scoretype", 'S', "RGF", "The method used to compute a similarity score from a pair of graphlet vectors. Valid options are: RGF=relative graphlet frequency distance, CUS=custom."));
                settings.add(new Setting("PTGLgraphComputation_B_compute_graphlet_similarities_pg", 'B', "false", "Whether to compute graphlet similarities for Protein Graphs. Only used if PTGLgraphComputation_B_compute_graphlet_similarities is true."));
                settings.add(new Setting("PTGLgraphComputation_B_compute_graphlet_similarities_cg", 'B', "false", "Whether to compute graphlet similarities for Complex Graphs. Only used if PTGLgraphComputation_B_compute_graphlet_similarities is true."));
                settings.add(new Setting("PTGLgraphComputation_B_compute_graphlet_similarities_aag", 'B', "false", "Whether to compute graphlet similarities for Amino Acid Graphs. Only used if PTGLgraphComputation_B_compute_graphlet_similarities is true."));
                settings.add(new Setting("PTGLgraphComputation_I_compute_all_graphlet_similarities_num_to_save_in_db", 'I', "25", "The number of the most similar protein chain to store in the database after graphlet similarity computation. Set to n to store the n most similar for each chain."));
                settings.add(new Setting("PTGLgraphComputation_I_compute_all_graphlet_similarities_start_graphlet_index", 'I', "0", "Determines the graphlets from the array in the DB which are considered for similarity computation. This is the index of the first (start) graphlet used. Do not forget to also set the end index properly. This is inclusive."));
                settings.add(new Setting("PTGLgraphComputation_I_compute_all_graphlet_similarities_end_graphlet_index", 'I', "29", "Determines the graphlets from the array in the DB which are considered for similarity computation. This is the index of the last (end) graphlet used. Do not forget to also set the start index properly. This is inclusive."));
                settings.add(new Setting("PTGLgraphComputation_I_number_of_graphlets", 'I', "30", "The length of the graphlet vector in the database (the PostgreSQL SQL array). This is the number of graphlets used to compute similarity."));
                settings.add(new Setting("PTGLgraphComputation_B_write_chains_file", 'B', "false", "Whether to write a chains file containing all chain names of the currently handled PDB file. Can be used by GraphletAnalyzer later to construct graph file names for all chains."));
                break;
                
            case "Random tree generation":
                settings.add(new Setting("PTGLgraphComputation_B_generate_random_binary_tree", 'B', "false", "Whether to generate a random binary tree."));
                settings.add(new Setting("PTGLgraphComputation_B_tree_rooted", 'B', "true", "Whether the generated tree should be rooted or not."));
                settings.add(new Setting("PTGLgraphComputation_I_tree_number_leafs", 'I', "4", "How many leafs the tree should contain."));
                settings.add(new Setting("PTGLgraphComputation_S_leaf_labels", 'S', "", "Comma-seperated labels for the leaves. Must exacly match the number of PTGLgraphComputation_I_tree_number_leafs. Will use IDs from 0..n if empty."));
                break;
                
            case "Debug":
                settings.add(new Setting("PTGLgraphComputation_I_debug_level", 'I', "0", "Debug level. Higher value means more output."));
                settings.add(new Setting("PTGLgraphComputation_B_debug_compareSSEContacts", 'B', "false", "Whether to compare the computed SSE level contacts to those in the geom_neo output file that is supplied."));
                settings.add(new Setting("PTGLgraphComputation_S_debug_compareSSEContactsFile", 'S', "geo.dat_ptgl", "The path to the geo.dat file to use for SSE level contact comparison."));
                settings.add(new Setting("PTGLgraphComputation_B_contact_debug_dysfunct", 'B', "false", "Atom level contact debugging mode. WARNING: When this is true, PTGLgraphComputation will abort after the first few residues and produce wrong overall results!"));
                settings.add(new Setting("PTGLgraphComputation_B_debug_only_parse", 'B', "false", "Exit after parsing. WARNING: When this is true, PTGLgraphComputation will abort after parsing and not produce results!"));
                settings.add(new Setting("PTGLgraphComputation_B_debug_only_contact_comp", 'B', "false", "Exit after contact computation. WARNING: When this is true, PTGLgraphComputation will abort after contact computation and not produce results!"));
                break;
                
            case "Disabled":
                settings.add(new Setting("PTGLgraphComputation_B_parse_binding_sites", 'B', "false", "[DISABLED] Whether to parse binding site data from the REMARK 800 and SITE lines of legacy PDB file."));
                settings.add(new Setting("PTGLgraphComputation_B_write_lig_geolig", 'B', "true", "[DISABLED] Determines whether ligand contacts are included in the <pdbid>.geolig file."));
                settings.add(new Setting("PTGLgraphComputation_B_consider_all_ligands_for_each_chain", 'B', "false", "[DISABLED] Whether to ignore the assignement of a ligand to a chain in the PDB file, and assign a ligand to each chain it has contacts with. WARNING: This setting is ignored and off atm."));
                settings.add(new Setting("PTGLgraphComputation_B_complex_graph_same", 'B', "false", "[DISABLED] Determines whether the complex graph is drawn with all nodes of the same type."));
                settings.add(new Setting("PTGLgraphComputation_B_complex_graph_mere", 'B', "false", "[DISABLED] Determines whether the complex graph is drawn with nodes of different type for each mere."));
                settings.add(new Setting("PTGLgraphComputation_I_ligSAS", 'I', "20", "[DISABLED] The solvent accessible surface value that is written to the dssplig file for ligands (not used atm)"));
                break;
                
            case "Alternate AAG":
                settings.add(new Setting("PTGLgraphComputation_B_handle_hydrogen_atoms_from_reduce", 'B', "false", "[EXPERIMENTAL] Whether to parse hydrogen atoms added to PDB files by the external Reduce software."));
                settings.add(new Setting("PTGLgraphComputation_B_alternate_aminoacid_contact_model", 'B', "false", "[EXPERIMENTAL] Whether to use alternate residue contact model by A. Scheck. Skips all computations except AA graphs."));
                settings.add(new Setting("PTGLgraphComputation_B_alternate_aminoacid_contact_model_with_ligands", 'B', "false", "[EXPERIMENTAL] Whether to use alternate residue contact model including ligands by A. Scheck. Skips all computations except AA graphs."));
                break;
                
            default:
                DP.getInstance().e(Settings.PACKAGE_TAG, "Settings section with name '" + name + "' not implemented. Please contact a developer.");
        }
    }
    
    
    /**
     * Creates a map for setting names to their index. Always call after settings are created!
     */
    private void initIndexMap() {
        mapSettingNameToSettingIndex = new HashMap<>();
        for (int i = 0; i < settings.size(); i++) {
            mapSettingNameToSettingIndex.put(settings.get(i).getName(), i);
        }
    }

    
    Setting getSetting(String name) {
        return settings.get(mapSettingNameToSettingIndex.get(name));
    }
    
    
    String asFormattedString() {
        String formattedString = "";
        formattedString += "## " + name + " ##\n\n";
        
        // settings
        for (Setting tmpSetting : settings) {
            formattedString += tmpSetting.asFormattedString() + "\n";
        }
        
        return formattedString;
    }
    
    
    /**
     * Returns each setting with its value in an own line.
     * @return 
     */
    String asShortString() {
        String shortString = "";
        for (Setting tmpSetting : settings) {
            shortString += tmpSetting.asKeyValueString() + "\n";
        }
        return shortString;
    }
}
