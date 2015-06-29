<?php
/** 
 * @author Tim Schäfer
 */

ini_set('display_errors', 0);
ini_set('display_startup_errors', 0);
ini_set('log_errors', TRUE);
error_reporting(E_ERROR);

$valid_values = FALSE;
// get config values
include('./backend/config.php'); 

if($DEBUG){
	ini_set('display_errors', 1);
	ini_set('display_startup_errors', 1);
	ini_set('log_errors', TRUE);
	error_reporting(E_ALL);
}


function get_protein_graph_file_name_no_ext($pdbid, $chain, $graphtype_string) {
  return $pdbid . "_" . $chain . "_" . $graphtype_string . "_PG";
}


function get_protein_graph_path_and_file_name_no_ext($pdbid, $chain, $graphtype_string) {
  $path = get_path_to($pdbid, $chain);
  $fname = get_protein_graph_file_name_no_ext($pdbid, $chain, $graphtype_string);
  return $path . $fname;
}


function get_graphtype_string($graphtype_int){
	switch ($graphtype_int){
		case 1:
			return "alpha";
			break;
		case 2:
			return "beta";
			break;
		case 3:
			return "albe";
			break;
		case 4:
			return "alphalig";
			break;
		case 5:
			return "betalig";
			break;
		case 6:
			return "albelig";
			break;
	}
}

function check_valid_pdbid($str) {
  if (preg_match('/^[A-Z0-9]{4}$/i', $str)) {
    return true;
  }  
  return false;
}

function check_valid_chainid($str) {
  if (preg_match('/^[A-Z0-9]{1}$/i', $str)) {
    return true;
  }  
  return false;
}

function get_folding_graph_file_name_no_ext($pdbid, $chain, $graphtype_string, $fg_number) {
  return $pdbid . "_" . $chain . "_" . $graphtype_string . "_FG_" . $fg_number;
}

function get_path_to($pdbid, $chain) {
  $mid2chars = substr($pdbid, 1, 2);
  return $mid2chars . "/" . $pdbid . "/". $chain . "/";
}

function get_folding_graph_path_and_file_name_no_ext($pdbid, $chain, $graphtype_string, $fg_number) {
  $path = get_path_to($pdbid, $chain);
  $fname = get_folding_graph_file_name_no_ext($pdbid, $chain, $graphtype_string, $fg_number);
  return $path . $fname;
}

$pageload_was_search = FALSE;
$valid_values = FALSE;

if(isset($_GET['first_pdbchain']) && isset($_GET['first_graphtype_int']) &&  isset($_GET['second_pdbchain']) && isset($_GET['second_graphtype_int'])){
        $pageload_was_search = TRUE;
	$valid_values = FALSE;
	$first_pdbchain = $_GET["first_pdbchain"];
	$second_pdbchain = $_GET["second_pdbchain"];
	$first_graphtype_int = $_GET["first_graphtype_int"];
	$second_graphtype_int = $_GET["second_graphtype_int"];

	$valid_values_first = FALSE;
	if($first_graphtype_int === "1" || $first_graphtype_int === "2" || $first_graphtype_int === "3" || 
	    $first_graphtype_int === "4" || $first_graphtype_int === "5" || $first_graphtype_int === "6") 
	    { 

			if(strlen($first_pdbchain) === 5) {
				$first_pdb_id = substr($first_pdbchain, 0, 4);
				$first_chain_name = substr($first_pdbchain, 4, 1);						
			
				if(check_valid_pdbid($first_pdb_id) && check_valid_chainid($first_chain_name)) {
					$valid_values_first = TRUE;				  
				}
			}					
	}
	
	$valid_values_second = FALSE;
	if($second_graphtype_int === "1" || $second_graphtype_int === "2" || $second_graphtype_int === "3" || 
	    $second_graphtype_int === "4" || $second_graphtype_int === "5" || $second_graphtype_int === "6") 
	    { 

			if(strlen($second_pdbchain) === 5) {
				$second_pdb_id = substr($second_pdbchain, 0, 4);
				$second_chain_name = substr($second_pdbchain, 4, 1);						
			
				if(check_valid_pdbid($second_pdb_id) && check_valid_chainid($second_chain_name)) {
					$valid_values_second = TRUE;				  
				}
			}					
	}
	
	$valid_values = ($valid_values_first && $valid_values_second);
}

$num_found = 0;
$tableString = "";
$gml_files_available = FALSE;

if($valid_values){    	
	// check for the GML files
	$first_graphtype_str = get_graphtype_string($first_graphtype_int);	
	$first_graph_file_name_no_ext = get_protein_graph_path_and_file_name_no_ext($first_pdb_id, $first_chain_name, $first_graphtype_str);
	$first_full_file = $IMG_ROOT_PATH . $first_graph_file_name_no_ext . ".gml";
	

	$second_graphtype_str = get_graphtype_string($second_graphtype_int);	
	$second_graph_file_name_no_ext = get_protein_graph_path_and_file_name_no_ext($second_pdb_id, $second_chain_name, $second_graphtype_str);
	$second_full_file = $IMG_ROOT_PATH . $second_graph_file_name_no_ext . ".gml";
	if(file_exists($first_full_file) && file_exists($second_full_file)){
	    $gml_files_available = TRUE;
	    //$tableString .= "GML files found, running bk.<br>";
		//$bk_random_tag = rand(10000, 99999);
		$bk_random_tag = 12345;
		$outfile_prefix = "bkrun" . $bk_random_tag . "_";
		
		//DEBUG:
		$outfile_prefix = "bkrun12345_";
		
	    $stdoutput = "";
		
		if(! is_dir('bk_web')) {
			echo "ERROR: The required bk_web directory was not found. Aborting.<br>";
			die();
		}
		if( ! chdir('bk_web')) {
		    die("Could not change into BK run directory. Aborting.");
		}
		if( ! is_writable("./tmp_output/")) {
			echo "ERROR: Temporary directory not writeable, exiting. The server admin needs to fix this. Aborting.<br>";
			die();
		}
		$first_full_file_changed_rel = '../' . $first_full_file;
		$second_full_file_changed_rel = '../' . $second_full_file;
		//exec("./bk_protsim $first_full_file_changed_rel $second_full_file_changed_rel -l -f --output-dir tmp_output --output-prefix $outfile_prefix", $stdoutput);
		if( ! chdir('..')) {
		    die("Could not leave BK run directory. Aborting.");
		}
		
		// now go and search for the result files generated by bk_protsim
		$results = array();
		$results_second = array();
        $directory = 'bk_web/tmp_output/';
        $handler = opendir($directory);
        while ($file = readdir($handler)) {
  
          if ($file != "." && $file != "..") {
		    // find file names like 'bkrunXXXXX_results_0_first.txt', where the 0 can be any number and XXXXX is a random number prefix.
            if(preg_match('#^' . $outfile_prefix . '(results_)[^\s]*_first(\.(txt))#', $file)) {
              // add to our file array for later use
              $results[] = $file;
            }
			
			// find file names like 'bkrunXXXXX_results_0_second.txt', where the 0 can be any number and XXXXX is a random number prefix.
			if(preg_match('#^' . $outfile_prefix . '(results_)[^\s]*_second(\.(txt))#', $file)) {
              // add to our file array for later use
              $results_second[] = $file;
            }
			
           }
        }
		$num_found = count($results);
		$num_found_second = count($results_second);
		
		// TODO: we could check whether both arrays have the same length, and throw an error if not
		
		$tableString .= "<form  method='POST' action='bkweb_vis.php' target='_blank'>";
		$tableString .= "<table id='bktable'><tr><th>Result ID</th><th>File names</th><th>Preview</th><th>Substructure size</th><th>Select substructure</th></tr>";

		for($i = 0; $i < count($results); $i++) {
		  $file_first = $results[$i];
		  $file_second = $results_second[$i];
		  $num_lines = count(file('bk_web/tmp_output/' . $file_first));
		  $result_id = $i; // TODO: should parse from file name instead of simply incrementing!
		  $checked_string = "";
		  $file_content_first = file_get_contents('bk_web/tmp_output/' . $file_first);
		  $file_content_second = file_get_contents('bk_web/tmp_output/' . $file_second);
		  $preview_first = trim(preg_replace('/\s+/', ' ', $file_content_first));  // replace newlines with space
		  $preview_second = trim(preg_replace('/\s+/', ' ', $file_content_second));  // replace newlines with space
		  $preview = $preview_first . "<br>" . $preview_second;
		  if($i === 0) { $checked_string = " checked='checked'"; }
		  $tableString .= "<tr><td>$result_id</td><td>$file_first <br> $file_second</td><td>$preview</td><td>$num_lines</td><td><input type='radio' name='result_id' value='$result_id' $checked_string><br></td></tr>";
		}
		$tableString .= "</table>";		
		//$tableString .= "<input type='hidden' name='gml_file_first' value='$file_first'>";
		//$tableString .= "<input type='hidden' name='gml_file_second' value='$file_second'>";
		$tableString .= "<input type='hidden' name='first_graphtype_int' value='$first_graphtype_int'>";
		$tableString .= "<input type='hidden' name='second_graphtype_int' value='$second_graphtype_int'>";
		$tableString .= "<input type='hidden' name='first_pdbchain' value='$first_pdbchain'>";
		$tableString .= "<input type='hidden' name='second_pdbchain' value='$second_pdbchain'>";
		$tableString .= "<input type='hidden' name='bk_random_tag' value='$bk_random_tag'>";		
		
		$tableString .= "<button type='submit' id='sendit_visualize' name='sendit_visualize' class='btn btn-default' style='margin-top:35px;'><span>Visualize selected substructure</span></button>";
		$tableString .= "</form>";
	}
	else {
	    $tableString .= "GML files for PDB chains not available ($first_full_file and $second_full_file).";
	}	

	

}

//EOF
?>