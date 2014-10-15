<?php
/** This file stores the configuration and the database connection information
 *  in an array in returns it.
 *   
 * @author Daniel Bruness <dbruness@gmail.com>
 * @author Andreas Scheck <andreas.scheck.home@googlemail.com>
 */

/**
 * Returns the configuration data in an associative array
 */

//---- Database settings ----//
$DB_HOST		= '127.0.0.1';
$DB_USER		= 'vplg';
$DB_PASSWORD		= 'vplg';
$DB_NAME		= 'vplg';
$DB_PORT		= '5432';
//---- End of database settings ----//

//---- Debug settings ----//
$SITE_TITLE = "PTGL";
$TITLE_SPACER = " -- ";
//---- End of debug settings ----//	

//---- Search Settings ----//
$ENABLE_COMPLEX_GRAPHS	= FALSE;
$ENABLE_MOTIF_SEARCH	= TRUE;
$ENABLE_BLAST_SEARCH	= TRUE;
$USE_LOGIC_OPERATORS	= FALSE;
//---- End of search settings ----//


//---- File settings ----//
$BUILD_FILE_PATH	= FALSE;
$IMG_ROOT_PATH		= './data/';
//---- End of file settings ----//


//---- Debug settings ----//
$DEBUG_MODE = TRUE;
//---- End of debug settings ----//

//EOF
?>