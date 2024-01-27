## Graph-theoretical prediction of biological modules

This folder is part of the Supporting Information for the paper '_Graph-theoretical prediction of biological modules in quaternary structures of large protein complexes_' by Gisdon et al.

### Approach
We modelled the quaternary structure of protein complexes as undirected, labelled graphs called complex graphs (Wolf et al. (2020), Bioinformatics, 37:1032-1034). The vertices represent protein chains and the edges spatial chain-chain contacts. We applied the Leiden clustering algorithm to compute clusters within the protein complex graph. We hypothesized that the computed clusters correspond to functional biological
modules.

### Folders for the discussed case studies
Jupyter notebooks with the workflow are available in the folder of the presented case studies
- ComplexI
- Phi29 (comment: please unzip the structure files 6qyd.cif.zip and annotated_6qyd.cif.zip before running the notebook)
- TRiC_CCT

### Provided Python code
Custom Python functions are available:    
- generalFunctions.py                 
- pymolRepresentationClustering.py   

### Additional files
The settings for Protein Topology Graph Computation (PTGL) software is provided in the file 'PTGLgraphComputation_settings.txt'. 
The relevant files for each case study are provided to run the examples without the software PTGL, DSSP, or PyMOL. Output files are provided.
The structure file and the DSSP-annotated structure file for the case study Phi29 are zipped to mach the allowed file size on GitHub. Please unzip before running the notebook.
