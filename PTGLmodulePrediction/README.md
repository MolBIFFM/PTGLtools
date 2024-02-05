## Graph-theoretical prediction of biological modules

This folder is part of the Supporting Information for the paper '_Graph-theoretical prediction of biological modules in quaternary structures of large protein complexes_' by Gisdon et al.

### Approach
We modelled the quaternary structure of protein complexes as undirected, labelled graphs called complex graphs (Wolf et al. (2020), _Bioinformatics_, 37:1032-1034). The vertices represent protein chains and the edges spatial chain-chain contacts. We applied the Leiden clustering algorithm (Traag et al. (2019), _Scientific Reports_, 9:5233) to compute clusters within the protein complex graph. We hypothesized that the computed clusters correspond to functional biological modules. We discussed three case studies with protein complexes of different size. The respiratory complex I has 45 protein chains (44 unique protein chains), the head of the bacteriophage Phi29 has 400 protein chains (2 unique protein chains), and the group II chaperonin TRiC/CCT has 17 unique protein chains with 16 protein chains being paralogs. The chaperonin TRiC/CCT was analyzed in an open and a closed conformation.

### Folders for the discussed case studies
Jupyter notebooks with the workflow are available in the folder of the presented case studies
- ComplexI
- Phi29 (comment: please unzip the structure files 6qyd.cif.zip and annotated_6qyd.cif.zip in case you want to use the provided structures)
- TRiC_CCT

### Provided Python code
Custom Python functions:    
- generalFunctions.py                 
- pymolRepresentationClustering.py   

### Additional files
The settings for Protein Topology Graph Computation (PTGL) software is provided in the file 'PTGLgraphComputation_settings.txt'. 
The relevant files for each case study are provided to run the examples without the software PTGL, DSSP, or PyMOL. Output files are provided.
The structure file and the DSSP-annotated structure file for the case study Phi29 are zipped to mach the allowed file size on GitHub. Please unzip before running the notebook.
