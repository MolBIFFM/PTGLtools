# PTGLtools (formerly labeled VPLG)

This package contains software for the analysis of topologies of protein structures.
The main part of this package is the computation of topology graphs based on three-dimensional structures from the Protein Data Bank (PDB) and the secondary structure assignment from Define Secondary Structure of Proteins (DSSP).
This software is able to read in legacy PDB files and macromolecular Crystallographic Information Files (mmCIFs) and therefore able to process large structures of > 62 chains or > 99,999 atoms.
The software is used to fill the Protein Topology Graph Library (PTGL).

## Contained software

Each top-level directory contains one software part of this package.

- 'PTGLgraphComputation' (formerly labeled PLCC) is the main part of PTGLtools computing graphs representing the topology of three-dimensional protein structures
- 'PTGLdynamics' [unpublished & WIP]: Pipeline to analyze snapshots of MD simulation with PTGLgraphComputation
- 'splitpdb': splits a multi-model legacy PDB file into single PDB files of one model
- 'graphletAnalyzer': finds graphlets in graph modeling language (GML) files
- 'bk': implementation of the Bron_Kerbosch algorithm to find cliques
- 'vpg' [closed]: graphical user interface for PTGLgraphComputation

This is just the README for the code repository, see the doc/ subdirectory for help with the software in this repo.
For example, find the documentation of PTGLgraphComputation in [/PTGLgraphComputation/doc](https://github.com/MolBIFFM/PTGLtools/tree/main/PTGLgraphComputation/doc).

## Installing

Both the 'splitpdb' and the 'PTGLgraphComputation' directories contain not only the source code but the full Netbeans projects.
You can directly open the project using the Netbeans IDE.
Simply build the project: button in the Netbeans taskbar or [F11] by default.
The ready to use JAR file and lib folder can be found in the newly created dist folder by default.
You can test the program and access the help by running `java -jar PTGLgraphComputation.jar --help`

## Websites

- [PTGL website](http://ptgl.uni-frankfurt.de/): the website PTGLgraphComputation fills with data
- [PTGLgraphComputation website](http://www.bioinformatik.uni-frankfurt.de/tools/PTGLgraphComputation/index.html)
- [Molbi group](https://www.uni-frankfurt.de/57211826/People)

## Cite us

If you use our software, please [cite us](http://ptgl.uni-frankfurt.de/citing.php):

Wolf JN, Keßler M, Ackermann J, Koch I (2020). PTGL: extension to graph-based topologies of cryo-EM data for large protein structures. *Bioinformatics*. ([DOI 10.1093/bta706 at Oxford University Press](https://academic.oup.com/bioinformatics/article/37/7/1032/5891172?guestAccessKey=c5d067b4-1e81-4893-81a4-55076e9c0e17))

## License

This project is licensed under the GPL 2 license - see the [license](https://github.com/MolBIFFM/PTGLtools/blob/master/LICENSE) file for details.

## Notes

The source is at GitHub, at https://github.com/dfsp-spirit/vplg/ .
The code has been forked at April 26, 2019, to and will be maintained on https://github.com/MolBIFFM/PTGLtools/ .

PTGLtools is free software and comes WITHOUT ANY WARRANTY. See the LICENSE file for details.

--
jnw
