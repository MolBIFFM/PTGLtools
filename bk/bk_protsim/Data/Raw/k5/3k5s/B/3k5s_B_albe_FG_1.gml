graph [
  id 1
  label "VPLG Protein Graph 3k5s-B-albe-FG1[7V,7E]"
  comment "pdbid=3k5s|graphtype=albe|chainid=B|graphclass=folding graph|foldinggraphnumber=2|"
  directed 0
  isplanar 0
  creator "PLCC version 0.98.1"
  pdb_id "3k5s"
  chain_id "B"
  graph_type "albe"
  is_protein_graph 0
  is_folding_graph 1
  is_SSE_graph 1
  is_AA_graph 0
  is_all_chains_graph 0
  node [
    id 0
    label "0-E"
    num_in_chain 8
    num_residues 6
    pdb_res_start "B-112- "
    pdb_res_end "B-117- "
    dssp_res_start 330
    dssp_res_end 335
    pdb_residues_full "B-112- ,B-113- ,B-114- ,B-115- ,B-116- ,B-117- "
    aa_sequence "YVGHVM"
    index_in_parent_pg "7"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 1
    label "1-E"
    num_in_chain 9
    num_residues 4
    pdb_res_start "B-125- "
    pdb_res_end "B-128- "
    dssp_res_start 343
    dssp_res_end 346
    pdb_residues_full "B-125- ,B-126- ,B-127- ,B-128- "
    aa_sequence "TVMR"
    index_in_parent_pg "8"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 2
    label "2-E"
    num_in_chain 10
    num_residues 8
    pdb_res_start "B-145- "
    pdb_res_end "B-152- "
    dssp_res_start 363
    dssp_res_end 370
    pdb_residues_full "B-145- ,B-146- ,B-147- ,B-148- ,B-149- ,B-150- ,B-151- ,B-152- "
    aa_sequence "RYNILKQT"
    index_in_parent_pg "9"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 3
    label "3-E"
    num_in_chain 11
    num_residues 3
    pdb_res_start "B-161- "
    pdb_res_end "B-163- "
    dssp_res_start 379
    dssp_res_end 381
    pdb_residues_full "B-161- ,B-162- ,B-163- "
    aa_sequence "FYI"
    index_in_parent_pg "10"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 4
    label "4-E"
    num_in_chain 12
    num_residues 4
    pdb_res_start "B-169- "
    pdb_res_end "B-172- "
    dssp_res_start 387
    dssp_res_end 390
    pdb_residues_full "B-169- ,B-170- ,B-171- ,B-172- "
    aa_sequence "DIVT"
    index_in_parent_pg "11"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 5
    label "5-E"
    num_in_chain 13
    num_residues 9
    pdb_res_start "B-188- "
    pdb_res_end "B-196- "
    dssp_res_start 406
    dssp_res_end 414
    pdb_residues_full "B-188- ,B-189- ,B-190- ,B-191- ,B-192- ,B-193- ,B-194- ,B-195- ,B-196- "
    aa_sequence "KYELVIEAK"
    index_in_parent_pg "12"
    sse_type "E"
    fg_notation_label "e"
  ]
  node [
    id 6
    label "6-E"
    num_in_chain 14
    num_residues 11
    pdb_res_start "B-206- "
    pdb_res_end "B-216- "
    dssp_res_start 424
    dssp_res_end 434
    pdb_residues_full "B-206- ,B-207- ,B-208- ,B-209- ,B-210- ,B-211- ,B-212- ,B-213- ,B-214- ,B-215- ,B-216- "
    aa_sequence "TGTATATILID"
    index_in_parent_pg "13"
    sse_type "E"
    fg_notation_label "e"
  ]
  edge [
    source 0
    target 1
    label "a"
    spatial "a"
  ]
  edge [
    source 0
    target 6
    label "p"
    spatial "p"
  ]
  edge [
    source 1
    target 4
    label "a"
    spatial "a"
  ]
  edge [
    source 2
    target 3
    label "a"
    spatial "a"
  ]
  edge [
    source 2
    target 5
    label "a"
    spatial "a"
  ]
  edge [
    source 3
    target 4
    label "a"
    spatial "a"
  ]
  edge [
    source 5
    target 6
    label "a"
    spatial "a"
  ]
]
