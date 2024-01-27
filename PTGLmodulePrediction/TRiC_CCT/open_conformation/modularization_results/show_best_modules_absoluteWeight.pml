
### PML file to represent the clusters in the structure.
# Python is by default case insensitive for chain IDs but for large clusters with upper and lower case letters we have to care!
set ignore_case, off

# Load the structure
load 7ylv.cif, complex

# First color everything white
color white, all
# Select and color the found clusters
select cluster0, chain A+G+D+Z
select cluster1, chain a+g+d+z+P
select cluster2, chain B+E+H+Q
select cluster3, chain b+e+h+q

set_color newColor0, [1.0, 0.0, 0.0] 
color newColor0, cluster0
set_color newColor1, [0.0, 1.0, 0.0] 
color newColor1, cluster1
set_color newColor2, [0.0, 0.0, 1.0] 
color newColor2, cluster2
set_color newColor3, [1.0, 1.0, 0.0] 
color newColor3, cluster3

# Create pseudoatoms at the geometric center of the chains
pseudoatom psA, chain A
pseudoatom psa, chain a
pseudoatom psB, chain B
pseudoatom psb, chain b
pseudoatom psG, chain G
pseudoatom psg, chain g
pseudoatom psD, chain D
pseudoatom psd, chain d
pseudoatom psE, chain E
pseudoatom pse, chain e
pseudoatom psZ, chain Z
pseudoatom psz, chain z
pseudoatom psH, chain H
pseudoatom psh, chain h
pseudoatom psQ, chain Q
pseudoatom psq, chain q
pseudoatom psP, chain P

group pseudoatoms, ps*

# Create all distance objects originating from each chain (reverse dublicates are needed!)
# Hide the distance label, generate a pseudoatom between both pseudoatoms for which the distance was calculated
# and label it according to the edge weight
dist dist_A-G, psA////PS1, psG////PS1
hide label, dist_A-G
dist dist_A-D, psA////PS1, psD////PS1
hide label, dist_A-D
dist dist_A-h, psA////PS1, psh////PS1
hide label, dist_A-h
dist dist_A-q, psA////PS1, psq////PS1
hide label, dist_A-q
dist dist_a-P, psa////PS1, psP////PS1
hide label, dist_a-P
dist dist_a-g, psa////PS1, psg////PS1
hide label, dist_a-g
dist dist_a-d, psa////PS1, psd////PS1
hide label, dist_a-d
dist dist_a-H, psa////PS1, psH////PS1
hide label, dist_a-H
dist dist_a-Q, psa////PS1, psQ////PS1
hide label, dist_a-Q
dist dist_B-b, psB////PS1, psb////PS1
hide label, dist_B-b
dist dist_B-D, psB////PS1, psD////PS1
hide label, dist_B-D
dist dist_B-d, psB////PS1, psd////PS1
hide label, dist_B-d
dist dist_B-E, psB////PS1, psE////PS1
hide label, dist_B-E
dist dist_b-B, psb////PS1, psB////PS1
hide label, dist_b-B
dist dist_b-D, psb////PS1, psD////PS1
hide label, dist_b-D
dist dist_b-d, psb////PS1, psd////PS1
hide label, dist_b-d
dist dist_b-e, psb////PS1, pse////PS1
hide label, dist_b-e
dist dist_G-A, psG////PS1, psA////PS1
hide label, dist_G-A
dist dist_G-Z, psG////PS1, psZ////PS1
hide label, dist_G-Z
dist dist_G-z, psG////PS1, psz////PS1
hide label, dist_G-z
dist dist_G-q, psG////PS1, psq////PS1
hide label, dist_G-q
dist dist_g-a, psg////PS1, psa////PS1
hide label, dist_g-a
dist dist_g-Z, psg////PS1, psZ////PS1
hide label, dist_g-Z
dist dist_g-z, psg////PS1, psz////PS1
hide label, dist_g-z
dist dist_g-Q, psg////PS1, psQ////PS1
hide label, dist_g-Q
dist dist_D-A, psD////PS1, psA////PS1
hide label, dist_D-A
dist dist_D-B, psD////PS1, psB////PS1
hide label, dist_D-B
dist dist_D-b, psD////PS1, psb////PS1
hide label, dist_D-b
dist dist_D-e, psD////PS1, pse////PS1
hide label, dist_D-e
dist dist_d-a, psd////PS1, psa////PS1
hide label, dist_d-a
dist dist_d-B, psd////PS1, psB////PS1
hide label, dist_d-B
dist dist_d-b, psd////PS1, psb////PS1
hide label, dist_d-b
dist dist_d-E, psd////PS1, psE////PS1
hide label, dist_d-E
dist dist_E-B, psE////PS1, psB////PS1
hide label, dist_E-B
dist dist_E-d, psE////PS1, psd////PS1
hide label, dist_E-d
dist dist_E-H, psE////PS1, psH////PS1
hide label, dist_E-H
dist dist_e-b, pse////PS1, psb////PS1
hide label, dist_e-b
dist dist_e-D, pse////PS1, psD////PS1
hide label, dist_e-D
dist dist_e-h, pse////PS1, psh////PS1
hide label, dist_e-h
dist dist_Z-G, psZ////PS1, psG////PS1
hide label, dist_Z-G
dist dist_Z-g, psZ////PS1, psg////PS1
hide label, dist_Z-g
dist dist_Z-z, psZ////PS1, psz////PS1
hide label, dist_Z-z
dist dist_Z-Q, psZ////PS1, psQ////PS1
hide label, dist_Z-Q
dist dist_z-G, psz////PS1, psG////PS1
hide label, dist_z-G
dist dist_z-g, psz////PS1, psg////PS1
hide label, dist_z-g
dist dist_z-Z, psz////PS1, psZ////PS1
hide label, dist_z-Z
dist dist_z-q, psz////PS1, psq////PS1
hide label, dist_z-q
dist dist_H-a, psH////PS1, psa////PS1
hide label, dist_H-a
dist dist_H-E, psH////PS1, psE////PS1
hide label, dist_H-E
dist dist_H-Q, psH////PS1, psQ////PS1
hide label, dist_H-Q
dist dist_h-A, psh////PS1, psA////PS1
hide label, dist_h-A
dist dist_h-e, psh////PS1, pse////PS1
hide label, dist_h-e
dist dist_h-q, psh////PS1, psq////PS1
hide label, dist_h-q
dist dist_Q-a, psQ////PS1, psa////PS1
hide label, dist_Q-a
dist dist_Q-g, psQ////PS1, psg////PS1
hide label, dist_Q-g
dist dist_Q-Z, psQ////PS1, psZ////PS1
hide label, dist_Q-Z
dist dist_Q-H, psQ////PS1, psH////PS1
hide label, dist_Q-H
dist dist_Q-q, psQ////PS1, psq////PS1
hide label, dist_Q-q
dist dist_q-A, psq////PS1, psA////PS1
hide label, dist_q-A
dist dist_q-G, psq////PS1, psG////PS1
hide label, dist_q-G
dist dist_q-z, psq////PS1, psz////PS1
hide label, dist_q-z
dist dist_q-h, psq////PS1, psh////PS1
hide label, dist_q-h
dist dist_q-Q, psq////PS1, psQ////PS1
hide label, dist_q-Q
dist dist_P-a, psP////PS1, psa////PS1
hide label, dist_P-a

# Group all distance measurements originating from one chain to all other connected ones
group chain-A_contacts, dist_A-*
group chain-a_contacts, dist_a-*
group chain-B_contacts, dist_B-*
group chain-b_contacts, dist_b-*
group chain-G_contacts, dist_G-*
group chain-g_contacts, dist_g-*
group chain-D_contacts, dist_D-*
group chain-d_contacts, dist_d-*
group chain-E_contacts, dist_E-*
group chain-e_contacts, dist_e-*
group chain-Z_contacts, dist_Z-*
group chain-z_contacts, dist_z-*
group chain-H_contacts, dist_H-*
group chain-h_contacts, dist_h-*
group chain-Q_contacts, dist_Q-*
group chain-q_contacts, dist_q-*
group chain-P_contacts, dist_P-*

# Disallow addition/removal from the group (allow with 'open')
group chain-A_contacts, close
group chain-a_contacts, close
group chain-B_contacts, close
group chain-b_contacts, close
group chain-G_contacts, close
group chain-g_contacts, close
group chain-D_contacts, close
group chain-d_contacts, close
group chain-E_contacts, close
group chain-e_contacts, close
group chain-Z_contacts, close
group chain-z_contacts, close
group chain-H_contacts, close
group chain-h_contacts, close
group chain-Q_contacts, close
group chain-q_contacts, close
group chain-P_contacts, close

# Set pseudoatoms along the distance measurement line and label them with the corresponding edge weights
pseudoatom ps_A-G, psA or psG, label=40
pseudoatom ps_A-D, psA or psD, label=37
pseudoatom ps_A-h, psA or psh, label=16
pseudoatom ps_A-q, psA or psq, label=2
pseudoatom ps_a-P, psa or psP, label=2
pseudoatom ps_a-g, psa or psg, label=51
pseudoatom ps_a-d, psa or psd, label=44
pseudoatom ps_a-H, psa or psH, label=18
pseudoatom ps_a-Q, psa or psQ, label=1
pseudoatom ps_B-b, psB or psb, label=22
pseudoatom ps_B-D, psB or psD, label=42
pseudoatom ps_B-d, psB or psd, label=14
pseudoatom ps_B-E, psB or psE, label=79
pseudoatom ps_b-B, psb or psB, label=22
pseudoatom ps_b-D, psb or psD, label=12
pseudoatom ps_b-d, psb or psd, label=42
pseudoatom ps_b-e, psb or pse, label=61
pseudoatom ps_G-A, psG or psA, label=40
pseudoatom ps_G-Z, psG or psZ, label=48
pseudoatom ps_G-z, psG or psz, label=2
pseudoatom ps_G-q, psG or psq, label=14
pseudoatom ps_g-a, psg or psa, label=51
pseudoatom ps_g-Z, psg or psZ, label=3
pseudoatom ps_g-z, psg or psz, label=46
pseudoatom ps_g-Q, psg or psQ, label=20
pseudoatom ps_D-A, psD or psA, label=37
pseudoatom ps_D-B, psD or psB, label=42
pseudoatom ps_D-b, psD or psb, label=12
pseudoatom ps_D-e, psD or pse, label=25
pseudoatom ps_d-a, psd or psa, label=44
pseudoatom ps_d-B, psd or psB, label=14
pseudoatom ps_d-b, psd or psb, label=42
pseudoatom ps_d-E, psd or psE, label=22
pseudoatom ps_E-B, psE or psB, label=79
pseudoatom ps_E-d, psE or psd, label=22
pseudoatom ps_E-H, psE or psH, label=76
pseudoatom ps_e-b, pse or psb, label=61
pseudoatom ps_e-D, pse or psD, label=25
pseudoatom ps_e-h, pse or psh, label=81
pseudoatom ps_Z-G, psZ or psG, label=48
pseudoatom ps_Z-g, psZ or psg, label=3
pseudoatom ps_Z-z, psZ or psz, label=13
pseudoatom ps_Z-Q, psZ or psQ, label=45
pseudoatom ps_z-G, psz or psG, label=2
pseudoatom ps_z-g, psz or psg, label=46
pseudoatom ps_z-Z, psz or psZ, label=13
pseudoatom ps_z-q, psz or psq, label=54
pseudoatom ps_H-a, psH or psa, label=18
pseudoatom ps_H-E, psH or psE, label=76
pseudoatom ps_H-Q, psH or psQ, label=68
pseudoatom ps_h-A, psh or psA, label=16
pseudoatom ps_h-e, psh or pse, label=81
pseudoatom ps_h-q, psh or psq, label=66
pseudoatom ps_Q-a, psQ or psa, label=1
pseudoatom ps_Q-g, psQ or psg, label=20
pseudoatom ps_Q-Z, psQ or psZ, label=45
pseudoatom ps_Q-H, psQ or psH, label=68
pseudoatom ps_Q-q, psQ or psq, label=4
pseudoatom ps_q-A, psq or psA, label=2
pseudoatom ps_q-G, psq or psG, label=14
pseudoatom ps_q-z, psq or psz, label=54
pseudoatom ps_q-h, psq or psh, label=66
pseudoatom ps_q-Q, psq or psQ, label=4
pseudoatom ps_P-a, psP or psa, label=2

group label_pseudoatoms, ps_*

# Set some parameters
set dash_gap, 0.7
set dash_radius, 0.2
set dash_color, black
set label_position,(1,1,1) # label offset
set label_color, black
bg_color white
set cartoon_highlight_color, grey40

## If coordinates of the pseudoatoms should be written
#xyz = cmd.get_coords("pseudoatoms", 1)
#python
#import numpy as np
#with open('pseudoatomCoords.npy', 'wb') as f:
#    np.save(f, xyz)
#python end

	