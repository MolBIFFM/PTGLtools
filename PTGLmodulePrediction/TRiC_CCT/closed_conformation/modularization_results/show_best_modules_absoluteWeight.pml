
### PML file to represent the clusters in the structure.
# Python is by default case insensitive for chain IDs but for large clusters with upper and lower case letters we have to care!
set ignore_case, off

# Load the structure
load 7ylx.cif, complex

# First color everything white
color white, all
# Select and color the found clusters
select cluster0, chain A+G+Z
select cluster1, chain B+D+E+H+Q
select cluster2, chain a+g+z
select cluster3, chain b+d+e+h+p+q

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
pseudoatom psB, chain B
pseudoatom psD, chain D
pseudoatom psE, chain E
pseudoatom psG, chain G
pseudoatom psH, chain H
pseudoatom psQ, chain Q
pseudoatom psZ, chain Z
pseudoatom psa, chain a
pseudoatom psb, chain b
pseudoatom psd, chain d
pseudoatom pse, chain e
pseudoatom psg, chain g
pseudoatom psh, chain h
pseudoatom psp, chain p
pseudoatom psq, chain q
pseudoatom psz, chain z

group pseudoatoms, ps*

# Create all distance objects originating from each chain (reverse dublicates are needed!)
# Hide the distance label, generate a pseudoatom between both pseudoatoms for which the distance was calculated
# and label it according to the edge weight
dist dist_A-B, psA////PS1, psB////PS1
hide label, dist_A-B
dist dist_A-D, psA////PS1, psD////PS1
hide label, dist_A-D
dist dist_A-G, psA////PS1, psG////PS1
hide label, dist_A-G
dist dist_A-Z, psA////PS1, psZ////PS1
hide label, dist_A-Z
dist dist_A-b, psA////PS1, psb////PS1
hide label, dist_A-b
dist dist_A-e, psA////PS1, pse////PS1
hide label, dist_A-e
dist dist_A-h, psA////PS1, psh////PS1
hide label, dist_A-h
dist dist_A-q, psA////PS1, psq////PS1
hide label, dist_A-q
dist dist_B-A, psB////PS1, psA////PS1
hide label, dist_B-A
dist dist_B-D, psB////PS1, psD////PS1
hide label, dist_B-D
dist dist_B-E, psB////PS1, psE////PS1
hide label, dist_B-E
dist dist_B-H, psB////PS1, psH////PS1
hide label, dist_B-H
dist dist_B-a, psB////PS1, psa////PS1
hide label, dist_B-a
dist dist_B-b, psB////PS1, psb////PS1
hide label, dist_B-b
dist dist_B-d, psB////PS1, psd////PS1
hide label, dist_B-d
dist dist_D-A, psD////PS1, psA////PS1
hide label, dist_D-A
dist dist_D-B, psD////PS1, psB////PS1
hide label, dist_D-B
dist dist_D-b, psD////PS1, psb////PS1
hide label, dist_D-b
dist dist_D-e, psD////PS1, pse////PS1
hide label, dist_D-e
dist dist_D-h, psD////PS1, psh////PS1
hide label, dist_D-h
dist dist_E-B, psE////PS1, psB////PS1
hide label, dist_E-B
dist dist_E-H, psE////PS1, psH////PS1
hide label, dist_E-H
dist dist_E-a, psE////PS1, psa////PS1
hide label, dist_E-a
dist dist_E-d, psE////PS1, psd////PS1
hide label, dist_E-d
dist dist_E-g, psE////PS1, psg////PS1
hide label, dist_E-g
dist dist_G-A, psG////PS1, psA////PS1
hide label, dist_G-A
dist dist_G-H, psG////PS1, psH////PS1
hide label, dist_G-H
dist dist_G-Z, psG////PS1, psZ////PS1
hide label, dist_G-Z
dist dist_G-e, psG////PS1, pse////PS1
hide label, dist_G-e
dist dist_G-q, psG////PS1, psq////PS1
hide label, dist_G-q
dist dist_H-B, psH////PS1, psB////PS1
hide label, dist_H-B
dist dist_H-E, psH////PS1, psE////PS1
hide label, dist_H-E
dist dist_H-G, psH////PS1, psG////PS1
hide label, dist_H-G
dist dist_H-z, psH////PS1, psz////PS1
hide label, dist_H-z
dist dist_H-Q, psH////PS1, psQ////PS1
hide label, dist_H-Q
dist dist_H-Z, psH////PS1, psZ////PS1
hide label, dist_H-Z
dist dist_H-a, psH////PS1, psa////PS1
hide label, dist_H-a
dist dist_H-g, psH////PS1, psg////PS1
hide label, dist_H-g
dist dist_Q-H, psQ////PS1, psH////PS1
hide label, dist_Q-H
dist dist_Q-z, psQ////PS1, psz////PS1
hide label, dist_Q-z
dist dist_Q-Z, psQ////PS1, psZ////PS1
hide label, dist_Q-Z
dist dist_Q-a, psQ////PS1, psa////PS1
hide label, dist_Q-a
dist dist_Q-g, psQ////PS1, psg////PS1
hide label, dist_Q-g
dist dist_Z-A, psZ////PS1, psA////PS1
hide label, dist_Z-A
dist dist_Z-G, psZ////PS1, psG////PS1
hide label, dist_Z-G
dist dist_Z-H, psZ////PS1, psH////PS1
hide label, dist_Z-H
dist dist_Z-Q, psZ////PS1, psQ////PS1
hide label, dist_Z-Q
dist dist_Z-z, psZ////PS1, psz////PS1
hide label, dist_Z-z
dist dist_Z-h, psZ////PS1, psh////PS1
hide label, dist_Z-h
dist dist_Z-q, psZ////PS1, psq////PS1
hide label, dist_Z-q
dist dist_a-B, psa////PS1, psB////PS1
hide label, dist_a-B
dist dist_a-E, psa////PS1, psE////PS1
hide label, dist_a-E
dist dist_a-H, psa////PS1, psH////PS1
hide label, dist_a-H
dist dist_a-Q, psa////PS1, psQ////PS1
hide label, dist_a-Q
dist dist_a-z, psa////PS1, psz////PS1
hide label, dist_a-z
dist dist_a-d, psa////PS1, psd////PS1
hide label, dist_a-d
dist dist_a-g, psa////PS1, psg////PS1
hide label, dist_a-g
dist dist_a-p, psa////PS1, psp////PS1
hide label, dist_a-p
dist dist_b-A, psb////PS1, psA////PS1
hide label, dist_b-A
dist dist_b-B, psb////PS1, psB////PS1
hide label, dist_b-B
dist dist_b-D, psb////PS1, psD////PS1
hide label, dist_b-D
dist dist_b-d, psb////PS1, psd////PS1
hide label, dist_b-d
dist dist_b-e, psb////PS1, pse////PS1
hide label, dist_b-e
dist dist_b-p, psb////PS1, psp////PS1
hide label, dist_b-p
dist dist_d-B, psd////PS1, psB////PS1
hide label, dist_d-B
dist dist_d-E, psd////PS1, psE////PS1
hide label, dist_d-E
dist dist_d-a, psd////PS1, psa////PS1
hide label, dist_d-a
dist dist_d-b, psd////PS1, psb////PS1
hide label, dist_d-b
dist dist_d-p, psd////PS1, psp////PS1
hide label, dist_d-p
dist dist_e-A, pse////PS1, psA////PS1
hide label, dist_e-A
dist dist_e-D, pse////PS1, psD////PS1
hide label, dist_e-D
dist dist_e-G, pse////PS1, psG////PS1
hide label, dist_e-G
dist dist_e-b, pse////PS1, psb////PS1
hide label, dist_e-b
dist dist_e-h, pse////PS1, psh////PS1
hide label, dist_e-h
dist dist_e-p, pse////PS1, psp////PS1
hide label, dist_e-p
dist dist_g-E, psg////PS1, psE////PS1
hide label, dist_g-E
dist dist_g-H, psg////PS1, psH////PS1
hide label, dist_g-H
dist dist_g-Q, psg////PS1, psQ////PS1
hide label, dist_g-Q
dist dist_g-a, psg////PS1, psa////PS1
hide label, dist_g-a
dist dist_g-z, psg////PS1, psz////PS1
hide label, dist_g-z
dist dist_g-p, psg////PS1, psp////PS1
hide label, dist_g-p
dist dist_h-A, psh////PS1, psA////PS1
hide label, dist_h-A
dist dist_h-D, psh////PS1, psD////PS1
hide label, dist_h-D
dist dist_h-Z, psh////PS1, psZ////PS1
hide label, dist_h-Z
dist dist_h-e, psh////PS1, pse////PS1
hide label, dist_h-e
dist dist_h-z, psh////PS1, psz////PS1
hide label, dist_h-z
dist dist_h-p, psh////PS1, psp////PS1
hide label, dist_h-p
dist dist_h-q, psh////PS1, psq////PS1
hide label, dist_h-q
dist dist_p-a, psp////PS1, psa////PS1
hide label, dist_p-a
dist dist_p-b, psp////PS1, psb////PS1
hide label, dist_p-b
dist dist_p-d, psp////PS1, psd////PS1
hide label, dist_p-d
dist dist_p-e, psp////PS1, pse////PS1
hide label, dist_p-e
dist dist_p-g, psp////PS1, psg////PS1
hide label, dist_p-g
dist dist_p-h, psp////PS1, psh////PS1
hide label, dist_p-h
dist dist_p-z, psp////PS1, psz////PS1
hide label, dist_p-z
dist dist_p-q, psp////PS1, psq////PS1
hide label, dist_p-q
dist dist_q-A, psq////PS1, psA////PS1
hide label, dist_q-A
dist dist_q-G, psq////PS1, psG////PS1
hide label, dist_q-G
dist dist_q-Z, psq////PS1, psZ////PS1
hide label, dist_q-Z
dist dist_q-h, psq////PS1, psh////PS1
hide label, dist_q-h
dist dist_q-p, psq////PS1, psp////PS1
hide label, dist_q-p
dist dist_q-z, psq////PS1, psz////PS1
hide label, dist_q-z
dist dist_z-H, psz////PS1, psH////PS1
hide label, dist_z-H
dist dist_z-Q, psz////PS1, psQ////PS1
hide label, dist_z-Q
dist dist_z-Z, psz////PS1, psZ////PS1
hide label, dist_z-Z
dist dist_z-a, psz////PS1, psa////PS1
hide label, dist_z-a
dist dist_z-g, psz////PS1, psg////PS1
hide label, dist_z-g
dist dist_z-h, psz////PS1, psh////PS1
hide label, dist_z-h
dist dist_z-p, psz////PS1, psp////PS1
hide label, dist_z-p
dist dist_z-q, psz////PS1, psq////PS1
hide label, dist_z-q

# Group all distance measurements originating from one chain to all other connected ones
group chain-A_contacts, dist_A-*
group chain-B_contacts, dist_B-*
group chain-D_contacts, dist_D-*
group chain-E_contacts, dist_E-*
group chain-G_contacts, dist_G-*
group chain-H_contacts, dist_H-*
group chain-Q_contacts, dist_Q-*
group chain-Z_contacts, dist_Z-*
group chain-a_contacts, dist_a-*
group chain-b_contacts, dist_b-*
group chain-d_contacts, dist_d-*
group chain-e_contacts, dist_e-*
group chain-g_contacts, dist_g-*
group chain-h_contacts, dist_h-*
group chain-p_contacts, dist_p-*
group chain-q_contacts, dist_q-*
group chain-z_contacts, dist_z-*

# Disallow addition/removal from the group (allow with 'open')
group chain-A_contacts, close
group chain-B_contacts, close
group chain-D_contacts, close
group chain-E_contacts, close
group chain-G_contacts, close
group chain-H_contacts, close
group chain-Q_contacts, close
group chain-Z_contacts, close
group chain-a_contacts, close
group chain-b_contacts, close
group chain-d_contacts, close
group chain-e_contacts, close
group chain-g_contacts, close
group chain-h_contacts, close
group chain-p_contacts, close
group chain-q_contacts, close
group chain-z_contacts, close

# Set pseudoatoms along the distance measurement line and label them with the corresponding edge weights
pseudoatom ps_A-B, psA or psB, label=1
pseudoatom ps_A-D, psA or psD, label=93
pseudoatom ps_A-G, psA or psG, label=112
pseudoatom ps_A-Z, psA or psZ, label=1
pseudoatom ps_A-b, psA or psb, label=1
pseudoatom ps_A-e, psA or pse, label=22
pseudoatom ps_A-h, psA or psh, label=20
pseudoatom ps_A-q, psA or psq, label=1
pseudoatom ps_B-A, psB or psA, label=1
pseudoatom ps_B-D, psB or psD, label=106
pseudoatom ps_B-E, psB or psE, label=117
pseudoatom ps_B-H, psB or psH, label=1
pseudoatom ps_B-a, psB or psa, label=1
pseudoatom ps_B-b, psB or psb, label=7
pseudoatom ps_B-d, psB or psd, label=16
pseudoatom ps_D-A, psD or psA, label=93
pseudoatom ps_D-B, psD or psB, label=106
pseudoatom ps_D-b, psD or psb, label=13
pseudoatom ps_D-e, psD or pse, label=18
pseudoatom ps_D-h, psD or psh, label=1
pseudoatom ps_E-B, psE or psB, label=117
pseudoatom ps_E-H, psE or psH, label=112
pseudoatom ps_E-a, psE or psa, label=18
pseudoatom ps_E-d, psE or psd, label=21
pseudoatom ps_E-g, psE or psg, label=4
pseudoatom ps_G-A, psG or psA, label=112
pseudoatom ps_G-H, psG or psH, label=1
pseudoatom ps_G-Z, psG or psZ, label=128
pseudoatom ps_G-e, psG or pse, label=4
pseudoatom ps_G-q, psG or psq, label=17
pseudoatom ps_H-B, psH or psB, label=1
pseudoatom ps_H-E, psH or psE, label=112
pseudoatom ps_H-G, psH or psG, label=1
pseudoatom ps_H-z, psH or psz, label=6
pseudoatom ps_H-Q, psH or psQ, label=95
pseudoatom ps_H-Z, psH or psZ, label=3
pseudoatom ps_H-a, psH or psa, label=17
pseudoatom ps_H-g, psH or psg, label=5
pseudoatom ps_Q-H, psQ or psH, label=95
pseudoatom ps_Q-z, psQ or psz, label=14
pseudoatom ps_Q-Z, psQ or psZ, label=82
pseudoatom ps_Q-a, psQ or psa, label=1
pseudoatom ps_Q-g, psQ or psg, label=16
pseudoatom ps_Z-A, psZ or psA, label=1
pseudoatom ps_Z-G, psZ or psG, label=128
pseudoatom ps_Z-H, psZ or psH, label=3
pseudoatom ps_Z-Q, psZ or psQ, label=82
pseudoatom ps_Z-z, psZ or psz, label=18
pseudoatom ps_Z-h, psZ or psh, label=6
pseudoatom ps_Z-q, psZ or psq, label=12
pseudoatom ps_a-B, psa or psB, label=1
pseudoatom ps_a-E, psa or psE, label=18
pseudoatom ps_a-H, psa or psH, label=17
pseudoatom ps_a-Q, psa or psQ, label=1
pseudoatom ps_a-z, psa or psz, label=1
pseudoatom ps_a-d, psa or psd, label=74
pseudoatom ps_a-g, psa or psg, label=115
pseudoatom ps_a-p, psa or psp, label=5
pseudoatom ps_b-A, psb or psA, label=1
pseudoatom ps_b-B, psb or psB, label=7
pseudoatom ps_b-D, psb or psD, label=13
pseudoatom ps_b-d, psb or psd, label=97
pseudoatom ps_b-e, psb or pse, label=117
pseudoatom ps_b-p, psb or psp, label=9
pseudoatom ps_d-B, psd or psB, label=16
pseudoatom ps_d-E, psd or psE, label=21
pseudoatom ps_d-a, psd or psa, label=74
pseudoatom ps_d-b, psd or psb, label=97
pseudoatom ps_d-p, psd or psp, label=6
pseudoatom ps_e-A, pse or psA, label=22
pseudoatom ps_e-D, pse or psD, label=18
pseudoatom ps_e-G, pse or psG, label=4
pseudoatom ps_e-b, pse or psb, label=117
pseudoatom ps_e-h, pse or psh, label=112
pseudoatom ps_e-p, pse or psp, label=14
pseudoatom ps_g-E, psg or psE, label=4
pseudoatom ps_g-H, psg or psH, label=5
pseudoatom ps_g-Q, psg or psQ, label=16
pseudoatom ps_g-a, psg or psa, label=115
pseudoatom ps_g-z, psg or psz, label=110
pseudoatom ps_g-p, psg or psp, label=8
pseudoatom ps_h-A, psh or psA, label=20
pseudoatom ps_h-D, psh or psD, label=1
pseudoatom ps_h-Z, psh or psZ, label=6
pseudoatom ps_h-e, psh or pse, label=112
pseudoatom ps_h-z, psh or psz, label=1
pseudoatom ps_h-p, psh or psp, label=4
pseudoatom ps_h-q, psh or psq, label=102
pseudoatom ps_p-a, psp or psa, label=5
pseudoatom ps_p-b, psp or psb, label=9
pseudoatom ps_p-d, psp or psd, label=6
pseudoatom ps_p-e, psp or pse, label=14
pseudoatom ps_p-g, psp or psg, label=8
pseudoatom ps_p-h, psp or psh, label=4
pseudoatom ps_p-z, psp or psz, label=5
pseudoatom ps_p-q, psp or psq, label=2
pseudoatom ps_q-A, psq or psA, label=1
pseudoatom ps_q-G, psq or psG, label=17
pseudoatom ps_q-Z, psq or psZ, label=12
pseudoatom ps_q-h, psq or psh, label=102
pseudoatom ps_q-p, psq or psp, label=2
pseudoatom ps_q-z, psq or psz, label=75
pseudoatom ps_z-H, psz or psH, label=6
pseudoatom ps_z-Q, psz or psQ, label=14
pseudoatom ps_z-Z, psz or psZ, label=18
pseudoatom ps_z-a, psz or psa, label=1
pseudoatom ps_z-g, psz or psg, label=110
pseudoatom ps_z-h, psz or psh, label=1
pseudoatom ps_z-p, psz or psp, label=5
pseudoatom ps_z-q, psz or psq, label=75

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

	