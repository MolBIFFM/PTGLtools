
### PML file to represent the clusters in the structure.
# Python is by default case insensitive for chain IDs but for large clusters with upper and lower case letters we have to care!
set ignore_case, off

# Load the structure
load 5xtd.cif, complex

# First color everything white
color white, all
# Select and color the found clusters
select cluster0, chain A+F+K+L+M+O
select cluster1, chain B+C+E+G+H+I+J+N+P+Q+T
select cluster2, chain S+U+W+h+j+k+m+s+u
select cluster3, chain X+Y+Z+a+b+c+d+e+l+n+o+p+r+v
select cluster4, chain V+f+g+i+w

set_color newColor0, [1.0, 0.0, 0.0] 
color newColor0, cluster0
set_color newColor1, [0.0, 1.0, 0.0] 
color newColor1, cluster1
set_color newColor2, [0.0, 0.0, 1.0] 
color newColor2, cluster2
set_color newColor3, [1.0, 1.0, 0.0] 
color newColor3, cluster3
set_color newColor4, [1.0, 0.0, 1.0] 
color newColor4, cluster4

# Create pseudoatoms at the geometric center of the chains
pseudoatom psA, chain A
pseudoatom psB, chain B
pseudoatom psC, chain C
pseudoatom psE, chain E
pseudoatom psF, chain F
pseudoatom psG, chain G
pseudoatom psH, chain H
pseudoatom psI, chain I
pseudoatom psJ, chain J
pseudoatom psK, chain K
pseudoatom psL, chain L
pseudoatom psM, chain M
pseudoatom psN, chain N
pseudoatom psO, chain O
pseudoatom psP, chain P
pseudoatom psQ, chain Q
pseudoatom psS, chain S
pseudoatom psT, chain T
pseudoatom psU, chain U
pseudoatom psV, chain V
pseudoatom psW, chain W
pseudoatom psX, chain X
pseudoatom psY, chain Y
pseudoatom psZ, chain Z
pseudoatom psa, chain a
pseudoatom psb, chain b
pseudoatom psc, chain c
pseudoatom psd, chain d
pseudoatom pse, chain e
pseudoatom psf, chain f
pseudoatom psg, chain g
pseudoatom psh, chain h
pseudoatom psi, chain i
pseudoatom psj, chain j
pseudoatom psk, chain k
pseudoatom psl, chain l
pseudoatom psm, chain m
pseudoatom psn, chain n
pseudoatom pso, chain o
pseudoatom psp, chain p
pseudoatom psr, chain r
pseudoatom pss, chain s
pseudoatom psu, chain u
pseudoatom psv, chain v
pseudoatom psw, chain w

group pseudoatoms, ps*

# Create all distance objects originating from each chain (reverse dublicates are needed!)
# Hide the distance label, generate a pseudoatom between both pseudoatoms for which the distance was calculated
# and label it according to the edge weight
dist dist_A-I, psA////PS1, psI////PS1
hide label, dist_A-I
dist dist_A-K, psA////PS1, psK////PS1
hide label, dist_A-K
dist dist_A-L, psA////PS1, psL////PS1
hide label, dist_A-L
dist dist_A-M, psA////PS1, psM////PS1
hide label, dist_A-M
dist dist_A-O, psA////PS1, psO////PS1
hide label, dist_A-O
dist dist_B-C, psB////PS1, psC////PS1
hide label, dist_B-C
dist dist_B-H, psB////PS1, psH////PS1
hide label, dist_B-H
dist dist_B-I, psB////PS1, psI////PS1
hide label, dist_B-I
dist dist_B-J, psB////PS1, psJ////PS1
hide label, dist_B-J
dist dist_B-s, psB////PS1, pss////PS1
hide label, dist_B-s
dist dist_B-L, psB////PS1, psL////PS1
hide label, dist_B-L
dist dist_B-M, psB////PS1, psM////PS1
hide label, dist_B-M
dist dist_B-N, psB////PS1, psN////PS1
hide label, dist_B-N
dist dist_B-P, psB////PS1, psP////PS1
hide label, dist_B-P
dist dist_B-Q, psB////PS1, psQ////PS1
hide label, dist_B-Q
dist dist_B-S, psB////PS1, psS////PS1
hide label, dist_B-S
dist dist_B-T, psB////PS1, psT////PS1
hide label, dist_B-T
dist dist_B-U, psB////PS1, psU////PS1
hide label, dist_B-U
dist dist_B-W, psB////PS1, psW////PS1
hide label, dist_B-W
dist dist_C-B, psC////PS1, psB////PS1
hide label, dist_C-B
dist dist_C-j, psC////PS1, psj////PS1
hide label, dist_C-j
dist dist_C-E, psC////PS1, psE////PS1
hide label, dist_C-E
dist dist_C-J, psC////PS1, psJ////PS1
hide label, dist_C-J
dist dist_C-s, psC////PS1, pss////PS1
hide label, dist_C-s
dist dist_C-N, psC////PS1, psN////PS1
hide label, dist_C-N
dist dist_C-P, psC////PS1, psP////PS1
hide label, dist_C-P
dist dist_C-Q, psC////PS1, psQ////PS1
hide label, dist_C-Q
dist dist_E-C, psE////PS1, psC////PS1
hide label, dist_E-C
dist dist_E-j, psE////PS1, psj////PS1
hide label, dist_E-j
dist dist_E-G, psE////PS1, psG////PS1
hide label, dist_E-G
dist dist_E-J, psE////PS1, psJ////PS1
hide label, dist_E-J
dist dist_E-L, psE////PS1, psL////PS1
hide label, dist_E-L
dist dist_E-M, psE////PS1, psM////PS1
hide label, dist_E-M
dist dist_E-P, psE////PS1, psP////PS1
hide label, dist_E-P
dist dist_E-Q, psE////PS1, psQ////PS1
hide label, dist_E-Q
dist dist_F-M, psF////PS1, psM////PS1
hide label, dist_F-M
dist dist_G-E, psG////PS1, psE////PS1
hide label, dist_G-E
dist dist_H-B, psH////PS1, psB////PS1
hide label, dist_H-B
dist dist_H-I, psH////PS1, psI////PS1
hide label, dist_H-I
dist dist_H-w, psH////PS1, psw////PS1
hide label, dist_H-w
dist dist_H-P, psH////PS1, psP////PS1
hide label, dist_H-P
dist dist_H-Q, psH////PS1, psQ////PS1
hide label, dist_H-Q
dist dist_I-A, psI////PS1, psA////PS1
hide label, dist_I-A
dist dist_I-B, psI////PS1, psB////PS1
hide label, dist_I-B
dist dist_I-H, psI////PS1, psH////PS1
hide label, dist_I-H
dist dist_I-W, psI////PS1, psW////PS1
hide label, dist_I-W
dist dist_I-M, psI////PS1, psM////PS1
hide label, dist_I-M
dist dist_I-N, psI////PS1, psN////PS1
hide label, dist_I-N
dist dist_I-P, psI////PS1, psP////PS1
hide label, dist_I-P
dist dist_I-Q, psI////PS1, psQ////PS1
hide label, dist_I-Q
dist dist_J-B, psJ////PS1, psB////PS1
hide label, dist_J-B
dist dist_J-C, psJ////PS1, psC////PS1
hide label, dist_J-C
dist dist_J-E, psJ////PS1, psE////PS1
hide label, dist_J-E
dist dist_J-T, psJ////PS1, psT////PS1
hide label, dist_J-T
dist dist_J-j, psJ////PS1, psj////PS1
hide label, dist_J-j
dist dist_J-m, psJ////PS1, psm////PS1
hide label, dist_J-m
dist dist_J-s, psJ////PS1, pss////PS1
hide label, dist_J-s
dist dist_J-L, psJ////PS1, psL////PS1
hide label, dist_J-L
dist dist_J-M, psJ////PS1, psM////PS1
hide label, dist_J-M
dist dist_J-P, psJ////PS1, psP////PS1
hide label, dist_J-P
dist dist_K-A, psK////PS1, psA////PS1
hide label, dist_K-A
dist dist_K-L, psK////PS1, psL////PS1
hide label, dist_K-L
dist dist_K-M, psK////PS1, psM////PS1
hide label, dist_K-M
dist dist_K-O, psK////PS1, psO////PS1
hide label, dist_K-O
dist dist_L-A, psL////PS1, psA////PS1
hide label, dist_L-A
dist dist_L-B, psL////PS1, psB////PS1
hide label, dist_L-B
dist dist_L-E, psL////PS1, psE////PS1
hide label, dist_L-E
dist dist_L-J, psL////PS1, psJ////PS1
hide label, dist_L-J
dist dist_L-K, psL////PS1, psK////PS1
hide label, dist_L-K
dist dist_L-T, psL////PS1, psT////PS1
hide label, dist_L-T
dist dist_L-M, psL////PS1, psM////PS1
hide label, dist_L-M
dist dist_L-N, psL////PS1, psN////PS1
hide label, dist_L-N
dist dist_L-O, psL////PS1, psO////PS1
hide label, dist_L-O
dist dist_L-P, psL////PS1, psP////PS1
hide label, dist_L-P
dist dist_M-A, psM////PS1, psA////PS1
hide label, dist_M-A
dist dist_M-B, psM////PS1, psB////PS1
hide label, dist_M-B
dist dist_M-E, psM////PS1, psE////PS1
hide label, dist_M-E
dist dist_M-F, psM////PS1, psF////PS1
hide label, dist_M-F
dist dist_M-I, psM////PS1, psI////PS1
hide label, dist_M-I
dist dist_M-J, psM////PS1, psJ////PS1
hide label, dist_M-J
dist dist_M-K, psM////PS1, psK////PS1
hide label, dist_M-K
dist dist_M-L, psM////PS1, psL////PS1
hide label, dist_M-L
dist dist_M-N, psM////PS1, psN////PS1
hide label, dist_M-N
dist dist_M-O, psM////PS1, psO////PS1
hide label, dist_M-O
dist dist_M-P, psM////PS1, psP////PS1
hide label, dist_M-P
dist dist_M-Q, psM////PS1, psQ////PS1
hide label, dist_M-Q
dist dist_M-T, psM////PS1, psT////PS1
hide label, dist_M-T
dist dist_N-B, psN////PS1, psB////PS1
hide label, dist_N-B
dist dist_N-C, psN////PS1, psC////PS1
hide label, dist_N-C
dist dist_N-I, psN////PS1, psI////PS1
hide label, dist_N-I
dist dist_N-L, psN////PS1, psL////PS1
hide label, dist_N-L
dist dist_N-M, psN////PS1, psM////PS1
hide label, dist_N-M
dist dist_N-T, psN////PS1, psT////PS1
hide label, dist_N-T
dist dist_N-s, psN////PS1, pss////PS1
hide label, dist_N-s
dist dist_N-Q, psN////PS1, psQ////PS1
hide label, dist_N-Q
dist dist_O-A, psO////PS1, psA////PS1
hide label, dist_O-A
dist dist_O-K, psO////PS1, psK////PS1
hide label, dist_O-K
dist dist_O-L, psO////PS1, psL////PS1
hide label, dist_O-L
dist dist_O-M, psO////PS1, psM////PS1
hide label, dist_O-M
dist dist_O-T, psO////PS1, psT////PS1
hide label, dist_O-T
dist dist_P-B, psP////PS1, psB////PS1
hide label, dist_P-B
dist dist_P-C, psP////PS1, psC////PS1
hide label, dist_P-C
dist dist_P-E, psP////PS1, psE////PS1
hide label, dist_P-E
dist dist_P-H, psP////PS1, psH////PS1
hide label, dist_P-H
dist dist_P-I, psP////PS1, psI////PS1
hide label, dist_P-I
dist dist_P-J, psP////PS1, psJ////PS1
hide label, dist_P-J
dist dist_P-L, psP////PS1, psL////PS1
hide label, dist_P-L
dist dist_P-M, psP////PS1, psM////PS1
hide label, dist_P-M
dist dist_P-Q, psP////PS1, psQ////PS1
hide label, dist_P-Q
dist dist_Q-B, psQ////PS1, psB////PS1
hide label, dist_Q-B
dist dist_Q-C, psQ////PS1, psC////PS1
hide label, dist_Q-C
dist dist_Q-E, psQ////PS1, psE////PS1
hide label, dist_Q-E
dist dist_Q-H, psQ////PS1, psH////PS1
hide label, dist_Q-H
dist dist_Q-I, psQ////PS1, psI////PS1
hide label, dist_Q-I
dist dist_Q-M, psQ////PS1, psM////PS1
hide label, dist_Q-M
dist dist_Q-N, psQ////PS1, psN////PS1
hide label, dist_Q-N
dist dist_Q-P, psQ////PS1, psP////PS1
hide label, dist_Q-P
dist dist_Q-i, psQ////PS1, psi////PS1
hide label, dist_Q-i
dist dist_Q-j, psQ////PS1, psj////PS1
hide label, dist_Q-j
dist dist_Q-k, psQ////PS1, psk////PS1
hide label, dist_Q-k
dist dist_Q-l, psQ////PS1, psl////PS1
hide label, dist_Q-l
dist dist_Q-r, psQ////PS1, psr////PS1
hide label, dist_Q-r
dist dist_Q-s, psQ////PS1, pss////PS1
hide label, dist_Q-s
dist dist_Q-w, psQ////PS1, psw////PS1
hide label, dist_Q-w
dist dist_Q-T, psQ////PS1, psT////PS1
hide label, dist_Q-T
dist dist_Q-W, psQ////PS1, psW////PS1
hide label, dist_Q-W
dist dist_Q-e, psQ////PS1, pse////PS1
hide label, dist_Q-e
dist dist_S-B, psS////PS1, psB////PS1
hide label, dist_S-B
dist dist_S-W, psS////PS1, psW////PS1
hide label, dist_S-W
dist dist_S-m, psS////PS1, psm////PS1
hide label, dist_S-m
dist dist_S-s, psS////PS1, pss////PS1
hide label, dist_S-s
dist dist_S-u, psS////PS1, psu////PS1
hide label, dist_S-u
dist dist_S-h, psS////PS1, psh////PS1
hide label, dist_S-h
dist dist_T-B, psT////PS1, psB////PS1
hide label, dist_T-B
dist dist_T-J, psT////PS1, psJ////PS1
hide label, dist_T-J
dist dist_T-L, psT////PS1, psL////PS1
hide label, dist_T-L
dist dist_T-M, psT////PS1, psM////PS1
hide label, dist_T-M
dist dist_T-N, psT////PS1, psN////PS1
hide label, dist_T-N
dist dist_T-O, psT////PS1, psO////PS1
hide label, dist_T-O
dist dist_T-Q, psT////PS1, psQ////PS1
hide label, dist_T-Q
dist dist_U-B, psU////PS1, psB////PS1
hide label, dist_U-B
dist dist_U-j, psU////PS1, psj////PS1
hide label, dist_U-j
dist dist_U-W, psU////PS1, psW////PS1
hide label, dist_U-W
dist dist_U-m, psU////PS1, psm////PS1
hide label, dist_U-m
dist dist_U-s, psU////PS1, pss////PS1
hide label, dist_U-s
dist dist_U-u, psU////PS1, psu////PS1
hide label, dist_U-u
dist dist_U-h, psU////PS1, psh////PS1
hide label, dist_U-h
dist dist_V-i, psV////PS1, psi////PS1
hide label, dist_V-i
dist dist_V-k, psV////PS1, psk////PS1
hide label, dist_V-k
dist dist_V-l, psV////PS1, psl////PS1
hide label, dist_V-l
dist dist_V-a, psV////PS1, psa////PS1
hide label, dist_V-a
dist dist_V-r, psV////PS1, psr////PS1
hide label, dist_V-r
dist dist_W-B, psW////PS1, psB////PS1
hide label, dist_W-B
dist dist_W-I, psW////PS1, psI////PS1
hide label, dist_W-I
dist dist_W-Q, psW////PS1, psQ////PS1
hide label, dist_W-Q
dist dist_W-S, psW////PS1, psS////PS1
hide label, dist_W-S
dist dist_W-U, psW////PS1, psU////PS1
hide label, dist_W-U
dist dist_W-j, psW////PS1, psj////PS1
hide label, dist_W-j
dist dist_W-m, psW////PS1, psm////PS1
hide label, dist_W-m
dist dist_W-a, psW////PS1, psa////PS1
hide label, dist_W-a
dist dist_W-s, psW////PS1, pss////PS1
hide label, dist_W-s
dist dist_W-u, psW////PS1, psu////PS1
hide label, dist_W-u
dist dist_W-h, psW////PS1, psh////PS1
hide label, dist_W-h
dist dist_X-l, psX////PS1, psl////PS1
hide label, dist_X-l
dist dist_X-Y, psX////PS1, psY////PS1
hide label, dist_X-Y
dist dist_X-Z, psX////PS1, psZ////PS1
hide label, dist_X-Z
dist dist_X-p, psX////PS1, psp////PS1
hide label, dist_X-p
dist dist_X-b, psX////PS1, psb////PS1
hide label, dist_X-b
dist dist_X-c, psX////PS1, psc////PS1
hide label, dist_X-c
dist dist_Y-X, psY////PS1, psX////PS1
hide label, dist_Y-X
dist dist_Y-l, psY////PS1, psl////PS1
hide label, dist_Y-l
dist dist_Y-Z, psY////PS1, psZ////PS1
hide label, dist_Y-Z
dist dist_Y-c, psY////PS1, psc////PS1
hide label, dist_Y-c
dist dist_Y-v, psY////PS1, psv////PS1
hide label, dist_Y-v
dist dist_Z-X, psZ////PS1, psX////PS1
hide label, dist_Z-X
dist dist_Z-Y, psZ////PS1, psY////PS1
hide label, dist_Z-Y
dist dist_Z-l, psZ////PS1, psl////PS1
hide label, dist_Z-l
dist dist_Z-p, psZ////PS1, psp////PS1
hide label, dist_Z-p
dist dist_a-V, psa////PS1, psV////PS1
hide label, dist_a-V
dist dist_a-W, psa////PS1, psW////PS1
hide label, dist_a-W
dist dist_a-i, psa////PS1, psi////PS1
hide label, dist_a-i
dist dist_a-l, psa////PS1, psl////PS1
hide label, dist_a-l
dist dist_a-n, psa////PS1, psn////PS1
hide label, dist_a-n
dist dist_a-p, psa////PS1, psp////PS1
hide label, dist_a-p
dist dist_a-r, psa////PS1, psr////PS1
hide label, dist_a-r
dist dist_a-u, psa////PS1, psu////PS1
hide label, dist_a-u
dist dist_a-b, psa////PS1, psb////PS1
hide label, dist_a-b
dist dist_a-d, psa////PS1, psd////PS1
hide label, dist_a-d
dist dist_a-e, psa////PS1, pse////PS1
hide label, dist_a-e
dist dist_a-g, psa////PS1, psg////PS1
hide label, dist_a-g
dist dist_a-h, psa////PS1, psh////PS1
hide label, dist_a-h
dist dist_b-X, psb////PS1, psX////PS1
hide label, dist_b-X
dist dist_b-a, psb////PS1, psa////PS1
hide label, dist_b-a
dist dist_b-l, psb////PS1, psl////PS1
hide label, dist_b-l
dist dist_b-p, psb////PS1, psp////PS1
hide label, dist_b-p
dist dist_b-d, psb////PS1, psd////PS1
hide label, dist_b-d
dist dist_b-v, psb////PS1, psv////PS1
hide label, dist_b-v
dist dist_c-X, psc////PS1, psX////PS1
hide label, dist_c-X
dist dist_c-Y, psc////PS1, psY////PS1
hide label, dist_c-Y
dist dist_c-l, psc////PS1, psl////PS1
hide label, dist_c-l
dist dist_c-o, psc////PS1, pso////PS1
hide label, dist_c-o
dist dist_c-p, psc////PS1, psp////PS1
hide label, dist_c-p
dist dist_c-r, psc////PS1, psr////PS1
hide label, dist_c-r
dist dist_c-d, psc////PS1, psd////PS1
hide label, dist_c-d
dist dist_c-v, psc////PS1, psv////PS1
hide label, dist_c-v
dist dist_d-a, psd////PS1, psa////PS1
hide label, dist_d-a
dist dist_d-b, psd////PS1, psb////PS1
hide label, dist_d-b
dist dist_d-c, psd////PS1, psc////PS1
hide label, dist_d-c
dist dist_d-l, psd////PS1, psl////PS1
hide label, dist_d-l
dist dist_d-n, psd////PS1, psn////PS1
hide label, dist_d-n
dist dist_d-o, psd////PS1, pso////PS1
hide label, dist_d-o
dist dist_d-r, psd////PS1, psr////PS1
hide label, dist_d-r
dist dist_d-v, psd////PS1, psv////PS1
hide label, dist_d-v
dist dist_d-e, psd////PS1, pse////PS1
hide label, dist_d-e
dist dist_d-g, psd////PS1, psg////PS1
hide label, dist_d-g
dist dist_e-Q, pse////PS1, psQ////PS1
hide label, dist_e-Q
dist dist_e-a, pse////PS1, psa////PS1
hide label, dist_e-a
dist dist_e-d, pse////PS1, psd////PS1
hide label, dist_e-d
dist dist_e-l, pse////PS1, psl////PS1
hide label, dist_e-l
dist dist_e-o, pse////PS1, pso////PS1
hide label, dist_e-o
dist dist_e-p, pse////PS1, psp////PS1
hide label, dist_e-p
dist dist_e-r, pse////PS1, psr////PS1
hide label, dist_e-r
dist dist_e-w, pse////PS1, psw////PS1
hide label, dist_e-w
dist dist_e-g, pse////PS1, psg////PS1
hide label, dist_e-g
dist dist_f-w, psf////PS1, psw////PS1
hide label, dist_f-w
dist dist_f-g, psf////PS1, psg////PS1
hide label, dist_f-g
dist dist_g-a, psg////PS1, psa////PS1
hide label, dist_g-a
dist dist_g-d, psg////PS1, psd////PS1
hide label, dist_g-d
dist dist_g-e, psg////PS1, pse////PS1
hide label, dist_g-e
dist dist_g-f, psg////PS1, psf////PS1
hide label, dist_g-f
dist dist_g-i, psg////PS1, psi////PS1
hide label, dist_g-i
dist dist_g-o, psg////PS1, pso////PS1
hide label, dist_g-o
dist dist_g-r, psg////PS1, psr////PS1
hide label, dist_g-r
dist dist_g-u, psg////PS1, psu////PS1
hide label, dist_g-u
dist dist_g-w, psg////PS1, psw////PS1
hide label, dist_g-w
dist dist_g-h, psg////PS1, psh////PS1
hide label, dist_g-h
dist dist_h-S, psh////PS1, psS////PS1
hide label, dist_h-S
dist dist_h-U, psh////PS1, psU////PS1
hide label, dist_h-U
dist dist_h-W, psh////PS1, psW////PS1
hide label, dist_h-W
dist dist_h-a, psh////PS1, psa////PS1
hide label, dist_h-a
dist dist_h-g, psh////PS1, psg////PS1
hide label, dist_h-g
dist dist_h-i, psh////PS1, psi////PS1
hide label, dist_h-i
dist dist_h-k, psh////PS1, psk////PS1
hide label, dist_h-k
dist dist_h-m, psh////PS1, psm////PS1
hide label, dist_h-m
dist dist_h-s, psh////PS1, pss////PS1
hide label, dist_h-s
dist dist_h-u, psh////PS1, psu////PS1
hide label, dist_h-u
dist dist_i-Q, psi////PS1, psQ////PS1
hide label, dist_i-Q
dist dist_i-V, psi////PS1, psV////PS1
hide label, dist_i-V
dist dist_i-a, psi////PS1, psa////PS1
hide label, dist_i-a
dist dist_i-g, psi////PS1, psg////PS1
hide label, dist_i-g
dist dist_i-h, psi////PS1, psh////PS1
hide label, dist_i-h
dist dist_i-j, psi////PS1, psj////PS1
hide label, dist_i-j
dist dist_i-k, psi////PS1, psk////PS1
hide label, dist_i-k
dist dist_i-l, psi////PS1, psl////PS1
hide label, dist_i-l
dist dist_i-m, psi////PS1, psm////PS1
hide label, dist_i-m
dist dist_i-r, psi////PS1, psr////PS1
hide label, dist_i-r
dist dist_i-u, psi////PS1, psu////PS1
hide label, dist_i-u
dist dist_i-w, psi////PS1, psw////PS1
hide label, dist_i-w
dist dist_j-C, psj////PS1, psC////PS1
hide label, dist_j-C
dist dist_j-E, psj////PS1, psE////PS1
hide label, dist_j-E
dist dist_j-J, psj////PS1, psJ////PS1
hide label, dist_j-J
dist dist_j-Q, psj////PS1, psQ////PS1
hide label, dist_j-Q
dist dist_j-U, psj////PS1, psU////PS1
hide label, dist_j-U
dist dist_j-W, psj////PS1, psW////PS1
hide label, dist_j-W
dist dist_j-i, psj////PS1, psi////PS1
hide label, dist_j-i
dist dist_j-k, psj////PS1, psk////PS1
hide label, dist_j-k
dist dist_j-m, psj////PS1, psm////PS1
hide label, dist_j-m
dist dist_j-s, psj////PS1, pss////PS1
hide label, dist_j-s
dist dist_k-Q, psk////PS1, psQ////PS1
hide label, dist_k-Q
dist dist_k-V, psk////PS1, psV////PS1
hide label, dist_k-V
dist dist_k-h, psk////PS1, psh////PS1
hide label, dist_k-h
dist dist_k-i, psk////PS1, psi////PS1
hide label, dist_k-i
dist dist_k-j, psk////PS1, psj////PS1
hide label, dist_k-j
dist dist_k-l, psk////PS1, psl////PS1
hide label, dist_k-l
dist dist_k-m, psk////PS1, psm////PS1
hide label, dist_k-m
dist dist_l-Q, psl////PS1, psQ////PS1
hide label, dist_l-Q
dist dist_l-V, psl////PS1, psV////PS1
hide label, dist_l-V
dist dist_l-X, psl////PS1, psX////PS1
hide label, dist_l-X
dist dist_l-Y, psl////PS1, psY////PS1
hide label, dist_l-Y
dist dist_l-Z, psl////PS1, psZ////PS1
hide label, dist_l-Z
dist dist_l-a, psl////PS1, psa////PS1
hide label, dist_l-a
dist dist_l-b, psl////PS1, psb////PS1
hide label, dist_l-b
dist dist_l-c, psl////PS1, psc////PS1
hide label, dist_l-c
dist dist_l-d, psl////PS1, psd////PS1
hide label, dist_l-d
dist dist_l-e, psl////PS1, pse////PS1
hide label, dist_l-e
dist dist_l-i, psl////PS1, psi////PS1
hide label, dist_l-i
dist dist_l-k, psl////PS1, psk////PS1
hide label, dist_l-k
dist dist_l-o, psl////PS1, pso////PS1
hide label, dist_l-o
dist dist_l-p, psl////PS1, psp////PS1
hide label, dist_l-p
dist dist_l-r, psl////PS1, psr////PS1
hide label, dist_l-r
dist dist_l-v, psl////PS1, psv////PS1
hide label, dist_l-v
dist dist_m-J, psm////PS1, psJ////PS1
hide label, dist_m-J
dist dist_m-S, psm////PS1, psS////PS1
hide label, dist_m-S
dist dist_m-U, psm////PS1, psU////PS1
hide label, dist_m-U
dist dist_m-W, psm////PS1, psW////PS1
hide label, dist_m-W
dist dist_m-h, psm////PS1, psh////PS1
hide label, dist_m-h
dist dist_m-i, psm////PS1, psi////PS1
hide label, dist_m-i
dist dist_m-j, psm////PS1, psj////PS1
hide label, dist_m-j
dist dist_m-k, psm////PS1, psk////PS1
hide label, dist_m-k
dist dist_m-s, psm////PS1, pss////PS1
hide label, dist_m-s
dist dist_n-a, psn////PS1, psa////PS1
hide label, dist_n-a
dist dist_n-d, psn////PS1, psd////PS1
hide label, dist_n-d
dist dist_n-r, psn////PS1, psr////PS1
hide label, dist_n-r
dist dist_o-c, pso////PS1, psc////PS1
hide label, dist_o-c
dist dist_o-d, pso////PS1, psd////PS1
hide label, dist_o-d
dist dist_o-e, pso////PS1, pse////PS1
hide label, dist_o-e
dist dist_o-g, pso////PS1, psg////PS1
hide label, dist_o-g
dist dist_o-l, pso////PS1, psl////PS1
hide label, dist_o-l
dist dist_o-p, pso////PS1, psp////PS1
hide label, dist_o-p
dist dist_o-r, pso////PS1, psr////PS1
hide label, dist_o-r
dist dist_p-X, psp////PS1, psX////PS1
hide label, dist_p-X
dist dist_p-Z, psp////PS1, psZ////PS1
hide label, dist_p-Z
dist dist_p-a, psp////PS1, psa////PS1
hide label, dist_p-a
dist dist_p-b, psp////PS1, psb////PS1
hide label, dist_p-b
dist dist_p-c, psp////PS1, psc////PS1
hide label, dist_p-c
dist dist_p-e, psp////PS1, pse////PS1
hide label, dist_p-e
dist dist_p-l, psp////PS1, psl////PS1
hide label, dist_p-l
dist dist_p-o, psp////PS1, pso////PS1
hide label, dist_p-o
dist dist_p-r, psp////PS1, psr////PS1
hide label, dist_p-r
dist dist_r-Q, psr////PS1, psQ////PS1
hide label, dist_r-Q
dist dist_r-V, psr////PS1, psV////PS1
hide label, dist_r-V
dist dist_r-a, psr////PS1, psa////PS1
hide label, dist_r-a
dist dist_r-c, psr////PS1, psc////PS1
hide label, dist_r-c
dist dist_r-d, psr////PS1, psd////PS1
hide label, dist_r-d
dist dist_r-e, psr////PS1, pse////PS1
hide label, dist_r-e
dist dist_r-g, psr////PS1, psg////PS1
hide label, dist_r-g
dist dist_r-i, psr////PS1, psi////PS1
hide label, dist_r-i
dist dist_r-l, psr////PS1, psl////PS1
hide label, dist_r-l
dist dist_r-n, psr////PS1, psn////PS1
hide label, dist_r-n
dist dist_r-o, psr////PS1, pso////PS1
hide label, dist_r-o
dist dist_r-p, psr////PS1, psp////PS1
hide label, dist_r-p
dist dist_r-u, psr////PS1, psu////PS1
hide label, dist_r-u
dist dist_r-w, psr////PS1, psw////PS1
hide label, dist_r-w
dist dist_s-B, pss////PS1, psB////PS1
hide label, dist_s-B
dist dist_s-C, pss////PS1, psC////PS1
hide label, dist_s-C
dist dist_s-J, pss////PS1, psJ////PS1
hide label, dist_s-J
dist dist_s-N, pss////PS1, psN////PS1
hide label, dist_s-N
dist dist_s-Q, pss////PS1, psQ////PS1
hide label, dist_s-Q
dist dist_s-S, pss////PS1, psS////PS1
hide label, dist_s-S
dist dist_s-U, pss////PS1, psU////PS1
hide label, dist_s-U
dist dist_s-W, pss////PS1, psW////PS1
hide label, dist_s-W
dist dist_s-h, pss////PS1, psh////PS1
hide label, dist_s-h
dist dist_s-j, pss////PS1, psj////PS1
hide label, dist_s-j
dist dist_s-m, pss////PS1, psm////PS1
hide label, dist_s-m
dist dist_s-u, pss////PS1, psu////PS1
hide label, dist_s-u
dist dist_u-S, psu////PS1, psS////PS1
hide label, dist_u-S
dist dist_u-U, psu////PS1, psU////PS1
hide label, dist_u-U
dist dist_u-W, psu////PS1, psW////PS1
hide label, dist_u-W
dist dist_u-a, psu////PS1, psa////PS1
hide label, dist_u-a
dist dist_u-g, psu////PS1, psg////PS1
hide label, dist_u-g
dist dist_u-h, psu////PS1, psh////PS1
hide label, dist_u-h
dist dist_u-i, psu////PS1, psi////PS1
hide label, dist_u-i
dist dist_u-r, psu////PS1, psr////PS1
hide label, dist_u-r
dist dist_u-s, psu////PS1, pss////PS1
hide label, dist_u-s
dist dist_v-Y, psv////PS1, psY////PS1
hide label, dist_v-Y
dist dist_v-b, psv////PS1, psb////PS1
hide label, dist_v-b
dist dist_v-c, psv////PS1, psc////PS1
hide label, dist_v-c
dist dist_v-d, psv////PS1, psd////PS1
hide label, dist_v-d
dist dist_v-l, psv////PS1, psl////PS1
hide label, dist_v-l
dist dist_w-H, psw////PS1, psH////PS1
hide label, dist_w-H
dist dist_w-Q, psw////PS1, psQ////PS1
hide label, dist_w-Q
dist dist_w-e, psw////PS1, pse////PS1
hide label, dist_w-e
dist dist_w-f, psw////PS1, psf////PS1
hide label, dist_w-f
dist dist_w-g, psw////PS1, psg////PS1
hide label, dist_w-g
dist dist_w-i, psw////PS1, psi////PS1
hide label, dist_w-i
dist dist_w-r, psw////PS1, psr////PS1
hide label, dist_w-r

# Group all distance measurements originating from one chain to all other connected ones
group chain-A_contacts, dist_A-*
group chain-B_contacts, dist_B-*
group chain-C_contacts, dist_C-*
group chain-E_contacts, dist_E-*
group chain-F_contacts, dist_F-*
group chain-G_contacts, dist_G-*
group chain-H_contacts, dist_H-*
group chain-I_contacts, dist_I-*
group chain-J_contacts, dist_J-*
group chain-K_contacts, dist_K-*
group chain-L_contacts, dist_L-*
group chain-M_contacts, dist_M-*
group chain-N_contacts, dist_N-*
group chain-O_contacts, dist_O-*
group chain-P_contacts, dist_P-*
group chain-Q_contacts, dist_Q-*
group chain-S_contacts, dist_S-*
group chain-T_contacts, dist_T-*
group chain-U_contacts, dist_U-*
group chain-V_contacts, dist_V-*
group chain-W_contacts, dist_W-*
group chain-X_contacts, dist_X-*
group chain-Y_contacts, dist_Y-*
group chain-Z_contacts, dist_Z-*
group chain-a_contacts, dist_a-*
group chain-b_contacts, dist_b-*
group chain-c_contacts, dist_c-*
group chain-d_contacts, dist_d-*
group chain-e_contacts, dist_e-*
group chain-f_contacts, dist_f-*
group chain-g_contacts, dist_g-*
group chain-h_contacts, dist_h-*
group chain-i_contacts, dist_i-*
group chain-j_contacts, dist_j-*
group chain-k_contacts, dist_k-*
group chain-l_contacts, dist_l-*
group chain-m_contacts, dist_m-*
group chain-n_contacts, dist_n-*
group chain-o_contacts, dist_o-*
group chain-p_contacts, dist_p-*
group chain-r_contacts, dist_r-*
group chain-s_contacts, dist_s-*
group chain-u_contacts, dist_u-*
group chain-v_contacts, dist_v-*
group chain-w_contacts, dist_w-*

# Disallow addition/removal from the group (allow with 'open')
group chain-A_contacts, close
group chain-B_contacts, close
group chain-C_contacts, close
group chain-E_contacts, close
group chain-F_contacts, close
group chain-G_contacts, close
group chain-H_contacts, close
group chain-I_contacts, close
group chain-J_contacts, close
group chain-K_contacts, close
group chain-L_contacts, close
group chain-M_contacts, close
group chain-N_contacts, close
group chain-O_contacts, close
group chain-P_contacts, close
group chain-Q_contacts, close
group chain-S_contacts, close
group chain-T_contacts, close
group chain-U_contacts, close
group chain-V_contacts, close
group chain-W_contacts, close
group chain-X_contacts, close
group chain-Y_contacts, close
group chain-Z_contacts, close
group chain-a_contacts, close
group chain-b_contacts, close
group chain-c_contacts, close
group chain-d_contacts, close
group chain-e_contacts, close
group chain-f_contacts, close
group chain-g_contacts, close
group chain-h_contacts, close
group chain-i_contacts, close
group chain-j_contacts, close
group chain-k_contacts, close
group chain-l_contacts, close
group chain-m_contacts, close
group chain-n_contacts, close
group chain-o_contacts, close
group chain-p_contacts, close
group chain-r_contacts, close
group chain-s_contacts, close
group chain-u_contacts, close
group chain-v_contacts, close
group chain-w_contacts, close

# Set pseudoatoms along the distance measurement line and label them with the corresponding edge weights
pseudoatom ps_A-I, psA or psI, label=6
pseudoatom ps_A-K, psA or psK, label=39
pseudoatom ps_A-L, psA or psL, label=19
pseudoatom ps_A-M, psA or psM, label=48
pseudoatom ps_A-O, psA or psO, label=105
pseudoatom ps_B-C, psB or psC, label=44
pseudoatom ps_B-H, psB or psH, label=1
pseudoatom ps_B-I, psB or psI, label=27
pseudoatom ps_B-J, psB or psJ, label=2
pseudoatom ps_B-s, psB or pss, label=18
pseudoatom ps_B-L, psB or psL, label=9
pseudoatom ps_B-M, psB or psM, label=16
pseudoatom ps_B-N, psB or psN, label=77
pseudoatom ps_B-P, psB or psP, label=10
pseudoatom ps_B-Q, psB or psQ, label=66
pseudoatom ps_B-S, psB or psS, label=4
pseudoatom ps_B-T, psB or psT, label=25
pseudoatom ps_B-U, psB or psU, label=7
pseudoatom ps_B-W, psB or psW, label=7
pseudoatom ps_C-B, psC or psB, label=44
pseudoatom ps_C-j, psC or psj, label=13
pseudoatom ps_C-E, psC or psE, label=1
pseudoatom ps_C-J, psC or psJ, label=17
pseudoatom ps_C-s, psC or pss, label=42
pseudoatom ps_C-N, psC or psN, label=7
pseudoatom ps_C-P, psC or psP, label=21
pseudoatom ps_C-Q, psC or psQ, label=51
pseudoatom ps_E-C, psE or psC, label=1
pseudoatom ps_E-j, psE or psj, label=8
pseudoatom ps_E-G, psE or psG, label=21
pseudoatom ps_E-J, psE or psJ, label=13
pseudoatom ps_E-L, psE or psL, label=10
pseudoatom ps_E-M, psE or psM, label=14
pseudoatom ps_E-P, psE or psP, label=34
pseudoatom ps_E-Q, psE or psQ, label=2
pseudoatom ps_F-M, psF or psM, label=53
pseudoatom ps_G-E, psG or psE, label=21
pseudoatom ps_H-B, psH or psB, label=1
pseudoatom ps_H-I, psH or psI, label=9
pseudoatom ps_H-w, psH or psw, label=2
pseudoatom ps_H-P, psH or psP, label=50
pseudoatom ps_H-Q, psH or psQ, label=16
pseudoatom ps_I-A, psI or psA, label=6
pseudoatom ps_I-B, psI or psB, label=27
pseudoatom ps_I-H, psI or psH, label=9
pseudoatom ps_I-W, psI or psW, label=17
pseudoatom ps_I-M, psI or psM, label=14
pseudoatom ps_I-N, psI or psN, label=14
pseudoatom ps_I-P, psI or psP, label=50
pseudoatom ps_I-Q, psI or psQ, label=47
pseudoatom ps_J-B, psJ or psB, label=2
pseudoatom ps_J-C, psJ or psC, label=17
pseudoatom ps_J-E, psJ or psE, label=13
pseudoatom ps_J-T, psJ or psT, label=10
pseudoatom ps_J-j, psJ or psj, label=4
pseudoatom ps_J-m, psJ or psm, label=3
pseudoatom ps_J-s, psJ or pss, label=2
pseudoatom ps_J-L, psJ or psL, label=7
pseudoatom ps_J-M, psJ or psM, label=1
pseudoatom ps_J-P, psJ or psP, label=17
pseudoatom ps_K-A, psK or psA, label=39
pseudoatom ps_K-L, psK or psL, label=5
pseudoatom ps_K-M, psK or psM, label=2
pseudoatom ps_K-O, psK or psO, label=27
pseudoatom ps_L-A, psL or psA, label=19
pseudoatom ps_L-B, psL or psB, label=9
pseudoatom ps_L-E, psL or psE, label=10
pseudoatom ps_L-J, psL or psJ, label=7
pseudoatom ps_L-K, psL or psK, label=5
pseudoatom ps_L-T, psL or psT, label=1
pseudoatom ps_L-M, psL or psM, label=43
pseudoatom ps_L-N, psL or psN, label=3
pseudoatom ps_L-O, psL or psO, label=3
pseudoatom ps_L-P, psL or psP, label=43
pseudoatom ps_M-A, psM or psA, label=48
pseudoatom ps_M-B, psM or psB, label=16
pseudoatom ps_M-E, psM or psE, label=14
pseudoatom ps_M-F, psM or psF, label=53
pseudoatom ps_M-I, psM or psI, label=14
pseudoatom ps_M-J, psM or psJ, label=1
pseudoatom ps_M-K, psM or psK, label=2
pseudoatom ps_M-L, psM or psL, label=43
pseudoatom ps_M-N, psM or psN, label=23
pseudoatom ps_M-O, psM or psO, label=32
pseudoatom ps_M-P, psM or psP, label=16
pseudoatom ps_M-Q, psM or psQ, label=35
pseudoatom ps_M-T, psM or psT, label=14
pseudoatom ps_N-B, psN or psB, label=77
pseudoatom ps_N-C, psN or psC, label=7
pseudoatom ps_N-I, psN or psI, label=14
pseudoatom ps_N-L, psN or psL, label=3
pseudoatom ps_N-M, psN or psM, label=23
pseudoatom ps_N-T, psN or psT, label=25
pseudoatom ps_N-s, psN or pss, label=5
pseudoatom ps_N-Q, psN or psQ, label=1
pseudoatom ps_O-A, psO or psA, label=105
pseudoatom ps_O-K, psO or psK, label=27
pseudoatom ps_O-L, psO or psL, label=3
pseudoatom ps_O-M, psO or psM, label=32
pseudoatom ps_O-T, psO or psT, label=1
pseudoatom ps_P-B, psP or psB, label=10
pseudoatom ps_P-C, psP or psC, label=21
pseudoatom ps_P-E, psP or psE, label=34
pseudoatom ps_P-H, psP or psH, label=50
pseudoatom ps_P-I, psP or psI, label=50
pseudoatom ps_P-J, psP or psJ, label=17
pseudoatom ps_P-L, psP or psL, label=43
pseudoatom ps_P-M, psP or psM, label=16
pseudoatom ps_P-Q, psP or psQ, label=70
pseudoatom ps_Q-B, psQ or psB, label=66
pseudoatom ps_Q-C, psQ or psC, label=51
pseudoatom ps_Q-E, psQ or psE, label=2
pseudoatom ps_Q-H, psQ or psH, label=16
pseudoatom ps_Q-I, psQ or psI, label=47
pseudoatom ps_Q-M, psQ or psM, label=35
pseudoatom ps_Q-N, psQ or psN, label=1
pseudoatom ps_Q-P, psQ or psP, label=70
pseudoatom ps_Q-i, psQ or psi, label=24
pseudoatom ps_Q-j, psQ or psj, label=19
pseudoatom ps_Q-k, psQ or psk, label=5
pseudoatom ps_Q-l, psQ or psl, label=6
pseudoatom ps_Q-r, psQ or psr, label=12
pseudoatom ps_Q-s, psQ or pss, label=35
pseudoatom ps_Q-w, psQ or psw, label=8
pseudoatom ps_Q-T, psQ or psT, label=7
pseudoatom ps_Q-W, psQ or psW, label=30
pseudoatom ps_Q-e, psQ or pse, label=3
pseudoatom ps_S-B, psS or psB, label=4
pseudoatom ps_S-W, psS or psW, label=12
pseudoatom ps_S-m, psS or psm, label=8
pseudoatom ps_S-s, psS or pss, label=36
pseudoatom ps_S-u, psS or psu, label=34
pseudoatom ps_S-h, psS or psh, label=9
pseudoatom ps_T-B, psT or psB, label=25
pseudoatom ps_T-J, psT or psJ, label=10
pseudoatom ps_T-L, psT or psL, label=1
pseudoatom ps_T-M, psT or psM, label=14
pseudoatom ps_T-N, psT or psN, label=25
pseudoatom ps_T-O, psT or psO, label=1
pseudoatom ps_T-Q, psT or psQ, label=7
pseudoatom ps_U-B, psU or psB, label=7
pseudoatom ps_U-j, psU or psj, label=9
pseudoatom ps_U-W, psU or psW, label=13
pseudoatom ps_U-m, psU or psm, label=2
pseudoatom ps_U-s, psU or pss, label=16
pseudoatom ps_U-u, psU or psu, label=42
pseudoatom ps_U-h, psU or psh, label=3
pseudoatom ps_V-i, psV or psi, label=13
pseudoatom ps_V-k, psV or psk, label=1
pseudoatom ps_V-l, psV or psl, label=10
pseudoatom ps_V-a, psV or psa, label=4
pseudoatom ps_V-r, psV or psr, label=5
pseudoatom ps_W-B, psW or psB, label=7
pseudoatom ps_W-I, psW or psI, label=17
pseudoatom ps_W-Q, psW or psQ, label=30
pseudoatom ps_W-S, psW or psS, label=12
pseudoatom ps_W-U, psW or psU, label=13
pseudoatom ps_W-j, psW or psj, label=1
pseudoatom ps_W-m, psW or psm, label=17
pseudoatom ps_W-a, psW or psa, label=2
pseudoatom ps_W-s, psW or pss, label=30
pseudoatom ps_W-u, psW or psu, label=46
pseudoatom ps_W-h, psW or psh, label=36
pseudoatom ps_X-l, psX or psl, label=3
pseudoatom ps_X-Y, psX or psY, label=3
pseudoatom ps_X-Z, psX or psZ, label=22
pseudoatom ps_X-p, psX or psp, label=26
pseudoatom ps_X-b, psX or psb, label=13
pseudoatom ps_X-c, psX or psc, label=2
pseudoatom ps_Y-X, psY or psX, label=3
pseudoatom ps_Y-l, psY or psl, label=34
pseudoatom ps_Y-Z, psY or psZ, label=26
pseudoatom ps_Y-c, psY or psc, label=1
pseudoatom ps_Y-v, psY or psv, label=20
pseudoatom ps_Z-X, psZ or psX, label=22
pseudoatom ps_Z-Y, psZ or psY, label=26
pseudoatom ps_Z-l, psZ or psl, label=22
pseudoatom ps_Z-p, psZ or psp, label=13
pseudoatom ps_a-V, psa or psV, label=4
pseudoatom ps_a-W, psa or psW, label=2
pseudoatom ps_a-i, psa or psi, label=5
pseudoatom ps_a-l, psa or psl, label=3
pseudoatom ps_a-n, psa or psn, label=17
pseudoatom ps_a-p, psa or psp, label=16
pseudoatom ps_a-r, psa or psr, label=27
pseudoatom ps_a-u, psa or psu, label=16
pseudoatom ps_a-b, psa or psb, label=10
pseudoatom ps_a-d, psa or psd, label=45
pseudoatom ps_a-e, psa or pse, label=25
pseudoatom ps_a-g, psa or psg, label=26
pseudoatom ps_a-h, psa or psh, label=38
pseudoatom ps_b-X, psb or psX, label=13
pseudoatom ps_b-a, psb or psa, label=10
pseudoatom ps_b-l, psb or psl, label=27
pseudoatom ps_b-p, psb or psp, label=31
pseudoatom ps_b-d, psb or psd, label=37
pseudoatom ps_b-v, psb or psv, label=20
pseudoatom ps_c-X, psc or psX, label=2
pseudoatom ps_c-Y, psc or psY, label=1
pseudoatom ps_c-l, psc or psl, label=38
pseudoatom ps_c-o, psc or pso, label=65
pseudoatom ps_c-p, psc or psp, label=8
pseudoatom ps_c-r, psc or psr, label=7
pseudoatom ps_c-d, psc or psd, label=4
pseudoatom ps_c-v, psc or psv, label=39
pseudoatom ps_d-a, psd or psa, label=45
pseudoatom ps_d-b, psd or psb, label=37
pseudoatom ps_d-c, psd or psc, label=4
pseudoatom ps_d-l, psd or psl, label=26
pseudoatom ps_d-n, psd or psn, label=7
pseudoatom ps_d-o, psd or pso, label=12
pseudoatom ps_d-r, psd or psr, label=23
pseudoatom ps_d-v, psd or psv, label=9
pseudoatom ps_d-e, psd or pse, label=43
pseudoatom ps_d-g, psd or psg, label=32
pseudoatom ps_e-Q, pse or psQ, label=3
pseudoatom ps_e-a, pse or psa, label=25
pseudoatom ps_e-d, pse or psd, label=43
pseudoatom ps_e-l, pse or psl, label=1
pseudoatom ps_e-o, pse or pso, label=1
pseudoatom ps_e-p, pse or psp, label=3
pseudoatom ps_e-r, pse or psr, label=49
pseudoatom ps_e-w, pse or psw, label=5
pseudoatom ps_e-g, pse or psg, label=5
pseudoatom ps_f-w, psf or psw, label=13
pseudoatom ps_f-g, psf or psg, label=22
pseudoatom ps_g-a, psg or psa, label=26
pseudoatom ps_g-d, psg or psd, label=32
pseudoatom ps_g-e, psg or pse, label=5
pseudoatom ps_g-f, psg or psf, label=22
pseudoatom ps_g-i, psg or psi, label=29
pseudoatom ps_g-o, psg or pso, label=2
pseudoatom ps_g-r, psg or psr, label=4
pseudoatom ps_g-u, psg or psu, label=24
pseudoatom ps_g-w, psg or psw, label=10
pseudoatom ps_g-h, psg or psh, label=5
pseudoatom ps_h-S, psh or psS, label=9
pseudoatom ps_h-U, psh or psU, label=3
pseudoatom ps_h-W, psh or psW, label=36
pseudoatom ps_h-a, psh or psa, label=38
pseudoatom ps_h-g, psh or psg, label=5
pseudoatom ps_h-i, psh or psi, label=25
pseudoatom ps_h-k, psh or psk, label=14
pseudoatom ps_h-m, psh or psm, label=15
pseudoatom ps_h-s, psh or pss, label=1
pseudoatom ps_h-u, psh or psu, label=16
pseudoatom ps_i-Q, psi or psQ, label=24
pseudoatom ps_i-V, psi or psV, label=13
pseudoatom ps_i-a, psi or psa, label=5
pseudoatom ps_i-g, psi or psg, label=29
pseudoatom ps_i-h, psi or psh, label=25
pseudoatom ps_i-j, psi or psj, label=1
pseudoatom ps_i-k, psi or psk, label=38
pseudoatom ps_i-l, psi or psl, label=21
pseudoatom ps_i-m, psi or psm, label=14
pseudoatom ps_i-r, psi or psr, label=25
pseudoatom ps_i-u, psi or psu, label=8
pseudoatom ps_i-w, psi or psw, label=55
pseudoatom ps_j-C, psj or psC, label=13
pseudoatom ps_j-E, psj or psE, label=8
pseudoatom ps_j-J, psj or psJ, label=4
pseudoatom ps_j-Q, psj or psQ, label=19
pseudoatom ps_j-U, psj or psU, label=9
pseudoatom ps_j-W, psj or psW, label=1
pseudoatom ps_j-i, psj or psi, label=1
pseudoatom ps_j-k, psj or psk, label=6
pseudoatom ps_j-m, psj or psm, label=47
pseudoatom ps_j-s, psj or pss, label=66
pseudoatom ps_k-Q, psk or psQ, label=5
pseudoatom ps_k-V, psk or psV, label=1
pseudoatom ps_k-h, psk or psh, label=14
pseudoatom ps_k-i, psk or psi, label=38
pseudoatom ps_k-j, psk or psj, label=6
pseudoatom ps_k-l, psk or psl, label=13
pseudoatom ps_k-m, psk or psm, label=87
pseudoatom ps_l-Q, psl or psQ, label=6
pseudoatom ps_l-V, psl or psV, label=10
pseudoatom ps_l-X, psl or psX, label=3
pseudoatom ps_l-Y, psl or psY, label=34
pseudoatom ps_l-Z, psl or psZ, label=22
pseudoatom ps_l-a, psl or psa, label=3
pseudoatom ps_l-b, psl or psb, label=27
pseudoatom ps_l-c, psl or psc, label=38
pseudoatom ps_l-d, psl or psd, label=26
pseudoatom ps_l-e, psl or pse, label=1
pseudoatom ps_l-i, psl or psi, label=21
pseudoatom ps_l-k, psl or psk, label=13
pseudoatom ps_l-o, psl or pso, label=25
pseudoatom ps_l-p, psl or psp, label=43
pseudoatom ps_l-r, psl or psr, label=52
pseudoatom ps_l-v, psl or psv, label=22
pseudoatom ps_m-J, psm or psJ, label=3
pseudoatom ps_m-S, psm or psS, label=8
pseudoatom ps_m-U, psm or psU, label=2
pseudoatom ps_m-W, psm or psW, label=17
pseudoatom ps_m-h, psm or psh, label=15
pseudoatom ps_m-i, psm or psi, label=14
pseudoatom ps_m-j, psm or psj, label=47
pseudoatom ps_m-k, psm or psk, label=87
pseudoatom ps_m-s, psm or pss, label=23
pseudoatom ps_n-a, psn or psa, label=17
pseudoatom ps_n-d, psn or psd, label=7
pseudoatom ps_n-r, psn or psr, label=17
pseudoatom ps_o-c, pso or psc, label=65
pseudoatom ps_o-d, pso or psd, label=12
pseudoatom ps_o-e, pso or pse, label=1
pseudoatom ps_o-g, pso or psg, label=2
pseudoatom ps_o-l, pso or psl, label=25
pseudoatom ps_o-p, pso or psp, label=37
pseudoatom ps_o-r, pso or psr, label=33
pseudoatom ps_p-X, psp or psX, label=26
pseudoatom ps_p-Z, psp or psZ, label=13
pseudoatom ps_p-a, psp or psa, label=16
pseudoatom ps_p-b, psp or psb, label=31
pseudoatom ps_p-c, psp or psc, label=8
pseudoatom ps_p-e, psp or pse, label=3
pseudoatom ps_p-l, psp or psl, label=43
pseudoatom ps_p-o, psp or pso, label=37
pseudoatom ps_p-r, psp or psr, label=10
pseudoatom ps_r-Q, psr or psQ, label=12
pseudoatom ps_r-V, psr or psV, label=5
pseudoatom ps_r-a, psr or psa, label=27
pseudoatom ps_r-c, psr or psc, label=7
pseudoatom ps_r-d, psr or psd, label=23
pseudoatom ps_r-e, psr or pse, label=49
pseudoatom ps_r-g, psr or psg, label=4
pseudoatom ps_r-i, psr or psi, label=25
pseudoatom ps_r-l, psr or psl, label=52
pseudoatom ps_r-n, psr or psn, label=17
pseudoatom ps_r-o, psr or pso, label=33
pseudoatom ps_r-p, psr or psp, label=10
pseudoatom ps_r-u, psr or psu, label=7
pseudoatom ps_r-w, psr or psw, label=3
pseudoatom ps_s-B, pss or psB, label=18
pseudoatom ps_s-C, pss or psC, label=42
pseudoatom ps_s-J, pss or psJ, label=2
pseudoatom ps_s-N, pss or psN, label=5
pseudoatom ps_s-Q, pss or psQ, label=35
pseudoatom ps_s-S, pss or psS, label=36
pseudoatom ps_s-U, pss or psU, label=16
pseudoatom ps_s-W, pss or psW, label=30
pseudoatom ps_s-h, pss or psh, label=1
pseudoatom ps_s-j, pss or psj, label=66
pseudoatom ps_s-m, pss or psm, label=23
pseudoatom ps_s-u, pss or psu, label=3
pseudoatom ps_u-S, psu or psS, label=34
pseudoatom ps_u-U, psu or psU, label=42
pseudoatom ps_u-W, psu or psW, label=46
pseudoatom ps_u-a, psu or psa, label=16
pseudoatom ps_u-g, psu or psg, label=24
pseudoatom ps_u-h, psu or psh, label=16
pseudoatom ps_u-i, psu or psi, label=8
pseudoatom ps_u-r, psu or psr, label=7
pseudoatom ps_u-s, psu or pss, label=3
pseudoatom ps_v-Y, psv or psY, label=20
pseudoatom ps_v-b, psv or psb, label=20
pseudoatom ps_v-c, psv or psc, label=39
pseudoatom ps_v-d, psv or psd, label=9
pseudoatom ps_v-l, psv or psl, label=22
pseudoatom ps_w-H, psw or psH, label=2
pseudoatom ps_w-Q, psw or psQ, label=8
pseudoatom ps_w-e, psw or pse, label=5
pseudoatom ps_w-f, psw or psf, label=13
pseudoatom ps_w-g, psw or psg, label=10
pseudoatom ps_w-i, psw or psi, label=55
pseudoatom ps_w-r, psw or psr, label=3

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

	