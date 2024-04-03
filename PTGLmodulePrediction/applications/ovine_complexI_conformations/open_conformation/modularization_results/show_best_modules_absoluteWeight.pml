
### PML file to represent the clusters in the structure.
# Python is by default case insensitive for chain IDs but for large clusters with upper and lower case letters we have to care!
set ignore_case, off

# Load the structure
load 7zd6.cif, complex

# First color everything white
color white, all
# Select and color the found clusters
select cluster0, chain A+H+J+K+N+Y+k+l+m+o+q+x+z
select cluster1, chain L+M+V+W+X+Z+n+p+r+s+t+u+v+w+y
select cluster2, chain 1+2+a
select cluster3, chain 3+4+5+6+9+b+c+d+e+f+g+h+i+j

set_color newColor0, [0.0, 1.0, 0.0] 
color newColor0, cluster0
set_color newColor1, [0.0, 0.0, 1.0] 
color newColor1, cluster1
set_color newColor2, [1.0, 0.0, 0.0] 
color newColor2, cluster2
set_color newColor3, [1.0, 1.0, 0.0] 
color newColor3, cluster3

# Create pseudoatoms at the geometric center of the chains
pseudoatom psA, chain A
pseudoatom psH, chain H
pseudoatom psJ, chain J
pseudoatom psK, chain K
pseudoatom psL, chain L
pseudoatom psM, chain M
pseudoatom psN, chain N
pseudoatom psV, chain V
pseudoatom psW, chain W
pseudoatom psX, chain X
pseudoatom psY, chain Y
pseudoatom psZ, chain Z
pseudoatom psk, chain k
pseudoatom psl, chain l
pseudoatom psm, chain m
pseudoatom psn, chain n
pseudoatom pso, chain o
pseudoatom psp, chain p
pseudoatom psq, chain q
pseudoatom psr, chain r
pseudoatom pss, chain s
pseudoatom pst, chain t
pseudoatom psu, chain u
pseudoatom psv, chain v
pseudoatom psw, chain w
pseudoatom psx, chain x
pseudoatom psy, chain y
pseudoatom psz, chain z
pseudoatom ps1, chain 1
pseudoatom ps2, chain 2
pseudoatom ps3, chain 3
pseudoatom ps4, chain 4
pseudoatom ps5, chain 5
pseudoatom ps6, chain 6
pseudoatom ps9, chain 9
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

group pseudoatoms, ps*

# Create all distance objects originating from each chain (reverse dublicates are needed!)
# Hide the distance label, generate a pseudoatom between both pseudoatoms for which the distance was calculated
# and label it according to the edge weight
dist dist_A-H, psA////PS1, psH////PS1
hide label, dist_A-H
dist dist_A-6, psA////PS1, ps6////PS1
hide label, dist_A-6
dist dist_A-J, psA////PS1, psJ////PS1
hide label, dist_A-J
dist dist_A-q, psA////PS1, psq////PS1
hide label, dist_A-q
dist dist_A-K, psA////PS1, psK////PS1
hide label, dist_A-K
dist dist_A-N, psA////PS1, psN////PS1
hide label, dist_A-N
dist dist_A-d, psA////PS1, psd////PS1
hide label, dist_A-d
dist dist_A-g, psA////PS1, psg////PS1
hide label, dist_A-g
dist dist_A-m, psA////PS1, psm////PS1
hide label, dist_A-m
dist dist_A-4, psA////PS1, ps4////PS1
hide label, dist_A-4
dist dist_H-A, psH////PS1, psA////PS1
hide label, dist_H-A
dist dist_H-6, psH////PS1, ps6////PS1
hide label, dist_H-6
dist dist_H-J, psH////PS1, psJ////PS1
hide label, dist_H-J
dist dist_H-q, psH////PS1, psq////PS1
hide label, dist_H-q
dist dist_H-9, psH////PS1, ps9////PS1
hide label, dist_H-9
dist dist_H-Y, psH////PS1, psY////PS1
hide label, dist_H-Y
dist dist_H-h, psH////PS1, psh////PS1
hide label, dist_H-h
dist dist_H-z, psH////PS1, psz////PS1
hide label, dist_H-z
dist dist_H-i, psH////PS1, psi////PS1
hide label, dist_H-i
dist dist_H-m, psH////PS1, psm////PS1
hide label, dist_H-m
dist dist_H-4, psH////PS1, ps4////PS1
hide label, dist_H-4
dist dist_J-A, psJ////PS1, psA////PS1
hide label, dist_J-A
dist dist_J-H, psJ////PS1, psH////PS1
hide label, dist_J-H
dist dist_J-q, psJ////PS1, psq////PS1
hide label, dist_J-q
dist dist_J-K, psJ////PS1, psK////PS1
hide label, dist_J-K
dist dist_J-N, psJ////PS1, psN////PS1
hide label, dist_J-N
dist dist_J-z, psJ////PS1, psz////PS1
hide label, dist_J-z
dist dist_J-k, psJ////PS1, psk////PS1
hide label, dist_J-k
dist dist_J-l, psJ////PS1, psl////PS1
hide label, dist_J-l
dist dist_J-m, psJ////PS1, psm////PS1
hide label, dist_J-m
dist dist_K-A, psK////PS1, psA////PS1
hide label, dist_K-A
dist dist_K-J, psK////PS1, psJ////PS1
hide label, dist_K-J
dist dist_K-L, psK////PS1, psL////PS1
hide label, dist_K-L
dist dist_K-N, psK////PS1, psN////PS1
hide label, dist_K-N
dist dist_K-l, psK////PS1, psl////PS1
hide label, dist_K-l
dist dist_K-4, psK////PS1, ps4////PS1
hide label, dist_K-4
dist dist_L-K, psL////PS1, psK////PS1
hide label, dist_L-K
dist dist_L-M, psL////PS1, psM////PS1
hide label, dist_L-M
dist dist_L-N, psL////PS1, psN////PS1
hide label, dist_L-N
dist dist_L-V, psL////PS1, psV////PS1
hide label, dist_L-V
dist dist_L-W, psL////PS1, psW////PS1
hide label, dist_L-W
dist dist_L-X, psL////PS1, psX////PS1
hide label, dist_L-X
dist dist_L-Z, psL////PS1, psZ////PS1
hide label, dist_L-Z
dist dist_L-n, psL////PS1, psn////PS1
hide label, dist_L-n
dist dist_L-p, psL////PS1, psp////PS1
hide label, dist_L-p
dist dist_L-r, psL////PS1, psr////PS1
hide label, dist_L-r
dist dist_L-s, psL////PS1, pss////PS1
hide label, dist_L-s
dist dist_L-t, psL////PS1, pst////PS1
hide label, dist_L-t
dist dist_L-u, psL////PS1, psu////PS1
hide label, dist_L-u
dist dist_L-v, psL////PS1, psv////PS1
hide label, dist_L-v
dist dist_L-4, psL////PS1, ps4////PS1
hide label, dist_L-4
dist dist_M-L, psM////PS1, psL////PS1
hide label, dist_M-L
dist dist_M-N, psM////PS1, psN////PS1
hide label, dist_M-N
dist dist_M-V, psM////PS1, psV////PS1
hide label, dist_M-V
dist dist_M-W, psM////PS1, psW////PS1
hide label, dist_M-W
dist dist_M-Y, psM////PS1, psY////PS1
hide label, dist_M-Y
dist dist_M-Z, psM////PS1, psZ////PS1
hide label, dist_M-Z
dist dist_M-k, psM////PS1, psk////PS1
hide label, dist_M-k
dist dist_M-o, psM////PS1, pso////PS1
hide label, dist_M-o
dist dist_M-p, psM////PS1, psp////PS1
hide label, dist_M-p
dist dist_M-t, psM////PS1, pst////PS1
hide label, dist_M-t
dist dist_M-v, psM////PS1, psv////PS1
hide label, dist_M-v
dist dist_M-w, psM////PS1, psw////PS1
hide label, dist_M-w
dist dist_M-y, psM////PS1, psy////PS1
hide label, dist_M-y
dist dist_M-4, psM////PS1, ps4////PS1
hide label, dist_M-4
dist dist_N-A, psN////PS1, psA////PS1
hide label, dist_N-A
dist dist_N-J, psN////PS1, psJ////PS1
hide label, dist_N-J
dist dist_N-K, psN////PS1, psK////PS1
hide label, dist_N-K
dist dist_N-L, psN////PS1, psL////PS1
hide label, dist_N-L
dist dist_N-M, psN////PS1, psM////PS1
hide label, dist_N-M
dist dist_N-o, psN////PS1, pso////PS1
hide label, dist_N-o
dist dist_N-V, psN////PS1, psV////PS1
hide label, dist_N-V
dist dist_N-W, psN////PS1, psW////PS1
hide label, dist_N-W
dist dist_N-Y, psN////PS1, psY////PS1
hide label, dist_N-Y
dist dist_N-k, psN////PS1, psk////PS1
hide label, dist_N-k
dist dist_N-l, psN////PS1, psl////PS1
hide label, dist_N-l
dist dist_N-4, psN////PS1, ps4////PS1
hide label, dist_N-4
dist dist_V-L, psV////PS1, psL////PS1
hide label, dist_V-L
dist dist_V-M, psV////PS1, psM////PS1
hide label, dist_V-M
dist dist_V-N, psV////PS1, psN////PS1
hide label, dist_V-N
dist dist_V-W, psV////PS1, psW////PS1
hide label, dist_V-W
dist dist_W-L, psW////PS1, psL////PS1
hide label, dist_W-L
dist dist_W-M, psW////PS1, psM////PS1
hide label, dist_W-M
dist dist_W-N, psW////PS1, psN////PS1
hide label, dist_W-N
dist dist_W-V, psW////PS1, psV////PS1
hide label, dist_W-V
dist dist_W-X, psW////PS1, psX////PS1
hide label, dist_W-X
dist dist_W-Y, psW////PS1, psY////PS1
hide label, dist_W-Y
dist dist_W-Z, psW////PS1, psZ////PS1
hide label, dist_W-Z
dist dist_W-l, psW////PS1, psl////PS1
hide label, dist_W-l
dist dist_W-o, psW////PS1, pso////PS1
hide label, dist_W-o
dist dist_W-q, psW////PS1, psq////PS1
hide label, dist_W-q
dist dist_W-r, psW////PS1, psr////PS1
hide label, dist_W-r
dist dist_W-t, psW////PS1, pst////PS1
hide label, dist_W-t
dist dist_W-w, psW////PS1, psw////PS1
hide label, dist_W-w
dist dist_W-y, psW////PS1, psy////PS1
hide label, dist_W-y
dist dist_X-L, psX////PS1, psL////PS1
hide label, dist_X-L
dist dist_X-W, psX////PS1, psW////PS1
hide label, dist_X-W
dist dist_X-r, psX////PS1, psr////PS1
hide label, dist_X-r
dist dist_X-t, psX////PS1, pst////PS1
hide label, dist_X-t
dist dist_X-u, psX////PS1, psu////PS1
hide label, dist_X-u
dist dist_X-v, psX////PS1, psv////PS1
hide label, dist_X-v
dist dist_X-n, psX////PS1, psn////PS1
hide label, dist_X-n
dist dist_Y-H, psY////PS1, psH////PS1
hide label, dist_Y-H
dist dist_Y-M, psY////PS1, psM////PS1
hide label, dist_Y-M
dist dist_Y-N, psY////PS1, psN////PS1
hide label, dist_Y-N
dist dist_Y-W, psY////PS1, psW////PS1
hide label, dist_Y-W
dist dist_Y-o, psY////PS1, pso////PS1
hide label, dist_Y-o
dist dist_Y-q, psY////PS1, psq////PS1
hide label, dist_Y-q
dist dist_Y-z, psY////PS1, psz////PS1
hide label, dist_Y-z
dist dist_Y-l, psY////PS1, psl////PS1
hide label, dist_Y-l
dist dist_Y-m, psY////PS1, psm////PS1
hide label, dist_Y-m
dist dist_Z-L, psZ////PS1, psL////PS1
hide label, dist_Z-L
dist dist_Z-M, psZ////PS1, psM////PS1
hide label, dist_Z-M
dist dist_Z-W, psZ////PS1, psW////PS1
hide label, dist_Z-W
dist dist_Z-o, psZ////PS1, pso////PS1
hide label, dist_Z-o
dist dist_Z-p, psZ////PS1, psp////PS1
hide label, dist_Z-p
dist dist_Z-r, psZ////PS1, psr////PS1
hide label, dist_Z-r
dist dist_Z-s, psZ////PS1, pss////PS1
hide label, dist_Z-s
dist dist_Z-v, psZ////PS1, psv////PS1
hide label, dist_Z-v
dist dist_Z-w, psZ////PS1, psw////PS1
hide label, dist_Z-w
dist dist_Z-y, psZ////PS1, psy////PS1
hide label, dist_Z-y
dist dist_k-J, psk////PS1, psJ////PS1
hide label, dist_k-J
dist dist_k-M, psk////PS1, psM////PS1
hide label, dist_k-M
dist dist_k-N, psk////PS1, psN////PS1
hide label, dist_k-N
dist dist_k-o, psk////PS1, pso////PS1
hide label, dist_k-o
dist dist_k-w, psk////PS1, psw////PS1
hide label, dist_k-w
dist dist_k-f, psk////PS1, psf////PS1
hide label, dist_k-f
dist dist_k-x, psk////PS1, psx////PS1
hide label, dist_k-x
dist dist_k-4, psk////PS1, ps4////PS1
hide label, dist_k-4
dist dist_l-J, psl////PS1, psJ////PS1
hide label, dist_l-J
dist dist_l-K, psl////PS1, psK////PS1
hide label, dist_l-K
dist dist_l-N, psl////PS1, psN////PS1
hide label, dist_l-N
dist dist_l-W, psl////PS1, psW////PS1
hide label, dist_l-W
dist dist_l-Y, psl////PS1, psY////PS1
hide label, dist_l-Y
dist dist_l-o, psl////PS1, pso////PS1
hide label, dist_l-o
dist dist_l-q, psl////PS1, psq////PS1
hide label, dist_l-q
dist dist_l-z, psl////PS1, psz////PS1
hide label, dist_l-z
dist dist_l-m, psl////PS1, psm////PS1
hide label, dist_l-m
dist dist_m-A, psm////PS1, psA////PS1
hide label, dist_m-A
dist dist_m-H, psm////PS1, psH////PS1
hide label, dist_m-H
dist dist_m-J, psm////PS1, psJ////PS1
hide label, dist_m-J
dist dist_m-Y, psm////PS1, psY////PS1
hide label, dist_m-Y
dist dist_m-l, psm////PS1, psl////PS1
hide label, dist_m-l
dist dist_m-q, psm////PS1, psq////PS1
hide label, dist_m-q
dist dist_m-9, psm////PS1, ps9////PS1
hide label, dist_m-9
dist dist_n-L, psn////PS1, psL////PS1
hide label, dist_n-L
dist dist_n-X, psn////PS1, psX////PS1
hide label, dist_n-X
dist dist_n-t, psn////PS1, pst////PS1
hide label, dist_n-t
dist dist_n-u, psn////PS1, psu////PS1
hide label, dist_n-u
dist dist_o-M, pso////PS1, psM////PS1
hide label, dist_o-M
dist dist_o-N, pso////PS1, psN////PS1
hide label, dist_o-N
dist dist_o-W, pso////PS1, psW////PS1
hide label, dist_o-W
dist dist_o-Y, pso////PS1, psY////PS1
hide label, dist_o-Y
dist dist_o-Z, pso////PS1, psZ////PS1
hide label, dist_o-Z
dist dist_o-k, pso////PS1, psk////PS1
hide label, dist_o-k
dist dist_o-l, pso////PS1, psl////PS1
hide label, dist_o-l
dist dist_o-p, pso////PS1, psp////PS1
hide label, dist_o-p
dist dist_o-w, pso////PS1, psw////PS1
hide label, dist_o-w
dist dist_o-x, pso////PS1, psx////PS1
hide label, dist_o-x
dist dist_p-L, psp////PS1, psL////PS1
hide label, dist_p-L
dist dist_p-M, psp////PS1, psM////PS1
hide label, dist_p-M
dist dist_p-Z, psp////PS1, psZ////PS1
hide label, dist_p-Z
dist dist_p-o, psp////PS1, pso////PS1
hide label, dist_p-o
dist dist_p-t, psp////PS1, pst////PS1
hide label, dist_p-t
dist dist_p-v, psp////PS1, psv////PS1
hide label, dist_p-v
dist dist_q-A, psq////PS1, psA////PS1
hide label, dist_q-A
dist dist_q-H, psq////PS1, psH////PS1
hide label, dist_q-H
dist dist_q-J, psq////PS1, psJ////PS1
hide label, dist_q-J
dist dist_q-W, psq////PS1, psW////PS1
hide label, dist_q-W
dist dist_q-Y, psq////PS1, psY////PS1
hide label, dist_q-Y
dist dist_q-l, psq////PS1, psl////PS1
hide label, dist_q-l
dist dist_q-m, psq////PS1, psm////PS1
hide label, dist_q-m
dist dist_q-9, psq////PS1, ps9////PS1
hide label, dist_q-9
dist dist_q-h, psq////PS1, psh////PS1
hide label, dist_q-h
dist dist_q-z, psq////PS1, psz////PS1
hide label, dist_q-z
dist dist_q-4, psq////PS1, ps4////PS1
hide label, dist_q-4
dist dist_r-L, psr////PS1, psL////PS1
hide label, dist_r-L
dist dist_r-W, psr////PS1, psW////PS1
hide label, dist_r-W
dist dist_r-X, psr////PS1, psX////PS1
hide label, dist_r-X
dist dist_r-Z, psr////PS1, psZ////PS1
hide label, dist_r-Z
dist dist_r-s, psr////PS1, pss////PS1
hide label, dist_r-s
dist dist_r-t, psr////PS1, pst////PS1
hide label, dist_r-t
dist dist_s-L, pss////PS1, psL////PS1
hide label, dist_s-L
dist dist_s-Z, pss////PS1, psZ////PS1
hide label, dist_s-Z
dist dist_s-r, pss////PS1, psr////PS1
hide label, dist_s-r
dist dist_s-u, pss////PS1, psu////PS1
hide label, dist_s-u
dist dist_s-v, pss////PS1, psv////PS1
hide label, dist_s-v
dist dist_t-L, pst////PS1, psL////PS1
hide label, dist_t-L
dist dist_t-M, pst////PS1, psM////PS1
hide label, dist_t-M
dist dist_t-W, pst////PS1, psW////PS1
hide label, dist_t-W
dist dist_t-X, pst////PS1, psX////PS1
hide label, dist_t-X
dist dist_t-n, pst////PS1, psn////PS1
hide label, dist_t-n
dist dist_t-p, pst////PS1, psp////PS1
hide label, dist_t-p
dist dist_t-r, pst////PS1, psr////PS1
hide label, dist_t-r
dist dist_t-v, pst////PS1, psv////PS1
hide label, dist_t-v
dist dist_t-w, pst////PS1, psw////PS1
hide label, dist_t-w
dist dist_u-L, psu////PS1, psL////PS1
hide label, dist_u-L
dist dist_u-X, psu////PS1, psX////PS1
hide label, dist_u-X
dist dist_u-n, psu////PS1, psn////PS1
hide label, dist_u-n
dist dist_u-s, psu////PS1, pss////PS1
hide label, dist_u-s
dist dist_v-L, psv////PS1, psL////PS1
hide label, dist_v-L
dist dist_v-M, psv////PS1, psM////PS1
hide label, dist_v-M
dist dist_v-X, psv////PS1, psX////PS1
hide label, dist_v-X
dist dist_v-Z, psv////PS1, psZ////PS1
hide label, dist_v-Z
dist dist_v-p, psv////PS1, psp////PS1
hide label, dist_v-p
dist dist_v-s, psv////PS1, pss////PS1
hide label, dist_v-s
dist dist_v-t, psv////PS1, pst////PS1
hide label, dist_v-t
dist dist_w-M, psw////PS1, psM////PS1
hide label, dist_w-M
dist dist_w-W, psw////PS1, psW////PS1
hide label, dist_w-W
dist dist_w-Z, psw////PS1, psZ////PS1
hide label, dist_w-Z
dist dist_w-k, psw////PS1, psk////PS1
hide label, dist_w-k
dist dist_w-o, psw////PS1, pso////PS1
hide label, dist_w-o
dist dist_w-t, psw////PS1, pst////PS1
hide label, dist_w-t
dist dist_w-4, psw////PS1, ps4////PS1
hide label, dist_w-4
dist dist_x-k, psx////PS1, psk////PS1
hide label, dist_x-k
dist dist_x-o, psx////PS1, pso////PS1
hide label, dist_x-o
dist dist_y-M, psy////PS1, psM////PS1
hide label, dist_y-M
dist dist_y-W, psy////PS1, psW////PS1
hide label, dist_y-W
dist dist_y-Z, psy////PS1, psZ////PS1
hide label, dist_y-Z
dist dist_z-H, psz////PS1, psH////PS1
hide label, dist_z-H
dist dist_z-J, psz////PS1, psJ////PS1
hide label, dist_z-J
dist dist_z-Y, psz////PS1, psY////PS1
hide label, dist_z-Y
dist dist_z-l, psz////PS1, psl////PS1
hide label, dist_z-l
dist dist_z-q, psz////PS1, psq////PS1
hide label, dist_z-q
dist dist_z-9, psz////PS1, ps9////PS1
hide label, dist_z-9
dist dist_z-h, psz////PS1, psh////PS1
hide label, dist_z-h
dist dist_1-a, ps1////PS1, psa////PS1
hide label, dist_1-a
dist dist_1-c, ps1////PS1, psc////PS1
hide label, dist_1-c
dist dist_1-h, ps1////PS1, psh////PS1
hide label, dist_1-h
dist dist_1-2, ps1////PS1, ps2////PS1
hide label, dist_1-2
dist dist_1-3, ps1////PS1, ps3////PS1
hide label, dist_1-3
dist dist_2-1, ps2////PS1, ps1////PS1
hide label, dist_2-1
dist dist_2-a, ps2////PS1, psa////PS1
hide label, dist_2-a
dist dist_2-b, ps2////PS1, psb////PS1
hide label, dist_2-b
dist dist_2-c, ps2////PS1, psc////PS1
hide label, dist_2-c
dist dist_2-3, ps2////PS1, ps3////PS1
hide label, dist_2-3
dist dist_3-1, ps3////PS1, ps1////PS1
hide label, dist_3-1
dist dist_3-2, ps3////PS1, ps2////PS1
hide label, dist_3-2
dist dist_3-5, ps3////PS1, ps5////PS1
hide label, dist_3-5
dist dist_3-9, ps3////PS1, ps9////PS1
hide label, dist_3-9
dist dist_3-a, ps3////PS1, psa////PS1
hide label, dist_3-a
dist dist_3-b, ps3////PS1, psb////PS1
hide label, dist_3-b
dist dist_3-c, ps3////PS1, psc////PS1
hide label, dist_3-c
dist dist_3-d, ps3////PS1, psd////PS1
hide label, dist_3-d
dist dist_3-e, ps3////PS1, pse////PS1
hide label, dist_3-e
dist dist_3-g, ps3////PS1, psg////PS1
hide label, dist_3-g
dist dist_3-h, ps3////PS1, psh////PS1
hide label, dist_3-h
dist dist_3-i, ps3////PS1, psi////PS1
hide label, dist_3-i
dist dist_3-4, ps3////PS1, ps4////PS1
hide label, dist_3-4
dist dist_4-A, ps4////PS1, psA////PS1
hide label, dist_4-A
dist dist_4-H, ps4////PS1, psH////PS1
hide label, dist_4-H
dist dist_4-K, ps4////PS1, psK////PS1
hide label, dist_4-K
dist dist_4-L, ps4////PS1, psL////PS1
hide label, dist_4-L
dist dist_4-M, ps4////PS1, psM////PS1
hide label, dist_4-M
dist dist_4-N, ps4////PS1, psN////PS1
hide label, dist_4-N
dist dist_4-k, ps4////PS1, psk////PS1
hide label, dist_4-k
dist dist_4-q, ps4////PS1, psq////PS1
hide label, dist_4-q
dist dist_4-w, ps4////PS1, psw////PS1
hide label, dist_4-w
dist dist_4-3, ps4////PS1, ps3////PS1
hide label, dist_4-3
dist dist_4-5, ps4////PS1, ps5////PS1
hide label, dist_4-5
dist dist_4-6, ps4////PS1, ps6////PS1
hide label, dist_4-6
dist dist_4-9, ps4////PS1, ps9////PS1
hide label, dist_4-9
dist dist_4-b, ps4////PS1, psb////PS1
hide label, dist_4-b
dist dist_4-f, ps4////PS1, psf////PS1
hide label, dist_4-f
dist dist_4-g, ps4////PS1, psg////PS1
hide label, dist_4-g
dist dist_4-h, ps4////PS1, psh////PS1
hide label, dist_4-h
dist dist_5-3, ps5////PS1, ps3////PS1
hide label, dist_5-3
dist dist_5-4, ps5////PS1, ps4////PS1
hide label, dist_5-4
dist dist_5-6, ps5////PS1, ps6////PS1
hide label, dist_5-6
dist dist_5-9, ps5////PS1, ps9////PS1
hide label, dist_5-9
dist dist_5-c, ps5////PS1, psc////PS1
hide label, dist_5-c
dist dist_5-d, ps5////PS1, psd////PS1
hide label, dist_5-d
dist dist_5-f, ps5////PS1, psf////PS1
hide label, dist_5-f
dist dist_5-g, ps5////PS1, psg////PS1
hide label, dist_5-g
dist dist_5-h, ps5////PS1, psh////PS1
hide label, dist_5-h
dist dist_6-A, ps6////PS1, psA////PS1
hide label, dist_6-A
dist dist_6-H, ps6////PS1, psH////PS1
hide label, dist_6-H
dist dist_6-4, ps6////PS1, ps4////PS1
hide label, dist_6-4
dist dist_6-5, ps6////PS1, ps5////PS1
hide label, dist_6-5
dist dist_6-9, ps6////PS1, ps9////PS1
hide label, dist_6-9
dist dist_6-d, ps6////PS1, psd////PS1
hide label, dist_6-d
dist dist_6-h, ps6////PS1, psh////PS1
hide label, dist_6-h
dist dist_6-i, ps6////PS1, psi////PS1
hide label, dist_6-i
dist dist_9-H, ps9////PS1, psH////PS1
hide label, dist_9-H
dist dist_9-m, ps9////PS1, psm////PS1
hide label, dist_9-m
dist dist_9-q, ps9////PS1, psq////PS1
hide label, dist_9-q
dist dist_9-z, ps9////PS1, psz////PS1
hide label, dist_9-z
dist dist_9-3, ps9////PS1, ps3////PS1
hide label, dist_9-3
dist dist_9-4, ps9////PS1, ps4////PS1
hide label, dist_9-4
dist dist_9-5, ps9////PS1, ps5////PS1
hide label, dist_9-5
dist dist_9-6, ps9////PS1, ps6////PS1
hide label, dist_9-6
dist dist_9-b, ps9////PS1, psb////PS1
hide label, dist_9-b
dist dist_9-c, ps9////PS1, psc////PS1
hide label, dist_9-c
dist dist_9-d, ps9////PS1, psd////PS1
hide label, dist_9-d
dist dist_9-h, ps9////PS1, psh////PS1
hide label, dist_9-h
dist dist_9-i, ps9////PS1, psi////PS1
hide label, dist_9-i
dist dist_a-1, psa////PS1, ps1////PS1
hide label, dist_a-1
dist dist_a-2, psa////PS1, ps2////PS1
hide label, dist_a-2
dist dist_a-3, psa////PS1, ps3////PS1
hide label, dist_a-3
dist dist_a-c, psa////PS1, psc////PS1
hide label, dist_a-c
dist dist_b-2, psb////PS1, ps2////PS1
hide label, dist_b-2
dist dist_b-3, psb////PS1, ps3////PS1
hide label, dist_b-3
dist dist_b-4, psb////PS1, ps4////PS1
hide label, dist_b-4
dist dist_b-9, psb////PS1, ps9////PS1
hide label, dist_b-9
dist dist_b-c, psb////PS1, psc////PS1
hide label, dist_b-c
dist dist_b-d, psb////PS1, psd////PS1
hide label, dist_b-d
dist dist_b-i, psb////PS1, psi////PS1
hide label, dist_b-i
dist dist_c-1, psc////PS1, ps1////PS1
hide label, dist_c-1
dist dist_c-2, psc////PS1, ps2////PS1
hide label, dist_c-2
dist dist_c-3, psc////PS1, ps3////PS1
hide label, dist_c-3
dist dist_c-5, psc////PS1, ps5////PS1
hide label, dist_c-5
dist dist_c-9, psc////PS1, ps9////PS1
hide label, dist_c-9
dist dist_c-a, psc////PS1, psa////PS1
hide label, dist_c-a
dist dist_c-b, psc////PS1, psb////PS1
hide label, dist_c-b
dist dist_c-d, psc////PS1, psd////PS1
hide label, dist_c-d
dist dist_c-f, psc////PS1, psf////PS1
hide label, dist_c-f
dist dist_c-g, psc////PS1, psg////PS1
hide label, dist_c-g
dist dist_c-i, psc////PS1, psi////PS1
hide label, dist_c-i
dist dist_d-A, psd////PS1, psA////PS1
hide label, dist_d-A
dist dist_d-3, psd////PS1, ps3////PS1
hide label, dist_d-3
dist dist_d-5, psd////PS1, ps5////PS1
hide label, dist_d-5
dist dist_d-6, psd////PS1, ps6////PS1
hide label, dist_d-6
dist dist_d-9, psd////PS1, ps9////PS1
hide label, dist_d-9
dist dist_d-b, psd////PS1, psb////PS1
hide label, dist_d-b
dist dist_d-c, psd////PS1, psc////PS1
hide label, dist_d-c
dist dist_d-g, psd////PS1, psg////PS1
hide label, dist_d-g
dist dist_d-i, psd////PS1, psi////PS1
hide label, dist_d-i
dist dist_e-3, pse////PS1, ps3////PS1
hide label, dist_e-3
dist dist_f-k, psf////PS1, psk////PS1
hide label, dist_f-k
dist dist_f-4, psf////PS1, ps4////PS1
hide label, dist_f-4
dist dist_f-5, psf////PS1, ps5////PS1
hide label, dist_f-5
dist dist_f-c, psf////PS1, psc////PS1
hide label, dist_f-c
dist dist_f-h, psf////PS1, psh////PS1
hide label, dist_f-h
dist dist_g-A, psg////PS1, psA////PS1
hide label, dist_g-A
dist dist_g-3, psg////PS1, ps3////PS1
hide label, dist_g-3
dist dist_g-4, psg////PS1, ps4////PS1
hide label, dist_g-4
dist dist_g-5, psg////PS1, ps5////PS1
hide label, dist_g-5
dist dist_g-c, psg////PS1, psc////PS1
hide label, dist_g-c
dist dist_g-d, psg////PS1, psd////PS1
hide label, dist_g-d
dist dist_g-j, psg////PS1, psj////PS1
hide label, dist_g-j
dist dist_h-H, psh////PS1, psH////PS1
hide label, dist_h-H
dist dist_h-q, psh////PS1, psq////PS1
hide label, dist_h-q
dist dist_h-z, psh////PS1, psz////PS1
hide label, dist_h-z
dist dist_h-1, psh////PS1, ps1////PS1
hide label, dist_h-1
dist dist_h-3, psh////PS1, ps3////PS1
hide label, dist_h-3
dist dist_h-4, psh////PS1, ps4////PS1
hide label, dist_h-4
dist dist_h-5, psh////PS1, ps5////PS1
hide label, dist_h-5
dist dist_h-6, psh////PS1, ps6////PS1
hide label, dist_h-6
dist dist_h-9, psh////PS1, ps9////PS1
hide label, dist_h-9
dist dist_h-f, psh////PS1, psf////PS1
hide label, dist_h-f
dist dist_h-i, psh////PS1, psi////PS1
hide label, dist_h-i
dist dist_i-H, psi////PS1, psH////PS1
hide label, dist_i-H
dist dist_i-3, psi////PS1, ps3////PS1
hide label, dist_i-3
dist dist_i-6, psi////PS1, ps6////PS1
hide label, dist_i-6
dist dist_i-9, psi////PS1, ps9////PS1
hide label, dist_i-9
dist dist_i-b, psi////PS1, psb////PS1
hide label, dist_i-b
dist dist_i-c, psi////PS1, psc////PS1
hide label, dist_i-c
dist dist_i-d, psi////PS1, psd////PS1
hide label, dist_i-d
dist dist_i-h, psi////PS1, psh////PS1
hide label, dist_i-h
dist dist_j-g, psj////PS1, psg////PS1
hide label, dist_j-g

# Group all distance measurements originating from one chain to all other connected ones
group chain-A_contacts, dist_A-*
group chain-H_contacts, dist_H-*
group chain-J_contacts, dist_J-*
group chain-K_contacts, dist_K-*
group chain-L_contacts, dist_L-*
group chain-M_contacts, dist_M-*
group chain-N_contacts, dist_N-*
group chain-V_contacts, dist_V-*
group chain-W_contacts, dist_W-*
group chain-X_contacts, dist_X-*
group chain-Y_contacts, dist_Y-*
group chain-Z_contacts, dist_Z-*
group chain-k_contacts, dist_k-*
group chain-l_contacts, dist_l-*
group chain-m_contacts, dist_m-*
group chain-n_contacts, dist_n-*
group chain-o_contacts, dist_o-*
group chain-p_contacts, dist_p-*
group chain-q_contacts, dist_q-*
group chain-r_contacts, dist_r-*
group chain-s_contacts, dist_s-*
group chain-t_contacts, dist_t-*
group chain-u_contacts, dist_u-*
group chain-v_contacts, dist_v-*
group chain-w_contacts, dist_w-*
group chain-x_contacts, dist_x-*
group chain-y_contacts, dist_y-*
group chain-z_contacts, dist_z-*
group chain-1_contacts, dist_1-*
group chain-2_contacts, dist_2-*
group chain-3_contacts, dist_3-*
group chain-4_contacts, dist_4-*
group chain-5_contacts, dist_5-*
group chain-6_contacts, dist_6-*
group chain-9_contacts, dist_9-*
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

# Disallow addition/removal from the group (allow with 'open')
group chain-A_contacts, close
group chain-H_contacts, close
group chain-J_contacts, close
group chain-K_contacts, close
group chain-L_contacts, close
group chain-M_contacts, close
group chain-N_contacts, close
group chain-V_contacts, close
group chain-W_contacts, close
group chain-X_contacts, close
group chain-Y_contacts, close
group chain-Z_contacts, close
group chain-k_contacts, close
group chain-l_contacts, close
group chain-m_contacts, close
group chain-n_contacts, close
group chain-o_contacts, close
group chain-p_contacts, close
group chain-q_contacts, close
group chain-r_contacts, close
group chain-s_contacts, close
group chain-t_contacts, close
group chain-u_contacts, close
group chain-v_contacts, close
group chain-w_contacts, close
group chain-x_contacts, close
group chain-y_contacts, close
group chain-z_contacts, close
group chain-1_contacts, close
group chain-2_contacts, close
group chain-3_contacts, close
group chain-4_contacts, close
group chain-5_contacts, close
group chain-6_contacts, close
group chain-9_contacts, close
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

# Set pseudoatoms along the distance measurement line and label them with the corresponding edge weights
pseudoatom ps_A-H, psA or psH, label=75
pseudoatom ps_A-6, psA or ps6, label=6
pseudoatom ps_A-J, psA or psJ, label=53
pseudoatom ps_A-q, psA or psq, label=2
pseudoatom ps_A-K, psA or psK, label=13
pseudoatom ps_A-N, psA or psN, label=1
pseudoatom ps_A-d, psA or psd, label=3
pseudoatom ps_A-g, psA or psg, label=2
pseudoatom ps_A-m, psA or psm, label=9
pseudoatom ps_A-4, psA or ps4, label=4
pseudoatom ps_H-A, psH or psA, label=75
pseudoatom ps_H-6, psH or ps6, label=38
pseudoatom ps_H-J, psH or psJ, label=37
pseudoatom ps_H-q, psH or psq, label=33
pseudoatom ps_H-9, psH or ps9, label=24
pseudoatom ps_H-Y, psH or psY, label=6
pseudoatom ps_H-h, psH or psh, label=1
pseudoatom ps_H-z, psH or psz, label=46
pseudoatom ps_H-i, psH or psi, label=2
pseudoatom ps_H-m, psH or psm, label=27
pseudoatom ps_H-4, psH or ps4, label=33
pseudoatom ps_J-A, psJ or psA, label=53
pseudoatom ps_J-H, psJ or psH, label=37
pseudoatom ps_J-q, psJ or psq, label=26
pseudoatom ps_J-K, psJ or psK, label=92
pseudoatom ps_J-N, psJ or psN, label=23
pseudoatom ps_J-z, psJ or psz, label=8
pseudoatom ps_J-k, psJ or psk, label=3
pseudoatom ps_J-l, psJ or psl, label=19
pseudoatom ps_J-m, psJ or psm, label=1
pseudoatom ps_K-A, psK or psA, label=13
pseudoatom ps_K-J, psK or psJ, label=92
pseudoatom ps_K-L, psK or psL, label=12
pseudoatom ps_K-N, psK or psN, label=56
pseudoatom ps_K-l, psK or psl, label=20
pseudoatom ps_K-4, psK or ps4, label=7
pseudoatom ps_L-K, psL or psK, label=12
pseudoatom ps_L-M, psL or psM, label=68
pseudoatom ps_L-N, psL or psN, label=27
pseudoatom ps_L-V, psL or psV, label=19
pseudoatom ps_L-W, psL or psW, label=3
pseudoatom ps_L-X, psL or psX, label=3
pseudoatom ps_L-Z, psL or psZ, label=29
pseudoatom ps_L-n, psL or psn, label=29
pseudoatom ps_L-p, psL or psp, label=29
pseudoatom ps_L-r, psL or psr, label=29
pseudoatom ps_L-s, psL or pss, label=27
pseudoatom ps_L-t, psL or pst, label=47
pseudoatom ps_L-u, psL or psu, label=50
pseudoatom ps_L-v, psL or psv, label=52
pseudoatom ps_L-4, psL or ps4, label=5
pseudoatom ps_M-L, psM or psL, label=68
pseudoatom ps_M-N, psM or psN, label=33
pseudoatom ps_M-V, psM or psV, label=7
pseudoatom ps_M-W, psM or psW, label=41
pseudoatom ps_M-Y, psM or psY, label=8
pseudoatom ps_M-Z, psM or psZ, label=24
pseudoatom ps_M-k, psM or psk, label=2
pseudoatom ps_M-o, psM or pso, label=3
pseudoatom ps_M-p, psM or psp, label=41
pseudoatom ps_M-t, psM or pst, label=8
pseudoatom ps_M-v, psM or psv, label=11
pseudoatom ps_M-w, psM or psw, label=46
pseudoatom ps_M-y, psM or psy, label=26
pseudoatom ps_M-4, psM or ps4, label=9
pseudoatom ps_N-A, psN or psA, label=1
pseudoatom ps_N-J, psN or psJ, label=23
pseudoatom ps_N-K, psN or psK, label=56
pseudoatom ps_N-L, psN or psL, label=27
pseudoatom ps_N-M, psN or psM, label=33
pseudoatom ps_N-o, psN or pso, label=31
pseudoatom ps_N-V, psN or psV, label=19
pseudoatom ps_N-W, psN or psW, label=4
pseudoatom ps_N-Y, psN or psY, label=8
pseudoatom ps_N-k, psN or psk, label=62
pseudoatom ps_N-l, psN or psl, label=37
pseudoatom ps_N-4, psN or ps4, label=31
pseudoatom ps_V-L, psV or psL, label=19
pseudoatom ps_V-M, psV or psM, label=7
pseudoatom ps_V-N, psV or psN, label=19
pseudoatom ps_V-W, psV or psW, label=5
pseudoatom ps_W-L, psW or psL, label=3
pseudoatom ps_W-M, psW or psM, label=41
pseudoatom ps_W-N, psW or psN, label=4
pseudoatom ps_W-V, psW or psV, label=5
pseudoatom ps_W-X, psW or psX, label=2
pseudoatom ps_W-Y, psW or psY, label=19
pseudoatom ps_W-Z, psW or psZ, label=45
pseudoatom ps_W-l, psW or psl, label=34
pseudoatom ps_W-o, psW or pso, label=22
pseudoatom ps_W-q, psW or psq, label=1
pseudoatom ps_W-r, psW or psr, label=7
pseudoatom ps_W-t, psW or pst, label=18
pseudoatom ps_W-w, psW or psw, label=30
pseudoatom ps_W-y, psW or psy, label=18
pseudoatom ps_X-L, psX or psL, label=3
pseudoatom ps_X-W, psX or psW, label=2
pseudoatom ps_X-r, psX or psr, label=13
pseudoatom ps_X-t, psX or pst, label=30
pseudoatom ps_X-u, psX or psu, label=5
pseudoatom ps_X-v, psX or psv, label=1
pseudoatom ps_X-n, psX or psn, label=25
pseudoatom ps_Y-H, psY or psH, label=6
pseudoatom ps_Y-M, psY or psM, label=8
pseudoatom ps_Y-N, psY or psN, label=8
pseudoatom ps_Y-W, psY or psW, label=19
pseudoatom ps_Y-o, psY or pso, label=22
pseudoatom ps_Y-q, psY or psq, label=56
pseudoatom ps_Y-z, psY or psz, label=36
pseudoatom ps_Y-l, psY or psl, label=21
pseudoatom ps_Y-m, psY or psm, label=45
pseudoatom ps_Z-L, psZ or psL, label=29
pseudoatom ps_Z-M, psZ or psM, label=24
pseudoatom ps_Z-W, psZ or psW, label=45
pseudoatom ps_Z-o, psZ or pso, label=33
pseudoatom ps_Z-p, psZ or psp, label=15
pseudoatom ps_Z-r, psZ or psr, label=50
pseudoatom ps_Z-s, psZ or pss, label=15
pseudoatom ps_Z-v, psZ or psv, label=2
pseudoatom ps_Z-w, psZ or psw, label=64
pseudoatom ps_Z-y, psZ or psy, label=11
pseudoatom ps_k-J, psk or psJ, label=3
pseudoatom ps_k-M, psk or psM, label=2
pseudoatom ps_k-N, psk or psN, label=62
pseudoatom ps_k-o, psk or pso, label=16
pseudoatom ps_k-w, psk or psw, label=20
pseudoatom ps_k-f, psk or psf, label=2
pseudoatom ps_k-x, psk or psx, label=22
pseudoatom ps_k-4, psk or ps4, label=14
pseudoatom ps_l-J, psl or psJ, label=19
pseudoatom ps_l-K, psl or psK, label=20
pseudoatom ps_l-N, psl or psN, label=37
pseudoatom ps_l-W, psl or psW, label=34
pseudoatom ps_l-Y, psl or psY, label=21
pseudoatom ps_l-o, psl or pso, label=8
pseudoatom ps_l-q, psl or psq, label=53
pseudoatom ps_l-z, psl or psz, label=3
pseudoatom ps_l-m, psl or psm, label=3
pseudoatom ps_m-A, psm or psA, label=9
pseudoatom ps_m-H, psm or psH, label=27
pseudoatom ps_m-J, psm or psJ, label=1
pseudoatom ps_m-Y, psm or psY, label=45
pseudoatom ps_m-l, psm or psl, label=3
pseudoatom ps_m-q, psm or psq, label=17
pseudoatom ps_m-9, psm or ps9, label=7
pseudoatom ps_n-L, psn or psL, label=29
pseudoatom ps_n-X, psn or psX, label=25
pseudoatom ps_n-t, psn or pst, label=18
pseudoatom ps_n-u, psn or psu, label=28
pseudoatom ps_o-M, pso or psM, label=3
pseudoatom ps_o-N, pso or psN, label=31
pseudoatom ps_o-W, pso or psW, label=22
pseudoatom ps_o-Y, pso or psY, label=22
pseudoatom ps_o-Z, pso or psZ, label=33
pseudoatom ps_o-k, pso or psk, label=16
pseudoatom ps_o-l, pso or psl, label=8
pseudoatom ps_o-p, pso or psp, label=5
pseudoatom ps_o-w, pso or psw, label=2
pseudoatom ps_o-x, pso or psx, label=32
pseudoatom ps_p-L, psp or psL, label=29
pseudoatom ps_p-M, psp or psM, label=41
pseudoatom ps_p-Z, psp or psZ, label=15
pseudoatom ps_p-o, psp or pso, label=5
pseudoatom ps_p-t, psp or pst, label=47
pseudoatom ps_p-v, psp or psv, label=68
pseudoatom ps_q-A, psq or psA, label=2
pseudoatom ps_q-H, psq or psH, label=33
pseudoatom ps_q-J, psq or psJ, label=26
pseudoatom ps_q-W, psq or psW, label=1
pseudoatom ps_q-Y, psq or psY, label=56
pseudoatom ps_q-l, psq or psl, label=53
pseudoatom ps_q-m, psq or psm, label=17
pseudoatom ps_q-9, psq or ps9, label=13
pseudoatom ps_q-h, psq or psh, label=17
pseudoatom ps_q-z, psq or psz, label=20
pseudoatom ps_q-4, psq or ps4, label=44
pseudoatom ps_r-L, psr or psL, label=29
pseudoatom ps_r-W, psr or psW, label=7
pseudoatom ps_r-X, psr or psX, label=13
pseudoatom ps_r-Z, psr or psZ, label=50
pseudoatom ps_r-s, psr or pss, label=23
pseudoatom ps_r-t, psr or pst, label=40
pseudoatom ps_s-L, pss or psL, label=27
pseudoatom ps_s-Z, pss or psZ, label=15
pseudoatom ps_s-r, pss or psr, label=23
pseudoatom ps_s-u, pss or psu, label=29
pseudoatom ps_s-v, pss or psv, label=61
pseudoatom ps_t-L, pst or psL, label=47
pseudoatom ps_t-M, pst or psM, label=8
pseudoatom ps_t-W, pst or psW, label=18
pseudoatom ps_t-X, pst or psX, label=30
pseudoatom ps_t-n, pst or psn, label=18
pseudoatom ps_t-p, pst or psp, label=47
pseudoatom ps_t-r, pst or psr, label=40
pseudoatom ps_t-v, pst or psv, label=11
pseudoatom ps_t-w, pst or psw, label=3
pseudoatom ps_u-L, psu or psL, label=50
pseudoatom ps_u-X, psu or psX, label=5
pseudoatom ps_u-n, psu or psn, label=28
pseudoatom ps_u-s, psu or pss, label=29
pseudoatom ps_v-L, psv or psL, label=52
pseudoatom ps_v-M, psv or psM, label=11
pseudoatom ps_v-X, psv or psX, label=1
pseudoatom ps_v-Z, psv or psZ, label=2
pseudoatom ps_v-p, psv or psp, label=68
pseudoatom ps_v-s, psv or pss, label=61
pseudoatom ps_v-t, psv or pst, label=11
pseudoatom ps_w-M, psw or psM, label=46
pseudoatom ps_w-W, psw or psW, label=30
pseudoatom ps_w-Z, psw or psZ, label=64
pseudoatom ps_w-k, psw or psk, label=20
pseudoatom ps_w-o, psw or pso, label=2
pseudoatom ps_w-t, psw or pst, label=3
pseudoatom ps_w-4, psw or ps4, label=5
pseudoatom ps_x-k, psx or psk, label=22
pseudoatom ps_x-o, psx or pso, label=32
pseudoatom ps_y-M, psy or psM, label=26
pseudoatom ps_y-W, psy or psW, label=18
pseudoatom ps_y-Z, psy or psZ, label=11
pseudoatom ps_z-H, psz or psH, label=46
pseudoatom ps_z-J, psz or psJ, label=8
pseudoatom ps_z-Y, psz or psY, label=36
pseudoatom ps_z-l, psz or psl, label=3
pseudoatom ps_z-q, psz or psq, label=20
pseudoatom ps_z-9, psz or ps9, label=2
pseudoatom ps_z-h, psz or psh, label=1
pseudoatom ps_1-a, ps1 or psa, label=52
pseudoatom ps_1-c, ps1 or psc, label=19
pseudoatom ps_1-h, ps1 or psh, label=9
pseudoatom ps_1-2, ps1 or ps2, label=123
pseudoatom ps_1-3, ps1 or ps3, label=50
pseudoatom ps_2-1, ps2 or ps1, label=123
pseudoatom ps_2-a, ps2 or psa, label=37
pseudoatom ps_2-b, ps2 or psb, label=2
pseudoatom ps_2-c, ps2 or psc, label=4
pseudoatom ps_2-3, ps2 or ps3, label=27
pseudoatom ps_3-1, ps3 or ps1, label=50
pseudoatom ps_3-2, ps3 or ps2, label=27
pseudoatom ps_3-5, ps3 or ps5, label=22
pseudoatom ps_3-9, ps3 or ps9, label=21
pseudoatom ps_3-a, ps3 or psa, label=4
pseudoatom ps_3-b, ps3 or psb, label=19
pseudoatom ps_3-c, ps3 or psc, label=62
pseudoatom ps_3-d, ps3 or psd, label=4
pseudoatom ps_3-e, ps3 or pse, label=46
pseudoatom ps_3-g, ps3 or psg, label=20
pseudoatom ps_3-h, ps3 or psh, label=13
pseudoatom ps_3-i, ps3 or psi, label=23
pseudoatom ps_3-4, ps3 or ps4, label=35
pseudoatom ps_4-A, ps4 or psA, label=4
pseudoatom ps_4-H, ps4 or psH, label=33
pseudoatom ps_4-K, ps4 or psK, label=7
pseudoatom ps_4-L, ps4 or psL, label=5
pseudoatom ps_4-M, ps4 or psM, label=9
pseudoatom ps_4-N, ps4 or psN, label=31
pseudoatom ps_4-k, ps4 or psk, label=14
pseudoatom ps_4-q, ps4 or psq, label=44
pseudoatom ps_4-w, ps4 or psw, label=5
pseudoatom ps_4-3, ps4 or ps3, label=35
pseudoatom ps_4-5, ps4 or ps5, label=93
pseudoatom ps_4-6, ps4 or ps6, label=58
pseudoatom ps_4-9, ps4 or ps9, label=78
pseudoatom ps_4-b, ps4 or psb, label=5
pseudoatom ps_4-f, ps4 or psf, label=19
pseudoatom ps_4-g, ps4 or psg, label=4
pseudoatom ps_4-h, ps4 or psh, label=41
pseudoatom ps_5-3, ps5 or ps3, label=22
pseudoatom ps_5-4, ps5 or ps4, label=93
pseudoatom ps_5-6, ps5 or ps6, label=21
pseudoatom ps_5-9, ps5 or ps9, label=8
pseudoatom ps_5-c, ps5 or psc, label=47
pseudoatom ps_5-d, ps5 or psd, label=21
pseudoatom ps_5-f, ps5 or psf, label=66
pseudoatom ps_5-g, ps5 or psg, label=30
pseudoatom ps_5-h, ps5 or psh, label=47
pseudoatom ps_6-A, ps6 or psA, label=6
pseudoatom ps_6-H, ps6 or psH, label=38
pseudoatom ps_6-4, ps6 or ps4, label=58
pseudoatom ps_6-5, ps6 or ps5, label=21
pseudoatom ps_6-9, ps6 or ps9, label=54
pseudoatom ps_6-d, ps6 or psd, label=21
pseudoatom ps_6-h, ps6 or psh, label=1
pseudoatom ps_6-i, ps6 or psi, label=10
pseudoatom ps_9-H, ps9 or psH, label=24
pseudoatom ps_9-m, ps9 or psm, label=7
pseudoatom ps_9-q, ps9 or psq, label=13
pseudoatom ps_9-z, ps9 or psz, label=2
pseudoatom ps_9-3, ps9 or ps3, label=21
pseudoatom ps_9-4, ps9 or ps4, label=78
pseudoatom ps_9-5, ps9 or ps5, label=8
pseudoatom ps_9-6, ps9 or ps6, label=54
pseudoatom ps_9-b, ps9 or psb, label=41
pseudoatom ps_9-c, ps9 or psc, label=11
pseudoatom ps_9-d, ps9 or psd, label=3
pseudoatom ps_9-h, ps9 or psh, label=35
pseudoatom ps_9-i, ps9 or psi, label=74
pseudoatom ps_a-1, psa or ps1, label=52
pseudoatom ps_a-2, psa or ps2, label=37
pseudoatom ps_a-3, psa or ps3, label=4
pseudoatom ps_a-c, psa or psc, label=7
pseudoatom ps_b-2, psb or ps2, label=2
pseudoatom ps_b-3, psb or ps3, label=19
pseudoatom ps_b-4, psb or ps4, label=5
pseudoatom ps_b-9, psb or ps9, label=41
pseudoatom ps_b-c, psb or psc, label=4
pseudoatom ps_b-d, psb or psd, label=16
pseudoatom ps_b-i, psb or psi, label=26
pseudoatom ps_c-1, psc or ps1, label=19
pseudoatom ps_c-2, psc or ps2, label=4
pseudoatom ps_c-3, psc or ps3, label=62
pseudoatom ps_c-5, psc or ps5, label=47
pseudoatom ps_c-9, psc or ps9, label=11
pseudoatom ps_c-a, psc or psa, label=7
pseudoatom ps_c-b, psc or psb, label=4
pseudoatom ps_c-d, psc or psd, label=9
pseudoatom ps_c-f, psc or psf, label=3
pseudoatom ps_c-g, psc or psg, label=30
pseudoatom ps_c-i, psc or psi, label=4
pseudoatom ps_d-A, psd or psA, label=3
pseudoatom ps_d-3, psd or ps3, label=4
pseudoatom ps_d-5, psd or ps5, label=21
pseudoatom ps_d-6, psd or ps6, label=21
pseudoatom ps_d-9, psd or ps9, label=3
pseudoatom ps_d-b, psd or psb, label=16
pseudoatom ps_d-c, psd or psc, label=9
pseudoatom ps_d-g, psd or psg, label=8
pseudoatom ps_d-i, psd or psi, label=2
pseudoatom ps_e-3, pse or ps3, label=46
pseudoatom ps_f-k, psf or psk, label=2
pseudoatom ps_f-4, psf or ps4, label=19
pseudoatom ps_f-5, psf or ps5, label=66
pseudoatom ps_f-c, psf or psc, label=3
pseudoatom ps_f-h, psf or psh, label=6
pseudoatom ps_g-A, psg or psA, label=2
pseudoatom ps_g-3, psg or ps3, label=20
pseudoatom ps_g-4, psg or ps4, label=4
pseudoatom ps_g-5, psg or ps5, label=30
pseudoatom ps_g-c, psg or psc, label=30
pseudoatom ps_g-d, psg or psd, label=8
pseudoatom ps_g-j, psg or psj, label=22
pseudoatom ps_h-H, psh or psH, label=1
pseudoatom ps_h-q, psh or psq, label=17
pseudoatom ps_h-z, psh or psz, label=1
pseudoatom ps_h-1, psh or ps1, label=9
pseudoatom ps_h-3, psh or ps3, label=13
pseudoatom ps_h-4, psh or ps4, label=41
pseudoatom ps_h-5, psh or ps5, label=47
pseudoatom ps_h-6, psh or ps6, label=1
pseudoatom ps_h-9, psh or ps9, label=35
pseudoatom ps_h-f, psh or psf, label=6
pseudoatom ps_h-i, psh or psi, label=17
pseudoatom ps_i-H, psi or psH, label=2
pseudoatom ps_i-3, psi or ps3, label=23
pseudoatom ps_i-6, psi or ps6, label=10
pseudoatom ps_i-9, psi or ps9, label=74
pseudoatom ps_i-b, psi or psb, label=26
pseudoatom ps_i-c, psi or psc, label=4
pseudoatom ps_i-d, psi or psd, label=2
pseudoatom ps_i-h, psi or psh, label=17
pseudoatom ps_j-g, psj or psg, label=22

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

	