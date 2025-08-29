
## convert RNA structure from PDB to SEQ

NUC_To_SEQ 3hjw_D.pdb 3hjw_D.seq

calculate the distance matrix of a given RNA structure + it's sequence

NUC_To_DistMat 3hjw_D.seq 3hjw_D.pdb 3hjw_D.dist_mat 0

dump the key atoms of a given RNA structure + it's sequence

NUC_To_KeyAtom 3hjw_D.seq 3hjw_D.pdb 3hjw_D.key_pdb 7 1

./NUC_To_AllAtom ${id}.seq example/${id}.pdb ${id}.full_pdb
