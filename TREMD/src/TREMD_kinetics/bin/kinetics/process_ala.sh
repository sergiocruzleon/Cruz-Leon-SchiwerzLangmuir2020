#!/bin/bash

xtc_fn=$1
pdb=$2
rmsd=$3

process_ala.py "$xtc_fn" -pdb $pdb -rmsd=$rmsd

mv *_f.xtc rmsd
mv *rmsd.txt* rmsd

mv *psi_phi*txt phi_psi
mv *.png phi_psi 
