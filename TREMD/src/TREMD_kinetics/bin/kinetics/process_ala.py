#!/usr/bin/env python
import glob
from ala_kinetics import *
#from kinetics.ala_kinetics import  *
import argparse
from rmsd_st import *
import shutil
from mkdir_p import mkdir_p

def process_ala_dih(trj_l, pdb):
    trj_d = calc_phi_psi_multi(trj_l, pdb)
    phi_psi_l = glob.glob("*_phi_psi.txt")
    process_dih_sign(phi_psi_l, lower_upper_col1=[150, 180],lower_upper_col2=[-180.0, -100.0])
    return trj_d

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xtc_fn", help="xtc trajectory")
    parser.add_argument("-pdb", help="reference pdb")
    parser.add_argument("-rmsd", default=None, help="Set to true to RMSD fit the trajectories.")
    args = parser.parse_args()

    print args

    #for a in args.:
    #    print a

    trj_l = glob.glob(args.xtc_fn)

    print trj_l

    if args.rmsd:
        rmsd_file = calc_rmsd_st(args.pdb, trj_l)

    mkdir_p("rmsd")
    # shutil.move(xtc, "rmsd") for xtc in glob.glob("*_f.xtc")
    #print x
    #shutil.move(x, "rmsd")
    #shutil.move([str(xtc) for xtc in glob.glob("*_rmsd.txt")], "rmsd")
    #shutil.move([str(xtc) for xtc in glob.glob("*_rmsd.txt.png")], "rmsd")

    rmsd_l = glob.glob("*_f.xtc")
    process_ala_dih(rmsd_l, args.pdb)
    plot_all_rama(trj_d)

    mkdir_p("phi_psi")
    #shutil.move("phi_psi*.txt", "phi_psi")
    #shutil.move([str(xtc) for xtc in glob.glob("*phi_psi*.txt")], "phi_psi")

if __name__ == "__main__":
    main()