# Introduction

This project contains the software developed and required to reproduce the 
results published in the paper:

**Hofmeister Series for Metal-Cation–RNA Interactions: 
The Interplay of Binding Affinity and Exchange Kinetics**

Sergio Cruz-León and Nadine Schwierz*

DOI: https://doi.org/10.1021/acs.langmuir.0c00851
[see the paper online](https://pubs.acs.org/doi/10.1021/acs.langmuir.0c00851).

# What's contained in this project?

- **Extended-PB-cylindrical** - Fortran code to solve the extended Poisson-Boltzmann
equation using the ion-site effective potentials as described 
[here](https://pubs.acs.org/doi/10.1021/acs.langmuir.0c00851). It contains all
ion-binding site effective potentials and a working example.


- **TREMD** - Code to calculate the rates from Temperature Replica Exchange
MD (TREMD) by using the method published [here](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00372).
It includes a working example with all required input files.


# Extended-PB-cylindrical

The code is contained in the file: `pbe_cyl-mesoscale.f`
To compile the fortran code type:

```diff
      gfortran -O3 -o run pbe_cyl-mesoscale.f
```
To run the program just type:
```diff
       ./run
```
The concentrations, potentials as well as the excess of ions as a function of distance is stored in the file `pbe.out2`

### Dependencies

Note that to run the program it requires multiple files (they need to be in the same folder)
```shell
pbe.parm
E1_pmf.dat
E2_pmf.dat
E3_pmf.dat
E4_pmf.dat
E5_pmf.dat
E6_pmf.dat
E7_pmf.dat
```
The `pbe.parm` contains all the input parameters of the sytem while `E*_pmf.dat`
contains the effective potentials that describe the ion-site interactions.

- `pbe.parm` - 
    - q1: Charge of the anion 
    - q2: Charge of the buffer cation
    - q3: Charge of the competitor cation 
    - c1_01: Buffer charge of the buffer cation (mol/Liter)
    - c3_01: Buffer charge of the competitor cation (mol/Liter)
    - sigma (in e/nm^2): Surface charge of the cylinder.
    - radius: of the cylinder
    - a1 a2 a3: Defined for the Fermionic distribution (see eq 12 [here](https://pubs.acs.org/doi/10.1021/acs.langmuir.0c00851))
    - zeta: proportion surface between sites OP, N7 and O6
- `E1_pmf.dat` PMF of the Cl-to RNA
- `E2_pmf.dat` PMF of the buffer cation to RNA backbone
- `E3_pmf.dat` PMF of the titrated cation to RNA backbone
- `E4_pmf.dat` PMF of the buffer cation to RNA nucleobase   (N7)
- `E5_pmf.dat` PMF of the titrated cation to RNA nucleobase (N7)
- `E6_pmf.dat` PMF of the buffer cation to RNA nucleobase   (O6)
- `E7_pmf.dat` PMF of the titrated cation to RNA nucleobase (O6)
- Other files:
    - Epsilon_z.dat: Dielectric constant as a function of distance. 


### Example: 50mM of NaCl and 100mM of KCl

We need to provide the PMF files:

- `cp PMFs/CXY_pmf.dat example/50mM_NaCl-100mM_KCl/E1_pmf.dat` PMF of the Cl-to RNA
- `cp PMFs/NIO-Phos.dat example/50mM_NaCl-100mM_KCl/E2_pmf.dat` PMF of the buffer cation to RNA backbone
- `cp PMFs/KIO-Phos.dat example/50mM_NaCl-100mM_KCl/E3_pmf.dat` PMF of the titrated cation to RNA backbone
- `cp PMFs/NIO-N7.dat example/50mM_NaCl-100mM_KCl/E4_pmf.dat` PMF of the buffer cation to RNA nucleobase   (N7)
- `cp PMFs/KIO-N7.dat example/50mM_NaCl-100mM_KCl/E5_pmf.dat` PMF of the titrated cation to RNA nucleobase (N7)
- `cp PMFs/NIO-O6.dat example/50mM_NaCl-100mM_KCl/E6_pmf.dat` PMF of the buffer cation to RNA nucleobase   (O6)
- `cp PMFs/KIO-O6.dat example/50mM_NaCl-100mM_KCl/E7_pmf.dat` PMF of the titrated cation to RNA nucleobase (O6)


# TREMD

The trajectory has to be continuous in space. For `gromacs` use the command `-demux`.

### Dependencies

This version of the rates code requires to use anaconda/2_4.3.0 (at mpibp:`module load anaconda/2_4.3.0`)


## Indicator function and substates

The input for the code are time series of a discrete set of substates. In this case we use the states [0,1] for unbound and bound respectively.

## Temperature index

A file that contains the index of the time series of the temperatures for all the replicas. 

## Example Li-RNA dinucleotide
The example described in the following corresponds to Li interacting with the OP binding site of the dinucloetide. 
We used 26 replicas between (300-400) K, and each trajectory is 250 ns long. From the time series of the distances, we evaluated the indicator function (eq. 2 of the manuscript).
All the details can be found [here](https://pubs.acs.org/doi/10.1021/acs.langmuir.0c00851).

To run the kinetics analysis type:

```diff
      sh run_kinetics_TREMD.sh 
```

- The output files are in `LIO/kinetics/OP/RESULTS/`, and it contains:
    - Lifetimes: `lifetime-binding.dat `,	`lifetime-unbinding.dat`
    - Populations:  `populations.dat`, 
    - Rates: `rates_binding.dat`, `rates_unbinding.dat`, `s-rates_binding.dat`, `s-rates_unbinding.dat`
    - Binding and unbinding times: `tp_binding.dat`, `tp_unbinding.dat`


# Reference links

- [Group webpage](https://www.biophys.mpg.de/de/schwierz.html)
- [Force-field cations TIP3P-water](https://github.com/bio-phys/Force-fields-for-metal-cations)
