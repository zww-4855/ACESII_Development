Contains the relevant folder for ACESII property calculations. The primary modification was made to prop.F via calling

1) built_scf_orbdens.F
2) built_CC_orbdens.F

These routines further call newly created routines orbdens.F (orbdensCC.F), mkden.F, and analyze.F to perform the relevant property analysis.

The purpose of these auxillary softwares is intended to expand current ACESII capabilities involving the calculation of isotropic hyperfine coupling. Essentially, these files take the SCF/CC alpha and beta density matrix, and decomposes it according to occupied Molecular Orbital (MO) contributions. After this process is complete, the spin density (and consequently the Fermi Contact) is calculated on an individual MO basis. The software then writes these values, alongside the spin density, to an external output file(s) CCspinDen.txt and hfspinDen.txt.

To facilitate this procedure, the existing data structure for the Unrestricted Hartree Fock and Coupled Cluster density is diagonalized - the objective being to use the corresponding eigenvectors to build natural orbital density matrices for each alpha/beta MO. For example, the UHF alpha density matrix for MO 1 would take the first column of eigenvector matrix and matrix multiply it with its transpose. 

In order to use alongside existing ACESII binaries, modify PATH to include the location of this directory. 



TO DO: 
* clean up unnecessary printing
* combine or clean redundancy of built_XX_orbdens.F and orbdensXX.F
