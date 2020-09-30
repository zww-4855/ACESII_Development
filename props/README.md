# props
Contains modified ACESII source code devoted to electronic property calculations. 


The primary modification was made to **prop.F** via calling

* built_scf_orbdens.F
* built_CC_orbdens.F

These subroutines further call newly created subroutines **orbdens.F** (orbdensCC.F), **mkden.F**, and **analyze.F** to perform the relevant property analysis.

The purpose of these auxillary softwares is intended to expand current ACESII capabilities involving the calculation of isotropic hyperfine coupling. Essentially, these files take the SCF/CC alpha and beta density matrix, and decomposes it according to occupied Molecular Orbital (MO) contributions. After this process is complete, the spin density (and consequently the Fermi Contact) is calculated on an individual MO basis. The software then writes these values, alongside the spin density, to an external output file(s) **CCspinDen.txt** and **hfspinDen.txt**.

To facilitate this procedure, the existing ACESII data structure for the Unrestricted Hartree Fock and Coupled Cluster density is diagonalized - the objective being to use the corresponding eigenvectors to build natural orbital density matrices for each alpha/beta MO. For example, the UHF alpha density matrix for MO 1 could be described by the outer product between the first column of the eigenvector matrix and its transpose. 

In order to use alongside existing ACESII binaries, modify the PATH environment variable to include the location of this directory, for example:

```
export PATH=/home/z.windom/ACESII/props:/apps/shared/bartlett/ACESII-2.12.0/bin/:$PATH
```
 Then ACESII will reference this "xprops" executable upon completing the relevant electronic structure calculation.  


#TO DO: 
* clean up unnecessary printing
* combine or clean redundancy of built_XX_orbdens.F and orbdensXX.F
* results for first-order property calculations besides Fermi Contact are unverified, but should still be legit => would be interesting to verify
