# ACESII_Development
Additions and/or modifications to the ACESII computational chemistry software, focused on Coupled Cluster Theory. Each sub-directory is
meant to be independent.



## Sub-Directories:

* **NLSCCcodes** - Contains the relevant software needed to perform ground and excited state Natural Linear Scaling Coupled Cluster 
                   (NLS-CC) calculations, as well as custom software designed to automate the fragmentation process
               
* **driver** - Utilizes the ACESII ZMAT and GENBAS files, as well as the ACESII executables xjoda, xvmol, xvmol2ja to harvest one and two 
               electron quantities. Custom software that performs Restricted Hartree-Fock theory(HF), Configuration Interaction Singles (CIS), 
               Time-Dependent Hartree Fock(TDHF), Coupled Perturbed Hartree Fock (CPHF), and Many Body Perturbation Theory (MBPT2) for a water molecule 
               =>>>>>>> MUST GENERALIZE TO ANY GEOMETRY & TO UHF REFERENCE IN THE FUTURE
               
* **props** - Extension of existing ACESII source code containing property calculation subroutines. Contains custom additions to 
              decompose the isotropic hyperfine coupling constant (Fermi Contact) into contributions from natural orbitals. 
              General software that is applicable to any electronic structure calculation in ACESII.
              
              
              
## Dependencies

  * intel/2019 compilers
  * GNU gmake 3.82
