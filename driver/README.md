# driver
Custom electronic structure software that interfaces to ACESII to obtain relevant one and two center integrals. Performs RHF, CIS, MBPT(2), CPHF, and
TDHF calculations, given an ACESII ZMAT input file and associated GENBAS basis set file. Prior to running **xdriver**, the following ACESII executables 
must be called to generate the required integrals:

* xjoda
* xvmol
* xvmol2ja

This can then be followed by **~/driver/xdriver**
Thusfar, I have only been able to clean and generalize the RHF code.

## UNDER CONSTRUCTION TO DO:
  * finalize the generalization of the remaining programs
  * adapt all methods to an UHF reference
