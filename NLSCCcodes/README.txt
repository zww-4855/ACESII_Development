Directory containing all work dedicated towards the Natural Linear Scaling Coupled Cluster theory. 

Thus far there are 2 directories:

1) NLSCCenergy:
      Contains all subroutines used to calculate NLSCCSD energy. Inspiration was taken from multiple existing ACES2 routines in */libr including cmpeng.F, tener.F, ftau.F.
      The prerequistites for running the executable are a file specifying the QM1 region - called QMcenter - used together with standard ACES2 output from the NLMO module
      called "nbocenters". The formatting for QMcenter is:
      
      # of atoms in QM1
      space
      list of atom index within QM1 wrt to ZMAT input file
      space
      # of atoms in QM2
      space
      list of atom index within QM2 wrt to ZMAT input file
      
      The formatting for nbocenters is 2 columns, where the first column in NLMO number, and the second column is the atom index it belongs to. Note, a NLMO can belong
      to several atoms. Once these two files are appropriately read, we load up the appropriate T1,fia,T2, and 2 e- integrals. The formatting of these data files was 
      obtained using ACES2v1.0.0 manual as a reference. They are stored like <AB||IJ>. This software has been tested and works using a RHF reference. UHF reference
      loads appropriate data structures, but I have not added NLSCC capabilities on top of this. 
      
2) PartitionQMs:
    This directory contains relevant subroutines intended to take a list of QM1 index files, and a primary txt file containing their names, as well as appropriate 
    cutoff values defining the borders of QM2. The input is a .xyz file containing the entire molecule of interest. The output, based on the QM1 index files, is
    a new set of .xyz files which include only QM1 and QM2 - nothing else of the original molecule remains. Whichever bond is cut, a terminal hydrogen atom is added
    in the direction of that bond. Very crude program. It works for linear molecules, but has difficulty cutting rings and double bonds where its necessary to add 2
    hydrogen atom caps. Will need to add an additional layer of sophistication to this module. 
    
**** Work in Progress *****

Add EOM-NLSCCSD capabilities to calculate spectra 
