# NLSCCcodes

Custom software focused on performing ground state Natural Linear Scaling Coupled Cluster (NLS-CC) calculations. 
Completed prototype extending the NLS formalism to an excited state framework => (NLS-CIS)
Future work intends to generalize the NLS-CIS software to include NLS-EOM-CC calculations

## Sub-Directories

* **NLSCCenergy** - Contains the custom software to perform NLS-CCSD and NLS-MBPT(2)

* **eomccs** - Contains the custom software to perform NLS-CIS calculations for a UHF reference. 

* **include** - Contains ACESII header files needed in order to successfully compile the NLS executables

* **lib** - Contains compiled static libraries that are needed in order to successfully compile the NLS executables 

* partitionQMs - Custom software to automate fragmentation of a molecule into logical subunits given an input geometry.



## Using the NLS executable(s)

In order to use the NLS executables **~NLSCCenergy/xcalculateCCenergy** or **~eomccs/xeomccs**, a prior ACESII calculation is required. 
These executables harvest integrals, T1/T2 amplitudes, etc from the output files written at the end of a "parent" ACESII calculation that employs 
the Natural Localized Molecular Orbital (NLMO) localization scheme. Thus, the speed of the NLS algorithm depends on the size of the subunit
being studied in the "parent" ZMAT. **Note** SYMMETRY MUST BE TURNED OFF

Example of a "parent" ZMAT employing the NLMO scheme (with the intention of eventually running the NLS executables):
```
H2 geom opt
 O                 -1.65574779   -0.12310111    0.00000000
 O                  1.24513516    0.10232880    0.00041621
 H                 -0.70663626    0.04646664   -0.00004418
 H                 -2.05988123    0.74421549    0.00000108
 H                  1.57490516   -0.38288120   -0.75814479
 H                  1.57188214   -0.37843325    0.76244793

*ACES2
REF=UHF
EXCITE=CIS
BASIS=STO-3G
ESTATE_SYM=10
NON-HF=1
SYMMETRY=NONE
MEM=1GB
```

and an example runscript calling Norbert Flocke's NLMO transformation procedure might look like:

```
#!/bin/bash
#SBATCH --job-name=h2o
#SBATCH --output=zzz.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=99:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --qos=bartlett

echo $SLURM_NODELIST
echo $SLURM_JOBID

export TESTROOT=$SLURM_SUBMIT_DIR
export WORKDIR=$TESTROOT/temp

cd $TESTROOT
mkdir -p $WORKDIR
cd $WORKDIR
cp $TESTROOT/ZMAT .
#cp $TESTROOT/QMcenter .
cp $TESTROOT/GENBAS .
#ln -s ~/GENBAS

echo "`date`" > out

xjoda >> out
xvmol >> out
xvmol2ja >> out
xvscf >> out
/home/auybv/NLSCC4/nlorb/xnlorb >> nlscc.out
xvtran >> out
xintprc >> out
xvee >> out
```
where "xnlorb" is the executable performing the NLMO localization and originates in the ACESII source code under the "nlorb" directory

Once this "parent" ACESII calculation has been performed, an additional file called "QMcenter" must be created that separates the constituent atoms
of the ZMAT into QM1 or QM2 regions. An example of "QMcenter" for our water dimer example looks like (nothing that the # must be deleted in a real calc):

```
3     # Number of atoms in QM1

1     # Indices of the atoms in QM1, ordered according to ZMAT input file
3     # So water molecule 1 - having atom numbers 1,3,4 - are inside QM1
4

3     # Number of atoms in QM2

2     # Indices of the atoms in QM2
5
6
```

The NLS software references the atom labels in this file, in conjunction with output from Norbert Flocke's NLMO procedure into "nbocenters" to determine which atoms "own" which
localized molecular orbital. 

An example of output from the NLMO procedure looks like: 

```
   1     1
   2     2
   3     1
   4     1
   5     2
   6     2
   7     1
   7     4
   8     1
   8     3
   9     2
   9     5
  10     2
  10     6
  11     1
  11     4
  12     1
  12     3
  13     2
  13     5
  14     2
  14     6

```

where the first column refers to the NLMO #, and the second column refers to the atom which "owns" the NLMO. **NOTE** It is possible for NLMO's 
to be share across QM regions; the NLS adjusts for this eventuality. 


In summary, once the ACESII "parent" calculation utilizing the NLMO procedure completes and you are inside this directory, create "QMcenter" based on
your chemical intuition and simply run the executable.





