










      Subroutine Built_scf_orbdens(Work,Maxcor,Nbfns,Naobfns,Nirrep,
     &                             Occ,Vrt,Nocc,Nvrt,Iuhf)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Work(Maxcor)
      Integer Occ(8,2),Vrt(8,2),Nocc(2),Nvrt(2)
      Character*8 Vecs(2) 

      Data Vecs /"SCFEVCA0","SCFEVCB0"/
      Data Done,Dzero /1.0D0,0.0D0/



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



     
      Nbfns2    = Nbfns*Nbfns 

      I000      = 1
      Iscfvec_a = I000
      Iscfden_a = Iscfvec_a + Nbfns2
      Iscfvec_b = Iscfden_a + Nbfns2*Nocc(1)
      If (Iuhf .Ne. 0) Then 
          Iscfden_b = Iscfvec_b + Nbfns2
          Iend      = Iscfden_b + Nbfns2*Nocc(2)
          If (Iend .GT. Maxcor) Call Insmem("built_scforb_dens",
     &                                       Iend,Maxcor)
      Else
         Write(6,"(a)") " The spin-densities are zero for closed"
     &                  " molecules" 
         Call Errex
      Endif 

C Read in alpha and beta SCF vectors from JOBARC. The SCFEVECA0
C and SCFEVECB0 are (AO,MO) quantities. 

      Call Getrec(20,"JOBARC",Vecs(1),Nbfns2*Iintfp,Work(Iscfvec_a))
      Call Getrec(20,"JOBARC",Vecs(2),Nbfns2*Iintfp,Work(Iscfvec_b))

      Write(6,"(a)") " The alpha eigenvectors"
      call output(Work(Iscfvec_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
      Write(6,"(a)") " The beta eigenvectors"
      call output(Work(Iscfvec_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

C Built the density matrix from eigenvectors and compare with what 
C is on JOBRAC records, SCFDENSA. I do this only for alpha density.

      Call Xgemm("N","T",Nbfns,Nbfns,Nocc,Done,Work(Iscfvec_a),
     &            Nbfns,Work(Iscfvec_a),Nbfns,Dzero,Work(Iscfden_a),
     &            Nbfns)
      Write(6,"(a)") " The alpha density"
      call output(Work(Iscfden_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

C check with the density matrix from JOBARC 

      Call Getrec(20,"JOBARC","SCFDENSA",Nbfns2*Iintfp,Work(Iscfden_b))
      print*, 'first element of scfdensa is; ', Work(Iscfden_b)
      Call Daxpy(Nbfns2,-Done,Work(Iscfden_b),1,Work(Iscfden_a),1)
      Write(6,"(a)") "(Dens_a(prop)-Densa(jarc)),should be zero!!"
      call output(Work(Iscfden_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

C Now proceed to built alpha and beta orbital density matrices
      
      Call Orbdens(Work(Iscfvec_a),Work(Iscfden_a),Nocc(1),Nbfns,1)
      Call Orbdens(Work(Iscfvec_b),Work(Iscfden_b),Nocc(2),Nbfns,2)

C First sum the alpha and beta orbital density matrices to verify
C that they match with JOBRAC records, SCFDENSA and SCFDENSB. 
C Then proceed to do the spin-densities for each orbital. Make
C sure to pass AO integrals. 

      Memleft = Maxcor - Iend 
      Call Analyze(Work(Iscfden_a),Work(Iscfden_b),Work(Iend),Nocc(1),
     &             Nocc(2),Memleft,Nbfns,0)
     
      Return
      End 
