      Subroutine Built_CC_orbdens(Work,Maxcor,
     &                             Nbfns,Naobfns,Nirrep,
     &                             Occ,Vrt,Nocc,Nvrt,Iuhf)

      Implicit Double Precision(A-H,O-Z)
     
      Dimension Work(Maxcor),eigval_a(Nbfns*Nbfns),eigval_b(Nbfns*Nbfns)
      Integer Occ(8,2),Vrt(8,2),Nocc(2),Nvrt(2)
      Character*8 Vecs(2) 
      Dimension scftransA(Nbfns*Nbfns),scftransB(Nbfns*Nbfns)
      Data Vecs /"SCFEVCA0","SCFEVCB0"/
      Data Done,Dzero /1.0D0,0.0D0/

#include "machsp.com"
     
      Nbfns2    = Nbfns*Nbfns 
        ! relDens is for RELDENSA/B correlated density matrix
        ! den_a/b is for the density matrices associated with 
        !               alpha/beta MO 1,2,..., etc
        ! eig_a/b is meant to store eigenvectors of relDens

      I000      = 1
      IscfrelDens_a = I000
      Iscfden_a = IscfrelDens_a + Nbfns2
      IscfrelDens_b = Iscfden_a + Nbfns2*Nocc(1)
      If (Iuhf .Ne. 0) Then 
          Iscfden_b = IscfrelDens_b + Nbfns2
         
          Iscfeig_a = Iscfden_b + Nbfns2*Nocc(2)
          Iscfeig_b = Iscfeig_a + Nbfns2

          Iend      = Iscfeig_b + Nbfns2
          If (Iend .GT. Maxcor) Call Insmem("built_scforb_dens",
     &                                       Iend,Maxcor)
      Else
         Write(6,"(a)") " The spin-densities are zero for closed"
     &                  " molecules" 
         Call Errex
      Endif 

C Read in alpha and beta correlated density matrices from JOBARC. The SCFEVECA0
C and SCFEVECB0 are (AO,MO) quantities. 


        CALL GETREC(20,'JOBARC','RELDENSA',IINTFP*Nbfns*Nbfns,!NCOMP*NCOMP,
     &              Work(IscfrelDens_a))
         CALL GETREC(20,'JOBARC','RELDENSB',IINTFP*Nbfns*Nbfns,
     &               Work(IscfrelDens_b))
c      Call Getrec(20,"JOBARC",Vecs(1),Nbfns2*Iintfp,Work(IscfrelDens_a))
c      Call Getrec(20,"JOBARC",Vecs(2),Nbfns2*Iintfp,Work(IscfrelDens_b))

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The correlated alpha Density Matrix"
      call output(Work(IscfrelDens_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
      Write(6,"(a)") " The correlated beta Density Matrix"
      call output(Work(IscfrelDens_b),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif 

      call dcopy(nbfns*nbfns,Work(IscfrelDens_a),1,eigval_a,1)
      call dcopy(nbfns*nbfns,Work(IscfrelDens_b),1,eigval_b,1)

      call eig(eigval_a,Work(Iscfeig_a),1,Nbfns,-1)
      Write(6,"(a)") " The correlated alpha eigenvalues"
      call output(eigval_a,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

      call eig(eigval_b,Work(Iscfeig_b),1,Nbfns,-1)
      Write(6,"(a)") " The correlated beta eigenvalues"
      call output(eigval_b,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
      ! check the trace to verify it is the total number of alpha/beta
      ! electrons
      total=0.0d0
      totalb=0.0d0
      do i=1,Nbfns
        total=total+Work(IscfrelDens_a+Nbfns*(i-1)+(i-1))
        totalb=totalb+Work(IscfrelDens_b+Nbfns*(i-1)+(i-1))
        total1=0.0d0
        total2=0.0d0
        do j=1,Nbfns*Nbfns
        total1=total1+eigval_a(j)
        total2=total2+eigval_b(j)
        enddo
      enddo
      print*,'correlated trace is: ', total,totalb
      print*,'corrected** correlated trace is: ', total1,total2
      print*,'eiga cc chekc:', Work(Iscfeig_b),Work(Iscfeig_a)
c      call dcopy(nbfns*nbfns,Iscfeig_a,1,Work(IscfrelDens_a),1)
c      call dcopy(nbfns*nbfns,Iscfeig_b,1,Work(IscfrelDens_b),1)

            Call Getrec(20,"JOBARC","SCFEVCA0",Nbfns*Nbfns*Iintfp,
     &                  scftransA)

            Call Getrec(20,"JOBARC","SCFEVCB0",Nbfns*Nbfns*Iintfp,
     &                  scftransB)

            print*,'transform',scftransA(1),scftransB(1)



        print*, 'before matmul inlocal', Work(Iscfeig_a),Work(Iscfeig_b)
        call xgemm('N','N',Nbfns,Nbfns,Nbfns,1.0D0,scftransA,
     &          Nbfns,Work(Iscfeig_a),Nbfns,0.0D0,Work(Iscfeig_a),Nbfns)
        call xgemm('N','N',Nbfns,Nbfns,Nbfns,1.0D0,scftransB,
     &          Nbfns,Work(Iscfeig_b),Nbfns,0.0D0,Work(Iscfeig_b),Nbfns)
        print*, 'after matmul inlocal:', Work(Iscfeig_a),Work(Iscfeig_b)
      Call OrbdensCC(eigval_a,Work(Iscfeig_a),Work(Iscfden_a)
     &                          ,Nocc(1),Nbfns,1)
      Call OrbdensCC(eigval_b,Work(Iscfeig_b),Work(Iscfden_b)
     &                          ,Nocc(2),Nbfns,2)


      Memleft = Maxcor - Iend
      Call Analyze(Work(Iscfden_a),Work(Iscfden_b),Work(Iend),Nocc(1),
     &             Nocc(2),Memleft,Nbfns,1)
c#ifdef _DEBUG_LVL0
C Built the density matrix from eigenvectors and compare with what 
C is on JOBRAC records, SCFDENSA. I do this only for alpha density.

c      Call Xgemm("N","T",Nbfns,Nbfns,Nocc,Done,Work(IscfrelDens_a),
c     &            Nbfns,Work(IscfrelDens_a),Nbfns,Dzero,Work(Iscfden_a),
c     &            Nbfns)
c      Write(6,"(a)") " The alpha density"
c      call output(Work(Iscfden_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
c
cC check with the density matrix from JOBARC 
c:q

c      Call Getrec(20,"JOBARC","SCFDENSA",Nbfns2*Iintfp,Work(Iscfden_b))
c      print*, 'first element of scfdensa is; ', Work(Iscfden_b)
c      Call Daxpy(Nbfns2,-Done,Work(Iscfden_b),1,Work(Iscfden_a),1)
c      Write(6,"(a)") "(Dens_a(prop)-Densa(jarc)),should be zero!!"
c      call output(Work(Iscfden_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
c#endif
c
cC Now proceed to built alpha and beta orbital density matrices
c      
c      Call Orbdens(Work(IscfrelDens_a),Work(Iscfden_a),Nocc(1),Nbfns,1)
c      Call Orbdens(Work(IscfrelDens_b),Work(Iscfden_b),Nocc(2),Nbfns,2)
c
cC First sum the alpha and beta orbital density matrices to verify
cC that they match with JOBRAC records, SCFDENSA and SCFDENSB. 
cC Then proceed to do the spin-densities for each orbital. Make
cC sure to pass AO integrals. 
c
c      Memleft = Maxcor - Iend 
c      Call Analyze(Work(Iscfden_a),Work(Iscfden_b),Work(Iend),Nocc(1),
c     &             Nocc(2),Memleft,Nbfns)
     
      Return
      End 
