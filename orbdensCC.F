      Subroutine OrbdensCC(Evals,Evecs,Dens,Nocc,Nbfns,Ispin)

      Implicit Double Precision(A-H,O-Z)

      Dimension Evecs(Nbfns,Nbfns), Dens(Nbfns*Nbfns,Nocc)
      Dimension Evals(Nbfns*Nbfns)
      dimension scftransA(Nbfns,Nbfns),scftransB(Nbfns,Nbfns)
      dimension tempDens(Nbfns,Nbfns),B(Nbfns,Nbfns)
      Data Done,Dnull /1.0D0,0.0D0/


#include "machsp.com"


      Call Dzero(Dens,Nbfns*Nbfns*Nocc)
      print*,'eiga cccheck',Evecs(1,1)

      Do Iocc = 1, Nocc
         call Dzero(B,Nbfns*Nbfns)
         Call Mkden(Evecs(1,Iocc),Dens(1,Iocc),Nbfns)

         print*,'b4 dens:',Dens(1,Iocc)
         Dens(1,Iocc)=Dens(1,Iocc)*Evals(Nbfns*(Iocc-1)+Iocc)
        do j=1,Nbfns*Nbfns
         Dens(j,Iocc)=Dens(j,Iocc)*Evals(Nbfns*(Iocc-1)+Iocc)
        enddo
 
        print*,'after dens:',Dens(1,Iocc),Evals(Nbfns*(Iocc-1)+Iocc)
        
         CALL DCOPY(Nbfns*Nbfns, Dens(1,Iocc),1,B,1)
         print*, 'ispin is: ', Ispin
c         if (Ispin.eq.1) then
c            Call Getrec(20,"JOBARC","SCFEVCA0",Nbfns*Nbfns*Iintfp,
c     &                  scftransA(1,1))
c            print*,'transform',scftransA(1,1),Dens(1,Iocc),B(1,1)
c            call xgemm('N','N',Nbfns,Nbfns,Nbfns,1.0D0,scftransA(1,1),
c     &               Nbfns,B,Nbfns,0.0D0,Dens(1,Iocc),Nbfns)
c        else
c            Call Getrec(20,"JOBARC","SCFEVCB0",Nbfns*Nbfns*Iintfp,
c     &                  scftransB)
c            call xgemm('N','N',Nbfns,Nbfns,Nbfns,1.0D0,B,Nbfns,
c     &                 scftransB,Nbfns,0.0D0,Dens(1,Iocc),Nbfns)
c        endif
         print*,'after matmul:',Dens(1,Iocc)
c        print*,'tempDens', tempDens(1,1)
#ifdef _DEBUG_LVL0
        Write(6,*) 
        If (Ispin .EQ. 1) Then
        write(6,"(a)") "The alpha orbital density matrices"
        call output(Dens(1,Iocc),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
        else 
        write(6,"(a)") "The beta orbital density matrices"
        call output(Dens(1,Iocc),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
        endif 
#endif

      Enddo 

      Return
      End 


 
