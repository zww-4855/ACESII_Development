      Subroutine Analyze(Dens_a,Dens_b,Work,Nocc_a,Nocc_b,Maxcor,
     &                   Nbfns,writeFile)

      Implicit Double Precision(A-H,O-Z)

      Dimension Dens_a(Nbfns*Nbfns,Nocc_a)
      Dimension Dens_b(Nbfns*Nbfns,Nocc_b)
      dimension propInt(Nbfns*Nbfns)
      Dimension Work(Maxcor)
      Dimension checka(Nbfns*Nbfns)
      Dimension checkb(Nbfns*Nbfns)
      double precision spens ! spin density from FC calculation
      Data Done /1.0D0/
      COMMON /FILES/ LUOUT,MOINTS
      dimension diag(Nbfns*Nbfns),diagb(Nbfns*Nbfns)
      dimension mata(Nbfns,Nbfns)
      real(kind=8),allocatable:: spdena(:),spdenb(:),spdenALL(:)
      integer NATOM,closestatus,openstatus,writeFile

#include "machsp.com"
C       SUBROUTINE takes the decomposed alpha&beta density matrix
C       and, alongside the property integrals, calls compp.F. The emphasis
C       of the property integrals focused on the evaluation of the 
C       electronic density at the nucleus.
C       compp.F then calculates the spin density, in this case per
C       alpha&beta MO, and returns the spin density to SUBROUTINE.
C       SUBROUTINE then writes this to output file.
 
      IONE=1
c       Sort spin density according to atom center
      CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOM)
      print*, 'number of atoms is: ', NATOM

        Nbfns2 = Nbfns*Nbfns 

      Idens_a = I000
      Idens_b = Idens_a + Nbfns2 
      Iend    = Idens_b + Nbfns2 

      allocate(spdena(Nocc_a),spdenb(Nocc_b),spdenall(Nocc_a))
      
      if (writeFile.eq.0) then
         open(unit=159, file='hfspinDen.txt',status='unknown',!old',
     &  action='write')!,position='rewind')
      else
           open(unit=159, file='CCspinDen.txt',status='unknown',!'old',
     &  action='write')!,position='rewind')
      endif
C       Do loop runs for every atom in ZMAT
C       Work(Idens_a/b) is intended to hold a running summation
C       ,using LAPACK daxpy routine,
C       of the MO density matrices for later comparison against
C       the total alpha/beta Density matrix. 
      do n=1,NATOM 
              Call Dzero(Work(Idens_a),Nbfns2)
              Call Dzero(Work(Idens_b),Nbfns2) 
C             2 Do loops for alpha&beta quantities
              Do Iocc = 1, Nocc_a
C           SEEKLB offsets JOBARC to the appropriate density
                 if (n.eq.1) then
                     IRWND=0
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                     IF(IERR.NE.0)RETURN
                 else
                     IRWND=0
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                     IF(IERR.NE.0)RETURN
                     IRWND=1
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                 endif
                 spens=0.0d0

                Call Daxpy(Nbfns2,Done,Dens_a(1,Iocc),1,Work(Idens_a),1)
                call COMPPR(spens, Dens_a(1,Iocc), propInt, Nbfns, 
     &                                                  .False.)
                spdena(Iocc)=spens
#ifdef _DEBUG_LVL0
                print*, 'spin Density for alpha e: ', Iocc, spens
#endif
              Enddo 

              Do Iocc = 1, Nocc_b
                spens=0.0d0
                 if (n.eq.1) then
                     IRWND=0
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                     IF(IERR.NE.0)RETURN
                 else
                     IRWND=0
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                     IF(IERR.NE.0)RETURN
                     IRWND=1
                     CALL SEEKLB('   DEN  ',IERR,IRWND)
                 endif
                Call Daxpy(Nbfns2,Done,Dens_b(1,Iocc),1,Work(Idens_b),1)
                call COMPPR(spens, Dens_b(1,Iocc), propInt, Nbfns,
     &                                                   .False.)
                spdenb(Iocc)=spens
#ifdef _DEBUG_LVL0
                print*, 'spin Density for beta e: ', Iocc, spens
#endif
              Enddo 
C       Determines the spin density, defined as the difference
C       between alpha and beta density, in this case at the nucleus
              do i=1,Nocc_a
                if (i.eq.Nocc_a) then
                  spdenall(i)=spdena(i)
                  exit
                endif
                spdenall(i)=spdena(i)-spdenb(i)
              enddo
C        Write to file
         do i=1,Nocc_a
            write(159,"(A14,I1,A9,F14.7,3X,F14.7,3X,F14.7)")
     &      'a/b MO ',i,' density:' , spdena(i),spdenb(i),spdenall(i)
            if (i.eq.Nocc_a) then
              write(159,*) '  '
              exit
            endif 
         enddo
         iatom=iatom+1
      enddo

      close(159,iostat=closestatus)
      IRWND=0 ! necessary to reset this for subsequent calcs
      deallocate(spdena,spdenb,spdenall)

C                       *********DEBUG**********
C       Ensure the density matrices per MO summed equals the total
C       Density matrix. Also verify diagonalization of density matrix
C       yields eigenvector matrix whose trace is number of electrons
#ifdef _DEBUG_LVL0

      if (writeFile.eq.0) then
         Call Getrec(20,"JOBARC","SCFDENSA",Nbfns2*Iintfp,checka(1))
         Call Getrec(20,"JOBARC","SCFDENSB",Nbfns2*Iintfp,checkb(1))
      else
         Call Getrec(20,"JOBARC","RELDENSA",Nbfns2*Iintfp,checka(1))
         Call Getrec(20,"JOBARC","RELDENSB",Nbfns2*Iintfp,checkb(1)) 
      endif
      Call Daxpy(Nbfns2,-Done,Work(Idens_a),1,checka(1),1)
      print*, 'Difference corr between JOBarc/Calc alpha Density==0',
     &          checka(1),Work(Idens_a)
      call output(checka,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)




      counter=1
      total=0.0d0
      do i=1,Nbfns
        do j=1,Nbfns
          mata(i,j)=checka(counter)
          counter=counter+1
          if (i.eq.j) total=total+mata(i,i)
        enddo
      enddo
      print*,'verify is: ', total
      
      call eig(mata,diag,1,Nbfns,1)


      totalb=0.0d0
      total=0.0d0
      totalmat=0.0d0
      do i=1,Nbfns !13!Nbfns
        total=total+checka(Nbfns*(i-1)+i)
        totalb=totalb+checkb(Nbfns*(i-1)+i)
        totalmat=totalmat+mata(i,i)
      enddo

      print*, 'trace alphaDens matrix',checka(1), total,totalb
      print*,'verify trace', Work(Idens_a),totalmat



      call eig(checka,diag,1,Nbfns,1)
      call eig(checkb,diagb,1,Nbfns,1)
      print*,'diagonal element a',checka(1),diag(1)
      print*,'diagonal element b',checkb(1),diagb(1)

      total=0.0d0
      totalb=0.0d0
      do i=1,Nbfns*Nbfns
        total=total+diag(i)!(Nbfns*(i-1)+i)
      enddo
      print*, 'trace alpha matrix', total
      call output(checka(1),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

#endif

      Return
      End 


 
