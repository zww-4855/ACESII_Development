










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



 
c       Sort spin density according to atom center
      CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOM)
c      CALL GETREC(20,'JOBARC','MAP2ZMAT',2,NATOM)
c      CALL GETREC(20,'JOBARC','ATOMCHRG',NUMATOM,IATNUM)
      print*, 'number of atoms is: ', NATOM

        Nbfns2 = Nbfns*Nbfns 

      Idens_a = I000
      Idens_b = Idens_a + Nbfns2 
      Iend    = Idens_b + Nbfns2 

c      Call Dzero(Work(Idens_a),Nbfns2)
c      Call Dzero(Work(Idens_b),Nbfns2)
c      CALL SEEKLB('   DEN  ',IERR,IRWND)

      allocate(spdena(Nocc_a),spdenb(Nocc_b),spdenall(Nocc_a))
      
      if (writeFile.eq.0) then
         open(unit=159, file='hfspinDen.txt',status='unknown',!old',
     & ! form='FORMATTED',iostat=openstatus,
     &  action='write')!,position='rewind')
      else
           open(unit=159, file='CCspinDen.txt',status='unknown',!'old',
     & ! form='FORMATTED',iostat=openstatus,
     &  action='write')!,position='rewind')
      endif

      do n=1,2 !number of atoms in ZMAT!!!
              Call Dzero(Work(Idens_a),Nbfns2)
              Call Dzero(Work(Idens_b),Nbfns2) 
c               CALL SEEKLB('   DEN  ',IERR,IRWND)
c                IF(IERR.NE.0)RETURN
c              IRWND=1
              Do Iocc = 1, Nocc_a
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
c                CALL SEEKLB('   DEN  ',IERR,IRWND)
c                IF(IERR.NE.0)RETURN
c                 if (n.eq.1) then
c                     IRWND=0
c                 else
c                     IRWND=1
c                 endif

                Call Daxpy(Nbfns2,Done,Dens_a(1,Iocc),1,Work(Idens_a),1)
             call COMPPR(spens, Dens_a(1,Iocc), propInt, Nbfns, .False.)
                spdena(Iocc)=spens
                print*, 'spin Density for alpha e: ', Iocc, spens
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
c                IRWND=0
c                CALL SEEKLB('   DEN  ',IERR,IRWND)
                Call Daxpy(Nbfns2,Done,Dens_b(1,Iocc),1,Work(Idens_b),1)
             call COMPPR(spens, Dens_b(1,Iocc), propInt, Nbfns, .False.)
                spdenb(Iocc)=spens
                print*, 'spin Density for beta e: ', Iocc, spens
              Enddo 

              do i=1,Nocc_a
                if (i.eq.Nocc_a) then
                  spdenall(i)=spdena(i)
                  exit
                endif
                spdenall(i)=spdena(i)-spdenb(i)
              enddo
      
c      if (writeFile.eq.0) then
c            open(unit=159, file='hfspinDen.txt',status='unknown',!old',
c     & ! form='FORMATTED',iostat=openstatus,
c     &  action='write')!,position='rewind')
c      else
c           open(unit=159, file='CCspinDen.txt',status='unknown',!'old',
c     & ! form='FORMATTED',iostat=openstatus,
c     &  action='write')!,position='rewind')
c      endif
      do i=1,Nocc_a
        write(159,"(A14,I1,A9,F14.7,3X,F14.7,3X,F14.7)")
     &      'a/b MO ',i,' density:' , spdena(i),spdenb(i),spdenall(i)
        if (i.eq.Nocc_a) then
          write(159,*) '  '
          exit
        endif 
      enddo
c      close(159,iostat=closestatus)
c      IRWND=1
c      CALL SEEKLB('   DEN  ',IERR,1)
      iatom=iatom+1
      enddo
      close(159,iostat=closestatus)
      IRWND=0
      deallocate(spdena,spdenb,spdenall)
c        write(6,"(a)") "The alpha density matrix"
c        call output(Work(Idens_a),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)

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
        !totalb=totalb+checkb(Nbfns*(i-1)+i)
      enddo
      print*, 'trace alpha matrix', total
      call output(checka(1),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
      !print*, 'trace beta matrix', totalb
c      Call Daxpy(Nbfns2,Done,checkb(Idens_a),1,checka(1),1)!total DM



c      print*, 'checka(1)',checka(1),checkb(1)
c
c      Call Daxpy(Nbfns2,-Done,Work(Idens_a),1,checka(1),1)
c      print*, 'Difference between JOBarc/Calc alpha Density==0',
c     &          checka(1),Work(Idens_a)
c      call output(checka,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
c        write(6,"(a)") "The beta density matrix"
c        call output(Work(Idens_b),1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
c
c
c      Call Daxpy(Nbfns2,-Done,Work(Idens_b),1,checkb(1),1)
c      print*, 'Difference between JOBarc/Calc beta Density==0',
c     &          checkb(1),Work(Idens_b)
c
c
c      CALL SEEKLB('   DEN  ',IERR,IRWND)
c      Call Daxpy(Nbfns2,-Done,Work(Idens_b),1,Work(Idens_a),1)
c      call COMPPR(spens,Work(Idens_a) , propInt, Nbfns,.False.)
c         print*, 'spin Density for alpha e: ', Iocc, spens



      Return
      End 


 
