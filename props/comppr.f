










      SUBROUTINE COMPPR(PROP,DENS,PRPINT,NSIZ,NUCLEAR)
C
C THIS ROUTINE READS A LIST OF PROPERTY INTEGRALS AND THE DENSITY
C  MATRIX AND COMPUTES THE PARTICULAR PROPERTY.
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ICORE
      PARAMETER (MAXCOR=10 000 000)
      DIMENSION ICORE(MAXCOR)
      LOGICAL NUCLEAR
      CHARACTER*3  TYPE(3),FRONT
      CHARACTER*32 CRAP
      DIMENSION DENS(NSIZ,NSIZ),PRPINT(NSIZ,NSIZ),BUF(600),IBUF(600)
      DIMENSION IXX(2)
      COMMON /SYMINF/ NSTART,NIRREP
      EQUIVALENCE (PRPNUC,IXX(1))
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT
      real(kind=8)::summed,summedB,diff,total, sumMO
      integer :: na,nb,nelec,noccorb(2)
      integer :: a,b,NSO,NAO,counted
      real(kind=8),allocatable::tad(:,:),tbd(:,:),scfevca0(:)
      real(kind=8),allocatable::scfA(:,:),scfB(:,:),scfevcb0(:)
      real(kind=8),allocatable::tempDen(:,:),scfdensa(:),scfdensb(:)
      DATA ONE /1.0D0/,TWO/2.D0/
        NNM1O2(IX)=(IX*(IX-1))/2
      IEXTI(IX)=1+(-1+INT(DSQRT(8.D0*IX+0.999D0)))/2
      IEXTJ(IX)=IX-NNM1O2(IEXTI(IX))
C
C FIRST GET NUCLEAR CONTRIBUTION BY BACKSPACING FILE.
C
c      print*, 'inside commpp.f'
      BACKSPACE(10)
      READ(10)CRAP,PRPNUC
      IF(.NOT.NUCLEAR) PRPNUC=0.D0
c      print*, 'reading prop ints'
C
C READ IN THE PROPERTY INTEGRALS
C
      CALL ZERO(PRPINT,NSIZ*NSIZ)
1     READ(10)BUF,IBUF,NUT
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 10 I=1,NUT
c       print*, 'is i ', I
       INDI=IEXTI(IBUF(I))
       INDJ=IEXTJ(IBUF(I))
       PRPINT(INDI,INDJ)=BUF(I)
       PRPINT(INDJ,INDI)=BUF(I)
10    CONTINUE
      IF(NUT.EQ.600)GOTO 1
      print*, 'first element of prpint', PRPINT(1,1)
      print*,'first element of density atrix;',DENS(1,1)
C
C COMPUTE THE PROPERTY
C
c      total=0.0d0
c      do i=1,nsiz
c        sumMO=0.0d0
c        do j=1,nsiz
c          sumMO=sumMO+DENS(i,j)*PRPINT(j,i)
c        enddo
cc        write(*,*) 'contribution from MO:', i, sumMO
c        total=total+sumMO
c      enddo
c      write(*,*) 'Trace of DENS matrix:', total+PRPNUC
c       courtesy of pop.f line 306
c      IDMOA=1
c      IUHF=1
c      call getrec(20,"JOBARC","NBASTOT",1,NAO)
c      NCOMP=NAO
c      IDMOB=1+IUHF*NCOMP*NCOMP
c      NORBS2=NAO*NAO
c       IDAO=IDMOB+NAO*NAO
c       ISCR=IDAO+NAO*NAO
c
c         INBASIR=ISCR
c         IOCCA=INBASIR+NIRREP
c         IOCCB=IOCCA
c         SCFOCC=ONE
c         CALL ZERO(ICORE(IDMOA),NCOMP*NCOMP)
c         CALL GETREC(20,'JOBARC','NUMBASIR',NIRREP,ICORE(INBASIR))
c         CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,ICORE(IOCCA))
c         CALL FILLDHF2(SCFOCC,ICORE(IDMOA),NCOMP,ICORE(IOCCA),
c     &                 ICORE(INBASIR))
c
c          CALL ZERO(ICORE(IDMOB),NCOMP*NCOMP)
c          CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,ICORE(IOCCB))
c          CALL FILLDHF2(SCFOCC,ICORE(IDMOB),NCOMP,ICORE(IOCCB),
c     &                  ICORE(INBASIR))
c      print*,'new comparison:', ICORE(1),ICORE(IDMOB)
c      z=0
c      k=0
c      m=0
c      do i=1,NAO
c        do j=1,NAO
c        print*, 'diagonal of alpha/beta',ICORE(1+z),ICORE(IDMOB+z)
c        k=k+ICORE(1+z)
c        m=m+ICORE(IDMOB+z)
c        z=z+1
c        enddo
c      enddo
c prints number of alpha/beta occ
c      print*, 'trace of diag', k,m
c     Harvest alpha/beta eigenvectors 
c      call getrec(20,"JOBARC","NAOBFORB",1,NSO)
c      write(*,*) 'number of spin orbs', NSO
c      write(*,*) 'number of atomic orbs', NAO
c      allocate(scfevca0(NAO*NAO),scfevcb0(NAO*NAO))
c      allocate(scfdensa(NSO*NSO),scfdensb(NSO*NSO))
cc      write(*,*)'alpha molecular orbs:', scfevca0
c      call getrec(20,"JOBARC","SCFEVCA0",NAO*NAO,scfevca0)
c      call getrec(20,"JOBARC","SCFEVCB0",NAO*NAO,scfevcb0)
c      call getrec(20,"JOBARC", "NOCCORB",2, noccorb)
c      call getrec(20,"JOBARC", "SCFDENSA",NSO*NSO,scfdensa)
c      call getrec(20,"JOBARC", "SCFDENSB",NSO*NSO,scfdensb)
c      nelec=noccorb(1)+noccorb(2)
c      print*, "number of electrons is: ", nelec
c      do i=1,NAO
c          write(*,*) i,scfevca0(i)
c      enddo
c
cc   Check alpha density + beta density matrix is same as overall Den
c      allocate(tad(NAO,NAO),tbd(NAO,NAO),tempDen(NAO,NAO))
c      allocate(scfA(NAO,NAO),scfB(NAO,NAO))
c      counted=1
cc       alpha/beta coeff vectors in matrix form--scfa.scfb
cc       scfa/scfb are alpha/beta density matrix
c      do i=1,NAO
c        do j=1,NAO
c          tad(i,j)=scfevca0(counted)
c          tbd(i,j)=scfevcb0(counted)
c          counted=counted+1
cc          scfA(i,j)=tad(j,i)*tad(i,j)
cc          scfB(i,j)=tbd(j,i)*tbd(i,j)
cc          tempDen(i,j)=scfA(i,j)-scfB(i,j)
c          write(*,*) 'nad ij',  tad(i,j)
c        enddo
c      enddo
c
cc    Add alpha/beta density matrix and compare with real AO
cc    Density Matrix
c      na=noccorb(1)
c      nb=noccorb(2)
c      do mu=1,NAO
c        do nu=1,NAO
c          summed=0.0
c          summedB=0.0
c          do a=1,NAO
c            summed=summed+tad(nu,a)*tad(mu,a)
c          enddo 
c          do b=1,NAO
c            summedB=summedB+tbd(nu,b)*tbd(mu,b)
c          enddo
c          scfA(mu,nu)=summed
c          scfB(mu,nu)=summedB
c          tempDen(mu,nu)=scfA(mu,nu)-scfB(mu,nu)
c        enddo
c      enddo
c
c      do i=1,NAO
c        do j=1,NAO
c          diff=abs(tempDen(i,j) - DENS(i,j))
c          if (diff .gt. 0.0001) then
c            write(*,*)'element i,j gt crit',i,j,diff
c          endif
c        enddo
c      enddo
c      summed=0.0
c      summedB=0.0
c      do i=1,NAO
c        summed=summed+DENS(i,i)
c        summed=summedB+scfA(i,i)+scfB(i,i)
c      enddo
c      print*,'trace of alpha density mat', summed,summedB
c      scfA=0.0d0
c      do i=1,NSO
c        scfA(i,1)=scfdensA(i)
c      enddo
c      do i=1,NSO
c        diff=abs(scfA(i,1)-scfdensA(i))
c        print*, 'compare dens and scfdensA',scfA(i,1),scfdensA(i),
c     $                          scfdensB(i), DENS(i,1)
c        if (diff .gt. 0.00001) then
c          print*, 'error in mo basis setup'
c        endif
c      enddo
c      print*, 'compare',scfdensa(1),scfdensb(1),scfA(1,1)
c      print*, 'write nsize is',NSIZ,DENS(1,1),scfA(1,1),scfB(1,1)
c      deallocate(tempDen,scfdensa,scfdensb)
c      deallocate(scfevca0,scfevcb0,tad,tbd,scfA,scfB)  
      PROP=SDOT(NSIZ*NSIZ,DENS,1,PRPINT,1)+PRPNUC
      RETURN
      END
