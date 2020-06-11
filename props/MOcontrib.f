










      SUBROUTINE MOcontrib(NAO,DENS,NSO,scfdensa,scfdensb,na,nb,PRPINT)
C       TAKES THE OVERALL DENSITY MATRIX
C       RETURNS THE DENSITY MATRIX SPECIFIC TO A PARTICULAR ORBITAL      
C       SEE SZABO PG213
C       PURPOSE: BREAK OVERALL DENS MATRIX INTO ORBITAL CONTRIBUTIONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer :: na,nb,nelec,noccorb(2)
      integer :: a,b,NSO,NAO
      DIMENSION::DENS(NAO,NAO),cvecDENS(NSO,NSO),cDENS(NAO,NAO)
      real(kind=8)::scfdensa(NSO*NSO),scfdensb(NSO*NSO)
      real(kind=8)::cbDENS(NSO,NSO),scfMOa(NSO,NSO),scfMOb(NSO,NSO)
      real(kind=8)::scfMOtot(NSO,NSO),MOspens
        COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
            DIMENSION BUF(600),IBUF(600),PRPINT(NAO,NAO)
        COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT
      COMMON /SYMINF/ NSTART,NIRREP
        real(kind=8)::total
        integer::junk
        NNM1O2(IX)=(IX*(IX-1))/2
      IEXTI(IX)=1+(-1+INT(DSQRT(8.D0*IX+0.999D0)))/2
      IEXTJ(IX)=IX-NNM1O2(IEXTI(IX))

1     READ(10)BUF,IBUF,NUT
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 10 I=1,NUT
       INDI=IEXTI(IBUF(I))
       INDJ=IEXTJ(IBUF(I))
       PRPINT(INDI,INDJ)=BUF(I)
       PRPINT(INDJ,INDI)=BUF(I)
10    CONTINUE
      IF(NUT.EQ.600)GOTO 1
C

      cvecDENS=0.0d0
c      call DCOPY(NAO*NAO,DENS,1,cDENS,1)
      z=1
      do i=1,NSO
        do j=1,NSO
          cDENS(i,j)=scfdensa(z)
          cbDENS(i,j)=scfdensb(z)
          z=z+1
        enddo
      enddo
c      cDens=scfdensa
c      call eig(cDENS,cvecDENS,junk,NSO,0)
      print*, 'cvecs dens', cDENS(1,1),DENS(1,1)
      MOspens=0.0d0
      print*,'compare',DENS(1,1),cDENS(1,1),cbDENS(1,1)
      do i=1,NSO
        do j=1,NSO
          MOspens=abs(DENS(i,j))-abs(cDENS(i,j)-cbDENS(i,j))
          if (abs(MOspens) .gt. 0.0001 ) then
            print*, 'error in element',i,j,MOspens
        endif
      enddo
        enddo
      print*, 'trace is: ', MOspens
      total=0.0d0
      do j=1,NSO
        z=0.0d0
        do i=1,NSO
          total=total+cDENS(i,j)*PRPINT(j,i)
          z=z+cDENS(i,j)*PRPINT(j,i)
        enddo
        print*, 'MO contribution of: ', j, z
      enddo
      print*,'final realdens mat mo1',total
        
      z=0
      total=0.0d0
      do i=1,NSO
        do j=1,NSO
          scfMOb(i,j)=scfdensb(z)
          z=z+1
             total=total+scfMOb(i,j)
          if (i .eq. j) then
             print*, "diagonal", scfMOb(i,j)
          endif
        enddo
      enddo
      print*, 'sum of diagonal b',total 
      do i=1,NSO
        scfMOa(i,1)=scfdensA(i)
        scfMOb(i,1)=scfdensB(i)
      enddo
      scfMOtot=0.0d0
      scfMOtot=scfMOa-scfMOb
      print*,'difference between files:',scfMOa(1,1)-scfMOb(1,1)
c      call comppr(MOspens,scfMOtot,SCR,NSO,.FALSE.)
      do i=1,NSO
        z=0.0
        do j=1,NSO
          print*,'prpint',i,j, PRPINT(j,i),scfMOtot(i,j)
          z=z+scfMOtot(i,j)*PRPINT(j,i)
        enddo
        total=total+z
      enddo
      c=0.0d0
      do i=1,NSO
        c=c+scfMOtot(i,1)*PRPINT(i,1)
      enddo
      
      WRITE(*,*) 'mo1', c,total
      RETURN
      END
