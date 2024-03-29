

















c Macros beginning with M_ are machine dependant macros.
c Macros beginning with B_ are blas/lapack routines.

c  M_REAL         set to either "real" or "double precision" depending on
c                 whether single or double precision math is done on this
c                 machine

c  M_IMPLICITNONE set iff the fortran compiler supports "implicit none"
c  M_TRACEBACK    set to the call which gets a traceback if supported by
c                 the compiler























cYAU - ACES3 stuff . . . we hope - #include <aces.par>





      SUBROUTINE EDENPRC(DDIFF,SECMOM,SCR,NMO,NAO,IUHF,ISPIN,ITYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD
      CHARACTER*8 LBLVEC(2),LBLORB(2),LBLOCC(2),LBL2ND(6)
      CHARACTER*5 SPTYPE(2)
      LOGICAL PRINT1,PRINT2,VPROP
      DIMENSION DDIFF(NMO,NMO),SECMOM(6,NMO),SCR(*),MAP(6)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/  IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /VDINTLEN/ LENS(8),LENT(8)
      COMMON/AOSYM/IAOPOP(8),IOFFAOX(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON /INTPROG/ VPROP
C
      DATA MAP /1,4,2,5,6,3/
      DATA ONE,ZILCH /1.0D0,0.0D0/
      DATA LBLVEC /'SCFEVECA','SCFEVECB'/
      DATA LBLORB /'DIFFORBA','DIFFORBB'/
      DATA LBLOCC /'OCCNUM_A','OCCNUM_B'/
      DATA LBL2ND /'2NDMO_XX','2NDMO_YY','2NDMO_ZZ',
     &             '2NDMO_XY','2NDMO_XZ','2NDMO_YZ'/
      DATA SPTYPE /'alpha','beta '/
C
      NNP1O2(I)=(I*(I+1))/2
C
      PRINT1=IFLAGS(1).GE.100
      PRINT2=.NOT.PRINT1.AND.IFLAGS(1).GE.2
C
      IOFFMO=1
      IOFFAO=IOFFMO+NMO*NMO
      IOFFE =IOFFAO+NAO*NAO
      IOFFS1=IOFFE+NAO*NAO
      IOFFS2=IOFFS1+NAO*NAO
      CALL SSCAL(NMO*NMO,2.0D0-DFLOAT(IUHF),DDIFF,1)
C
C DIAGONALIZE DENSITY DIFFERENCE MATRIX
C
cYAU - old
c      CALL EIG2(SCR,DDIFF,NMO,NMO,-1)
cYAU - new
      IF (NMO.EQ.1) THEN
         SCR(1)     = DDIFF(1,1)
         DDIFF(1,1) = 1.0d0
      ELSE
         CALL dsyev('V','L',NMO,DDIFF,NMO,SCR(1),SCR(1+NMO),3*NMO,ierr)
      END IF
cYAU - end
C
C TRANSFORM EIGENVECTORS TO THE AO BASIS AND WRITE TO JOBARC FILE
C
      CALL GETREC(20,'JOBARC',LBLVEC(ISPIN),NMO*NAO*IINTFP,SCR(IOFFE))
      CALL XGEMM('N','N',NAO,NMO,NMO,ONE,SCR(IOFFE),NAO,DDIFF,NMO,
     &           ZILCH,SCR(IOFFAO),NAO)
      CALL PUTREC(20,'JOBARC',LBLORB(ISPIN),NAO*NMO*IINTFP,SCR(IOFFAO))
      CALL PUTREC(20,'JOBARC',LBLOCC(ISPIN),NMO*IINTFP,SCR)
C
C CALCULATE SECOND MOMENT OF CHARGE DISTRIBUTION FOR THESE ORBITALS
C
      IDIR=0
      DO 20 IXYZ1=1,3
       DO 21 IXYZ2=1,IXYZ1
        IDIR=IDIR+1
        IF(VPROP)THEN
         LENGTH=NNP1O2(NAO)
         CALL GETREC(20,'JOBARC',LBL2ND(MAP(IDIR)),IINTFP*LENGTH,
     &               SCR(IOFFS1))
         CALL EXPND2(SCR(IOFFS1),SCR(IOFFE),NAO) 
        ELSE
         IRRPRT=DIRPRD(IDIPSM(IXYZ1),IDIPSM(IXYZ2))
         LENGTH=LENT(IRRPRT)
         CALL GETREC(20,'JOBARC',LBL2ND(MAP(IDIR)),IINTFP*LENGTH,
     &               SCR(IOFFS1))
         CALL VMINUS(SCR(IOFFS1),LENGTH)
         CALL ZERO(SCR(IOFFE),NAO*NAO)
         CALL MATEXP(IRRPRT,IAOPOP,SCR(IOFFS1),SCR(IOFFE))
         CALL MATEXP2(IRRPRT,SCR(IOFFE),SCR(IOFFS1),NAO)
         CALL SCOPY  (NAO*NAO,SCR(IOFFS1),1,SCR(IOFFE),1)
        ENDIF
        ILOC=IOFFAO
        DO 22 IMO=1,NMO
         CALL XGEMM ('N','T',NAO,NAO,1,ONE,SCR(ILOC),NAO,SCR(ILOC),
     &               NAO,ZILCH,SCR(IOFFS1),NAO)
         SECMOM(MAP(IDIR),IMO)=SDOT(NAO*NAO,SCR(IOFFS1),1,SCR(IOFFE),1)
         ILOC=ILOC+NAO
22      CONTINUE
21     CONTINUE
20    CONTINUE      
C
      IF(PRINT1)THEN
       IF(ITYPE.EQ.1)THEN
        WRITE(6,1005)SPTYPE(ISPIN)
       ELSE
        WRITE(6,1006)SPTYPE(ISPIN)
       ENDIF
       DO 13 I=1,NMO
        WRITE(6,*)' Eigenvalue : ',SCR(I)
        WRITE(6,*)' Eigenvector : '
        WRITE(6,'((6F12.8))')(DDIFF(J,I),J=1,NMO)
13     CONTINUE
      ELSE
       IF(ITYPE.EQ.1)THEN
        WRITE(6,1005)SPTYPE(ISPIN)
       ELSE
        WRITE(6,1006)SPTYPE(ISPIN)
       ENDIF
       WRITE(6,1000)
       WRITE(6,1010)
       WRITE(6,1000)
       DO 14 I=1,NMO
        IX=ISAMAX(NMO,DDIFF(1,I),1)
        IF(ABS(SCR(I)).GT.0.075.OR.PRINT2)THEN
         WRITE(6,1011)I,SCR(I),IX,ABS(DDIFF(IX,I)),
     &                (SECMOM(IK,I),IK=1,6)
        ENDIF
14     CONTINUE
       WRITE(6,1000)
      ENDIF
C
      RETURN
1000  FORMAT(72('-'))
1001  FORMAT(T3,'Orbital',T19,'Eigenvalue',T35,'Major Component',
     &       T58,'Contribution')
1002  FORMAT(T5,I3,T19,F8.5,T40,I3,T60,F8.5)
1003  format(6f10.6)
1005  FORMAT(T3,'Spin case ',A5,' natural difference orbitals:')
1006  FORMAT(T3,'Spin case ',A5,' occupation numbers.')
1010  FORMAT(T2,'Orbital',T11,'Occ. Major',T43,'Second Moment',/,
     &       T11,'Num. Cont.',/,
     &       T29,'xx',T37,'yy',T45,'zz',T53,'xy',T61,'xz',T69,'yz')
1011  FORMAT(T3,I4,T9,F6.3,T16,I3,T20,'[',F4.2,']',6(1X,F7.3))
      END
