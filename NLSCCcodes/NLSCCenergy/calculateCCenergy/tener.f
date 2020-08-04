










      SUBROUTINE TENER(NLIST,ET2,EAA,NUMSYT,DISSYT,T2,W,T1A,
     &                 T1B,ISPIN,TAU,IRREP,POP1,POP2,VRT1,
     &                 VRT2,TMP)
C
C THIS ROUTINE CALCULATES THE CORRELATION ENERGY FOR A GIVEN
C SET OF AMPLITUDES (ALWAYS CALLED FOM CMPENG) FOR ONE
C SPECIFIC SPIN CASE
C
C    NLIST ..... OFFSET OF THE T2 LIST    
C    ET2 ....... CORRELATION ENERGY CALCULATED FROM T2 ONLY
C    EAA ....... CORRELATION ENERGY CALCULATED FROM TAU
C    NUMSYT .... NUMBER OF DISTRUBUTION IN T2 LIST
C    DISSYT .... DISTRIBUTION SIZE OF T2 LIST
C    T2 ........ HOLDS T2 AMPLITUDES
C    W ......... HOLDS INTEGRALS <AB//IJ>
C    T1A ....... HOLDS T1 AMPLITUDES (ALPHA)
C    T1B ....... HOLDS T1 AMPLITUDES (BETA)
C       FOR ISPIN=1,2 BOTH T1 ARRAYS ARE IDENTICAL
C    ISPIN ..... SPIN CASE
C    TAU ....... LOGICAl FLAG WHICH TELLS IF TAU AMPLITUDES ARE REQUIRED 
C    IRREP ..... IRREP OF T2 LIST
C    POP1 ...... POPULATION VECTOR OF I (OCCUPIED ORBITALS)
C    POP2 ...... POPULATIOn VECTOR OF J (OCCUPIED ORBITALS)
C    VRT1 ...... POPULATION VECTOR OF A (VIRTUAL ORBITALS)
C    TMP ....... SCRATCH ARRAY
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TAU
      INTEGER i,j,DISSYT,POP1,POP2,VRT1,VRT2,DIRPRD
      DIMENSION T1A(1),T1B(1),POP1(8),POP2(8),VRT1(8),VRT2(8)
      DIMENSION T2(DISSYT,NUMSYT),W(DISSYT,NUMSYT),TMP(1),
     &       NUMIRT(1),IBUFT(1)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      DATA AZERO,ONE /0.0D0,1.0D0/
      CALL GETLST(W, 1,NUMSYT,2,IRREP,13+ISPIN)
      CALL GETLST(T2,1,NUMSYT,1,IRREP,ISPIN+NLIST)


      ET2=DDOT(NUMSYT*DISSYT,T2,1,W,1)
      IF(TAU) THEN
       CALL FTAU(T2,T1A,T1B,DISSYT,NUMSYT,POP1,POP2,VRT1,VRT2,
     &           IRREP,ISPIN,ONE)
      ENDIF
      EAA=DDOT(NUMSYT*DISSYT,T2,1,W,1)
      write(*,"(a)") 
      Call checksum("T2", T2, NUMSYT*DISSYT)
      Call checksum("W ", W, NUMSYT*DISSYT)
      Write(*, "(a,1x,F15.10)") "W*T2 =", EAA
      RETURN
      END
