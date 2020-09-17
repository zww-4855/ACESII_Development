










      SUBROUTINE AFTERRG(NORDER,EVALR,EVALI,EVEC,SCR,ISCR,IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EVALR(NORDER),EVEC(NORDER,NORDER),SCR(NORDER,NORDER)
      DIMENSION ISCR(NORDER),EVALI(NORDER)
C
      DATA ONE/1.D0/
C
      DO 10 I=1,NORDER
       ISCR(I)=I
10    CONTINUE
      CALL PIKSR2(NORDER,EVALR,ISCR)
C
      DO 11 I=1,NORDER
       J=ISCR(I)
       X=ONE/SNRM2(NORDER,EVEC(1,J),1)
       CALL SSCAL(NORDER,X,EVEC(1,J),1)
       CALL SCOPY(NORDER,EVEC(1,J),1,SCR(1,I),1)
11    CONTINUE
      CALL SCOPY(NORDER*NORDER,SCR,1,EVEC,1)
C
      DO 13 I=1,NORDER
       J=ISCR(I)
       CALL SCOPY(1,EVALI(J),1,SCR(I,1),1)
13    CONTINUE
      CALL SCOPY(NORDER,SCR,1,EVALI,1)
C
CSSS      IF(IPRINT.GE.10)THEN
CSSS      ENDIF

      IF(IPRINT.GT.20)THEN
CSSS       DO 12 I=1,NORDER
       DO 12 I=1,1
        WRITE(6,1001)I
        WRITE(6,'((8F10.7))')(EVEC(K,I),K=1,NORDER)
1001   FORMAT(T3,'Eigenvector #',I3,':') 
12     CONTINUE
      ENDIF
      RETURN
      END
