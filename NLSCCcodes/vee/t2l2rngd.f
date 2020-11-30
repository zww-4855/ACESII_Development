      SUBROUTINE T2L2RNGD(ICORE,MAXCOR,IUHF,R0)
C
C   THIS ROUTINE CALCULATES THE FOLLOWING INTERMEDIATE 
C
C   H(ME,IA) = - SUM N,F  L(MN,EF) T(NI,FA)
C
C  THIS TERM IS VERY SIMILAR TO THE T1(IJ,AB) CONTRIBUTION TO
C  THE W-RING INTERMEDIATE
C
C  THE SPIN CASES ARE
C
C    AAAA : =  - SUM M,E L(IM,BE) T(MJ,EA) - SUM m,e L(Im,Be) T(Jm,Ae)
C
C    ABAB : =  - SUM m,E L(Im,Eb) T(Jm,Ea)
C
C    ABBA : =  - SUM M,E L(IM,BE) T(Mj,Ea) - SUM m,e L(Im,Be) T(mj,eb)
C
C FOR UHF IN ADDITION THE BBBB BABA AND BABA SPIN CASES HAS TO CALCULATES
C
C
C  AAAA : LIST 414  (UHF ONLY)
C  BBBB : LIST 415  (UHF ONLY)
C  ABBA : LIST 418
C  BAAB : LIST 419  (UHF ONLY)
C  ABAB : LIST 416
C  BABA : LIST 417  (UHF ONLY)
C
CEND
C 
C  CODED SEPTEMBER/93 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ICORE(MAXCOR)
C
      COMMON/STATSYM/IRREPX
C
C CREATE LIST AREAS (414 -- 419)
C
      CALL SETLST(ICORE,MAXCOR,IUHF,IRREPX)
C
C CALL T2L2RNG TO CALCULATE THE VARIOUS H INTERMEDIATES
C
      CALL T2L2RNG(ICORE,MAXCOR,'ABBA',IUHF,IRREPX,R0)
      CALL T2L2RNG(ICORE,MAXCOR,'ABAB',IUHF,IRREPX,R0)
      IF(IUHF.EQ.1) THEN
       CALL T2L2RNG(ICORE,MAXCOR,'AAAA',IUHF,IRREPX,R0)
      ELSE
       CALL QUIKAA3(ICORE,MAXCOR,IRREPX,414)
      ENDIF
      IF(IUHF.NE.0) THEN
       CALL T2L2RNG(ICORE,MAXCOR,'BBBB',IUHF,IRREPX,R0)
       CALL T2L2RNG(ICORE,MAXCOR,'BABA',IUHF,IRREPX,R0)
C
C  BAAB HAS BEEN CALCULATED ALREADY IN THE CALL WITH ABBA
C
      ENDIF
C
C RESORT H INTERMEDIATES
C
      CALL SORTH(ICORE,MAXCOR,IUHF,414,IRREPX)
C
      RETURN
      END
