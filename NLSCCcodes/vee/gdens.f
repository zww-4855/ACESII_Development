










      SUBROUTINE GDENS(IRREPXL,IRREPXR,IRREPX,DOO,DVV,DVO,DOV,
     &                 SCR,MXCOR,IUHF,FACT,R0,L0,ZNORM,
     &                 LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &                 LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &                 LSTR2RS,LSTL2RS,LSTT1,LSTT1OFF)
C
C DRIVER FOR FORMATION OF GENERAL ONE-PARTICLE REDUCED DENSITY MATRIX
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L0
      DIMENSION SCR(MXCOR),LENVV(2),LENOO(2),LENVO(2)
      DIMENSION DOO(*),DVV(*),DVO(*),DOV(*)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONEM /-1.0D0/
C
C DO SOME ACCOUNTING AND SET UP ADDRESSES
C
      CALL IZERO(LENVV,2)
      CALL IZERO(LENOO,2)
      CALL IZERO(LENVO,2)
      DO 10 ISPIN=1,1+IUHF
       LENVV(ISPIN)=IRPDPD(IRREPX,18+ISPIN)
       LENOO(ISPIN)=IRPDPD(IRREPX,20+ISPIN)
       LENVO(ISPIN)=IRPDPD(IRREPX, 8+ISPIN)
10    CONTINUE
C
C INITIALIZE THE DENSITY MATRIX.  
C
C IF THE SYMMETRY OF THE RH STATE IS NOT TOTALLY SYMMETRIC, ZERO IT OUT.
C
C IF THE RH STATE IS TOTALLY SYMMETRIC, INITIALIZE WITH R0 CONTRIBUTION.
C
      CALL ZERO(DOO,LENOO(1)+IUHF*LENOO(2))
      CALL ZERO(DVV,LENVV(1)+IUHF*LENVV(2))
      CALL ZERO(DVO,LENVO(1)+IUHF*LENVO(2))
      CALL ZERO(DOV,LENVO(1)+IUHF*LENVO(2))
      IF(IRREPXR.EQ.1)THEN
       CALL GGSDEN(IRREPXL,DOO,DVV,DVO,DOV,
     &            LSTL1OFF,LSTL1,0,90,LISTL2,44,LSTGTL,LSTGTLOF,
     &            L0,SCR,MXCOR,IUHF)
       CALL SSCAL(LENOO(1)+IUHF*LENOO(2),R0,DOO,1)
       CALL SSCAL(LENVV(1)+IUHF*LENVV(2),R0,DVV,1)
       CALL SSCAL(LENVO(1)+IUHF*LENVO(2),R0,DVO,1)
       CALL SSCAL(LENVO(1)+IUHF*LENVO(2),R0,DOV,1)
      ENDIF
C
C                      +
C NOW CALCULATE <0|L {p q exp(T)} R_n|0> CONTRIBUTIONS TO THE
C DENSITY, WHERE n=1 AND 2.
C
C EVALUATE L2*R2 G-TYPE INTERMEDIATES AND PLACE ON LISTS 91 AND 92
C
      IF(LSTGRL.NE.-1)THEN
       CALL GFORMG(IRREPXL,IRREPXR,LISTL2,LISTR2,LSTGRL,SCR,
     &             MXCOR,LSTGRLOF,ONEM,ONEM,IUHF)
      ENDIF
C
C EVALUATE L2*T2 G-TYPE INTERMEDIATES AND PLACE ON LISTS 491 AND 492
C
      CALL GFORMG(IRREPXL,1,LISTL2,44,LSTGTL,SCR,
     &            MXCOR,LSTGTLOF,ONEM,ONEM,IUHF)
C
C OCCUPIED-VIRTUAL TERMS
C
      CALL GDENSOV(IRREPXL,IRREPXR,IRREPX,DOV,R0,SCR,
     &             MXCOR,IUHF,
     &             LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &             LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &             LSTR2RS,LSTL2RS)
C
C VIRTUAL-VIRTUAL AND OCCUPIED-OCCUPIED TERMS
C
      CALL GDENS1 (IRREPXL,IRREPXR,IRREPX,DOO,DVV,
     &             DOV,R0,SCR,MXCOR,IUHF,
     &             LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &             LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &             LSTR2RS,LSTL2RS,LSTT1,LSTT1OFF)
C
C VIRTUAL-OCCUPIED TERMS
C
      CALL GDENSVO(IRREPXL,IRREPXR,IRREPX,DVO,DOV,R0,L0,ZNORM,
     &             SCR,MXCOR,IUHF,
     &             LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &             LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &             LSTR2RS,LSTL2RS)
      
C
      RETURN
      END
