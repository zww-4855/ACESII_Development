










      SUBROUTINE AUGMENTP(SCR, MAXCOR, ISIDE, IUHF, IRREPX)
C
C  AUGMENT OVERLAP MATRIX 'P' OF VECTORS PROJECTED ON
C      EXCITATION PATTERN (ON ICOLPR1 LISTH0)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SCR(MAXCOR)
C
      PARAMETER (MAXORD=100)
C
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/PMAT/ P(10000)
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
C
      INDXF(I,J,N)=I+(J-1)*N
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      LISTC = 470
C
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      I030=I020+NSIZEC
      CALL GETLST(SCR(I010),IGET(1),1,1,ISIDE,LISTC)
      CALL GETLST(SCR(I020),ICOLPR1,1,1,1,LISTH0)
      CALL VECPRD(SCR(I020),SCR(I010),SCR(I000),NSIZEC)
      IF(IUHF.EQ.0)THEN
       CALL SPNTSING(IRREPX,SCR(I000),SCR(I010),MAXCOR-I010+1)
      ENDIF
      DO 5 I=1,NDIMR
       CALL GETLST(SCR(I010),IGET(I),1,1,ISIDE,LISTC)
       P(INDXF(1,I,MAXORD))=SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
       P(INDXF(I,1,MAXORD)) =  P(INDXF(1,I,MAXORD))
 5    CONTINUE 
C
      RETURN
      END