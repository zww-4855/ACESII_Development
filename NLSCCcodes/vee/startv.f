      SUBROUTINE STARTV(Z, GUESSV, NSIZEC, NDIMR, IOLDEST, MAXORD, SCR,
     &   IRREPX, IUHF, ISIDE, IPRINT, IREAD, IJUNK)
C
C THE EXPANSION COEFFICIENTS OF GUESSV IN TERMS OF THE BASIS VECTORS
C IS CALCULATED
C
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION Z(100), GUESSV(NSIZEC), SCR(*), SDOT
      DOUBLE PRECISION ROOT 
      DIMENSION IJUNK(100)
      LOGICAL PRINT
C
      IGET(I)=1+MOD(IOLDEST+MAXORD-I,MAXORD+1)
C
      PRINT = IPRINT .GE. 50
C
      CALL ZERO(Z,100)
      CALL GETGES(GUESSV,NSIZEC,IREAD,IRREPX,ROOT,IJUNK)
      IF (PRINT) THEN
         WRITE(6,*) ' ORIGINAL  GUESSV'
         CALL OUTPUT(GUESSV,1,1,1,NSIZEC,1,1,1)
       ENDIF
       IF(IUHF.EQ.0) CALL SPNTSING(IRREPX,GUESSV,SCR,2*NSIZEC)
C
       DO 10 I=1,NDIMR
        CALL GETLST(SCR,IGET(I),1,1,ISIDE,470) 
        Z(I)=SDOT(NSIZEC,GUESSV,1,SCR,1)
   10  CONTINUE
C
       IF (PRINT) THEN
C
C  Reconstruct GUESSV from Z-coefficients
C
         I000 = 1
         I010 = I000 + NSIZEC
         CALL ZERO(SCR(I000),NSIZEC)
         DO 20 I=1,NDIMR
           CALL GETLST(SCR(I010),IGET(I),1,1,ISIDE,470) 
           CALL VADD(SCR(I000),SCR(I000),SCR(I010),NSIZEC, Z(I))
   20    CONTINUE
         WRITE(6,*) ' RECONSTRUCTED GUESSV'
         CALL OUTPUT(SCR(I000),1,1,1,NSIZEC,1,1,1)
C
         WRITE(6,*) ' NEW EXPANSION VECTOR IN STARTV '
         CALL OUTPUT(Z,1,ndimr,1,1,1,ndimr,1)
       ENDIF
C
       RETURN
       END