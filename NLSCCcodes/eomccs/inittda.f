      SUBROUTINE INITTDA(ICORE,IUHF)
C
C  THIS WAS PART OF SETMET, BUT DID NOT WORK WITH THE NEW LIST STRUCTURE.
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,EOMCC,CISD,FULDIAG,INCORE,READGUES,DOUBLE,NONSTD
      LOGICAL RPA
      DIMENSION ICORE(*)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /CALCINFO/ NROOT(8)
      COMMON /GUESS/ DOUBLE,NONSTD
      COMMON /GUESS2/ IMAP(100,8)
      COMMON/TDALIST/LISTETDA, LISTVTDA
C
      DO 33 IRREP=1,NIRREP  
C        DO 32 I=1,NROOT(IRREP)
C         IMAX=MAX(IMAX,IMAP(I,IRREP))
C  32      CONTINUE
C        IF(.NOT.FULDIAG)THEN
C         IF(.NOT.NONSTD)THEN
C          CALL UPDMOI(IMAX+3,1,IRREP,LISTETDA,0,0)
         IF(CIS.AND.NONSTD)THEN
          CALL UPDMOI(NROOT(IRREP),1,IRREP,LISTETDA,0,0)
         ENDIF
C        ENDIF
33    CONTINUE
C
      RETURN
      END



