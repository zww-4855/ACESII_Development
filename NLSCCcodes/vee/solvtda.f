      SUBROUTINE SOLVTDA(ALARGE, NDIML, ASMALL, NDIMS, EVECL, NUMSOL,
     &   EVECS, EVALS, INDEX,TRIPLET, IRREPX, IUHF, SCR, MAXCOR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL TRIPLET, PRINT
C
      DIMENSION ALARGE(NDIML,NDIML), ASMALL(NDIMS, NDIMS),
     &   EVECL(NDIML, NUMSOL), EVECS(NDIMS, NDIMS), EVALS(NDIMS),
     &   INDEX(NDIMS), SCR(MAXCOR)
c
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH /0.0D0/
C
      PRINT = IFLAGS(1) .GE. 5
      ONE = 1.0D0
C
C
C  READ IN THE EXCITATION PATTERNS
C
       IF(IUHF.EQ.0) THEN
        NTDA=IRPDPD(IRREPX,9)
        CALL GETLST(SCR,1,1,1, 1, 490)
      ELSE
        NUMAA=IRPDPD(IRREPX,ISYTYP(1,19))
        NUMBB=IRPDPD(IRREPX,ISYTYP(1,20))
        NTDA = NUMAA + NUMBB
        I000 = 1
        I010 = I000 + NUMAA
        I020 = I010 + NUMBB
        CALL GETLST(SCR(I000),1,1,1, 1, 490)
        CALL GETLST(SCR(I010),1,1,1, 2, 490)
      ENDIF
C
      IF (NTDA .NE. NDIML) THEN
        WRITE(6,*) ' SOMETHING WRONG IN SOLVTDA', NDIML, NTDA
        CALL ERREX
      ENDIF
C
      ICOUNT = 0
      DO I = 1, NDIML
        IF (SCR(I) .GT. 0.5) THEN
          ICOUNT = ICOUNT + 1
          INDEX(ICOUNT) = I
        ENDIF
      ENDDO
C
      IF (ICOUNT .NE. NDIMS) THEN
        WRITE(6,*) ' SOMETHING WRONG IN SOLVTDA', ICOUNT, NDIMS
        CALL ERREX
      ENDIF
C
C FILL PROJECTED TDA MATRIX
C
      DO I = 1, NDIMS
        DO J = 1, NDIMS
          ASMALL(I,J) = ALARGE(INDEX(I),INDEX(J))
        ENDDO
      ENDDO
C
      IF (PRINT) THEN
        WRITE(6,*) ' PROJECTED TDA MATRIX '
        CALL OUTPUT(ASMALL, 1, NDIMS, 1, NDIMS, NDIMS, NDIMS, 1)
      ENDIF
C
C DIAGONALIZE ASMALL
C
       CALL EIG(ASMALL,EVECS,NDIMS,NDIMS,0)
       DO I = 1, NDIMS
         EVALS(I) = ASMALL(I,I)
       ENDDO
C
       IF (PRINT) THEN
         WRITE(6,*) ' EIGENVALUES IN PROJECTED TDA '
         CALL OUTPUT(EVALS, 1, 1, 1, NDIMS, 1, NDIMS, 1)
         WRITE(6,*) ' EIGENVECTORS IN PROJECTED TDA '
         CALL OUTPUT(EVECS, 1, NDIMS, 1, NDIMS, NDIMS, NDIMS, 1)
       ENDIF
C
C     SELECT EIGENVALUES
C
       IF (TRIPLET) THEN
C
C  SET EIGENVALUES CORRESPONDING TO SINGLETS TO HUGE
C
         HUGE = 1.0D8
         IF (MOD(NDIMS,2) .NE. 0) THEN
           WRITE(6,*)' SOMETHING WRONG IN SOLVTDA', NDIMS
           CALL ERREX
         ENDIF
         NDIMA = NDIMS/2
         DO I = 1, NDIMS
           DIFF = ZILCH
           DO J = 1, NDIMA
             DIFF = DIFF + ABS(EVECS(J,I) - EVECS(J+NDIMA,I))
           ENDDO
           IF (DIFF .LT. 0.01) THEN
             EVALS(I) = HUGE
           ENDIF
         ENDDO
C
       IF (PRINT) THEN
         WRITE(6,*) ' EIGENVALUES IN PROJECTED TDA AFTER TRIPLET'
         CALL OUTPUT(EVALS, 1, 1, 1, NDIMS, 1, NDIMS, 1)
       ENDIF
C
       ENDIF
       I000 = 1
       I010 = I000 + NDIMS * NDIMS
C
       CALL ORDERTDA(NDIMS, NDIMS, EVALS, EVECS, SCR(I000),
     &    SCR(I010))
C
       IF (PRINT) THEN
         WRITE(6,*) ' EIGENVALUES IN TDA AFTER ORDERTDA'
         CALL OUTPUT(EVALS, 1, 1, 1, NDIMS, 1, NDIMS, 1)
         WRITE(6,*) ' EIGENVECTORS IN PROJECTED TDA '
         CALL OUTPUT(EVECS, 1, NDIMS, 1, NDIMS, NDIMS, NDIMS, 1)
       ENDIF
C
C  CONSTRUCT EVECL
C
       CALL ZERO(EVECL, NUMSOL*NDIML)
       DO J = 1, NUMSOL
         DO I = 1, NDIMS
           EVECL(INDEX(I),J) = EVECS(I,J)
         ENDDO
       ENDDO
C
C  NORMALIZE EIGENVECTORS
C
       DO J = 1, NUMSOL
        Z=SNRM2(NDIML,EVECL(1,J),1)
        CALL SSCAL(NDIML,ONE/Z,EVECL(1,J),1)       
      ENDDO
C
       IF (PRINT) THEN
         WRITE(6,*) ' FULL EIGENVECTORS IN PROJECTED TDA '
         CALL OUTPUT(EVECL, 1, NDIML, 1, NUMSOL, NDIML, NUMSOL, 1)
       ENDIF
C
       RETURN
       END