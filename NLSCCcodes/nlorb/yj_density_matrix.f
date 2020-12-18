      SUBROUTINE YJ_DENSITY_MATRIX (COEFFS,NBAS)

      IMPLICIT NONE

      INTEGER NBAS,NPROT,NOCC,I,J,K


      DOUBLE PRECISION COEFFS(NBAS,NBAS)
      DOUBLE PRECISION P(NBAS,NBAS)


      CALL GETREC(10,"JOBARC","NMPROTON",1,NPROT)

      NOCC=NPROT/2

      DO I=1,NBAS
       DO J=1,NBAS
        P(I,J)=0
        DO K=1,NOCC
         P(I,J)=P(I,J)+2*COEFFS(I,K)*COEFFS(J,K)
        END DO
       END DO
      END DO

         CALL  PRINTM
     +
     +                     ("Current density matrix",
     +                      NBAS,NBAS,
     +                      NBAS,NBAS,
     +                      P )
     +




      RETURN
      END


