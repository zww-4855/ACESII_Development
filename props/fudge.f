
      SUBROUTINE FUDGE(DENS,NORB,NOCC,IUHF)
C
C CONSTRUCTS THE SCF DENSITY IN THE MO BASIS.  REAL TOUGH STUFF.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DENS(NORB,NORB)
      Z=1.0
      IF(IUHF.EQ.0)Z=2.0
      CALL ZERO(DENS,NORB*NORB)
      DO 10 I=1,NOCC
       DENS(I,I)=Z
10    CONTINUE
      RETURN
      END