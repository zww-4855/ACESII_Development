      SUBROUTINE GSCHMIDT(VEC,VORTH,NSIZE,NDIM,TMP,RESID)
C
C THIS PROJECTS OUT ALL PARTS OF AN INPUT VECTOR (VEC)
C WHICH LIE IN THE SPACE SPANNED BY THE ORTHOGONAL BASIS
C VORTH.
C
C   |v'> = |v> - SUM <i|v> |i>
C                 i 
C
C WHERE THE |i> ARE NORMALIZED BASIS VECTORS FOR THE SPACE VORTH
C
C INPUT:
C       VEC : THE VECTOR WHICH IS TO BE ORTHOGONALIZED TO
C             THE EXISTING BASIS.  *IT IS ASSUMED THAT VEC
C             IS NORMALIZED ON INPUT)
C     VORTH : THE BASIS VECTORS FOR THE EXISTING ORTHOGONAL
C             BASIS
C     NSIZE : THE LENGTH OF THE BASIS VECTORS
C      NDIM : THE DIMENSION OF THE ORTHOGONAL SPACE
C       TMP : A SCRATCH VECTOR OF LENGTH NSIZE
C     RESID : THE NORM OF VEC, AFTER ORTHONALIZATION AND
C             BEFORE NORMALIZATION
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(NSIZE),VORTH(NSIZE,NDIM),TMP(NSIZE)
      logical print
C
      DATA ONE /1.0D0/
      DATA TOL /1.D-10/
C
      IF(NDIM.EQ.0)THEN
       RESID=SNRM2(NSIZE,VEC,1)
       X = ONE / RESID
       CALL SSCAL(NSIZE,X,VEC,1)
       RETURN
      ENDIF
C
      print = .false.
      if (print) then
        write(6,*) ' old vector in gschmidt'
        call output(vec,1,1,1,nsize,1,nsize,2)
        write(6,*) ' orthogonal vectors in gschmidt'
        call output(vorth,1,nsize,1,ndim,nsize,ndim,1)
        do i = 1, ndim
          TMP(i)=SDOT(NSIZE,VORTH(1,I),1,VEC,1)
        enddo
        write(6,*) ' old overlap in gschmidt'
        call output(tmp,1,1,1,ndim,1,ndim,1)
      endif
C
      CALL SCOPY(NSIZE,VEC,1,TMP,1)
      DO 10 I=NDIM ,1 ,-1
       FACT=  SDOT(NSIZE,VORTH(1,I),1,VEC,1)
       if (print) write(6,*) i, fact
       CALL SAXPY(NSIZE,- FACT,VORTH(1,I),1,TMP,1)
10    CONTINUE
      CALL SCOPY(NSIZE,TMP,1,VEC,1)
C
C RENORMALIZE THE RESIDUAL
C
      RESID=SNRM2(NSIZE,VEC,1)
      IF(RESID.GT.TOL)THEN
       X=ONE/RESID
       CALL SSCAL(NSIZE,X,VEC,1)
      ENDIF 
C
      if (print) then
        write(6,*) ' new vector in gschmidt'
        call output(vec,1,1,1,nsize,1,nsize,1)
        write(6,*) ' orthogonal vectors in gschmidt'
        call output(vorth,1,nsize,1,ndim,nsize,ndim,1)
        do i = 1, ndim
          TMP(i)=SDOT(NSIZE,VORTH(1,I),1,VEC,1)
        enddo
        write(6,*) ' current overlap in gschmidt'
        call output(tmp,1,1,1,ndim,1,ndim,1)
      endif
      RETURN
      END

