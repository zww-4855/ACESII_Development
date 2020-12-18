         SUBROUTINE  NLO__NORMALIZE_OVERLAP_DENSITY
     +
     +                    ( NBAS,
     +                      X,
     +                      LABEL,
     +
     +                             DENORM,
     +                             D,S )
     +
C------------------------------------------------------------------------
C  OPERATION   : NLO__NORMALIZE_OVERLAP_DENSITY
C  MODULE      : Natural Localized Orbitals
C  MODULE-ID   : NLO
C  DESCRIPTION : This routine normalizes the overlap matrix and also
C                updates the density matrix. This is necessary to get
C                the proper weighting in the averaging of the density
C                matrix blocks over the m-quantum number.
C
C                The density matrix elements D (i,j) are defined as the
C                expansion coefficients of the one electron density
C                rho (r,r) in terms of atomic orbital basis products:
C
C
C                          rho (r,r) = sum  D (i,j) |i> <j|
C                                       ij
C
C                Lets assume the AO basis is not normalized, that is:
C
C                             S (i,i) = <i|i> neq 1
C                             S (j,j) = <j|j> neq 1
C
C                We obtain new normalized AO's via:
C
C                                 ___________
C                        |i'> = \/ (1/S(i,i)) * |i>
C                                 ___________
C                        |j'> = \/ (1/S(j,j)) * |j>
C
C                which we can use to write a new rho (r,r) expansion:
C
C
C                          rho (r,r) = sum  D'(i,j) |i'> <j'|
C                                       ij
C
C                from which we deduce that:
C
C                                   _______     _______
C                       D'(i,j) = \/ S(i,i) * \/ S(j,j) * D(i,j)
C
C
C                The normalized matrix elements are then given by:
C
C                                ___________     ___________
C                    S'(i,j) = \/ (1/S(i,i)) * \/ (1/S(j,j)) * S(i,j)
C
C
C                Since the square roots of the diagonal overlap matrix
C                elements will be needed in order to express the
C                natural orbitals in the original unnormalized AO basis,
C                they are stored in a special array and returned.
C
C
C
C                Input:
C
C                     NBAS = size of AO basis.
C                        X = storage array of floating type.
C                    LABEL = storage array of integer type of size
C                            NBAS.
C                      D,S = original density and overlap matrix.
C
C                Output:
C
C                   DENORM = denormalization factors.
C                      D,S = normalized density and overlap matrix.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    I,J,K,N
         INTEGER    NBAS

         INTEGER    LABEL (1:NBAS)

         DOUBLE PRECISION  SDIAG
         DOUBLE PRECISION  SSQR
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  TINY
         DOUBLE PRECISION  XK,YK

         DOUBLE PRECISION  DENORM (1:NBAS)
         DOUBLE PRECISION  X      (1:NBAS)

         DOUBLE PRECISION  D (1:NBAS,1:NBAS)
         DOUBLE PRECISION  S (1:NBAS,1:NBAS)

         PARAMETER    (ONE   =  1.D0 )
         PARAMETER    (TINY  =  1.D-10)
C
C
C------------------------------------------------------------------------
C
C
C             ...search through the diagonal elements of the
C                overlap matrix S and identify those which are
C                different from 1 within a tolerance limit.
C                Save the square root and the inverse of the
C                square root for further use.
C
C

C----------------------------------------------------------
C Print initial overlap matrix
C         write(*,*) ((S(I,J),J=1,NBAS),I=1,NBAS)
C         write(*,*) (S(I,I),I=1,NBAS)
C----------------------------------------------------------

C----------------------------------------------------------------
C    Print original density matrix
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Density matrix before normalization',
C     +                  NBAS,NBAS,
C     +                  NBAS,NBAS,
C     +                  D )
C     +
C     +
C----------------------------------------------------------------

         N = 0
         DO 10 I = 1,NBAS
            SDIAG = S (I,I)
            IF ( DABS (SDIAG - ONE).GT.TINY ) THEN
                 SSQR = DSQRT (SDIAG)
                 N = N + 1
                 X (N) = SSQR
                 DENORM (N) = ONE / SSQR
                 LABEL (N) = I
            END IF
   10    CONTINUE
C
C
C             ...update the overlap and density matrix with inner
C                loops over rows to minimize cache misses.
C
C
         IF (N.NE.0) THEN
             DO 20 J = 1,NBAS
             DO 20 K = 1,N
                I = LABEL (K)
                D (I,J) =      X (K) * D (I,J)
                S (I,J) = DENORM (K) * S (I,J)
   20        CONTINUE

C----------------------------------------------------------
C Print intermediate overlap matrix
C         write(*,*) ((S(I,J),J=1,NBAS),I=1,NBAS)
C         write(*,*) "intermediate overlap matrix"
C         write(*,*) (S(I,I),I=1,NBAS)
C----------------------------------------------------------

             DO 30 K = 1,N
                J = LABEL (K)
                XK = X (K)
                YK = DENORM (K)
                DO 32 I = 1,NBAS
                   D (I,J) = XK * D (I,J)
                   S (I,J) = YK * S (I,J)
   32           CONTINUE
   30        CONTINUE
         END IF

C----------------------------------------------------------
C Print final overlap matrix
C         write(*,*) ((S(I,J),J=1,NBAS),I=1,NBAS)
C         write(*,*) "final overlap matrix"
C         write(*,*) (S(I,I),I=1,NBAS)
C----------------------------------------------------------

C
C
C             ...generate the denormalization array.
C
C
         DO 40 I = 1,NBAS
            DENORM (I) = ONE
   40    CONTINUE

         DO 50 K = 1,N
            I = LABEL (K)
            DENORM (I) = X (K)
   50    CONTINUE
C
C
C             ...print for checking.
C
C
C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  ' Overlap matrix after normalization',
C     +                  NBAS,NBAS,
C     +                  NBAS,NBAS,
C     +                  S )
C     +
C     +
         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
     +
     +                ( 1,
     +                  ' Density matrix after normalization',
     +                  NBAS,NBAS,
     +                  NBAS,NBAS,
     +                  D )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
