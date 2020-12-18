         SUBROUTINE  MAT__C_EQ_ORTHOTRAN_DIAG
     +
     +                    ( DDROWA,DDCOLA,
     +                      DDROWB,DDCOLB,
     +                      DDROWC,DDCOLC,
     +                      DDWORK,
     +                      ROW,SUM,
     +                      ROWOFF,
     +                      SAVEB,LTRG,UTRG,
     +                      A,B,
     +
     +                              WORK,
     +                              C )
     +
C------------------------------------------------------------------------
C  OPERATION   : MAT__C_EQ_ORTHOTRAN_DIAG
C  MODULE      : Matrix
C  MODULE-ID   : MAT
C  DESCRIPTION : This routine generates a diagonal block of the
C                orthogonally transformed matrix B:
C
C                              C = A(T) * B * A
C
C                for which we can write for each element:
C
C                          C(ij)  =  sum  A(ki) * B(kl) * A(lj)
C                                     kl
C
C                where the index ranges are:
C
C                           i index    ->   1 to ROW
C                           j index    ->   1 to ROW
C                         k=l indices  ->   1 to SUM
C
C
C                There are two possible paths the routine can take:
C
C
C                   i) If SAVEB = .false. is passed, the original B
C                      matrix can be overwritten. Note, that this is
C                      only possible, if DDCOLB >= ROW, because only
C                      then can the intermediate product B*A be copied
C                      to B. This option is provided at the moment
C                      for those cases in which the full product
C                      result is wanted and B can be wasted, as then
C                      we can invoke efficient full matrix multiply
C                      routines. Hence if at all possible, especially
C                      for large ROW and SUM values, try to use this
C                      case. Note also that here one does not need
C                      the working array WORK.
C
C                  ii) If SAVEB = .true. is passed, we have to save
C                      the original B matrix. This path is chosen
C                      whenever B contains data that should not be
C                      destroyed or any subset of the result is
C                      wanted (lower/upper triangle or only diagonals).
C
C
C                There is the option of placing the resulting square
C                matrix inside the transmitted C array. The final
C                transformed ROW x ROW matrix can be placed starting at
C                a specified row index with offset ROWOFF. Pictorially
C                this looks as follows inside the C array:
C
C
C                                   1          ROW
C
C                                    -----------
C                                   |           |
C                                   |           |
C                                   |           |
C                    IOFF + 1   ->  |-----------|
C                                   |/ / / / / /|
C                                   |/ / / / / /|
C                                   |/ C block /|
C                                   |/ / / / / /|
C                                   |/ / / / / /|
C                  ROWOFF + ROW ->  |-----------|
C                                   |           |
C                                   |           |
C                                   |           |
C                                    -----------
C
C
C                If for some reason one wants the C block to start at
C                the first row, simply pass ROWOFF = 0 in the argument.
C                Since the C matrix has to accomodate the intermediate
C                B*A result and the final ROW x ROW result at the
C                appropriate row index, the dimension requirements
C                for the C matrix are:
C
C                           MAX0 (SUM,ROW+ROWOFF) x ROW
C
C                If LTRG / UTRG are set true or false, the following
C                actions are being taken:
C
C
C                        LTRG | UTRG |    return on exit
C                      ------------------------------------------
C                          T  |   T  |        full matrix
C                          T  |   F  |  lower triangle of matrix
C                          F  |   T  |  upper triangle of matrix
C                          F  |   F  |      diagonal of matrix
C
C
C                In case only the diagonal part of the result is
C                wanted, the ROW elements will be placed directly
C                into the working array, hence WORK will contain
C                the desired result.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     SAVEB,LTRG,UTRG

         INTEGER     DDROWA,DDCOLA
         INTEGER     DDROWB,DDCOLB
         INTEGER     DDROWC,DDCOLC
         INTEGER     DDWORK
         INTEGER     I,J,K
         INTEGER     ROW,SUM
         INTEGER     ROWOFF

         DOUBLE PRECISION  X,ZERO

         DOUBLE PRECISION  WORK (1:DDWORK)

         DOUBLE PRECISION  A (1:DDROWA,1:DDCOLA)
         DOUBLE PRECISION  B (1:DDROWB,1:DDCOLB)
         DOUBLE PRECISION  C (1:DDROWC,1:DDCOLC)
         DOUBLE PRECISION  temp (ROW,ROW)

         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------

C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  'Density matrix before diag',
C     +                  DDROWB,DDCOLB,
C     +                  DDROWB,DDCOLB,
C     +                  B )
C     +


C
C
C             ...check dimensions for matrix B.
C
C
         IF  (SUM.GT.DDROWB .OR. SUM.GT.DDCOLB) THEN
              WRITE (*,*) ' Dimensions of matrix B too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_diag '
              WRITE (*,*) ' DDROWB,DDCOLB,SUM = ',DDROWB,DDCOLB,SUM
              WRITE (1,*) ' Dimensions of matrix B too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_diag '
              WRITE (1,*) ' DDROWB,DDCOLB,SUM = ',DDROWB,DDCOLB,SUM
              STOP
         END IF

         IF  (.NOT.SAVEB .AND. LTRG .AND. UTRG) THEN
              IF (ROW.GT.DDCOLB) THEN
                  WRITE (*,*) ' Cannot store B*A product into B! '
                  WRITE (*,*) ' mat__c_eq_orthotran_diag '
                  WRITE (*,*) ' DDCOLB,ROW = ',DDCOLB,ROW
                  WRITE (1,*) ' Cannot store B*A product into B! '
                  WRITE (1,*) ' mat__c_eq_orthotran_diag '
                  WRITE (1,*) ' DDCOLB,ROW = ',DDCOLB,ROW
                  STOP
              END IF
         END IF
C
C
C             ...check dimensions for matrix A.
C
C
         IF  (SUM.GT.DDROWA .OR. ROW.GT.DDCOLA) THEN
              WRITE (*,*) ' Dimensions of matrix A too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_diag '
              WRITE (*,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                      DDROWA,DDCOLA,ROW,SUM
              WRITE (1,*) ' Dimensions of matrix A too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_diag '
              WRITE (1,*) ' DDROWA,DDCOLA,ROW,SUM = ',
     +                      DDROWA,DDCOLA,ROW,SUM
              STOP
         END IF
C
C
C             ...check dimensions for matrix C.
C
C
         IF  (MAX0 (SUM,ROW+ROWOFF).GT.DDROWC .OR. ROW.GT.DDCOLC) THEN
              WRITE (*,*) ' Dimensions of matrix C too small! '
              WRITE (*,*) ' mat__c_eq_orthotran_diag '
              WRITE (*,*) ' DDROWC,DDCOLC,ROW,ROWOFF,SUM = ',
     +                      DDROWC,DDCOLC,ROW,ROWOFF,SUM
              WRITE (1,*) ' Dimensions of matrix C too small! '
              WRITE (1,*) ' mat__c_eq_orthotran_diag '
              WRITE (1,*) ' DDROWC,DDCOLC,ROW,ROWOFF,SUM = ',
     +                      DDROWC,DDCOLC,ROW,ROWOFF,SUM
              STOP
         END IF
C
C
C             ...check dimensions for working array WORK (if needed).
C
C
         IF  (SAVEB .OR. (.NOT.LTRG) .OR. (.NOT.UTRG)) THEN
              IF  (ROW.GT.DDWORK) THEN
                   WRITE (*,*) ' Dimensions of array WORK too small! '
                   WRITE (*,*) ' mat__c_eq_orthotran_diag '
                   WRITE (*,*) ' DDWORK,ROW = ',DDWORK,ROW
                   WRITE (1,*) ' Dimensions of array WORK too small! '
                   WRITE (1,*) ' mat__c_eq_orthotran_diag '
                   WRITE (1,*) ' DDWORK,ROW = ',DDWORK,ROW
                   STOP
              END IF
         END IF
C
C
C             ...perform intermediate multiplication B * A for all
C                cases and store result in C.
C
C
         CALL    MAT__C_EQ_A_TIMES_B_FLOAT
     +
     +                ( DDROWB,DDCOLB,
     +                  DDROWA,DDCOLA,
     +                  DDROWC,DDCOLC,
     +                  SUM,SUM,ROW,
     +                  B,A,
     +
     +                          C )
     +
     +
         IF (LTRG.AND.UTRG) THEN

             IF (.NOT.SAVEB) THEN
C
C
C             ...copy intermediate result B*A sitting in C to B
C                and calculate full product A(T)*B*A.
C
C
                 CALL    MAT__C_EQ_A_FLOAT
     +
     +                        ( DDROWC,DDCOLC,
     +                          DDROWB,DDCOLB,
     +                          SUM,ROW,
     +                          C,
     +
     +                                  B )
     +
     +
                 CALL    MAT__C_EQ_AT_TIMES_B_FLOAT
     +
     +                        ( DDROWA,DDCOLA,
     +                          DDROWB,DDCOLB,
     +                          DDROWC,DDCOLC,
     +                          ROW,SUM,ROW,
     +                          A,B,
     +
     +                                  C )
     +
     +
                 IF (ROWOFF.NE.0) THEN
                     DO 100 J = 1,ROW
                     DO 100 I = 1,ROW
                        C (ROWOFF+I,J) = C (I,J)
  100                CONTINUE
                 END IF
             ELSE
C
C
C             ...calculate full product A(T)*B*A directly into C.
C
C
                 DO 200 J = 1,ROW
                    DO 210 I = 1,ROW
                       X = ZERO
                       DO 220 K = 1,SUM
                          X = X + A (K,I) * C (K,J)
  220                  CONTINUE
                       WORK (I) = X
  210               CONTINUE
                    DO 230 I = 1,ROW
                       C (ROWOFF+I,J) = WORK (I)
  230               CONTINUE
  200            CONTINUE
             END IF

         ELSE IF (LTRG) THEN
C
C
C             ...calculate lower triangle of product A(T)*B*A directly
C                into C.
C
C
             DO 300 J = 1,ROW
                DO 310 I = J,ROW
                   X = ZERO
                   DO 320 K = 1,SUM
                      X = X + A (K,I) * C (K,J)
  320              CONTINUE
                   WORK (I) = X
  310           CONTINUE
                DO 330 I = J,ROW
                   C (ROWOFF+I,J) = WORK (I)
  330           CONTINUE
  300        CONTINUE

         ELSE IF (UTRG) THEN
C
C
C             ...calculate upper triangle of product A(T)*B*A directly
C                into C.
C
C
             DO 400 J = 1,ROW
                DO 410 I = 1,J
                   X = ZERO
                   DO 420 K = 1,SUM
                      X = X + A (K,I) * C (K,J)
  420              CONTINUE
                   WORK (I) = X
  410           CONTINUE
                DO 430 I = 1,J
                   C (ROWOFF+I,J) = WORK (I)
  430           CONTINUE
  400        CONTINUE

         ELSE
C
C
C             ...calculate diagonal of product A(T)*B*A and place
C                result into WORK.
C
C
             DO 500 I = 1,ROW
                X = ZERO
                DO 510 K = 1,SUM
                   X = X + A (K,I) * C (K,I)
  510           CONTINUE
                WORK (I) = X
C                write(1,*) WORK(I)
  500        CONTINUE

         END IF

C      temp=ZERO
C      DO I=1,ROW
C       DO J=1,ROW
C        DO K=1,SUM
C         temp(I,J)=temp(I,J)+A(K,I)*C(K,J)
C        END DO
C       END DO
C      END DO

C         CALL    MAT__PRINT_A_FLOAT_5_NOZEROS
C     +
C     +                ( 1,
C     +                  'Density matrix after diag',
C     +                  ROW,ROW,
C     +                  ROW,ROW,
C     +                  temp )
C     +

C
C
C             ...ready!
C
C
         RETURN
         END
