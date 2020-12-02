      SUBROUTINE PHBARXC(ICORE,MAXCOR,IUHF,ISIDE,IRREPX,
     &   SS, SD, DS, DD, ROOT, IPOWER)
C
C THIS ROUTINE DRIVES THE CALCULATION OF THE MATRIX-VECTOR
C PRODUCT BETWEEN A SINGLES AND DOUBLES HAMILTONIAN MATRIX 
C AND A VECTOR USING EXPLICITLY A PARTITIONED EOM SCHEME
C
C CALCULATE SINGLES TO SINGLES
C CALCULATE SINGLES TO DOUBLES PART
C DIVIDE DOUBLES BY - DIAGONAL IN DD-BLOCK
C COPY THESE DOUBLES IN INPUT VECTOR
C MULTIPLY TO GET DOUBLES IN SINGLES PART
C EFFECTIVELY WE CALCULATE [ A +  C (ROOT - D)^-IPOWER B ] SS
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL SS, SD, DS, DD
      DOUBLE PRECISION TIN, TOUT, TIMDUM, TDAVID, TMULT
      DOUBLE PRECISION ROOT
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTIJKA(2),LISTAIBJ(2),LISTABCI(2)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA LISTABCI/127, 27/
      DATA LISTIJKA/107,  7/
      DATA LISTAIBJ/ 54, 54/
      DATA LISTT1,LISTT2,LISTT2IN,LISTT2RS  /490,444,461,440/
      DATA LSTT2RNG /434/
C
      IF (DD) THEN
        WRITE(6,*) ' PHBARXC  SHOULDN"T BE CALLED IF FULL',
     &     ' DOUBLE-DOUBLE BLOCK IS TO BE INCLUDED'
        CALL ERREX
      ENDIF
      IF (SS .AND. IPOWER .NE. 1) THEN
        WRITE(6,*) ' IN PHBARXC BOTH SS AND IPOWER .NE. 1', IPOWER
      ENDIF
C
      CALL TIMER(1)
C
C  ZERO OUT INCREMENT LISTS
C
      IF (IPOWER .EQ. 1 .OR. SD) THEN
        DO 10 ISPIN=3,3-2*IUHF,-1
          CALL ZEROLIST(ICORE,MAXCOR,LISTT2IN-1+ISPIN)
   10   CONTINUE
      ENDIF
C
      CALL ZERO(ICORE, NSIZEC)
      DO ISPIN = 1, 1+IUHF
        CALL PUTLST(ICORE, 1,1,1,2+ISPIN, 490)
      ENDDO
C
C CALCULATE H x C1 -> C1
C
      IF (SS) THEN
C
        CALL DT1INT1(ICORE,MAXCOR,IUHF,IRREPX,IUHF.EQ.0,LISTAIBJ(1),
     &     LISTT1,ISIDE)
        CALL DFT1INT1(ICORE,MAXCOR,IUHF,IRREPX,LISTT1,ISIDE)
C
      ENDIF
C
C CALCULATE H x C1 -> C2
C
        IF (SD) THEN
C
          CALL DT1INT2A(ICORE,MAXCOR,IUHF,IRREPX,LISTIJKA(ISIDE),
     &       LISTT2IN,LISTT1)
          CALL DT1INT2B(ICORE,MAXCOR,IUHF,IRREPX,LISTABCI(ISIDE),
     &       LISTT2IN,LISTT1)
          CALL RSC1INC2(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,LISTT2IN,
     &                   LISTT2RS)
C
        ENDIF
C
C  DIVIDE INCREMENT DOUBLES BY DENOMINATOR AND PUT ON LISTT2
C
      CALL HDIVIDEDD(ICORE, MAXCOR, IUHF, IRREPX, NSIZEC, ROOT,
     &     IPOWER)
C
        IF (DS) THEN
C
C CALCULATE H x C2 -> C1
C
          CALL DT2INT1A(ICORE,MAXCOR,IUHF,IRREPX,LISTABCI(3-ISIDE),
     &       LISTT2,LISTT1)
          CALL DT2INT1B(ICORE,MAXCOR,IUHF,IRREPX,LISTIJKA(3-ISIDE),
     &       LISTT2,LISTT1)
          IF (ISIDE.EQ.1) THEN
            CALL RESORT(ICORE,MAXCOR,IUHF,IRREPX,LISTT2,LSTT2RNG)
            CALL DFT2INT1(ICORE,MAXCOR,IUHF,1,IRREPX)
          ELSEIF (ISIDE .EQ. 2) THEN
            CALL RSC2INC1(ICORE,MAXCOR,IUHF,IRREPX,ISIDE)
          ENDIF
C
        ENDIF
C
      CALL TIMER(1)      
      TMULT = TMULT + TIMENEW
C
      RETURN
      END