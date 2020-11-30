










      SUBROUTINE HBARXC(ICORE,MAXCOR,IUHF,ISIDE,IRREPX)
C
C THIS ROUTINE DRIVES THE CALCULATION OF THE MATRIX-VECTOR
C PRODUCT BETWEEN A SINGLES AND DOUBLES HAMILTONIAN MATRIX 
C AND A VECTOR.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL DOUBLE, SS, SD, DS, DD
      DOUBLE PRECISION TIN, TOUT, TIMDUM, TDAVID, TMULT
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTIJKL(2),LISTIJKA(2),LISTAIBJ(2),LISTABCI(2)
      DIMENSION LISTABCD(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
      DATA LISTABCD/231,231/
      DATA LISTABCI/127, 27/
      DATA LISTIJKA/107,  7/
      DATA LISTIJKL/ 51, 51/
      DATA LISTAIBJ/ 54, 54/
      DATA LISTT1,LISTT2,LISTT2IN,LISTT2RS  /490,444,461,440/
      DATA LSTT2RNG /434/


      CALL TIMER(1)
C
C INITIALIZE T2 INCREMENT LIST
C
      DO 10 ISPIN=3,3-2*IUHF,-1
        CALL ZEROLIST(ICORE,MAXCOR,LISTT2IN-1+ISPIN)
 10   CONTINUE
C
C CALCULATE H x C1 -> C1
C
      IF (SS) THEN

        If (Ispar) Then 
           Call RESTORE_CC_WMBEJ(ICORE,MAXCOR,IUHF)
           if (iuhf .eq. 0) CALL MAKESS(ICORE,MAXCOR,IUHF)
        Endif 

        CALL DT1INT1(ICORE,MAXCOR,IUHF,IRREPX,IUHF.EQ.0,LISTAIBJ(1),
     &     LISTT1,ISIDE)

        If (Ispar) Then 
           Call RESTORE_PDCC_WMBEJ(ICORE,MAXCOR,IUHF)
           if (iuhf .eq. 0) CALL MAKESS(ICORE,MAXCOR,IUHF)
        Endif 

        CALL DFT1INT1(ICORE,MAXCOR,IUHF,IRREPX,LISTT1,ISIDE)
C
      ENDIF
C
      double=.true.

      IF (DOUBLE) THEN
C
C CALCULATE H x C2 -> C1
C
        IF (DS) THEN

          CALL DT2INT1A(ICORE,MAXCOR,IUHF,IRREPX,LISTABCI(3-ISIDE),
     &       LISTT2,LISTT1)
          CALL DT2INT1B(ICORE,MAXCOR,IUHF,IRREPX,LISTIJKA(3-ISIDE),
     &       LISTT2,LISTT1)
         CALL RESORT(ICORE,MAXCOR,IUHF,IRREPX,LISTT2,LSTT2RNG)
         IF (ISIDE.EQ.1) THEN
            CALL DFT2INT1(ICORE,MAXCOR,IUHF,1,IRREPX)
         ENDIF
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
C
        ENDIF

C
C CALCULATE H x C2 -> C2
C
        IF (DD) THEN
C
          CALL DT2INT2 (ICORE,MAXCOR,IUHF,IRREPX,ISIDE,LISTIJKL,
     &       LISTABCD,LISTAIBJ,LISTT2,LISTT2IN,LISTT2RS,LSTT2RNG)
C
        ELSE
C
          CALL RSDSNSD (ICORE,MAXCOR,IUHF,IRREPX,ISIDE,LISTT2IN,
     &       LISTT2RS)
          CALL HDIAGXDD (ICORE,MAXCOR/IINTFP,IUHF,IRREPX,NSIZEC)
C 
        ENDIF
C
      ENDIF 
C
      CALL TIMER(1)      
      TMULT = TMULT + TIMENEW
C
      RETURN
      END
