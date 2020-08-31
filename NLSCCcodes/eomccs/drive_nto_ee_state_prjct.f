




































































































































































































       SUBROUTINE DRIVE_NTO_EE_STATE_PRJCT(SCR,MAXCOR,IRREPX,IUHF,IROOT)
C
C The core EEs are obtained by setting the R(ij,ab) = 0 when i and
C j do  not belong to the core region. The present version works
C state specifc fashion.
C
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C 
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR)
      DOUBLE PRECISION VALUE,O_SCALE,V_SCALE 
      DOUBLE PRECISION VAL1,VAL2,DIFF,THRES2 
      INTEGER INTOO(2),INTOV(2)
C
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
C
      LOGICAL PROJECT_SINGLES, NTOS_EXIST
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      CHARACTER*1 BLANK,STRING3
      CHARACTER*4 STRING
      CHARACTER*31 STRING1 
      CHARACTER*26 STRING2

      DIMENSION LS2OUT(2,2), IMAP_A(MAXBASFN),IMAP_B(MAXBASFN) 
      DIMENSION OCA_OFF(8),OCB_OFF(8)
      DIMENSION VRA_OFF(8),VRB_OFF(8)
      DIMENSION BGN_ORB(MAXBASFN,8),END_ORB(MAXBASFN,8)
      DIMENSION O_SCALE(MAXBASFN,8),V_SCALE(MAXBASFN,8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /PROJECT/ IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      DATA BLANK /" "/

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
C
      PROJECT_SINGLES = (Iflags2(173) .NE. 0)

      LUNITN = 30
      INQUIRE(FILE='NTO_EXCIT',EXIST=NTOS_EXIST)
      IF (NTOS_EXIST) THEN
       OPEN(LUNITN,FILE='NTO_EXCIT',STATUS='OLD',FORM='FORMATTED',
     &     POSITION='REWIND')
C          REWIND(LUNITN)
C
       READ(LUNITN,'(A,1X,I1)') STRING1, IRREPC
C      READ(LUNITN,'(A)') BLANK
       IRREPX = IRREPC
       ISPIN  = 1

       VAL1 = 0.0D0 
       DO IRREP = 1, NIRREP
          IF (POP(IRREP,ISPIN).GE.1) THEN
             READ(LUNITN,'(A)') BLANK
             READ(LUNITN,'(A,1x,I1)') STRING2,IRREPC
             READ(LUNITN,'(A)') BLANK
             IF (STRING2(13:16) .EQ. "occ") IRREPI = IRREPC
       
             DO I = 1, POP(IRREP,ISPIN)
                READ(LUNITN,'(A,I7,1X,F15.10)') STRING,INDEX,VALUE
                IF (STRING .EQ. "OCC") THEN
                    BGN_ORB(I,IRREP) = INDEX
                    O_SCALE(I,IRREP) = VALUE
                    IF (DABS(VALUE).GE.VAL1) THEN
                       VAL1=DABS(VALUE)
                       INTOO(1)=IRREP
                       INTOO(2)=I
                    ENDIF 
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       WRITE(6,'(A,I1,A,I3)') 'Largest excitation (OCC): Irrep.',
     &      INTOO(1),',  No.',INTOO(2)
       
       VAL2 = 0.0D0 
       DO IRREP = 1, NIRREP
          IF (VRT(IRREP,ISPIN).GE.1) THEN
             READ(LUNITN,'(A)') BLANK
             READ(LUNITN,'(A,1x,I1)') STRING2,IRREPC
             READ(LUNITN,'(A)') BLANK
             IF (STRING(13:16) .EQ. "vrt") IRREPA = IRREPC
       
             DO I = 1, VRT(IRREP,ISPIN)
                READ(LUNITN,'(A,I7,1X,F15.10)') STRING,INDEX,VALUE
                IF (STRING .EQ. "VRT") THEN
                   END_ORB(I,IRREP) = INDEX
                   V_SCALE(I,IRREP) = VALUE
                   IF (DABS(VALUE).GE.VAL2) THEN
                       VAL2=DABS(VALUE)
                       INTOV(1)=IRREP
                       INTOV(2)=I
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO 
       WRITE(6,'(A,I1,A,I3)') 'Largest excitation (OCC): Irrep.',
     &      INTOV(1),',  No.',INTOV(2)
      ENDIF

      IF (NTOS_EXIST) THEN
       Noccs = nocco(1)
       Nvrts = nvrto(1)
       Write(6,"(a,15(1x,I2))") "Begin: ", ((BGN_ORB(i,j),i=1,Noccs),
     &                                                j=1,Nirrep)
       Write(6,*)
       Write(6,"(a,15(1x,I2))") "End  : ", ((End_ORB(i,j),i=1,Nvrts),
     &                                                j=1,Nirrep)
      ENDIF
C
      LS1OUT      = 490
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      STRINGI(2)='i  [i_SYM]  '
      STRINGJ(2)='j  [j_SYM]  '
      STRINGA(2)='a  [a_SYM]  '
      STRINGB(2)='b  [b_SYM]  '
      STRINGI(1)='I  [I_SYM]  '
      STRINGJ(1)='J  [J_SYM]  '
      STRINGA(1)='A  [A_SYM]  '
      STRINGB(1)='B  [B_SYM]  '
     
      NBAS   = NOCCO(1) + NVRTO(1)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_A)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_B)
      IF (IUHF.NE.0) CALL GETREC(20,"JOBARC","ORBMAP_B",NBAS,
     &                           IMAP_B)
      
      OCA_OFF(1) = 0
      OCB_OFF(1) = 0

      DO I = 1, (NIRREP-1)
         OCA_OFF(I+1) = OCA_OFF(I) + POP(I,1)
         OCB_OFF(I+1) = OCB_OFF(I) + POP(I,2)
      ENDDO

      VRA_OFF(1) = NOCCO(1)
      VRA_OFF(1) = NOCCO(2)

      DO I = 1, (NIRREP-1)
         VRA_OFF(I+1) = VRA_OFF(I) + VRT(I,1)
         VRB_OFF(I+1) = VRB_OFF(I) + VRT(I,2)
      ENDDO

      CALL ZERO(SCR,NSIZEC)
c
C  FIRST CONSIDER SINGLE EXCITATION COEFFICIENTS
C
      IF (PROJECT_SINGLES) THEN

      DO 500 SSPIN = 1, 1+IUHF
        CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
        ICOUNT = 1
        DO IIRREP = 1, NIRREP
          AIRREP=DIRPRD(IIRREP,IRREPX)
C         WRITE(6,'(A,3(X,I1))') 'IIRREP*AIRREP=IRREPX',
C    &         IIRREP,AIRREP,IRREPX
          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)
              IF (SSPIN .EQ. 1) THEN
                 II = IMAP_A(I+OCA_OFF(IIRREP))
CYCP             AA = IMAP_A(A+VRA_OFF(IIRREP))
                 AA = IMAP_A(A+VRA_OFF(AIRREP))
              ELSE
CYCP             II = IMAP_B(I+OCB_OFF(AIRREP))
CYCP             AA = IMAP_A(A+VRB_OFF(IIRREP))
                 II = IMAP_B(I+OCB_OFF(IIRREP))
                 AA = IMAP_B(A+VRB_OFF(AIRREP))
              ENDIF

              IF (NTOS_EXIST) THEN

                IF (IIRREP.EQ.INTOO(1) .AND. AIRREP.EQ.INTOV(1)) THEN
                   IF (I.EQ.INTOO(2) .AND. A.EQ.INTOV(2)) THEN
                      WRITE(6,'(2(A,I3,A,I1,A))')'R1:',I,'[',IIRREP,']',
     &                                           ' ->',A,'[',AIRREP,']'
                      SCR(ICOUNT) = 1.0D0
                   ELSE 
                      SCR(ICOUNT) = 0.0D0
                   ENDIF 
                ELSE 
                   SCR(ICOUNT) = 0.0D0
                ENDIF 
              ENDIF

CSSS                 IF (DABS(O_SCALE(I,IIRREP)) .GT. THRES .AND.
CSSS     +               DABS(V_SCALE(A,AIRREP)) .GT. THRES) THEN
CSSS                     IF (IORB_A(IROOT) .EQ. II) THEN

C               IF (I .EQ. 1 .AND. A .EQ. 1) THEN
C                   WRITE(6,*) 'IIRREP/AIRREP',IIRREP,AIRREP
C                   SCR(ICOUNT) = 1.0D0
C               ELSE 
C                   SCR(ICOUNT) = 0.0D0
C               ENDIF 

              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO

        CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
  500 CONTINUE
C
      ENDIF 

C
      DO ICASE = 1, 1 + 2 * IUHF
        IF (ICASE .EQ. 1) THEN
          ASPIN = 1
          BSPIN = 2
        ELSEIF (ICASE .EQ. 2) THEN
          ASPIN = 1
          BSPIN = 1
        ELSEIF (ICASE .EQ. 3) THEN
          ASPIN = 2
          BSPIN = 2
        ENDIF

      LISTS2EX = LS2OUT(ASPIN, BSPIN)

      DO 10 RIRREP = 1, NIRREP
        LIRREP = DIRPRD(RIRREP,IRREPX)
        IF (ICASE .EQ. 1) THEN
          DISSYS = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSS = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYA = DISSYS
          NUMDSA = NUMDSS
        ELSE
          DISSYA = IRPDPD(LIRREP, ISYTYP(1,LISTS2EX))
          NUMDSA = IRPDPD(RIRREP, ISYTYP(2,LISTS2EX))
          DISSYS = IRPDPD(LIRREP, 18 + ASPIN)
          NUMDSS = IRPDPD(RIRREP, 20 + ASPIN)
        ENDIF

       IF (DISSYA * NUMDSA .GT. 0) THEN
          I000 = 1
          I010 = I000 + DISSYS*NUMDSS *IINTFP
          IF (ICASE .EQ. 1) THEN
            CALL GETLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $       LISTS2EX)
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ELSE
            CALL GETLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $       LISTS2EX)
            CALL SYMEXP(RIRREP, POP(1,ASPIN),DISSYA,SCR(I000))
            CALL SYMEXP2(LIRREP,VRT(1,ASPIN),DISSYS, DISSYA,
     &         NUMDSS, SCR(I000), SCR(I000))
            CALL ZERO(SCR(I000),DISSYS*NUMDSS)
          ENDIF

          ICOUNT = I000
          DO 20 JIRREP = 1, NIRREP
            IIRREP = DIRPRD(JIRREP, RIRREP)
            DO 30 J= 1, POP(JIRREP,BSPIN)
              DO 40 I = 1, POP(IIRREP, ASPIN)
                DO 50 BIRREP = 1, NIRREP
                  AIRREP = DIRPRD(LIRREP, BIRREP)
                  DO 60 B = 1, VRT(BIRREP, BSPIN)
                    DO 70 A = 1, VRT(AIRREP, ASPIN)

                       IF (ICASE .EQ. 1) THEN

                          SCR(ICOUNT) = 1.0D0
                          II = IMAP_A(I+OCA_OFF(IIRREP))
                          JJ = IMAP_A(J+OCB_OFF(JIRREP))

                       ELSE IF (ICASE .EQ. 2) THEN 

                          SCR(ICOUNT) = 1.0D0

                          II = IMAP_A(I+OCA_OFF(IIRREP))
                          JJ = IMAP_A(J+OCA_OFF(JIRREP))

                       ELSE IF (ICASE .EQ. 3) THEN 

                          SCR(ICOUNT) = 1.0D0
                          II = IMAP_B(I+OCB_OFF(IIRREP))
                          JJ = IMAP_B(J+OCB_OFF(JIRREP))

                       ENDIF 

                    ICOUNT = ICOUNT + 1
   70             CONTINUE
   60           CONTINUE
   50         CONTINUE
   40       CONTINUE
   30     CONTINUE
   20   CONTINUE

      ENDIF

      IF (ICASE .EQ. 1) THEN
         CALL PUTLST(SCR(I000), 1, NUMDSS, 1, RIRREP,
     $               LISTS2EX)
      ELSE
C
C SQUEEZE ARRAY IN PROPER FORM A<B, I<J
C
         CALL SQSYM(LIRREP, VRT(1,ASPIN),DISSYA,DISSYS,
     &              NUMDSS,SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),NUMDSS,DISSYA)
         CALL SQSYM(RIRREP,POP(1,ASPIN), NUMDSA, NUMDSS,
     &              DISSYA, SCR(I010), SCR(I000))
         CALL TRANSP(SCR(I010), SCR(I000),DISSYA,NUMDSA)
         CALL PUTLST(SCR(I000), 1, NUMDSA, 1, RIRREP,
     $               LISTS2EX)
      ENDIF

   10 CONTINUE
      ENDDO
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
C
C  PUT NORMALIZED PROJECTION VECTORS BACK ON LIST
C     
      CALL UPDATES(IRREPX,SCR,444,0,490,IUHF)
C
      RETURN
      END
                






