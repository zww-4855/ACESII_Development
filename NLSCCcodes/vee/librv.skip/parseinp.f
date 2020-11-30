      SUBROUTINE PARSEINP(NBAS,SCR,ISCRAI,ISCRIJ,ISCRAB,
     &                    IMAPMO,IRREPORB,ILOC)
C
C READS EXCITATION INPUT GIVEN AS DOMINANT SINGLE EXCITATION
C IN GENERAL, THESE ARE IN THE FORM
C
C  2*
C    A 3 15 
C    B 11 103 
C
C WHICH MEANS THAT TWO ROOTS WILL BE SEARCHED FOR - ONE PRINCIPALLY
C DESCRIBED AS AN EXCITATION FROM ALPHA ORBITAL 3 TO ALPHA ORBITAL 15 
C WITH THE BEING BETA ORBITAL 11 TO BETA ORBITAL 103.  THE ORBITAL 
C NUMBERING USED HERE CORRESPONDS TO THE NUMBERING SCHEME OF ORBITALS 
C IN WHICH THE EIGENVALUES ARE SORTED FROM SMALLEST TO LARGEST.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*4 SPTYPE
      DIMENSION SCR(NBAS,*),ISCRAI(NBAS*NBAS,3),IMAPMO(NBAS,*)
      DIMENSION ISCRAB(NBAS*NBAS,3),ISCRIJ(NBAS*NBAS,3)
      DIMENSION IRREPORB(NBAS,*),IOFFIJ(8,3),IOFFAB(8,3)
      DIMENSION NUMBIJ(8,3),NUMBAB(8,3),IOFFABIJ(8,8,3)
      DIMENSION IOFFO(8,2),IOFFV(8,2),ILOC(*)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/CALCINFO/NROOT(8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/GUESS2/IMAP(100,8)
      COMMON/GUESS3/ZMAP2(10,100,8),IMAP2(10,100,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      INDX(I,J,N)=I+(J-1)*N
C
      IUHF=MIN(IFLAGS(11),1)
      CALL IZERO(NROOT,8)
      CALL IZERO(IMAP,800)
      CALL IZERO(IOFFAB,24)
      CALL IZERO(IOFFIJ,24)
      CALL IZERO(IOFFABIJ,8*8*3)
C
      WRITE(6,1000)
1000  FORMAT(T3,'@PARSEINP-I, Input particle-hole excitations used ',
     &          'as initial guesses.')
  
C
C READ IN EIGENVALUES, SORT THEM AND HOLD MAPPING VECTOR RELATING
C SYMMETRY ORDER AND UNPACKED ORDER
C
      DO 100 ISPIN=1,1+IUHF
       DO 101 I=1,NBAS
        IMAPMO(I,ISPIN)=I
101    CONTINUE
100   CONTINUE
      CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*NBAS,SCR)
      CALL GETREC(20,'JOBARC','IRREPALP',NBAS,IRREPORB(1,1))
      CALL PIKSR2(NBAS,SCR,IMAPMO(1,1))
      IF(IUHF.NE.0)THEN
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*NBAS,SCR(1,2))
       CALL GETREC(20,'JOBARC','IRREPBET',NBAS,IRREPORB(1,2))
       CALL PIKSR2(NBAS,SCR(1,2),IMAPMO(1,2))
      ELSE
       CALL ICOPY(NBAS,IRREPORB,1,IRREPORB(1,2),1)
       CALL ICOPY(NBAS,IMAPMO,  1,IMAPMO(1,2),  1)
      ENDIF
C
C CONSTRUCT SYMMETRY VECTOR FOR ALL IRREPS SO THAT WE CAN ULTIMATELY
C GET THE OFFSET INTO THE AI OR ai VECTOR FROM THE INDIVIDUAL VALUES
C OF A AND I 
C
C
      DO 1 ISPIN=1,2
       IOFFO(1,ISPIN)=0
       IOFFV(1,ISPIN)=0
       DO 2 IRREP=1,NIRREP-1
        IOFFO(IRREP+1,ISPIN)=IOFFO(IRREP,ISPIN)+POP(IRREP,ISPIN)
        IOFFV(IRREP+1,ISPIN)=IOFFV(IRREP,ISPIN)+VRT(IRREP,ISPIN)
2      CONTINUE
1     CONTINUE
C
       DO 10 IRREPX=1,NIRREP
        ITHRU=0
        DO 15 ISPIN=1,1+IUHF
         DO 11 IRREPI=1,NIRREP
          IRREPA=DIRPRD(IRREPI,IRREPX)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          DO 12 INDI=1,NUMI
           DO 13 INDA=1,NUMA
            ITHRU=ITHRU+1
            INDABSI=INDI+IOFFO(IRREPI,ISPIN)
            INDABSA=INDA+IOFFV(IRREPA,ISPIN)
            ISCRAI(INDX(INDABSA,INDABSI,NVRTO(ISPIN)),ISPIN)=ITHRU
13         CONTINUE
12        CONTINUE
11       CONTINUE
15      CONTINUE
10     CONTINUE
C
C IJ SYMMETRY VECTOR
C
       DO 29 ISPIN=3,3-2*IUHF,-1
       ITHRU=0
       DO 20 IRREPX=1,NIRREP
        IOFFIJ(IRREPX,ISPIN)=IOFFIJ(MAX(IRREPX-1,1),ISPIN)+ITHRU
        ITHRU=0
         DO 21 IRREPJ=1,NIRREP
          IRREPI=DIRPRD(IRREPJ,IRREPX)
          IF(ISPIN.EQ.3)THEN
           NUMI=POP(IRREPI,1)
           NUMJ=POP(IRREPJ,2)
           DO 22 INDJ=1,NUMJ
            DO 23 INDI=1,NUMI
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFO(IRREPI,1)
             INDABSJ=INDJ+IOFFO(IRREPJ,2)
             ISCRIJ(INDX(INDABSI,INDABSJ,NOCCO(2)),ISPIN)=ITHRU
23          CONTINUE
22         CONTINUE
          ELSEIF(IRREPI.LT.IRREPJ.AND.ISPIN.LE.2)THEN
           NUMI=POP(IRREPI,ISPIN)
           NUMJ=POP(IRREPJ,ISPIN)
           DO 24 INDJ=1,NUMJ
            DO 25 INDI=1,NUMI
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFO(IRREPI,ISPIN)
             INDABSJ=INDJ+IOFFO(IRREPJ,ISPIN)
             ISCRIJ(INDX(INDABSI,INDABSJ,NOCCO(ISPIN)),ISPIN)=ITHRU
25          CONTINUE
24         CONTINUE
          ELSEIF(IRREPI.EQ.IRREPJ.AND.ISPIN.LE.2)THEN
           NUMI=POP(IRREPI,ISPIN)
           NUMJ=POP(IRREPJ,ISPIN)
           DO 26 INDJ=2,NUMJ
            DO 27 INDI=1,INDJ-1
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFO(IRREPI,ISPIN)
             INDABSJ=INDJ+IOFFO(IRREPJ,ISPIN)
             ISCRIJ(INDX(INDABSI,INDABSJ,NOCCO(ISPIN)),ISPIN)=ITHRU
27          CONTINUE
26         CONTINUE
          ENDIF
21       CONTINUE
         NUMBIJ(IRREPX,ISPIN)=ITHRU
20      CONTINUE
29     CONTINUE
C
C AB SYMMETRY VECTOR
C
       ITHRU=0
       DO 39 ISPIN=3,3-2*IUHF,-1
        DO 30 IRREPX=1,NIRREP
         IOFFAB(IRREPX,ISPIN)=IOFFAB(MAX(IRREPX-1,1),ISPIN)+
     &                   ITHRU
         ITHRU=0
         DO 31 IRREPJ=1,NIRREP
          IRREPI=DIRPRD(IRREPJ,IRREPX)
          IF(ISPIN.EQ.3)THEN
           NUMI=VRT(IRREPI,1)
           NUMJ=VRT(IRREPJ,2)
           DO 32 INDJ=1,NUMJ
            DO 33 INDI=1,NUMI
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFV(IRREPI,1)
             INDABSJ=INDJ+IOFFV(IRREPJ,2)
             ISCRAB(INDX(INDABSI,INDABSJ,NVRTO(2)),ISPIN)=ITHRU
33          CONTINUE
32         CONTINUE
          ELSEIF(IRREPI.LT.IRREPJ.AND.ISPIN.LE.2)THEN
           NUMI=VRT(IRREPI,ISPIN)
           NUMJ=VRT(IRREPJ,ISPIN)
           DO 34 INDJ=1,NUMJ
            DO 35 INDI=1,NUMI
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFV(IRREPI,ISPIN)
             INDABSJ=INDJ+IOFFV(IRREPJ,ISPIN)
             ISCRAB(INDX(INDABSI,INDABSJ,NVRTO(ISPIN)),ISPIN)=ITHRU
35          CONTINUE
34         CONTINUE
          ELSEIF(IRREPI.EQ.IRREPJ.AND.ISPIN.LE.2)THEN
           NUMI=VRT(IRREPI,ISPIN)
           NUMJ=VRT(IRREPJ,ISPIN)
           DO 36 INDJ=2,NUMJ
            DO 37 INDI=1,INDJ-1
             ITHRU=ITHRU+1
             INDABSI=INDI+IOFFV(IRREPI,ISPIN)
             INDABSJ=INDJ+IOFFV(IRREPJ,ISPIN)
             ISCRAB(INDX(INDABSI,INDABSJ,NVRTO(ISPIN)),ISPIN)=ITHRU
37          CONTINUE
36         CONTINUE
          ENDIF
31       CONTINUE
         NUMBAB(IRREPX,ISPIN)=ITHRU
30      CONTINUE
39     CONTINUE
C
       DO 40 ISPIN=3,3-2*IUHF,-1 
        DO 41 IRREPX=1,NIRREP
         DO 42 IRREPIJ=1,NIRREP-1
          IRREPAB=DIRPRD(IRREPIJ,IRREPX)
          IOFFABIJ(IRREPIJ+1,IRREPX,ISPIN)
     &     =IOFFABIJ(IRREPIJ,IRREPX,ISPIN)
     &                              +NUMBAB(IRREPAB,ISPIN)
     &                              *NUMBIJ(IRREPIJ,ISPIN)
42       CONTINUE
41      CONTINUE
40     CONTINUE
C
C NOW PARSE THE INPUT
C
      READ(30,*)NROOTS
      WRITE(6,5000)NROOTS
      WRITE(6,2000)
      WRITE(6,8000)
      WRITE(6,2000)
      DO 1001 III=1,NROOTS
       READ(30,*,END=992)NCONT
       IF(NCONT.EQ.0)GOTO 993
       DO 1002 JJJ=1,NCONT
500     READ(30,*,END=991)ISPIN,INDEXI,INDEXJ,INDEXA,INDEXB,FACT
        IF(ISPIN.EQ.3)THEN
         ISP1=1
         ISP2=2
        ELSE
         ISP1=ISPIN
         ISP2=ISPIN
        ENDIF
        IF(INDEXB.EQ.0)THEN
C
C SINGLE EXCITATION 
C
         IF(ISPIN.EQ.0.AND.INDEXI.EQ.0.AND.INDEXA.EQ.0)GOTO 999
         IF(ISPIN.EQ.1)THEN
          SPTYPE=' AA '
          ELSEIF(ISPIN.EQ.2)THEN
          SPTYPE=' BB '
         ENDIF
         IPACKI=IMAPMO(INDEXI,ISPIN)
         IPACKJ=0
         IPACKA=IMAPMO(INDEXA,ISPIN)-NOCCO(ISPIN)
         IPACKB=0-NOCCO(ISPIN)
         ISYMI =IRREPORB(IPACKI,ISPIN)
         ISYMA =IRREPORB(IPACKA+NOCCO(ISPIN),ISPIN)
         ISYMM=DIRPRD(ISYMA,ISYMI)
         IF(JJJ.NE.1.AND.ISYMM.NE.ISYMAI)THEN
          WRITE(6,8010)
          CALL ERREX
         ENDIF
         ISYMAI=DIRPRD(ISYMA,ISYMI)
         IRREPX=ISYMAI
         INDXAI=INDX(IPACKA,IPACKI,NVRTO(ISPIN))
         IF(JJJ.EQ.1)NROOT(ISYMAI)=NROOT(ISYMAI)+1
         IMAP2(JJJ,NROOT(ISYMAI),ISYMAI)=ISCRAI(INDXAI,ISPIN)
         ZMAP2(JJJ,NROOT(ISYMAI),ISYMAI)=FACT
C
        ELSE
C
C DOUBLE EXCITATION 
C
         IF(ISPIN.EQ.1)THEN
          SPTYPE='AAAA'
         ELSEIF(ISPIN.EQ.2)THEN
          SPTYPE='BBBB'
         ELSE
          SPTYPE='ABAB'
         ENDIF
         IPACKI=IMAPMO(INDEXI,ISP1)
         IPACKJ=IMAPMO(INDEXJ,ISP2)
         IPACKA=IMAPMO(INDEXA,ISP1)-NOCCO(ISP1)
         IPACKB=IMAPMO(INDEXB,ISP2)-NOCCO(ISP2)
         ISYMI =IRREPORB(IPACKI,ISP1)
         ISYMJ =IRREPORB(IPACKJ,ISP2)
         ISYMA =IRREPORB(IPACKA+NOCCO(ISP1),ISP1)
         ISYMB =IRREPORB(IPACKB+NOCCO(ISP2),ISP2)
         ISYMAB=DIRPRD(ISYMA,ISYMB)
         ISYMIJ=DIRPRD(ISYMI,ISYMJ)
         INDXAB=INDX(IPACKA,IPACKB,NVRTO(ISP1))
         INDXIJ=INDX(IPACKI,IPACKJ,NOCCO(ISP1))
         IPOSAB=ISCRAB(INDXAB,ISPIN)
         IPOSIJ=ISCRIJ(INDXIJ,ISPIN)
         ISYMM =DIRPRD(ISYMAB,ISYMIJ)
         IF(JJJ.NE.1.AND.ISYMM.NE.IRREPX)THEN
          WRITE(6,8010)
          CALL ERREX
         ENDIF
         IRREPX=DIRPRD(ISYMAB,ISYMIJ)
         IOFF0 =IOFFABIJ(ISYMIJ,IRREPX,ISPIN)
         IOFF1 =INDX(IPOSAB,IPOSIJ,NUMBAB(ISYMAB,ISPIN))
         LENT1 =IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10)
         IF(ISPIN.LT.3)THEN
          LENT1=LENT1+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
          IF(ISPIN.EQ.2)THEN
           LENT1=LENT1+IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
          ENDIF
         ENDIF
         IF(JJJ.EQ.1)NROOT(IRREPX)=NROOT(IRREPX)+1
c       write(6,*)' positions '
c       write(6,*)' i,j,a,b           ',ipacki,ipackj,ipacka,ipackb
c       write(6,*)' irrepij           ',isymij
c       write(6,*)' irrepab           ',isymab
c       write(6,*)' irrepx            ',irrepx
c       write(6,*)' iposab and iposij ',iposab,iposij
c       write(6,*)' absolute offset   ',iabspos
c       write(6,*)' ioff0             ',ioff0
c       write(6,*)' ioff1             ',ioff1
c       write(6,*)' absolute offset   ',ioff0+ioff1+lent1
         IMAP2(JJJ,NROOT(IRREPX),IRREPX)=IOFF0+IOFF1+LENT1
         ZMAP2(JJJ,NROOT(IRREPX),IRREPX)=FACT
C
        ENDIF
        WRITE(6,8001)III,JJJ,IPACKI,IPACKJ,IPACKA+NOCCO(ISP1),
     &               IPACKB+NOCCO(ISP2),SPTYPE,FACT
999     CONTINUE
1002   CONTINUE
       WRITE(6,2004)IRREPX 
1001  CONTINUE
      WRITE(6,2000)
C
      RETURN
2000  FORMAT(71('-'))
2001  FORMAT(T8,'Hole orbital',T32,'Particle orbital',T53,'Transition ',
     &       /,
     &       T5,'Offset',T14,'Eigenvalue',T30,'Offset',T39,'Eigenvalue',
     &       T54,'Symmetry',T66,'Spin') 
2002  FORMAT(T7,I3,T13,F12.6,T32,I3,T38,F12.6,T57,I1,T68,A1)
2003  FORMAT(T7,4I3,F20.10)
2004  FORMAT(T3,'Guess vectors transform as symmetry ',I1,'.')
8000  FORMAT('Root',T7,'Element',T21,'I',T25,'J',T30,'A',T34,
     &       'B',T40,'Spin Case',T55,'Weight')
8001  FORMAT(T2,I2,T8,I3,T20,I3,T24,I3,T29,I3,T33,I3,T43,A4,T51,F10.4)
8010  FORMAT(T3,'@PARSEINP-F, Input vector has mixed symmetry.')
5000  FORMAT(T3,'@PARSEINP-I, Starting vectors supplied for ',I5,
     &       ' roots.')
8011  FORMAT(T3,'@PARSEINP-F, ',I5,' contributions not specified ',
     &       'for root ',I5,'.')
8012  FORMAT(T3,'@PARSEINP-F, Input for root ',I5,' not found.')
8013  FORMAT(T3,'@PARSEINP-F, At least one contribution required for ',
     &          'root ',I5,'.')
C
991   WRITE(6,8011)NCONT,III
      CALL ERREX
      RETURN
C
992   WRITE(6,8012)III
      CALL ERREX
      RETURN
C
993   WRITE(6,8013)III
      CALL ERREX
      RETURN
C
      END
