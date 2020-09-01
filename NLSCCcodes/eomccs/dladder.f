










      SUBROUTINE DLADDER(ICORE,MAXCOR,IUHF,IRREPX,ITYPE,
     &                   LISTT0,LISTW0,LISTZ0,ISIDE)
C
C THIS SUBROUTINE CALCULATES THE HOLE-HOLE LADDER CONTRIBUTION
C
C     Z(ab,ij) = SUM T(ab,mn) * <mn||ij>  [ITYPE=1] 
C                m,n
C
C     Z(ab,ij) = SUM W(ab,ef) * T(ef,ij)  [ITYPE=6]
C                e,f
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH,FACT
      CHARACTER*1 MATTYP(2)
C
      DIMENSION ICORE(MAXCOR),I0L(2)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE  /1.0D0/
      DATA ZILCH/0.0D0/
      DATA MATTYP/'N','T'/


      DO 50 ISPIN=3,3-2*IUHF,-1
C
C ALPHA-BETA SPIN CASE
C
       DO 100 IRREPZR=1,NIRREP
C
C LOOP OVER KET IRREPS OF *TARGET*.  THIS IS NOT THE SAME AS THE
C  IRREPS OF THE INTEGRALS AND AMPLITUDES UNLESS IRREPX=1.
C
        IRREPZL=DIRPRD(IRREPZR,IRREPX)
        LISTW=LISTW0+ISPIN-1
        LISTT=LISTT0+ISPIN-1
        LISTZ=LISTZ0+ISPIN-1
        IF(ITYPE.EQ.1)THEN
         IRREPW=IRREPZR
        ELSE
         IRREPW=IRREPZL
        ENDIF
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        DISSYT=DISSYZ
        NUMDST=NUMDSZ
        I000=1
        I010=I000+IINTFP*DISSYZ*NUMDSZ
        I020=I010+IINTFP*DISSYT*NUMDST
C
C USE GENERAL ALGORITHM ALLOWING BOTH IN-CORE AND OUT-OF-CORE
C SOLUTIONS
C
        FACT=ZILCH
        CORLFT=MAXCOR-I020+1
        IF(DISSYW.NE.0)THEN
         NINCOR=CORLFT/(DISSYW*IINTFP)
         IF(NINCOR.EQ.0)THEN
          WRITE(6,1000)
1000      FORMAT(T3,'@DLADDER-F, Not enough memory for ladders.')
          CALL INSMEM('DLADDER',DISSYW*IINTFP,CORLFT)
         ENDIF
        ELSE
         NINCOR=DISSYW
        ENDIF
        IOFFT=I010
        IOFFZ=I000
        NLEFT=NUMDSW
        NFIRST=1
        CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPZR,LISTT)
1       NREAD =MIN(NLEFT,NINCOR)
        CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREPW,LISTW)
        IF(ITYPE.EQ.1)THEN
         CALL XGEMM('N',MATTYP(ISIDE),DISSYZ,NREAD,NUMDSZ,ONE,
     &              ICORE(IOFFT),DISSYT,ICORE(I020),DISSYW,FACT,
     &              ICORE(IOFFZ),DISSYZ)
         IOFFZ=IOFFZ+DISSYZ*NREAD*IINTFP
        ELSE
         IF(ISIDE.EQ.1)THEN
          CALL XGEMM ('N','N',DISSYZ,NUMDSZ,NREAD,ONE,ICORE(I020),
     &                DISSYW,ICORE(IOFFT),DISSYT,FACT,ICORE(IOFFZ),
     &                DISSYZ)
          IOFFT=IOFFT+IINTFP*NREAD
          FACT =ONE
         ELSE
          CALL XGEMM ('T','N',NREAD,NUMDSZ,DISSYZ,ONE,ICORE(I020),
     &                DISSYW,ICORE(IOFFT),DISSYT,FACT,ICORE(IOFFZ),
     &                DISSYZ)
          IOFFZ=IOFFZ+NREAD*IINTFP
         ENDIF
        ENDIF
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT -NREAD
        IF(NLEFT.NE.0)GOTO 1
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONE,ICORE(I000),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
C
C FORM CONTRIBUTION TO L1 INCREMENT FOR LEFT-HAND SIDE SOLUTION.
C THIS IS NECESSITATED BY THE FACT THAT THE CONTRIBUTION TO THE
C ABCI HBAR INTERMEDIATE ARISING FROM THE ABCD INTEGRALS IS NOT
C CALCULATED.  THAT IS, THIS ROUTINE CALCULATES A MISSING CONTRIBUTION
C OF HBAR(AIBC)-> L1 (A PART OF L2INL1).  THE CURRENT INCREMENT IS
C HELD AT ICORE(I000)
C
        IF(ISIDE.EQ.2.AND.ITYPE.EQ.6)THEN
         IRREPIJ=IRREPZR
         IRREPAB=IRREPZL
         IF(IUHF.NE.0)THEN
          NUMABX=IRPDPD(IRREPAB,18+MIN(ISPIN,2))
          NUMIJX=IRPDPD(IRREPIJ,20+MIN(ISPIN,2))
         ELSE
          NUMABX=DISSYZ
          NUMIJX=NUMDSZ
         ENDIF
         MAXBUF=MAX(NUMABX,NUMIJX,NT(1),NT(2),IRPDPD(IRREPX,9),
     &              IRPDPD(IRREPX,10),DISSYZ,NUMDSZ)
         I0L(1)=I000+IINTFP*MAX(DISSYZ*NUMDSZ,NUMABX*NUMIJX)
         CALL GETLST(ICORE(I0L(1)),1,1,1,3,490)
         IF(IUHF.NE.0)THEN
          I0L(2)=I0L(1)+IINTFP*IRPDPD(IRREPX,9)
          CALL GETLST(ICORE(I0L(2)),1,1,1,4,490)
         ELSE
          I0L(2)=I0L(1)
         ENDIF
         I010=I0L(2)+IINTFP*IRPDPD(IRREPX,10)
         IF(ISPIN.EQ.3.AND.IUHF.EQ.0)THEN
          I020=I010+IINTFP*MAXBUF
          I030=I020+IINTFP*MAXBUF
          CALL SPINAD1(IRREPIJ,POP(1,1),DISSYZ,ICORE(I000),
     &                 ICORE(I010),ICORE(I020))
          CALL GETLST(ICORE(I010),1,1,1,1,90)
          CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &                 ICORE(I0L(1)),ICORE(I010),ICORE(I000),
     &                 ICORE(I020),DISSYZ,VRT(1,1),POP(1,1),
     &                 VRT(1,1),VRT(1,1),POP(1,1),POP(1,1),
     &                 'STST')
          CALL PUTLST(ICORE(I0L(1)),1,1,1,3,490)
         ELSEIF(ISPIN.LE.2.AND.ISIDE.EQ.2)THEN
          I020=I010+IINTFP*NUMABX*NUMIJX
          I030=I020+IINTFP*MAXBUF
          I040=I030+IINTFP*MAXBUF
          CALL SYMEXP (IRREPIJ,POP(1,ISPIN),DISSYZ,ICORE(I000))
          CALL SYMEXP2(IRREPAB,VRT(1,ISPIN),NUMABX,DISSYZ,NUMIJX,
     &                 ICORE(I010),ICORE(I000))
          CALL GETLST(ICORE(I020),1,1,1,ISPIN,90)
          CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &                 ICORE(I0L(ISPIN)),ICORE(I020),ICORE(I010),
     &                 ICORE(I030),NUMABX,VRT(1,ISPIN),POP(1,ISPIN),
     &                 VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &                 POP(1,ISPIN),'TSTS')
          CALL PUTLST(ICORE(I0L(ISPIN)),1,1,1,2+ISPIN,490)
         ELSEIF(ISPIN.EQ.3.AND.IUHF.NE.0)THEN
          I020=I010+IINTFP*MAXBUF
          I030=I020+IINTFP*MAXBUF
          CALL GETLST(ICORE(I010),1,1,1,2,90)
          CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &                 ICORE(I0L(1)),ICORE(I010),ICORE(I000),
     &                 ICORE(I020),DISSYZ,VRT(1,1),POP(1,1),
     &                 VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &                 'TSTS')
          CALL GETLST(ICORE(I010),1,1,1,1,90)
          CALL DDDOT24(IRREPX,1,IRREPAB,IRREPIJ,
     &                 ICORE(I0L(2)),ICORE(I010),ICORE(I000),
     &                 ICORE(I020),DISSYZ,VRT(1,2),POP(1,2),
     &                 VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &                 'STST')
          CALL PUTLST(ICORE(I0L(1)),1,1,1,3,490)
          CALL PUTLST(ICORE(I0L(2)),1,1,1,4,490)
         ENDIF
        ENDIF
100    CONTINUE
50    CONTINUE
C
      RETURN
      END