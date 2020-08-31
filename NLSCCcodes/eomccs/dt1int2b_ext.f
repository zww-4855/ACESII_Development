










      SUBROUTINE DT1INT2B_EXT(ICORE,MAXCOR,IUHF,IRREPT,IRREPW,
     &                    LISTT1,IOFFT1,LISTW0,LISTZ0,ANTI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD,DISSYW,DISSYZ,E
      LOGICAL ANTI
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0T(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      DATA ZILCH,ONE,ONEM/0.D0,1.D0,-1.D0/
      FACT=ONE
      IF(ANTI) FACT=ONEM
      I0T(1)=1
      I0T(2)=I0T(1)+IUHF*IINTFP*IRPDPD(IRREPT,9)
      I000=I0T(2)+IINTFP*IRPDPD(IRREPT,10)
      CALL GETLST(ICORE(I0T(1)),1,1,1,1+IOFFT1,LISTT1)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0T(2)),1,1,1,2+IOFFT1,LISTT1)
      IRREPX=DIRPRD(IRREPW,IRREPT)
      DO 100 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPX,IRREPZR)
       IRREPWL=IRREPZL
       IRREPWR=DIRPRD(IRREPW,IRREPWL)
       LISTW=LISTW0+3
       LISTZ=LISTZ0+2
       DISSYW=IRPDPD(IRREPWL,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPWR,ISYTYP(2,LISTW))
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
       MAXBUF=MAX(NUMDSZ,NUMDSW,DISSYW,DISSYZ)
       I010=I000+IINTFP*DISSYZ*NUMDSZ
       I020=I010+IINTFP*MAX(DISSYW*NUMDSW,NUMDSZ*DISSYZ)
       ITMP1=I020
       ITMP2=ITMP1+IINTFP*MAXBUF
       ITMP3=ITMP2+IINTFP*MAXBUF
       IEND=ITMP3+IINTFP*MAXBUF
       IF(IEND.LT.MAXCOR) THEN
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPWR,LISTW)
        CALL SYMTR1(IRREPWR,VRT(1,1),POP(1,2),DISSYW,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        DO 120 IRREPE=1,NIRREP
         IRREPI=DIRPRD(IRREPE,IRREPT)
         IRREPJ=DIRPRD(IRREPE,IRREPWR)
         NROW=DISSYW*POP(IRREPJ,2)
         NCOL=POP(IRREPI,1)
         NSUM=VRT(IRREPE,1)
         IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,24)-1)
         IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPWR,25)-1)
         IT=I0T(1)+IINTFP*(ISYMOFF(IRREPI,IRREPT,9)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,FACT,ICORE(IW),NROW,
     &              ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
120     CONTINUE
        CALL SYMTR1(IRREPZR,POP(1,2),POP(1,1),DISSYZ,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
       ELSE
        IZ=I000
        IDISW1=1
        DO 130 IRREPJ=1,NIRREP
         IRREPE=DIRPRD(IRREPJ,IRREPWR)
         IRREPI=DIRPRD(IRREPJ,IRREPZR)
         NUMJ=POP(IRREPJ,2)
         NUMI=POP(IRREPI,1) 
         NUME=VRT(IRREPE,1)
         ITOP=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYW*NUME) 
         IF(ITOP.GT.MAXCOR) CALL INSMEM('DT1INT2B',ITOP,MAXCOR)
         DO 131 J=1,NUMJ
          CALL GETLST(ICORE(I010),IDISW1,NUME,1,IRREPWR,LISTW)
          NROW=DISSYW
          NCOL=NUMI
          NSUM=NUME
          IT=I0T(1)+IINTFP*(ISYMOFF(IRREPI,IRREPT,9)-1)
          CALL XGEMM('N','N',NROW,NCOL,NSUM,FACT,ICORE(I010),NROW,
     &               ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW) 
          IZ=IZ+IINTFP*DISSYZ*NUMI
          IDISW1=IDISW1+NUME
131      CONTINUE
130     CONTINUE
       ENDIF
       IF(IUHF.EQ.0) THEN
        ITMP1=I010
        ITMP2=ITMP1+IINTFP*MAXBUF
        ITMP3=ITMP2+IINTFP*MAXBUF
        CALL SYMRHF3(IRREPZL,IRREPZR,VRT(1,1),POP(1,1),DISSYZ,
     &               ICORE(I000),ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(ITMP3))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY(NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
       ELSE
        LISTW=LISTW0+2
        DISSYW=IRPDPD(IRREPWL,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPWR,ISYTYP(2,LISTW))
        IEND=I010+IINTFP*MAX(NUMDSZ*DISSYZ,NUMDSW*DISSYW)
        IF(IEND.LT.MAXCOR) THEN
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPWR,LISTW)
         DO 150 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPT)
          IRREPJ=DIRPRD(IRREPE,IRREPWR)
          NROW=DISSYW*POP(IRREPJ,1)
          NCOL=POP(IRREPI,2)
          NSUM=VRT(IRREPE,2)
          IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,14)-1)
          IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPWR,18)-1)
          IT=I0T(2)+IINTFP*(ISYMOFF(IRREPI,IRREPT,10)-1)
          CALL XGEMM('N','N',NROW,NCOL,NSUM,FACT,ICORE(IW),NROW,
     &               ICORE(IT),NSUM,ONE,ICORE(IZ),NROW)
150      CONTINUE
        ELSE
         IDISW1=1
         DO 151 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPT)
          IRREPJ=DIRPRD(IRREPE,IRREPZR)
          IZ=I000+IINTFP*(ISYMOFF(IRREPI,IRREPZR,14)-1) 
          NUMJ=POP(IRREPJ,1)
          NUMI=POP(IRREPI,2)
          NUME=VRT(IRREPE,2)
          ITOP=I010+IINTFP*MAX(NUMDSZ*DISSYZ,DISSYW*NUMJ)
          IF(ITOP.GE.MAXCOR) CALL INSMEM('DT1INT2B',ITOP,MAXCOR)
          DO 152 E=1,NUME
           CALL GETLST(ICORE(I010),IDISW1,NUMJ,1,IRREPWR,LISTW)
           NROW=DISSYZ*NUMJ
           NCOL=NUMI
           NSUM=1
           IT=I0T(2)+IINTFP*(ISYMOFF(IRREPI,IRREPT,10)+E-2)
           CALL XGEMM('N','N',NROW,NCOL,NSUM,FACT,ICORE(I010),NROW,
     &                ICORE(IT),NUME,ONE,ICORE(IZ),NROW)
           IDISW1=IDISW1+NUMJ
152       CONTINUE
151      CONTINUE
        ENDIF
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY(NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
       ENDIF
100   CONTINUE
      IF(IUHF.EQ.0) RETURN
      DO 200 ISPIN=1,1+IUHF
       DO 210 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZR,IRREPX) 
        IRREPWL=IRREPZL
        IRREPWR=DIRPRD(IRREPW,IRREPWL)
        LISTW=LISTW0-1+ISPIN
        LISTZ=LISTZ0-1+ISPIN
        DISSYW=IRPDPD(IRREPWL,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPWR,ISYTYP(2,LISTW))
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        NUMDSZF=IRPDPD(IRREPZR,20+ISPIN)
        MAXBUF=MAX(DISSYW,NUMDSW,DISSYZ,NUMDSZ,NUMDSZF)
        I010=I000+IINTFP*DISSYZ*NUMDSZF
        I020=I010+IINTFP*MAX(NUMDSZ*DISSYZ,NUMDSW*DISSYW)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXBUF 
        ITMP3=ITMP2+IINTFP*MAXBUF
        IEND=ITMP3+IINTFP*MAXBUF
        IF(IEND.LT.MAXCOR) THEN
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPWR,LISTW) 
         CALL SYMTR1(IRREPWR,VRT(1,ISPIN),POP(1,ISPIN),DISSYW, 
     &               ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     &               ICORE(ITMP3))
         DO 230 IRREPE=1,NIRREP
          IRREPI=DIRPRD(IRREPE,IRREPT)
          IRREPJ=DIRPRD(IRREPE,IRREPWR)
          NROW=DISSYW*POP(IRREPJ,ISPIN)
          NCOL=POP(IRREPI,ISPIN)
          NSUM=VRT(IRREPE,ISPIN)
          IZ=I000+IINTFP*DISSYZ*(ISYMOFF(IRREPI,IRREPZR,20+ISPIN)-1)
          IW=I010+IINTFP*DISSYW*(ISYMOFF(IRREPE,IRREPWR,15+ISPIN)-1)
          IT=I0T(ISPIN)+IINTFP*(ISYMOFF(IRREPI,IRREPT,8+ISPIN)-1)
          CALL XGEMM('N','N',NROW,NCOL,NSUM,-FACT,ICORE(IW),NROW,
     &               ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
230      CONTINUE
        ELSE
         IZ=I000
         IDSIW1=1
         DO 231 IRREPJ=1,NIRREP
          IRREPE=DIRPRD(IRREPJ,IRREPWR) 
          IRREPI=DIRPRD(IRREPJ,IRREPZR)
          NUMJ=POP(IRREPJ,ISPIN) 
          NUMI=POP(IRREPI,ISPIN) 
          NUME=VRT(IRREPE,ISPIN)
          ITOP=I010+IINTFP*MAX(NUMDSZ*DISSYZ,DISSYW*NUME)
          IF(ITOP.GE.MAXCOR) CALL INSMEM('DT1INT2B',ITOP,MAXCOR)
          DO 232 J=1,NUMJ
           CALL GETLST(ICORE(I010),IDISW1,NUME,1,IRREPWR,LISTW)
           NROW=DISSYW
           NCOL=NUMI
           NSUM=NUME
           IT=I0T(ISPIN)+IINTFP*(ISYMOFF(IRREPI,IRREPT,8+ISPIN)-1)
           CALL XGEMM('N','N',NROW,NCOL,NSUM,-FACT,ICORE(I020),NROW,
     &                ICORE(IT),NSUM,ZILCH,ICORE(IZ),NROW)
           IZ=IZ+IINTFP*DISSYZ*NUMI
           IDISW1=IDISW1+NUME
232       CONTINUE
231      CONTINUE
        ENDIF
        CALL ASSYM2(IRREPZR,POP(1,ISPIN),DISSYZ,ICORE(I000))
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL SAXPY(NUMDSZ*DISSYZ,ONE,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
210    CONTINUE
200   CONTINUE
      RETURN
      END 
