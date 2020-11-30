      SUBROUTINE DRHFRNG_MODF(ICORE,MAXCOR,IUHF,IRREPX,DISCO,INCREM,
     &                        WSPIN,LISTT0,LISTZ0,LSTTMP0,ISIDE)
C
C THIS SUBROUTINE CALCULATES THE RING CONTRIBUTION FOR RHF
C REFERENCES, USING THE SPIN-ADAPTED SINGLET FORMULA:
C
C                 _        _ 
C  Z(ab,ij) = 1/2 T(ae,im) W(mb,ej) - 1/2 T(Im,Ea) W(Mb,jE) 
C
C           - T(Im,Eb) W(Ma,Je)
C
C       _
C WHERE X IS SHORTHAND FOR 2*X(ij,kl)-X(ji,kl)
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*2 MATTYP(2)
      DOUBLE PRECISION ONE,ONEM,ZILCH,HALF,TWO
      LOGICAL DISCO,INCREM,WSPIN
      LOGICAL  MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
      DIMENSION ICORE(MAXCOR),IDID(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
C
      DATA ONE  /1.0D0/
      DATA ONEM /-1.0D0/
      DATA ZILCH/0.0D0/
      DATA HALF /0.5D0/
      DATA TWO  /2.0D0/
      DATA MATTYP /'N','T'/
C
      LISTZ =LISTZ0-1 +3
      LSTTMP=LSTTMP0-1+3
      DO 100 IRREPZR=1,NIRREP
       IRREPZL=DIRPRD(IRREPZR,IRREPX)
       DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LSTTMP))
       NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LSTTMP))
       I000=1
       I010=I000+IINTFP*DISSYZ*NUMDSZ
       IRREPW=IRREPZR
       IRREPTR=IRREPZR
       IRREPTL=IRREPZL
C
C FORM SPIN-ADAPTED T INTERMEDIATE.  SPIN-ADAPTED INTEGRALS ARE ALREADY
C ON DISK
C                _          _
C     Z(AI,BJ) = T(AI,EM) * W(EM,BJ)
C
C 56 -> 118
C 
       LISTW=118
       LISTT1=LISTT0-1+4
       LISTT2=LISTT0-1+6
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT1))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT1))
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW))
       I020=I010+IINTFP*DISSYT*NUMDST
       I030=I020+IINTFP*DISSYW*NUMDSW

       IF (CC2) THEN

       CALL ZERO(ICORE(I000),DISSYZ*NUMDSZ)

       ELSE
       CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPTR,LISTT1)
       CALL SSCAL (NUMDST*DISSYT,TWO,ICORE(I010),1)
       CALL GETLST(ICORE(I020),1,NUMDST,1,IRREPTR,LISTT2)
       CALL SAXPY (NUMDST*DISSYT,ONEM,ICORE(I020),1,ICORE(I010),1)
       IF(.NOT.WSPIN)THEN
C
C 54,58 -> 123, 125
C
        LISTWA=123
        LISTWB=125
        CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPW,LISTWA)
        CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTWB)
        CALL SAXPY (NUMDSW*DISSYW,-TWO,ICORE(I000),1,ICORE(I020),1)
       ELSE
        CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
       ENDIF
       CALL XGEMM('N',MATTYP(ISIDE),DISSYT,NUMDSW,DISSYW,-HALF,
     &            ICORE(I010),DISSYT,ICORE(I020),DISSYW,ZILCH,
     &            ICORE(I000),DISSYZ)

       ENDIF 

       IF(INCREM)THEN
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPZR,LSTTMP)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEm,ICORE(I010),1,ICORE(I000),1)
       ENDIF
C
C DISCO CONTRIBUTION: Z(AI,BJ) = T(AI) * F(BJ)
C
       IF(DISCO.AND.IRREPZR.EQ.1)THEN
        I020=I010+IINTFP*NT(1)
        I030=I020+IINTFP*IRPDPD(IRREPX,9)
        CALL GETLST(ICORE(I010),1,1,1,1,93)
        CALL GETLST(ICORE(I020),1,1,1,1,490)
        NROW=IRPDPD(IRREPX,9)
        NCOL=NT(1)
        CALL XGEMM('N','N',NROW,NCOL,1,ONEM,ICORE(I020),NROW,
     &             ICORE(I010),1,ONE,ICORE(I000),NROW)
       ENDIF
       CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LSTTMP)
C
C    Z(ai,bj) = T(Im,Ea) * W(Mb,Je)
C
C 58 -> 125
C
       LISTW=125
       LISTT=LISTT0-1+6
       DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
       NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
       DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
       NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW))
       I020=I010+IINTFP*DISSYT*NUMDST
       I030=I020+IINTFP*DISSYW*NUMDSW

       IF (CC2) THEN

       CALL ZERO(ICORE(I000),NUMDSZ*DISSYZ)

       ELSE

       CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPTR,LISTT)
       CALL GETLST(ICORE(I020),1,NUMDSW,1,IRREPW,LISTW)
       CALL XGEMM('N',MATTYP(ISIDE),DISSYT,NUMDSW,DISSYW,ONE,
     &            ICORE(I010),DISSYT,ICORE(I020),DISSYW,ZILCH,
     &            ICORE(I000),DISSYZ)
      ENDIF 

       CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LSTTMP+1)

100   CONTINUE

C
C NOW READ IN aI,bJ INCREMENTS AND SSTRNG THEM
C
      ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LSTTMP+1),ISYTYP(2,LSTTMP+1))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      CALL GETALL(ICORE(I010),ISIZE,IRREPX,LSTTMP+1)
      CALL SSTGEN(ICORE(I010),ICORE(I000),ISIZE,VRT(1,1),POP(1,1),
     &            VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1432')
      CALL SAXPY (ISIZE,HALF,ICORE(I000),1,ICORE(I010),1)
c      IF(INCREM)THEN
c       CALL GETALL(ICORE(I000),ISIZE,IRREPX,441)
c       CALL SAXPY (ISIZE,ONE,ICORE(I000),1,ICORE(I010),1)
c      ENDIF
      CALL SSTGEN(ICORE(I010),ICORE(I000),ISIZE,VRT(1,1),POP(1,1),
     &            VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1432')
      CALL GETALL(ICORE(I010),ISIZE,IRREPX,LSTTMP)
      CALL SAXPY (ISIZE,ONE,ICORE(I010),1,ICORE(I000),1)
      CALL PUTALL(ICORE(I000),ISIZE,IRREPX,LSTTMP)
C
C SYMMETRIZE INCREMENTS
C
      CALL IZERO(IDID,NIRREP)
      DO 200 IRREPR=1,NIRREP
       IF(IDID(IRREPR).EQ.0)THEN
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LSTTMP))
        DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LSTTMP))
        I000=1
        I010=I000+IINTFP*NUMDIS*DISSIZ
        IF(IRREPX.EQ.1)THEN
         CALL GETLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
         CALL MPMT  (ICORE(I000),DISSIZ)
         CALL PUTLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
        ELSE
         I020=I010+IINTFP*NUMDIS*DISSIZ
         CALL GETLST(ICORE(I000),1,DISSIZ,1,IRREPL,LSTTMP)
         CALL TRANSP(ICORE(I000),ICORE(I010),DISSIZ,NUMDIS)
         CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREPR,LSTTMP)
         CALL SAXPY (NUMDIS*DISSIZ,ONE,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMDIS,1,IRREPR,LSTTMP)
         CALL TRANSP(ICORE(I000),ICORE(I010),NUMDIS,DISSIZ)
         CALL PUTLST(ICORE(I010),1,DISSIZ,1,IRREPL,LSTTMP)
        ENDIF
        IDID(IRREPR)=1
        IDID(IRREPL)=1
       ENDIF
200   CONTINUE
C
      ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTZ),ISYTYP(2,LISTZ))
      I000=1
      I010=I000+IINTFP*ISIZE
      I020=I010+IINTFP*ISIZE
      CALL GETALL(ICORE(I000),ISIZE,IRREPX,LSTTMP)
      CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),POP(1,1),
     &            VRT(1,1),POP(1,1),ICORE(I020),IRREPX,'1324')
      CALL GETALL(ICORE(I000),ISIZE,IRREPX,LISTZ)
      CALL SAXPY (ISIZE,ONEM,ICORE(I010),1,ICORE(I000),1)
      CALL PUTALL(ICORE(I000),ISIZE,IRREPX,LISTZ)
C
      RETURN
      END
