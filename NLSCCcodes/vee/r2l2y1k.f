      SUBROUTINE R2L2Y1K(Y1,ICORE,MAXCOR,IUHF)
C
C Y2(ai) = -1/4 R(ef,mn)*L(ef,oi)*W(mn,oa)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH
      DIMENSION ICORE(MAXCOR),Y1(*)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C          
C Y2(AI)=  [2R(Ef,Mn)-R(Fe,Mn)]*L(Ef,Oi)*W(Mn,Oa)
C
C
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        LISTL=446
        LISTR=463
        LISTW=10
        DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYQ=NUMDSR
        NUMDSQ=NUMDSL
        MAXT=MAX(NUMDSR,DISSYR,DISSYW,NUMDSW,NUMDSQ,DISSYQ,
     &           NUMDSL,DISSYL)
        I000=1
        I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYL,
     &                       NUMDSW*DISSYW,NUMDSQ*DISSYQ)
        I020=I010+IINTFP*MAX(NUMDSW*DISSYW,NUMDSR*DISSYR,
     &                       NUMDSL*DISSYL,NUMDSQ*DISSYQ)
        ITMP1=I010
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT 
        CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
        CALL SPINAD3(IRREPL,VRT(1,1),DISSYR,NUMDSR,ICORE(I000),
     &               ICORE(ITMP1),ICORE(ITMP2))
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
C                              _       +
C FORM INTERMEDIATE Q(Mn,Oi) = R(Ef,Mn) * L(Ef,Oi)
C
        CALL XGEMM('T','N',NUMDSR,NUMDSL,DISSYR,ONE,ICORE(I000),
     &             DISSYR,ICORE(I010),DISSYL,ZILCH,ICORE(I020),NUMDSR)
C
C NOW READ IN W(Mn,Oa)
C
        CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPR,LISTW)
C
C FORM PRODUCT Y(ai) = W(Mn,Oa) * Q(Mn,Oi)
C
        IOFFW=I000
        IOFFQ=I020
        IOFFY=1 
        DO 11 IRREPI=1,NIRREP
         IRREPA=IRREPI
         IRREPO=DIRPRD(IRREPA,IRREPR)
         NUMI=POP(IRREPI,1)
         NUMA=VRT(IRREPA,1)
         NUMO=POP(IRREPO,1)
         NROW=NUMA
         NCOL=NUMI
         NSUM=DISSYW*NUMO
         CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &              ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
         IOFFW=IOFFW+IINTFP*NROW*NSUM
         IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
         IOFFY=IOFFY+IINTFP*NROW*NCOL
11      CONTINUE
C
10     CONTINUE
C
      ELSE
C
C Y2(AI) <= -1/4 R(EF,MN)*L(EF,OI)*W(MN,OA)
C
       IOFFY0=1 
       DO 20 ISPIN=1,2
        DO 21 IRREPR=1,NIRREP
         IRREPL=DIRPRD(IRREPR,IRREPX)
         LISTL=443+ISPIN
         LISTR=460+ISPIN
         LISTW=6+ISPIN
         DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
         NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
         NUMDSW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
         DISSYQ=NUMDSR
         NUMDSQ=NUMDSL
         I000=1
         I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYL,
     &                        NUMDSW*DISSYW)
         I020=I010+IINTFP*MAX(NUMDSW*DISSYW,NUMDSR*DISSYR,
     &                        NUMDSL*DISSYL)
         CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
         CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
C                                         +
C FORM INTERMEDIATE Q(M<N,O<I) = R(E<F,M<N) * L(E<F,O<I)
C
         CALL XGEMM('T','N',NUMDSR,NUMDSL,DISSYR,ONE,ICORE(I000),
     &              DISSYR,ICORE(I010),DISSYL,ZILCH,ICORE(I020),DISSYQ)
C
C EXPAND TO Q(M<N,OI)
C
         CALL SYMEXP(IRREPR,POP(1,ISPIN),NUMDSL,ICORE(I020))
C
C NOW READ IN W(M<N,OA)
C
         CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPR,LISTW)
C
C FORM PRODUCT Y(AI) <= W(M<N,OA) * Q(M<N,OI)
C
         IOFFW=I000
         IOFFQ=I020
         IOFFY=IOFFY0
         DO 22 IRREPI=1,NIRREP
          IRREPA=IRREPI
          IRREPO=DIRPRD(IRREPA,IRREPR)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMO=POP(IRREPO,ISPIN)
          NROW=NUMA
          NCOL=NUMI
          NSUM=DISSYQ*NUMO
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &               ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
          IOFFW=IOFFW+IINTFP*NROW*NSUM
          IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
          IOFFY=IOFFY+IINTFP*NROW*NCOL
22       CONTINUE
C
C Y2(AI) <=  R(Ef,Mn)*L(Ef,Io)*W(Mn,Ao) [ISPIN=1]
C Y2(AI) <=  R(Ef,Mn)*L(Ef,Oi)*W(Mn,Oa) [ISPIN=2]
C
         LISTL=446
         LISTR=463
         LISTW=8+ISPIN
         DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
         NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
         DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
         NUMDSW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYQ=NUMDSR
         NUMDSQ=NUMDSL
         MAXT=MAX(NUMDSR,DISSYR,DISSYW,NUMDSW,NUMDSQ,DISSYQ,
     &            NUMDSL,DISSYL)
         I000=1
         I010=I000+IINTFP*MAX(NUMDSR*DISSYR,NUMDSL*DISSYL,
     &                        NUMDSW*DISSYW,MAXT*3)
         I020=I010+IINTFP*MAX(NUMDSW*DISSYW,NUMDSR*DISSYR,
     &                        NUMDSL*DISSYL,3*MAXT)
         ITMP1=I010
         ITMP2=ITMP1+IINTFP*MAXT
         ITMP3=ITMP2+IINTFP*MAXT 
         CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPR,LISTR)
         CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
C                                       +
C FORM INTERMEDIATE Q(Mn,Io)  = R(Ef,Mn) * L(Ef,Io) [ISPIN=1] 
C                   Q(Mn,Oi)  = R(Ef,Mn) * L(Ef,Oi) [ISPIN=2] 
C
         CALL XGEMM('T','N',DISSYQ,NUMDSQ,DISSYR,ONE,ICORE(I000),
     &              DISSYR,ICORE(I010),DISSYL,ZILCH,ICORE(I020),DISSYQ)
C
C NOW READ IN W(Mn,Ao)  [ISPIN=1]
C NOW READ IN W(Mn,Oa)  [ISPIN=2]
C
         CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPR,LISTW)
C
C W(Mn,Ao) -> W(Mn,oA) [ISPIN=1 ONLY]
C Q(Mn,Io) -> Q(Mn,oI) [ISPIN=1 ONLY]
C
         IF(ISPIN.EQ.1)THEN
          CALL SYMTR1(IRREPR,POP(1,1),POP(1,2),DISSYQ,ICORE(I020),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
          CALL SYMTR1(IRREPR,VRT(1,1),POP(1,2),DISSYW,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
         ENDIF
C
C FORM PRODUCT Y(AI) <= W(Mn,oA)*Q(Mn,oI) [ISPIN=1]
C FORM PRODUCT Y(AI) <= W(Mn,Oa)*Q(Mn,Oi) [ISPIN=2]
C
         IOFFW=I000
         IOFFQ=I020
         IOFFY=IOFFY0
         DO 23 IRREPI=1,NIRREP
          IRREPA=IRREPI
          IRREPO=DIRPRD(IRREPA,IRREPR)
          NUMI=POP(IRREPI,ISPIN)
          NUMA=VRT(IRREPA,ISPIN)
          NUMO=POP(IRREPO,3-ISPIN)
          NROW=NUMA
          NCOL=NUMI
          NSUM=DISSYW*NUMO
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &               ICORE(IOFFQ),NSUM,ONE,Y1(IOFFY),NROW)
          IOFFW=IOFFW+IINTFP*NROW*NSUM
          IOFFQ=IOFFQ+IINTFP*NCOL*NSUM
          IOFFY=IOFFY+IINTFP*NROW*NCOL
23       CONTINUE
C
21      continue
        IOFFY0=IOFFY0+IINTFP*NT(ISPIN)
C
20     CONTINUE
C
      ENDIF
C
      RETURN
      END
