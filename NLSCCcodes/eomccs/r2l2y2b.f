      SUBROUTINE R2L2Y2B(ICORE,MAXCOR,IUHF)
C
C Y2(ab,ij) = R(mn,ef)*L(ij,ef)*W(ab,mn)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH
      DIMENSION ICORE(MAXCOR)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
      DO 1 ISPIN=3,3-2*IUHF,-1 
       LISTL=443+ISPIN
       LISTR=460+ISPIN
       LISTW=13+ISPIN
       LISTZ=60+ISPIN
       DO 10 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        DISSYW=DISSYZ
        NUMDSW=NUMDSZ
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYR=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        NUMDSR=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        MAXSIZ=MAX(DISSYZ*NUMDSZ,DISSYL*NUMDSL,DISSYR*NUMDSR)
        I000=1
        I010=I000+IINTFP*MAXSIZ
        I020=I010+IINTFP*MAXSIZ
        I030=I020+IINTFP*MAXSIZ
        CALL GETLST(ICORE(I000),1,NUMDSL,1,IRREPR,LISTL)
        CALL GETLST(ICORE(I010),1,NUMDSR,1,IRREPR,LISTR)
C
C                    +
C Q(mn,ij) = R(ef,mn) * L(ef,ij)
C
        CALL XGEMM('T','N',NUMDSR,NUMDSL,DISSYL,ONE,ICORE(I010),
     &             DISSYR,ICORE(I000),DISSYL,ZILCH,ICORE(I020),
     &             NUMDSR)
C
        CALL GETLST(ICORE(I000),1,NUMDSW,1,IRREPR,LISTW)
C                    
C Z(ab,ij) = W(ab,mn) * Q(mn,ij)
C
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL XGEMM('N','N',DISSYZ,NUMDSZ,NUMDSW,ONE,ICORE(I000),
     &             DISSYW,ICORE(I020),NUMDSW,ONE,ICORE(I010),DISSYZ)
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
10     CONTINUE
1     CONTINUE
C
      RETURN
      END
