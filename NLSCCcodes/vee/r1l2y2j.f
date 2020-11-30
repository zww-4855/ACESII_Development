      SUBROUTINE R1L2Y2J(F1,R1,ICORE,MAXCOR,IUHF)
C
C Y1(ab,ij) = R(me)*F(ma)*L(eb,ij)
C  
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR),R1(*),F1(*),I0X(2),I0F(2),I0R(2)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      DATA ONE,ONEM,ZILCH/1.0D0,-1.0D0,0.0D0/
C
C CALCULATE X(ea) = R(me)*F(ma) INTERMEDIATES
C
      I0X(1)=1
      I0F(1)=1
      I0R(1)=1
      IF(IUHF.NE.0)THEN
       I0X(2)=I0X(1)+IINTFP*IRPDPD(IRREPX,19)
       I0F(2)=I0F(1)+IINTFP*NT(1)
       I0R(2)=I0R(1)+IINTFP*IRPDPD(IRREPX,9)
      ELSE
       I0X(2)=I0X(1)
       I0F(2)=I0F(1)
       I0R(2)=I0R(1)
      ENDIF
      I000=I0X(2)+IINTFP*IRPDPD(IRREPX,20)
      DO 5 ISPIN=1,1+IUHF
       DO 6 IRREPA=1,NIRREP
        IRREPE=DIRPRD(IRREPA,IRREPX)
        IRREPM=IRREPA
        NUMA=VRT(IRREPA,ISPIN)
        NUME=VRT(IRREPE,ISPIN)
        NUMM=POP(IRREPM,ISPIN)
        NROW=NUME
        NCOL=NUMA
        NSUM=NUMM
        IOFFR=I0R(ISPIN)+IINTFP*(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
        IOFFF=I0F(ISPIN)+IINTFP*(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IOFFX=I0X(ISPIN)+IINTFP*(ISYMOFF(IRREPA,IRREPX,18+ISPIN)-1)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONE,R1(IOFFR),NROW,
     &             F1(IOFFF),NCOL,ZILCH,ICORE(IOFFX),NROW)
6      CONTINUE
5     CONTINUE
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE.
C
C   Y(Ab,Ij) = Q(Ab,Ij) + Q(Ba,Ji)
C
C WHERE
C
C   Q(Ab,Ij) = R(ME)*F(MA)*L(Eb,Ij)
C
C
C READ IN L(Eb,Ij) = L(Be,Ji) AND TRANSPOSE TO L(Ji,Be)
C
       DO 10 IRREPR=1,NIRREP
        LISTL=446
        LISTZ=63 
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        MAXT  =MAX(DISSYL,NUMDSL,DISSYZ,NUMDSZ)
        I010=I000+IINTFP*MAX(DISSYL*NUMDSL,3*MAXT,DISSYZ*NUMDSZ)
        I020=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYL*NUMDSL)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
C
C  Q(Ji,Ba) = L(Ji,Be) * X(ea)
C
        DO 20 IRREPA=1,NIRREP
         IRREPE=DIRPRD(IRREPX,IRREPA)
         IRREPB=DIRPRD(IRREPE,IRREPL)
         NUMA=VRT(IRREPA,1)
         NUME=VRT(IRREPE,1)
         NUMB=VRT(IRREPB,1)
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUME
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPE,IRREPL,19)-1)*NUMDSL
         IOFFX=I0X(1)+IINTFP*(ISYMOFF(IRREPA,IRREPX,19)-1)
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPA,IRREPR,19)-1)*NUMDSZ
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NSUM,ZILCH,ICORE(IOFFZ),NROW)
20      CONTINUE
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
        CALL SCOPY (NUMDSZ*DISSYZ,ICORE(I000),1,ICORE(I010),1)
        CALL SYMRHF(IRREPR,VRT(1,1),POP(1,1),DISSYZ,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPR,LISTZ)
C
10     CONTINUE
C
      ELSE
C
C UHF
C
C AAAA AND BBBB
C
C       Y1(AB,IJ) = R(ME)*F(MA)*L(EB,IJ) + ANTISYMMETRIZATION
C

       DO 100 ISPIN=1,2
        DO 110 IRREPR=1,NIRREP
         LISTL=443+ISPIN
         LISTZ=60+ISPIN
         IRREPL=DIRPRD(IRREPR,IRREPX)
         DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
         DISSYLX=IRPDPD(IRREPL,18+ISPIN)
         NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
         DISSYZX=IRPDPD(IRREPR,18+ISPIN)
         DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
         NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
         MAXT  =MAX(DISSYLX,NUMDSL,DISSYZX,NUMDSZ)
         I010=I000+IINTFP*MAX(DISSYLX*NUMDSL,3*MAXT,DISSYZX*NUMDSZ)
         I020=I010+IINTFP*MAX(DISSYZX*NUMDSZ,DISSYLX*NUMDSL)
         ITMP1=I020
         ITMP2=ITMP1+IINTFP*MAXT
         ITMP3=ITMP2+IINTFP*MAXT
         CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
         CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
         CALL SYMEXP(IRREPL,VRT(1,ISPIN),NUMDSL,ICORE(I000))
C
C  -Y1(I<J,AB) = Y1(I<J,BA) = L(I<J,BE)*X(EA)
C
        DO 120 IRREPA=1,NIRREP
         IRREPE=DIRPRD(IRREPX,IRREPA)
         IRREPB=DIRPRD(IRREPE,IRREPL)
         NUMA=VRT(IRREPA,ISPIN)
         NUME=VRT(IRREPE,ISPIN)
         NUMB=VRT(IRREPB,ISPIN)
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUME
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPE,IRREPL,18+ISPIN)-1)*NUMDSL
         IOFFX=I0X(ISPIN)+IINTFP*(ISYMOFF(IRREPA,IRREPX,18+ISPIN)-1)
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPA,IRREPR,18+ISPIN)-1)*NUMDSZ
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NSUM,ZILCH,ICORE(IOFFZ),NROW)
120      CONTINUE
         CALL ASSYM2(IRREPR,VRT(1,ISPIN),NUMDSZ,ICORE(I010))
         CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
         CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
         CALL SAXPY (NUMDSZ*DISSYZ,ONEm,ICORE(I000),1,ICORE(I010),1)
         CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
C
110     CONTINUE
C
100    CONTINUE
C
C ABAB : Y1(Ab,Ij) = L(Eb,Ij)*X(EA) + L(Ae,Ij)*X(eb)
C
       DO 210 IRREPR=1,NIRREP
        LISTL=446
        LISTZ=63
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DISSYL=IRPDPD(IRREPL,ISYTYP(1,LISTL))
        NUMDSL=IRPDPD(IRREPR,ISYTYP(2,LISTL))
        DISSYZ=IRPDPD(IRREPR,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPR,ISYTYP(2,LISTZ))
        MAXT  =MAX(DISSYL,NUMDSL,DISSYZ,NUMDSZ)
        I010=I000+IINTFP*MAX(DISSYL*NUMDSL,3*MAXT,DISSYZ*NUMDSZ)
        I020=I010+IINTFP*MAX(DISSYZ*NUMDSZ,DISSYL*NUMDSL)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXT
        ITMP3=ITMP2+IINTFP*MAXT
        CALL GETLST(ICORE(I010),1,NUMDSL,1,IRREPR,LISTL)
        CALL TRANSP(ICORE(I010),ICORE(I000),NUMDSL,DISSYL)
        CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMDSL,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C  Y1(Ij,bA) <= L(Ij,bE)*X(EA)
C
        DO 220 IRREPA=1,NIRREP
         IRREPE=DIRPRD(IRREPX,IRREPA)
         IRREPB=DIRPRD(IRREPE,IRREPL)
         NUMA=VRT(IRREPA,1)
         NUME=VRT(IRREPE,1)
         NUMB=VRT(IRREPB,2)
         NROW=NUMDSL*NUMB
         NCOL=NUMA
         NSUM=NUME
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPE,IRREPL,23)-1)*NUMDSL
         IOFFX=I0X(1)+IINTFP*(ISYMOFF(IRREPA,IRREPX,19)-1)
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPA,IRREPR,23)-1)*NUMDSZ
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NSUM,ZILCH,ICORE(IOFFZ),NROW)
220      CONTINUE
C
C Y1(Ij,bA) -> Y1(Ij,Ab)
C L (Ij,eA) -> L (Ij,Ae)
C
        CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMDSL,ICORE(I000),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
        CALL SYMTR1(IRREPR,VRT(1,2),VRT(1,1),NUMDSZ,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C  Y1(Ij,Ab) <= L(Ij,Ae)*X(eb)
C
        DO 230 IRREPB=1,NIRREP
         IRREPE=DIRPRD(IRREPX,IRREPB)
         IRREPA=DIRPRD(IRREPE,IRREPL)
         NUMA=VRT(IRREPA,1)
         NUME=VRT(IRREPE,2)
         NUMB=VRT(IRREPB,2)
         NROW=NUMDSL*NUMA
         NCOL=NUMB
         NSUM=NUME
         IOFFL=I000+IINTFP*(ISYMOFF(IRREPE,IRREPL,13)-1)*NUMDSL
         IOFFX=I0X(2)+IINTFP*(ISYMOFF(IRREPB,IRREPX,20)-1)
         IOFFZ=I010+IINTFP*(ISYMOFF(IRREPB,IRREPR,13)-1)*NUMDSZ
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFX),NSUM,ONE,ICORE(IOFFZ),NROW)
230     CONTINUE
C
        CALL TRANSP(ICORE(I010),ICORE(I000),DISSYZ,NUMDSZ)
        CALL GETLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
        CALL SAXPY (DISSYZ*NUMDSZ,ONEM,ICORE(I000),1,ICORE(I010),1)
        CALL PUTLST(ICORE(I010),1,NUMDSZ,1,IRREPR,LISTZ)
C
210    CONTINUE 
C
      ENDIF
C
      RETURN
      END
