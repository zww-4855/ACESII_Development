      SUBROUTINE S2EXPECT(ICORE,MAXCOR,IUHF,RNGAFTER,READDEN,
     &                    IRD)
C
C CALCULATES THE GENERALIZED EXPECTATION VALUE OF S**2 BY CONTRACTING
C MATRIX ELEMENTS OF S**2 WITH ELEMENTS OF THE ONE AND TWO-PARTICLE
C DENSITY MATRICES.  
C
C <A+|S**2|B> = SAMP(PQ)*D(PQ) + SAMP(PQRS)*D(PQRS)
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL READDEN
      CHARACTER*1 RNGAFTER
      DOUBLE PRECISION ONE,ZILCH,TWO,SDOT,FACT1,HALF
      DOUBLE PRECISION CIJKL,CIJKA,CIJAB,CAIBJ,CABCI,CABCD,TOT2PDM
      DOUBLE PRECISION CIJ,CAB,CAI,TOT1PDM,DIFFO,S200,SUMO
      DOUBLE PRECISION FABCD,FIJKL,FABCI,FIJKA,FABIJ,FAIBJ
      LOGICAL PRINT
      DIMENSION ICORE(MAXCOR),IDDOO(2),IDDVV(2),IDDVO(2)
      DIMENSION IDNOO(2),IDNVV(2),IDNVO(2),IDNOV(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA HALF,ONE,ZILCH,TWO /0.5D0,1.0D0,0.0D0,2.0D0/
C
      PRINT=IFLAGS(1).GE.4
C
      WRITE(6,900)
900   FORMAT(T3,'@S2EXPECT-I, S**2 calculated as expectation value.')
C
C CREATE LISTS
C
      DO 3 ISPIN=1,1+IUHF
       CALL UPDMOI(1,NFMI(ISPIN),2+ISPIN,91,0,0)
       CALL UPDMOI(1,NFEA(ISPIN),2+ISPIN,92,0,0)
3     CONTINUE
C
C INITIALIZE S**2
C
      FACT1=ONE
      FIJKL=ONE
      FABCD=ONE
      FABIJ=TWO
      FAIBJ=ONE
      FIJKA=TWO
      FABCI=TWO
C
      CALL CHKRNG(ICORE,MAXCOR,IUHF,'N')
C
C DO ONE-PARTICLE CONTRIBUTIONS
C
C   D(AI)*DEL(Ae)*DEL(eI) + D(ai)*DEL(aM)*DEL(Mi) 
C
C + D(IJ)*DEL(Ie)*DEL(eJ) + D(ij)*DEL(iM)*DEL(Mi)
C
C + D(AB)*DEL(Ae)*DEL(eB) + D(ab)*DEL(aM)*DEL(Mb)
C
      IDDOO(1)=1
      IDDOO(2)=IDDOO(1)+IINTFP*NFMI(1)
      IDDVV(1)=IDDOO(2)+IINTFP*NFMI(2)
      IDDVV(2)=IDDVV(1)+IINTFP*NFEA(1)
      IDDVO(1)=IDDVV(2)+IINTFP*NFEA(2)
      IDDVO(2)=IDDVO(1)+IINTFP*NT(1)
      ITOP    =IDDVO(2)+IINTFP*NT(2)
      LENOOAB =IRPDPD(1,14)
      LENVVAB =IRPDPD(1,13)
      LENVOAB =IRPDPD(1,11)
      LENOVAB =IRPDPD(1,18)
      IDOOAB  =ITOP
      IDVVAB  =IDOOAB+IINTFP*LENOOAB
      IDVOAB  =IDVVAB+IINTFP*LENVVAB
      IDOVAB  =IDVOAB+IINTFP*LENVOAB
      ITOP    =IDOVAB+IINTFP*LENOVAB
      IDNOO(1)=ITOP
      IDNOO(2)=IDNOO(1)+IINTFP*NFMI(1)
      IDNVV(1)=IDNOO(2)+IINTFP*NFMI(2)
      IDNVV(2)=IDNVV(1)+IINTFP*NFEA(1)
      IDNVO(1)=IDNVV(2)+IINTFP*NFEA(2)
      IDNVO(2)=IDNVO(1)+IINTFP*NT(1)
      IDNOV(1)=IDNVO(2)+IINTFP*NT(2)
      IDNOV(2)=IDNOV(1)+IINTFP*NT(1)
      ITOP    =IDNOV(2)+IINTFP*NT(2)
      MXCOR   =MAXCOR-ITOP+1
      CALL GETLST(ICORE(IDOOAB),1,1,1,5,90)
      CALL GETLST(ICORE(IDVVAB),1,1,1,8,90)
      CALL GETLST(ICORE(IDVOAB),1,1,1,6,90)
      CALL GETLST(ICORE(IDOVAB),1,1,1,7,90)
C
C EVALUATE GROUND STATE DENSITY
C
      IF(READDEN)THEN
       DO 2 ISPIN=1,2
        CALL GETLST(ICORE(IDNOO(ISPIN)),IRD,1,1,ISPIN,160)
        CALL GETLST(ICORE(IDNVV(ISPIN)),IRD,1,1,ISPIN+2,160)
        CALL GETLST(ICORE(IDNVO(ISPIN)),IRD,1,1,ISPIN+4,160)
2      CONTINUE
      ELSE
       CALL GSDEN(ICORE(IDNOO(1)),ICORE(IDNVV(1)),ICORE(IDNVO(1)),
     &            ICORE(IDNOV(1)),0,190,0,90,144,44,ONE,ICORE(ITOP),
     &            MXCOR,IUHF)
       DO 1 ISPIN=1,2
        LAI=NT(ISPIN)
        CALL SAXPY(LAI,ONE,ICORE(IDNOV(ISPIN)),1,ICORE(IDNVO(ISPIN)),1)
        CALL SSCAL(LAI,HALF,ICORE(IDNVO(ISPIN)),1)
1      CONTINUE
      ENDIF
C
C EVALUATE REFERENCE STATE MULTIPLICITY
C
      NALPHA=NOCCO(1)
      NBETA =NOCCO(2)
      SUMO  =DFLOAT(NALPHA+NBETA)
      DIFFO =DFLOAT(NALPHA-NBETA)
      S200  =HALF*(SUMO+HALF*DIFFO*DIFFO)
     &        -SDOT(LENOOAB,ICORE(IDOOAB),1,ICORE(IDOOAB),1)
C
C CALCULATE DEL(Ae)*DEL(Ie) = DD(AI)
C
      IOFF1=IDVVAB
      IOFF2=IDOVAB
      IOFF3=IDDVO(1)
      DO 10 IRREPI=1,NIRREP
       IRREPA=IRREPI
       IRREPE=IRREPI
       NUMA=VRT(IRREPA,1)
       NUMI=POP(IRREPI,1)
       NUME=VRT(IRREPE,2)
       NROW=NUMA
       NCOL=NUMI
       NSUM=NUME
       CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT1,ICORE(IOFF1),NROW,
     &            ICORE(IOFF2),NCOL,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
10    CONTINUE
C
C CALCULATE - DEL(Ma)*DEL(Mi) = DD(ai)
C
      IOFF1=IDOVAB
      IOFF2=IDOOAB
      IOFF3=IDDVO(2)
      DO 20 IRREPI=1,NIRREP
       IRREPA=IRREPI
       IRREPM=IRREPI
       NUMA=VRT(IRREPA,2)
       NUMI=POP(IRREPI,2)
       NUMM=POP(IRREPM,1)
       NROW=NUMA
       NCOL=NUMI
       NSUM=NUMM
       CALL XGEMM('T','N',NROW,NCOL,NSUM,-FACT1,ICORE(IOFF1),NSUM,
     &            ICORE(IOFF2),NSUM,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
20    CONTINUE
C
C CALCULATE DEL(Ae)*DEL(Be) = DD(AB)
C
      IOFF1=IDVVAB
      IOFF2=IDVVAB
      IOFF3=IDDVV(1)
      DO 30 IRREPA=1,NIRREP
       IRREPB=IRREPA
       IRREPE=IRREPA
       NUMA=VRT(IRREPA,1)
       NUMB=VRT(IRREPB,1)
       NUME=VRT(IRREPE,2)
       NROW=NUMA
       NCOL=NUMB
       NSUM=NUME
       CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT1,ICORE(IOFF1),NROW,
     &            ICORE(IOFF2),NCOL,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
30    CONTINUE
C
C CALCULATE DEL(Ma)*DEL(Mb) = D(ab)
C
      IOFF1=IDOVAB
      IOFF2=IDOVAB
      IOFF3=IDDVV(2)
      DO 40 IRREPB=1,NIRREP
       IRREPA=IRREPB
       IRREPE=IRREPA
       NUMA=VRT(IRREPA,2)
       NUMB=VRT(IRREPB,2)
       NUME=POP(IRREPE,1)
       NROW=NUMA
       NCOL=NUMB
       NSUM=NUME
       CALL XGEMM('T','N',NROW,NCOL,NSUM,-FACT1,ICORE(IOFF1),NSUM,
     &            ICORE(IOFF2),NSUM,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
40    CONTINUE
C
C CALCULATE DEL(Ie)*DEL(Je) = DD(IJ)
C
      IOFF1=IDOVAB
      IOFF2=IDOVAB
      IOFF3=IDDOO(1)
      DO 50 IRREPA=1,NIRREP
       IRREPB=IRREPA
       IRREPE=IRREPA
       NUMA=POP(IRREPA,1)
       NUMB=POP(IRREPB,1)
       NUME=VRT(IRREPE,2)
       NROW=NUMA
       NCOL=NUMB
       NSUM=NUME
       CALL XGEMM('N','T',NROW,NCOL,NSUM,FACT1,ICORE(IOFF1),NROW,
     &            ICORE(IOFF2),NCOL,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
50    CONTINUE
C
C CALCULATE DEL(Mi)*DEL(Mj) = DD(ij)
C
      IOFF1=IDOOAB
      IOFF2=IDOOAB
      IOFF3=IDDOO(2)
      DO 60 IRREPB=1,NIRREP
       IRREPA=IRREPB
       IRREPE=IRREPA
       NUMA=POP(IRREPA,2)
       NUMB=POP(IRREPB,2)
       NUME=POP(IRREPE,1)
       NROW=NUMA
       NCOL=NUMB
       NSUM=NUME
       CALL XGEMM('T','N',NROW,NCOL,NSUM,-FACT1,ICORE(IOFF1),NSUM,
     &            ICORE(IOFF2),NSUM,ZILCH,ICORE(IOFF3),NROW)
       IOFF1=IOFF1+IINTFP*NROW*NSUM
       IOFF2=IOFF2+IINTFP*NSUM*NCOL
       IOFF3=IOFF3+IINTFP*NROW*NCOL
60    CONTINUE
C
C CALCULATE ONE PARTICLE CONTRIBUTION TO S**2
C
      CIJ=ZILCH
      CAB=ZILCH
      CAI=ZILCH
      DO 5 ISPIN=1,2
       LIJ=NFMI(ISPIN)
       LAB=NFEA(ISPIN)
       LAI=NT(ISPIN)
       CIJ=CIJ+SDOT(LIJ,ICORE(IDNOO(ISPIN)),1,ICORE(IDDOO(ISPIN)),1)
       CAB=CAB+SDOT(LAB,ICORE(IDNVV(ISPIN)),1,ICORE(IDDVV(ISPIN)),1)
       CAI=CAI+TWO*SDOT(LAI,ICORE(IDNVO(ISPIN)),1,ICORE(IDDVO(ISPIN)),1)
5     CONTINUE
      TOT1PDM=CAI+CIJ+CAB
C
C NOW DO TWO-PARTICLE CONTRIBUTIONS.   
C
C IJKL 2PDM PART
C
      LISTW=113
      CALL ZERO(ICORE(ITOP),LENOOAB)
      DO 70 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENOOAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL DOT24 (IRREP,ICORE(I000),ICORE(IDOOAB),ICORE(I010),
     &             ICORE(ITMP),DISSIZ,POP(1,2),POP(1,1),
     &             POP(1,1),POP(1,2),POP(1,1),POP(1,2),'STTS')
70    CONTINUE
      CALL SYMTRA(1,POP(1,2),POP(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CIJKL=-FIJKL*SDOT(LENOOAB,ICORE(ITMP),1,ICORE(IDOOAB),1)
C
C IJKA 2PDM PART
C
C IjKa
C
      LISTW=110
      CALL ZERO(ICORE(ITOP),LENOOAB)
      DO 80 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENOOAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL DOT24 (IRREP,ICORE(I000),ICORE(IDOVAB),ICORE(I010),
     &             ICORE(ITMP),DISSIZ,POP(1,2),POP(1,1),
     &             POP(1,1),POP(1,2),POP(1,1),VRT(1,2),'STTS')
80    CONTINUE
      CALL SYMTRA(1,POP(1,2),POP(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CIJKA=-FIJKA*SDOT(LENOOAB,ICORE(ITMP),1,ICORE(IDOOAB),1)
C
C IjAk
C
      LISTW=109
      CALL ZERO(ICORE(ITOP),LENVOAB)
      DO 90 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENVOAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL DOT24 (IRREP,ICORE(I000),ICORE(IDOOAB),ICORE(I010),
     &             ICORE(ITMP),DISSIZ,POP(1,2),VRT(1,1),
     &             POP(1,1),POP(1,2),VRT(1,1),POP(1,2),'STTS')
90    CONTINUE
      CALL SYMTRA(1,POP(1,2),VRT(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CIJKA=CIJKA-FIJKA*SDOT(LENVOAB,ICORE(ITMP),1,ICORE(IDVOAB),1)  
C
C ABIJ 2PDM PART
C
      LISTW=116
      CALL ZERO(ICORE(ITOP),LENOVAB)
      DO 100 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENOVAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL DOT24 (IRREP,ICORE(I000),ICORE(IDVOAB),ICORE(I010),
     &             ICORE(ITMP),DISSIZ,VRT(1,2),POP(1,1),
     &             VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),'STTS')
100   CONTINUE
      CALL SYMTRA(1,VRT(1,2),POP(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CIJAB=-FABIJ*SDOT(LENOVAB,ICORE(ITMP),1,ICORE(IDOVAB),1)  
C
C ABCI 2PDM PART
C
C AbCi
C
      LISTW=130
      CALL ZERO(ICORE(ITOP),LENVVAB)
      DO 110 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENVVAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       ITOP2=ITMP+IINTFP*MAX(LENOVAB,LENVVAB)
       IF(ITOP2.LE.MXCOR)THEN
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        CALL DOT24 (IRREP,ICORE(I000),ICORE(IDVOAB),ICORE(I010),
     &              ICORE(ITMP),DISSIZ,VRT(1,2),VRT(1,1),
     &              VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),'STTS')
       ELSE  
        CALL DOT24X(IRREP,ICORE(I000),ICORE(IDVOAB),ICORE(I010),
     &              VRT(1,1),VRT(1,2),VRT(1,1),POP(1,2),LISTW,
     &              23,11,13,'STTS',1,1,1)
       ENDIF
110   CONTINUE
      ITMP=I000+IINTFP*LENVVAB
      CALL SYMTRA(1,VRT(1,2),VRT(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CABCI=-FABCI*SDOT(LENVVAB,ICORE(ITMP),1,ICORE(IDVVAB),1)  
C
C AbIc
C
      LISTW=129
      CALL ZERO(ICORE(ITOP),LENOVAB)
      DO 120 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENOVAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       ITOP2=ITMP+IINTFP*LENOVAB
       IF(ITOP2.LE.MXCOR)THEN
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        CALL DOT24 (IRREP,ICORE(I000),ICORE(IDVVAB),ICORE(I010),
     &              ICORE(ITMP),DISSIZ,VRT(1,2),POP(1,1),
     &              VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),'STTS')
       ELSE
        CALL DOT24X(IRREP,ICORE(I000),ICORE(IDVVAB),ICORE(I010),
     &              VRT(1,1),VRT(1,2),POP(1,1),VRT(1,2),LISTW,
     &              12,13,13,'STTS',1,1,1)
       ENDIF
120   CONTINUE
      ITMP=I000+IINTFP*LENOVAB
      CALL SYMTRA(1,VRT(1,2),POP(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CABCI=CABCI-FABCI*SDOT(LENOVAB,ICORE(ITMP),1,ICORE(IDOVAB),1)  
C
C ABCD 2PDM PART
C
      LISTW=133
      CALL ZERO(ICORE(ITOP),LENVVAB)
      DO 130 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENVVAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       ITOP2=ITMP+IINTFP*LENVVAB
       IF(ITOP2.LE.MXCOR)THEN
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        CALL DOT24 (IRREP,ICORE(I000),ICORE(IDVVAB),ICORE(I010),
     &              ICORE(ITMP),DISSIZ,VRT(1,2),VRT(1,1),
     &              VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),'STTS')
       ELSE
        CALL DOT24X(IRREP,ICORE(I000),ICORE(IDVVAB),ICORE(I010),
     &              VRT(1,1),VRT(1,2),VRT(1,1),VRT(1,2),LISTW,
     &              23,13,13,'STTS',1,1,1)
       ENDIF
130   CONTINUE
      ITMP=I000+IINTFP*LENVVAB
      CALL SYMTRA(1,VRT(1,2),VRT(1,1),1,ICORE(I000),ICORE(ITMP)) 
      CABCD=-FABCD*SDOT(LENVVAB,ICORE(ITMP),1,ICORE(IDVVAB),1)  
C
C RING TERMS.  THERE ARE TWO CONTRIBUTIONS
C
C D(Jb)*D(Ia)*GAMMA(Ib,aJ) - D(Jb)*D(Ia)*G(aI,bJ)
C
      LISTW=126
      DISSIZ=IRPDPD(1,ISYTYP(1,LISTW))
      NUMDIS=IRPDPD(1,ISYTYP(2,LISTW))
      I000=ITOP
      I010=I000+IINTFP*LENOVAB
      ITMP=I010+IINTFP*DISSIZ*NUMDIS
      CALL SYMTRA(1,POP(1,1),VRT(1,2),1,ICORE(IDOVAB),ICORE(ITMP)) 
      CALL GETLST(ICORE(I010),1,NUMDIS,1,1,LISTW)
      CALL XGEMM('N','N',1,NUMDIS,DISSIZ,ONE,ICORE(ITMP),1,
     &           ICORE(I010),DISSIZ,ZILCH,ICORE(I000),1)
      CAIBJ=-FAIBJ*SDOT(LENOVAB,ICORE(ITMP),1,ICORE(I000),1)
C
C D(Ai)*D(Bj)*GAMMA(iB,Aj) - D(Ai)*D(Bj)*G(Ai,Bj)
C
      LISTW=125
      DISSIZ=IRPDPD(1,ISYTYP(1,LISTW))
      NUMDIS=IRPDPD(1,ISYTYP(2,LISTW))
      I000=ITOP
      I010=I000+IINTFP*LENVOAB
      ITMP=I010+IINTFP*DISSIZ*NUMDIS
      CALL GETLST(ICORE(I010),1,NUMDIS,1,1,LISTW)
      CALL XGEMM('N','N',1,NUMDIS,DISSIZ,ONE,ICORE(IDVOAB),1,
     &           ICORE(I010),DISSIZ,ZILCH,ICORE(I000),1)
      CAIBJ=CAIBJ-FAIBJ*SDOT(LENVOAB,ICORE(IDVOAB),1,ICORE(I000),1)
C
C D(Ab)*D(Ij)*GAMMA(Ib,Aj) - D(Ab)*D(Ij)*G(AI,bj)
C
      LISTW=118
      CALL ZERO(ICORE(ITOP),LENOOAB)
      DO 140 IRREP=1,NIRREP
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       I000=ITOP
       I010=I000+IINTFP*LENOOAB
       ITMP=I010+IINTFP*DISSIZ*NUMDIS
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL DOT24 (IRREP,ICORE(I000),ICORE(IDVVAB),ICORE(I010),
     &             ICORE(ITMP),DISSIZ,POP(1,1),POP(1,2),
     &             VRT(1,1),POP(1,1),VRT(1,2),POP(1,2),'STST')
140   CONTINUE
      CAIBJ=CAIBJ-TWO*FAIBJ*SDOT(LENOOAB,ICORE(I000),1,ICORE(IDOOAB),1)
C
      TOT2PDM=CABCD+CABCI+CIJAB+CAIBJ+CIJKA+CIJKL
C
      IF(PRINT)THEN
       WRITE(6,1000)
1000   FORMAT(T3,'@S2EXPECT-I, Contributions to <S**2> : ')
       WRITE(6,1001)
1001   FORMAT(T3,'From one particle density : ')
       WRITE(6,1002)CAI,CIJ,CAB
1002   FORMAT(T5,'D(ai)   : ',F15.10,/,T5,'D(ij)   : ',F15.10,/,
     &        T5,'D(ab)   : ',F15.10)  
       WRITE(6,1003)
1003   FORMAT(T3,'From two particle density : ')
       WRITE(6,1004)CABCD,CABCI,CIJAB,CAIBJ,CIJKA,CIJKL
1004   FORMAT(T5,'G(abcd) : ',F15.10,/,T5,'G(abci) : ',F15.10,/,
     &        T5,'G(abij) : ',F15.10,/,T5,'G(aibj) : ',F15.10,/,
     &        T5,'G(ijka) : ',F15.10,/,T5,'G(ijkl) : ',F15.10)
      ENDIF
C
      WRITE(6,1005)S200
1005  FORMAT(T15,'Reference state <S**2>        : ',F15.10,'.')
      WRITE(6,1006)TOT1PDM,TOT2PDM
1006  FORMAT(T15,'Total one particle correction : ',F15.10,'.',
     &       /,T15,'Total two particle correction : ',F15.10,'.')
      WRITE(6,1007)TOT1PDM+TOT2PDM+S200
1007  FORMAT(T15,'Expectation value of <S**2>   : ',F15.10,'.')
C
      CALL CHKRNG(ICORE,MAXCOR,IUHF,RNGAFTER)
C
      RETURN
      END