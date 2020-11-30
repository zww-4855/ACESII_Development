      SUBROUTINE GW5INRHF(ICORE,MAXCOR,IUHF,IRREPT,LZOFF,LWOFF,LTOFF)
C
C THIS ROUTINE CALCULATES THE INITIAL T2*W CONTRIBUTION TO THE
C W5 HBAR INTERMEDIATES FOR RHF CASES ONLY.
C
C Modified version to allow T2 (hence target) not to be totally symmetric.
C W is currently assumed to be totally symmetric.
C
C New arguments :
C                IRREPT,LZOFF,LWOFF,LTOFF
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL TERM1,TERM2
      DOUBLE PRECISION ONE,ONEM,TWO,ZILCH,HALF,snrm2
      CHARACTER*3 GETTYP
      DIMENSION ICORE(MAXCOR),ILOCT(8),ILOCW(8,8),NDIMW(8)
      DIMENSION POPDUM(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      DATA ONE,ONEM,ZILCH,TWO,HALF /1.0D0,-1.0D0,0.0D0,2.0D0,0.5D0/
C
      LISTZ = LZOFF + 4
C
C BEGIN BY REORDERING <Ce|Am> => <Ae|Cm> ON DISK
C
C     LISTW=30
      LISTW=LWOFF + 4
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,LISTW,LISTW,1,
     &               .TRUE.,.FALSE.,.FALSE.)
      CALL IZERO(POPDUM,8)
      POPDUM(1)=1
C
C CALCULATE SIZES OF W(Ae,m) ARRAYS FOR EACH SYMMETRY OF C   
C AND OFFSETS INTO REORDERED W(em,A) ARRAY
C
      DO 5 IRREPC=1,NIRREP
       NDIMW(IRREPC)=0
       INCREM2=0
       DO 6 IRREPM=1,NIRREP
        IRREPAE=DIRPRD(IRREPM,IRREPC)
        INCREM=POP(IRREPM,1)*IRPDPD(IRREPAE,19)
        ILOCW(IRREPM,IRREPC)=IINTFP*INCREM2
        INCREM2=VRT(IRREPM,1)*IRPDPD(IRREPAE,9)+INCREM2
        NDIMW(IRREPC)=NDIMW(IRREPC)+INCREM
6      CONTINUE
5     CONTINUE
C
C READ IN T VECTOR AS T(em,bi) AND SPIN ADAPT IT 
C
      NSIZET=IDSYMSZ(IRREPT,ISYTYP(1,LTOFF + 37),ISYTYP(2,LTOFF + 37))
      I000=1
      I010=I000+IINTFP*NSIZET
      IOFFT=1
      DO 7 IRREP=1,NIRREP
       ILOCT(IRREP)=IOFFT
       DISSYT=IRPDPD(DIRPRD(IRREPT,IRREP),ISYTYP(1,LTOFF + 37))
       NUMDST=IRPDPD(       IRREP        ,ISYTYP(2,LTOFF + 37))
C      NUMDST=DISSYT
       CALL GETLST(ICORE(IOFFT),1,NUMDST,1,IRREP,LTOFF + 37)
       CALL GETLST(ICORE(I010),1,NUMDST,1,IRREP,LTOFF + 39)
       CALL SSCAL (NUMDST*DISSYT,TWO,ICORE(IOFFT),1)
       CALL SAXPY (NUMDST*DISSYT,ONEM,ICORE(I010),1,ICORE(IOFFT),1)
       IOFFT=IOFFT+IINTFP*NUMDST*DISSYT 
7     CONTINUE
C      
C FIRST TERM:
C
C   Z(Ab,Ci) <= [2 * T(im,eb) - T(mi,eb)] * <Ce|Am>
C
C INTEGRALS STORED NOW AS W(Ae,Cm)
C T VECTOR HELD IN CORE AS T(em,bi)
C
C LOOP OVER IRREPS OF TARGET
C
      DO 10 IRREPCI=1,NIRREP
       DO 11 IRREPI=1,NIRREP
        IRREPC=DIRPRD(IRREPI,IRREPCI)
        NUMI=POP(IRREPI,1)
        NUMC=VRT(IRREPC,1)
        I020=I010+IINTFP*NDIMW(IRREPC)
        I030=I020+IINTFP*NDIMW(IRREPC)
        DO 12 C=1,NUMC
C
C READ IN ALL W(Ae;m) FOR THIS VALUE OF c
C
         CALL GET3(ICORE(I020),LISTW,1,C,IRREPC,VRT,VRT,
     &             VRT,POP,19,9,'124',.FALSE.,.FALSE.,ICORE(I020))
C
C REORDER TO W(em;A)
C
         CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),VRT,
     &               VRT,POP,POPDUM,ICORE(I030),IRREPC,'2314') 
C
         DO 13 I=1,NUMI
C
          IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,9)-1
          CALL GETLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
C
          DO 14 IRREPB=1,NIRREP
           IRREPBI=DIRPRD(IRREPB,IRREPI)
           IRREPEM=DIRPRD(IRREPT,IRREPBI)
Cold       IRREPA=DIRPRD(IRREPB,IRREPCI)
Cnew
           IRREPAB=DIRPRD(IRREPT,IRREPCI)
           IRREPA=DIRPRD(IRREPB,IRREPAB)
C
           NUMB=VRT(IRREPB,1)
           NUMA=VRT(IRREPA,1)
           NROW=NUMA
           NCOL=NUMB
           NSUM=IRPDPD(IRREPEM,9)
C           IOFFT=ILOCT(IRREPEM)+IINTFP*
C     &           (NSUM*(ISYMOFF(IRREPI,IRREPEM,9)-1)+NSUM*NUMB*(I-1))
           IOFFT=ILOCT(IRREPBI)+IINTFP*
     &           (NSUM*(ISYMOFF(IRREPI,IRREPBI,9)-1)+NSUM*NUMB*(I-1))
           IOFFW=I010+ILOCW(IRREPA,IRREPC)
           IOFFZ=I020+IINTFP*(ISYMOFF(IRREPB,IRREPAB,19)-1)
C
C FORM Z(ab) = W(em,a) * T(em,bi) FOR FIXED i
C
           CALL XGEMM('T','N',NROW,NCOL,NSUM,HALF,ICORE(IOFFW),NSUM,
     &                ICORE(IOFFT),NSUM,ONE,ICORE(IOFFZ),NUMA)
14        CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD
C          
          IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,9)-1
          CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
13       CONTINUE
12      CONTINUE
11     CONTINUE
10    CONTINUE 
C
C END OF FIRST CONTRACTION
C
C REORDER INTEGRALS ON DISK AGAIN W(Ae,Cm)=>W(Ce,Am) AND UNSPIN-ADAPT
C THEM
C
C     LISTW=30
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,LISTW,LISTW,1,
     &               .FALSE.,.FALSE.,.FALSE.)
c
c temporary code to unspinadapt integrals
c
      DO 1000 IRREPCI=1,NIRREP
       DISSYW=IRPDPD(IRREPCI,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPCI,ISYTYP(2,LISTW))
       DO 1001 IDIS=1,NUMDSW
        CALL GETLST(ICORE(I000),IDIS,1,1,IRREPCI,LISTW)
        ITMP=I000+IINTFP*DISSYW
        CALL SYMTRA(IRREPCI,VRT,VRT,1,ICORE(I000),ICORE(ITMP))
        CALL SAXPY (DISSYW,TWO,ICORE(I000),1,ICORE(ITMP),1)
        CALL SSCAL (DISSYW,1.0D0/3.0D0,ICORE(ITMP),1)
        CALL PUTLST(ICORE(ITMP),IDIS,1,1,IRREPCI,LISTW)
1001   CONTINUE
1000  CONTINUE
C
C NOW REORDER INTEGRALS W(Ae,Cm) => W(AC,em)
C      
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,LISTW,LISTW,2,
     &               .FALSE.,.FALSE.,.FALSE.)
C   
C SECOND CONTRACTION: T(Im,Ea)*<Ec|Bm>     
C
C READ IN T VECTOR AS T(Em,aI) 
C
      NSIZET=IDSYMSZ(IRREPT,ISYTYP(1,LTOFF + 39),ISYTYP(2,LTOFF + 39))
      I000=1
      I010=I000+IINTFP*NSIZET
      CALL GETALL(ICORE(I000),NSIZET,IRREPT,LTOFF + 39)
C      
C
C   Z(Ab,Ci) <= T(Em,aI) * <Ec|Bm>
C
C INTEGRALS STORED NOW AS W(EB,cm)
C T VECTOR HELD IN CORE AS T(Em,aI)
C
C LOOP OVER IRREPS OF TARGET
C
      DO 110 IRREPCI=1,NIRREP
       DO 111 IRREPI=1,NIRREP
        IRREPC=DIRPRD(IRREPI,IRREPCI)
        NUMI=POP(IRREPI,1)
        NUMC=VRT(IRREPC,1)
        I020=I010+IINTFP*NDIMW(IRREPC)
        I030=I020+IINTFP*NDIMW(IRREPC)
        DO 112 C=1,NUMC
C
C READ IN ALL W(EB;m) FOR THIS VALUE OF c
C
         CALL GET3(ICORE(I020),LISTW,1,C,IRREPC,VRT,VRT,
     &             VRT,POP,19,9,'124',.FALSE.,.FALSE.,ICORE(I020))
C
C REORDER TO W(Em;B)
C
         CALL SSTGEN(ICORE(I020),ICORE(I010),NDIMW(IRREPC),VRT,
     &               VRT,POP,POPDUM,ICORE(I030),IRREPC,'1324')
C
         DO 113 I=1,NUMI
          LENAB=IRPDPD(DIRPRD(IRREPT,IRREPCI),19)
          DO 114 IRREPB=1,NIRREP
           IRREPAB=DIRPRD(IRREPT,IRREPCI)
           IRREPA=DIRPRD(IRREPB,IRREPAB)
           IRREPBI=DIRPRD(IRREPB,IRREPI)
           IRREPAI=DIRPRD(IRREPA,IRREPI)
           IRREPEM=DIRPRD(IRREPT,IRREPAI)
           NUMB=VRT(IRREPB,1)
           NUMA=VRT(IRREPA,1)
           NROW=NUMB
           NCOL=NUMA
           NSUM=IRPDPD(IRREPEM,9)
           IOFFT=ILOCT(IRREPAI)+IINTFP*
     &           (NSUM*(ISYMOFF(IRREPI,IRREPAI,9)-1)+NSUM*NUMA*(I-1))
           IOFFW=I010+ILOCW(IRREPB,IRREPC)
           IOFFZ=I020+IINTFP*(ISYMOFF(IRREPA,IRREPAB,19)-1)
C
C FORM Z(ba) = W(em,b) * T(em,ai) FOR FIXED i
C
           CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFW),NSUM,
     &                ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFZ),NUMB)
114       CONTINUE
C
C CALCULATE ADDRESS OF CI RECORD AND FORM 1/2 Z(ba)+Z(ab)
C          
          ITMP=I020+IINTFP*LENAB
          CALL SYMTRA(IRREPAB,VRT,VRT,1,ICORE(I020),ICORE(ITMP))
          CALL SAXPY (LENAB,HALF,ICORE(I020),1,ICORE(ITMP),1)
          IREC=C+NUMC*(I-1)+ISYMOFF(IRREPI,IRREPCI,9)-1
          CALL GETLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
C
C         CALL SAXPY (LENAB,ONEM,ICORE(I020),1,ICORE(ITMP),1)
C try this to get signs straight --- see negative sign in commented
C code below.
          CALL SAXPY (LENAB,ONEM,ICORE(ITMP),1,ICORE(I020),1)
C
          CALL PUTLST(ICORE(I020),IREC,1,1,IRREPCI,LISTZ)
113      CONTINUE
112     CONTINUE
111    CONTINUE
c
c for debugging ...
c
c       numdis=irpdpd(irrepci,isytyp(2,30))
c       dissiz=irpdpd(irrepci,isytyp(1,30))
c       call getlst(icore(i010),1,numdis,1,irrepci,130+iofflst)
c       write(6,*)' checksum for irrep ',irrepci,' is ',
c     &            snrm2(numdis*dissiz,icore(i010),1)
110   CONTINUE 
C
C NOW REORDER INTEGRALS W(AC,em) => W(Ae,Cm)
C      
      CALL SYMTRLST(ICORE,MAXCOR,IUHF,VRT,VRT,VRT,POP,
     &               1,19,19,19,9,9,LISTW,LISTW,2,
     &               .FALSE.,.FALSE.,.FALSE.)
C
C NOW ADD LISTS 30 AND 130
C
c      DO 2000 IRREP=1,NIRREP
c       I000=1
c       DISSIZ=IRPDPD(IRREP,ISYTYP(1,30))
c       NUMDIS=IRPDPD(IRREP,ISYTYP(2,30))
c       IF(DISSIZ.NE.0)THEN
c        NINCOR=MAXCOR/(DISSIZ*IINTFP)
c       ELSE
c        NINCOR=2*NUMDIS
c       ENDIF
c       NBUNCH=MIN(NINCOR/2,NUMDIS)
c       I010=I000+IINTFP*DISSIZ*NBUNCH
c       ISTART=1
c       NLEFT=NUMDIS
c1      NREAD=MIN(NLEFT,NBUNCH)
c       IF(TERM1)THEN
c        CALL GETLST(ICORE(I000),ISTART,NREAD,1,IRREP,30)
c       ELSE
c        CALL ZERO  (ICORE(I000),NREAD*DISSIZ)
c       ENDIF
c       IF(TERM2)THEN
c        CALL GETLST(ICORE(I010),ISTART,NREAD,1,IRREP,130+IOFFLIST)
c       ELSE
c        CALL ZERO  (ICORE(I010),NREAD*DISSIZ)
c       ENDIF
c       CALL SAXPY (NREAD*DISSIZ,ONEM,ICORE(I010),1,ICORE(I000),1)
c       CALL PUTLST(ICORE(I000),ISTART,NREAD,1,IRREP,130+IOFFLIST)
c       NLEFT=NLEFT-NREAD
c       ISTART=ISTART+NREAD
c       IF(NLEFT.NE.0)GOTO 1
c2000  CONTINUE
C
      RETURN
      END
