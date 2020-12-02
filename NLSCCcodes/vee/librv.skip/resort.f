      SUBROUTINE RESORT(ICORE,MAXCOR,IUHF,IRREPX,LISTT2,LISTT2RS)
C
C FORM RESORTED T2 LISTS
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 SPCASE(2),SPCAS2
      DIMENSION ICORE(MAXCOR),LIST(6)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      LISTT2AB=LISTT2-1+3
C
      IF(IUHF.EQ.0)THEN
C
C RHF RESORTS
C
       ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT2AB),ISYTYP(2,LISTT2AB))
       I000=1
       I010=I000+ISIZE*IINTFP
       I020=I010+ISIZE*IINTFP
       CALL GETALL(ICORE(I000),ISIZE,IRREPX,LISTT2AB)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),VRT(1,2),
     &             POP(1,1),POP(1,2),ICORE(I020),IRREPX,'1324')
       CALL PUTALL(ICORE(I010),ISIZE,IRREPX,LISTT2RS-1+4)
       CALL SSTGEN(ICORE(I000),ICORE(I010),ISIZE,VRT(1,1),VRT(1,2),
     &             POP(1,1),POP(1,2),ICORE(I020),IRREPX,'1423')
       CALL PUTALL(ICORE(I010),ISIZE,IRREPX,LISTT2RS-1+6)
C
      ELSE
C
C UHF RESORTS 
C
       LIST(1)=LISTT2RS-1+4
       LIST(2)=LISTT2RS-1+3
       LIST(3)=LISTT2RS-1+6
       LIST(4)=LISTT2RS-1+5
       LIST(5)=LISTT2RS-1+1
       LIST(6)=LISTT2RS-1+2
       SPCASE(1)='AABB'
       SPCASE(2)='BBAA'
       ISIZE=IDSYMSZ(IRREPX,13,14)
       I000=1
       I010=I000+ISIZE*IINTFP
       I020=I010+ISIZE*IINTFP
       DO 10 I=1,2
        CALL GETALL(ICORE(I000),ISIZE,IRREPX,LISTT2AB)
        CALL GSST002(ICORE(I000),ICORE(I010),ISIZE,ISIZE,
     &               ICORE(I020),SPCASE(I),IRREPX)
        CALL PUTALL(ICORE(I010),ISIZE,IRREPX,LIST(I))
        CALL GSSTRNG(ICORE(I010),ICORE(I000),ISIZE,ISIZE,
     &               ICORE(I020),SPCASE(I),IRREPX)
        CALL PUTALL(ICORE(I000),ISIZE,IRREPX,LIST(2+I))
10     CONTINUE
       SPCAS2='AAAA'
       DO 30 ISPIN=1,2
        ISZTOT=IDSYMSZ(IRREPX,ISPIN,ISPIN+2)
        ISZTAR=IDSYMSZ(IRREPX,8+ISPIN,8+ISPIN)
        I000=1
        I010=I000+IINTFP*MAX(ISZTOT,ISZTAR)
        I020=I010+IINTFP*MAX(ISZTOT,ISZTAR)
        CALL GETALL(ICORE(I000),ISZTOT,IRREPX,LISTT2-1+ISPIN)
        CALL GSST003(ICORE(I000),ICORE(I010),ISZTOT,ISZTAR,ICORE(I020),
     &               SPCAS2,'AJBI',IRREPX)
        CALL PUTALL(ICORE(I010),ISZTAR,IRREPX,LIST(4+ISPIN))
        SPCAS2='BBBB'
30     CONTINUE
C
      ENDIF
C
      RETURN
      END