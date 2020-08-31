      SUBROUTINE R2L2Y2(ICORE,MAXCOR,IUHF)
C
C DRIVER FOR CALCULATING CONTRIBUTIONS OF R2 AND L2 VECTORS
C TO THE Y1 INTERMEDIATE.
C  
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      DOUBLE PRECISION ONE,ZILCH,X,SDOT
      DIMENSION ICORE(MAXCOR)
      COMMON/STATSYM/IRREPX
      COMMON/REFTYPE/MBPT2,CC, CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
C ADD R(EF,MN) * W(EF,MN) * L(AB,IJ) TERM FOR IRREPX=1 ONLY
C
      I000=1
      IF(.NOT. (MBPT2 .OR.CC2) .AND.IRREPX.EQ.1)THEN
       X=ZILCH 
       DO 10 ISPIN=3,3-2*IUHF,-1
        LISTW=13+ISPIN
        LISTR=460+ISPIN
        DO 20 IRREP=1,NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTR))
         DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTR))
         MAXT=MAX(NUMDIS,DISSIZ)
         I010=I000+IINTFP*NUMDIS*DISSIZ
         I020=I010+IINTFP*NUMDIS*DISSIZ
         CALL GETLST (ICORE(I000),1,NUMDIS,1,IRREP,LISTW)
         CALL GETLST (ICORE(I010),1,NUMDIS,1,IRREP,LISTR)
         IF(IUHF.EQ.0)THEN
          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXT  
          ITMP3=ITMP2+IINTFP*MAXT  
          CALL SPINAD1(IRREP,POP(1,1),DISSIZ,ICORE(I000),
     &                 ICORE(ITMP1),ICORE(ITMP2))
         ENDIF
         X=X+SDOT(NUMDIS*DISSIZ,ICORE(I000),1,ICORE(I010),1)
20      CONTINUE
10     CONTINUE
C
       DO 11 ISPIN=3,3-2*IUHF,-1
        ISIZE=ISYMSZ(ISYTYP(1,60+ISPIN),ISYTYP(2,60+ISPIN))
        I010=I000+IINTFP*ISIZE 
        CALL GETALL(ICORE(I000),ISIZE,1,443+ISPIN)
        CALL GETALL(ICORE(I010),ISIZE,1,60+ISPIN)
        CALL SAXPY (ISIZE,X,ICORE(I000),1,ICORE(I010),1)
        CALL PUTALL(ICORE(I010),ISIZE,1,60+ISPIN)
11     CONTINUE
      ENDIF
C
c CALCULATE CONTRACTIONS
C
       CALL R2L2Y2A(ICORE,MAXCOR,IUHF)
       CALL R2L2Y2B(ICORE,MAXCOR,IUHF)
C
       CALL R2L2Y2C(ICORE,MAXCOR,IUHF)
C
C
C
c SG 2/96 lowercase c's are excess I/O that accomplish nothing
c      IOFF=1
c      DO 112 ISPIN=1,1+IUHF
c       CALL GETLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
c       IOFF=IOFF+NT(ISPIN)*IINTFP
c112    CONTINUE
c      DO 113 ISPIN=3,3-2*IUHF,-1
c       LIST=60+ISPIN
c       LENGTH=ISYMSZ(ISYTYP(1,LIST),ISYTYP(2,LIST))
c       CALL GETALL(ICORE(IOFF),LENGTH,1,LIST)
c       IOFF=IOFF+LENGTH*IINTFP
c113    CONTINUE
c
C
C PLACE G INTERMEDIATES ON DISK
C
      CALL GFORMG(1,IRREPX,14,461,400,ICORE(I000),MAXCOR,0,ONE,ONE,
     &            IUHF)
      CALL GFORMG(IRREPX,IRREPX,444,461,100,ICORE(I000),
     &            MAXCOR,0,ONE,ONE,IUHF)
       CALL R2L2Y2D(ICORE,MAXCOR,IUHF)
c      IOFF=1
c      DO 1112 ISPIN=1,1+IUHF
c       CALL GETLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
c       IOFF=IOFF+NT(ISPIN)*IINTFP
c1112    CONTINUE
c      DO 1113 ISPIN=3,3-2*IUHF,-1
c       LIST=60+ISPIN
c       LENGTH=ISYMSZ(ISYTYP(1,LIST),ISYTYP(2,LIST))
c       CALL GETALL(ICORE(IOFF),LENGTH,1,LIST)
c       IOFF=IOFF+LENGTH*IINTFP
c1113    CONTINUE
c
       CALL R2L2Y2E(ICORE,MAXCOR,IUHF)
c      IOFF=1
c      DO 2112 ISPIN=1,1+IUHF
c       CALL GETLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
c       IOFF=IOFF+NT(ISPIN)*IINTFP
c2112    CONTINUE
c      DO 2113 ISPIN=3,3-2*IUHF,-1
c       LIST=60+ISPIN
c       LENGTH=ISYMSZ(ISYTYP(1,LIST),ISYTYP(2,LIST))
c       CALL GETALL(ICORE(IOFF),LENGTH,1,LIST)
c       IOFF=IOFF+LENGTH*IINTFP
c2113    CONTINUE
c
       CALL R2L2Y2F(ICORE,MAXCOR,IUHF)
c      IOFF=1
c      DO 3112 ISPIN=1,1+IUHF
c       CALL GETLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
c       IOFF=IOFF+NT(ISPIN)*IINTFP
c3112    CONTINUE
c      DO 3113 ISPIN=3,3-2*IUHF,-1
c       LIST=60+ISPIN
c       LENGTH=ISYMSZ(ISYTYP(1,LIST),ISYTYP(2,LIST))
c       CALL GETALL(ICORE(IOFF),LENGTH,1,LIST)
c       IOFF=IOFF+LENGTH*IINTFP
c3113    CONTINUE
c
       CALL R2L2Y2G(ICORE,MAXCOR,IUHF)
c      IOFF=1
c      DO 4112 ISPIN=1,1+IUHF
c       CALL GETLST(ICORE(IOFF),1,1,1,2+ISPIN,90)
c       IOFF=IOFF+NT(ISPIN)*IINTFP
c4112    CONTINUE
c      DO 4113 ISPIN=3,3-2*IUHF,-1
c       LIST=60+ISPIN
c       LENGTH=ISYMSZ(ISYTYP(1,LIST),ISYTYP(2,LIST))
c       CALL GETALL(ICORE(IOFF),LENGTH,1,LIST)
c       IOFF=IOFF+LENGTH*IINTFP
c4113    CONTINUE
C
      RETURN
      END
