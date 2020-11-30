      SUBROUTINE GETDT1(ICORE,MAXCOR,MXCOR,IUHF,IOFFT1,IRROMEGA,L)
C
C THIS ROUTINE LOADS THE T1AA AND T1BB VECTORS INTO THE VERY TOP OF
C  CORE, AND RETURNS A POINTER ARRAY GIVING THE OFFSETS WHERE EACH
C  IRREP OF THE TWO SPIN CASES BEGINS.  FOR RHF, THE ADDRESSES OF THE
C  AA AND BB T1 VECTORS ARE IDENTICAL AND ONLY ONE IS HELD.
C
C THIS ROUTINE IS A GENERALIZATION OF GETT1 AND DOES NOT ASSUME
C  THAT T1 IS TOTALLY SYMMETRIC, BUT RATHER TRANSFORMS AS IRROMEGA
C
C  PARAMETERS:
C              ICORE - THE CORE VECTOR (T1 RETURNED AT TOP)
C             MAXCOR - THE TOTAL CORE SIZE 
C              MXCOR - THE SIZE OF CORE BELOW THE START OF THE T1 VECTORS.
C               IUHF - THE UHF/RHF FLAG
C             IOFFT1 - A TWO DIMENSIONAL ARRAY GIVING THE ADDRESS OF
C                       THE BEGINNING OF EACH IRREP IN THE T1 VECTOR.
C                       FOR EXAMPLE, IOFFT1(3,2) GIVES THE ADDRESS OF
C                       THE FIRST ELEMENT OF THE THIRD IRREP FOR SPIN
C                       CASE 2 (BB).
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),IOFFT1(8,2),ITSTART(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB
      MXCOR=MAXCOR
C
C COMPUTE OFFSETS FOR BEGINNING OF T1AA AND T1BB (1000 LOOP) AND OFFSETS FOR
C  BEGINNING OF IRREPS (2000 LOOP).
C
      TLIST=L
      DO 1000 ISPIN=2,2-IUHF,-1
       TLIST2=ISPIN
       IF(IUHF.EQ.0)TLIST2=1
       NLIST=22+ISPIN
       IF(IUHF.EQ.0)THEN
        NLIST=23
        TLIST=L
       ENDIF
       T1SIZ=IRPDPD(IRROMEGA,ISYTYP(1,NLIST))
       ITSTART(ISPIN)=MXCOR-T1SIZ*IINTFP+1
       IF(IUHF.EQ.0)ITSTART(1)=ITSTART(2)
       MXCOR=MXCOR-T1SIZ*IINTFP
       IOFF=ITSTART(ISPIN)
       DO 2000 IRREPI=1,NIRREP
        IRREPA=DIRPRD(IRREPI,IRROMEGA)
        IOFFT1(IRREPI,ISPIN)=IOFF
        IF(IUHF.EQ.0)IOFFT1(IRREPI,1)=IOFF
        IOFF=IOFF+VRT(IRREPA,ISPIN)*POP(IRREPI,ISPIN)*IINTFP
2000   CONTINUE
       CALL GETLST(ICORE(ITSTART(ISPIN)),1,1,1,TLIST2,TLIST)
1000  CONTINUE
      RETURN
      END 
