      SUBROUTINE R2L2Y1E(Y1,ICORE,MAXCOR,IUHF)
C
C Y1(ai) = - 1/2 R(mn,ef)*L(mo,ef)*W(ni,oa)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),Y1(*)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      LENGTH=IRPDPD(1,21)+IUHF*IRPDPD(1,22)
C
      I0G=1
      I000=I0G+IINTFP*LENGTH
      MXCOR = MAXCOR - I000 + 1
C
C FIRST PICK UP G(no) FROM LIST 191
C
      CALL GETLST(ICORE(I0G),1,1,1,1,191)
      IF(IUHF.NE.0)THEN
       CALL GETLST(ICORE(I0G+IINTFP*IRPDPD(1,21)),1,1,1,2,191)
      ENDIF
C
C FORM PRODUCT Y(ai) = - W(Ni,Oa) * G(no)
C
      CALL VMINUS (ICORE(I0G),LENGTH)
      IOFFY=1
      DO 10 ISPIN=1,1+IUHF
       CALL SYMTRA (1,VRT(1,ISPIN),POP(1,ISPIN),1,Y1(IOFFY),ICORE(I000))
       CALL SCOPY  (NT(ISPIN),ICORE(I000),1,Y1(IOFFY),1)
       IOFFY=IOFFY+NT(ISPIN)*IINTFP
10    CONTINUE
      CALL DOOINVO(1,Y1,ICORE(I0G),ICORE(I000),MXCOR,IUHF,7)
      IOFFY=1
      DO 11 ISPIN=1,1+IUHF
       CALL SYMTRA (1,POP(1,ISPIN),VRT(1,ISPIN),1,Y1(IOFFY),ICORE(I000))
       CALL SCOPY  (NT(ISPIN),ICORE(I000),1,Y1(IOFFY),1)
       IOFFY=IOFFY+NT(ISPIN)*IINTFP
11    CONTINUE
C
      RETURN
      END