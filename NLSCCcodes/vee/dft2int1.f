










      SUBROUTINE DFT2INT1(ICORE,MAXCOR,IUHF,ISPIN,IRREPX)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTION OF THE Fme INTERMEDIATE
C  IN THE T1 EQUATION.
C
C         t1(i,a) = t2(im,ae) * F(me) 
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DATA ONE /1.0D0/


      LENA=IRPDPD(IRREPX,9)
      LENB=IRPDPD(IRREPX,10)
      IF(IUHF.EQ.0)THEN
       I0TA=1
       I000=I0TA+IINTFP*LENA
      ELSE
       I0TA=1
       I0TB=I0TA+IINTFP*LENA
       I000=I0TB+IINTFP*LENB
      ENDIF
      MXCOR=(MAXCOR-I000+1)/2
      CALL GETLST(ICORE(I000),1,1,1,3,490)
C
      CALL GT2XF(ICORE(I0TA),IRREPX,1,434,93,0,ICORE(I000),MXCOR,IUHF,0)
C
      CALL GETLST(ICORE(I000),1,1,1,3,490)
      CALL SAXPY (LENA,ONE,ICORE(I000),1,ICORE(I0TA),1)
      CALL PUTLST(ICORE(I0TA) ,1,1,1,3,490)
      IF(IUHF.NE.0)THEN
       CALL GETLST(ICORE(I000),1,1,1,4,490)
       CALL SAXPY (LENB,ONE,ICORE(I000),1,ICORE(I0TB),1)
       CALL PUTLST(ICORE(I0TB) ,1,1,1,4,490)
      ENDIF
C
      RETURN
      END
