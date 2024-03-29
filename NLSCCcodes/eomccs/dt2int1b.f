










      SUBROUTINE DT2INT1B(ICORE,MAXCOR,IUHF,IRROMEGA,LISTW0,LISTT0,
     &                    LISTT1)
C
C THIS SUBROUTINE CALCULATES THE CONTRIBUTION OF T2 TO
C   T1 (T1<-T2).  THIS CODE IS GENERAL IN THE SENSE THAT
C   IT DOES NOT ASSUME THAT T (AND Z) ARE TOTALLY SYMMETRIC
C   QUANTITIES, BUT RATHER TRANSFORM AS IRROMEGA.
C
C CONTRACTION 1:
C
C     Z(A,I) = SUM  T(AE,MN) * <MN||IE> + T(Ae,Mn) * <Mn|Ie>  [AA]
C              men
C
C     Z(a,i) = SUM  T(ae,mn) * <mn||ie> + T(aE,mN) * <mN|iE>  [BB]
C              men
C
C     Z(A,I) = SUM  [2 T(Ae,Mn) - T(Ea,Mn)] <Mn|Ie> [SPIN ADAPTED RHF] 
C              men
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM
      LOGICAL RHF
C
      DIMENSION ICORE(MAXCOR),IOFFW(8),IOFFZ(8),IOFFT(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE  /1.0D0/
      DATA ONEM/-1.0D0/


      RHF=.FALSE.
      IF(IUHF.EQ.0)RHF=.TRUE.
C
C LOOP OVER SPIN CASES 
C
      DO 50 ISPIN=1,1+IUHF
       LENTAR=IRPDPD(IRROMEGA,8+ISPIN)
       I000=1
       I010=I000+LENTAR*IINTFP
       CALL ZERO(ICORE(I000),LENTAR)
C
C LOOP OVER IRREPS OF RHS OF T VECTOR
C
       DO 100 IRREPTR=1,NIRREP
        IRREPTL=DIRPRD(IRREPTR,IRROMEGA)
        IRREPW=IRREPTR
        LISTW=LISTW0+4-ISPIN
        LISTT=LISTT0+2
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT))
        NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT))
        MAXW=MAX(NUMDSW,DISSYW)
        MAXT=MAX(DISSYT,NUMDST)
        I020=I010+IINTFP*DISSYW*NUMDSW
        I030=I020+IINTFP*MAX(DISSYT*NUMDST,3*MAXW)
        I040=I030+IINTFP*MAX(3*MAXT,DISSYT*NUMDST)
        IF(I040.LT.MAXCOR)THEN
C
C DO IN-CORE ALGORITHM
C
C
C READ W INTO 
C
C        W(Mn,Ie) AND TRANSPOSE KET INDICES TO W(Mne,I)   [ISPIN=1]
C
C        W(Nm,Ei)                                         [ISPIN=2]
C
         CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
         IF(ISPIN.EQ.1)THEN
          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXW
          ITMP3=ITMP2+IINTFP*MAXW
          ITMP4=ITMP3+IINTFP*MAXW
          CALL SYMTR1(IRREPW,POP(1,1),VRT(1,2),DISSYW,ICORE(I010),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
         ENDIF
C
C READ T2 VECTOR INTO 
C
C           T(Ae,Mn) AND REORDER TO T(Mn,eA)           [ISPIN=1]
C
C           T(Ea,Nm) AND REORDER TO T(NM,Ea)           [ISPIN=2]
C
         CALL GETLST(ICORE(I030),1,NUMDST,1,IRREPTR,LISTT)
         CALL TRANSP(ICORE(I030),ICORE(I020),NUMDST,DISSYT)
         IF(ISPIN.EQ.1)THEN
          ITMP1=I030
          ITMP2=ITMP1+IINTFP*MAXT
          ITMP3=ITMP2+IINTFP*MAXT
          ITMP4=ITMP3+IINTFP*MAXT
          CALL SYMTR1(IRREPTL,VRT(1,1),VRT(1,2),NUMDST,ICORE(I020),
     &                ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C SPIN ADAPT FOR RHF
C
          IF(RHF)THEN
           CALL SPINAD3(IRREPTR,POP(1,1),NUMDST,DISSYT,
     &                  ICORE(I020),ICORE(ITMP1),ICORE(ITMP2))
          ENDIF
         ENDIF
C
C COMPUTE OFFSETS FOR W, T AND Z VECTORS
C
         IW=0
         IT=0
         IZ=0
         DO 110 IRR2=1,NIRREP
          IOFFZ(IRR2)=I000+IZ
          IOFFW(IRR2)=I010+IW
          IOFFT(IRR2)=I020+IT
          IRR1Z=DIRPRD(IRR2,IRROMEGA)
          IRR1W=DIRPRD(IRR2,IRREPW)
          IRR1T=DIRPRD(IRR2,IRREPTL)
          IZ=IZ+VRT(IRR1Z,ISPIN)*POP(IRR2,ISPIN)*IINTFP
          IW=IW+DISSYW*VRT(IRR1W,3-ISPIN)*POP(IRR2,ISPIN)*IINTFP
          IT=IT+NUMDST*VRT(IRR1T,3-ISPIN)*VRT(IRR2,ISPIN)*IINTFP
110      CONTINUE
C
C PERFORM MATRIX MULTIPLICATION
C
C                              +          
C         Z(A,I) = SUM T(Mne,A)  * W(Mne,I)             [ISPIN=1]
C                  Efm
C  
C                              +          
C         Z(a,i) = SUM W(Fe,Ma) * T(Fe,Mi)              [ISPIN=2]
C                  FeM
C
         DO 120 IRREPI=1,NIRREP
          IRREPE=DIRPRD(IRREPI,IRREPW)
          IRREPA=DIRPRD(IRREPE,IRREPTL)
          NROW=VRT(IRREPA,ISPIN)
          NCOL=POP(IRREPI,ISPIN)
          NSUM=DISSYW*VRT(IRREPE,3-ISPIN)
          IT=IOFFT(IRREPA)
          IW=IOFFW(IRREPI)
          IZ=IOFFZ(IRREPI) 
          CALL XGEMM('T','N',NROW,NCOL,NSUM,ONEM,ICORE(IT),NSUM,
     &               ICORE(IW),NSUM,ONE,ICORE(IZ),NROW)
120      CONTINUE
        ELSE
C
C OUT-OF-CORE ALGORITHM
C
         WRITE(6,*)' out-of-core AB not coded '
         call errex
C
        ENDIF
C
        IF(.NOT.RHF)THEN
C
C NOW DO THE OTHER SPIN CASE FOR UHF CALCS
C
C     Z(A,I) = Z(A,I) + SUM  T(MN,AE) * <MN||IE>    [ISPIN=1]
C                       MEF
C
C     Z(a,i) = Z(a,i) + SUM  T(mn,ae) * <mn||ie>    [ISPIN=2]
C                       MEF
C
         LISTW=LISTW0-1+ISPIN
         LISTT=LISTT0-1+ISPIN
         DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW)) 
         NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
         DISSYT=IRPDPD(IRREPTL,ISYTYP(1,LISTT)) 
         DISSYT0=IRPDPD(IRREPTL,18+ISPIN)
         NUMDST=IRPDPD(IRREPTR,ISYTYP(2,LISTT)) 
         MAXT=MAX(DISSYT,NUMDST)
         MAXW=MAX(DISSYW,NUMDSW)
         I020=I010+IINTFP*DISSYW*NUMDSW
         I030=I020+IINTFP*MAX(DISSYT0*NUMDST,3*DISSYW)
         I040=I030+IINTFP*MAX(DISSYT0*NUMDST,3*MAXT)
         IF(I040.LT.MAXCOR)THEN
C
C READ W INTO W(M<N,IE) AND TRANSPOSE KET INDICES TO W(M<NE,I)
C
          CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXW
          ITMP3=ITMP2+IINTFP*MAXW
          ITMP4=ITMP3+IINTFP*MAXW
          CALL SYMTR1(IRREPW,POP(1,ISPIN),VRT(1,ISPIN),DISSYW,
     &                ICORE(I010),ICORE(ITMP1),ICORE(ITMP2),
     &                ICORE(ITMP3))
C
C READ T INTO T(A<E,M<N) AND REORDER TO T(M<NE,A)
C
          CALL GETLST(ICORE(I030),1,NUMDST,1,IRREPTR,LISTT)
          CALL TRANSP(ICORE(I030),ICORE(I020),NUMDST,DISSYT)
          CALL SYMEXP(IRREPTL,VRT(1,ISPIN),NUMDST,ICORE(I020))
C
C COMPUTE OFFSETS FOR W, T AND Z VECTORS
C
          IW=0
          IT=0
          IZ=0
          DO 111 IRR2=1,NIRREP
           IOFFZ(IRR2)=I000+IZ
           IOFFW(IRR2)=I010+IW
           IOFFT(IRR2)=I020+IT
           IRR1Z=DIRPRD(IRR2,IRROMEGA)
           IRR1W=DIRPRD(IRR2,IRREPW)
           IRR1T=DIRPRD(IRR2,IRREPTL)
           IZ=IZ+VRT(IRR1Z,ISPIN)*POP(IRR2,ISPIN)*IINTFP
           IW=IW+DISSYW*VRT(IRR1W,ISPIN)*POP(IRR2,ISPIN)*IINTFP
           IT=IT+NUMDST*VRT(IRR1T,ISPIN)*VRT(IRR2,ISPIN)*IINTFP
111       CONTINUE
C
C PERFORM MATRIX MULTIPLICATION
C
C                                     +
C         Z(A,I)   =   SUM   T(M<NE,A) *  W(M<NE,I) 
C                      MEF
C  
          DO 150 IRREPI=1,NIRREP
           IRREPE=DIRPRD(IRREPI,IRREPW)
           IRREPA=DIRPRD(IRREPE,IRREPTL)
           NROW=VRT(IRREPA,ISPIN)
           NCOL=POP(IRREPI,ISPIN)
           NSUM=DISSYW*VRT(IRREPE,ISPIN)
           IZ=IOFFZ(IRREPI)
           IW=IOFFW(IRREPI)
           IT=IOFFT(IRREPA)
           CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,ICORE(IT),NSUM,
     &                ICORE(IW),NSUM,ONE,ICORE(IZ),NROW)
150       CONTINUE
C
         ELSE
C
          WRITE(6,*)' out-of-core not coded yet '
          call errex
C
C ENDS OUT-OF-CORE/INCORE IF
C
         ENDIF
C
C ENDS RHF/UHF IF
C
        ENDIF
C
100    CONTINUE
C
       CALL GETLST(ICORE(I010),1,1,1,2+ISPIN,LISTT1)
       CALL SAXPY (LENTAR,ONE,ICORE(I000),1,ICORE(I010),1)
       CALL PUTLST(ICORE(I010),1,1,1,2+ISPIN,LISTT1)
C
50    CONTINUE
C
      RETURN
      END
