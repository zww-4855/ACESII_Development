      SUBROUTINE L1W1(ICORE,MAXCOR,IUHF,IRREPX,LISTL0,LISTT,
     &                LISTW0,LISTZ0)
C
C THIS ROUTINE CARRIES OUT THE CONTRACTION:
C
C  Z(AB,IJ) = P(EF) L(EF,IJ) * T(FM) * <ME||AB>
C
C WHILE T(FM)*<ME||AB> IS FORMALLY A PART OF THE ABEF
C HBAR MATRIX, EVALUATING THIS CONTRIBUTION IN THE
C MANNER ABOVE ALLOWS THE USE OF AO-BASED EVALUATION
C OF THE PARTICLE-PARTICLE LADDER TERM.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH
      DIMENSION ICORE(MAXCOR)
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/   
C
C AAAA AND BBBB SPIN CASES FOR UHF
C
      IF(IUHF.NE.0)THEN
C
       DO 10 ISPIN=1,2
C
        LISTW=LISTW0+ISPIN
        LISTL=LISTL0+ISPIN
        LISTZ=LISTZ0+ISPIN
C
C FIRST FORM INTERMEDIATE Q(I<J,EM) = L(I<J,E<F) * T(FM)
C
        DO 20 IRREPIJ=1,NIRREP
         IRREPEF=DIRPRD(IRREPIJ,IRREPX)
         DISSYL=IRPDPD(IRREPEF,ISYTYP(1,43+ISPIN))
         DISSYLX=IRPDPD(IRREPEF,18+ISPIN)
         NUMDSL=IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
         DISSYQ=IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
         NUMDSQ=IRPDPD(IRREPEF,ISYTYP(1,18+ISPIN))
         I000=1
         I010=I000+NT(ISPIN)*IINTFP
         I020=I010+MAX(DISSYLX*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
         I030=I020+MAX(DISSYL*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
         CALL GETLST(ICORE(I000),1,1,1,ISPIN,LISTT)
         CALL GETLST(ICORE(I020),1,NUMDSL,1,IRREPIJ,LISTL)
         CALL TRANSP(ICORE(I020),ICORE(I010),NUMDSL,DISSYL)
         CALL SYMEXP(IRREPEF,VRT(1,ISPIN),NUMDSL,ICORE(I010)) 
C
C FORM PRODUCT 
C
C      Q(I<J,EM) = L(I<JE,F) * T(FM)
C
         IOFFQ=I020
         IOFFL=I010
         IOFFT=I000
         DO 21 IRREPM=1,NIRREP
          IRREPE=DIRPRD(IRREPM,IRREPEF)
          IRREPF=IRREPM
          NUME=VRT(IRREPE,ISPIN)
          NUMF=VRT(IRREPF,ISPIN)
          NUMM=POP(IRREPM,ISPIN)
          NROW=DISSYQ*NUME
          NCOL=NUMM
          NSUM=NUMF
          CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &               ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFQ),NROW)
          IOFFQ=IOFFQ+NROW*NCOL*IINTFP
          IOFFL=IOFFL+NROW*NSUM*IINTFP
          IOFFT=IOFFT+NUMF*NUMM*IINTFP
21       CONTINUE
C
C NOW FORM CONTRACTION Z(A<B,I<J) = W(A<B,EM) * Q(EM,I<J)
C
         CALL TRANSP(ICORE(I020),ICORE(I000),NUMDSQ,DISSYQ)
         NUMDSZ=IRPDPD(IRREPIJ,ISYTYP(2,43+ISPIN))
         DISSYZ=IRPDPD(IRREPEF,ISYTYP(1,43+ISPIN))
         NUMDSW=IRPDPD(IRREPEF,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPEF,ISYTYP(1,LISTW))
         I010=I000+MAX(NUMDSQ*DISSYQ,NUMDSZ*DISSYZ)*IINTFP 
         I020=I010+NUMDSZ*DISSYZ*IINTFP
         MAXSIZE=MAXCOR-I020+1
         NINCOR =MAXSIZE/(MAX(1,DISSYW)*IINTFP)
         NLEFT=NUMDSW
         NFIRST=1
         NPASS=0
         CALL ZERO(ICORE(I010),NUMDSZ*DISSYZ)
         IOFFQ=I000
1        NREAD=MIN(NLEFT,NINCOR)
         CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREPEF,LISTW)
C
C PERFORM MULTIPLICATION
C
C             Z(A<B,I<J) = W(A<B,EM) * Q(EM,I<J)
C
         CALL XGEMM('N','N',DISSYZ,NUMDSZ,NREAD,ONE,ICORE(I020),
     &              DISSYW,ICORE(IOFFQ),NUMDSQ,ONE,ICORE(I010),DISSYZ)
         IOFFQ=IOFFQ+IINTFP*NREAD
         NFIRST=NFIRST+NREAD
         NLEFT =NLEFT-NREAD
         NPASS=NPASS+1
         IF(NLEFT.NE.0)GOTO 1
C
         CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
         CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(I010),1,ICORE(I000),1)
         CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
C
20      CONTINUE
C
10     CONTINUE
C
      ENDIF
C
C NOW HANDLE ABAB SPIN CASE
C
C
C Z(Ab,Ij) = L(Ef,Ij) * T(fm) * <Ab|Em> + L(Ef,Ij) * T(EM) * <Ab|Mf>
C
      LISTZ=LISTZ0+3
      LISTL=LISTL0+3
C
C DO FIRST TERM - L(Ef,Ij) * T(fm) * <Ab|Em>
C
C
C FIRST FORM INTERMEDIATE Q(Ij,Em) = L(Ij,Ef) * T(fm)
C
      DO 120 IRREPIJ=1,NIRREP
       IRREPEF=DIRPRD(IRREPIJ,IRREPX)
       DISSYL=IRPDPD(IRREPEF,ISYTYP(1,46))
       NUMDSL=IRPDPD(IRREPIJ,ISYTYP(2,46))
       DISSYQ=IRPDPD(IRREPIJ,ISYTYP(2,46))
       NUMDSQ=IRPDPD(IRREPEF,ISYTYP(1,25))
       I000=1
       I010=I000+NT(2)*IINTFP
       I020=I010+MAX(DISSYL*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
       I030=I020+MAX(DISSYL*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
       CALL GETLST(ICORE(I000),1,1,1,1+IUHF,LISTT)
       CALL GETLST(ICORE(I020),1,NUMDSL,1,IRREPIJ,LISTL)
       CALL TRANSP(ICORE(I020),ICORE(I010),NUMDSL,DISSYL)
C
C FORM PRODUCT 
C
C      Q(Ij,Em) = L(Ij,Ef) * T(fm)
C
       IOFFQ=I020
       IOFFL=I010
       IOFFT=I000
       DO 121 IRREPM=1,NIRREP
        IRREPE=DIRPRD(IRREPM,IRREPEF)
        IRREPF=IRREPM
        NUME=VRT(IRREPE,1)
        NUMF=VRT(IRREPF,2)
        NUMM=POP(IRREPM,2)
        NROW=DISSYQ*NUME
        NCOL=NUMM
        NSUM=NUMF
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &             ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFQ),NROW)
        IOFFQ=IOFFQ+NROW*NCOL*IINTFP
        IOFFL=IOFFL+NROW*NSUM*IINTFP
        IOFFT=IOFFT+NUMF*NUMM*IINTFP
121    CONTINUE
C
C NOW FORM CONTRACTION Z(Ab,Ij) = W(Ab,Em) * Q(Em,Ij)
C
       LISTW=LISTW0+4
       CALL TRANSP(ICORE(I020),ICORE(I000),NUMDSQ,DISSYQ)
       NUMDSZ=IRPDPD(IRREPIJ,ISYTYP(2,46))
       DISSYZ=IRPDPD(IRREPEF,ISYTYP(1,46))
       NUMDSW=IRPDPD(IRREPEF,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPEF,ISYTYP(1,LISTW))
       I010=I000+MAX(NUMDSQ*DISSYQ,NUMDSZ*DISSYZ)*IINTFP 
       I020=I010+NUMDSZ*DISSYZ*IINTFP
       MAXSIZE=MAXCOR-I020+1
       NINCOR =MAXSIZE/(MAX(1,DISSYW)*IINTFP)
       NLEFT=NUMDSW
       NFIRST=1
       NPASS=0
       CALL ZERO(ICORE(I010),NUMDSZ*DISSYZ)
       IOFFQ=I000
2      NREAD=MIN(NLEFT,NINCOR)
       CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREPEF,LISTW)
C
C PERFORM MULTIPLICATION
C
C             Z(Ab,Ij) = W(Ab,Em) * Q(Em,Ij)
C
       CALL XGEMM('N','N',DISSYZ,NUMDSZ,NREAD,ONE,ICORE(I020),
     &            DISSYW,ICORE(IOFFQ),NUMDSQ,ONE,ICORE(I010),DISSYZ)
       IOFFQ=IOFFQ+IINTFP*NREAD
       NFIRST=NFIRST+NREAD
       NLEFT =NLEFT-NREAD
       NPASS=NPASS+1
       IF(NLEFT.NE.0)GOTO 2
C
C SYMMETRIZE AND WRITE TO DISK FOR RHF.
C
       IF(IUHF.EQ.0)THEN
        MAXZ=MAX(DISSYZ,NUMDSZ)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXZ
        ITMP3=ITMP2+IINTFP*MAXZ
        CALL SYMRHF3(IRREPEF,IRREPIJ,VRT(1,1),POP(1,1),
     &               DISSYZ,ICORE(I010),ICORE(ITMP1),
     &               ICORE(ITMP2),ICORE(ITMP3))
        CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
       ELSE
        CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
C
C HAVE TO DO TERM 2 FOR UHF CASES
C
C
C Z(Ab,Ij) = L(Ef,Ij) * T(EM) * <Ab|Mf>
C
C
        LISTW=LISTW0+3
C
C FIRST FORM INTERMEDIATE Q(Ij,Mf) = L(Ij,Ef) * T(EM)
C
        IRREPEF=DIRPRD(IRREPIJ,IRREPX)
        DISSYL=IRPDPD(IRREPEF,ISYTYP(1,46))
        NUMDSL=IRPDPD(IRREPIJ,ISYTYP(2,46))
        DISSYQ=IRPDPD(IRREPIJ,ISYTYP(2,46))
        NUMDSQ=IRPDPD(IRREPEF,ISYTYP(2,21))
        MAXL=MAX(NUMDSL,DISSYL)
        MAXQ=MAX(NUMDSQ,DISSYQ)
        I000=1
        I010=I000+NT(1)*IINTFP
        I020=I010+MAX(DISSYL*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
        I030=I020+MAX(DISSYL*NUMDSL,DISSYQ*NUMDSQ)*IINTFP
        CALL GETLST(ICORE(I000),1,1,1,1,LISTT)
        CALL GETLST(ICORE(I020),1,NUMDSL,1,IRREPIJ,LISTL)
        CALL TRANSP(ICORE(I020),ICORE(I010),NUMDSL,DISSYL)
        ITMP1=I020
        ITMP2=ITMP1+IINTFP*MAXL
        ITMP3=ITMP2+IINTFP*MAXL
        CALL SYMTR1(IRREPEF,VRT(1,1),VRT(1,2),NUMDSL,ICORE(I010),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C FORM PRODUCT 
C
C      Q(Ij,fM) = L(Ij,fE) * T(EM)
C
        IOFFQ=I020
        IOFFL=I010
        IOFFT=I000
        DO 122 IRREPM=1,NIRREP
         IRREPE=IRREPM
         IRREPF=DIRPRD(IRREPE,IRREPEF)
         NUME=VRT(IRREPE,1)
         NUMF=VRT(IRREPF,2)
         NUMM=POP(IRREPM,1)
         NROW=DISSYQ*NUMF
         NCOL=NUMM
         NSUM=NUME
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,ICORE(IOFFL),NROW,
     &              ICORE(IOFFT),NSUM,ZILCH,ICORE(IOFFQ),NROW)
         IOFFQ=IOFFQ+NROW*NCOL*IINTFP
         IOFFL=IOFFL+NROW*NSUM*IINTFP
         IOFFT=IOFFT+NUME*NUMM*IINTFP
122     CONTINUE
        ITMP1=I030
        ITMP2=ITMP1+IINTFP*MAXQ
        ITMP3=ITMP2+IINTFP*MAXQ
        CALL SYMTR1(IRREPEF,VRT(1,2),POP(1,1),DISSYQ,ICORE(I020),
     &              ICORE(ITMP1),ICORE(ITMP2),ICORE(ITMP3))
C
C NOW FORM CONTRACTION Z(Ab,Ij) = W(Ab,Mf) * Q(Mf,Ij)
C
        CALL TRANSP(ICORE(I020),ICORE(I000),NUMDSQ,DISSYQ)
        NUMDSZ=IRPDPD(IRREPIJ,ISYTYP(2,46))
        DISSYZ=IRPDPD(IRREPEF,ISYTYP(1,46))
        NUMDSW=IRPDPD(IRREPEF,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPEF,ISYTYP(1,LISTW))
        I010=I000+MAX(NUMDSQ*DISSYQ,NUMDSZ*DISSYZ)*IINTFP 
        I020=I010+NUMDSZ*DISSYZ*IINTFP
        MAXSIZE=MAXCOR-I020+1
        NINCOR =MAXSIZE/(MAX(1,DISSYW)*IINTFP)
        NLEFT=NUMDSW
        NFIRST=1
        NPASS=0
        CALL ZERO(ICORE(I010),NUMDSZ*DISSYZ)
        IOFFQ=I000
3       NREAD=MIN(NLEFT,NINCOR)
        CALL GETLST(ICORE(I020),NFIRST,NREAD,1,IRREPEF,LISTW)
C
C PERFORM MULTIPLICATION
C
C             Z(Ab,Ij) = W(Ab,Mf) * Q(Mf,Ij)
C
        CALL XGEMM('N','N',DISSYZ,NUMDSZ,NREAD,ONE,ICORE(I020),
     &             DISSYW,ICORE(IOFFQ),NUMDSQ,ONE,ICORE(I010),DISSYZ)
        IOFFQ=IOFFQ+IINTFP*NREAD
        NFIRST=NFIRST+NREAD
        NLEFT =NLEFT-NREAD
        NPASS=NPASS+1
        IF(NLEFT.NE.0)GOTO 3
C
        CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
        CALL SAXPY (NUMDSZ*DISSYZ,ONEM,ICORE(I010),1,ICORE(I000),1)
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPIJ,LISTZ)
C
       ENDIF
C
120   CONTINUE
C
      RETURN
      END
