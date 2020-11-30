      SUBROUTINE DOOINVO(IRREPX,IOV,DOO,ICORE,MXCOR,IUHF,LISTW0)
C
C THIS ROUTINE COMPUTES THE CONTRIBUTIONS:
C
C      IOV(I,A) =  [SUM D(M,N) * <MI||NA> + D(m,n) * <Im|An>] (ISPIN=1)
C
C      IOV(i,a) =  [SUM D(m,n) * <mi||na> + D(M,N) * <Mi|Na>] (ISPIN=2)
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION IOV,DOO
      DIMENSION DOO(1),IOV(1),ICORE(MXCOR)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      IF(IUHF.EQ.0)THEN
C
C SPIN ADAPTED RHF CODE
C
       LISTW=LISTW0-1+4
       DO 120 IRREP=1,NIRREP
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        I000=1
        I010=I000+IINTFP*NUMDIS*DISSIZ
        I020=I010+IINTFP*NUMDIS
        I030=I020+IINTFP*NUMDIS
        CALL GETLST(ICORE(I000),1,NUMDIS,1,IRREP,LISTW)
        CALL SPINAD3(IRREP,POP(1,1),DISSIZ,NUMDIS,ICORE(I000),
     &               ICORE(I010),ICORE(I020))
        CALL DDOT24(IRREPX,IRREP,IOV,DOO,ICORE(I000),ICORE(I010),DISSIZ,
     &              POP(1,1),VRT(1,1),POP(1,1),POP(1,1),POP(1,1),
     &              VRT(1,1),'STST')
120    CONTINUE
      ELSEIF(IUHF.NE.0)THEN
       DO 10 ISPIN=1,1+IUHF
        LISTW1=LISTW0-1+ISPIN
        LISTW2=LISTW0-1+2+ISPIN
        DO 20 IRREP=1,NIRREP
         NDSZ1=IRPDPD(IRREP,ISYTYP(1,LISTW1))
         NDSZF=IRPDPD(IRREP,20+ISPIN)
         NDIS1=IRPDPD(IRREP,ISYTYP(2,LISTW1))
         NDSZ2=IRPDPD(IRREP,ISYTYP(1,LISTW2))
         NDIS2=IRPDPD(IRREP,ISYTYP(2,LISTW2))
         IOFFI=1+(ISPIN-1)*IRPDPD(IRREPX,9)
C
C DO FIRST PART - D(M,N) * <MI||NA> [D(m,n) * <mi||na>]
C
         I000=1
         I010=I000+IINTFP*NDSZF*NDIS1
         IOFFD=1+(ISPIN-1)*IRPDPD(IRREPX,21)
         CALL GETLST(ICORE(I000),1,NDIS1,2,IRREP,LISTW1)
         CALL SYMEXP2(IRREP,POP(1,ISPIN),NDSZF,NDSZ1,NDIS1,ICORE(I000),
     &                ICORE(I000))
         CALL DDOT24(IRREPX,IRREP,IOV(IOFFI),DOO(IOFFD),ICORE(I000),
     &              ICORE(I010),NDSZF,POP(1,ISPIN),VRT(1,ISPIN),
     &              POP(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &              'STST')
C
C NOW DO SECOND PART 
C
C          D(m,n) * <Im|An>   [ISPIN=1]
C
C          D(M,N) * <Mi|Na>   [ISPIN=2]
C
         I000=1
         I010=I000+IINTFP*NDSZ2*NDIS2
         IOFFD=1+(2-ISPIN)*IRPDPD(IRREPX,21)
         CALL GETLST(ICORE(I000),1,NDIS2,2,IRREP,LISTW2)
         IF(ISPIN.EQ.1)THEN
          CALL DDOT24(IRREPX,IRREP,IOV(IOFFI),DOO(IOFFD),ICORE(I000),
     &               ICORE(I010),NDSZ2,POP(1,ISPIN),VRT(1,ISPIN),
     &               POP(1,1),POP(1,2),VRT(1,1),POP(1,2),'TSTS')
         ELSE
          CALL DDOT24(IRREPX,IRREP,IOV(IOFFI),DOO(IOFFD),ICORE(I000),
     &               ICORE(I010),NDSZ2,POP(1,ISPIN),VRT(1,ISPIN),
     &               POP(1,1),POP(1,2),POP(1,1),VRT(1,2),'STST')
         ENDIF
20      CONTINUE
10     CONTINUE
      ENDIF
      RETURN
      END
