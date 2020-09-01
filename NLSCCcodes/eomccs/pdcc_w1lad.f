










      SUBROUTINE PDCC_W1LAD(ICORE,MAXCOR,IUHF,ISIDE,IRROMEGA)
C
C THIS ROUTINE CALCULATES THE CONTRIBUTION
C
C Z(ab,ij) = [R(ef,ij)*<mn||ef>]*T(ab,mn)] [RIGHT EIGENPROBLEM]
C
C Z(ab,ij) = <mn||ab> * [T(ef,mn)*L(ef,ij)] [LEFT EIGENPROBLEM ]
C
C AND IS USED IN CASES WHERE THE W(abcd) HBAR MATRIX ELEMENTS
C ARE NOT STORED ON DISK [LIST 233 EITHER CONTAINS ONLY THE
C CORRESPONDING MO INTEGRALS OR IS ALTOGETHER ABSENT, SUCH
C AS IS THE CASE IN AO-BASED ALGORITHMS.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR),I0T(2)
      DOUBLE PRECISION ONE,ZILCH,ALPHA
      LOGICAL HBAR_4LCCSD
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS/ IFLAGS2(500)
C

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
      DATA ONE  /1.0D0/
      DATA ZILCH/0.0D0/
C
      LISTZ0=461
      LISTR0=444
      IF(ISIDE.EQ.1)THEN
       LISTW0=14
       LISTT0=44
      ELSEIF(ISIDE.EQ.2)THEN
       LISTW0=44
       LISTT0=14
      ENDIF 

      HBAR_4LCCSD=.FALSE.
      IF (IFLAGS(2) .EQ. 6 .OR. IFLAGS2(117) .EQ. 7) HBAR_4LCCSD=.TRUE.
C
      I0T(1)=1
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      IF(IUHF.EQ.1)THEN
       I0T(2)=I0T(1)+NT(1)*IINTFP
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
      ELSE
       I0T(2)=I0T(1)
      ENDIF
      I000=I0T(2)+NT(2)*IINTFP
C
C Z(ab,ij) = [R(ef,ij)*W(ef,mn)]*TAU(ab,mn)  [ISIDE=1]
C
C Z(ab,ij) = [L(ef,ij)*TAU(ef,mn)]*W(ab,mn)  [ISIDE=2]
C
      DO 10 ISPIN=3,3-2*IUHF,-1
       IF(ISPIN.LE.2)THEN
        I1=ISPIN
        I2=ISPIN
       ELSE
        I1=1
        I2=2
       ENDIF
C
C LOOP OVER KET IRREPS OF *TARGET*.
C
       DO 100 IRREPZR=1,NIRREP
        IRREPZL=DIRPRD(IRREPZR,IRROMEGA)
        IRREPW=IRREPZL
        IRREPT=IRREPW
        LISTW=LISTW0+ISPIN-1
        LISTR=LISTR0+ISPIN-1
        LISTT=LISTT0+ISPIN-1
        LISTZ=LISTZ0+ISPIN-1
        DISSYW=IRPDPD(IRREPW,ISYTYP(1,LISTW))
        NUMDSW=IRPDPD(IRREPW,ISYTYP(2,LISTW)) 
        DISSYT=IRPDPD(IRREPT,ISYTYP(1,LISTW))
        NUMDST=IRPDPD(IRREPT,ISYTYP(2,LISTW)) 
        DISSYZ=IRPDPD(IRREPZL,ISYTYP(1,LISTZ))
        NUMDSZ=IRPDPD(IRREPZR,ISYTYP(2,LISTZ))
        DISSYR=DISSYZ
        NUMDSR=NUMDSZ
        DISSYQ=NUMDSW
        NUMDSQ=NUMDSR
        SIZEQ=NUMDSQ*DISSYQ
        SIZEW=NUMDSW*DISSYW
        SIZET=NUMDST*DISSYT
        SIZEZ=NUMDSZ*DISSYZ
        MAXSIZ=MAX(SIZEW,SIZET,SIZEZ)
        I010=I000+IINTFP*MAXSIZ
        I020=I010+IINTFP*MAXSIZ
        I030=I020+IINTFP*SIZEQ
        CALL GETLST(ICORE(I000),1,NUMDSR,1,IRREPZR,LISTR)
        CALL GETLST(ICORE(I010),1,NUMDSW,1,IRREPW,LISTW)
        
CSSS#ifdef _DCC_FLAG
        If (Ispar) then
        IF (ISIDE .EQ. 2) THEN
           Write(6, "(a,F5.2)") " In PDCC_W1LAD Scaling by",
     &                            Wmnij_scale
           CALL DSCAL(NUMDST*DISSYT, Wmnij_scale, ICORE(I010), 1)
        ENDIF 
        Endif 
CSSS#endif
        IF(ISIDE.EQ.2)THEN
         IF (.NOT. HBAR_4LCCSD) THEN
         CALL FTAU(ICORE(I010),ICORE(I0T(I1)),ICORE(I0T(I2)),
     &             DISSYW,NUMDSW,POP(1,I1),POP(1,I2),VRT(1,I1),
     &             VRT(1,I2),IRREPW,ISPIN,ONE)
         ENDIF 
       ENDIF
C
C                          +
C       Q(mn,ij) = W(ef,mn) * R(ef,ij)
C
        CALL XGEMM ('T','N',DISSYQ,NUMDSQ,DISSYW,ONE,ICORE(I010),
     &              DISSYW,ICORE(I000),DISSYR,ZILCH,ICORE(I020),
     &              DISSYQ)
C
        CALL GETLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
        CALL GETLST(ICORE(I010),1,NUMDST,1,IRREPT ,LISTT)

CSSS#ifdef _DCC_FLAG
        If (Ispar) Then
        IF (ISIDE .EQ. 1) THEN
           Write(6, "(a,F5.2)") " In PDCC_W1LAD Scaling by",
     &                            Wmnij_scale
           CALL DSCAL(NUMDST*DISSYT, Wmnij_scale, ICORE(I010), 1)
        ENDIF 
        Endif 
CSSS#endif
        IF(ISIDE.EQ.1)THEN
         IF (.NOT. HBAR_4LCCSD) THEN
         CALL FTAU(ICORE(I010),ICORE(I0T(I1)),ICORE(I0T(I2)),
     &             DISSYW,NUMDSW,POP(1,I1),POP(1,I2),VRT(1,I1),
     &             VRT(1,I2),IRREPW,ISPIN,ONE)
         ENDIF 
       ENDIF
C
C   Z(ab,ij) = T(ab,mn) * Q(mn,ij)
C
        CALL XGEMM('N','N',DISSYZ,NUMDSZ,DISSYQ,ONE,ICORE(I010),
     &             DISSYT,ICORE(I020),DISSYQ,ONE,ICORE(I000),
     &             DISSYZ)
C
        CALL PUTLST(ICORE(I000),1,NUMDSZ,1,IRREPZR,LISTZ)
C
100    CONTINUE                    
C
10    CONTINUE
C
      RETURN
      END