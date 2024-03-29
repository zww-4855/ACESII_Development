










      SUBROUTINE DRVTDA(W,MAXCOR,IUHF)
C
C THIS ROUTINE DRIVES THE TDA CALCULATION.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSIZ
      LOGICAL CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES,PRINT
      LOGICAL DOUBLE,NONSTD, TRIPLET
      DIMENSION W(MAXCOR),KROOT(8),NROOT_ORG(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON/CALCINFO/NROOT(8)
      COMMON/STATSYM/IRREP
      COMMON/GUESS/DOUBLE,NONSTD
      COMMON/GUESS2/IMAP(100,8)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/TDALIST/LISTETDA, LISTVTDA
      COMMON /SPINSTATE/TRIPLET
C
      TOL=10.0D0**(-IFLAGS(98))
      NBAS=NOCCO(1)+NVRTO(1)
      NDIMS = 0
C
      PRINT=IFLAGS(1).GE.0
      CALL IZERO(KROOT,8)
C
      LISTDIAG=42
C
C COMPUTE AND STORE DIPOLE MOMENT INTEGRALS
C
      IONE=1
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NBAST)
      I000=1
      I010=I000+NBAST*NBAST
      CALL DIPLIST(W,W(I010),NBAS,NBAST,IUHF,IRREP)
      CALL ICOPY(8,NROOT,1,NROOT_ORG,1)
C
C WRITE OUT SOME STUFF
C
       IF (PRINT) THEN
       IF(FULDIAG)THEN
        WRITE(6,1010)
1010    FORMAT(T2,' @DRVTDA-I, Roots obtained by full diagonalization.')
       ELSE
        WRITE(6,1020)
1020    FORMAT(T2,' @DRVTDA-I, Roots obtained by Davidson method.')
       ENDIF
       IF (IPATTERN .NE. 0) THEN
         WRITE(6,*) ' Projection is used in TDA calculation'
       ENDIF
       ENDIF 
C    
        IBLOCK = 1
        IF (TRIPLET) IBLOCK = 2
C
C LOOP OVER SYMMETRY BLOCKS, READ THE MATRIX IN AND DIAGONALIZE IT
C
      DO 10 IRREP=1,NIRREP
       IF(NROOT(IRREP).EQ.0)GOTO 10
       CALL NEWLST(IRREP,W,MAXCOR,IUHF)
C This is borrowed from John to debug RPA(D). Commented for standard runs.
CSSS       CALL DDENOMTDA(W,MAXCOR,IUHF,IRREP,447)
C 
       IF(IUHF.EQ.0)THEN
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,18))
        MATDIM = NUMDIS
        I000=1
        I010=I000+NUMDIS*NUMDIS
        CALL GETLST(W(I000),1,NUMDIS,1,IRREP,LISTDIAG)
       ELSE
C
C FORM UHF MATRIX
C
        NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
        NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
        MATDIM=NUMAA+NUMBB
        I000=1
        I010=I000+MATDIM*MATDIM
        IOFFAAAA=1
        IOFFBBBB=IOFFAAAA+MATDIM*NUMAA+NUMAA
        IOFFBBAA=IOFFAAAA+NUMAA
        IOFFAABB=IOFFAAAA+MATDIM*NUMAA
        LISTW=23
        NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
        IOFF=IOFFAAAA
        DO 11 IDIS=1,NUMDIS
         CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
         IOFF=IOFF+MATDIM
11      CONTINUE
        LISTW=24
        NUMDIS=IRPDPD(IRREP,ISYTYP(1,LISTW))        
        IOFF=IOFFBBBB
        DO 12 IDIS=1,NUMDIS
         CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
         IOFF=IOFF+MATDIM
12      CONTINUE
        LISTW=17
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))        
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))        
        IOFF=IOFFBBAA
        IOFF2=IOFFAABB
        DO 13 IDIS=1,NUMDIS
         CALL GETLST(W(IOFF),IDIS,1,1,IRREP,LISTW)
         CALL SCOPY (DISSIZ,W(IOFF),1,W(IOFF2),MATDIM)
         IOFF=IOFF+MATDIM
         IOFF2=IOFF2+1
13      CONTINUE
        NUMDIS=MATDIM
       ENDIF
C       

C  IF IPATTERN .GT. 0 SPECIAL TREATMENT 
C
       IF (IPATTERN .GT. 0) THEN
C
         CALL cALCEXCP(IUHF, W(I010), MAXCOR-I010+1, IRREP,
     &      .FALSE., .FALSE.,NSIZEC, 1)
         CALL FNDNDIMS(IUHF, W(I010), MAXCOR-I010+1, NDIMS,
     &      IRREP)
C
      
      IF (NDIMS .GT. 5000) THEN
         Write(6,"(3a)") " The Dimension of the CIS matrix is greater", 
     &                   " than 5000. Full diagonalization"
         Write(6,"(a)")  " is not used"
      ENDIF 
      IF (NDIMS .LE. 5000) THEN
C
C DIAGONALIZE PROJECTED TDA MATRIX IN CORE
C
        WRITE(6,*) ' Projected TDA matrix is diagonalize in core'
        NUMSOL = MIN(NROOT(IRREP)+3, NDIMS/IBLOCK)

CSSS        Write(*,*) "NDIMS, BLOCK, NROOT(IRREP)+3",NDIMS BLOCK,
CSSS     &               NROOT(IRREP)+3
        IALARGE = 1
        IASMALL = IALARGE + MATDIM * MATDIM
        IEVECL = IASMALL + NDIMS*NDIMS
        IEVECS = IEVECL + MATDIM * NUMSOL
        IEVALS = IEVECS + NDIMS * NDIMS
        IINDEX = IEVALS + NDIMS
        ISCR = IINDEX + NDIMS
        KROOT(IRREP) = NUMSOL
C
        CALL SOLVTDA(W(IALARGE), MATDIM, W(IASMALL),
     &     NDIMS, W(IEVECL), NUMSOL,W(IEVECS), W(IEVALS), W(IINDEX),
     &     TRIPLET, IRREP, IUHF, W(ISCR), MAXCOR-ISCR+1)
C
        NROOT(IRREP) = MIN(NROOT(IRREP), KROOT(IRREP))
        NUMSOL = NROOT(IRREP)
C
        CALL UPDMOI(NUMSOL,MATDIM,IRREP,LISTVTDA,0,0)
        CALL UPDMOI(NUMSOL,1,IRREP,LISTETDA,0,0)
        CALL PUTLST(W(IEVECL),1,NUMSOL,1,IRREP,LISTVTDA)
        CALL PUTLST(W(IEVALS),1,NUMSOL,1,IRREP,LISTETDA)
C
        IF(IFLAGS(91).GE.1)THEN
          CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
          CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,W(IEVALS)+ESCF) 
          CALL CALCDEN(IFLAGS(88),IRREP,W(I010),
     &       (MAXCOR-I010+1)*IINTFP,IUHF)
        ENDIF
      ELSE
        WRITE(6,*) ' IMPLICIT PROJECTION IS USED IN TDA'
        CALL MODTDA(W, NUMDIS, W(I010), MAXCOR-I010+1, IUHF,
     &     NSIZEC, IRREP)
      ENDIF
      ENDIF

C If the IPATTERN (for core projection) is zero, then NDIMS is zero, Test is 
C trivial, Ajith Perera, 09/2020 
C
CSSS      IF (IPATTERN .EQ. 0 .OR. NDIMS .GT. 5000) THEN
C
      IF (IPATTERN .EQ. 0) THEN

      IF (.NOT. (FULDIAG .OR. NROOT(IRREP).GT.15)) THEN
         Write(6,"(3a)") " The Dimension of the CIS matrix is greater",
     &                   " than 5000. Full diagonalization"
         Write(6,"(a)")  " is not used"
      ENDIF
C
C  OBTAIN RELEVANT EIGENVALUES BY DAVIDTDA
C
       IF (FULDIAG .OR. NROOT(IRREP).GT.15) THEN
        NUMSOL=NUMDIS
        ILOCEVC=I000
        ILOCEVL=I010
cYAU 20001127
c old:
c        CALL EIG2(W(ILOCEVL),W,NUMDIS,NUMDIS,0)
c new:
        I020 = I010 + ( NUMDIS * IINTFP )
        SCR_SIZE = MAXCOR + 1 - I020
        CALL DSYEV('V','L',NUMDIS,W(ILOCEVC),NUMDIS,W(ILOCEVL),
     &             W(I020),(SCR_SIZE/IINTFP)-1,INFO)
        IF (INFO.NE.0) THEN
           WRITE(*,*)
     & '@DRVTDA: There was a problem diagonalizing the matrix.'
           CALL ERREX
        END IF
c :end
cYAU
        KROOT(IRREP)=NUMDIS
        IROOT=NUMDIS
        I040=I010+NUMDIS
        I050=I040+NUMDIS*NUMDIS
        ITOP=I040+NUMDIS
       ELSE

      Write(6,"(2a)") " Roots are searched using the Davidson",
     &                " algorithm:"
        MAXOFF1=0
        DO 14 J=1,NROOT(IRREP)
         MAXOFF1=MAX(MAXOFF1,IMAP(J,IRREP))
14      CONTINUE
        maxoff=nroot(irrep)
        IF(MAXOFF.NE.0.AND.MAXOFF1.EQ.MAXOFF.AND..NOT.NONSTD)THEN
         NUMSOL=MIN(MAX(10,MAXOFF+5),NUMDIS/IBLOCK)
        ELSEIF((NONSTD.AND.MAXOFF.NE.0).OR.
     &         (MAXOFF.NE.0.AND.MAXOFF1.NE.MAXOFF))THEN
         NUMSOL=MAXOFF
        ELSE
         NUMSOL=0
        ENDIF
        IF(NUMSOL.EQ.0)THEN
         KROOT(IRREP)=0
         GOTO 10
        ENDIF
        MAXIT=50
        I020=I010+MAXIT*MAXIT
        I030=I020+NUMSOL
        I040=I030+NUMSOL*NUMDIS
        I050=I040+NUMDIS*(MAXIT+NUMSOL)
        I060=I050+NUMDIS
        I070=I060+NUMDIS
        I080=I070+NUMDIS*MAXIT
        I090=I080+NUMDIS
        IF (I090 .GE. MAXCOR) CALL INSMEM("@-DRVTD:",I090,MAXCOR)

        CALL DAVIDTDA(W(I000),W(I010),W(I020),W(I030),W(I040),
     &              W(I050),W(I060),W(I070),W(I080),
     &              MAXIT,NUMSOL,IROOT,NUMDIS,NUMDIS,TOL,IERR,
     &              IMAP(1,IRREP),IRREP, TRIPLET)

        CALL SCOPY(NUMSOL,W(I020),1,W(I000),1)
        ILOCEVC=I030
        ILOCEVL=I000
        KROOT(IRREP)=MIN(IROOT,MAX(NROOT(IRREP),MAXOFF))
        NROOT(IRREP) = MIN(NROOT(IRREP), KROOT(IRREP))
       ENDIF
C
C RESORT BY EIGENVALUE IF NECESSARY
C
       CALL ORDERTDA(IROOT,NUMDIS,W(ILOCEVL),W(ILOCEVC),W(I040),W(I050))

C Reset the IROOT. IROOT=NUMDIS if DAVIDTDA is not used and that
C is correct. It is used to indicate something else in Davidtda. 
C Lets change the old calls so that this can work for both routes.`
CSSS       CALL UPDMOI(IROOT,NUMDIS,IRREP,LISTVTDA,0,0)
CSSS       CALL UPDMOI(IROOT,1,IRREP,LISTETDA,0,0)
CSSS       CALL PUTLST(W(ILOCEVC),1,IROOT,1,IRREP,LISTVTDA)
CSSS       CALL PUTLST(W(ILOCEVL),1,IROOT,1,IRREP,LISTETDA)
C 08/2017, Ajith Perera
C 
       CALL UPDMOI(IROOT,NUMDIS,IRREP,LISTVTDA,0,0)
       CALL UPDMOI(IROOT,1,IRREP,LISTETDA,0,0)
       CALL PUTLST(W(ILOCEVC),1,IROOT,1,IRREP,LISTVTDA)
       CALL PUTLST(W(ILOCEVL),1,IROOT,1,IRREP,LISTETDA)
C
C IF THIS IS A PROPERTY/OPTIMIZATION CALCULATION, THEN COMPUTE
C THE DENSITY AND STORE IT ON DISK.
C
       IF(IFLAGS(91).GE.1)THEN
        CALL GETREC(20,'JOBARC','SCFENEG ',IINTFP,ESCF)
C        WRITE(6,1001)W(ILOCEVL)
C        WRITE(6,1002)W(ILOCEVL)+ESCF
        CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,W(ILOCEVL)+ESCF) 

C Note that Calcden assume a closed shell molecules. 

        CALL CALCDEN(IFLAGS(88),IRREP,W(I010),
     &               (MAXCOR-I010+1)*IINTFP,IUHF)
       ENDIF
C
      ENDIF
C
10    CONTINUE
C
CMN        IF(IUHF.NE.0)THEN
CMN         I000=1
CMN         I010=I000+IINTFP*NBAST*2
CMN         CALL NEW2324(W(I010),W(I000),(MAXCOR-I010)*IINTFP,IUHF,1)
CMN        ENDIF
C
C CALCULATE TDA TRANSITION MOMENTS
C
C DETERMINE MAXIMUM POSSIBLE SINGLE EXCITATION LENGTH
C
      MAXLEN = 0
      DO IRREP = 1, NIRREP
        LEN = IRPDPD(IRREP, 9)
        IF (IUHF .NE. 0) LEN = LEN + IRPDPD(IRREP,10)
        MAXLEN = MAX(MAXLEN,LEN)
      ENDDO
      I000=1
      I010=I000+3*MAXLEN
c     I020 = I000 + MAXLEN
      I020 = I010 + MAXLEN
      CALL TRNMOM1(W(I000),W(I010),W(I020),MAXCOR-I020+1,
     &             NROOT,NROOT_ORG,IUHF)
C
C Lets try to do CIS(D) as debugging mechanism of RPA(D). If this works
C nicely then we will make this option available.A
 
CSSS      IF (CISD) THEN
CSSS          CALL FORM_DBLS_CORRECTNS_2CIS(W,MAXCOR,IUHF)
CSSS      ENDIF 

      RETURN
1001  FORMAT(T3,'Total CIS/TDA excitation energy ',F22.14,' a.u.')
1002  FORMAT(T3,'Total CIS/TDA energy ',F22.14,' a.u.')
      END
