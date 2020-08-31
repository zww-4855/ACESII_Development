




































































































































































































      SUBROUTINE NEXTDAV(BUF,BUF3,SCR,MAXCOR,CONVRG,
     &                    LISTT2IN,LISTT2,NROOT,IRREPX,ISIDE,
     &                    IUHF,ITER,MAXITER)
C
C DAVIDSON EXTRAPOLATION ROUTINE FOR UNSYMMETRIC MATRICES.  THIS
C APPROACH FOLLOWS THAT OF HIRAO AND NAKATSUJI (J. COMP. PHYS. 45, 246).
C
C THE VECTORS STORED ON DISK ARE AS FOLLOWS:
C
C LIST 470 - OLD C VECTORS, ONE VECTOR PER LOGICAL RECORD
C LIST 471 - OLD HC VECTORS, ONE VECTOR PER LOGICAL RECORD
C LIST 472 - DIAGONAL PART OF H, STORED ON ONE LOGICAL RECORD
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVRG,DOUBLE,ESPROP,LAMTHERE,MBPT2,CC,CCD,RCCD,DRCCD
      LOGICAL NEWVEC,LOCK,NONSTD,EMINFOL,EVECFOL,PRINT,EOM_TRPS 
      LOGICAL LCCD,LCCSD,CC2,GRAD_CALC
      LOGICAL ESTATE_GEOM_OPT,CORE_PROJECT,CORE_SEARCH
      CHARACTER*1 NATURE(100,8)
      INTEGER BGN,END,BGN_IRP,END_IRP
C
      PARAMETER (MAXORD=100)
C
      DIMENSION BUF(*),SCR(MAXCOR)
      DIMENSION BUF3(*),sjunk(100),ijunk(100)
C
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF2/ROOT,OVRLAP,RESID
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/GUESS/DOUBLE,NONSTD
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/PROPGRAD/ESPROP,IOPTROOT,IOPTSYM
      COMMON/RMAT/ R(10000)
      COMMON/FLAGS/IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE  

      COMMON/CNVRGE/EMINFOL,EVECFOL
      COMMON/LAMSTATE/LAMTHERE
CMN
      COMMON/PMAT/ P(10000)
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
      DATA FACT /27.2113957D0/
      SAVE LOCK,Z0, NEWVEC
C
      INDXF(I,J,N)=I+(J-1)*N
C
      CALL TIMER(1)      
      PRINT=IFLAGS(1).GE.10
C
      IF (IPATTERN .LE. 1) THEN
        EXCPTHRS = 0.90D0
      ELSEIF (IPATTERN .GE. 2) THEN
        EXCPTHRS = 0.70D0
      ENDIF
C
C   THRESHOLD TO DECIDE IF VECTOR IS OF INTEREST
C
CMN END
C   
      ESTATE_GEOM_OPT=(Iflags2(5)  .GT. 0)
      CORE_PROJECT   =(Iflags2(120) .GT. 0)
      CORE_SEARCH    =(Iflags2(119) .GT. 0)

      DEGTOL=10.0D0**(-NTOL)
      IF(NDIMR.EQ.1)THEN
       IOLDEST=1
       IF (ISIDE .EQ. 1) LOCK=.FALSE.
      ENDIF
      IONE=1
C
C MOVE MOST RECENT H*C VECTOR TO MOST CURRENT LOCATION.  DO SAME
C FOR C IF THIS IS THE FIRST ITERATION
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &              .FALSE.)
      CALL PUTLST(SCR,IOLDEST,1,1,ISIDE,470)
      Write(6,"(a)")" The checksum of the 490,443-445"
      Call checksum(" @-NEXTDAV       :",SCR(I000),NSIZEC,S1)
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,2,460,NSIZEC,
     &              .FALSE.)
      Write(6,"(a)")" The checksum of the 490,460-463"
      Call checksum(" @-NEXTDAV       :",SCR(I000),NSIZEC,S1)
C
      IF (IPROJECT .GE. 2) THEN
C
C PROJECT OUT UNINTERESTING COMPONENTS.
C
        I000 = 1
        I010 = 1 + NSIZEC
        CALL GETLST(SCR(I010),ICOLPR2,1,1,1,LISTH0)  

C If the core projection is turned on with %EXCITE, lets start the
C projection after first iteration. Otherwise, it is quite 
C possible to end up with a zero eigenvalue at first iteration and
C quit. This can also be controlled with LOCK and the answers are
C different. I do not know which one is better. Ajith Perera, 12/2019.
 
       IF (NONSTD .AND. CORE_PROJECT) THEN
          IF (NDIM .GT. 1) CALL VECPRD(SCR(I000),SCR(I010),SCR(I000),
     &                                 NSIZEC)
        ELSE 
          CALL VECPRD(SCR(I000),SCR(I010),SCR(I000),NSIZEC) 
        ENDIF 

      Write(6,"(a)")" The checksum of the PR1+PR2"
      Call checksum(" @-NEXTDAV       :",SCR(I000),NSIZEC,S1)
      Write(6,*)
      ENDIF
C
      CALL PUTLST(SCR,IOLDEST,1,1,ISIDE,471)
C
      IOLDEST=1+MOD(IOLDEST,MAXORD+1)

      Write(6,"(a,i3)") " Oldest:", ioldest 
C
C AUGMENT EXISTING SUBSPACE MATRIX
C
      CALL AUGMENTR(SCR, MAXCOR, ISIDE, IUHF, IRREPX)
CMN
      IF (IPROJECT .GE. 1) THEN
         CALL AUGMENTP(SCR, MAXCOR, ISIDE, IUHF, IRREPX)
      ENDIF
C
C DIAGONALIZE SUBSPACE MATRIX
C
      CALL DIAGR(BUF3, BUF, MAXORD, IEVAL, IEVEC, IEXCP,
     $   IEVALSEL,EXCPTHRS, PRINT)
C
      IBUFTMP = IEVALSEL + NDIMR*NDIMR
C
C FOLLOW EIGENVECTOR ASSOCIATED WITH LOWEST UNCONVERGED EIGENVALUE!
C
      IF (IPROJECT .GE. 1) THEN
        IEVALACT = IEVALSEL
      ELSE
        IEVALACT = IEVAL
      ENDIF
C
C  IF THIS IS A NEW VECTOR, UPDATE SJUNK
C
      IF(NDIMR.NE.1) THEN
        IF (NEWVEC) THEN
          I000 = 1
          I010 = I000 + NSIZEC
          IREAD = IROOT + 1
          IPRINT = IFLAGS(1)
          CALL STARTV(SJUNK, SCR(I000), NSIZEC, NDIMR, IOLDEST,
     &       MAXORD, SCR(I010), IRREPX, IUHF, ISIDE, IPRINT, IREAD,
     &       IJUNK)
        ELSE
          CALL GETREC(20,'JOBARC','LASTVECT',100*IINTFP,SJUNK)
        ENDIF
      ENDIF
CMN END
      IF(.NOT.LOCK)THEN
       IF(EMINFOL)THEN
         CALL FNDMINE(NDIMR,BUF(IEVALACT),EIGVAL(1,IRREPX),BUF(IBUFTMP),
     &      IROOT,Z,ILOC,1.D-3) 
       ELSEIF(EVECFOL.AND.NDIMR.NE.1)THEN
         CALL FNDMAXS(NDIMR,SJUNK,BUF(IEVEC),NDIMR,
     &      BUF(IBUFTMP),OVRLAP,ILOC,DUMMY)
       ENDIF
C
C LOCK ON ROOT IF RELATIVE CHANGE IN EIGENVALUE IS LESS THAN TOLERANCE
C
       IF(ROOT.NE.ZILCH)THEN
        XTEST=ABS((BUF(IEVAL -1 + ILOC)-ROOT)/ROOT)
       ELSE
        XTEST=ABS(BUF(IEVAL - 1 + ILOC)-ROOT)
       ENDIF
       IF(ABS(XTEST).LT.1.D-3.AND.NDIMR.NE.1)THEN
        LOCK=.TRUE.
        WRITE(6,*)' L-O-C-K-I-N-G O-N R-O-O-T '
       ENDIF
      ELSE
        IF (ISIDE .EQ. 2) ROOT = EIGVAL(IROOT+1,IRREPX)
         CALL FNDNEAR(NDIMR,BUF(IEVALACT),EIGVAL(1,IRREPX),
     &     ROOT,BUF(IBUFTMP),IROOT,XCLOSE,ILOC,1.D-3) 
C
C       CALL FNDCLOSE(NDIMR,BUF(IEVALACT),ROOT,XCLOSE,ILOC)
c        CALL FNDMINE(NDIMR,BUF,EIGVAL(1,IRREPX),BUF(IBUFTMP),IROOT,
c     &               Z,ILOC,1.D-2) 
      ENDIF
C
      IF(NDIMR.EQ.1)THEN
       CALL ZERO(SJUNK,NDIMR+1)
       ILOC=1
       CALL SCOPY(NDIMR,BUF(IEVEC+(ILOC-1)*NDIMR),1,SJUNK(2),1)
       CALL PUTREC(20,'JOBARC','LASTVECT',100*IINTFP,SJUNK)
      ENDIF
C
      IBUFLC=IEVEC-1+INDXF(1,ILOC,NDIMR)
      OVRLAP=SDOT(NDIMR,SJUNK,1,BUF(IBUFLC),1)
      ROOT=BUF(IEVAL+ILOC-1)

      Write(6,"(a,F15.8)") "Root = ", Root
C
      CALL ZERO(SJUNK,NDIMR+1)
      CALL SCOPY(NDIMR,BUF(IBUFLC),1,SJUNK(2),1)
      CALL PUTREC(20,'JOBARC','LASTVECT',100*IINTFP,SJUNK)
C
C CALCULATE RESIDUAL : [H C(new) - E C(new)
C
      I000=1
      I010=I000+NSIZEC
      I020=I010+NSIZEC
      ITOP=I020+NSIZEC
      CALL FORMS(NSIZEC,NDIMR,SCR(I000),SCR(I020),BUF(IBUFLC),
     &           ISIDE,471,IOLDEST,MAXORD)
      Write(6,*) 
      Write(6,"(a)")" The checksum of 471:HC"
      Call checksum(" @-NEXTDAV       :",SCR(I000),NSIZEC,S1)
      CALL FORMS(NSIZEC,NDIMR,SCR(I010),SCR(I020),BUF(IBUFLC),
     &           ISIDE,470,IOLDEST,MAXORD)

      Write(6,"(a)")" The checksum of 470:C"
      Call checksum(" @-NEXTDAV       :",SCR(I010),NSIZEC,S1)
      CALL SAXPY(NSIZEC,ROOT*ONEM,SCR(I010),1,SCR(I000),1)

      Write(6,"(a)")" The checksum of the residual"
      Call checksum(" @-NEXTDAV",SCR(I000),NSIZEC,S1)
C
      IF (IPROJECT .GE. 2) THEN
C
C PROJECT OUT UNINTERESTING COMPONENTS.
C
        CALL GETLST(SCR(I010),ICOLPR2,1,1,1,LISTH0)
        CALL VECPRD(SCR(I000),SCR(I010),SCR(I000),NSIZEC)
      ENDIF
C
      CALL FNDMAXD(NSIZEC, SCR(I000), DIFF, I)
C
      CALL GETLST(SCR(I010),1,1,1,1,472)
      Write(6,"(a)")" The checksum of the denominator"
      Call checksum(" @-NEXTDAV       :",SCR(I010),NSIZEC,S1)
C
C SAVE THIS VECTOR SO THAT WE CAN USE IT IF THE DAVIDSON
C CORRECTION GENERATES A NEAR LINEAR DEPENDENCE
C
      CALL PUTLST(SCR(I000),ICOLPR2+1,1,1,1,472)
C
C NOW APPLY THE DAVIDSON CORRECTION
C
      CALL VECDIV2(ROOT,SCR(I000),SCR(I010),SCR(I000),NSIZEC)
      IF(IUHF.EQ.0)THEN
       CALL SYMT2AB(IRREPX,SCR,SCR(I010),MAXCOR-I010+1)
      ENDIF
      Write(6,"(a)")" The checksum of the correction vector"
      Call checksum(" @-NEXTDAV       :",SCR(I000),NSIZEC,S1)
      Write(6,*)
C
C SCHMIDT ORTHOGONALIZE NEW BASIS VECTOR TO EXISTING SPACE.
C
      CALL GSORTHOG(IRREPX,NSIZEC,NDIMR,SCR(I000),SCR(I010),
     &              SCR(I020),ISIDE,470,IUHF.EQ.0,IOLDEST,MAXORD,
     &              BUF(IBUFTMP),RESID) 
      WRITE(6,5000)ITER,NDIMR,ROOT,ROOT*FACT,ABS(OVRLAP),DIFF
C
      IF(ABS(ROOT).LT.1.D-5)THEN
       WRITE(6,*)' SUBSPACE EXHAUSTED! '
       CONVRG=.TRUE.
       IROOT=NROOT
       RETURN
      ENDIF
C
      IF(DIFF.GT.DEGTOL .AND. ITER .LE. MAXITER)THEN
       CONVRG=.FALSE.
       NEWVEC = .FALSE.
      ELSE
        IF (DIFF .GT. 1.E-2) THEN
          CONVRG = .FALSE.
          NEWVEC = .FALSE.
          WRITE(6,*)
          WRITE(6,"(a,a,a)")" Maximum allowed iterations is reached",
     &                    " and the convergence criteria is still"
          WRITE(6,"(a)")  " above 0.01."
          WRITE(6,"(a,a)")  " Completely giving up on this symmetry",
     &                      " block!!"
          WRITE(6,*)
          IROOT = NROOT
        ELSE
          CONVRG=.TRUE.
        ENDIF
      ENDIF
      IF (CONVRG) THEN
       IF (DIFF .GT. DEGTOL) THEN
         WRITE(6,"(a,a,ES8.2E2)")" Convergence criteria is below"
     &                          " 0.01 but not below the convergence",
     &                          " tolerance ", DEGTOL
         WRITE(6,"(a,a)")" Not strictly converged!! The excitation",
     &                    " energy is assigned zero (not converged)."
         WRITE(6,"(a)")  " Moving to the next root in this symmetry",
     &                   " block.!!"
         WRITE(6,*)
         IROOT=IROOT+1
         EIGVAL(IROOT,IRREPX)=ZILCH 
         OSCSTR(IROOT,IRREPX)=ZILCH
       ELSE
         IROOT=IROOT+1
         WRITE(6,*)' Converged eigenvalue: ',root,' a.u.'
         WRITE(6,*)'                       ',root*fact,' eV'
         EIGVAL(IROOT,IRREPX)=ROOT
         OSCSTR(IROOT,IRREPX)=ZILCH
       ENDIF
       CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
C 
C Here IFLAGS(87)=7 is P-EOM. 
C
       IF (IFLAGS(87) .EQ. 7) THEN
         IF(CC)THEN
           WRITE(6,98)ECC+ROOT
         ELSE
           WRITE(6,198)ECC+ROOT
         ENDIF
       ELSE
         IF(CC)THEN
           WRITE(6,99)ECC+ROOT
         ELSE IF (MBPT2) THEN
           WRITE(6,199)ECC+ROOT
         ELSE IF (CCD) THEN
           WRITE(6,200)ECC+ROOT
         ELSE IF (RCCD) THEN
           WRITE(6,201)ECC+ROOT
         ELSE IF (DRCCD) THEN
           WRITE(6,202)ECC+ROOT
         ELSE IF (LCCD) THEN
           WRITE(6,203)ECC+ROOT
         ELSE IF (LCCSD) THEN
           WRITE(6,204)ECC+ROOT
         ELSE IF (CC2) THEN
           WRITE(6,205)ECC+ROOT
         ENDIF
       ENDIF
 98    FORMAT(T3,' Total P-EOM-CCSD electronic energy ',
     &    F20.12,' a.u.')
 99    FORMAT(T3,' Total EOM-CCSD electronic energy ',
     &    F20.12,' a.u.')
 198   FORMAT(T3,' Total P-EOM-MBPT(2) electronic energy ',
     &    F20.12,' a.u.')
 199   FORMAT(T3,' Total EOM-MBPT(2) electronic energy ',
     &    F20.12,' a.u.')
 200   FORMAT(T3,' Total CCD electronic energy ', F20.12,' a.u.')
 201   FORMAT(T3,' Total rCCD electronic energy ', F20.12,' a.u.')
 202   FORMAT(T3,' Total drCCD electronic energy ', F20.12,' a.u.')
 203   FORMAT(T3,' Total LCCD electronic energy ', F20.12,' a.u.')
 204   FORMAT(T3,' Total LCCSD electronic energy ', F20.12,' a.u.')
 205   FORMAT(T3,' Total CC2 electronic energy ', F20.12,' a.u.')

       CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,ECC+ROOT)
C
C PUT CONVERGED VECTOR ON LIST 472
C
       CALL FORMS(NSIZEC,NDIMR,SCR(I000),SCR(I010),BUF(IBUFLC),
     &            ISIDE,470,IOLDEST,MAXORD)
       CALL SCOPY(NSIZEC,SCR(I000),1,SCR(I010),1)
       IF(IUHF.EQ.0)THEN
        CALL SPNTSING(IRREPX,SCR(I010),SCR(I020),MAXCOR-I010+1)
       ENDIF
       Z=SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
       X=ONE/SQRT(Z)
       CALL SSCAL (NSIZEC,X,SCR(I000),1)
C
C PRINT MAIN COMPONENTS OF EIGENVECTOR
C
       CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
       CALL PRMAINX(IUHF, SCR(I010), MAXCOR, IRREPX,IROOT)
C
C Save the R1 and R2 to a unformatted file. They may be used as a
C starting guess for certain CC calculations. In orde to do this
C the EOM calculation must be root specific (and asked for a only
C one root). Otherwise, RGUESS file corrspond to the last root 
C of the last irrep. 09/2015, Ajith Perera
     
       IF (ISIDE .EQ. 1) CALL DDMPRGSS(SCR(ITOP),MAXCOR-ITOP,IUHF,
     &                                 400,'RGUESS  ')
C
C Write a file in the %EXCITE format. This can be used as a starting
C guess of another EOM calculation to obtain a specific root.
C Lets not do this when EOM gradients are computed.

      GRAD_CALC = (IFLAGS2(h_IFLAGS2_grad_cal) .GT. 0 .OR.
     &             IFLAGS2(5) .GT. 0)

       IF (ISIDE .EQ. 1 .AND. .NOT. GRAD_CALC) 

     &     CALL FORM_PERCENT_EXCITE(SCR(ITOP),
     &                              MAXCOR-ITOP,IRREPX,IUHF,ROOT,
     &                              ECC,IROOT)
C
       IF(ESPROP)THEN
         CALL TIMER(1)      
         TDAVID = TDAVID + TIMENEW
        IF(ISIDE.EQ.1)THEN
         CALL RNORM(IRREPX,NSIZEC,ROOT,SCR(I000),SCR(I010),
     &              MAXCOR-I010+1,IUHF,Z0)
         CALL PUTLST(SCR(I000),1+ISIDE,1,1,1,472)
         IF(LAMTHERE)THEN
          CALL TDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &               Z*(ONE+DFLOAT(1-IUHF)),Z0,ROOT,CORE_SEARCH)
         ENDIF
         CALL PUTREC(20,'JOBARC','RIGHTDIM',IONE,NDIMR)
         CALL PUTREC(20,'JOBARC','RIGHTOLD',IONE,IOLDEST)
         CALL PUTREC(20,'JOBARC','RIGHTVEC',MAXEXP*MAXEXP*IINTFP,
     &               BUF(IEVEC))
        ELSEIF(ISIDE.EQ.2)THEN

C Use L=T^t approximation 
CSSS
CSSS         CALL GETLST(SCR(I000), 2, 1, 1, 1, 472)
CSSS
         CALL PUTLST(SCR(I000),1+ISIDE,1,1,1,472)
         CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
         CALL LNORM(IRREPX,ROOT,SCR,MAXCOR,IUHF,Z,Z0)

         IF(LAMTHERE)THEN
          CALL TDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &               Z*(ONE+DFLOAT(1-IUHF)),Z0,ROOT,CORE_SEARCH)
          CALL PRINTSUM(ROOT,FSTR,CORE_SEARCH)
          OSCSTR(IROOT,IRREPX) = FSTR
         ENDIF
         CALL GETLST(SCR(I000),3,1,1,1,472)
         CALL SSCAL (NSIZEC,Z,SCR(I000),1)
         CALL PUTLST(SCR(I000),3,1,1,1,472)
         ITROOT=ITROOT+1

         IF (IFLAGS(91) .LT. 2) THEN

           CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &        Z*(ONE+DFLOAT(1-IUHF)),Z0)

         ELSE 

           IF (ESTATE_GEOM_OPT .AND. IROOT .EQ. NROOT) THEN

               CALL INIGAM(IUHF)
               CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                    Z*(ONE+DFLOAT(1-IUHF)),Z0)
C
C CALCULATE EOM-CCSD TWO-PARTICLE DENSITY
C
               CALL TPDENS(SCR,IINTFP*MAXCOR,IUHF,Z0)
               CALL ACES_FIN
               STOP
           ELSE 

               IF (IROOT .EQ. NROOT) THEN
                   CALL INIGAM(IUHF)
                   CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                        Z*(ONE+DFLOAT(1-IUHF)),Z0)

C CALCULATE EOM-CCSD TWO-PARTICLE DENSITY

                   CALL TPDENS(SCR,IINTFP*MAXCOR,IUHF,Z0)
                   CALL ACES_FIN
                   STOP

               ELSE 
C
C This block is added to compute the S^2 for many states in a 
C single run
                   CALL SAVE_HBAR(SCR,MAXCOR,IUHF)
                   CALL INIGAM(IUHF)
                   CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                        Z*(ONE+DFLOAT(1-IUHF)),Z0)
C
C CALCULATE EOM-CCSD TWO-PARTICLE DENSITY
C
                   CALL TPDENS(SCR,IINTFP*MAXCOR,IUHF,Z0)
                   CALL COPY_HBAR(SCR,MAXCOR,IUHF)
                   IF (IROOT .EQ. NROOT) THEN
                      CALL ACES_IO_REMOVE(54,"DERGAM")  
                      CALL MAKLST(SCR,MAXCOR,IUHF)
                      RETURN
                   ENDIF 

               ENDIF 
           ENDIF
         ENDIF
        ENDIF
         CALL TIMER(1)      
       ENDIF
C
C The EOM triples are added 02/2014. Most of the routines that are
C used here originated from John D. Watts contributions. Ajith Perera.

       EOM_TRPS = (iflags(87) .EQ. 11)
       If (EOM_TRPS .AND. ISIDE .EQ. 2)  
     &     CALL ADD_TRPS(SCR, MAXCOR, IRREPX, ROOT, IUHF, ISIDE,
     &     IROOT) 
C
C ADD NEW VECTOR TO SPACE IF WE NEED TO GET MORE ROOTS
C
       NEWVEC=(IROOT.LT.NROOT).AND.
     &        (.NOT.ESPROP.OR.ESPROP.AND.ISIDE.EQ.2)
 
       IF(NEWVEC)THEN
        WRITE(6,*)' Adding new vector to expansion space.'
        LOCK=.FALSE.
        IF(ESPROP)THEN
C
C RESTORE EXISTING RIGHT VECTOR EXPANSION SPACE FOR NEXT ROOT
C
         CALL GETREC(20,'JOBARC','RIGHTDIM',IONE,NDIMR)
         CALL GETREC(20,'JOBARC','RIGHTOLD',IONE,IOLDEST)
         CALL GETREC(20,'JOBARC','RIGHTVEC',NDIMR*NDIMR*IINTFP,
     &               BUF(IEVEC))
         CALL REFORM(R,SCR,NDIMR,NSIZEC,IRREPX,MAXCOR,1,IUHF)
         IF (IPROJECT .GE. 1) THEN
           CALL REFORMP(P,SCR,NDIMR,NSIZEC,IRREPX,MAXCOR,1,IUHF)
         ENDIF
         IF (NDIMR .EQ. MAXEXP) THEN
           CALL DIAGR(BUF3, BUF, MAXORD,IEVAL,IEVEC,IEXCP,
     $        IEVALSEL,EXCPTHRS,PRINT)
         ENDIF
        ENDIF
C
        IREAD=IROOT+1

        CALL GETGES(SCR, NSIZEC, IREAD, IRREPX, ROOT,IJUNK)
C
CMN        IF(IUHF.EQ.0)THEN
CMN       CALL SYMT2AB(IRREPX,SCR,SCR(I010),MAXCOR-I010+1)
CMN        ENDIF
C
        IF (EVECFOL) THEN
           NDIMR  = 0
           IOLDEST= 1
        ENDIF 

        IF(.NOT.NONSTD)WRITE(6,*)' Guess for next eigenvalue ',root
        CALL GSORTHOG(IRREPX,NSIZEC,NDIMR,SCR(I000),SCR(I010),
     &                SCR(I020),1,470,IUHF.EQ.0,IOLDEST,MAXORD,
     &                SJUNK(2),RESID) 
        SJUNK(1)=ONEM
        Z=SNRM2(NDIMR,SJUNK,1)
        CALL SSCAL(NDIMR,ONE/Z,SJUNK,1)
        CALL PUTREC(20,'JOBARC','LASTVECT',100*IINTFP,SJUNK)
        IF(ISIDE.EQ.2)ISIDE=1
C     
       ELSEIF(ESPROP.AND.ISIDE.EQ.1)THEN
        CALL ZERO(R,MAXORD*MAXORD)
        IROOT=IROOT-1
        ISIDE=2
        NDIMR=0
        CALL GETLST(SCR,2,1,1,1,472)
       ENDIF
      ENDIF
C
C UPDATE T VECTOR ON LISTS
C
      IF(IROOT.NE.NROOT)THEN 
       CALL UPDATES(IRREPX,SCR,444,0,490,IUHF)
      ENDIF
C
C TRUNCATE EXPANSION SPACE IF WE HAVE REACHED THAT POINT
C
      IF(NDIMR.GE.(MAXEXP-1) .AND. .NOT. NEWVEC)THEN
       CALL ZERO(R,MAXORD*MAXORD)
       I000 = 1
       I010 = I000 + NREDUCE
       I020 = I010 + NDIMR * NREDUCE
       IF (IPATTERN .NE. 1) THEN
         IE = IEVAL
       ELSE
         IE = IEVALSEL
       ENDIF
       NDIMNEW = NREDUCE
C
C  IN PRINCIPLE NDIMNEW WILL EQUAL NREDUCE, HOWEVER IF TOO FEW
C  REAL, INTERESTING EIGENVALUES ARE AVAILABLE THEN NDIMNEW WILL BE
C  EQUAL TO THE NUMBER OF SUCH EIGENVALUES.
C
       CALL TRUNCATE(BUF(IE),BUF(IEVEC),SCR(I000), SCR(I010),
     &    SCR(I020),MAXCOR-I020+1,NSIZEC,
     &   NDIMR,NDIMNEW,IRREPX,ISIDE,IUHF,ICOLPR2)
       IOLDEST=NDIMNEW + 1
       NDIMR=NDIMNEW
       CALL REFORM(R,SCR,NDIMR,NSIZEC,IRREPX,MAXCOR,ISIDE,IUHF)
       IF (IPROJECT .GE. 1) THEN
         CALL REFORMP(P,SCR,NDIMR,NSIZEC,IRREPX,MAXCOR,1,IUHF)
       ENDIF
      ENDIF
C
C MOVE R AND P MATRICES
C
       DO 50 I=NDIMR,1,-1
        DO 51 J=NDIMR,1,-1
         R(INDXF(I+1,J+1,MAXORD))=R(INDXF(I,J,MAXORD))
         P(INDXF(I+1,J+1,MAXORD))=P(INDXF(I,J,MAXORD))
51      CONTINUE
50     CONTINUE
C     
      CALL TIMER(1)      
      TDAVID = TDAVID + TIMENEW
      RETURN
5000  FORMAT(T4,I4,T13,I4,T21,D14.7,T36,D14.7,T51,D11.6,T63,D8.3)
      END

