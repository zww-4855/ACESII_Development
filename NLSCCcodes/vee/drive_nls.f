




































































































































































































        subroutine DRIVE_NLS(SCR,MAXCOR,IRREPX,IUHF,IROOT)
C
C The core EEs are obtained by setting the R(ij,ab) = 0 when i and
C j do  not belong to the core region. The present version works
C state specifc fashion.
C
C                           all; each
C LIST 444:    C(IJ,AB )       A<B ; I<J     AA AA
C      445:    C(ij,ab )       a<b ; i<j     BB BB
C      446:    C(Ij,Ab )       A,b ; I,j     AB A
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR(MAXCOR)
      DOUBLE PRECISION VALUE,O_SCALE,V_SCALE
      DOUBLE PRECISION VAL1,VAL2,DIFF,THRES2
      INTEGER INTOO(2),INTOV(2)
C
c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN
      PARAMETER (MAXBASFN=1000)
c maxbasfn.par : end
C
C
      CHARACTER*12 STRINGI(2), STRINGJ(2), STRINGA(2),STRINGB(2)
      DIMENSION LS2OUT(2,2)
      INTEGER LS1OUT
      LOGICAL NLS_EXIST,NBO_EXIST,PROJECT_SINGLES
      integer ierr,natm,nbas
      INTEGER QM1num,NLMOnum,NLMOnum2,NLMOorigin,NLMOorigin2,indx
      INTEGER,allocatable::QM1atoms(:),QM2atoms(:),NLMO(:,:)
      integer,allocatable:: NLMOQM1(:),NLMOQM2(:)
      integer::compareIA(2),QMregIA(2)


      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /PROJECT/ IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
C
      DATA BLANK /" "/

c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
C
      PROJECT_SINGLES = (Iflags2(173) .NE. 0)
      print*, 'value of proj singles: ', PROJECT_SINGLES
      call getrec(1,'JOBARC','NBASTOT',1,nbas)
      call getrec(1,'JOBARC','NREALATM',1,natm)
      allocate(NLMOQM1(nbas),NLMOQM2(nbas))
      NLMOQM1=0
      NLMOQM2=0
!******************************************************************
! * Read NLS file detailing QM regions 'QMcenter'
!******************************************************************
      LUNITN = 30
      INQUIRE(FILE='QMcenter',EXIST=NLS_EXIST)
      if (NLS_EXIST) THEN
        open(unit=LUNITN,file='QMcenter',status='old',iostat=ierr)
        if (ierr.eq.0) then
          rewind LUNITN
          read(LUNITN,*) QM1num
          read(LUNITN,*)
          if (natm-QM1num .gt. 0) then
            allocate(QM1atoms(QM1num),QM2atoms(natm-QM1num))
          else 
            allocate(QM1atoms(QM1num))
          endif
          do i=1,QM1num
            read(LUNITN,*) site
            QM1atoms(i)=site
          enddo
          read(LUNITN,*)
          read(LUNITN,*) QM2num
          if (QM2num.gt.0) then
            do i=1,QM2num
              read(LUNITN,*) site
              QM2atoms(i)=site
            enddo
          endif
        else
          print*, "error opening QMcenter"
        endif
        close(LUNITN)  
      else
        print*,'FATAL ERROR: QMcenter-NLS input file-can not be found'
        call exit(1)
      endif

!******************************************************************
! * Read output of N. Flocke NLMO code 'nbocenters'
! * 'NLMO(x,y)' stores linked list of NLMOs associated with atom site
!               where x is the NLMO number and y are the atoms sharing
!               the NLMO
!******************************************************************
      INQUIRE(FILE='nbocenters',EXIST=NBO_EXIST)
      if (NBO_EXIST) THEN
        NBOfile=300
        allocate(NLMO(nbas,5))
        open(unit=NBOfile,file='nbocenters',status='old',iostat=ierr)
        rewind NBOfile
        read(NBOfile,30,IOSTAT=ierr) NLMOnum,NLMOorigin
        NLMO=100
        NLMO(NLMOnum,1)=NLMOorigin
        indx=1
        do while (ierr.eq.0)
           read(NBOfile,30,IOSTAT=ierr) NLMOnum2,NLMOorigin2
           if (ierr.ne.0) exit
           if (NLMOnum.eq.NLMOnum2) then
              indx=indx+1
           else
              indx=1
           endif
           NLMO(NLMOnum2,indx)=NLMOorigin2
           NLMOnum=NLMOnum2
           NLMOorigin2=NLMOorigin
        enddo
        close(NBOfile)
      ELSE
        print*,'FATAL ERROR: nbocenters-output of N. Flocke NLMO code-
     &               can not be found'
        call exit(1)
      endif     

!******************************************************************
        print*, 'INSIDE DRIVE_NLS.F'

!******************************************************************
! * Determine if the NLMO read from nbocenters is entirely within
! * QM1, QM2, or both QM1&QM2
!******************************************************************
        do i=1,nbas
          terminalIndex=findloc(NLMO(i,1:5),value=100,dim=1)
          qm1switch=0
          qm2switch=0
          do j=1,terminalIndex-1
                if (any(NLMO(i,j).eq.QM1atoms)) then
                   qm1switch=1
                else if (any(NLMO(i,j).eq.QM2atoms)) then
                   qm2switch=1
                endif
          enddo
!         * Handles if NLMO is shared between QM1 and QM2
          if ((qm1switch.eq.1) .and. (qm2switch.eq.1)) then
                NLMOQM1(counti)=i
                counti=counti+1

                NLMOQM2(countj)=i
                countj=countj+1
!         * Only if NLMO is in QM1
          else if ((qm1switch.eq.1) .and. (qm2switch.eq.0)) then
                NLMOQM1(counti)=i
                counti=counti+1
!         * Only if NLMO is in QM2
          else if ((qm1switch.eq.0) .and. (qm2switch.eq.1)) then
                NLMOQM2(countj)=i
                countj=countj+1
          endif
        enddo
!#ifdef _DEBUG_LVL0
        print*
        print*,'NLMOQM1', NLMOQM1
        print*,'NLMOQM2',NLMOQM2
!#endif
!******************************************************************
!******************************************************************
!       *       DONE OBTAINING AUXILLARY NLS INFO       *
!       *       NOW MAKE MODIFCATIONS TO GUESS VECTOR   *
!******************************************************************
!******************************************************************

      LS1OUT      = 490
      LS2OUT(1,1) = 444
      LS2OUT(1,2) = 446
      LS2OUT(2,1) = 446
      LS2OUT(2,2) = 445
C
      STRINGI(2)='i  [i_SYM]  '
      STRINGJ(2)='j  [j_SYM]  '
      STRINGA(2)='a  [a_SYM]  '
      STRINGB(2)='b  [b_SYM]  '
      STRINGI(1)='I  [I_SYM]  '
      STRINGJ(1)='J  [J_SYM]  '
      STRINGA(1)='A  [A_SYM]  '
      STRINGB(1)='B  [B_SYM]  '
      CALL ZERO(SCR,NSIZEC)


!******************************************************************
!******************************************************************
!******************************************************************
! * Modified version of 'drive_nto_ee_state_prjct.F' lines 165-233
! * Works on the single excitation coefficient C(i,a) by zeroing out
! * the i not in QM1
! *** FIRST TEST FOR JUST THE UHF CASE, THEN TRY RHF WHICH HAS NEVER 
! *** BEEN DONE IN THE NLS SCHEME 12/1/2020


      IF (PROJECT_SINGLES) THEN
        DO 500 SSPIN = 1, 1+IUHF
          CALL GETLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
          ICOUNT = 1
          DO IIRREP = 1, NIRREP
            AIRREP=DIRPRD(IIRREP,IRREPX)

          DO 1 I = 1, POP(IIRREP, SSPIN)
            DO 2 A = 1, VRT(AIRREP,SSPIN)

            compareIA=0
            QMregIA=0
            compareIA=(/ i,a /)
            call findQMregion(compareIA,size(compareIA),NLMOQM1,
     &           size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIA)

            if (QMregIA(1).eq.2) then
              SCR(ICOUNT) = 0.0D0
            else
              SCR(ICOUNT) =1.0d0
            endif
            print*, SCR(ICOUNT)
!                IF (IIRREP.EQ.INTOO(1) .AND. AIRREP.EQ.INTOV(1)) THEN
!                   IF (I.EQ.INTOO(2) .AND. A.EQ.INTOV(2)) THEN
!                      WRITE(6,'(2(A,I3,A,I1,A))')'R1:',I,'[',IIRREP,']',
!     &                                           ' ->',A,'[',AIRREP,']'
!                      SCR(ICOUNT) = 1.0D0
!                   ELSE
!                      SCR(ICOUNT) = 0.0D0
!                   ENDIF
!                ELSE
!                   SCR(ICOUNT) = 0.0D0
!                ENDIF
              ICOUNT = ICOUNT + 1
    2       CONTINUE
    1     CONTINUE
        ENDDO
        CALL PUTLST(SCR, 1, 1, 1, SSPIN, LS1OUT)
  500 CONTINUE
      ENDIF
      print*, SCR(1:ICOUNT)
30    FORMAT (I4,2X,I4)
!******************************************************************
!******************************************************************
!       * DONE MODIFYING CIS COEFF      *
!******************************************************************
!******************************************************************







!      deallocate(NLMOQM1,NLMOQM2)
          return
          end subroutine

