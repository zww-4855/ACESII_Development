




































































































































































































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
      call getrec(1,'JOBARC','NBASTOT',1,nbas)
      call getrec(1,'JOBARC','NREALATM',1,natm)
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
          allocate(QM1atoms(QM1num),QM2atoms(natm-QM1num))
          do i=1,QM1num
            read(LUNITN,*) site
            QM1atoms(i)=site
          enddo
          read(LUNITN,*)
          read(LUNITN,*) QM2num
          do i=1,QM2num
            read(LUNITN,*) site
            QM2atoms(i)=site
          enddo
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
! lines 203-214 of driveNTO code:

30      FORMAT (I4,2X,I4)
      end subroutine
