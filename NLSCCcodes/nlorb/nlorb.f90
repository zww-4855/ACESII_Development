      PROGRAM NLORB

      IMPLICIT NONE

!      integer :: icore(1)
      INTEGER :: i,j,k,l,a,b,c,d
      INTEGER :: core
      INTEGER :: memscr, iiicrtop, iuhf, iintln, ifltln, iintfp, ialone, ibitwd, jjjcrtop
      INTEGER :: i0, iiicoefs, iiidtran, iiiinvtr, iiiovrlp, iiidenmx, iiidscrs, iiiachrg, iiishell
      INTEGER :: iiinbasa, iiinchoo, iiichoos, iiicoord, iiimomto, iiicbfto
      INTEGER :: iiinshel, iiifockm, iiiordvc, iiiiodvc, iiievals
      INTEGER :: jjj00000, jjjcbfnt, jjjmomnt, jjjcbfto, jjjmomto
      INTEGER :: jjjsjunk, jjjcntsc, jjjcnts2, jjjiscrs, jjjtrans, jjjrowss, jjjscrsp, jjjpivot

      INTEGER :: nbas, nbastot, nelec, natoms, nocc(2), nvir(2)
      INTEGER :: bondsize, mxchoose, mxnal, largel, angmax, switch12
      INTEGER :: istat

      DOUBLE PRECISION :: a2u, enuc, maxocc, eelec, det, acesener
      DOUBLE PRECISION :: occthresh, virthresh

      LOGICAL :: ryd2hyb, mjump, spheric, onlynao, onlynho, onlynbo, onlynlmo, verprn, putaces, setchoose
      LOGICAL :: creatlist, yesno

      COMMON /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /istart/ i0, jjjcrtop 
      COMMON /iending/ jjj00000, iiicrtop

!---------------------------------------------------------------------------------
! iintln : the byte-length of a default integer
! ifltln : the byte-length of a double precision float
! iintfp : the number of integers in a double precision float
! ialone : the bitmask used to filter out the lowest fourth bits in aninteger
! ibitwd : the number of bits in one-fourth of an integer
!---------------------------------------------------------------------------------

      integer,allocatable,dimension(:)::icore

      memscr   = 1000000
      angmax   = 20
      ryd2hyb  = .FALSE.
      mjump    = .TRUE.
      spheric  = .TRUE.
      onlynao  = .FALSE.
      onlynho  = .FALSE.
      onlynbo  = .FALSE.
      onlynlmo = .TRUE.
      verprn   = .FALSE.
      mxchoose = 1
      a2u      = 1.88972599D0
      putaces  = .TRUE.
      setchoose= .FALSE.
      creatlist= .FALSE.

      iuhf = 0
      core = 0
      CALL aces_init(core,i0,jjjcrtop,iuhf,.TRUE.)


!******************************
      i0=1
!******************************

      print*, i0, jjjcrtop, iuhf

      iiicrtop = jjjcrtop - memscr
      jjj00000 = iiicrtop

      write(*,*)
      write(*,*) "========================================"
      write(*,*) "Natural Localized Orbital ACES Interface"
      write(*,*) "      written by Thomas F. Hughes       "
      write(*,*) "========================================"
      write(*,*)

      open(unit=100,file='bondsize',status='old',iostat=istat)
      if (istat.eq.0) then
       read(unit=100,fmt=*) bondsize
      else
       bondsize=4
      end if
      close(100)


!      open(unit=200,file='threshold')
!         read(unit=200,fmt=*) occthresh
!         read(unit=200,fmt=*) virthresh
!      close(200)
!      write(*,*) 'occ threshold = ', occthresh
!      write(*,*) 'vir threshold = ', virthresh
 
      CALL GETREC(0,'JOBARC','12SWITCH',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         go to 950
      end if
      nbas = 0
      CALL GETREC(20,'JOBARC','12SWITCH',1,switch12)
      if (switch12.eq.1) then
         write(*,*)
         write(*,*) "THE FIRST TWO ATOMS OF ZMAT HAVE BEEN SWAPPED"
         write(*,*) "FIX ZMAT ORDERING"
         write(*,*) "--- Localization will exit ---"
         write(*,*)
         stop
      end if

 950  continue

      CALL initdim(nbas,nbastot,nocc,nvir,natoms,nelec)
 
      iiicoord = i0
      CALL CKMEM(iiicoord,iiicrtop)
      iiiachrg = iiicoord+natoms*3*iintfp
      CALL CKMEM(iiiachrg,iiicrtop)
      jjjcbfnt = jjj00000
      CALL CKMEM(jjjcbfnt,jjjcrtop)
      jjjmomnt = jjjcbfnt+nbas
      CALL CKMEM(jjjmomnt,jjjcrtop)

!************************************************
      allocate (icore(jjjcrtop))
!      icore(iiicoord)=0
!      icore(iiiachrg)=0
!      icore(jjjmomnt)=0
!      icore(jjjcbfnt)=0
      print *, iiicoord, iiiachrg, jjjmomnt, jjjcbfnt
      print *, icore(iiicoord), icore(iiiachrg)
      print *, icore(jjjmomnt), icore(jjjcbfnt)
!************************************************

      CALL initize(icore(iiicoord),icore(iiiachrg),icore(jjjmomnt),icore(jjjcbfnt),natoms,nbas)

!-------------------------------------------------------------------------
!  icore(iiicoord): atomic coordinates
!  icore(iiiachrg): atomic charge 
!  icore(jjjmomnt): l-quantum number of each AO basis function
!  icore(jjjcbfnt): atomic center to which each AO basis function belongs
!-------------------------------------------------------------------------

      jjjcbfto = jjjmomnt+nbas
      CALL CKMEM(jjjcbfto,jjjcrtop)

      jjjmomto = jjjcbfto+nbas
      CALL CKMEM(jjjmomto,jjjcrtop)

      CALL iseccopy(icore(jjjcbfnt),icore(jjjcbfto),nbas)
      CALL iseccopy(icore(jjjmomnt),icore(jjjmomto),nbas)
!-----------------------------------------------------------------------------
! copy icore(jjjcbfnt) to icore(jjjcbfto), icore(jjjmomnt) to icore(jjjmomto)
!-----------------------------------------------------------------------------

      if (nbastot.ne.nbas) then
         CALL izerosec(icore(jjjcbfto),nbas)
         CALL izerosec(icore(jjjmomto),nbas)

         jjjcbfto = jjjmomnt+nbas
          CALL CKMEM(jjjcbfto,jjjcrtop)
         jjjmomto = jjjcbfto+nbastot
          CALL CKMEM(jjjmomto,jjjcrtop)
         jjjrowss = jjjmomto+nbastot
          CALL CKMEM(jjjrowss,jjjcrtop)
         jjjtrans = jjjrowss+nbas-nbastot
          CALL CKMEM(jjjtrans,jjjcrtop)
         jjjscrsp = jjjtrans+nbastot*nbas
          CALL CKMEM(jjjscrsp,jjjcrtop)

         CALL doitran(icore(jjjrowss),icore(jjjtrans), &
      icore(jjjscrsp),icore(jjjcbfnt),icore(jjjmomnt), &
      icore(jjjcbfto),icore(jjjmomto),nbas,nbastot,angmax)
      end if

      iiicbfto = iiiachrg+natoms
      CALL CKMEM(iiicbfto,iiicrtop)

      iiimomto = iiicbfto+nbastot
      CALL CKMEM(iiimomto,iiicrtop)

      CALL iseccopy(icore(jjjcbfto),icore(iiicbfto),nbastot)
      CALL iseccopy(icore(jjjmomto),icore(iiimomto),nbastot)

      CALL prncoor(natoms,icore(iiicoord),a2u,icore(iiiachrg))
!---------------------------------------------------------------------
!  prncoor: print the cartesian coordinates with the unit of Bohr
!---------------------------------------------------------------------

      CALL prnang(largel,nbastot,icore(iiimomto),icore(iiicbfto))
!---------------------------------------------------------------------
!  prnang: print the l-quantum number for each basis function
!  largel is the largest l (s=0, p=1, d=2,...)
!---------------------------------------------------------------------
     
      iiinshel = iiimomto+nbastot
       CALL CKMEM(iiinshel,iiicrtop)
      iiishell = iiinshel+natoms
       CALL CKMEM(iiishell,iiicrtop)
      iiinbasa = iiishell+(largel+1)*natoms
       CALL CKMEM(iiinbasa,iiicrtop)

      CALL izerosec(icore(jjj00000),nbas*4+nbas-nbastot+ &
      nbastot*nbas+nbas*nbas)

      jjjsjunk = jjj00000
       CALL CKMEM(jjjsjunk,jjjcrtop)
      jjjcntsc = jjjsjunk+nbastot
       CALL CKMEM(jjjcntsc,jjjcrtop)

      jjjcnts2 = jjjcntsc+natoms*angmax
!                              (angmax=20)
       CALL CKMEM(jjjcnts2,jjjcrtop)

      CALL sortmom(icore(iiicbfto),icore(iiimomto),icore(iiinshel), &
      icore(iiishell),icore(iiinbasa),largel,mxnal,natoms,nbastot, &
      angmax,icore(jjjsjunk),icore(jjjcntsc),icore(jjjcnts2))
!-------------------------------------------------------------------------
!  sortmom: summarize the number of shells on each atom and the number of
!           different type of angular momentum on each shell
!-------------------------------------------------------------------------

       iiidtran = iiinbasa+(largel+1)*natoms
      CALL CKMEM(iiidtran,iiicrtop)
       iiiinvtr = iiidtran+nbastot*nbastot*iintfp
      CALL CKMEM(iiiinvtr,iiicrtop)
       iiiordvc = iiiinvtr+nbastot*nbastot*iintfp
      CALL CKMEM(iiiordvc,iiicrtop)
       iiiiodvc = iiiordvc+nbastot
      CALL CKMEM(iiiiodvc,iiicrtop)

      CALL dodtrans(icore(iiidtran),icore(iiiinvtr), &
      icore(iiimomto),icore(iiicbfto),icore(iiinbasa), &
      icore(iiinshel),icore(iiishell), &
      largel,natoms,angmax,nbastot)

      call zpri_den("Owntrans",NBastot,nbastot,nbastot,nbastot,icore(iiidtran))
      call zpri_den("Owntrans inv",NBastot,nbastot,nbastot,nbastot,icore(iiiinvtr))

      CALL PRINT_MATRIX('Owntrans            ',icore(iiidtran),nbastot, &
      nbastot,verprn)
      CALL PRINT_MATRIX('Owntrans Inverse    ',icore(iiiinvtr),nbastot, &
      nbastot,verprn)

      CALL detord(icore(iiidtran),icore(iiiordvc),nbastot)
      CALL prnord(icore(iiiordvc),nbastot)
      CALL detord(icore(iiiinvtr),icore(iiiiodvc),nbastot)
      CALL prnord(icore(iiiiodvc),nbastot)


      write(*,*)
      write(*,*) "ACES: Initializing NBO-type Localization"
      write(*,*)
 
       iiicoefs = iiiiodvc+nbastot
      CALL CKMEM(iiicoefs,iiicrtop)
       iiiovrlp = iiicoefs+nbastot*nbastot*iintfp*(iuhf+1)
      CALL CKMEM(iiiovrlp,iiicrtop)
       iiidenmx = iiiovrlp+nbastot*nbastot*iintfp
      CALL CKMEM(iiidenmx,iiicrtop)

      CALL GETREC(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp, &
      icore(iiicoefs))
      CALL PRINT_MATRIX('Canonical Coef Alpha',icore(iiicoefs),nbastot, &
      nbastot,.TRUE.)

      CALL MAKE_DENSITY_MATRIX(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc), &
      nbastot,nocc(1),.FALSE.)
      if (iuhf.eq.0) then
         CALL scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      CALL PRINT_MATRIX('Density Alpha       ',icore(iiidenmx),nbastot, &
      nbastot,.TRUE.)

      if (iuhf.eq.1) then
         CALL GETREC(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp, &
      icore(iiicoefs+nbastot*nbastot*iintfp))
         CALL PRINT_MATRIX('Canonical Coef Beta ',icore(iiicoefs+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,.TRUE.)

         CALL MAKE_DENSITY_MATRIX(icore(iiicoefs+nbastot*nbastot*iintfp), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc), &
      nbastot,nocc(2),.FALSE.)
         CALL PRINT_MATRIX('Density Beta        ',icore(iiidenmx+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,.TRUE.)
      end if

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp)) 
      CALL GETREC(20,'JOBARC','AOOVRLAP',nbastot*nbastot*iintfp, &
      icore(iiiovrlp))
      CALL PRINT_MATRIX('Overlap             ',icore(iiiovrlp),nbastot, &
      nbastot,.TRUE.)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL transform2(icore(iiidtran),icore(iiidenmx), &
      icore(iiicoefs),nbastot)
      CALL PRINT_MATRIX('Canon Coef New Alpha',icore(iiicoefs),nbastot, &
      nbastot,verprn)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL MAKE_DENSITY_MATRIX(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc), &
      nbastot,nocc(1),.FALSE.)
      if (iuhf.eq.0) then
         CALL scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      CALL PRINT_MATRIX('Density New Alpha   ',icore(iiidenmx),nbastot, &
      nbastot,verprn)

      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+ &
      nbastot*nbastot*iintfp))
        CALL transform2(icore(iiidtran),icore(iiidenmx+ &
      nbastot*nbastot*iintfp),icore(iiicoefs+nbastot*nbastot* &
      iintfp),nbastot)
        CALL PRINT_MATRIX('Canon Coef New Beta ',icore(iiicoefs+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+ &
      nbastot*nbastot*iintfp))
        CALL MAKE_DENSITY_MATRIX(icore(iiicoefs+nbastot*nbastot*iintfp), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc), &
      nbastot,nocc(2),.FALSE.)
        CALL PRINT_MATRIX('Density New Beta    ',icore(iiidenmx+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if

      CALL GETREC(0,'JOBARC','CORRDENA',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         switch12 = 0
         go to 951
      end if
      nbas = 0
      switch12 = 0
      write(*,*)
      write(*,*) "Using coupled cluster response density alpha"
      write(*,*)
      CALL dzerosec(nbastot*nbastot,icore(iiidenmx))
      CALL GETREC(21,'JOBARC','CORRDENA',nbastot*nbastot, &
      icore(iiidenmx))
!
      CALL PRINT_MATRIX('CC Response AO densA',icore(iiidenmx),nbastot, &
      nbastot,.TRUE.)
      CALL dzerosec(nbastot*nbastot,icore(jjj00000))
      CALL transform(icore(iiidtran),icore(jjj00000), &
      icore(iiidenmx),nbastot)
      CALL PRINT_MATRIX('New CC Resp AO densA',icore(iiidenmx),nbastot, &
      nbastot,verprn)
      CALL dzerosec(nbastot*nbastot,icore(jjj00000))
      if (iuhf.eq.1) then
         write(*,*)
         write(*,*) "Using coupled cluster response density beta"
         write(*,*)
         CALL dzerosec(nbastot*nbastot,icore(iiidenmx+ &
      nbastot*nbastot))
        CALL GETREC(20,'JOBARC','CORRDENB',nbastot*nbastot, &
      icore(iiidenmx+nbastot*nbastot))
        CALL PRINT_MATRIX('CC Response AO densB',icore(iiidenmx+ &
      nbastot*nbastot),nbastot,nbastot,.TRUE.)
        CALL dzerosec(nbastot*nbastot,icore(jjj00000))
        CALL transform(icore(iiidtran),icore(jjj00000), &
      icore(iiidenmx+nbastot*nbastot),nbastot)
        CALL PRINT_MATRIX('New CC Resp AO densB',icore(iiidenmx+ &
      nbastot*nbastot),nbastot,nbastot,verprn)
         CALL dzerosec(nbastot*nbastot,icore(jjj00000))
      end if

 951  continue

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
      CALL transform(icore(iiidtran),icore(iiicoefs), &
      icore(iiiovrlp),nbastot)
      CALL PRINT_MATRIX('Overlap New         ',icore(iiiovrlp),nbastot, &
      nbastot,verprn)

       iiinchoo = iiidenmx+nbastot*nbastot*iintfp*(iuhf+1)
      CALL CKMEM(iiinchoo,iiicrtop)
       iiichoos= iiinchoo+bondsize
      CALL CKMEM(iiichoos,iiicrtop)
       iiidscrs = iiichoos+bondsize*bondsize*mxchoose
      CALL CKMEM(iiidscrs,iiicrtop)

      CALL izerosec(icore(jjj00000),nbastot+2*angmax*natoms)

      jjjiscrs = jjj00000
      CALL CKMEM(jjjiscrs,jjjcrtop)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs+ &
      nbastot*nbastot*iintfp))
      end if

      CALL initmax(maxocc,iuhf)

      write(*,*)
      write(*,*) "Entering NBO-type Localization written by Norbert Flocke"
      write(*,*)

!      CALL PUTREC(20,'JOBARC','IIINSHEL',natoms,icore(iiinshel))
!      CALL PUTREC(20,'JOBARC','IIISHELL',(largel+1)*natoms, &
!      icore(iiishell))
!      CALL PUTREC(20,'JOBARC','IIINBASA',(largel+1)*natoms, &
!       icore(iiinbasa))
!      CALL PUTREC(20,'JOBARC','IIIACHRG',natoms,icore(iiiachrg))
!      CALL PUTREC(20,'JOBARC','IIINCHOO',bondsize,icore(iiinchoo))
!      CALL PUTREC(20,'JOBARC','IIICHOOS',bondsize*mxchoose*bondsize, &
!       icore(iiichoos))
      CALL PUTREC(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp, &
      icore(iiiovrlp))

      if (iuhf.eq.0) then
         write(*,*)
         write(*,*)
         write(*,*) "--- LOCALIZING WITH RHF REFERENCE ---"
         write(*,*)
         write(*,*)
      end if
      if (iuhf.eq.1) then
         write(*,*)
         write(*,*)
         write(*,*) "--- LOCALIZING WITH UHF REFERENCE ALPHA ---"
         write(*,*)
         write(*,*)
      end if
 
      if (setchoose) then
         CALL setchoo(icore(iiinchoo),icore(iiichoos),mxchoose, &
      bondsize)
      end if


!      open (unit=40, file='nlscc4.dat')
!      write (40,*) jjjcrtop-jjjiscrs
!      write (40,*) iiicrtop-iiidscrs
!      write (40,*) nbastot
!      write (40,*) natoms
!      write (40,*) largel
!      write (40,*) mxnal
!      write (40,*) maxocc
!      write (40,*) bondsize
!      write (40,*) ryd2hyb
!      write (40,*) icore(iiinshel)
!      write (40,*) icore(iiishell)
!      write (40,*) icore(iiinbasa)
!      write (40,*) icore(iiiachrg)
!      write (40,*) icore(iiidenmx)
!      write (40,*) icore(iiiovrlp)
!      write (40,*) mjump
!      write (40,*) spheric
!      write (40,*) onlynao
!      write (40,*) onlynho
!      write (40,*) onlynbo
!      write (40,*) onlynlmo
!      write (40,*) mxchoose
!      write (40,*) icore(iiinchoo)
!      write (40,*) icore(iiichoos)
!      write (40,*) icore(jjjiscrs)
!      write (40,*) icore(iiidscrs)
!      write (40,*) icore(iiicoefs)
!      close (unit=40)


      CALL nlo__gener_nlo_orbitals(1,jjjcrtop-jjjiscrs, & 
      iiicrtop-iiidscrs,nbastot,natoms,largel,mxnal,maxocc, &
      bondsize,ryd2hyb,icore(iiinshel),icore(iiishell), &
      icore(iiinbasa),icore(iiiachrg),icore(iiidenmx), &
      icore(iiiovrlp),mjump,spheric,onlynao,onlynho,onlynbo, &
      onlynlmo,mxchoose,icore(iiinchoo),icore(iiichoos), &
      icore(jjjiscrs),icore(iiidscrs),icore(iiicoefs))


! yjin
!      call zpri_den("CoeF MATRIX 1",NBastot,nbastot,nbastot,nbastot,icore(iiicoefs))

!  1                      first
!  jjjcrtop-jjjiscrs      IMAX
!  iiicrtop-iiidscrs      ZMAX
!  nbastot                NBAS
!  natoms                 NATOM
!  largel                 MAXSHELL
!  mxnal                  MXNAL
!  icore(iiinshel)        NSHELLS
!  icore(iiishell)        SHELLS
!  icore(iiinbasa)        NBASAL
!  icore(iiiachrg)        ZATOM
!  icore(iiidenmx)        DENSITY
!  icore(iiiovrlp)        OVERLAP
!  icore(iiinchoo)        NCHOOSE
!  icore(iiichoos)        CHOOSE
!  icore(jjjiscrs)        ICORE(1:IMAX)
!  icore(iiidscrs)        ZCORE(1:ZMAX)
!  icore(iiicoefs)        COEFFS

      if (iuhf.eq.1) then

!         CALL GETREC(20,'JOBARC','IIINSHEL',natoms,icore(iiinshel))
!         CALL GETREC(20,'JOBARC','IIISHELL',(largel+1)*natoms, &
!      icore(iiishell))
!        CALL GETREC(20,'JOBARC','IIINBASA',(largel+1)*natoms, &
!      icore(iiinbasa))
!        CALL GETREC(20,'JOBARC','IIIACHRG',natoms,icore(iiiachrg))
!        CALL GETREC(20,'JOBARC','IIINCHOO',bondsize,icore(iiinchoo))
!        CALL GETREC(20,'JOBARC','IIICHOOS',bondsize*mxchoose*bondsize, &
!      icore(iiichoos))
        CALL GETREC(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp, &
      icore(iiiovrlp))

         write(*,*)
         write(*,*)
         write(*,*) "--- LOCALIZING WITH UHF REFERENCE BETA  ---"
         write(*,*)
         write(*,*)

         CALL nlo__gener_nlo_orbitals(0,jjjcrtop-jjjiscrs, &
      iiicrtop-iiidscrs,nbastot,natoms,largel,mxnal,maxocc, &
      bondsize,ryd2hyb,icore(iiinshel),icore(iiishell), &
      icore(iiinbasa),icore(iiiachrg), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiovrlp), &
      mjump,spheric,onlynao,onlynho,onlynbo,onlynlmo,mxchoose, &
      icore(iiinchoo),icore(iiichoos),icore(jjjiscrs), &
      icore(iiidscrs),icore(iiicoefs+nbastot*nbastot*iintfp))

      end if

      write(*,*)
      write(*,*) "Exiting NBO-type Localization"
      write(*,*)

!      CALL PRINT_MATRIX('Local New Alpha     ',icore(iiicoefs),nbastot, &
!      nbastot,verprn)

      call zpri_den("Local New Alpha",nbastot,nbastot,nbastot,nbastot,icore(iiicoefs))

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL MAKE_DENSITY_MATRIX(icore(iiicoefs),icore(iiidenmx),icore(iiiordvc), &
      nbastot,nocc(1),.FALSE.)
      if (iuhf.eq.0) then
         CALL scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      CALL PRINT_MATRIX('Local Dens New Alpha',icore(iiidenmx),nbastot, &
      nbastot,verprn)

      if (iuhf.eq.1) then
         CALL PRINT_MATRIX('Local New Beta      ',icore(iiicoefs+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

        CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+ &
      nbastot*nbastot*iintfp))
        CALL MAKE_DENSITY_MATRIX(icore(iiicoefs+nbastot*nbastot*iintfp), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiordvc), &
      nbastot,nocc(2),.FALSE.)
        CALL PRINT_MATRIX('Local Dens New Beta ',icore(iiidenmx+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp))
      CALL GETREC(20,'JOBARC','MYOVRLP ',nbastot*nbastot*iintfp, &
      icore(iiiovrlp))
      CALL PRINT_MATRIX('Overlap New         ',icore(iiiovrlp),nbastot, &
      nbastot,verprn)

! yjin
!      call zpri_den("CoeF MATRIX 3",NBastot,nbastot,nbastot,nbastot,icore(iiicoefs))

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL transform2(icore(iiiinvtr),icore(iiidenmx), &
      icore(iiicoefs),nbastot)
      CALL chksgn(icore(iiicoefs),nbastot) 


      inquire(file='orbsfile',exist=yesno)
      if (yesno) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiicoefs))
         CALL readorbfile(nbastot,icore(iiicoefs))
      end if

      CALL PRINT_MATRIX('Localized Alpha     ',icore(iiicoefs),nbastot, &
      nbastot,.TRUE.)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL GETREC(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp, &
      icore(iiidenmx))
      write(*,31) putaces
      if (.not.putaces) then
         write(*,*) "Non-canonical NLMOS will _NOT_ be put into JOBARC"
      end if
      if (putaces) then
         write(*,*)
         write(*,*)
         write(*,*) "Putting non-canonical NLMOS into JOBARC"
         write(*,*)
         write(*,*)
         CALL PUTREC(20,'JOBARC','SCFEVCA0',nbastot*nbastot*iintfp, &
      icore(iiicoefs))
! yjin
!      call zpri_den("CoeF MATRIX 2",NBastot,nbastot,nbastot,nbastot,icore(iiicoefs))

         CALL GETREC(0,'JOBARC','CORRDENA',nbas,switch12)
         if (nbas.eq.-1) then
            nbas = 0
            switch12 = 0
            go to 953
         end if
         nbas = 0
         switch12 = 0
         if (creatlist) then
            CALL orbfile(nbastot,icore(iiicoefs))
         end if
 953  continue

!         CALL PUTREC(20,'JOBARC','SCFCANA0',nbastot*nbastot*iintfp, &
!       icore(iiidenmx))
      end if
!      CALL PUTREC(20,'JOBARC','SCFLOCA0',nbastot*nbastot*iintfp, &
!       icore(iiicoefs))

      jjjpivot = jjj00000
      CALL CKMEM(jjjpivot,jjjcrtop)

      CALL izerosec(icore(jjjpivot),nbastot)
      CALL minv(icore(iiidenmx),nbastot,nbastot,icore(jjjpivot), &
      det,1.0D-8,0,0)
      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      CALL xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,icore(iiidenmx), &
      nbastot,icore(iiicoefs),nbastot,0.0D0,icore(iiidscrs),nbastot)
      CALL PRINT_MATRIX('Can Non Trans Alpha ',icore(iiidscrs),nbastot, &
      nbastot,.TRUE.)
!      CALL PUTREC(20,'JOBARC','CANNONAL',nbastot*nbastot*iintfp, &
!      icore(iiidscrs))

! yjin
!      call zpri_den("Fock matrix 2",nbastot,nbastot,nbastot,nbastot,icore(iiidscrs))

      CALL GETREC(0,'JOBARC','CORRDENA',nbas,switch12)
      if (nbas.eq.-1) then
         nbas = 0
         switch12 = 0
         go to 952
      end if
      nbas = 0
      switch12 = 0
      if (creatlist) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
         CALL droplist(nocc(1),nvir(1),nbastot,icore(iiidscrs), &
      occthresh,virthresh,icore(iiidenmx),icore(iiidenmx+nbastot))
      end if

 952  continue

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL izerosec(icore(jjjpivot),nbastot)
      CALL dtranspose(icore(iiidscrs),icore(iiidenmx),nbastot)
      CALL xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,icore(iiidenmx), &
      nbastot,icore(iiidscrs),nbastot,0.0D0,icore(jjjpivot),nbastot)
      CALL PRINT_MATRIX('Trans(T)Trans Alpha ',icore(jjjpivot),nbastot, &
      nbastot,.TRUE.)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
      CALL MAKE_DENSITY_MATRIX(icore(iiicoefs),icore(iiidenmx),icore(iiiiodvc), &
      nbastot,nocc(1),.FALSE.)
      if (iuhf.eq.0) then
         CALL scalmat(icore(iiidenmx),nbastot,2.0D0)
      end if
      CALL PRINT_MATRIX('Localized Dens Alpha',icore(iiidenmx),nbastot, &
      nbastot,verprn)

      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp, &
      icore(iiidenmx+nbastot*nbastot*iintfp))
         CALL transform2(icore(iiiinvtr), &
      icore(iiidenmx+nbastot*nbastot*iintfp), &
      icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         CALL chksgn(icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         CALL PRINT_MATRIX('Localized Beta      ',icore(iiicoefs+nbastot* &
      nbastot*iintfp),nbastot,nbastot,.TRUE.)

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+ &
      nbastot*nbastot*iintfp))
         CALL GETREC(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp, &
      icore(iiidenmx+nbastot*nbastot*iintfp))
         write(*,31) putaces
         if (.not.putaces) then
            write(*,*) "Non-canonical NLMOS will _NOT_ be put into JOBARC"
         end if
         if (putaces) then
            write(*,*)
            write(*,*)
            write(*,*) "Putting non-canonical NLMOS into JOBARC"
            write(*,*)
            write(*,*)
            CALL PUTREC(20,'JOBARC','SCFEVCB0',nbastot*nbastot*iintfp, &
      icore(iiicoefs+nbastot*nbastot*iintfp))
            CALL PUTREC(20,'JOBARC','SCFCANB0',nbastot*nbastot*iintfp, &
      icore(iiidenmx+nbastot*nbastot*iintfp))
         end if
      CALL PUTREC(20,'JOBARC','SCFLOCB0',nbastot*nbastot*iintfp, &
      icore(iiicoefs+nbastot*nbastot*iintfp))

         CALL izerosec(icore(jjjpivot),nbastot)
         CALL minv(icore(iiidenmx+nbastot*nbastot*iintfp), &
      nbastot,nbastot,icore(jjjpivot),det,1.0D-8,0,0)
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         CALL xgemm('N','N',nbastot,nbastot,nbastot,1.0D0, &
      icore(iiidenmx+nbastot*nbastot*iintfp),nbastot, &
      icore(iiicoefs+nbastot*nbastot*iintfp),nbastot, &
      0.0D0,icore(iiidscrs),nbastot)
         CALL PRINT_MATRIX('Can Non Trans Beta  ',icore(iiidscrs),nbastot, &
      nbastot,.TRUE.)
         CALL PUTREC(20,'JOBARC','CANNONBE',nbastot*nbastot*iintfp, &
      icore(iiidscrs))

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx+ &
      nbastot*nbastot*iintfp))
         CALL izerosec(icore(jjjpivot),nbastot)
         CALL dtranspose(icore(iiidscrs),icore(iiidenmx+ &
      nbastot*nbastot*iintfp),nbastot)
         CALL xgemm('N','N',nbastot,nbastot,nbastot,1.0D0, &
      icore(iiidenmx+nbastot*nbastot*iintfp),nbastot, &
      icore(iiidscrs),nbastot,0.0D0,icore(jjjpivot),nbastot)
         CALL PRINT_MATRIX('Trans(T)Trans Beta  ',icore(jjjpivot),nbastot, &
      nbastot,.TRUE.)

         CALL dzerosec(nbastot*nbastot*iintfp, &
      icore(iiidenmx+nbastot*nbastot*iintfp))
         CALL MAKE_DENSITY_MATRIX(icore(iiicoefs+nbastot*nbastot*iintfp), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiiiodvc), &
      nbastot,nocc(2),.FALSE.)
         CALL PRINT_MATRIX('Localized Dens Beta ',icore(iiidenmx+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)
      end if

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      CALL transform(icore(iiiinvtr),icore(iiidscrs), &
      icore(iiiovrlp),nbastot)
      CALL PRINT_MATRIX('Overlap             ',icore(iiiovrlp),nbastot, &
      nbastot,verprn)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      CALL ckorth(icore(iiiovrlp),icore(iiidscrs),icore(iiicoefs), &
      nbastot)
      CALL PRINT_MATRIX('Overlap in MOS Alpha ',icore(iiidscrs),nbastot, &
      nbastot,.TRUE.)

      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         CALL ckorth(icore(iiiovrlp),icore(iiidscrs), &
      icore(iiicoefs+nbastot*nbastot*iintfp),nbastot)
         CALL PRINT_MATRIX('Overlap in MOS Beta  ',icore(iiidscrs), &
      nbastot,nbastot,.TRUE.)
      end if

      iiifockm = iiidscrs+nbastot*nbastot*iintfp
      CALL CKMEM(iiifockm,iiicrtop)

      iiievals = iiifockm+nbastot*nbastot*iintfp*(iuhf+1)
      CALL CKMEM(iiievals,iiicrtop)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiifockm))

      natoms = 0
      CALL GETREC(0,'JOBARC','ONEHAO  ',natoms, &
      icore(iiifockm))
      CALL GETREC(20,'JOBARC','ONEHAO  ',natoms, &
      icore(iiifockm))
      natoms = 0
      CALL unpkoneh(icore(iiifockm),icore(iiidscrs),nbastot)
      CALL PRINT_MATRIX('Core Hamiltonian    ',icore(iiidscrs), &
      nbastot,nbastot,verprn)

      CALL dzerosec(nbastot*nbastot*iintfp,icore(iiifockm))
      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiifockm+ &
      nbastot*nbastot*iintfp))
      end if

      if (iuhf.eq.0) then
         CALL GETREC(20,'JOBARC','FOCKA   ',nbastot*nbastot*iintfp, &
      icore(iiifockm))
         CALL PRINT_MATRIX('Fock Alpha          ',icore(iiifockm), &
      nbastot,nbastot,verprn)
!      call zpri_den("Fock0",nbastot,nbastot,nbastot,nbastot,icore(iiifockm))

         CALL eneval(icore(iiidenmx),icore(iiidscrs), &
      icore(iiifockm),nbastot,eelec)
         
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidenmx))
         CALL orthock(icore(iiidenmx),icore(iiidscrs), &
      icore(iiicoefs),nbastot)
         CALL PRINT_MATRIX('Loc Core in MO Alpha',icore(iiidenmx), &
      nbastot,nbastot,.TRUE.)
!         CALL PUTREC(20,'JOBARC','LOCCOREA',nbastot*nbastot*iintfp, &
!      icore(iiidenmx))

! yjin
!      call zpri_den("Fock matrix 3",nbastot,nbastot,nbastot,nbastot,icore(iiidscrs))

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         CALL orthock(icore(iiidscrs),icore(iiifockm), &
      icore(iiicoefs),nbastot)
         CALL PRINT_MATRIX('Loc Fock in MO Alpha',icore(iiidscrs), &
      nbastot,nbastot,.TRUE.)
         CALL PUTREC(20,'JOBARC','LOCFOCKA',nbastot*nbastot*iintfp, &
      icore(iiidscrs))

! yjin
!      call zpri_den("Fock matrix 4",nbastot,nbastot,nbastot,nbastot,icore(iiidscrs))

         CALL dzerosec(2*nocc(1)*iintfp,icore(jjj00000))
         CALL GETREC(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         CALL hfperbond(icore(iiidenmx),icore(iiidscrs), &
      nbastot,nocc(1),icore(jjj00000),icore(jjj00000+nocc(1)),enuc)
!         CALL PUTREC(20,'JOBARC','MFBONDEN',nocc(1)*iintfp, &
!      icore(jjj00000))

         CALL dzerosec(nbastot*iintfp,icore(iiievals))
         CALL doeval(icore(iiidscrs),icore(iiievals),nbastot)
         CALL PRINT_MATRIX('Loc Fock Diag Alpha ',icore(iiievals), &
      nbastot,1,.TRUE.)

! yjin
         call putrec(20,'JOBARC','FOCKAD',nbastot,icore(iiievals))

         acesener = 0.0D0
         CALL getfia(nbastot,nocc(1),nvir(1),icore(iiidscrs), &
      acesener)
         acesener = 0.0D0

         CALL GETREC(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         write(*,*)
         write(*,19) eelec
         write(*,20) enuc
         write(*,21) enuc + eelec
         write(*,*)
         CALL GETREC(20,'JOBARC','SCFENEG ',iintfp,acesener)
         if (abs(acesener-(enuc+eelec)).gt.1.0D-6) then
            write(*,48) acesener
            write(*,*)
            write(*,*) "Localized Orbitals May NOT be Non-Canonical RHF"
!            stop
         end if
      end if

      if (iuhf.eq.1) then
         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiiovrlp))
         CALL totalden(icore(iiidenmx),icore(iiidenmx+ &
      nbastot*nbastot*iintfp),icore(iiiovrlp),nbastot)

         CALL GETREC(20,'JOBARC','FOCKA   ',nbastot*nbastot*iintfp, &
      icore(iiifockm))
         CALL PRINT_MATRIX('Fock Alpha          ',icore(iiifockm), &
      nbastot,nbastot,verprn)

         CALL GETREC(20,'JOBARC','FOCKB   ',nbastot*nbastot*iintfp, &
      icore(iiifockm+nbastot*nbastot*iintfp))
         CALL PRINT_MATRIX('Fock Beta           ',icore(iiifockm+ &
      nbastot*nbastot*iintfp),nbastot,nbastot,verprn)

         CALL eneval2(icore(iiiovrlp),icore(iiidenmx), &
      icore(iiidenmx+nbastot*nbastot*iintfp),icore(iiidscrs), &
      icore(iiifockm),icore(iiifockm+nbastot*nbastot*iintfp), &
      nbastot,eelec)

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         CALL orthock(icore(iiidscrs),icore(iiifockm), &
      icore(iiicoefs),nbastot)
         CALL PRINT_MATRIX('Loc Fock in MO Alpha',icore(iiidscrs), &
      nbastot,nbastot,.TRUE.)
         CALL PUTREC(20,'JOBARC','LOCFOCKA',nbastot*nbastot*iintfp, &
      icore(iiidscrs))

         CALL dzerosec(nbastot*iintfp,icore(iiievals))
         CALL doeval(icore(iiidscrs),icore(iiievals),nbastot)
         CALL PRINT_MATRIX('Loc Fock Diag Alpha ',icore(iiievals), &
      nbastot,1,.TRUE.)

         acesener = 0.0D0
         CALL getfia(nbastot,nocc(1),nvir(1),icore(iiidscrs), &
      acesener)
         acesener = 0.0D0

         CALL dzerosec(nbastot*nbastot*iintfp,icore(iiidscrs))
         CALL orthock(icore(iiidscrs),icore(iiifockm+ &
      nbastot*nbastot*iintfp),icore(iiicoefs+nbastot*nbastot*iintfp), &
      nbastot)
         CALL PRINT_MATRIX('Loc Fock in MO Beta ',icore(iiidscrs), &
      nbastot,nbastot,.TRUE.)
         CALL PUTREC(20,'JOBARC','LOCFOCKB',nbastot*nbastot*iintfp, &
      icore(iiidscrs))

         CALL dzerosec(nbastot*iintfp,icore(iiievals))
         CALL doeval(icore(iiidscrs),icore(iiievals),nbastot)
         CALL PRINT_MATRIX('Loc Fock Diag Beta  ',icore(iiievals), &
      nbastot,1,.TRUE.)

         acesener = 0.0D0
         CALL getfia(nbastot,nocc(2),nvir(2),icore(iiidscrs), &
      acesener)
         acesener = 0.0D0

         CALL GETREC(20,'JOBARC','NUCREP  ',1*iintfp,enuc)
         write(*,*)
         write(*,19) eelec
         write(*,20) enuc
         write(*,21) enuc + eelec
         write(*,*)
         CALL GETREC(20,'JOBARC','SCFENEG ',iintfp,acesener)
         if (abs(acesener-(enuc+eelec)).gt.1.0D-6) then
            write(*,48) acesener
            write(*,*)
            write(*,*) "Localized Orbitals May NOT be Non-Canonical UHF"
         end if
      end if

      write(*,*)
      write(*,*) "ACES: Successful NBO-type Localization - Exiting"
      write(*,*)

 16   format(A20, ' Matrix')
 19   format('Electronic Energy                 = ', F18.8)
 20   format('Nuclear Energy                    = ', F18.8)
 21   format('Total Energy                      = ', F18.8)
 31   format('Logical variable putaces = ',L5)
 48   format('HF energy from ACES = ',F20.15)

      CALL aces_fin

      end

