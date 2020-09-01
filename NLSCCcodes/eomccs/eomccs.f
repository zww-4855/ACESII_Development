




































































































































































































        program eomccs 
        IMPLICIT INTEGER (A-Z)
        integer,allocatable:: NLMOQM1(:),NLMOQM2(:)
        integer::counti,vector(5),indx,iuhf,natm,nbas,nocc,nvirt
        integer::countj
        integer*8, dimension(2) :: scr
        integer::t1Size,t2abSize,t2aaSize
c        real, allocatable :: tia(:),fia(:)
        DOUBLE PRECISION tia, fia
        DIMENSION tia(2*1060),fia(2*1060)
        DIMENSION I0T(2),I0F(2)
        DOUBLE PRECISION E,ETOT,FACTOR,ECORR(3),ESPIN,ET2,ETOTT2
     &                          ESING,SDOT
        LOGICAL TAU,NONHF
      LOGICAL CIS,EOMCC,CISD,FULDIAG,INCORE,READGUES,DOUBLE,NONSTD
      LOGICAL ESPROP,RPA,VPROP,LAMTHERE, NODAVID, TRIPLET
      LOGICAL SS,SD,DS,DD,CC,MBPT2,TRANABCI,CCD,RCCD,DRCCD
      LOGICAL LCCD,LCCSD,CC2,ADC2
      LOGICAL EOM_exite_exist
      DOUBLE PRECISION TDAVID, TMULT, POLTOT, R
!      COMMON / / ICORE(1)
!      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO /  NOCCO(2),NVRTO(2)
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
!      COMMON /SYM/ NT(2),NFMI(2),NFEA(2)!,POP(8,2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /LAMSTATE/ LAMTHERE
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL ,NSIZEC
      COMMON /RESTART_COM/ IRES
      COMMON /POLAR/ POLTOT(3,3)
      COMMON /RMAT/ R(10000)
      COMMON /GUESS/ DOUBLE,NONSTD

      COMMON /METH/CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES
      COMMON /VDINTPRT/NTPERT,NPERT(8),KPERT(8),IDIPSM(3),
     &                 IYZPERT,IXZPERT,IXYPERT,ITRANSX,
     &                 ITRANSY,ITRANSZ,NUCIND
      COMMON/EIGPROB/ISIDE
      COMMON /INTPROG/ VPROP
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON/PARTEOM/NODAVID
      COMMON /SPINSTATE/TRIPLET
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/CALCINFO/NROOT(8)
      COMMON/STATSYM/IRREP
!      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
!     &                DIRPRD(8,8),IRREP0(255,2)
!      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
!        COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
!     &                QCISD
!        COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
        integer::TOTAL,NLMOnum,NLMOnum2,NLMOorigin2,NLMOorigin
        integer:: offset,countinQM1,countinQM2,ierr,QM1num,QM2num,site
        double precision::Enlsccsd
        integer:: NUMDIS
        CHARACTER*4 input
        integer, allocatable::scrat(:),NLMO(:,:),QM1atoms(:),QM2atoms(:)
        double precision, allocatable::CISmatCOPY(:,:)
        double precision, allocatable::CISmat(:,:),CISevec(:,:)
        double precision, allocatable::CISmat0(:,:),CIS0vec(:,:)
        double precision, allocatable ::Waanew(:),
     &                                          Waa(:),diaFockA(:)
        double precision, allocatable::Wbb(:),Wab(:),Wabnew(:,:,:,:)
        double precision, allocatable:: WaanewOut(:)
        double precision, allocatable::wr(:),wi(:),vl(:,:),
     &                                  vr(:,:),work(:)
        double precision :: tempor(6)
        double precision:: revCIS(5,5)
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end


c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


        CALL ACES_INIT(icore,i0,icrsiz,iuhf,.true.)
        CALL SETMET(ICORE(I0),IUHF) ! sets NROOTS line$47 vee.F
        MAXCOR=icrsiz
        MXCOR=MAXCOR - MOD(MAXCOR,2)
c        allocate(space(MAXCORE))
!        CC=.TRUE. ! Global variable for cmpeng.F
!        SING1=.TRUE.
!        NONHF=.TRUE.
        CIS=.TRUE.
        call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NBASTOT',1,nbas)
        call getrec(1,'JOBARC','NOCCORB',2,scr)
! Applies only to even numbers of electrons
        nocc=scr(1)
        nvirt=nbas-nocc
        print*, 'number of atoms,#occ,#virt',natm,nocc,nvirt
        print*, 'nroots', NROOT


        t1Size=nocc*nvirt
        t2abSize=t1Size*t1Size
        t2aaSize=nocc*(nocc-1)*nvirt*(nvirt-1)/4

        NLIST1=0
        NLIST2=44
        indx=i0
        print*,'icore',MAXCOR,IINTFP,iuhf
        call makess(icore(i0),MAXCOR,iuhf)
        call DRVTDA(icore(i0),MAXCOR/IINTFP,iuhf)
        allocate(diaFockA(nocc+nvirt),Waa(t2abSize),CISmat(2*nocc*nvirt
     &          ,2*nocc*nvirt),CISevec(2*nocc*nvirt,2*nocc*nvirt),
     &          Waanew(t2abSize),scrat(nvirt*nvirt),
     &          WaanewOut(t2aaSize),Wbb(t2abSize),
     &          Wab(t1Size**2),Wabnew(nvirt,nvirt,nocc,nocc))

        allocate(NLMOQM1(nbas),NLMOQM2(nbas))
        allocate(CISmat0(nocc*nvirt,nocc*nvirt),
     &                  CIS0vec(nocc*nvirt,nocc*nvirt),
     &           CISmatCOPY(nocc*nvirt,nocc*nvirt))
        CALL GETREC(20,'JOBARC','SCFEVALA',(nocc+nvirt)*IINTFP,diaFockA)
        NLMOQM1=0
        NLMOQM2=0
        print*,'diafock elements:', diaFockA
        Waa=0.0d0
        Wab=0.0d0
        Waanew=0.0d0
        CALL GETALL(WaanewOut,t2aaSize,1,14)
!        CALL GETALL(Wab,t2abSize,1,16)
        CALL GETALL(Waanew,t2abSize,1,19)
        CALL GETALL(Wab,t2abSize,1,18)
        LISTW=23
        NUMDIS=IRPDPD(1,ISYTYP(1,LISTW))
        print*, 'listw&numdis', LISTW,NUMDIS
        NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
        NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
        NUMSYT=IRPDPD(1,ISYTYP(2,LISTW))
        MATDIM=NUMAA+NUMBB
        IOFF=1
        do idis=1,NUMDIS
          call GETLST(Waa(IOFF),idis,1,1,1,LISTW)
          IOFF=IOFF+nocc*nvirt
!        call GETALL(Wbb,t2abSize,1,24)
        enddo
        call GETALL(Wbb,t2abSize,1,24)
!        CALL GETALL(Waa,t2abSize,1,23)
!        CALL GETALL(Waanew,t2aaSize,1,14)

        print*,'waa/wab',WaanewOut(1),Wab(1)
!******************************************************************
! * Read QMcenter
!******************************************************************
        open(unit=250,file='QMcenter',status='old',iostat=ierr)
        if (ierr.eq.0) then
          rewind 250
          read(250,*) QM1num
          read(250,*)
          allocate(QM1atoms(QM1num),QM2atoms(natm-QM1num))
          do i=1,QM1num
            read(250,*) site
            QM1atoms(i)=site
          enddo
          read(250,*)
          read(250,*) QM2num
          do i=1,QM2num
            read(250,*) site
            QM2atoms(i)=site
          enddo
        else
          print*, "error opening QMcenter"
        endif
        close(250)
c        print*, 'qm sites', QM1atoms, QM2atoms
!******************************************************************
! * Read nbocenters
!******************************************************************
        allocate(NLMO(nbas,5))
        open(unit=150,file='nbocenters',status='old',iostat=ierr)
        rewind 150
        read(150,30,IOSTAT=ierr) NLMOnum,NLMOorigin
        NLMO=100
        NLMO(NLMOnum,1)=NLMOorigin
        indx=1
        do while (ierr.eq.0)
           read(150,30,IOSTAT=ierr) NLMOnum2,NLMOorigin2
           if (NLMOnum.eq.NLMOnum2) then
              indx=indx+1
           else
              indx=1
           endif
           NLMO(NLMOnum2,indx)=NLMOorigin2
           NLMOnum=NLMOnum2
           NLMOorigin2=NLMOorigin
        enddo 
        close(150)
        print*, 'printing QMNLMO ownershop'
        do i=1,nbas
          do j=1,5
                if (NLMO(i,j).eq.100) exit
                print*,i,NLMO(i,j)
          enddo
        enddo
        print*
        counti=1
        countj=1
        print*,'QM1/2atoms',QM1atoms
        print*,'QM2atoms',QM2atoms
        print*,any(QM2atoms(1:).eq.NLMO(1,:))
        do i=1,nbas
          terminalIndex=findloc(NLMO(i,1:5),value=100,dim=1)
          print*,'t index', NLMO(i,1:terminalIndex-1)
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
        print*
        print*,'NLMOQM1', NLMOQM1
        print*,'NLMOQM2',NLMOQM2
!******************************************************************
! * Calculate NLSCCSD energy per QM1 region
!****************************************************************** 
        total=(occ*virt)**2
        input='2134'
        print*,'*****************************************************'
        print*, '********        2 e- ints Wab           ************'
        call output(Wab,1,nocc*nvirt,1,nocc*nvirt,nocc*nvirt,
     &          nocc*nvirt,1)
        print*,'*****************************************************'

        print*,'*****************************************************'
        print*, '********        2 e- ints Waa           ************'
        call output(Waa,1,nocc*nvirt,1,nocc*nvirt,nocc*nvirt,
     &          nocc*nvirt,1)
        print*,'*****************************************************'
        print*
        print*,'*****************************************************'
        print*, '********        2 e- ints Wbb           ************'
        call output(Wbb,1,nocc*nvirt,1,nocc*nvirt,nocc*nvirt,
     &          nocc*nvirt,1)
        print*,'*****************************************************'
        print*
        iter=1
        Waanew=0.0d0
        do i=1,nocc
          do j=1,nocc
            do a=1,nvirt
              do b=1,nvirt
               ! if (j<i.and.b<a) then
                        print*,'waanew',Waa(iter)
                        iter=iter+1
                !endif
        enddo
        enddo
        enddo
        enddo 
c        CISmat=0.0d0
c        offset=nocc*nvirt
c        iter=1
c        do i=1,nocc*nvirt
c          do a=1,nocc*nvirt
c           print*,Wab(iter)
c           CISmat(i+offset,a)=Wab(iter)
c           CISmat(i,offset+a)=Wab(iter)
c           iter=iter+1
c          enddo
c           print*
c        enddo    
        call createCISmat(diaFockA,Waa,Wab,nocc,nvirt,CISmat)
        call dcopy((nocc*nvirt)**2,CISmat,1,CISmatCOPY,1)
        print*,'*****************************************************'
        print*,'******               CIS matrix                ******'
        print*,'*****************************************************'
        NUMAA=IRPDPD(1,ISYTYP(1,19))
        NUMBB=IRPDPD(1,ISYTYP(1,20))
        MATDIM=NUMAA+NUMBB
        Write(6,*) "The CIS matrix"
        call output(CISmat,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
        print*,'*****************************************************'
        print*, 'Calling Eigen for FULL CIS matrix'
        print*,'*****************************************************'
!        call eig(CISmat,CISevec,100,2*nocc*nvirt,1)
        print*
        print*,'*****************************************************'
        print*,'** Eigenvalues for full CIS matrix **'
        print*,'*****************************************************'
        print*
!        print*,'     Excite E (au)           Excite E (eV) '
!        do i=1,2*nocc*nvirt
!           print*, CISmat(i,i),CISmat(i,i)*27.2114
!        enddo
        print*
        print*
        print*,'*****************************************************'
        Print*,'*** End of eigenvalues for full CIS matrix ***'
        print*,'*****************************************************'
        print*,'*****************************************************'
        print*
        print*
        print*
        print*,'*****************************************************'
        print*,'******        0th order Approx. NLS-CIS        ******'
        print*,'*****************************************************'
        call nlscisZ(CISmat,CISmat0,nocc,nvirt,QM1atoms,QM1num,
     &                  QM2atoms,QM2num,NLMOQM1,NLMOQM2,nbas)
        call output(CISmat0,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
        print*,'after nlscisZERO call'
        skip=0
        skipC=0
        tempor=(/ (0.0d0,i=1,2*nocc*nvirt) /)
        print*, CISmat0(2,1:2*nocc*nvirt)
        print*
        print*
        print*,tempor
        do i=1,2*nocc*nvirt
           do j=1,2*nocc*nvirt
              print*, CISmat0(i,j)
              revCIS(i-skip,j-skipC)=CISmat0(i,j)
        enddo
        print*
        enddo
        call eig(CISmat,CISevec,100,2*nocc*nvirt,1)
        CALL output(CISmat,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
        CALL output(CISevec,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
!        where (CISmat0(i,i) .ne. 0.0d0) revCIS = CISmat0
!        call output(revCIS,1,MATDIM-1,1,MATDIM-1,MATDIM-1,MATDIM-1,1)
!        o=2*nocc*nvirt
        allocate(wr(o),wi(o),vl(o,o),vr(o,o),work(1))
!        call dgeev('V','V',o,revCIS,o,wr,wi,vl,o,vr,o,work,-1,info)
!        call eig(CISmat0,CIS0vec,100,2*nocc*nvirt,1)
!        print*,'info', info
!        do i=1,2*nocc*nvirt
!           print*, wr(i),wr(i)*27.2114
!        enddo


30      FORMAT (I4,2X,I4)
        deallocate(wr,wi,vl,vr,work)
        ! Deallocate memory for T1, fia, T2, W for RHF and UHF cases
        deallocate(CISmat0,CIS0vec,QM1atoms,QM2atoms,Wab,Wabnew,scrat)
        deallocate(NLMO,Waa,Waanew,Waanewout,diaFockA,CISmat,CISevec)
        deallocate(NLMOQM1,NLMOQM2,Wbb,CISmatCOPY)
        call aces_fin
        end program