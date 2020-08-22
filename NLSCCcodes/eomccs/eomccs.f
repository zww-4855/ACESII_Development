










        program eomccs 
        IMPLICIT INTEGER (A-Z)
        integer::indx,iuhf,natm,nbas,nocc,nvirt
        integer*8, dimension(2) :: scr
        integer::t1Size,t2abSize,t2aaSize
c        real, allocatable :: tia(:),fia(:)
        DOUBLE PRECISION tia, fia
        DIMENSION tia(2*1060),fia(2*1060)
        DIMENSION I0T(2),I0F(2)
        DOUBLE PRECISION E,ETOT,FACTOR,ECORR(3),ESPIN,ET2,ETOTT2
     &                          ESING,SDOT
        LOGICAL TAU,NONHF
        LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1

      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
        COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD
        COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
        integer::TOTAL,NLMOnum,NLMOnum2,NLMOorigin2,NLMOorigin
        integer:: offset,countinQM1,countinQM2,ierr,QM1num,QM2num,site
        double precision::Enlsccsd

        CHARACTER*4 input
        integer, allocatable::scrat(:),NLMO(:,:),QM1atoms(:),QM2atoms(:)
        double precision, allocatable::CISmat(:,:),CISevec(:,:)
        double precision, allocatable ::Waanew(:,:,:,:),
     &                                          Waa(:),diaFockA(:)
        double precision, allocatable::Wab(:),Wabnew(:,:,:,:)
        double precision, allocatable:: WaanewOut(:,:,:,:)
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

        CALL ACES_INIT(icore,i0,icrsiz,iuhf,.true.)
        MAXCOR=icrsiz
        MXCOR=MAXCOR - MOD(MAXCOR,2)
c        allocate(space(MAXCORE))
!        CC=.TRUE. ! Global variable for cmpeng.F
!        SING1=.TRUE.
!        NONHF=.TRUE.

        call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NBASTOT',1,nbas)
        call getrec(1,'JOBARC','NOCCORB',2,scr)
! Applies only to even numbers of electrons
        nocc=scr(1)
        nvirt=nbas-nocc
        print*, 'number of atoms,#occ,#virt',natm,nocc,nvirt



        t1Size=nocc*nvirt
        t2abSize=t1Size*t1Size

        NLIST1=0
        NLIST2=44
        indx=i0
        allocate(diaFockA(nocc+nvirt),Waa(t2abSize),CISmat(2*nocc*nvirt
     &          ,2*nocc*nvirt),CISevec(2*nocc*nvirt,2*nocc*nvirt),
     &          Waanew(nvirt,nvirt,nocc,nocc),scrat(nvirt*nvirt),
     &          WaanewOut(nvirt,nocc,nocc,nvirt),
     &          Wab(t1Size**2),Wabnew(nvirt,nvirt,nocc,nocc))
        CALL GETREC(20,'JOBARC','SCFEVALA',(nocc+nvirt)*IINTFP,diaFockA)
        print*,'diafock elements:', diaFockA
        Waa=0.0d0
        Wab=0.0d0
!        CALL GETALL(Waa,t2abSize,1,14)
!        CALL GETALL(Wab,t2abSize,1,16)
!        CALL GETALL(Waa,t2abSize,1,19)
        CALL GETALL(Wab,t2abSize,1,18)
        CALL GETALL(Waa,t2abSize,1,23)
!        CALL GETALL(Wab,t2abSize,1,25)

        print*,'waa/wab',Waa(1),Wab(1)
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
!******************************************************************
! * Calculate NLSCCSD energy per QM1 region
!****************************************************************** 
!        call expdAA(Waa,Waanew,T2aaSize,nocc,nvirt)
!        call expdAA(Wab,Wabnew,T2aaSize,nocc,nvirt)
        total=(occ*virt)**2
        input='2134'
!        call sstgen(Waa,WaanewOut,total,nocc,nvirt,nocc,nvirt,
!     &           scrat,1,input)
!        call sstgen(Wab,Wabnew,total,nvirt,nvirt,nocc,nocc,
!     &           scrat,1,input)

c        do i=1,nocc
c          do a =1,nvirt+1
c            do j=1,nocc
c                do b=1,nvirt+1
c                  print*,Wabnew(a,j,i,b),WaanewOut(a,j,i,b)
c        enddo
c        enddo
c        enddo   
c        enddo   
        print*, '******** 2 e- ints Wab ************'
        call output(Wab,1,nocc*nvirt,1,nocc*nvirt,nocc*nvirt,
     &          nocc*nvirt,1)
        CISmat=0.0d0
        offset=nocc*nvirt
        iter=1
        do i=1,nocc*nvirt
          do a=1,nocc*nvirt
           print*,Wab(iter)
           CISmat(i+offset,a)=Wab(iter)
           CISmat(i,offset+a)=Wab(iter)
           iter=iter+1
          enddo
           print*
        enddo    
        call createCISmat(diaFockA,Waa,Wab,nocc,nvirt,CISmat)
        print*,'*****************************************************'
        print*,'CIS matrix'
        print*,'*****************************************************'
        NUMAA=IRPDPD(1,ISYTYP(1,19))
        NUMBB=IRPDPD(1,ISYTYP(1,20))
        MATDIM=NUMAA+NUMBB
        Write(6,*) "The CIS matrix"
        call output(CISmat,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
c        do i=1,nocc
c          do a=1,nvirt
c             print*, CISmat(i,a),CISmat(i,a+1)
c             print*, CISmat(i+1,a), CISmat(i+1,a+1)
c          enddo
c        enddo
        print*, 'Calling Eigen'
        call eig(CISmat,CISevec,100,2*nocc*nvirt,1)
        print*
        print*,'*****************************************************'
        print*,'** Final NLSCCSD results for NLMO inside QM1 **'
        print*,'*****************************************************'
        print*
        do i=1,2*nocc*nvirt
           print*, CISmat(i,i),CISmat(i,i)*27.2114
        enddo
c       Enlsccsd=0.0d0
c        do i=1,nocc
c          shared=100 ! means there is nothing there
c          countinQM1=0
c          countinQM2=0
c          do j=1,5
c            if (NLMO(i,j).eq.100)then
c                 exit
c            endif
c            if (any(QM1atoms.eq.NLMO(i,j))) then
c                countinQM1=countinQM1+1
c            else
c                countinQM2=countinQM2+1
c            endif
c          enddo
c          if (countinQM2.gt.0) then
c                if (countinQM1.gt.0) then
c                  print*,"Orbital #", i,"is in QM1&QM2"
c                  Enlsccsd=Enlsccsd+0.5d0*Eccsd(i)
c                endif
c          else
c                print*,"Orbital #", i, "is in QM1"
c                Enlsccsd=Enlsccsd+Eccsd(i)
c          endif
c        enddo
c        print*
c        print*,"NLSCCSD fragmentation energy:", Enlsccsd
30      FORMAT (I4,2X,I4)
c        call testW(Wab,Wab,T2abSize, nocc,nvirt)
c        call testW(T2ab, T2ab, T2abSize,nocc,nvirt)
c        call testF(fia,fia,T1Size,nocc,nvirt)
c        call testF(T1aa,T1aa,T1Size,nocc,nvirt) 

        ! Deallocate memory for T1, fia, T2, W for RHF and UHF cases
        deallocate(QM1atoms,QM2atoms,Wab,Wabnew,scrat)
        deallocate(NLMO,Waa,Waanew,Waanewout,diaFockA,CISmat,CISevec)
        call aces_fin
        end program
