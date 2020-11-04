#include "flags.h"
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
        double precision, allocatable::CISreduce(:,:),CISreduceVec(:,:)
!        double precision ::CISreduce(10,10),CISreduceVec(10,10)
        integer::QM2NLMOcount
        double precision :: temp,tempor(6)
        double precision:: revCIS(5,5)
        logical :: incQM2,CTflag
        integer:: QMregIA(2), compareIA(2)
        double precision,allocatable:: intermedAC(:),acesCISevecs(:),
     &                  acesCISevecsbkup(:)
#include "sym.com"
#include "icore.com"
#include "istart.com"
#include "parcc.com"

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
        print*,'icore',NONHF,MAXCOR,IINTFP,iuhf
        call makess(icore(i0),MAXCOR,iuhf)
        call DRVTDA(icore(i0),MAXCOR/IINTFP,iuhf)
        allocate(diaFockA(nocc+nvirt),Waa(t2abSize),CISmat(2*nocc*nvirt
     &          ,2*nocc*nvirt),CISevec(2*nocc*nvirt,2*nocc*nvirt),
     &          Waanew(t2abSize),scrat(nvirt*nvirt),
     &          WaanewOut(t2aaSize),Wbb(t2abSize),
     &          Wab(t1Size**2),Wabnew(nvirt,nvirt,nocc,nocc),
     &          acesCISevecs(NROOT(1)*2*nocc*nvirt),
     &          acesCISevecsBKUP(NROOT(1)*2*nocc*nvirt),
     &          intermedAC(NROOT(1)*2*nocc*nvirt)) ! CAN HOLD 10 DAVIDSON
                                               !ROOTS

        allocate(NLMOQM1(nbas),NLMOQM2(nbas))
        allocate(CISmat0(2*nocc*nvirt,2*nocc*nvirt),
     &                  CIS0vec(2*nocc*nvirt,2*nocc*nvirt),
     &           CISmatCOPY(2*nocc*nvirt,2*nocc*nvirt))
        CALL GETREC(20,'JOBARC','SCFEVALA',(nocc+nvirt)*IINTFP,diaFockA)
        NLMOQM1=0
        NLMOQM2=0
        Waa=0.0d0
        Wab=0.0d0
        Waanew=0.0d0
        CALL GETALL(WaanewOut,t2aaSize,1,14)
        CALL GETALL(Waanew,t2abSize,1,19)
        CALL GETALL(Wab,t2abSize,1,18)
        LISTW=23
        NUMDIS=IRPDPD(1,ISYTYP(1,LISTW))
        NUMAA=IRPDPD(IRREP,ISYTYP(1,19))
        NUMBB=IRPDPD(IRREP,ISYTYP(1,20))
        NUMSYT=IRPDPD(1,ISYTYP(2,LISTW))
        MATDIM=NUMAA+NUMBB
        IOFF=1
        do idis=1,NUMDIS
          call GETLST(Waa(IOFF),idis,1,1,1,LISTW)
          IOFF=IOFF+nocc*nvirt
        enddo
        call GETALL(Wbb,t2abSize,1,24)
        acesCISevecsBKUP=0.0d0 
        acesCISevecs=0.0d0
!        Nirrep=1
!        Do Ispin=1,iuhf+1
!          Numdis=Irpdpd(Nirrep,8+Ispin)
!          call Getlst(acesCISevecs,1,Numdis,1,Nirrep,94)
!        enddo
        call Getlst(acesCISevecs,1,NROOT(1),1,1,94)
        CALL DCOPY(NROOT(1)*2*nocc*nvirt,acesCISevecs,1,acesCISevecsBKUP
     &                                  ,1)
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
        close(150)
        counti=1
        countj=1
!        print*,'QM1/2atoms',QM1atoms
!        print*,'QM2atoms',QM2atoms
!******************************************************************
! * Determine if the NLMO read from nbocenters is entirely within
! * QM1, QM2, or both QM1&QM2
!******************************************************************
        do i=1,nbas
          terminalIndex=findloc(NLMO(i,1:5),value=100,dim=1)
!          print*,'t index', NLMO(i,1:terminalIndex-1)
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
        call createCISmat(diaFockA,Waa,Wab,nocc,nvirt,CISmat)
        print*,'*****************************************************'
        print*,'******               CIS matrix                ******'
        print*,'*****************************************************'
        NUMAA=IRPDPD(1,ISYTYP(1,19))
        NUMBB=IRPDPD(1,ISYTYP(1,20))
        MATDIM=NUMAA+NUMBB
!#ifdef _DEBUG_LVL0
        print*,'*****************************************************'
        print*,'******           FullCIS matrix                ******'
        print*,'*****************************************************'
        call output(CISmat,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
        print*,'*****************************************************'
        print*,'******           FullCIS eigenvectors          ******'
        print*,'*****************************************************'
        call output(acesCISevecs,1,MATDIM,1,6,6,6,1)
        ! Re-Calculate the ACESII roots from Davidson 
        intermedAC=0.0d0
        iter=1
        do i=1,NROOT(1)!4!6
           call xgemm('N','N',2*nocc*nvirt,1,2*nocc*nvirt,1.0D0,CISmat,
     &          2*nocc*nvirt,acesCISevecs(iter),2*nocc*nvirt,0.0D0,
     &              intermedAC(iter),2*nocc*nvirt)
!           call output(acesCISevecs,1,2*nocc*nvirt,1,2*nocc*nvirt,
!     &          2*nocc*nvirt,2*nocc*nvirt,1)
           temp=sdot(2*nocc*nvirt,acesCISevecs(iter),1,intermedAC(iter)
     &                                  ,1) 
           print*,'recalculating ACESII vars',temp*27.2114
           iter=iter+2*nocc*nvirt
        enddo
!
!
!
!       Now omit relevant elements of the coefficient vectors

        evecSize=NROOT(1)*2*nocc*nvirt !multiplied by number of roots
        roots=NROOT(1)
        call nlsMultiply(acesCISevecs,evecSize,CISmat,nocc,nvirt
     &      ,nbas,roots,NLMOQM1,size(NLMOQM1),NLMOQM2,size(NLMOQM2))
!        iter=1
!        next=1
!        offset=nocc*nvirt
!        do roots=1,4!6!10
!         do i=1,nocc
!          do a=nocc+1,nbas
!            compareIA=0
!            QMregIA=0
!            compareIA=(/ i,a /)
!!            print*, compareIA
!!            print*, 'NLMO manip', NLMOQM1
!            call findQMregion(compareIA,size(compareIA),NLMOQM1,
!     &           size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIA)
!!            print*,QMregIA
!!            print*
!            if (QMregIA(1).eq.2) then
!              acesCISevecs(iter)=0.0d0
!              acesCISevecs(iter+offset)=0.0d0
!            endif
!!           print*, acesCISevecs(iter)
!          iter=iter+1
!          enddo
!        enddo
!           call xgemm('N','N',2*nocc*nvirt,1,2*nocc*nvirt,1.0D0,CISmat,
!     &          2*nocc*nvirt,acesCISevecs(next),2*nocc*nvirt,0.0D0,
!     &              intermedAC(next),2*nocc*nvirt)
!           temp=sdot(2*nocc*nvirt,acesCISevecs(next),1,intermedAC(next)
!     &                                  ,1)
!        iter=iter+nocc*nvirt
!        next=next+2*nocc*nvirt
!           print*,'Ajith suggested root', temp*27.2114
!        enddo
        print*,'*****************************************************'
        print*, 'Calling Eigen for FULL CIS matrix'
        print*,'*****************************************************'
!        call eig(CISmat,CISevec,100,2*nocc*nvirt,1)
!        print*
!        print*,'*****************************************************'
!        print*,'** Eigenvalues for full CIS matrix **'
!        print*,'*****************************************************'
!        print*
!        print*,'     Excite E (au)           Excite E (eV) '
!        do i=1,2*nocc*nvirt
!           print*, CISmat(i,i),CISmat(i,i)*27.2114
!        enddo
!        print*
        print*
        print*,'*****************************************************'
        Print*,'*** End of eigenvalues for full CIS matrix ***'
        print*,'*****************************************************'
        print*,'*****************************************************'
!#endif
        print*
        print*
        print*
        print*,'*****************************************************'
        print*,'******        0th order Approx. NLS-CIS        ******'
        print*,'*****************************************************'
        incQM2=.False.
        CTflag=.False.!.True.
        print*,'NLMOQM1', NLMOQM1
        call nlscisZ(CISmat,CISmat0,nocc,nvirt,QM1atoms,QM1num,
     &          QM2atoms,QM2num,NLMOQM1,NLMOQM2,nbas,incQM2,CTflag)
!
!
! ***** Call second method that zeros out elements of A, then multiplies
!                       C^t A C = E
        call nlsAmodify(CISmat0,nocc,nvirt,acesCISevecsBKUP,evecSize,
     &                          NROOT(1))
#ifdef _DEBUG_LVL0
        print*, "** Reformatted NLS-CIS matrix **"
        call output(CISmat0,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
#endif
        call ReduceCISmat(CISmat0,2*nocc*nvirt,
     &                    NLMOQM1,NLMOQM2,nbas,nocc,nvirt)

        print*
        print*
        print*,'*****************************************************'
        print*,'** END NLS-CIS **'
        print*,'*****************************************************'
30      FORMAT (I4,2X,I4)
        ! Deallocate memory for T1, fia, T2, W for RHF and UHF cases
        deallocate(CISmat0,CIS0vec,QM1atoms,QM2atoms,Wab,Wabnew,scrat)
        deallocate(NLMO,Waa,Waanew,Waanewout,diaFockA,CISmat,CISevec)
        deallocate(acesCISevecs,intermedAC,NLMOQM1,NLMOQM2,Wbb,
     &             acesCISevecsBKUP, CISmatCOPY)

        call aces_fin
        end program
