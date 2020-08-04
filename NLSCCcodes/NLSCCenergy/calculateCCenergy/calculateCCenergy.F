        program calculateCCenergy 
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
        integer::NLMOnum,NLMOnum2,NLMOorigin2,NLMOorigin
        integer:: countinQM1,countinQM2,ierr,QM1num,QM2num,site
        double precision::Enlsccsd
        integer, allocatable :: QM1atoms(:),QM2atoms(:)
        double precision, allocatable :: T1aa(:), T1bb(:)
        double precision, allocatable :: fia_aa(:),fia_bb(:)
        double precision, allocatable :: Waa(:),Wbb(:),
     &                    Wab(:),Waanew(:,:,:,:)
        double precision, allocatable :: T2aa(:),T2bb(:),
     &                    T2ab(:),T2aanew(:,:,:,:)
        double precision, allocatable::NLMO(:,:),Eccsd(:)
#include "sym.com"
#include "icore.com"
#include "istart.com"

        CALL ACES_INIT(icore,i0,icrsiz,iuhf,.true.)
        MAXCOR=icrsiz
        MXCOR=MAXCOR - MOD(MAXCOR,2)
c        allocate(space(MAXCORE))
        CC=.TRUE. ! Global variable for cmpeng.F
        SING1=.TRUE.
        NONHF=.TRUE.

        call getrec(1,'JOBARC','NREALATM',1,natm)
        call getrec(1,'JOBARC','NBASTOT',1,nbas)
        call getrec(1,'JOBARC','NOCCORB',2,scr)
! Applies only to even numbers of electrons
        nocc=scr(1)
        nvirt=nbas-nocc
        allocate(Eccsd(nocc))
        print*, 'number of atoms,#occ,#virt',natm,nocc,nvirt
c        print*,'NT(1) & NT(2)', NT(1),NT(2)
c        print*, 'NF1', NF1
c        print*, 'NF2', NF2
c        print*, iuhf

c        print*,'pop', POP
c        print*
c        print*,'virt',VRT

        ECORR=0.0d0
        ETOT=0.0d0
        ETOT2=0.0d0
        call CMPENG(ICORE(i0),MAXCOR, 43,0,ECORR,ETOT,ETOT2,iuhf,.TRUE.)



        t1Size=nocc*nvirt
        t2aaSize=nocc*(nocc-1)*nvirt*(nvirt-1)/4
        t2abSize=t1Size*t1Size

        NLIST1=0
        NLIST2=44
        indx=i0
!******************************************************************
! * this if statement loads all T1,fia,T2, and W data structures
! according to iuhf=1 (aka UHF) or rhf (aka RHF) reference
!******************************************************************
        if (iuhf.eq.1) then
           FACTOR=1.0D0
           allocate(T1aa(nocc*nvirt),T1bb(nocc*nvirt))
           allocate(Fia_aa(nocc*nvirt),Fia_bb(nocc*nvirt))
           allocate(T2aa(t2aaSize),T2aanew(nvirt,nvirt,nocc,nocc),
     &          T2bb(t2aaSize),T2ab(t2abSize))
           allocate(Waa(t2aaSize),Waanew(nvirt,nvirt,nocc,nocc),
     &          Wbb(t2aaSize),Wab(t2abSize))
           T1aa=0.0d0
           T1bb=0.0d0
           Fia_aa=0.0d0
           Fia_bb=0.0d0
           T2aa=0.0d0
           T2bb=0.0d0
           T2ab=0.0d0
           Waa=0.0d0
           Wbb=0.0d0
           Wab=0.0d0
           I0T(1)=indx ! T1aa
           CALL GETLST(ICORE(i0),1,1,1,1,90)
           indx=indx+NT(1)*IINTFP
           CALL GETLST(T1aa,1,1,1,1,90)
 
           print*, 'compare mem size:', t1Size, NT(1)*IINTFP
           I0T(2)=indx ! T1bb
           CALL GETLST(ICORE(indx),1,1,1,2,90)
           indx=indx+NT(2)*IINTFP
           CALL GETLST(T1bb,1,1,1,2,90)

           I0F(1)=indx ! fia_aa
           CALL GETLST(ICORE(indx),1,1,1,3,93)
           indx=indx+NT(1)*IINTFP
           CALL GETLST(fia_aa,1,1,1,3,93)
           
           I0F(2)=indx ! fia_bb
           CALL GETLST(ICORE(indx),1,1,1,4,93)
           indx=indx+NT(2)*IINTFP
           CALL GETLST(fia_bb,1,1,1,4,93)

!       * Load Waa, Wab, Wab terms
           CALL GETALL(Wab,t2abSize,1,16) ! Wab
           indx=indx+t2absize
           CALL GETALL(Waa,t2aaSize,1,14) ! Waa
           indx=indx+t2aaSize
           CALL GETALL(Wab,t2aaSize,1,15) ! Wbb
           indx=indx+t2aaSize

!       * Load T2aa, T2ab, T2bb terms
           CALL GETALL(T2ab,t2abSize,1,46) ! T2ab
           indx=indx+t2absize
           CALL GETALL(T2aa,t2aaSize,1,44) ! T2aa
           indx=indx+t2aaSize
           CALL GETALL(T2bb,t2aaSize,1,45) ! T2bb
           indx=indx+t2aaSize

        else
          FACTOR=2.0D0
           allocate(T1aa(nocc*nvirt))
           allocate(Fia_aa(nocc*nvirt))
           allocate(T2aa(t2aaSize),
     &          T2ab(t2abSize),T2aanew(nvirt,nvirt,nocc,nocc))
           allocate(Waa(t2aaSize),Waanew(nvirt,nvirt,nocc,nocc),
     &          Wab(t2abSize))
           T1aa=0.0d0
           Fia_aa=0.0d0
           T2aa=0.0d0
           T2ab=0.0d0
           Waa=0.0d0
           Wab=0.0d0
!       * Load T1 & fia terms for RHF
          CALL GETLST(T1aa,1,1,1,1,90)
          CALL GETLST(fia_aa,1,1,1,3,93)
!       * Load Waa, Wab, Wab terms
           CALL GETALL(Wab,t2abSize,1,16) ! Wab
           indx=indx+t2absize
           CALL GETALL(Waa,t2aaSize,1,14) ! Waa
           indx=indx+t2aaSize

!       * Load T2aa, T2ab, T2bb terms
           CALL GETALL(T2ab,t2abSize,1,46) ! T2ab
           indx=indx+t2absize
           CALL GETALL(T2aa,t2aaSize,1,44) ! T2aa
           indx=indx+t2aaSize
c!       * Load Waa, Wab terms for RHF
c           CALL GETALL(ICORE(indx),t2abSize,1,16) ! Wab
c           indx=indx+t2absize
c           CALL GETALL(ICORE(indx),t2aaSize,1,14) ! Waa
c           indx=indx+t2aaSize
c
c!       * Load T2aa, T2ab terms for RHF
c           CALL GETALL(ICORE(indx),t2abSize,1,46) ! T2ab
c           indx=indx+t2absize
c           CALL GETALL(ICORE(indx),t2aaSize,1,44) ! T2aa
c           indx=indx+t2aaSize
        endif

!******************************************************************
! * end of if statment to load in T1,fia, T2, W terms
!******************************************************************

        print*, 't1Size:', t1Size
c        do i=1, t1Size
c         print*,t1aa(i),t1bb(i)
c        enddo
        !print*, 'revised espin', espin
        print*,'ICORE comparison against ICORE and allocatable data 
     & struct for dot prod of T1aa&fia_aa'
        print*, SDOT(t1Size,ICORE(i0), 1,
     &              ICORE(I0F(1)),1)
        print*, SDOT(t1Size, T1aa,1,fia_aa,1)
        print*, "*****************************************"
! Prints the AA, BB, and AB spin contributions to the correlation energy 
! Can be used as a comparison reference to the corresponding values in
! ACESII output file
#ifdef _DEBUG_LVL0

        TAU=.TRUE.
        ECORR=0.0d0
        ETOT=0.0d0
        ETOT2=0.0d0
        ETOTT2=0.0D0
        DO ISPIN=1,IUHF+1
          ESPIN=0.0D0
          LISTT=43+ISPIN
          DO IRREP=1,NIRREP !W/O SYM => 1,1
            print*, IRREP, LISTT
            DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
            NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
            print*, 'dist and numsyt',DISSYT,NUMSYT
            I001=indx !W
            I002=I001+IINTFP*NUMSYT*DISSYT !T2 index
            I003=I002+IINTFP*NUMSYT*DISSYT
            I004=I003+NUMSYT

            
         CALL TENER(NLIST2,ET2,E,NUMSYT,DISSYT,ICORE(I001),
     &             ICORE(I002),ICORE(I0T(ISPIN)),ICORE(I0T(ISPIN)),
     &             ISPIN,TAU,IRREP,POP(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),ICORE(I003))
        print*,'test double prod', DDOT(NUMSYT*DISSYT,T2,1,W,1)       
c            CALL GETLST(ICORE(I001), 1,NUMSYT,2,IRREP
c     &                          ,13+ISPIN)
c            CALL GETLST(ICORE(I002),1,NUMSYT,1,IRREP,
c     &                          ISPIN+NLIST2)
c            ET2=DDOT(NUMSYT*DISSYT,ICORE(I002),1,ICORE(I001),1)
c            print*, 'ET2 is: ', ET2 
          ENDDO
       END DO
#endif

!******************************************************************
! * test on reformatting 1d T1 arrays using subroutine
!******************************************************************
        call expdAA(Waa,Waanew,T2aaSize, nocc,nvirt)
        call expdAA(T2aa,T2aanew,T2aaSize, nocc,nvirt)
        call totalEnergy(T1aa,T1bb,fia_aa,fia_bb,T2aanew,T2ab,T2bb,
     &                   Waanew,Wab,Wbb,nocc,nvirt,iuhf,Eccsd)


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
        print*
        print*,'*****************************************************'
        print*,'** Final NLSCCSD results for NLMO inside QM1 **'
        print*,'*****************************************************'
        print*
       Enlsccsd=0.0d0
        do i=1,nocc
          shared=100 ! means there is nothing there
          countinQM1=0
          countinQM2=0
          do j=1,5
            if (NLMO(i,j).eq.100)then
                 exit
            endif
            if (any(QM1atoms.eq.NLMO(i,j))) then
                countinQM1=countinQM1+1
            else
                countinQM2=countinQM2+1
            endif
          enddo
          if (countinQM2.gt.0) then
                if (countinQM1.gt.0) then
                  print*,"Orbital #", i,"is in QM1&QM2"
                  Enlsccsd=Enlsccsd+0.5d0*Eccsd(i)
                endif
          else
                print*,"Orbital #", i, "is in QM1"
                Enlsccsd=Enlsccsd+Eccsd(i)
          endif
        enddo
        print*
        print*,"NLSCCSD fragmentation energy:", Enlsccsd
30      FORMAT (I4,2X,I4)
c        call testW(Wab,Wab,T2abSize, nocc,nvirt)
c        call testW(T2ab, T2ab, T2abSize,nocc,nvirt)
c        call testF(fia,fia,T1Size,nocc,nvirt)
c        call testF(T1aa,T1aa,T1Size,nocc,nvirt) 

        ! Deallocate memory for T1, fia, T2, W for RHF and UHF cases
        if (iuhf.eq.1) then
          deallocate(T1aa,T1bb,fia_aa,fia_bb,QM1atoms,QM2atoms)
          deallocate(T2aa,T2ab,T2bb,Waa,Wab,Wbb,NLMO,Eccsd)
        else
          deallocate(T1aa,fia_aa,QM1atoms,QM2atoms,Eccsd)
          deallocate(T2aa,T2ab,Waa,Wab,Waanew,T2aanew,NLMO)
        endif
        call aces_fin
        end program
