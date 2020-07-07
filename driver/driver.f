










      program ringcc

c COMMON BLOCKS


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
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



      integer::ierr,imode,i,j,natom,nbas,mult,ibuf(600)
      real(kind=8),allocatable :: oneh(:),ovrlp(:),H(:,:),S(:,:)
      real(kind=8), allocatable :: G2(:,:,:,:), G3(:,:)
      real(kind=8) :: buf(600),FockSpinned(14,14),summed,denom
      real(kind=8) :: AO_dpx(7,7),AO_dpy(7,7),AO_dpz(7,7)
      real(kind=8) :: a,MO_dpx(10),MO_dpy(10),MO_dpz(10)
      real(kind=8) :: MSO_dpx(40),MSO_dpy(40),MSO_dpz(40)
      real(kind=8), parameter:: energyTol=10D-7
      real(kind=8), parameter :: densityTol=10D-7 
      real(kind=8), allocatable :: X(:,:),energy_den,total,num
      real(kind=8), allocatable :: MOEnergy(:),newDens(:,:),oldDens(:,:)
      real(kind=8), allocatable :: Gmatrix(:,:),Dens(:,:)
      real(kind=8), allocatable :: FockMat(:,:), FockPrime(:,:)
      real(kind=8), allocatable :: C(:,:), PrimeC(:,:),Q_trans(:,:,:,:) 
      real(kind=8) :: HFenergy, newHFenergy,rmsDens, energyError
      logical ::     converge_energy, converge_dens, converge
      integer::scfCount,counter, ierror,aa,bb,lamba
      integer::mu,nu,sigma,e,ff,m,n,ooo
      integer::nop,n_occ,n_virt,nalpha
        integer*8, dimension(2) :: scr
      integer :: counter_i, counter_j
      real(kind=8) :: mp2_energy,summed_y,summed_z
      real(kind=8):: fib,fab
      real(kind=8), allocatable :: MSO(:,:,:,:),Hcore_MO(:,:)
      real(kind=8), allocatable :: Hcore_MSO(:,:),FockSpin(:,:)
      real(kind=8), allocatable :: CIShamil(:,:),B(:,:),CISeigvec(:,:)
      real(kind=8), allocatable :: TDHF(:,:),TDHFeig(:,:)

      COMMON /SYMINF/ NSTART,NIRREP
      DOUBLE PRECISION ONE,TWO,SCFOCC
      INTEGER DIRPRD,JUNK(8),NOCCO(2),NVRTO(2)
      INTEGER NDRPOP(8),NDRVRT(8),POPFUL(8,2),VRTFUL(8,2),NBFIRR(8)
c      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
c      COMMON /FLAGS/  IFLAGS(100)
c      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
c      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
c      COMMON /FILES/  LUOUT,MOINTS
c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      call getrec(20,"JOBARC", "NBASTOT", 1, nbas)
      print*,nbas


      mult=(nbas*(nbas+1))/2
      allocate(oneh(mult), ovrlp(mult),H(nbas,nbas),S(nbas,nbas))
      allocate(G2(nbas,nbas,nbas,nbas), G3(mult,mult))
      allocate(X(nbas,nbas),MOEnergy(nbas),newDens(nbas,nbas))
      allocate(oldDens(nbas,nbas),Gmatrix(nbas,nbas),Dens(nbas,nbas))
      allocate(FockMat(nbas,nbas),FockPrime(nbas,nbas),C(nbas,nbas))
      allocate(PrimeC(nbas,nbas),Q_trans(nbas,nbas,nbas,nbas))
      allocate(MSO(nbas*2,nbas*2,nbas*2,nbas*2),FockSpin(2*nbas,2*nbas))
      allocate(Hcore_MO(nbas,nbas),Hcore_MSO(nbas*2,nbas*2))
      allocate(CIShamil(40,40),CISeigvec(40,40),B(40,40),TDHF(80,80))
      allocate(TDHFeig(80,80))

!     Harvest all 1 and 2 center AO integrals from ACES2
      call Get1EInt(oneh,ovrlp,buf,ibuf,mult)
      call EXPND2(oneh,H,nbas)
      call EXPND2(ovrlp,S,nbas)
      call Get2EInt(G2,G3,buf,ibuf,nbas,mult)

!      Primary SCF LOOP
      call GetTransform_Xmatrix(X,S,nbas)
      HFenergy=0.0d0
      scfCount=0
      Dens=0.0d0
      SCF:do !mu=1,1
        scfCount=scfCount+1
        call Get_GMatrix(Gmatrix,Dens,G2,nbas)
        call Get_FockMatrix(Gmatrix,H,FockMat,nbas)
        call Get_FockPrime(FockPrime,X,FockMat,nbas)
        call GetTransform_Cmatrix(FockPrime,PrimeC,MOEnergy,nbas)
        call Get_CMatrix(C,X,PrimeC,nbas)
        call Get_RevDensMatrix(newDens,C,nbas)
!        write(*,*) 'new density matrix HF:'
!        call output(newDens,1,Nbas,1,Nbas,Nbas,nbas,1)
        a=0.0d0
        do i=1,Nbas
              a=a+newDens(i,i)
        enddo
!        write(*,*) 'trace of density is: ', a
        call Get_newHFenergy(newDens,H,FockMat,newHFenergy,nbas)
        call Con(converge,HFenergy,newHFenergy,Dens,newDens,nbas)
        if (converge) EXIT SCF
        write(*,*) 'HF energy: ', newHFenergy, 'at iteration', scfCount
        Dens=newDens
        HFenergy=newHFenergy
        write(*,*)
        write(*,*)
      end do SCF

!      print*, MOEnergy
      write(*,*) 'Final converged SCF energy is: ', HFenergy
      write(*,*) 'SCF converged in: ',scfCount,' steps'
      write(*,*)
      !Find MP2 energy
c      call Get_QuarterTrans(Q_trans,C,G2,nbas)
c!      print*, Q_trans
c      call Get_MP2energy(Q_trans,FockPrime,nbas,mp2_energy)
c      print*, "MP2 correlation is:  ", mp2_energy
c      print*, "MP2 energy is: " , mp2_energy+HFenergy
c
c      !Transform Hcore from AO spatial -> MO spatial -> MSO
c      !Transform G2 from MO -> MSO
c      !Build Fock Spin matrix
c      !Lastly, build CIS Hamiltonian
c      call Get_MOhcore(H,Hcore_MO,C,nbas)
c      call Get_MSOhcore(Hcore_MO,Hcore_MSO,nbas)
c      call Get_G2MOtoMSO(Q_trans,MSO,nbas)
c      call Get_FockSpinMatrix(FockPrime,FockSpin,nbas)
c      
c      call CIShamiltonian(FockSpin,MSO,CIShamil,HFenergy,nbas)
c      !print*, 'CIS hamil',CIShamil
c      print*
c      print*
c      do i=1,40
c        print*, CIShamil(i,i)
c      enddo
c      call eig(CIShamil,CISeigvec,100,40,1)
c      print*, "CIS solution"
c      print*
c      do i=1,40
c        print*, CIShamil(i,i)
c      enddo
c      print*, "FockPrime(1,1)", FockPrime(1,1)
c      print*, "FockPrime(11,11)", FockPrime(6,6)
c      print*, "FockSpin(1,)", FockSpin(1,1)
c      print*, "2e spinint",MSO(1,11,11,1)!-MSO(1,1,11,11)
c      print*, "fock prime"
c      print*, FockPrime
c      term=0.0d0
c      print*, "SO mbpt2"
c      summed=0.0d0
c
c      call Get_Bmatrix(MSO,B,nbas)
c      call Get_TDHF(TDHF,CIShamil,B,nbas)
c      call eig(TDHF,TDHFeig,100,80,1)
c
c      do i=1,80
c        print*, TDHF(i,i)
c      enddo
c
c      do i=1,5
c        do a=6,7
c          summed=0.0d0
c          summed_y=0.0d0
c          summed_z=0.0d0
c          do mu=1,7
c            do nu=1,7
c              summed=summed+C(mu,i)*C(nu,a)*AO_dpx(mu,nu)
c              summed_y=summed_y+C(mu,i)*C(nu,a)*AO_dpy(mu,nu)
c              summed_z=summed_z+C(mu,i)*C(nu,a)*AO_dpz(mu,nu)
c            enddo
c          enddo
c          MO_dpx(i*(a-5))=summed
c          MO_dpy(i*(a-5))=summed_y
c          MO_dpz(i*(a-5))=summed_z
c         enddo
c      enddo
c      print*, "MO dp in x"
c      print*, MO_dpx 
c      print*, "MO dp in y"
c      print*, MO_dpy
c      print*, "MO dp in z"
c      print*, MO_dpz  
c      counter=1
c      MSO_dpx=0.0d0
c      do i=1,40,4
c        MSO_dpx(i)=MO_dpx(counter)
c        MSO_dpx(i+3)=MO_dpx(counter)
c        counter=counter+1
c      enddo
c      print*,"MSO dp"
c      print*, MSO_dpx

       
      deallocate(oneh,ovrlp,H,S,G2,G3,X,MOEnergy,newDens,Dens)
      deallocate(oldDens, Gmatrix,FockMat,FockPrime,C,PrimeC)
!      deallocate(MSO,Hcore_MO,Hcore_MSO,FockSpin,CIShamil,B,TDHF)
!      deallocate(TDHFeig)
      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

