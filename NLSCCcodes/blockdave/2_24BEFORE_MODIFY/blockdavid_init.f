










      Subroutine Block_david_driver(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
     +                              Tol,Maxiter,Nblocks)

      Implicit Double Precision(A-H,O-Z)
      integer jj
      Dimension Work(Maxcor)
      Logical Converged
      double precision thetaOld(Nblocks),diff(Nblocks)
      double precision temp(NSize)


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



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      
      Common /Calcinfo/Nroot(8)

      Parameter(Maxblocks=20)

C This code assume that the TDA vectors are availble for the 
C requested number of roots via estate_sym (and less <= 10 for
C each symmtery block)
      
      print*,'Entered Block David for EOM'
      print*,Irrepx,Maxcor,Nsize,Iuhf,Iside,
     +             Tol,Maxiter,Nblocks
      print*
     
      !Nblocks = Nroot(Irrepx)
      If (Nblocks .Gt. Maxblocks) Then
         Write(6,"(2a,i2,2a,i2)") " The requested number of roots",
     +                            " for irrep", Irrep, "greater",
     +                            " than the aximum allowed value", 
     +                               Maxblocks
         Call Errex
      Endif

C Pick the TDA eigenvectors for each root per symmetry block and
C carry out the multiplication and accmulate 

      Do Ispin=3,3-2*Iuhf,-1
         Call Zerolist(Work,Maxcor,443+Ispin)
      Enddo 

C This is the start of the block Davidson loop

       Converged = .False.
       Iter_count = 0

       Do while (.not.Converged) 
          Iter_count = Iter_count + 1
          Ioffc = 1
          print*,'********** Begin iteration:, ',Iter_count
!         ** Create list & put CIS vecs on oldRvec list 497
!         ** AFTER HBARXC::: go on 498
          if (Iter_count.eq.1) then
            Call Updmoi(Nblocks*Maxiter,nsize,Irrepx,497,0,0)
            Call Updmoi(Nblocks*Maxiter,nsize,Irrepx,498,0,0)
            Call Dzero(Work,Nsize)
            do k=1,Nblocks
                 Call Getlst(Work,k,1,1,Irrepx,94)
                 If (Iuhf .Eq. 0) Then
                    Fact = One/Dsqrt(Two)
                    Call Dscal(Irpdpd(Irrepx,9),Fact,Work,1)
                 Endif
      !           Work=0.0d0
      !           Work(k)=2.0d0
              call checksum("initVec:",Work,NSize,s)
	       Call Putlst(Work,k,1,1,Irrepx,497)
            enddo
           print*,'outside do'
          !stop
          else ! dump prev iterations vectors to 498 AFTER QR DECOMP
            call blockQR(Work(Ioffc),Iter_count,Nblocks,
     &                                          Nsize,IRREPX)
          endif

!! ** BEGIN PRIMARY INNER LOOP TO COMPUTE HBAR*C **
C It is assuemed that we need at least Nsize*Nblocks of memory to
C proceed. Note that this does not count any additional requirements 
C within blokdavison. Memory checks are essential since this is 
c very memory intensive process. 
           print*,'Nsize is: ', Nsize
           Ioff_hc = 1

           !call printHBAR(Work,Nsize,Maxcor,Iuhf,Irrepx)
           !stop

           Do Iblck = 1, Nblocks 
              I000 = Ione
              Iend = I000 + Nsize*Nblocks
              Memleft = Maxcor - Iend 
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)

              Call Getlst(Work(Ioff_hc),Iblck+Nblocks*(Iter_count-1),
     +                                                  1,0,Irrepx,497)
              Call LANCZOS_DUMP_VEC(Irrepx,Work(Ioff_hc),NSize,
     +                                 490,0,0,443,0,IUHF,.False.)

              Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)

              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
     +                      Nsize,.False.)
              Call Putlst(Work(Ioff_hc),Iblck+Nblocks*(Iter_count-1)
     +                                                  ,1,0,Irrepx,498)

              call checksum("H*C:",Work(Ioff_hc),NSize,s)
              Ioff_hc = Ioff_hc + Nsize
   
           Enddo 
          print*,"outside primary loop!!"

         Ioff_hc = 1
         Ioff_c=Ioff_hc+NSize*Nblocks*Iter_count
         Ioff_T=Ioff_c+NSize*Nblocks*Iter_count
         Ioff_Tvec=Ioff_T+(Nblocks*Iter_count)**2
         Ioff_Tval=Ioff_Tvec+(Nblocks*Iter_count)**2
         Iend=Ioff_Tval+Nblocks*Iter_count
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)



        call constructT(Work(Ioff_hc),Work(Ioff_c),Nsize,
     + Nblocks,Iter_count,Maxiter,Work(Ioff_T),Work(Ioff_Tvec),
     +                                  Work(Ioff_Tval))
        !stop
        print*,'evals:'
        print*,Work(Ioff_Tval:Iend-1)
        print*,'evals (eV):'
        print*,Work(Ioff_Tval:Iend-1)*27.2114
       ! print*,'checksum of eVecs:'
       !call checksum("eVecs:",Work(Ioff_Tvec),Nblocks*Nblocks,s) 
!! Call QR decomposition routine to orthogonalize records 497 w.r.t 498
!! ie orthogonalize Rnew w.r.t. Rold
!! **Add Putlst at end of blockQR to overwrite that in inner loop?
!! **Change k variable -- dependent upon # of householder reflectors --in blockQR to save time???
        print*,'before extrapBlock'
            call extrapBlock(Work(Ioff_hc),Work(Ioff_c),Work(Ioff_Tvec),
     +                Work(Ioff_Tval),Nblocks,NSize,Iter_count)

      if(Iter_count.eq.1) then
        thetaOld=1.0d0
      endif
      diff=abs(thetaOld-Work(Ioff_Tval:Ioff_Tval+Nblocks-1))
      print*,'diff:',diff
      print*,'TOL:',TOL
      print*,all(diff.lt.TOL)

       if (all(diff.lt.TOL)) then
        print*,"** Converged on iteration: ",Iter_count
        print*
        print*
        print*,'** List of converged roots: **'
        print*,thetaOld
        print*,'** in (eV) **'
        print*,thetaOld*27.2114
        Converged=.True.
!! Write final converged eVECS to 497, eVALS to 498
        print*,'** Writing eVecs to disk **'
        print*,'number of blocks:',Nblocks
        call storeFinalVec(Work(Ioff_c),Work(Ioff_Tvec),NSize,
     +                  Nblocks,Iter_count)
        do jj=1,Nblocks
          print*,'jj', jj
!          call getlst(Work,Nblocks*Iter_count+jj,1,1,1,497)
!          call checksum("FINAL EVec:",Work,NSize,s)
!        temp=Work(1:Nsize)
!        print*,'shape:',shape(temp),size(temp)
!        print*,'overlap',dot_product(temp,temp)
!!        print*,'overlap:',matmul(transpose(temp),temp)
!              Call LANCZOS_DUMP_VEC(Irrepx,Work,NSize,
!     +                                 490,0,0,443,0,IUHF,.False.)
!
!              Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)
!              Call Loadvec1(Irrepx,Work,Maxcor,Iuhf,490,2,460,
!     +                      Nsize,.False.)
!          run=0.0d0
!          do aa=1,NSize
!            run=run+temp(aa)*Work(aa)
!          enddo
!          print*,'Subspace Eval:', run !dot_product(transpose(temp),Work) 
!          print*,'final theta:', thetaOld(jj)
!          call putlst(Work,jj,1,1,1,497)
          call putlst(thetaOld(jj),jj,1,1,1,95)
          !call putlst(thetaOld(jj),jj,1,1,1,498)
          print*,'end do'
        !  stop
        enddo 
!        exit
       endif
       thetaOld=Work(Ioff_Tval:Ioff_Tval+Iter_count*Nblocks-1)
       print*,'theta old',thetaOld!*27.2114


!      if (Iter_count.gt.50) then
!        stop
!      endif
      Enddo 

!      print*,'Starting Single Vector David*** '
!      call singledavid_init(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
!     +                         Tol,Maxiter,Nblocks)
!      stop
      Return 
      End
   
