      Subroutine Block_david_driver(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
     +                              Tol,Maxiter,Nblocks)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Converged

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "flags2.com"
      
      Common /Calcinfo/Nroot(8)

      Parameter(Maxblocks=10)

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
!         ** Create list & put CIS vecs on oldRvec list 497
!         ** AFTER HBARXC::: go on 498
          if (Iter_count.eq.1) then
            Call Updmoi(Nblocks,nsize,Irrepx,497,0,0)
            Call Updmoi(Nblocks,nsize,Irrepx,498,0,0)
            Call Dzero(Work,Nsize)
            do k=1,Nblocks
                 Call Getlst(Work,k,1,1,Irrepx,94)
                 If (Iuhf .Eq. 0) Then
                    Fact = One/Dsqrt(Two)
                    Call Dscal(Irpdpd(Irrepx,9),Fact,Work,1)
                 Endif
		 Call Putlst(Work,k,1,1,Irrepx,497)
             enddo
          else ! dump prev iterations vectors to 498 AFTER QR DECOMP
            call blockQR(Work(Ioff_hc),Work(iend),Nblocks,Nsize,IRREPX)
            ioff=1
            do k=1,Nblocks
                call putlst(work(ioff),k,1,0,Irrepx,497)
                ioff=ioff+Nsize
            enddo
          endif

!! ** BEGIN PRIMARY INNER LOOP TO COMPUTE HBAR*C **
C It is assuemed that we need at least Nsize*Nblocks of memory to
C proceed. Note that this does not count any additional requirements 
C within blokdavison. Memory checks are essential since this is 
c very memory intensive process. 

           Ioff_hc = 1
           Do Iblck = 1, Nblocks 
              I000 = Ione
              Iend = I000 + Nsize*Nblocks
              Memleft = Maxcor - Iend 
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)

C This is start of the loop over No. blocks. First iteration starts 
C with the TDA vectors of (Nblocks). Subsequent interation will start
C use Davidson extrapolated vectors. These vectors, starting from CIS
C stored in 498. (lists 498: C_exp, 497: HC_exp) 
C 
              Call Getlst(Work(Ioff_hc),Iblck,1,0,Irrepx,497)
              Call LANCZOS_DUMP_VEC(Irrepx,Work(Ioff_hc),NSize,
     +                                 490,0,0,443,0,IUHF,.False.)

              Call Hbarxc(Work,Maxcor,Iuhf,1,Irrepx)
              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
     +                      Nsize,.False.)
              Call Putlst(Work(Ioff_hc),Iblck,1,0,Irrepx,498)

C Move this to blockdave and read-one vector at time in a loop.
CSSS              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
CSSS     +                      Nsize,.False.)
              Ioff_hc = Ioff_hc + Nsize
   
           Enddo 
          print*,"outside primary loop!!"
C After this loop, we have Nsize*Nblocks of HC vectors and we can call
C block davidson procedure and return Nsize*Nblocks of extraploted 
C vectors. Starting Ioff_hc=1, we have Nsize*Nblock new HC vectors. 

!! Build tridiagonal matrix T, diagonalize to find estimated
!!eigenvalues/vectors. T has maximal size of
!! (Nblocks*maxitr,Nblocks*maxitr)
         Ioff_hc = 1
         Ioff_c=Ioff_hc+NSize*Nblocks
         Ioff_T=Ioff_c+NSize*Nblocks
         Ioff_Tvec=Ioff_T+(Nblocks*maxiter)**2
         Ioff_Tval=Ioff_Tvec+(Nblocks*Iter_count)**2
         Iend=Ioff_Tval+Nblocks*Iter_count

        call constructT(Work(Ioff_hc),Work(Ioff_c),Nsize,
     + Nblocks,Iter_count,Maxiter,Work(Ioff_T),Work(Ioff_Tvec),
     +                                  Work(Ioff_Tval))
        print*,'evals:'
        print*,Work(Ioff_Tval:Iend-1)
        print*,'checksum of eVecs:'
       call checksum("eVecs:",Work(Ioff_Tvec),Nblocks*Nblocks,s) 
!! Call QR decomposition routine to orthogonalize records 497 w.r.t 498
!! ie orthogonalize Rnew w.r.t. Rold
!! **Add Putlst at end of blockQR to overwrite that in inner loop?
!! **Change k variable -- dependent upon # of householder reflectors --in blockQR to save time???

            call extrapBlock(Work(Ioff_hc),Work(Ioff_c),Work(Ioff_Tvec),Work(Ioff_Tval),thetaOld,Memleft,
     +                Nblocks,NSize,Iter_count,Irrepx,Tol,Converged)
        stop
CSSS       Call Blcokdave(Irepx,Work(Ioff_hc),Work(Iend),Memleft,
CSSS     +                Nbocks,Irrepx,Tol,Converged)

C When returned this call expects extrpolated vectors arrive in
C Work(ioff_hc) where Ioff_hc=1 

      Enddo 
  
      Return 
      End
   
