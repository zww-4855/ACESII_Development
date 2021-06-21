subroutine constructT(newRvecs,oldRvecs,Nsize,Nblocks,&
                     itr,maxitr,fullT,eVecs,eVals)
  integer, intent(in)::Nsize,Nblocks,itr,maxitr
  double precision, intent(inout)::newRvecs(Nsize*Nblocks*itr)
  double precision, intent(inout)::oldRvecs(Nsize*Nblocks*itr)
  double precision,intent(inout)::fullT(Nblocks*itr,Nblocks*itr)
  double precision,intent(inout)::eVecs(itr*Nblocks,itr*Nblocks)
  double precision,intent(inout)::eVals(Nblocks*itr)

  double precision:: ScrAlpha(Nblocks,Nblocks),ScrBeta(Nblocks,Nblocks)
  double precision::wr(Nblocks*itr),tempwr(Nblocks)
  double precision:: wi(Nblocks*itr),vl(1,itr*Nblocks)
  double precision:: vr(Nblocks*itr,itr*Nblocks),tempvr(Nblocks,Nblocks)
  double precision::work(4*Nblocks*itr)
  integer::i,info,aoff,aend,boff,ioff,const
  double precision::prevRvec(Nsize,Nblocks),total
  double precision::bkupFullT(Nblocks*itr,Nblocks*itr)
! **Variables:
! newRvecs - R vectors fresh out of hbarxc
! Nsize    - row dim. of newRvecs
! Nblocks  - col. dim. of newRvecs
! itr      - current iteration of the block algo
! maxitr   - maximum number of iterations in block algo
! fullT    - full (albeit mostly empty) space reserved for depositing relevant
!            elements of T. 
! ScrAlpha -  temporary scratch space to store alpha
! ScrBeta  - temporary scratch space to store beta
! eVecs    - Eigenvectors of matrix T
! eVals    - Eigenvalues of matrix T
!
! **Purpose:
! Performs the block vector-block vector multiplication to create the leements
! of T, namely alpha_i=Rnew^t(Hbar*Rnew),beta_i+1=Rold^t(Hbar*Rnew)
! beta_i+1^t=Rnew^t(Hbar*Rold). Each block is of size (Nblocks x Nblocks)
! Uses current itr to determine placement; must create tridiagonal T. 
! Then Diagonalizes this mini matrix and returns eVecs/values.

  ScrAlpha=0.0d0
  ScrBeta=0.0d0
  ioff=1
  print*,'inside constructT'
  eVals=0.0d0
  eVecs=0.0d0
  do i=1,Nblocks*itr
    call getlst(newRvecs(ioff),i,1,1,1,498)
    call getlst(oldRvecs(ioff),i,1,1,1,497)
    ioff=ioff+Nsize
  enddo
  fullT=0.0d0
!  call checksum("newRvecs      :",newRvecs,NSize*Nblocks,s)
!  call checksum("fullT zerod:",fullT,(maxiter*Nblocks)**2,s)
  const=Nblocks*itr
  call dgemm('T','N',const,const,Nsize,1.0d0,oldRvecs,Nsize,&
                     newRvecs,Nsize,0.0d0,fullT,const)
  print*,'maxiter:',maxitr
  print*,'Nblocks:',Nblocks
  print*,'size of eVals:',size(eVals)
  print*,'size of fullT:',size(fullT)
  print*,'shape of fullT:',shape(fullT)
     call checksum("fullT stored check:",fullT,(const*const),s)
     print*,'T matrix:'
!     call output(fullT,1,const,1,const,const,const,1)
!     call output(ScrAlpha,1,Nblocks,1,Nblocks,Nblocks,Nblocks,1)

     call dgeev("N","V",const,fullT,const,eVals,wi,vl,1,eVecs,const,work,&
                4*const,info)
     print*,'shape of Tvecs:',shape(eVecs),size(eVecs)
     do i=1,const
       call checksum("Tvecsb4sort:",eVecs(:,i),Nblocks,s)
     enddo
     if (info.ne.0) then
        print*,'Something went wrong diagonalizing T, in constructT.f90'
        stop
     endif
     print*,'eigenvals in T:', eVals!wr
     !call checksum("eVecs b4 sort:",eVecs,Nblocks*Nblocks,s)
     call sort(Evals,eVecs,itr*Nblocks)
     do i=1,const
       call checksum("TvecsAFTERsort:",eVecs(:,i),itr*Nblocks,s)
     enddo
     print*,'sorted evals in T:', eVals
     print*, 'sorted Tvecs:'
!     call output(eVecs,1,const,1,const,const,const,1)
     total=0.0d0
     do i=1,const
        total=total+fullT(i,i)
     enddo
     print*,'trace of fullT:', total
     print*,'sum of eigenvalues:',sum(Evals)
end subroutine
