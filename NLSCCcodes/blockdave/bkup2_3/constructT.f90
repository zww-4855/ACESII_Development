subroutine constructT(newRvecs,Nsize,Nblocks,&
                     itr,maxitr,fullT,eVecs,eVals)
  integer, intent(in)::Nsize,Nblocks,itr,maxitr
  double precision, intent(in)::newRvecs(Nsize,Nblocks)
  double precision,intent(inout)::fullT(Nblocks*maxitr,Nblocks*maxitr)
  double precision,intent(inout)::eVecs(itr*Nblocks,itr*Nblocks)
  double precision,intent(inout)::eVals(Nblocks)

  double precision:: ScrAlpha(Nblocks,Nblocks),ScrBeta(Nblocks,Nblocks)
  double precision::wr(Nblocks),wi(Nblocks),vl(1,Nblocks)
  double precision:: vr(Nblocks,Nblocks),work(4*Nblocks)
  integer::info,aoff,aend,boff
  double precision::prevRvec(Nsize,Nblocks)
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
  if (itr.eq.1) then !only perform alpha_1=newRvecs^t(Hbar*newRvecs)
     call dgemm('T','N',Nblocks,Nblocks,Nsize,1.0d0,newRvecs,Nsize,&
                        newRvecs,Nsize,0.0d0,ScrAlpha,Nblocks)
     fullT(1:Nblocks,1:Nblocks)=ScrAlpha
     call dgeev("N","V",Nblocks,ScrAlpha,Nblocks,wr,wi,vl,1,vr,Nblocks,work,&
                4*Nblocks,info)
     if (info.ne.0) then
        print*,'Something went wrong diagonalizing T, in constructT.f90'
        stop
     endif
     call sort(wr,vr,Nblocks)
  else ! determine alpha_{itr} & beta_{itr}
     aoff=(Nblocks*itr)+1
     aend=aoff+Nblocks
     call dgemm('T','N',Nblocks,Nblocks,Nsize,1.0d0,newRvecs,Nsize,&
                        newRvecs,Nsize,0.0d0,ScrAlpha,Nblocks)

     fullT(aoff:aend,aoff:aend)=ScrAlpha

     boff=aoff-Nblocks !alpha2 has same cols as beta2, and same rows as beta2^t
! Load in previous vectors
     do k=1,Nblocks
       call getlst(prevRvec,k,1,0,1,497)
     enddo

! Calculate overlap
     call dgemm('T','N',Nblocks,Nblocks,Nsize,1.0d0,prevRvecs,Nsize,&
                        newRvecs,Nsize,0.0d0,ScrBeta,Nblocks)
     
     fullT(boff:boff+Nblocks,aoff:aend)=ScrBeta
     fullT(aoff:aend,boff:boff+Nblocks)=transpose(ScrBeta)

     bkupfullT=fullT(1:aend,1:aend)


     call diagFULL(bkupfullT,aend,Nblocks,eVecs,eVals)

  endif

end subroutine
