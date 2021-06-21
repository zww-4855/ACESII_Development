subroutine partitionAlgo(CISmat,NSize,Nblock)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)

  double precision::Hbar(NSize,NSize),P(NSize,NBlock),Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  logical::sVecAlgo,Converged,BWPT,RSPT
  integer::i,j,order,itr,iter
 

  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize,NSize)
  double precision:: tse(NBlock,NBlock),E0(NBlock,NBlock)
  double precision::eps0(NBlock,NBlock),tseOld(NBlock,NBlock)
  double precision::pNew(NSize,NBlock)
  double precision::fullSpace(NSize,100*Nblock)
  
  tseNew=0.0d0
  pNew=0.0d0
  tseOld=0.0d0
  Hbar=CISmat
  BWPT=.TRUE.
  RSPT=.FALSE.
  Converged=.False.
  sVecAlgo=.True.
  order=0 
  itr=1
  call definePQ(CISmat,NSize,1,P,Q,projP,projQ,sVecAlgo)
  iter=1
  do while (.not.Converged)
    print*
    print*,'************************'
    print*,'************************'
    print*,'Beginning iteration:',iter
    if (iter.ne.1) then
      projP=matmul(pNew,transpose(pNew(:,:Nblock)))
      Hbar=matmul(matmul(projP,Hbar),projP)
      call definePQ(Hbar,NSize,1,P,Q,projP,projQ,sVecAlgo)
    endif

    call defineH0(Hbar,NSize,Nblock,P,tse,H0)
    call defineV(Hbar,NSize,V)
    print*,'tse',tse
    E0=tse
    if (iter.eq.1) print*,'Initial eval guess:',tse
    print*,'Going into R0def'
    call defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT)
    call extPSpace(p,R0,V,NSize,NBlock,order,itr,fullSpace,sVecAlgo)
    call buildRedH(Hbar,tse,fullSpace,NSize,Nblock,order,itr,sVecAlgo)
    print*,'Revised eigenvalue:',tse
    print*,'Correction:',tse-E0
!    call buildEps0(p,V,R0,NSize,NBlock,eps0)
!    tseNew=E0+eps0
!    print*,'Correction:',eps0
!    print*,'Revised eigenvalue:',tseNew


    if (abs(abs(tseOld(1,1))-abs(tse(1,1))).gt.1.0D-6) then
!      call redefineP(p,R0,V,pNew,NSize,NBlock)
      tseOld=tse
!      q=p
      pNew(:,1)=fullSpace(:,2)!+fullSpace(:,2)
!      p=pNew
    else
      print*,'************************'
      print*,'Converged on iteration:',iter
      print*,'Converged eigenvalue:',tseNew
      Converged=.true.

    endif

    if (iter.gt.100) then
        print*,'over 100 iterations.... stopping'
        exit
    endif

    iter=iter+1
  enddo
  
end subroutine
