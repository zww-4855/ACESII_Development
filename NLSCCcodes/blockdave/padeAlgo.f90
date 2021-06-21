subroutine padeAlgo(CISmat,NSize,Nblock,sVecAlgo)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)
  logical,intent(in)::sVecAlgo



  double precision::Hbar(NSize,NSize),P(NSize,NBlock),Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize*NSize*NBlock)
  logical::orderConverge,Converged,BWPT,RSPT
  integer::order
  double precision:: Hxx(Nblock,Nblock),corr(Nblock,Nblock)
  double precision:: HxxNew(NBlock,Nblock),corrOLD(Nblock,Nblock)
  double precision:: approxEval(Nblock,Nblock)
  double precision:: wiOrder(Nblock,Nblock),wiOrderOld(Nblock,Nblock)
  approxEval=0.0d0
  corr=0.0d0
  corrOLD=0.0d0
  BWPT=.TRUE.
  RSPT=.FALSE.
  Converged=.False.

  ! Determine initial approx. to root::  Hxx= <p|H|p>
  call definePQ(CISmat,NSize,Nblock,P,Q,projP,projQ,sVecAlgo) 
!  Hxx=matmul(matmul(transpose(p(:,:)),CISmat),p(:,:))
!  print*,'Initial approx.:', Hxx
  call defineH0PADE(CISmat,NSize,Nblock,P,Hxx,H0)
  print*,'Initial approx.:', Hxx
    

  ! Obtain correction to initial approx. by iterating through subsequent orders:
  ! 
  HxxNew=Hxx
!  do order=1,10
  order=1
  do while (.not.Converged)
!    call defineR0(R0,HxxNew,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,1)
    wiOrderOld=HxxNew
    orderConverge=.False.
!    do while (.not.orderConverge)
      call getPadeCorr(CISmat,HxxNew,NSize,Nblock,sVecAlgo,corr,order,projQ,p)!R0,p)



      HxxNew=Hxx+corr
      wiOrder=HxxNew
!      if (abs(abs(wiOrder(1,1))-abs(wiOrderOld(1,1))).lt.0.00001) then
!        orderConverge=.True.
!      else
!        wiOrderOld=wiOrder
!      endif
!    enddo
    if (abs(abs(corr(1,1))-abs(corrOld(1,1))) .lt.0.000001) then
     print*,'Converged at order:', order
     print*,'Converged eval:', HxxNew
     Converged=.True.
    else 
     corrOld=corr
     order=order+1
     print*,'iteration:',order-1
     print*,'current guess at eval:', HxxNew
    endif

    if (order>20) exit
  enddo
end subroutine
