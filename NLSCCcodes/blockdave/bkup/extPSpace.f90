subroutine extPSpace(p,R0,V,NSize,NBlock,order,itr,fullSpace,sVecAlgo)
  integer,intent(in)::NSize,Nblock,order,itr
  double precision,intent(in)::p(NSize,Nblock)
  double precision,intent(in)::V(NSize,NSize),R0(NSize,NSize)
  double precision,intent(inout)::fullSpace(NSize,100*Nblock)
  logical,intent(in)::sVecAlgo

  double precision::R0V(NSize,NSize),tmp(Nsize,NSize)
  double precision::tmpP(NSize,Nblock)
  integer::k

  R0V=matmul(R0,V)
  tmp=R0V
  if (sVecAlgo) then
    
    if (order.eq.0) then
      fullSpace(:,1:NBlock)=p
    endif
! build correction (R0V)^(order+1)|p>

    do k=0,order
      tmp=matmul(tmp,R0V)
    enddo

    tmpP=matmul(tmp,p)
    print*,'correction to p:',tmpP
!  Add correction to full space |R>
    fullSpace(:,(itr)*Nblock+1:(1+itr)*Nblock)=tmpP

! Orthogonalize space: |R>=|p phi0 phi1 ....phip>
! See eq. 
    call GramSchmidt(fullSpace,NSize,(1+itr)*Nblock)
    print*,'orthogonalize correction to p:',fullSpace(:,:(1+itr)*Nblock)
  else  

   print*,'havent coded block scheme yet'

  endif


end subroutine
