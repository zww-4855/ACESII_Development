subroutine defineH0(CISmat,NSize,P,initEval,H0)
  integer,intent(in)::NSize
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(in)::P(NSize,1)
  
  double precision,intent(inout)::H0(NSize,NSize)
  double precision,intent(inout):: initEval(1,1)
  integer::i
  initEval=matmul(matmul(transpose(P(:,1:1)),CISmat),P)
  print*,'initial guess to lowest Eval:',initEval




  H0=0.0d0
  do i=1,NSize
    H0(i,i)=CISmat(i,i)
  enddo

end subroutine
