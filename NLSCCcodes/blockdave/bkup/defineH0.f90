subroutine defineH0(CISmat,NSize,NBlock,P,initEval,H0)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(in)::P(NSize,NBlock)
  
  double precision,intent(inout)::H0(NSize,NSize)
  double precision,intent(inout):: initEval(NBlock,NBlock)
  integer::i
  print*,'p'
  print*,p
  initEval=matmul(matmul(transpose(P(:,1:NBlock)),CISmat),P)
  print*,'initial guess to lowest Eval:',initEval




  H0=0.0d0
  do i=1,NSize
    H0(i,i)=CISmat(i,i)
  enddo
print*,'H0'
call output(H0,1,NSize,1,NSize,NSize,NSize,1)
end subroutine
