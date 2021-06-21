subroutine defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::tse(NBlock,NBlock),H0(NSize,NSize)
  double precision,intent(in)::q(NSize,NSize-Nblock)
  double precision,intent(inout)::R0(NSize,NSize)
  logical, intent(inout):: RSPT,BWPT

  double precision::Denom(NSize-Nblock,NSize-Nblock)
  integer::MATDIM
MATDIM=Nblock
Denom=0.0d0
if (BWPT) then
  call defineD(Denom,tse,NSize,Nblock,H0,q,RSPT,BWPT)
  R0=matmul(matmul(q,Denom),transpose(q(:,:)))
!print*,R0
print*,'printing Denom'
call output(Denom,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
print*,'printing negative definite R0'
call output(R0,1,NSize,1,NSize,NSize,NSize,1)
else
  print*,'Havent written code for RSPT yet....'
endif
end subroutine
