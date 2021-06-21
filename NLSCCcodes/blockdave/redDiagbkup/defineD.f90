subroutine defineD(D,initEval,NSize,H0,q,guessE,RSPT,BWPT)
  integer,intent(in)::NSize
  double precision,intent(in)::initEval(1,1),H0(NSize,NSize)
  double precision,intent(in)::q(NSize,NSize-1)
  double precision,intent(inout)::D(NSize-1,NSize-1)
  logical, intent(in):: RSPT,BWPT
  double precision,intent(in)::guessE(1,1)

  double precision::QH0Q(NSize-1,NSize-1)
  double precision:: tmp,tmp2
  integer::i


if (RSPT) then
  QH0Q=0.0d0
  QH0Q=matmul(transpose(q(:,:)),matmul(H0,q))
!  tmp=1.0d0/initEval(1,1)
!  print*,'tmp in D:',tmp
  D=0.0d0
  print*,'QH0Q part of Resolvent:'
  do i=1,NSize-1
    print*,QH0Q(i,i)
    D(i,i)=1.0d0/(initEval(1,1)-QH0Q(i,i))
    print*,'value of D:',D(i,i)
  enddo

endif

if (BWPT) then
  QH0Q=0.0d0
  QH0Q=matmul(transpose(q(:,:)),matmul(H0,q))
  D=0.0d0
  do i=1,NSize-1
    tmp2=guessE(1,1)-QH0Q(i,i)
    if (abs(tmp2).lt.1.0D-2) then
      D(i,i)=sign(1.0D-2,tmp2)
    else
      D(i,i)=1.0d0/(guessE(1,1)-QH0Q(i,i))
    endif
  enddo

endif


end subroutine



