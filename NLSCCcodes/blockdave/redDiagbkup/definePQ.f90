subroutine definePQ(CISmat,NSize,P,Q,projP,projQ)
  integer,intent(in)::NSize
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision,intent(inout)::P(NSize,1),Q(Nsize,NSize-1)
  double precision, intent(inout)::projP(NSize,NSize),projQ(NSize,NSize)
  double precision::diagCIS(NSize),tmp,cmp
  integer::i,indx,j
  
  integer::MATDIM,resortedArr(NSize)
  logical,dimension(NSize)::mk

  P=0.0d0
  Q=0.0d0
  projP=0.0d0
  projQ=0.0d0
  resortedArr=0
  mk=.true.
  do i=1,NSize
   diagCIS(i)=CISmat(i,i)
  enddo

  print*,'DiagCIS:',diagCIS
    
  do i=1,NSize
    resortedArr(i)=minloc(diagCIS,1,mk)
    mk(minloc(diagCIS,mk))=.false.
  enddo

  P(resortedArr(1),1)=1.0d0

  print*,'lowest element of diagonal:',resortedArr(1)
  print*,'the rest of sorted arr',resortedArr(2:)

  do i=2,NSize
    Q(resortedArr(i),i-1)=1.0d0
  enddo

  MATDIM=NSize
  print*,'p'

  print*,p
  print*
  print*,'q'
  print*,q
  print*
!  call output(P,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
!  call output(Q,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
  projP=matmul(P,transpose(P))
  projQ=matmul(Q,transpose(Q))
 

  print*,'projP'

  print*
  print*,'projQ' 
!  call output(projP,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
!  call output(projQ,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)




end subroutine
