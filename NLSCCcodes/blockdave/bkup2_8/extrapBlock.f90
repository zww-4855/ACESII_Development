subroutine extrapBlock(revR,oldR,Tvecs,thetaOld, &
                     Nblocks,NSize,itr)
  integer, intent(in)::Nblocks,NSize,itr
  double precision,intent(inout)::revR(NSize,itr*Nblocks)
  double precision, intent(inout)::oldR(NSize,itr*Nblocks)
  double precision,intent(inout)::Tvecs(Nblocks*itr,Nblocks*itr)
  double precision,intent(in)::thetaOld(Nblocks*itr)

  integer::counter
  double precision::diff(Nblocks)
  double precision::HbarDiag(Nsize),r1(Nsize)
  double precision::SCR(Nsize,itr*Nblocks)
  double precision::factor
  double precision::revRSCR(Nsize,Nblocks)
! calculate (Hbar*C - eig*C)*Tvec
  print*,'Inside ExtrapBlock'
  counter=Nsize*(Nblocks*itr)
  HbarDiag=0.0d0
  call getlst(HbarDiag,1,1,1,1,472)
call checksum("HbarDiag:",HbarDiag,NSize,s)
  revRSCR=0.0d0
  r1=0.0d0
  do j=1,Nblocks
    SCR=revR
    call daxpy(counter,-thetaOld(j),oldR,1,SCR,1)
    r1=matmul(SCR,Tvecs(:,j))! r1
    factor=1.0d0/(thetaOld(j)-HbarDiag(j))
    print*,"HbarDiag:", HbarDiag(j)
    call dscal(Nsize,factor,r1,1)
    call putlst(r1,Nblocks*itr+j,1,1,1,497)
  enddo


end subroutine
