subroutine storeFinalVec(oldRvecs,Tvecs,NSize,Nblocks,itr)
  integer, intent(in)::NSize,Nblocks,itr
  double precision,intent(in)::oldRvecs(NSize,Nblocks*itr)
  double precision,intent(in)::Tvecs(itr*Nblocks,itr*Nblocks)

  double precision:: finalVec(NSize,Nblocks*itr),overlap
  integer::i

  finalVec=matmul(oldRvecs,Tvecs)


  do i=1,Nblocks
    overlap=0.0d0
    do j=1,Nsize
      overlap=overlap+finalVec(j,i)*finalVec(j,i)
    enddo
    print*,'finalVec:',i
    print*,'overlap:',overlap
  enddo

  do i=1,Nblocks
    call putlst(finalVec(:,i),i,1,1,1,497)
  enddo





end subroutine
