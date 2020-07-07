      subroutine Get_FockMatrix(Gmatrix,H,FockMat,nbas)
       implicit none
       integer :: nbas,i,j,k,l
       real(kind=8),intent(in) :: Gmatrix(nbas,nbas),H(nbas,nbas)
       real(kind=8),intent(inout) :: FockMat(nbas,nbas)
       do i=1,nbas
         do j=1,nbas
           FockMat(i,j)=H(i,j)+Gmatrix(i,j)
         end do
       end do
      end
