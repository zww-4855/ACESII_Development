










      subroutine Get_GMatrix(Gmatrix,Dens,G2,nbas)
       implicit none
       integer :: nbas,mu,nu,lambda,sigma
       real(kind=8),intent(inout) :: Gmatrix(nbas,nbas)
       real(kind=8),intent(in) :: Dens(nbas,nbas)
       real(kind=8),intent(in) :: G2(nbas,nbas,nbas,nbas)
       Gmatrix=0.d0
       do mu=1,nbas
         do nu=1,nbas
           do lambda=1,nbas
             do sigma=1,nbas
               !print*, "hi "
       Gmatrix(mu,nu)=Gmatrix(mu,nu)+Dens(lambda,sigma)*
     &          (G2(mu,nu,sigma,lambda)-(0.5d0)*G2(mu,lambda,sigma,nu))
             end do
           end do
         end do
       end do
      end
