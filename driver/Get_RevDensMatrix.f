










      subroutine Get_RevDensMatrix(newDens,C,nbas)
        implicit none
        real(kind=8)::newDens(nbas,nbas)
        real(kind=8)::C(nbas,nbas)
        integer ::nbas,nElec
        integer :: i,j,k
        real(kind=8) :: term
      nElec=10
      newDens=0.d0
      do i=1,nbas
        do j=1,nbas
          term = 0.d0
          do k=1,nElec/2
            term=term+C(i,k)*C(j,k)
          end do
         ! term = 2*term
          newDens(i,j)=2*term
        end do
      end do
      end
