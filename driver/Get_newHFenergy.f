










      subroutine Get_newHFenergy(newDens,H,FockMat,newHFenergy,nbas)
        implicit none
        real(kind=8)::newDens(nbas,nbas),H(nbas,nbas)
        real(kind=8)::FockMat(nbas,nbas)
        integer ::nbas, i,j
        real(kind=8)::newHFenergy
        newHFenergy=0.d0
        do i=1,nbas
          do j=1,nbas
            newHFenergy=newHFenergy+newDens(j,i)*(H(i,j)+FockMat(i,j))
          end do
        end do
        !Added in the nuclear-nuclear repulsion
        newHFenergy=0.5d0*newHFenergy+9.13173
        print*, "hf e form fxn: ", newHFenergy
      end
