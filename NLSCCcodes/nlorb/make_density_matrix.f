      subroutine make_density_matrix(coef,denmx,ordvc,nbas,nocc,ordered)
      
      implicit none

      integer nbas, nelec, nocc, i, j, k, l,
     & ordvc(nbas)

      double precision coef(nbas,nbas), denmx(nbas,nbas)

      logical ordered

      do i = 1, nbas
         do j = 1, nbas
            do k = 1, nocc 
               l = k
               if (ordered) then
                  l = ordvc(k)
               end if
               denmx(i,j) = denmx(i,j) +
     & coef(i,l)*coef(j,l)
            end do
         end do
      end do

      return
      end

