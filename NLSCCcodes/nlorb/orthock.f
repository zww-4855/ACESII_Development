      subroutine orthock(ovrmo,ovrao,coef,nbas)

      implicit none

      integer nbas, i, j, mu, nu
      double precision ovrmo(nbas,nbas), ovrao(nbas,nbas),
     &                 coef(nbas,nbas)

      call printm("AO matrix",nbas,nbas,nbas,nbas,ovrao)

      do i = 1, nbas
         do j = 1, nbas
            do mu = 1, nbas
               do nu = 1, nbas
                  ovrmo(i,j) = ovrmo(i,j) +
     & coef(mu,i)*coef(nu,j)*ovrao(mu,nu)
               end do
	    end do
         end do
      end do

      call printm("MO matrix",nbas,nbas,nbas,nbas,ovrmo)

      return
      end

