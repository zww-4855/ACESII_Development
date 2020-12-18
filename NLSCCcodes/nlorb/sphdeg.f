      subroutine sphdeg(degen,angmax)

      implicit none
 
      integer angmax, degen(angmax), i

      do i = 1, angmax
         degen(i) = 2*(i - 1) + 1
      end do

      return
      end

