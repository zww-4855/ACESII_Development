      subroutine izerosec(nchoose,bondsize)

      implicit none
  
      integer bondsize, nchoose(bondsize), i

      do i = 1, bondsize
         nchoose(i) = 0
      end do

      return
      end 

