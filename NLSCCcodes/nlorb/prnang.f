      subroutine prnang(largel,nbas,mombf,centbf)

C-----------------------------------------------------------
C  Print angular momentum quantum number (l) for each basis
C  function on each atom (s=0, p=1, ...)
C  largel: largest angular quantum number
C-----------------------------------------------------------

      implicit none

      integer i

      integer largel, nbas

      integer mombf(nbas), centbf(nbas)

      largel = 0

      do i = 1, nbas
         write(*,10) i, centbf(i), mombf(i)
         if (largel.lt.mombf(i)) then
            largel = mombf(i)
         end if
      end do 

      write(*,*)
      write(*,*)

 10   format ('Function',I5,' on center',I5,' with angular momentum',I5)

      return
      end

