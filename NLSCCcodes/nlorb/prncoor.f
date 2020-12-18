      subroutine prncoor(natoms,crd,a2u,ach)
      implicit none

C-----------------------------------------------------------------------
C  Print cartesian coordinates with units of Bohr
C  a2u=1.88972599 (transform Angstrom to Bohr)
C  ach: atomic number (H=1, He=2, ...)
C-----------------------------------------------------------------------

      integer natoms, i

      integer ach(natoms)

      double precision crd(3,natoms), a2u

      write(*,10)

      do i = 1, natoms
         write(*,99) i, ach(i), crd(1,i)/a2u, crd(2,i)/a2u,
     & crd(3,i)/a2u
      end do

      write(*,*)
      write(*,*)

 99   format(2I6,3F12.7)
 10   format('Cartesian coordinates')

      return
      end
