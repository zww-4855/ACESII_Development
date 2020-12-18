      subroutine hfperbond(core,fock,nbas,nocc,energy,
     & energyn,enuc)

C------------------------------------------------------------------
C     HF PER BOND
C     Print the energy of each localized orbitals
C     energy:  electronic energy / total energy
C     energyn: nuclear energy
C------------------------------------------------------------------

      implicit none

      integer i
      
      integer nbas, nocc

      double precision core(nbas,nbas), fock(nbas,nbas),
     & energy(nocc), energyn(nocc), enuc, total, totaln

      total = 0.0D0
      do i = 1, nocc
         energy(i) = core(i,i) + fock(i,i)
         write(*,10) i, energy(i)
         total = total + energy(i)
      end do
      write(*,20) 
      write(*,30) total
      write(*,*)
      
 10   format('Orbital ',I4,' has electronic energy ',F20.15)
 30   format('Electronic energy is ',F25.15)

      totaln = 0.0D0
      do i = 1, nocc
         energyn(i) = (enuc*energy(i))/total
         write(*,40) i, energyn(i)
         totaln = totaln + energyn(i)
      end do
      write(*,20) 
      write(*,50) totaln
      write(*,*)
     
 40   format('Orbital ',I4,' has nuclear energy ',F20.15)
 50   format('Nuclear energy is ',F25.15)

      total = 0.0D0 
      do i = 1, nocc
         energy(i) = energy(i) + energyn(i)
         write(*,60) i, energy(i)
         total = total + energy(i)
      end do
      write(*,20) 
      write(*,70) total
      write(*,*)

 60   format('Orbital ',I4,' has total energy ',F20.15)
 70   format('Total energy is ',F25.15)

 20   format('----------------------------------------')

      return
      end

