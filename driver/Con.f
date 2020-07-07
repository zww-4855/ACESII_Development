










      subroutine Con(converge,HFenergy,newHFenergy,Dens,newDens,nbas)
        implicit none
        real(kind=8)::newDens(nbas,nbas),Dens(nbas,nbas)
        integer ::nbas, i,j
        real(kind=8)::HFenergy,newHFenergy,rmsDens
        real(kind=8)::energyTol,densityTol,energyError
        logical :: converge_energy,converge_dens,converge
        energyTol=10D-7
        densityTol=10D-7

        converge_energy=.false.
        converge_dens=.false.
        converge=.false.

        energyError=abs(HFenergy-newHFenergy)
        if (energyError<energyTol) converge_energy=.true.

        rmsDens=0.d0
        do i=1,nbas
          do j=1,nbas
            rmsDens=rmsDens+(Dens(i,j)-newDens(i,j))**2
          end do
        end do
        rmsDens=dsqrt(rmsDens/nbas**2)

        if (rmsDens<densityTol) converge_dens=.true.

        if (converge_energy .and. converge_dens) converge=.true.
      end
