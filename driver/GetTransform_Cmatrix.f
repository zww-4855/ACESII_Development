










      subroutine GetTransform_Cmatrix(FockPrime,PrimeC,MOEnergy,nbas)
        implicit none
        real(kind=8)::FockPrime(nbas,nbas),PrimeC(nbas,nbas)
        real(kind=8)::MOEnergy(nbas)
        integer ::nbas, lwork, ok
        real(kind=8):: work(nbas*(3+nbas/2))
        real(kind=8)::vec(nbas,nbas)

        lwork=nbas*(3+nbas/2)
        call eig(FockPrime,vec,100,nbas,1)
        PrimeC=vec 
        print*, "eigenvalue of FockPrime"
        print*, FockPrime      
       ! call DSYEV('V','U',nbas,FockPrime,nbas,MOEnergy,work,lwork,ok)
       ! if (ok/=0) then
       !   write(*,*) 'error performing diagonalization'
       !   stop
       ! end if
       ! PrimeC=FockPrime
      end
