










      subroutine GetTransform_Xmatrix(X,S,nbas)
        implicit none
        real(kind=8)::X(nbas,nbas),S(nbas,nbas)
        integer ::nbas,i,lwork,ok
        real(kind=8)::work(nbas*(3+nbas/2))

        real(kind=8) :: diagInvS(nbas,nbas), eig(nbas),t(nbas,nbas)
        real(kind=8) :: eigVec(nbas,nbas)
        diagInvS=0.d0
        eig=0.d0
        lwork=nbas*(3+nbas/2)
        !call eig(S,eigVec,100,nbas,1)
        !do i=1,nbas
        !  diagInvS(i,i)=1.0d0/dsqrt(S(i,i))
        !enddo
        !print*, matmul(matmul(diagInvS,S),diagInvS)
        call DSYEV('V','U',nbas,S,nbas,eig,work,lwork,ok)
        if (ok/=0) then
           write(*,*) 'error performing diagonalization'
           stop
        end if
        do i=1,nbas
          diagInvS(i,i)=1.d0/dsqrt(eig(i))
          eigVec(i,i)=eig(i)
        end do
        ! The following verifies s^-1/2 S s^-1/2 =1 
        print*, matmul(matmul(diagInvS,eigVec),diagInvS)
        X=matmul(matmul(S,diagInvS),transpose(S))
      end
