        subroutine blockQR(newRvecs,allRvecs,Nblocks,Nsize,IRREPX)
        integer,intent(in)::Nblocks,Memleft,IRREPX
        double precision,intent(inout)::newRvecs(Nsize,Nblocks)
        double precision,intent(in)::allRvecs(Nsize,2*Nblocks)


        double precision:: tau(2*Nblocks),work(2*Nblocks)
        integer::lwork,info,j,k
! ** VARIABLES: **
! newRvecs - the R vectors recently compiled after hbarxc
! allRvecs - scratch space meant to house the prior and current
! iterations R vectors
! Nblocks - Number of columns per block in our Davidson scheme
! Nsize - row dimension of R vectors
!
! ** PURPOSE **
! The routine loads/stores the old R vectors in 'allRvecs(:,:Nblocks)'
! and places the current iterations R vectors in
! 'allRvecs(:,Nblocks:2*Nblocks)'. Then, allRvecs is sent to BLAS QR
! routine. Upon exit, the new orthogonalized R vectors will be placed at
! allRvecs(:,Nblocks:2*Nblocks); copy this over to 'newRvecs' and exit. 

! Load in old R vectors first.
        k=1
        do j=1,Nblocks
          call getlst(allRvecs(k),j,1,1,IRREPX,498)
          k=k+Nsize
        enddo

! Now copy over new R vectors -> allRvecs
        call dcopy(Nsize*Nblocks,newRvecs,1,allRvecs(k),1)

! Now perform Lapack QR
        lwork=Nblocks*2
        call dgeqrf(Nsize,Nblocks*2,allRvecs,Nsize,tau,work,lwork,info)
        if (info.lt.0) then
          print*,'QR decomposition failed; check blockQR.f90'
          stop
        endif
        call dorgqr(Nsize,Nblocks*2,Nblocks*2,allRvecs,Nsize,tau, &
                        work,lwork,info)
        if (info.lt.0) then
          print*,'QR decomposition failed; check blockQR.f90'
          stop
        endif
! Now overwrite newRvecs with allRvecs(:,Nblocks:Nblocks*2)
! allRvecs(:,Nblocks:Nblocks*2)  --->> newRvecs
        call dcopy(Nsize*Nblocks,allRvecs(k),1,newRvecs,1)

        end subroutine
