subroutine partitionAlgo(CISmat,NSize,Nblock)
  integer,intent(in)::NSize,Nblock
  double precision,intent(in)::CISmat(NSize,NSize)

  double precision::Hbar(NSize,NSize),P(NSize,NBlock),Q(Nsize,NSize-NBlock)
  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  logical::sVecAlgo,Converged,BWPT,RSPT
  integer::z,i,j,order,itr
 

  double precision::H0(NSize,NSize),V(NSize,NSize)
  double precision::R0(NSize,NSize)
  double precision:: tse(NBlock,NBlock),E0(NBlock,NBlock)
  double precision::eps0(NBlock,NBlock),tseOld(NBlock,NBlock)
  double precision::norm,pNew(NSize,2)
  double precision::fullSpace(NSize,100*Nblock)
  double precision::tmpProjP(NSize,NSize)!,evecs(NSize,NSize)
  double precision:: RHR(2,2),evecs(2,2),fatPhi(NSize,2)
  double precision::R0V(NSize,NSize),tmp(Nsize,NSize)
  double precision::tmpP(NSize,1),identity(Nsize,Nsize)
  double precision::orthoQ(NSize,NSize+1)
  tseNew=0.0d0
  pNew=0.0d0
  tseOld=0.0d0
  Hbar=CISmat
  BWPT=.TRUE.
  RSPT=.FALSE.
  Converged=.False.
  sVecAlgo=.True.
  order=0 
  call definePQ(CISmat,NSize,1,P,Q,projP,projQ,sVecAlgo)
  itr=1
  identity=0.0d0
  do i=1,NSize
    identity(i,i)=1.0d0
  enddo
  orthoQ=0.0d0
  orthoQ(:,2:NSize+1)=q
  do while (.not.Converged)
    print*
    print*,'************************'
    print*,'************************'
    print*,'Beginning itration:',itr
    if (itr.ne.1) then
!      projP=matmul(pNew,transpose(pNew(:,:Nblock)))
!      projP=matmul(fullSpace(:,:itr),transpose(fullSpace(:,:itr)))
!       projP=matmul(pNew,transpose(pNew(:,:1)))
!      tmpProjP=projP
!      call eig(tmpProjP,evecs,0,NSize,0)
!      print*,'eigvalues of projP:',tmpProjP
!      Hbar=matmul(matmul(projP,Hbar),projP)
!      call definePQ(Hbar,NSize,1,P,Q,projP,projQ,sVecAlgo)

    endif
    print*,'test on outerProj of H0 only'
    print*,matmul(matmul(projP,H0),projP)
    call defineH0(Hbar,NSize,Nblock,P,tse,H0)
    call defineV(Hbar,NSize,V)
    print*,'tse',tse
!    TSE=0.841916685500910
    E0=tse
    if (itr.eq.1) print*,'Initial eval guess:',tse
    print*,'Going into R0def'
    z=1
!    if (itr.eq.1) then
      call defineR0(R0,tse,NSize,Nblock,H0,q,RSPT,BWPT,sVecAlgo,z)
!    else

!    endif
!*****************************
! Build correction
  R0V=matmul(R0,V)!Hbar) !!! **MODIFIED 5/22/2021!V)
  tmp=R0V
  if (order.eq.0 .and. itr.eq.1) then ! change here 4.29.21
    fullSpace(:,1:NBlock)=p
  endif
! build correction (R0V)^(order+1)|p>

  do k=1,order!0,order ! change here 4.29.21
    tmp=matmul(tmp,R0V)
  enddo

    tmpP=matmul(tmp,p)
    print*,'correction to p:',tmpP
!  Add correction to full space |R>
!  if (sVecAlgo) then ! **SINGLE VECTOR ONLY
    fullSpace(:,2)=tmpP(:,1) !(itr)*Nblock+1:(1+itr)*Nblock)=tmpP
    call GramSchmidt(fullSpace,NSize,2)
!    norm=dot_product(fullSpace(:,1),fullSpace(:,1))
!    norm=1.0d0/sqrt(norm)
!    fullSpace=norm*fullSpace
!!! ******DONE


!************************
! Compute <R|H|R>, then diagonalize
    RHR=0.0d0
    RHR=matmul(matmul(transpose(fullSpace(:,:2)),Hbar),fullSpace)
  print*,'RHR',RHR
  call eig(RHR,evecs,0,2,0)
    print*,'new evalues:', RHR

!    call extPSpace(p,R0,V,NSize,NBlock,order,itr,fullSpace,sVecAlgo)
!    call computeRHR(Hbar,tse,fullSpace,NSize,Nblock,order,itr,sVecAlgo,fatPhi)
!    norm=dot_product(fatPhi(:,2),fatPhi(:,2))
!    print*,'norm of fatPhi:',norm
!    norm=1.0d0/sqrt(norm)
!    fatPhi=norm*fatPhi
!    tse=RHR(1,1)
    do i=1,2
      if (RHR(i,i).lt.0.0001) then
        cycle
      endif
      tse=RHR(i,i)
      exit
    enddo
    print*,'Revised eigenvalue:',tse
    print*,'Correction:',tse-E0
    fatPhi=matmul(fullSpace(:,:2),evecs)
    print*,'fatPhi|Hbar|fatPhi'
    print*,matmul(matmul(transpose(fatPhi(:,:2)),Hbar),fatPhi)


    print*,'Cobining 2 columns of fatPhi...'
    fatPhi(:,1)=fatPhi(:,1)!+fatPhi(:,2)
    fatPhi(:,2)=0.0d0

    norm=dot_product(fatPhi(:,1),fatPhi(:,1))
    print*,'norm of fatPhi:',norm
    norm=1.0d0/sqrt(norm)
    fatPhi=norm*fatPhi
    if (abs(abs(tseOld(1,1))-abs(tse(1,1))).gt.1.0D-7) then
      tseOld=tse
      print*,'fat phi:',fatPhi(:,1)
      print*,'fat phi2:',fatPhi(:,2)
      p=0.0d0
      p(:,1)=fatPhi(:,1)
      orthoQ(:,1)=p(:,1)
      call GramSchmidt(orthoQ,NSize,NSize+1)
      q=orthoQ(:,2:NSize+1)

      projP=matmul(fatPhi(:,:1),transpose(fatPhi(:,:1)))
      projQ=identity-projP
      Hbar=matmul(matmul(projP,Hbar),projP)
      print*,'printing new Hbar:'
      call output(Hbar,1,NSize,1,NSize,NSize,NSize,1)
!      pNew(:,1)=fullSpace(:,1)+fullSpace(:,itr+1)
!      p(:,1)=fullSpace(:,itr+1)!p(:,1)+fullSpace(:,itr+1)
!      q=p
!      pNew(:,1)=fullSpace(:,2)!+fullSpace(:,2)
!      p=pNew
!      exit
      print*,'test to first guess'
      print*,matmul(matmul(transpose(p(:,:)),Hbar),p)
    else
      print*,'************************'
      print*,'Converged on itration:',itr
      print*,'Converged eigenvalue:',tse
      Converged=.true.

    endif

    if (itr.gt.10) then
        print*,'over 10 itrations.... stopping'
        exit
    endif

    itr=itr+1
  enddo
  
end subroutine
