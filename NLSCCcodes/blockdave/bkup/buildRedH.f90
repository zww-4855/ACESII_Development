subroutine buildRedH(Hbar,tse,fullSpace,NSize,Nblock,order,itr,sVecAlgo)
  integer,intent(in)::NSize,Nblock,order,itr
  double precision,intent(in)::Hbar(NSize,NSize)
  double precision,intent(inout)::fullSpace(NSize,100*Nblock)
  double precision,intent(inout)::tse(NBlock,NBlock)
  logical,intent(in)::sVecAlgo


  double precision:: F(NSize,(1+itr)*NBlock)
  double precision:: Hfs((1+itr)*NBlock,(1+itr)*NBlock)
  double precision:: evecs((1+itr)*NBlock,(1+itr)*NBlock)
  integer::k,i,j
  double precision::HR(NSize,2)

  print*,'incoming H:', Hbar
  Hfs=0.0d0
  evecs=0.0d0
  F(:,:(1+itr)*NBlock)=fullSpace(:,:(1+itr)*NBlock)
  print*,'F is:',F

  HR=matmul(Hbar,F)
  print*,'H|R>', HR

!  do i=1,NSize
!    print*,dot_product(Hbar)
!  enddo
  Hfs=matmul(matmul(transpose(F(:,:)),Hbar),F)
  print*,'printing <R|H|R>',Hfs
! Diagonalize reducedSpace H to obtain new eigenvalues:

  call eig(Hfs,evecs,0,(1+itr)*NBlock,0)

! print new eigenvalues:
  print*,'Revised eigenvalue(s) at order:',order
  print*,HFs
  do k=1,(1+itr)*NBlock!Nblock
    print*,Hfs(k,k)
    if (Hfs(k,k).lt.0.0001d0) then
      cycle
    else
      tse(1,1)=Hfs(k,k)
      exit
    endif
  enddo


end subroutine
