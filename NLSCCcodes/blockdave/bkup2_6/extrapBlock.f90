subroutine extrapBlock(revR,oldR,Tvecs,TVals,thetaOld,Memleft, &
                     Nblocks,NSize,itr,Irrepx,Tol,Converged)
  integer, intent(in)::Nblocks,NSize,Irrepx,Memleft,itr
  double precision,intent(inout)::revR(NSize,Nblocks)
  double precision, intent(inout)::oldR(NSize,Nblocks)
  double precision,intent(inout)::Tvecs(Nblocks*itr,Nblocks*itr)
  double precision, intent(in)::TVals(Nblocks),thetaOld(Nblocks)
  double precision, intent(in)::Tol
  logical, intent(inout)::Converged

  integer::counter
  double precision::diff(Nblocks)
  double precision::HbarDiag(Nsize)
  double precision::SCR(Nsize,Nblocks)
  double precision::factor
  double precision::revRSCR(Nsize,Nblocks)
! calculate (Hbar*C - eig*C)*Tvec
  counter=(Nblocks*Nsize*itr)**2
  call getlst(HbarDiag,1,1,1,1,472)
  revRSCR=0.0d0
  do j=1,Nblocks
    !revR(:,j)=matmul(revR,Tvecs(:,j)) !colsV*X[:,j]
    SCR=revR
    call daxpy(counter,-thetaOld(j),oldR,1,SCR,1)
    SCR(:,j)=matmul(SCR,Tvecs(:,j))! r1
    factor=1.0/(thetaOld(j)-HbarDiag(j))
    call dscal(Nsize,factor,SCR(:,j),1)
    revRSCR(:,j)=SCR(:,j)
  enddo

  revR=revRSCR 






  diff=abs(thetaNew-thetaOld)
  print*,'List of Residual Norms for iteration ', i/NUMSOL,': '
  print*,diff
  if (all(diff.lt.TOL)) then
    print*,"** Converged on iteration: ",i/NUMSOL
    print*
    print*
    print*,'** List of converged roots: **'
    print*,thetaNew
    cVALS=thetaNew
    print*,'** in (eV) **'
    print*,thetaNew*27.2114
    cROOTStotal=NUMSOL
    Converged=.true.
   endif

end subroutine
