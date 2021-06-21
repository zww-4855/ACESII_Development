subroutine extrapBlock(revR,Tvecs,TVals,thetaOld,Memleft, &
                     Nblocks,NSize,itr,Irrepx,Tol,Converged)
  integer, intent(in)::Nblocks,NSize,Irrepx,Memleft,itr
  double precision,intent(inout)::revR(NSize,Nblocks)
  double precision,intent(inout)::Tvecs(Nblocks*itr,Nblocks*itr)
  double precision, intent(in)::TVals(Nblocks),thetaOld(Nblocks)
  double precision, intent(in)::Tol
  logical, intent(inout)::Converged

  double precision::diff(Nblocks)

  do j=1,Nblocks
    revR(:,j)=matmul(revR,Tvecs(:,j))






  enddo

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
