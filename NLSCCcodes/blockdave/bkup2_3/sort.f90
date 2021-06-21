subroutine sort(wr,vr,Nblocks)
  integer,intent(in)::Nblocks
  double precision,intent(inout)::wr(Nblocks),vr(Nblocks,Nblocks)


  integer::k,place,resortedArr(Nblocks)
  logical,dimension(Nblocks)::mk
  double precision:: thetaCPY(Nblocks),eVecCPY(Nblocks,Nblocks)
! **Variables::
! wr - eigenvalues of tridiagonal T
! vr - eigenvectors of tridiagonal T
! **Purpose::
! Sort the eigenvalues in increasing order; use this labeling to also
! sort the eigenvectors. Overwrite the data. 

      resortedArr=0
      mk=.true.
      do k=1,Nblocks
        resortedArr(k)=minloc(wr,1,mk)
        mk(minloc(wr,mk))=.false.
      enddo
      thetaCPY=0.0d0
      eVecCPY=0.0d0
      do k=1,Nblocks
        place=resortedArr(k)
        thetaCPY(k)=wr(place)
        eVecCPY(:,k)=vr(:,place)
      enddo
        wr=thetaCPY!(:Nblocks)
        vr=eVecCPY
end subroutine
