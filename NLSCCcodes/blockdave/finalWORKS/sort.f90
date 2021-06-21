subroutine sort(wr,vr,Nblocks)
  integer,intent(in)::Nblocks
  double precision,intent(inout)::wr(Nblocks),vr(Nblocks,Nblocks)


  integer::y,z,k,place,resortedArr(Nblocks)
  logical,dimension(Nblocks)::mk
  double precision:: thetaCPY(Nblocks),eVecCPY(Nblocks,Nblocks)
  double precision::total
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
      print*,'resorted Arr:', resortedArr
      thetaCPY=0.0d0
      eVecCPY=0.0d0
      do k=1,Nblocks
        place=resortedArr(k)
        thetaCPY(k)=wr(place)
        eVecCPY(:,k)=vr(:,place)
      enddo
        print*,'thetaCPY',thetaCPY
        wr=thetaCPY!(:Nblocks)
        print*,'wr after thetaCPY',wr
        vr=eVecCPY

! Find norm of vr
!      do y=1,Nblocks
!        total=0.0d0
!        do z=1,Nblocks
!          total=total+vr(z,y)**2
!          print*,'total:', total
!        enddo
!
!        print*,'vr scale factor:', sqrt(total)
!        print*
!        vr(:,y)=vr(:,y)*sqrt(total)
!      enddo



end subroutine
