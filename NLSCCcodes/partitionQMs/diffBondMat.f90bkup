subroutine diffBondMat(coords,fullBondMat,lineCount,&
        & bondMat,hcapList,QM1index,QM1count,QM2index)

        integer,intent(in)::lineCount,QM1count
        integer,intent(in)::QM1index(QM1count)
        integer, intent(in)::QM2index(lineCount)
        real, intent(in)::coords(lineCount,3)
        integer,intent(in)::fullBondMat(lineCount,4)
        integer,intent(in)::bondMat(lineCount,4)
        integer,intent(inout)::hcapList(lineCount,4)
        integer::union(lineCount,4)
        integer::i,j,k,scalarA,scalarB(4)
!!      Iterates through each element are visually compare
        do i=1,lineCount
          do j=1,4
!                  print*,"index ",i,"FB:",fullBondMat(i,j),&
!                &       "qmbond",bondMat(i,j)
          ! If an element of the total Bond Matrix, fullBondMat, 
          !             has the bonding
          ! information contained in QM2 bond matrix, bondMat, then
          ! continue to the next iteration of j loop; ie 
          ! Do not store anything in the HcapList
!                  if (any(fullBondMat(i,j) .eq.bondMat(i,:)))cycle

          ! If an element of fullBondMat is *NOT* in the QM bondMat
          !                     ****AND****
          !    the element (or originating,central atom from which we  
          !    define the bond in terms of) is in the QM2 region
          !                     ****THEN***
          !     That element has a severed bond and needs to be put in 
          !     hcapList
           !     if (bondMat(i,1).eq.100000) exit !dont iterate over QM1/2 
                scalarA=fullBondMat(i,j)
                scalarB=bondMat(i,:)
!                print*,'scalar AB',scalarA, scalarB
!                if (.not.any(scalarB.eq.scalarA)) then
!                    if (any(QM2index .eq. i))then
!                        hcapList(i,j)=fullBondMat(i,j)
!                    endif
!                endif
          k=1
          do
                print*,'ijk',i,j,k,bondMat(i,k),fullBondMat(i,j) 
                if (fullBondMat(i,j) .eq. bondMat(i,k)) then
                        !print*,'equal fbm/bm',fullBondMat(i,j),bondMat(i,k)
                        exit 
                endif
               ! print*,'fbm/bm',fullBondMat(i,j),bondMat(i,k)
                k=k+1
                if (k.gt.4)then 
               ! else if ((fullBondMat(i,j).ne.&
               !         &       bondMat(i,k))) then      
                        hcapList(i,j)=fullBondMat(i,j)
                print*,'exiting fbm/bm',hcapList(i,j)
                       exit 
                 endif
          enddo    
 
          enddo
        enddo
        print*,'hcapList:'
        do i=1,lineCount
           do j=1,4
             print*,i,j,bondMat(i,j),fullBondMat(i,j),hcapList(i,j)
           enddo
        enddo
end subroutine
