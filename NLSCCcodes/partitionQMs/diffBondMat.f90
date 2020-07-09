subroutine diffBondMat(coords,fullBondMat,lineCount,&
        & bondMat,hcapList,QM1index,QM1count,QM2index,outsideQMIndex)

        integer,intent(in)::lineCount,QM1count
        integer,intent(in)::QM1index(QM1count),outsideQMIndex(lineCount)
        integer, intent(in)::QM2index(lineCount)
        real, intent(in)::coords(lineCount,3)
        integer,intent(in)::fullBondMat(lineCount,4)
        integer,intent(in)::bondMat(lineCount,4)
        integer,intent(inout)::hcapList(lineCount,4)
        integer::indexOut,union(lineCount,4)
        integer::i,j,k,scalarA,scalarB(4)
!!      Iterates through each element are visually compare
!!      Notation of hcapList is:
!!              hcapList(atomIndex inside QM1/2)=[atomIndex of bonded site
!!                                                      outside QM1/2]
        do i=1,lineCount
          do j=1,4
             
             if (any(bondMat(i,j) .eq. outsideQMIndex)) then
                 !indexOut=findloc(outsideQMIndex,bondMat(i,j))
                 hcapList(i,j)=bondMat(i,j)
             endif     
          enddo    
 
        enddo
        print*,'hcapList:'
        do i=1,lineCount
           do j=1,4
             print*,i,j,bondMat(i,j),fullBondMat(i,j),hcapList(i,j)
           enddo
        enddo
end subroutine
