subroutine diffBondMat(coords,fullBondMat,lineCount,&
        & bondMat,hcapList,QM1index)

        integer,intent(in)::lineCount
        integer,intent(in)::QM1index(4)
        real, intent(in)::coords(lineCount,3)
        integer,intent(in)::fullBondMat(lineCount,4)
        integer,intent(in)::bondMat(lineCount,4)
        integer,intent(inout)::hcapList(lineCount,4)
        integer::union(lineCount,4)
        integer::i,j,k
!!      Iterates through each element are visually compare
        do i=1,lineCount
          if (i.eq.QM1index(i))cycle
          do j=1,4
          print*,"index ",i,"FB:",fullBondMat(i,j),&
                &       "qmbond",bondMat(i,j)

          do k=1,4
                if (fullBondMat(i,j) .eq. bondMat(i,k)) then
                  exit 
                else if ((k.eq.4).and.(fullBondMat(i,j).ne.&
                        &       bondMat(i,k))) then      
                        hcapList(i,j)=fullBondMat(i,j)
                endif
          enddo     
          enddo
        enddo
        print*,'hcapList:'
        do i=1,lineCount
           do j=1,4
             print*,i,j,hcapList(i,j)
           enddo
        enddo
end subroutine
