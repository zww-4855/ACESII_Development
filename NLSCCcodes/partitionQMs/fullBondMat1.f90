subroutine fullBondMat1(coords,lineCount,fullBondMat)
        integer,intent(in)::lineCount
        real, intent(in)::coords(lineCount,3)
        integer,intent(inout)::fullBondMat(lineCount,4)

        integer::i,j,bondCount,QM2index(lineCount)
        real::bondCut

     bondCut=3.35
     do i=1,lineCount
        QM2index(i)=i
     enddo
     do i=1,lineCount
        bondCount=0
        do j=1,lineCount
          if (i.eq.j) Cycle

          bondCalc=dist(coords,lineCount,QM2index,i,j)
          print*,'bond calc and iter i,j:',QM2index(i),&
        &               QM2index(j),bondCalc

          if (bondCalc.lt.bondCut) then
            bondCount=bondCount+1
            fullBondMat(QM2index(i),bondCount)=QM2index(j)
          endif
          print*,'bond count and mat',bondCount,&
        &       fullBondMat(QM2index(i),bondCount)
          if (bondCount.gt.4) then
            print*, "bondMatrix.f90 has calculated more than 4 bonds, which is&
         &      impossible. Please recheck bondCut in partition.f90"
            call EXIT(1)
          endif

        enddo
     enddo
    do i=1,lineCount
        do j=1,4
            print*,'i,j,bm',i,j, fullBondMat(i,j)
        enddo
        print*
    enddo
end subroutine 
