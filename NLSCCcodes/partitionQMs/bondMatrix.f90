subroutine bondMatrix(coords,lineCount,QM2index,bondMat,bondCut)
     integer,intent(in)::lineCount
     integer,intent(in)::QM2index(lineCount)
     real, intent(in)::coords(lineCount,3),bondCut
     integer,intent(inout)::bondMat(lineCount,4)
     real::bondCalc,dist
     integer::i,j
     do i=1,lineCount
        bondCount=0
        if (QM2index(i).eq.0) Cycle
        ! reference atom
        do j=1,lineCount
          if (QM2index(j).eq.0) Cycle
          if (i.eq.j) Cycle

          bondCalc=dist(coords,lineCount,QM2index,i,j)          
          print*,'bond calc and iter i,j:',QM2index(i),&
        &               QM2index(j),bondCalc

          if (bondCalc.lt.bondCut) then
            bondCount=bondCount+1
            bondMat(QM2index(i),bondCount)=QM2index(j)
          endif

          if (bondCount.gt.4) then
            print*, "bondMatrix.f90 has calculated more than 4 bonds, which is&
         &      impossible. Please recheck bondCut in partition.f90"
            call EXIT(1)
          endif

        enddo
     enddo  
    print*, 'bond matrix:'
    do i=1,lineCount
        do j=1,4
            print*,'i,j,bm',i,j, bondMat(i,j)
        enddo
        print*
    enddo
end subroutine 
