subroutine fullBondMat1(coords,lineCount,fullBondMat,atomName)
        integer,intent(in)::lineCount
        real, intent(in)::coords(lineCount,3)
     character,intent(inout)::atomName(lineCount)
        integer,intent(inout)::fullBondMat(lineCount,4)
     character(len=1)::site1,site2
        integer::i,j,sCount,test(lineCount)
        real::distance,bondCut
     logical :: checkbondDist
     do i=1,lineCount
        sCount=0
        do j=1,lineCount
          if (i.eq.j) Cycle

          bondCalc=distance(coords,lineCount,i,j)
          print*,'bond calc and iter i,j:',i,&
        &               j,bondCalc

          site1=atomName(i)
          site2=atomName(j)
         if (checkbondDist(site1,site2, bondCalc)) then
            sCount=sCount+1
            fullBondMat(i,sCount)=j
          print*,'bond count and mat',sCount,&
        &       fullBondMat(i,sCount)
          endif
          if (sCount.gt.4) then
            print*, "bondMatrix.f90 has calculated more than 4 bonds, which is&
         &      impossible. Please recheck bondCut in partition.f90"
            call EXIT(1)
          endif

        enddo
     enddo
    print*
    print*,'full bond matrix'
    print*
    do i=1,lineCount
        do j=1,4
            print*,'i,j,bm',i,j, fullBondMat(i,j)
        enddo
        print*
    enddo
end subroutine 
