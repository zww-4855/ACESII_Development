subroutine bondMatrix(coords,lineCount,QM1index,QM1count,&
                    &           QM2index,bondMat,atomName)
!     Determines the atoms that are bonded in the QM2 region of space *ONLY*
!     An error is thrown if there are more than 4 detectable bonds to a single
!     atom. Returned datastructure is bondMat, having an organization of
!       
!               bondMat( *atomIndex* ) = *list of at most 4 bonding partner*
!
     integer,intent(in)::lineCount,QM1count
     integer,intent(in)::QM2index(lineCount)
     character,intent(inout)::atomName(lineCount)
     real, intent(in)::coords(lineCount,3)
     integer,intent(inout)::bondMat(lineCount,4)
     integer,intent(in)::QM1index(QM1count)
     real::bondCalc,distance
     integer::i,j,absent
     logical :: checkbondDist
     character(len=1)::site1,site2 
     print*,'lincount is: ', lineCount 
     do i=1,lineCount
        bondCount=0
!       If statement ensures we only loop over elements within QM2
        if (.not.any(QM2index .eq. i)) then
                print*,'inside if',i,j,QM2index(i)
                bondMat(i,1)=100000
                cycle 
        endif

          print*,'i.j.',i,j,QM2index(i)
        do j=1,lineCount
          if (i.eq.j) cycle ! Cycle
          bondCalc=distance(coords,lineCount,i,j  )        
          print*,'i,j,dist,char',i,j,bondCalc,&
                &       atomName(i),coords(i,1),bondCalc


!       Determine if the distance is beneath threshhold of typical bond
!       function checkbondDist returns True if calculated distance signals
!       a bond, False otherwise
          site1=atomName(i)
          site2=atomName(j)
          if (checkbondDist(site1,site2, bondCalc)) then
            bondCount=bondCount+1
            bondMat(i,bondCount)=j
          endif

          if (bondCount.gt.4) then
            print*, "bondMatrix.f90 has calculated more than 4 bonds, which is&
         &      impossible. Please recheck bondCut in partition.f90"
            call EXIT(1)
          endif

        enddo
     enddo  
    print*
    print*, 'QM bond matrix'
    print* 

    do i=1,lineCount
        do j=1,4
            print*,'i,j,bm',i,j, bondMat(i,j)
        enddo
        print*
    enddo
end subroutine 
