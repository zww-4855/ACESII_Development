program partition
     character (len=15) :: geometry
     character (len=1) :: tmpStr
     character, allocatable :: atomName(:)
     real :: QMcutoff,bondCut 
     real, allocatable::xCoord(:),yCoord(:),zCoord(:),coords(:,:)     
     real, allocatable::QM1coords(:,:)
     integer,allocatable::QM1index(:),QM2index(:)
     integer::lineCount,QM1count
     integer::xmin,xmax,ymin,ymax,zmin,zmax
     integer,allocatable::bondMat(:,:),fullBondMat(:,:)
     integer,allocatable::hcapList(:,:)
     print*,"Please enter the name of the input file containing NLSCC geometry:"
     read*, geometry
     geometry=trim(geometry)
     print*,geometry
     open(unit=159, file=geometry)!,status=unknown)
     if ( ios /= 0 ) stop "Error opening file data.dat"
     lineCount=0
     do
        read(159, '(A)', iostat=ios) 
        if (ios /= 0) exit
        lineCount = lineCount + 1
     end do
     allocate(QM2index(lineCount),atomName(lineCount),&
        &coords(lineCount,3),fullBondMat(lineCount,4),&
        &       bondMat(lineCount,4),hcapList(lineCount,4))

     print*,'lineCount', lineCount   
     rewind(159)
     do i=1,lineCount
        read(159, fmt=*, iostat=ios) tmpStr, coords(i,1),&
                                &coords(i,2), coords(i,3)
        atomName(i)=tmpStr
!        print*,atomName(i),coords(i,1),coords(i,2),coords(i,3)
     enddo
     close(159)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       *** DETERMINE THE # OF ATOMS AND THEIR INDICES
!       *** W.R.T INPUT GEOMETRY FILE
!       *** QM1count-#ofatoms in QM1
!       *** QM1index-list of atom indexes in QM1
!     print*, atomName, coords
     print*, "Please enter the number of atoms you will declare to be in QM1"
     read*, QM1count
     allocate(QM1index(QM1count),QM1coords(QM1count,3))
     do i=1,QM1count
          print*, "Please enter numbering of the atoms you wish to be in QM1"
          read*,QM1index(i)
     enddo 
     print*, QM1index

     do i=1,QM1count
        ii=QM1index(i)
        print*,'ii and QM1 coords', ii, coords(ii,1)
        QM1coords(i,1)=coords(ii,1)
        QM1coords(i,2)=coords(ii,2)
        QM1coords(i,3)=coords(ii,3)
     enddo


     print*,'QM1 coords vs coords', QM1coords
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       User defines the boundary of QM1&QM2

    print*, "Please enter the cutoff value that partitions QM2"
    read*, QMcutoff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Find the boundaries of QM1 (-xmin,xmax,ymin,ymax,zmin,zmax)
!       in order to define the QM1&QM2 but adding the QMcutoff
! function: Returns the index of the coords() array that is least/largest
!           in magnitude in a particular direction
!           ie xmin,xmax,ymin,ymax,zmin,zmax

    call findQM1Bound(QM1coords,QM1index,QM1count,xmin,xmax,ymin,ymax,zmin,zmax) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Now find all atoms within the total QM1/2 boundaries
!       xmin-cutoff < atomXX < xmax+cutoff
!       ymin -cutoff< atomYY < ymax+cutoff
!       zmin-cutoff < atomZZ < zmax+cutoff
! function: Returns QM2index that lists the index of coords() 
!           that is within QM2 region. Double counts QM1 region!!
   
    call findQM2Bound(coords,lineCount,QM1coords,QM1count,QMcutoff,&
        &   QM1index,xmin,xmax,ymin,ymax,zmin,zmax,QM2index)
        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Calculate all distances with respect to the origin
!       of coords(QM2index)

!    call distance(coords,lineCount,QM2index,dist)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Create matrix holding the index of atoms involved in bonding
!       So it should look like the following where index is relevant 
!       site and B.N is the bonding neighbor - at most 4 bonds
!       Ex.
!
!       index#  | B.N.1 | B.N.2 | B.N.3 | B.N.4
!         1     |   3   |   5   | -     | -

     bondCut=2.45664d0
     call bondMatrix(coords,lineCount,QM2index,bondMat,bondCut)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Calculate full bond matrix for the entire system
!       Compare with bondMat() from prior function call
!       to see which atoms have bonds that have been cut
! ****IMPORTANT***
!       HAS A SENSITIVE VARIABLE USED TO DETERMINE A UNIVERSAL 
!       BOND LENGTH; MAY NEED TO PLAY AROUND WITH THIS PARAMETR
     call fullBondMat1(coords,lineCount,fullBondMat)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!       The difference between the fullBondMat and
!       bondMat will tell us which bonds we severed in defining 
!       the QM2 boundaries.

     print*, "qm2 inex", QM2index
     call diffBondMat(coords,fullBondMat,lineCount, bondMat, hcapList,QM1index)
     print*,'shape of QM2index:', shape(QM2index)
     deallocate(QM1index,atomName,bondMat,QM2index,&
        &       fullBondMat,coords,QM1coords,hcapList)!xCoord,yCoord,zCoord)


end program partition
