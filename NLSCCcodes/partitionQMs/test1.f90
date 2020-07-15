program test 
     character (len=15) :: geometry
     character (len=1) :: tmpStr
     character, allocatable :: atomName(:)
     real :: QMcutoff,bondCut 
     real, allocatable::xCoord(:),yCoord(:),zCoord(:),coords(:,:)     
     real, allocatable::QM1coords(:,:)
     integer,allocatable::QM1index(:),QM2index(:)
     integer::lineCount,QM1count
     integer::i,j,xmin,xmax,ymin,ymax,zmin,zmax
     integer,allocatable::bondMat(:,:),fullBondMat(:,:)
     integer,allocatable::hcapList(:,:),outsideQMIndex(:)
     real,allocatable:: revQMcoords(:,:)
     character,allocatable::revAtomName(:)
     real::scaleH
     character (len=15) :: geometryFile
     character (len=15) :: tmpString
     character (len=65), allocatable :: QMlist(:)
     real::QM2boundary

     open(unit=159, file='inputInfo.txt')!,status=unknown)
     if ( ios /= 0 ) stop "Error opening file inputInfo.dat"
     lineCount=0
     do
        read(159, '(A)', iostat=ios)
        if (ios /= 0) exit
        lineCount = lineCount + 1
     end do

     print*, lineCount
     rewind(159)
     allocate(QMlist(lineCount-2))

!    Reads Geometry file
     read(159,fmt=*,iostat=ios) geometryFile
     geometryFile=trim(geometryFile)
     print*, geometryFile

     read(159,fmt=*,iostat=ios) QM2boundary
     print*,QM2boundary

!   Reads QM1xx.txt files
! https://stackoverflow.com/questions/26920180/fortran-array-of-strings
    do i=1,lineCount-2
         read(159,fmt=*,iostat=ios) tmpString
         print*,trim(tmpString)
    enddo

    print*,QMlist
    deallocate(QMlist)
    close(159)
end program
