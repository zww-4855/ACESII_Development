subroutine rdInput(lineCount,coords,atomName)
    integer,intent(inout)::lineCount
     real,allocatable, intent(inout)::coords(:,:)
     character, allocatable,intent(inout) :: atomName(:)
     character (len=15) :: geometryFile
     character (len=15) :: tmpString
     character (len=65), allocatable :: QMlist(:)
     real::QM2boundary
     character (len=1) :: tmpStr
! Number of lines in inputInfo.txt
! Output: Geometry & atomName, as well as QMxx.txt listing
!          the indices of atoms in QM1
     open(unit=159, file='inputInfo.txt')!,status=unknown)
     if ( ios /= 0 ) stop "Error opening file data.dat"
     lineCount=0
     do
        read(159, '(A)', iostat=ios)
        if (ios /= 0) exit
        lineCount = lineCount + 1
     end do

     print*, lineCount
     rewind(159)
!     allocate(QMlist(lineCount-1))

!    Reads Geometry file
     read(159,fmt=*,iostat=ios) geometryFile
     geometryFile=trim(geometryFile)
     print*, geometryFile

!!! Reads Geometry into datastructure
     open(unit=179, file=geometryFile)!,status=unknown)
     if ( ios /= 0 ) stop "Error opening file data.dat"
     lineCount=0
     do
        read(179, '(A)', iostat=ios)
        if (ios /= 0) exit
        lineCount = lineCount + 1
     end do
     allocate(atomName(lineCount),coords(lineCount,3))

     print*,'lineCount', lineCount
     rewind(179)
     do i=1,lineCount
        read(179, fmt=*, iostat=ios) tmpStr, coords(i,1),&
                                &coords(i,2), coords(i,3)
        atomName(i)=tmpStr
!        print*,atomName(i),coords(i,1),coords(i,2),coords(i,3)
     enddo
     close(179)











     read(159,fmt=*,iostat=ios) QM2boundary
     print*,QM2boundary

!   Reads QM1xx.txt files
! https://stackoverflow.com/questions/26920180/fortran-array-of-strings
    do i=1,lineCount-2
         read(159,fmt=*,iostat=ios) tmpString
         QMlist(i)=trim(tmpString)
    enddo

    print*,QMlist


!    deallocate(QMlist)         
end subroutine 
