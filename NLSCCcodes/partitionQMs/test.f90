program test
        character(len=1):: alh 
        integer::z,a(3),b(3)
  integer :: x(4)
  integer ::za(8),testarr(8)  ,ix
  integer,allocatable ::zac(:)
      a(1)=2
        a(2)=5
        a(3)=7
        b(1)=0
        b(2)=2
        b(3)=0


        if ( any(a == b) ) print*,"true"
        z=2
        print*,any(a(1) == b)
        print*, any(z == a)
        do i=1,4
            do j=1,10
                
                if (j.eq.6) exit
                print*,i,j
            enddo
        enddo
        print*,'tet',(1==2)
testarr = [0,0,0,2,0,3,0,4]
print*, testarr
print*,'testarr after',pack([(ix,ix=1,size(testarr))],testarr.ne.0)
print*,'shape',size(pack([(ix,ix=1,size(testarr))],testarr.ne.0))
za=0
allocate(zac(size(pack([(ix,ix=1,size(testarr))],testarr.ne.0))))
za=pack([(ix,ix=1,size(testarr))],testarr.ne.0)
zac=pack([(ix,ix=1,size(testarr))],testarr.ne.0)
print*
print*,'za',za
print*
print*,'zac',zac
deallocate(zac)
!   = findloc([4,9,0,2,-9,9,1], value = 9)
x(1)=1
x(2)=2
x(3)=3
x(4)=4

!print*,shape(x(1:3))


end program
