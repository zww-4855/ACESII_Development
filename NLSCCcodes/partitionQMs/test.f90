program test
        character(len=1):: alh 
        integer::z,a(3),b(3)
  integer :: x(4)
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


!  x = findloc([4,9,0,2,-9,9,1], value = 9)
x(1)=1
x(2)=2
x(3)=3
x(4)=4

print*,shape(x(1:3))


end program
