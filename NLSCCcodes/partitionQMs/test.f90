program test

        integer::a(3),b(3)
        a(1)=2
        a(2)=5
        a(3)=7
        b(1)=0
        b(2)=2
        b(3)=0
        if ( any(a == b) ) print*,"true"
end program
