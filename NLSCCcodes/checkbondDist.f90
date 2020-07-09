function checkbondDist(atomStr1,atomStr2,bondCalc)
        character (len=*),intent(in) :: atomStr1, atomStr2
        real, intent(in):: bondCalc
        logical ::checkbondDist
!
!      
        print*,'str 1 and str2 are: ', atomStr1, atomStr2
        print* 
        if ( trim(atomStr1) == 'C' ) then
                if (trim(atomStr1) == 'C' ) then
                        checkbondDist=bondCalc.lt.1.6d0
                else if (trim(atomStr1) == 'H' ) then 
                        checkbondDist=bondCalc.lt.1.2d0
                else if (trim(atomStr1) == 'O' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else if (trim(atomStr1) == 'N' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else
                print*,"Atom distance not included in checkbondDist.f90"
                stop
                endif


        else if ( trim(atomStr1) == 'O' ) then
                if (trim(atomStr1) == 'C' ) then
                        checkbondDist=bondCalc.lt.1.6d0
                else if (trim(atomStr1) == 'H' ) then
                        checkbondDist=bondCalc.lt.1.2d0
                else if (trim(atomStr1) == 'O' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else if (trim(atomStr1) == 'N' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else
                print*,"Atom distance not included in checkbondDist.f90"
                stop
                endif

        else if ( trim(atomStr1) == 'H' ) then
                if (trim(atomStr1) == 'C' ) then
                        checkbondDist=bondCalc.lt.1.6d0
                else if (trim(atomStr1) == 'H' ) then
                        checkbondDist=bondCalc.lt.1.2d0
                else if (trim(atomStr1) == 'O' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else if (trim(atomStr1) == 'N' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else
                print*,"Atom distance not included in checkbondDist.f90"
                stop
                endif

        else if ( trim(atomStr1) == 'N' ) then
                if (trim(atomStr1) == 'C' ) then
                        checkbondDist=bondCalc.lt.1.6d0
                else if (trim(atomStr1) == 'H' ) then
                        checkbondDist=bondCalc.lt.1.2d0
                else if (trim(atomStr1) == 'O' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else if (trim(atomStr1) == 'N' ) then
                        checkbondDist=bondCalc.lt.1.5d0
                else
                print*,"Atom distance not included in checkbondDist.f90"
                stop
                endif
        endif

end function
