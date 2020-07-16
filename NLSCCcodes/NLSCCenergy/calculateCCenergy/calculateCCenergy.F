        program calculateCCenergy 
        integer::iuhf,n_atm,nbas
        integer*8, dimension(2) :: scr


#include "icore.com"
#include "istart.com"

        CALL ACES_INIT(icore,i0,icrsiz,iuhf,.true.)

        call getrec(1,'JOBARC','NREALATM',1,n_atm)
        call getrec(1,'JOBARC','NBASTOT',1,nbas)
        call getrec(1,'JOBARC','NOCCORB',2,scr)

        print*, 'number of atoms',n_atm





        call aces_fin
        end program
