










        program calculateCCenergy 
        integer::iuhf,n_atm,nbas
        integer*8, dimension(2) :: scr




c icore.com : begin

c icore(1) is an anchor in memory that allows subroutines to address memory
c allocated with malloc. This system will fail if the main memory is segmented
c or parallel processes are not careful in how they allocate memory.

      integer icore(1)
      common / / icore

c icore.com : end





c istart.com : begin
      integer         i0, icrsiz
      common /istart/ i0, icrsiz
      save   /istart/
c istart.com : end

        CALL ACES_INIT(icore,i0,icrsiz,iuhf,.true.)

        call getrec(1,'JOBARC','NREALATM',1,n_atm)
        call getrec(1,'JOBARC','NBASTOT',1,nbas)
        call getrec(1,'JOBARC','NOCCORB',2,scr)

        print*, 'number of atoms',n_atm





        call aces_fin
        end program
