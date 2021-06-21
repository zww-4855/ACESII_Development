










      Subroutine Block_david_driver(Irrepx,Work,Maxcor,Nsize,Iuhf,Iside,
     +                              Tol,Maxiter)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)
      Logical Converged



c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end



c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c flags.com : begin
      integer        iflags(100)
      common /flags/ iflags
c flags.com : end
c flags2.com : begin
      integer         iflags2(500)
      common /flags2/ iflags2
c flags2.com : end
      
      Common /Calcinfo/Nroot(8)

      Parameter(Maxblocks=10)

C This code assume that the TDA vectors are availble for the 
C requested number of roots via estate_sym (and less <= 10 for
C each symmtery block)

      Nblocks = Nroot(Irrepx)
      If (Nblocks .Gt. Maxblocks) Then
         Write(6,"(2a,i2,2a,i2)") " The requested number of roots",
     +                            " for irrep", Irrep, "greater",
     +                            " than the aximum allowed value", 
     +                               Maxblocks
         Call Errex
      Endif

C Pick the TDA eigenvectors for each root per symmetry block and
C carry out the multiplication and accmulate 

      Do Ispin=3,3-2*Iuhf,-1
         Call Zerolist(Scr,Maxcor,443+Ispin)
      Enddo 

C This is the start of the block Davidson loop

       Converged = .False.
       Iter_count = 0

       Do while (Converged) 
          Iter_count = Iter_count + 1
          Ioffc = 1

C It is assuemed that we need at least Nsize*Nblocks of memory to
C proceed. Note that this does not count any additional requirements 
C within blokdavison. Memory checks are essential since this is 
c very memory intensive process. 

           Ioff_hc = 1
           Do Iblck = 1, Nblocks 
              I000 = Ione
              Iend = I000 + Nsize*Nblocks
              Memleft = Maxcor - Iend 
              If (Iend .Ge. Maxcor) Call Insmem("block_david_driver",
     +                                            Iend,Maxcor)
C This is start of the loop over No. blocks. First iteration starts 
C with the TDA vectors of (Nblocks). Subsequent interation will start
C use Davidson extrapolated vectors.

              If (Iter_count .Eq. 1) Then
                 Call Getlst(Work,Iblck,1,1,Irrepx,94)
                 If (Iuhf .Eq. 0) Then
                    Fact = One/Dsqrt(Two)
                    Call Dscal(Irpdpd(Irrepx,9),Fact,Work,1)
                 Endif
                 Ioffc = Ione 
                 Do Ispin=1,Iuhf+1
                    Call Putlst(Work(Ioffc),1,1,Ispin,490)
                    Ioffc = Ioffc+Irpdpd(Irrepx,8+Ispin)
                 Enddo 
              Else 
                 Call Lanczos_dumpvec(Irrepx,Work(Ioff_hc),Maxcor,
     +                                 Iuhf,490,1,443,Nsize,.False.)
              Endif 
           
              Call Hbarxc(Work,Maxcor,Iuhf,Iside,Irrepx)

              Call Loadvec1(Irrepx,Work(Ioff_hc),Maxcor,Iuhf,490,2,460,
     +                      Nsize,.False.)
              Ioff_hc = Ioff_hc + Nsize
   
           Enddo 

C After this loop, we have Nsize*Nblocks of HC vectors and we can call
C block davidson procedure and return Nsize*Nblocks of extraploted 
C vectors. Starting Ioff_hc=1, we have Nsize*Nblock new HC vectors. 

           Ioff_hc = 1

CSSS       Call Blcokdave(Irepx,Work(Ioff_hc),Work(Iend),Memleft,
CSSS     +                Nbocks,Irrepx,Tol,Converged)

C When returned this call expects extrpolated vectors arrive in
C Work(ioff_hc) where Ioff_hc=1 

      Enddo 
  
      Return 
      End
   
