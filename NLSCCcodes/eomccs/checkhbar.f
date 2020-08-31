










      Subroutine CheckHbar(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     


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



c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 
































































































































































































c This common block contains the IFLAGS and IFLAGS2 arrays for JODA ROUTINES
c ONLY! The reason is that it contains both arrays back-to-back. If the
c preprocessor define MONSTER_FLAGS is set, then the arrays are compressed
c into one large (currently) 600 element long array; otherwise, they are
c split into IFLAGS(100) and IFLAGS2(500).

c iflags(100)  ASVs reserved for Stanton, Gauss, and Co.
c              (Our code is already irrevocably split, why bother anymore?)
c iflags2(500) ASVs for everyone else

      integer        iflags(100), iflags2(500)
      common /flags/ iflags,      iflags2
      save   /flags/





      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C First Hbar(A,E), Hbar(M,I) and Hbar(M,E) 

      IRREPX =1
      HHA_LENGTH = IRPDPD(IRREPX,21)
      HHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,1,91)
      Call Checksum("FMI", Work, HHA_LENGTH,S)
      IF (ISPAR) THEN
         Call Getlst(Work,1,1,1,9,91)
         Call Checksum("FMI", Work, HHA_LENGTH,S)
      ENDIF
      IF (UHF) Then
         Call Getlst(Work,1,1,1,2,91)
         Call Checksum("Fmi", Work, HHB_LENGTH,S)
         IF (ISPAR) THEN
            Call Getlst(Work,1,1,1,10,91)
            Call Checksum("Fmi", Work, HHB_LENGTH,S)
         ENDIF
      ENDIF

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,93)
      Call Checksum("FME", Work, HPA_LENGTH,s)
      IF (UHF) Call Getlst(Work,1,1,1,2,93)
      IF (UHF) Call Checksum("Fme", Work, HPB_LENGTH,S)

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,1,92)
      Call Checksum("FEA", Work, PPA_LENGTH, s)
      If (Ispar) THEN
         CALL Getlst(Work,1,1,1,9,92)
         Call Checksum("FEA", Work, PPA_LENGTH, s)
      Endif
      IF (UHF) THEN
        Call Getlst(Work,1,1,1,2,92)
        Call Checksum("Fea", Work, PPB_LENGTH,S)
        If (Ispar) Then
           Call Getlst(Work,1,1,1,10,92)
           Call Checksum("Fea", Work, PPB_LENGTH,S)
        Endif
      Endif

C Hbar(MN,IJ)
  
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,51),ISYTYP(2,51))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,52),ISYTYP(2,52))
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         Call Getall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Checksum("AAAA_MNIJ", Work, AAAA_LENGTH_MNIJ,S)
         Call Getall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Checksum("BBBB_MNIJ", Work, BBBB_LENGTH_MNIJ,S)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ,S)
      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ,S)
      Endif

C Hbar(AI,BC)
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         BBBB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
         ABAB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         Call Getall(Work, AAAA_LENGTH_AIBC, IRREPX, 27)
         Call Checksum("AAAA_AIBC", Work, AAAA_LENGTH_AIBC,S)
         Call Getall(Work, BBBB_LENGTH_AIBC, IRREPX, 28)
         Call Checksum("BBBB_AIBC", Work, BBBB_LENGTH_AIBC,S)
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 29)
         Call Checksum("ABAB_AIBC", Work, ABAB_LENGTH_AIBC,S)
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC,S)
      Else
CSSS         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
CSSS         Call Getall(Work, AAAA_LENGTH_AIBC, IRREPX, 27)
CSSS         Call Checksum("AAAA_AIBC", Work, AAAA_LENGTH_AIBC,S)
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC,S)
      Endif 

C Hbar(IJ,KA)
      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Checksum("AAAA_IJKA", Work, AAAA_LENGTH_IJKA,S)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 8)
         Call Checksum("BBBB_IJKA", Work, BBBB_LENGTH_IJKA,S)
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 9)
         Call Checksum("ABAB_IJKA", Work, ABAB_LENGTH_IJKA,S)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA,S)
      Else
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Checksum("AAAA_IJKA", Work, AAAA_LENGTH_IJKA,S)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA,S)
      Endif

C Hbar(MB,EJ)
      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
	 ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ,S)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ,S)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ,S)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ,S)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ,S)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ,S)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ,S)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ,S)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ,S)
      Endif

C Hbar(AB,CI)
      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Checksum("AAAA_ABCI", Work, AAAA_LENGTH_ABCI,S)
         Call Getall(Work, BBBB_LENGTH_ABCI, IRREPX, 128)
         Call Checksum("BBBB_ABCI", Work, BBBB_LENGTH_ABCI,S)
         Call Getall(Work, ABAB_LENGTH_ABCI, IRREPX, 129)
         Call Checksum("ABAB_ABCI", Work, ABAB_LENGTH_ABCI,S)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Checksum("ABBA_ABCI", Work, ABBA_LENGTH_ABCI,S)
      Else
CSSS         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
CSSS         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
CSSS         Call Checksum("AAAA_ABCI", Work, AAAA_LENGTH_ABCI,S)
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Checksum("ABBA_ABCI", Work, ABBA_LENGTH_ABCI,S)
      Endif

C Hbar(IA,JK)
      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Checksum("AAAA_IAJK", Work, AAAA_LENGTH_IJKA,S)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Checksum("BBBB_IAJK", Work, BBBB_LENGTH_IJKA,S)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 109)
         Call Checksum("ABAB_IAJK", Work, ABAB_LENGTH_IJKA,S)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Checksum("ABBA_IAJK", Work, ABBA_LENGTH_IJKA,S)
      Else
CSSS         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
CSSS         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
CSSS         Call Checksum("AAAA_IAJK", Work, AAAA_LENGTH_IJKA,S)
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Checksum("ABBA_IAJK", Work, ABBA_LENGTH_IJKA,S)
      Endif

C Hbar(AB,CD)
      Write(*,*)
      IF (IFLAGS(93) .EQ. 0) THEN
      IF (UHF) Then
         AAAA_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,231),ISYTYP(2,231))
         BBBB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,232),ISYTYP(2,232))
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         Call Getall(Work, AAAA_LENGTH_ABCD, IRREPX, 231)
         Call Checksum("AAAA_ABCD", Work, AAAA_LENGTH_ABCD,S)
         Call Getall(Work, BBBB_LENGTH_ABCD, IRREPX, 232)
         Call Checksum("BBBB_ABCD", Work, BBBB_LENGTH_ABCD,S)
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD,S)
      Else
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD,S)
      Endif
      Endif 


      Return
      End
