










      Subroutine dcc_hbar_fmi(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Double Precision Mone, One

      Dimension Work(Length)



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

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 

C
C Construct the DCC FMI intermediates (this is a duplicate from hbar)
C
      Irrepx = 1
      Write(6,"(a)") " Checksums of FMI : At Entry"
      Write(6,*)
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,9,91)
      Call Checksum("FMI", Work, IHHA_LENGTH,S)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,91)
      IF (IUHF.gt.0) Call Checksum("Fmi", Work, IHHB_LENGTH,S)

      Call pdcc_quad2(Work,Length,iuhf,Fmi_scale,8,1)

      Write(6,*)
      Write(6,"(a)") " Checksums of FMI : At Exit"
      Write(6,*)
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,9,91)
      Call Checksum("FMI", Work, IHHA_LENGTH,S)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,91)
      IF (IUHF.gt.0) Call Checksum("Fmi", Work, IHHB_LENGTH,S)
      Write(6,*)

      Return
      End
