










      Subroutine Check_vee_lists_init(Work,Maxcor,Irrepx,Iuhf,Iside)
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

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Write(6,*) 

      If (iside .Eq. 1) then 

      Call Getlst(Work,1,1,1,1,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif

      Elseif (Iside .Eq. 2) Then

      Call Getlst(Work,1,1,1,1,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif 

      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
      call getall(Work(1),length,irrepx,444)
      call checksum("List-444;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
      call getall(Work(1),length,irrepx,445)
      call checksum("List-445;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(ab,ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(Ai,Bj):",Work(1),length,s)
      Endif 

      elseif (iside .eq. 2) then

      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
      call getall(Work(1),length,irrepx,444)
      call checksum("List-444;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
      call getall(Work(1),length,irrepx,445)
      call checksum("List-445;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(ab,ij):",Work(1),length,s)
      Else 
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(Ai,Bj):",Work(1),length,s)
      Endif 
      Endif 

      Return 
      End
