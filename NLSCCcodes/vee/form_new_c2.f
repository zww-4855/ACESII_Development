










      Subroutine form_new_c2(Work,Length,Iuhf,Irrepx,Listz0,Listd0)
      
      Implicit Double Precision(A-H,O-Z)

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
c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end

C ABAB Block

      Ioff = 1 
      Do Irrepr = 1, Nirrep
         Irrepl = Dirprd(Irrepr,Irrepx)

         Listz = Listz0 + 2
         Listd = Listd0 + 3
         Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
         Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
         Nsize = Nrow * Ncol

         I000  = 1
         I010  = I000  + Nsize
         I020  = I010  + Nsize
         IScr1 = I020  + Nsize
         Iscr2 = Iscr1 + Nsize
         Iend  = Iscr2 + Nsize 
         If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                      Iend,Length)
         Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

         Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
C         If (Iuhf .Eq. 0) Then
C            Call Spinad1(Irrepr,Pop,Nrow,Work(Ioff),Work(Iscr1),
C     +                   Work(Iscr2))
C            Ioff = Ioff + Nrow*Ncol
C         Endif
         Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)

      Enddo
      If (Iuhf .eq. 0) Return 

C AAAA and BBBB blocks

      Do Ispin = 1, 1+Iuhf
         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)

            Listz = Listz0 + Ispin - 1
            Listd = Listd0 + Ispin 
            Nrow  = Irpdpd(Irrepl,Isytyp(1,Listz))
            Ncol  = Irpdpd(Irrepr,Isytyp(2,Listz))
            Nsize = Nrow * Ncol

            I000 = 1
            I010 = I000 + Nsize 
            I020 = I010 + Nsize
            Iend = I020
            If (Iend .Ge. Length) Call Insmem("form_new_c2",
     +                                         Iend,Length)

            Call Getlst(Work(I000),1,Ncol,1,Irrepr,Listz)
            Call Getlst(Work(I010),1,Ncol,1,Irrepr,Listd)

            Call Vecprd(Work(I000),Work(I010),Work(I000),Nsize)
            Call Putlst(Work(I000),1,Ncol,1,Irrepr,Listz)
         Enddo 
      Enddo 


      Return
      End
