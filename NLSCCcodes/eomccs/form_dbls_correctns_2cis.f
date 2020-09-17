










      Subroutine form_dbls_correctns_2cis(Work,Length,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Length)
      Dimension Nroot(8)
      Dimension Listijka(2),Listabci(2),Listmbej(3)
      Double Precision Ovlp 
      Dimension E2_rcc(100,8)
      Dimension Idoo(2),Idvv(2)

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c sympop.com : begin
      integer         irpdpd(8,22), isytyp(2,500), id(18)
      common /sympop/ irpdpd,       isytyp,        id
c sympop.com : end
c syminf.com : begin
      integer nstart, nirrep, irrepa(255), irrepb(255), dirprd(8,8)
      common /syminf/ nstart, nirrep, irrepa, irrepb, dirprd
c syminf.com : end

      Common /Tdalist/ Listetda,Listvtda 
      Data Listc1,Listt2,Listt2in,Listt2rs /490,444,461,440/
      Data Listijka /107,7/
      Data Listabci /127,27/
      Data Listd0   /447/
      Data Listmbej /54,54,254/
      Data Scale1,Imult,One,Onem,Dnull/1.0D0,1,1.0D0,-1.0D0,0.0D0/

      Call Getrec(20,'JOBARC','EESYMINF',NIRREP,NRoot)

      Scale2 = One + Dfloat(1-Iuhf)
      do Irrepx = 1, Nirrep

         Call Newlst(Irrepx,Work,Length,Iuhf)
         Len_ov = Irpdpd(Irrepx,9) + Irpdpd(Irrepx,10)*Iuhf 
         Itda_vec = 1
         It2c_vec = Itda_vec + Len_oV
         Iw2c_vec = It2c_vec + Len_ov 
         Iend     = Iw2c_vec + Len_ov
         Length   = Length - Iend
         If (Iend .Ge. Length) Call Insmem("dbls_correction",
     &                                      Iend,Length)

      Do My_root = 1, Nroot(Irrepx)

         Call Getlst(Work(Itda_vec),My_root,1,1,Irrepx,Listvtda)
         Call Getlst(Root,My_root,1,1,Irrepx,Listetda)

         Ioff = 1
         Call Form_d2_4rcc(Work(Iend),Root,Len_ov,Length,Iuhf,
     &                     Irrepx,Listd0,1)
         Do Ispin = 3,3-2*Iuhf,-1  
            Call Zerolist(Work(Iend),Length,Listt2in-1+Ispin)
         Enddo 
         Call Zero_490_lists(Work(Iend),Length,Irrepx,Iuhf)
         E2_rcc_aa = 0.0D0
         E2_rcc_bb = 0.0D0

         Do Ispin = 1, (Iuhf+1)
      
C Lets write the RPA vectors to the same lists that we used to
C store the CIS vectors. That will make the rest of the code's 
C work flow easier.

            Len_ph_aa = Irpdpd(Irrepx,9)
            Len_ph_bb = Irpdpd(Irrepx,10)
            Len_ph_pq = Len_ph_aa + Len_ph_bb * Iuhf

            Ioff = Ioff + (Ispin - 1) * Len_ph_aa
            Call Putlst(Work(Itda_vec-1+Ioff),1,1,1,Ispin,Listc1) 
               
         Enddo 

C The doubles correction to RPA requires doubles target lists
C (lists 461 have already been formed). Therefore, we can proceed
C to built the first part of what is required for doubles correction,

         Listz0 = Listt2in 
         Listw1 = Listijka(2)
         Listw2 = Listabci(2)

         Call C1rpaint2c2a(Work(Iend),Length,Iuhf,Irrepx,Listw1,Listz0,
     &                     Listc1,Imult)
         Call C1rpaint2c2b(Work(Iend),Length,Iuhf,Irrepx,Listw2,Listz0,
     &                     Listc1,Imult)
         Call form_new_c2(Work(Iend),Length,Iuhf,Irrepx,Listz0,Listd0)

         Call C2inc1a(Work(Iend),Length,Iuhf,Irrepx,Listw2,Listz0,
     &                Listc1,Imult)
         Call C2inc1b(Work(Iend),Length,Iuhf,Irrepx,Listw1,Listz0,
     &                Listc1,Imult)
         Call Dbls_correctns_2rpa(Work(Iend),Length,Iuhf,Work(Itda_vec),
     &                            Len_ph_aa,Len_ph_bb,Len_ph_pq,
     &                            E2_rcc_aa,E2_rcc_bb,Listc1)
         E2_rcc(My_root,Irrepx) = (E2_rcc_aa + E2_rcc_bb)*Scale1

C
C This term applies only to CIS. For RPA and DRPA this term is alrready 
C in HBAR

         Call Resort(Work(Iend),Length,Iuhf,1,44,34)
         Call Gt2xf(Work(It2c_vec),1,Irrepx,34,Listc1,0,Work(Iend),
     &              Length,Iuhf,0)
         Call Resort(Work(Iend),Length,Iuhf,1,14,34)
         Call Gt2xf(Work(Iw2c_vec),1,Irrepx,34,Listc1,0,Work(Iend),
     &              Length,Iuhf,0) 

         Do Ispin = 1, Iuhf+1
            Length_aa = Irpdpd(Irrepx,9)
            Length_bb = Irpdpd(Irrepx,10) 
            Itoff_aa  = It2c_vec
            Iwoff_aa  = Iw2c_vec  
            Itoff_bb  = It2c_vec + Length_aa 
            Iwoff_bb  = Iw2c_vec + Length_aa 
      
            if (Ispin .Eq. 1) E2_rcc_aa = Ddot(Length_aa,
     &                                         Work(Itoff_aa),
     &                                         1,Work(Iwoff_aa),1)
            if (Ispin .Eq. 2) E2_rcc_bb = Ddot(Length_bb,
     &                                         Work(Itoff_bb),
     &                                         1,Work(Iwoff_bb),1)
         Enddo 
   
         E2_rcc(My_root,Irrepx) = E2_rcc(My_root,Irrepx) + E2_rcc_aa + 
     &                            E2_rcc_bb
         E2_rcc_oo = Dnull
         E2_rcc_vv = Dnull
         Itda_vec_off = 1
         Do Ispin = 1, Iuhf+1
            Length_ai = Irpdpd(Irrepx,9)
            Idoo(1) = Iend 
            Idoo(2) = Idoo(1) + Nfmi(1)
            Idvv(1) = Idoo(2) + Iuhf*Nfmi(2) 
            Idvv(2) = Idvv(1) + Nfea(1) 
            Iend    = Idvv(2) + Iuhf*Nfea(1) 
            Max_scr = Max(Nfmi(1),Nfmi(2),Nfea(1),Nfea(2))
            Iend    = Iend + Max_scr
            Length  = Length - Iend 
            If (Iend .Ge. Length) Call Insmem("dbls_correction",
     &                                         Iend,Length)
            
            Call Dzero(Work(Iend),Nfmi(Ispin))
            Call Modgij(1,Irrepx,Irrepx,Work(Idoo(Ispin)),Work(Iend),
     &                  Work(Itda_vec_off),Work(Itda_vec_off),Ispin,
     &                  Onem)
            Call Dzero(Work(Iend),Nfea(Ispin))
            Call Modgab(1,Irrepx,Irrepx,Work(Idvv(Ispin)),Work(Iend),
     &                  Work(Itda_vec_off),Work(Itda_vec_off),Ispin,
     &                  One)
            Call Gformg(1,1,44,14,0,Work(Iend),Length,0,One,One,Iuhf)

            Length_oo = Nfmi(Ispin)
            Length_vv = Nfea(Ispin)
            Call Getlst(Work(Iend),1,1,1,Ispin,91)
            Write(6,*) Scale2
            E2_rcc_oo = E2_rcc_oo + Ddot(Length_oo,Work(Iend),1,
     &                                   Work(Idoo(Ispin)),1)
            Call Getlst(Work(Iend),1,1,1,Ispin,92)
            E2_rcc_vv = E2_rcc_vv + Ddot(Length_vv,Work(Iend),1,
     &                                   Work(Idvv(Ispin)),1)
            Itda_vec_off = Itda_vec_off +  Length_ai
         Enddo 

         E2_rcc(My_root,Irrepx) = E2_rcc(My_root,Irrepx) + E2_rcc_oo +
     &                            E2_rcc_vv
        
      Enddo     
      Enddo 

      Return
      End 

      
