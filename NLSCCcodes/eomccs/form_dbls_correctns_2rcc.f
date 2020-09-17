










      Subroutine form_dbls_correctns_2rcc(Rcc_eigs,Rcc_vecl,Rcc_vecr,
     &                                    E2_rcc,Work,Length,Iuhf,
     &                                    Len_ov,Irrepx,Imult,
     &                                    Eom_sf_rccd,Eom_sf_drccd)

      Implicit Double Precision(A-H,O-Z)

      Dimension Rcc_vecl(Len_ov,Len_ov),Rcc_vecr(Len_ov,Len_ov)
      Dimension Rcc_eigs(Len_ov)
      Dimension Work(Length),E2_rcc(Len_ov,8)
      Dimension Nroot(8)
      Dimension Idoo(2),Idvv(2)
      Dimension Listijka(2),Listabci(2),Listmbej(3)
      Double Precision Ovlp 
      Logical Eom_sf_rccd,Eom_sf_drccd

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
      Data Scale,One,Onem,Dnull/1.0D0,1.0D0,-1.0D0,0.0D0/

      Call Getrec(20,'JOBARC','EESYMINF',NIRREP,NRoot)


      Do My_root = 1, Nroot(Irrepx)
         Ioff = 1
         Call Form_d2_4rcc(Work,Rcc_eigs,Len_ov,Length,Iuhf,Irrepx,
     +                     Listd0,1)
         Do Ispin = 3,3-2*Iuhf,-1  
            Call Zerolist(Work,Length,Listt2in-1+Ispin)
         Enddo 
         Call Zero_490_lists(Work,Length,Irrepx,Iuhf)
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

C This was done to compare with the UHF for debugging purposes. 
C            If (iuhf .eq. 0) Call Dscal(Len_ph_aa,1.0D0/Dsqrt(2.0D0),
C     +                                  Rcc_vecr(ioff,My_root),1)

            Call Putlst(Rcc_vecr(ioff,My_root),1,1,1,Ispin,Listc1) 
               
         Enddo 

C The doubles correction to RPA requires doubles target lists
C (lists 461 have already been formed). Therefore, we can proceed
C to built the first part of what is required for doubles correction,

         Listz0 = Listt2in 
         Listw1 = Listijka(2)
         Listw2 = Listabci(2)

         Call C1rpaint2c2a(Work,Length,Iuhf,Irrepx,Listw1,Listz0,
     &                     Listc1,Imult)
         Call C1rpaint2c2b(Work,Length,Iuhf,Irrepx,Listw2,Listz0,
     &                     Listc1,Imult)
         Call form_new_c2(Work,Length,Iuhf,Irrepx,Listz0,Listd0)

         Call C2inc1a(Work,Length,Iuhf,Irrepx,Listw2,Listz0,
     &                Listc1,Imult)
         Call C2inc1b(Work,Length,Iuhf,Irrepx,Listw1,Listz0,
     &                Listc1,Imult)
         Call Dbls_correctns_2rpa(Work,Length,Iuhf,RCC_vecl(1,My_root),
     &                            Len_ph_aa,Len_ph_bb,Len_ph_pq,
     &                            E2_rcc_aa,E2_rcc_bb,Listc1)
         E2_rcc(My_root,Irrepx) = (E2_rcc_aa + E2_rcc_bb)*Scale 
C
      Enddo     

      Return
      End 

      
