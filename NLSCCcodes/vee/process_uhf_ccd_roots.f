










      Subroutine Process_uhf_ccd_roots(Zvec,Yvec,Eigval,Ndim,Nsizea,
     &                                 Nsizeb,Irrep,Iuhf,Rpa,
     &                                 Eom_sf_rccd,Eom_sf_drccd,
     &                                 Eom_s_rccd,Eom_s_drccd)

      Implicit Double Precision (A-H, O-Z)
      Integer Airrep, A 
      Logical Rpa,Rccd,DrcCd,Eom_s_rccd,Eom_s_drccd
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

      Dimension Zvec(Ndim,Ndim),Yvec(Ndim,Ndim)
      Dimension Eigval(Ndim),Nroot(8)
      Parameter(Thresh = 0.050D0)
      Character*7 Spin_state(40)

      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NRoot)
C
      If (Irrep .EQ. 1) then
      Write(6,*)
      If (Eom_sf_rccd) Then
      Write(6,"(a,a)") " EOM(SF)-RCCD (RPA) excitation energies",
     +                 " and state vectors"
      Write(6,"(a,a)")   " (only coefs. which are greater than", 
     +                   " 0.05 is printed.)"
      Elseif (Eom_sf_drccd) Then
      Write(6,"(a,a)") " EOM(SF)-DRCCD (RPA) excitation energies"
     +                 " and  state vectors"
      Write(6,"(a,a)")   " (only coefs. which are greater than", 
     +                   " 0.05 is printed.)"
      Elseif (Eom_s_rccd) Then
      Write(6,"(a,a)") " EOM(S)-RCCD excitation energies",
     +                 " and  state vectors"
      Write(6,"(a,a)")   " (only coefs. which are greater than", 
     +                   " 0.05 is printed.)"
      Elseif (Eom_s_drccd) Then
      Write(6,"(a,a)") " EOM(S)-DRCCD excitation energies",
     +                 " and  state vectors"
      Write(6,"(a,a)")   " (only coefs. which are greater than", 
     +                   " 0.05 is printed.)"
      Endif 
      Endif 

      Call Assign_rcc_spin_state(Zvec,Nsizea,Nsizeb,Spin_state,
     +                           Nroot,Irrep)

      Write(6,"(a,a)")  "---------------------------------------",
     &                  "-------------------"
      Write(6,*)
      Write(6,"(a,1x,i1)") "The excited states of symmetry block:",
     &                      irrep

      If  (Nirrep .EQ. 1) Then

          Do Iroot = 1, Nroot(Irrep)
          Do Ispin = 1, 2
             Icount = (Ispin - 1) * Nsizea  + 1

             If (Ispin .Eq. 1) Then
             Write(6,*)
             Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                                Eigval(Iroot)*27.2113961D0,
     &                             " eV"
             Write(6,"(a,7a)") " The spin-state        = ", 
     &                           Spin_state(Iroot)
             Endif 
             Write(6,*)
             If (Ispin .EQ. 1) Write(6,"(a)") " The AA coeffcients" 
             If (Ispin .EQ. 2) Write(6,"(a)") " The BB coeffcients" 
             Do i =1, Pop(Irrep, Ispin)
                Do a =1, Vrt(Irrep, Ispin)
 
                   If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                    
                      Write(6,1) I, Irrep, A, Irrep, 
     &                           Zvec(Icount,Iroot)
                    
                   Endif 
                   Icount = Icount + 1

                Enddo 
             Enddo 
          Enddo 
          Write(6,"(a,a)")  "---------------------------------------",
     &                      "-------------------"
          Enddo 
      Else 
          Do Iroot = 1, Nroot(Irrep)
          Do Ispin = 1, 2
             Icount = (Ispin - 1) * Nsizea + 1
 
             If (Ispin .EQ. 1) then
             Write(6,*)
             Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                                Eigval(Iroot)*27.2113961D0,
     &                             " eV"
             Write(6,"(a,7a)") " The spin-state        = ", 
     &                           Spin_state(Iroot)
             Endif 
             Write(6,*)
             If (Ispin .EQ. 1) Write(6,"(a)") " The AA coeffcients" 
             If (Ispin .EQ. 2) Write(6,"(a)") " The BB coeffcients" 
             Do Iirrep = 1, Nirrep
                Airrep = Dirprd(Iirrep, Irrep)

                Do i =1, Pop(Iirrep, Ispin)
                   Do a =1, Vrt(AIrrep, Ispin)
 
                      If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                       
                         Write(6,1) I, Iirrep, A, Airrep, 
     &                              Zvec(Icount,Iroot)
                    
                      Endif 
                      Icount = Icount + 1
                   Enddo
                Enddo 
             Enddo 
          Write(6,"(a,a)")  "---------------------------------------",
     &                      "-------------------"
          Enddo 
          Enddo
      Endif 

 1    Format(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     +       10X, F12.6)

      Return
      End
