      Subroutine Process_rhf_ccd_roots(Zvec, Yvec, Eigval, Ndim, Irrep,
     &                                 Ispin,Rpa,
     &                                 Eom_sf_rccd,Eom_sf_drccd,
     &                                 Eom_s_rccd,Eom_s_drccd,
     &                                 Imult)
 
      Implicit Double Precision (A-H, O-Z)
      Integer Airrep, A 
      Logical Rpa,Rccd,DrcCd,Eom_s_rccd,Eom_s_drccd
      Logical Eom_sf_rccd,Eom_sf_drccd

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      integer Nroot(8)
      Dimension Zvec(Ndim,Ndim),Yvec(Ndim,Ndim)
      Dimension Eigval(Ndim)
      Parameter(Thresh = 0.050D0)
C
      CALL GETREC(20,'JOBARC','EESYMINF',NIRREP,NRoot)
C     Write(6,*) 'NIRREP:',NIRREP

      Iroot = 1
      Iuhf  = 0
      If (Irrep .EQ. 1) then
      Write(6,*)
      If (Eom_sf_rccd) Then
      Write(6,"(a,a)") " EOM(SF)-RCCD (RPA) excitation energies",
     +                 " and state vectors" 
      Write(6,"(a,a)")   " (only coefs. which are greater than",
     +                   " 0.05 is printed.)"
      Elseif (Eom_sf_drccd) Then
      Write(6,"(a,a)") " EOM(SF)-DRCCD (DRPA) excitation energies"
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

      Write(6,"(a,a)")  "---------------------------------------",
     &                  "-------------------"
      Write(6,*)
      Write(6,"(a,1x,i1)") "The excited states of symmetry block:",irrep

      Do Iroot = 1, NRoot(Irrep)
         Write(6,*)
         Write(6,"(a,F12.6,a)") " The excitation energy = ",
     &                         Eigval(Iroot)*27.2113961D0,
     &                         "eV"
         If (Imult .EQ. 1) Then
            Write(6,"(a,a)") " The spin-state        = ", "Singlet"
         Elseif (Imult .EQ. 2) then 
            Write(6,"(a,a)") " The spin-state        = ", "Triplet"
         Endif 

         Write(6,*)
         Icount =  1
         Do Iirrep = 1, Nirrep
            Airrep = Dirprd(Iirrep, Irrep)

            Do i =1, Pop(Iirrep, Ispin)
               Do a =1, Vrt(AIrrep, Ispin)
 
                  If (Dabs(Zvec(Icount,Iroot)) .gt. Thresh) then
                   
                     Write(6,1) I, Iirrep, A, Airrep, 
     &                          Zvec(Icount,Iroot)
                
                  Endif 
                  Icount = Icount + 1
               Enddo
            Enddo 
         Enddo 
         Write(6,"(a,a)")  "---------------------------------------",
     &                     "-------------------"
      Enddo 
C
 1    Format(3X,I4,3X, ' [',I1,'],  ',1X, I4,3X, ' [',I1,']   ;',
     +       10X, F12.6)

      Return
      End
