










      Subroutine Assign_rcc_spin_state(Coefs,Nsizea,Nsizeb,Spin_state,
     &                                 Nroot,Irrep) 

      Implicit Double Precision (A-H, O-Z)
      
      Dimension Coefs(Nsizea+Nsizeb,Nsizea+Nsizeb)
      Dimension Nroot(8)
      Character*7 Spin_state(40)

      Thres1  = 1.0D-05
      Thres2  = 0.05D0


      If (Irrep .EQ. 1) Then
         Do Iroot = 1, Nroot(Irrep)
            Do I = 1, Nsizea

               If (Dabs(Coefs(I,Iroot)) .GT. Thres2  .AND. 
     &             Dabs(Coefs(I+Nsizeb,Iroot)) .GT. Thres2) Then

                    Diff = Coefs(I,Iroot)/Coefs(I+Nsizeb,Iroot)
                  
                    If (Dabs(1.0D0+Diff) .lt. Thres1) Then
                        Spin_state(Iroot) = "Triplet"
                    Else if (Dabs(1.0D0-Diff) .lt. Thres1) then
                        Spin_state(Iroot) = "Singlet"
                    Else 
                        Spin_state(Iroot) = "??????"
                    Endif

               Endif
            Enddo
         Enddo
      Else
         Do Iroot = 1, Nroot(Irrep)
            Do I = 1, Nsize

              If (Dabs(Coefs(I,Iroot)) .NE. Thres2 .AND.
     &            Dabs(Coefs(I+Nsizeb,Iroot)) .NE. Thres2) Then

                    Diff = Coefs(I,Iroot)/Coefs(I+Nsizeb,Iroot)

                    If (Dabs(1.0D0+Diff) .lt. Thres1) Then
                        Spin_state(Iroot) = "Triplet"
                    Else if (Dabs(1.0D0-Diff) .lt. Thres1) then
                        Spin_state(Iroot) = "Singlet"
                    Else 
                        Spin_state(Iroot) = "??????"
                    Endif 

CSSS                    If (Spin_state(Iroot) .EQ. "Singlet" .OR.
CSSS     &                  Spin_state(Iroot) .EQ. "Triplet")
CSSS     &                  Return
               Endif
            Enddo
         Enddo
      Endif 

      Return
      End 
