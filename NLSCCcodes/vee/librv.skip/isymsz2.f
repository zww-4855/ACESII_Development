      INTEGER FUNCTION ISYMSZ2(IrrepX, ITYPL,ITYPR)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: isymsz2.f,v 1.1.1.1 2003/04/02 19:22:07 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     isymsz2 -- Number of elements in a list of given overall symmetry
C
C SYNOPSIS
      Integer IrrepX, ITypL, ITypR
C
C ARGUMENTS
C     IrrepX  Overall irrep of list (IN)
C     ITypL   Symmetry type of LHS indices (IN)
C     ITypR   Symmetry type of RHS indices (IN)
C
C DESCRIPTION
C     Computes total size of all irreps of an MOIO list of overall
C     symmetry IRREPX with LHS indices of type ITYPL and RHS indices
C     of type ITYPR.  The symmetry type inputs can be interchanged
C     without changing the result.
C
C     Calling this routine with IRREPX=1 will give the same result as
C     ISYMSZ.
C
C COMMON BLOCKS
      Integer NStart, NIrrep, IrrepA, IrrepB, DirPrd
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                 DIRPRD(8,8)
      Integer IrpDPD, ISyTyp, NTot
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     IrrepR  Irrep of RHS
C     IrrepL  Irrep of LHS
C
      Integer IrrepR, IrrepL
C
      ISYMSZ2=0
      DO 10 IRREPR=1,NIRREP
         IrrepL = DirPrd(IrrepR, IrrepX)
         ISYMSZ2=ISYMSZ2+IRPDPD(IRREPL,ITYPL)*IRPDPD(IRREPR,ITYPR)
 10   CONTINUE
      RETURN
      END
