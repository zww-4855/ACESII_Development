      SUBROUTINE H4T2ALL(ICORE,MAXCOR,ISTART,IOFFT1,ISPIN,IUHF)
C
C DRIVER FOR RING CONTRIBUTIONS TO G(IJ,AB) AMPLITUDES.
C
C THIS CONTRIBUTION IS:
C
C        ab              ab
C       Z   = P(ij|ab) Q
C        ij              ij
C
C WHERE P(ij|ab) IS THE USUAL PERMUTATION OPERATOR, AND THE QUANTITY
C  Q IS DEFINED BY 
C
C        ab           ae
C       Q    =  SUM  T   H4(mejb)
C        ij     m,e   im                                         
c
C 
C  WITH H4 AS A TWO-PARTICLE INTERMEDIATE
C
C FOR THE VARIOUS SPIN CASES, Z IS DEFINED BY:
C
C       AB    AB  BA  AB  BA
C      Z  =  Q  +Q  -Q  -Q     (AAAA)
C       IJ    IJ  JI  IJ  JI
C
C       Ab    Ab  bA  Ab  bA
C      Z  =  Q  +Q  +Q  +Q     (ABAB)
C       Ij    Ij  jI  Ij  jI
C
C       aB    aB  Ba  aB  Ba
C      Z  =  Q  +Q  +Q  +Q     (ABBA)
C       Ij    Ij  jI  Ij  jI
C    
C
C  THE BBBB, BABA AND BAAB SPIN CASES CAN BE OBTAINED FROM THESE EQUATIONS
C   BY CHANGING ALL ALPHA LABELS TO BETA AND ALL BETA LABELS TO ALPHA.
C
C
C  THE SPIN CASES FOR Q ARE GIVEN BY:
C
C  AB             AE           Ae                
C Q  = Q(AAAA) = T  H4(EMJB) + T  H4(meJB)
C  IJ             IM           Im                   
C
C  ab             ae           aE                
C Q  = Q(BBBB) = T  H4(mejb) + T  H4(MEjb)
C  ij             im           iM                   
C
C  Ab             AE           Ae                
C Q  = Q(ABAB) = T  H4(MEjb) + T  H4(mejb)
C  Ij             IM           Im                   
C
C  aB             ae           aE
C Q  = Q(BABA) = T  H4(meJB) + T  H4(MEJB)
C  iJ             im           iM
C
C  Ab               eA
C Q  = Q(BAAB) = - T  H4(MeJb)
C  iJ               iM      
C
C  aB               Ae
C Q  = Q(ABBA) = - T  H4(mEjB)
C  Ij               Im
C
C (ALL QUANTITIES ARE SUMMED OVER m AND e).
C
C
C HAVING THE Q QUANTITIES, THE RING CONTRIBUTION TO THE G(IJ,AB) AMPLITUDES
C  (THE Z(ij,ab)) CAN BE COMPUTED FROM:
C
C                 AB             AB
C SPIN CASE AA : G   = P(IJ|AB) Q
C                 IJ             IJ
C
C                 ab             ab
C SPIN CASE BB : G   = P(ij|ab) Q
C                 ij             ij
C
C                 Ab    Ab    bA    bA    Ab
C SPIN CASE AB : G   = Q   + Q   - Q   - Q  
C                 Ij    Ij    jI    Ij    jI
C
C
C    OR G(AB)=Q(ABAB)+Q(BABA)-Q(ABBA)-Q(BAAB)
C
C NOTE THAT WE CALCULATE FOR THE SPIN CASES AAAA AND BBBB THE POSITIVE
C Q(IA,JB) WHILE FOR ALL OTHER SPIN CASES THE NEGATIVE Q(IA,JB) IS     
C EVALUATED.
C
C WE NEED ALSO THIS FUNNY TAU : T(IJ,AB)-2*T(I,B)*T(J,A)
C FOR THE FISRT TWO SPIN CASES. THE FACTOR IS TWO AND THIS CONTRIBUTION
C IS FORMED IN F2TAU.
C 
C THE FACTORS ARE CHOSEN SO THAT HERE THE OVERALL FACTOR IS ONE OR ONEM,
C THE FACTOR OF HALF IS TAKEN CARE OF IN HTRNGDRV  
C
CEND 
C
C  CODED OCTOBER/93  JG
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT
      CHARACTER*4 SPCASE
      DIMENSION ICORE(MAXCOR),IOFFT1(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFO/NOCCO(2),NVRTO(2) 
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE /1.0D0/,TWO/2.0D0/
C
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
C 
      MXCOR=MAXCOR-IINTFP*NT(1)-IUHF*IINTFP*NT(2)
C      
      INCREM=0
      IF(ISPIN.LT.3)THEN
       IF(ISPIN.EQ.1)SPCASE='AAAA'
       IF(ISPIN.EQ.2)SPCASE='BBBB'
C
C   Q(AAAA) AND Q(BBBB).  COMMENTS BELOW REFER TO THE AAAA CASE ONLY.
C                            FLIPPING THE SENSE OF A AND B WILL GIVE
C                            THE EQUATIONS FOR THE BBBB CASE.
C
C               AB         EA
C   SOLVE FOR  Q   = SUM  T   *[- H4(MEJB)]  [FIRST PART OF Q(AAAA)]
C               IJ   M,E   IM
C
       MAXSIZ=0
       LISTT=33+ISPIN
       LISTW=53+ISPIN
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'TxW',MAXSIZ,INCREM,39+ISPIN,.TRUE.,
     &               ICORE(IOFFT1(ISPIN)),ICORE(IOFFT1(ISPIN)),
     &               POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &               VRT(1,ISPIN),TWO,SPCASE,IUHF)
C
C               AB         Ae
C   SOLVE FOR  Q   = SUM  T    H(meJB)  [SECOND PART OF Q(AAAA)].
C               IJ   m,e   Im
C
       LISTW=58-ISPIN
       LISTT=38-ISPIN
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'TxW',MAXSIZ,1,39+ISPIN,.FALSE.,
     &               ICORE(IOFFT1(ISPIN)),ICORE(IOFFT1(ISPIN)),
     &               POP(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),
     &               VRT(1,ISPIN),TWO,SPCASE,IUHF) 
C
C RESORT AND ANTISYMMETRIZE X2 INCREMENT, AND THEN AUGMENT THE 
C  T2 INCREMENT LIST. 
C
       NOCC=NOCCO(ISPIN)
       NVRT=NVRTO(ISPIN)
       NVSIZ=(NVRT*(NVRT-1))/2
       NOSIZ=(NOCC*(NOCC-1))/2
       NTOTSZ=ISYMSZ(ISYTYP(1,39+ISPIN),ISYTYP(2,39+ISPIN))
       NTOTSZ2=ISYMSZ(ISYTYP(1,13+ISPIN),ISYTYP(2,13+ISPIN))
       ISCRSZ=NVSIZ+NOSIZ+NOCC*NVRT
       I000=ISTART
       I010=I000+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I020=I010+MAX(NTOTSZ,NTOTSZ2)*IINTFP
       I030=I020+ISCRSZ
       CALL GETALL(ICORE(I000),NTOTSZ,1,39+ISPIN)
       CALL SST03I(ICORE(I000),ICORE(I010),NTOTSZ,NTOTSZ2,ICORE(I020),
     &             SPCASE)
       CALL GETALL(ICORE(I000),NTOTSZ2,1,113+ISPIN)
       CALL SAXPY(NTOTSZ2,ONE,ICORE(I000),1,ICORE(I010),1)
       CALL PUTALL(ICORE(I010),NTOTSZ2,1,113+ISPIN)
      ELSEIF(ISPIN.EQ.3)THEN
C
C  SPIN CASES  ABAB, BABA, ABBA AND BAAB.
C
       MAXSIZ=0
C
C               Ab           EA
C   SOLVE FOR  Q   =  SUM  T    H(MEjb)  [PART 1 OF Q(ABAB)]
C               Ij     M,E   IM
C
        LISTT=34
        LISTW=56
        CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &                'TxW',MAXSIZ,INCREM,42,.TRUE.,
     &                ICORE(IOFFT1(1)),ICORE(IOFFT1(1)),      
     &                POP(1,1),POP(1,1),VRT(1,1),VRT(1,1),
     &                TWO,'AAAA',IUHF)
C
C               Ab         Ae
C   SOLVE FOR  Q   = SUM  T   *[-H(mejb)] [PART 2 OF Q(ABAB)]
C               Ij   m,e   Im
C
C USE AA LIST FOR RHF REFERENCE SINCE IT IS IDENTICAL WITH THE (NOT STOR
C  BB LIST.
C
       LISTW=54+IUHF
       LISTT=37
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'TxW',MAXSIZ,1,42,.FALSE.,
     &               ICORE(IOFFT1(1)),ICORE(IOFFT1(1)),
     &               POP(1,1),POP(1,1),VRT(1,1),VRT(1,1),
     &               TWO,'AAAA',IUHF)
C
C NEXT TWO TERMS DO NOT HAVE TO BE COMPUTED FOR RHF REFERENCE SINCE
C  THEY ARE SIMPLY THE TRANSPOSE OF THE Q(ABAB) EVALUATED ABOVE.
C
      IF(IUHF.EQ.0)GOTO 3001
C
C               aB         aE
C   SOLVE FOR  Q   = SUM  T   *[-H(MEJB)] [PART 1 OF Q(BABA)].
C               iJ   M,E   iM
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=36
       LISTW=54
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'WxT',MAXSIZ,1,42,.FALSE.,
     &               ICORE(IOFFT1(2)),ICORE(IOFFT1(2)),
     &               POP(1,2),POP(1,2),VRT(1,2),VRT(1,2),
     &               TWO,'BBBB',IUHF)
C
C               aB           ea
C   SOLVE FOR  Q   =   SUM  T    H(meJB) [PART 2 OF Q(BABA)].
C               iJ     m,e   im
C
C NEED TO CHANGE LISTW FOR W(mbej) INTERMEDIATE PICKUP!
C
       LISTT=35
       LISTW=57
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'WxT',MAXSIZ,1,42,.TRUE.,
     &               ICORE(IOFFT1(2)),ICORE(IOFFT1(2)),      
     &               POP(1,2),POP(1,2),VRT(1,2),VRT(1,2),
     &               TWO,'BBBB',IUHF)
C
C RHF REENTRY POINT - FORM Q(ABAB) + Q(ABAB) (t) IF RHF.
C
3001   IF(IUHF.EQ.0)THEN
        I000=ISTART
        I010=I000+MAXSIZ*IINTFP
        CALL RHFQ(ICORE(ISTART),MAXSIZ,42)
       ELSE
        CONTINUE
       ENDIF
C
C               Ab         eA
C   SOLVE FOR  Q   = SUM  T   *[- H(MeJb)]  [Q(BAAB)]
C               iJ   M,e   iM
C
       LISTW=58+IUHF
       LISTT=39
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'TxW',MAXSIZ,INCREM,43,.TRUE.,
     &               ICORE(IOFFT1(1)),ICORE(IOFFT1(2)),      
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &               TWO,'BABA',IUHF)
C
C SKIP NEXT TERM FOR RHF
C
      IF(IUHF.EQ.0)GOTO 3002
C
C               aB         Ea
C   SOLVE FOR  Q   = SUM  T    *[-H(mEjB]  [Q(ABBA)]
C               Ij   m,E   Im
C
C NEED TO CHANGE LISTW FOR W INTERMEDIATE PICKUP!
C
       LISTT=38
       LISTW=58
       CALL HTRNGDRV(ICORE(ISTART),MXCOR,LISTW,LISTT,
     &               'WxT',MAXSIZ,1,43,.TRUE.,
     &               ICORE(IOFFT1(2)),ICORE(IOFFT1(1)),      
     &               POP(1,2),POP(1,1),VRT(1,2),VRT(1,1),
     &               TWO,'ABAB',IUHF)
C
C RHF REENTRY POINT #2. FORM Q(BAAB) + Q(BAAB)(t).
C
3002   IF(IUHF.EQ.0)THEN
        I000=ISTART
        I010=I000+MAXSIZ*IINTFP
        CALL RHFQ(ICORE(ISTART),MAXSIZ,43)
       ELSE
        CONTINUE
       ENDIF
C
C SUM TOGETHER THE (AI,bj) AND (Aj,bI) RING INCREMENTS AND FORM THE
C FULL Z(AB) PIECE.  
C
C NOTE THAT THE INCREMENTS COMPUTED THUS FAR ARE NEGATED IN THE ROUTINE SUMRNG
C
       NSZAB=ISYMSZ(ISYTYP(1,42),ISYTYP(2,42))
       ISCRSZ=(NOCCA+NOCCB)*(NVRTA+NVRTB)
       I000=ISTART
       I010=I000+NSZAB*IINTFP
       I020=I010+NSZAB*IINTFP
       I030=I020+ISCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('H4T2ALL',I030,MAXCOR)
       CALL SUMRNG(ICORE(I000),ICORE(I010),ICORE(I020),NSZAB,NSZAB,
     &             ISCRSZ)
C
C NOW CONVERT THESE TO Ab-Ij INCREMENTS AND AUGMENT THE T2 INCREMENT LIST.
C
       I030=I020+ISCRSZ
       IF(I030.GT.MAXCOR)CALL INSMEM('H4X2ALL',I030,MAXCOR)
       CALL SST02I(ICORE(I000),ICORE(I010),NSZAB,NSZAB,ICORE(I020),
     &             'AABB')
C
C NOW SUM THIS WITH THE EXISTING AB INCREMENT AND OVERWRITE IT.
C
       CALL SUMSYM(ICORE(I010),ICORE(I000),NSZAB,116)
      ENDIF
      RETURN
      END
