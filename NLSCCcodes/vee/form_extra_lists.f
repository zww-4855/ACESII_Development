










      SUBROUTINE FORM_EXTRA_LISTS(WORK,MAXCOR,IUHF) 

      IMPLICIT INTEGER (A-Z)
      DIMENSION WORK(MAXCOR)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
C DRIVER FOR POINTER CREATION FOR LISTS 254-259 (RING INTERMEDIATES)
C

      DO ISPIN = 1, IUHF+1
         CALL UPDMOI(1,NF1(ISPIN),6+ISPIN,91,0,0)
         CALL UPDMOI(1,NF2(ISPIN),6+ISPIN,92,0,0)
         CALL ACES_LIST_MEMSET(6+ISPIN,91,0)
         CALL ACES_LIST_MEMSET(6+ISPIN,92,0)
      END DO
      CALL UPDMOI(1,NF1(1),9,91,0,0)
      CALL UPDMOI(1,NF2(1),9,92,0,0)
      CALL ACES_LIST_MEMSET(9,91,0)
      CALL ACES_LIST_MEMSET(9,92,0)

      RETURN
      END 

