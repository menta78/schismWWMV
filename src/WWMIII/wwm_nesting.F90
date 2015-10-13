#include "wwm_functions.h"
#ifdef NCDF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NESTING
      USE DATAPOOL
      IMPLICIT NONE
      real(rkind) :: XYTMP(2,MNP)
      real(rkind) :: eWI(3)
      integer iGrid, np_total_loc, IWBMNP_loc
      type(FILEDEF) eGRD, eBND
      integer IP, idx
      REAL(rkind) eX, eY
      integer eElt, NI(3), eIdx
      !
      IF (L_NESTING .eqv. .FALSE.) THEN
        RETURN
      END IF
      !
      ! First reading the grids
      !
      ALLOCATE(ListNestInfo(NB_GRID_NEST), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
      DO iGrid=1,NB_GRID_NEST
        IGRIDTYPE = ListIGRIDTYPE(iGrid)
        eGRD % FNAME = ListFILEGRID(iGrid)
        eGRD % FHNDL = 24037
        CALL READ_SPATIAL_GRID_TOTAL_KERNEL(ListNestInfo(iGrid) % eGrid, DIMMODE, LVAR1D, LSPHE, eGRD, IGRIDTYPE)
        !
!        np_total_loc = ListNestInfo(iGrid) % eGrid % np_total
!        allocate(ListNestInfo(iGrid) % NodeBelonging(np_total_loc), stat=istat)
!        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
        !
        allocate(ListNestInfo(iGrid) % IOBPtotal(np_total_loc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 3')
        eBND % FNAME = ListFILEBOUND(iGrid)
        eBND % FHNDL = 24977
        CALL SINGLE_READ_IOBP_TOTAL(IOBPtotal, IGRIDTYPE, eBND, np_total_loc)
        !
        IWBMNP_loc=0
        DO IP=1,np_total_loc
          IF ((IOBPtotal(IP) .eq. 2) .or. (IOBPtotal(IP) .eq. 3)) THEN
            IWBMNP_loc = IWBMNP_loc + 1
          END IF
        END DO
        ListNestInfo(iGrid) % IWBMNP = IWBMNP_loc
        !
        allocate(ListNestInfo(iGrid) % IWBNDLC(IWBMNP_loc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
        idx=0
        DO IP=1,np_total_loc
          IF ((IOBPtotal(IP) .eq. 2) .or. (IOBPtotal(IP) .eq. 3)) THEN
            idx=idx+1
            ListNestInfo(iGrid) % IWBNDLC(idx) = IP
          END IF
        END DO
        !
        
      END DO
      !
      ! Now we find the IE and weights for interpolation
      !
      XYTMP(1,:) = XP
      XYTMP(2,:) = YP
      DO iGrid=1,NB_GRID_NEST
        np_total_loc = ListNestInfo(iGrid) % eGrid % np_total
        IWBMNP_loc = ListNestInfo(iGrid) % IWBMNP
        IF (L_HOTFILE) THEN
          allocate(ListNestInfo(iGrid) % HOT_IE(np_total_loc),ListNestInfo(iGrid) % HOT_W(3, np_total_loc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
          DO IP=1,np_total_loc
            eX=ListNestInfo(iGrid) % eGrid % XPtotal(IP)
            eY=ListNestInfo(iGrid) % eGrid % YPtotal(IP)
            CALL FIND_ELE(MNE,MNP,INE,XYTMP,eX, eY, eElt)
            ListNestInfo(iGrid) % HOT_IE(IP) = eElt
            IF (eElt .gt. 0) THEN
              NI                 = INE(:,eElt)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), eX, eY, eWI)
              ListNestInfo(iGrid) % HOT_W(:, IP) = eWI
            END IF
          END DO
        END IF
        IF (L_BOUC_PARAM .or. L_BOUC_SPEC) THEN
          allocate(ListNestInfo(iGrid) % BOUC_IE(IWBMNP_loc),ListNestInfo(iGrid) % BOUC_W(3, IWBMNP_loc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 2')
          DO eIdx=1,IWBMNP_loc
            IP=ListNestInfo(iGrid) % IWBNDLC(eIdx)
            eX=ListNestInfo(iGrid) % eGrid % XPtotal(IP)
            eY=ListNestInfo(iGrid) % eGrid % YPtotal(IP)
            CALL FIND_ELE(MNE, MNP, INE, XYTMP, eX, eY, eElt)
            ListNestInfo(iGrid) % BOUC_IE(IP) = eElt
            IF (eElt .gt. 0) THEN
              NI                 = INE(:,eElt)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), eX, eY, eWI)
              ListNestInfo(iGrid) % BOUC_W(:, IP) = eWI
            END IF
          END DO
        END IF
      END DO
      !
      ! Now we timings needed by the model
      !
      DO iGRid=1,NB_GRID_NEST
        ListNestInfo(iGrid) % eTime % BEGT = ListBEGTC(iGrid)
        ListNestInfo(iGrid) % eTime % DELT = ListDELTC(iGrid)
        ListNestInfo(iGrid) % eTime % UNIT = ListUNITC(iGrid)
        ListNestInfo(iGrid) % eTime % ENDT = ListENDTC(iGrid)
        CALL CT2MJD(ListNestInfo(iGrid) % eTime % BEGT, ListNestInfo(iGrid) % eTime % BMJD)
        CALL CT2MJD(ListNestInfo(iGrid) % eTime % ENDT, ListNestInfo(iGrid) % eTime % EMJD)
        CALL CU2SEC(ListNestInfo(iGrid) % eTime % UNIT, ListNestInfo(iGrid) % eTime % DELT)

        ListNestInfo(iGrid) % eTime % TOTL = (ListNestInfo(iGrid) % eTime % EMJD - ListNestInfo(iGrid) % eTime % BMJD) * DAY2SEC
        ListNestInfo(iGrid) % eTime % ISTP = NINT(ListNestInfo(iGrid) % eTime % TOTL / ListNestInfo(iGrid) % eTime % DELT) + 1
        ListNestInfo(iGrid) % eTime % TMJD = ListNestInfo(iGrid) % eTime % BMJD
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NESTING_OUTPUT_HOTFILE(iGrid)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: iGrid
      !
      
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NESTING_BOUNDARY_CONDITION(iGrid)
      USE DATAPOOL
      IMPLICIT NONE
      integer, intent(in) :: iGrid
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DO_NESTING_OPERATION
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, SAVE :: INITDONE = .FALSE.
      integer iGrid
      real(rkind) DeltaTimeDiff
      !
      ! First the init
      !
      IF (INITDONE .eqv. .FALSE.) THEN
        CALL INIT_NESTING
      END IF
      INITDONE = .TRUE.
      !
      ! Then the HOTFILE
      !
      DO iGrid=1,NB_GRID_NEST
        IF ((MAIN%TMJD .GE. ListNestInfo(iGrid) % eTime % TMJD - 1.E-8) .AND. (MAIN%TMJD .LE. ListNestInfo(iGrid) % eTime % EMJD)) THEN
          DeltaTimeDiff = abs(MAIN % TMJD - ListNestInfo(iGrid) % eTime % BMJD)
          IF (L_HOTFILE .and. DeltaTimeDiff .le. 1.e-8) THEN
            CALL NESTING_OUTPUT_HOTFILE(iGrid)
          END IF
          IF (L_BOUC_PARAM .or. L_BOUC_SPEC) THEN
            CALL NESTING_BOUNDARY_CONDITION(iGrid)
          END IF
          ListNestInfo(iGrid) % eTime % TMJD = ListNestInfo(iGrid) % eTime % TMJD + ListNestInfo(iGrid) % eTime % DELT*SEC2DAY
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#endif
