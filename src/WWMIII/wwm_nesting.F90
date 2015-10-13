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
      USE WWM_HOTFILE_MOD
      IMPLICIT NONE
      integer, intent(in) :: iGrid
      character(len=140) FILERET
      real(rkind), allocatable :: ACwrite(:,:,:), VAR_ONEDwrite(:,:)
      real(rkind), allocatable :: ACsend(:,:,:), VAR_ONEDsend(:,:)
      integer, allocatable :: ListStatus(:)
      real(rkind) :: eVect(nbOned)
      !
      integer eInt(1)
      integer np_write, ne_write
      integer MULTIPLEOUT_W
      logical GRIDWRITE_W, IOBPD_HISTORY_W, WriteOutputProcess
      integer nbMatch, nbTime
      integer IP, IE, I, idx, IP2, iProc, nbMatchLoc
      real(rkind) eW
      integer, allocatable :: ListMatch(:)
      real(rkind) eTimeDay
      integer POS
      np_write=ListNestInfo(iGrid) % eGrid % np_total
      ne_write=ListNestInfo(iGrid) % eGrid % ne_total
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        FILERET = ListPrefix(iGrid) // '_hotfile.nc'
        MULTIPLEOUT_W = .FALSE.
        GRIDWRITE_W = .FALSE.
        IOBPD_HISTORY_W = .FALSE.
        WriteOutputProcess = .TRUE.
        nbTime=-1
        CALL WRITE_HOTFILE_PART_1(FILERET, nbTime, MULTIPLEOUT_W, GRIDWRITE_W, IOBPD_HISTORY_W, np_write, ne_write)
        CALL WRITE_NETCDF_HEADERS_2(FILERET, MULTIPLEOUT_W, WriteOutputProcess, GRIDWRITE_W, np_write, ne_write)
        allocate(ACwrite(MSC,MDC,np_write), VAR_ONEDwrite(nbOned, np_write), ListStatus(np_write), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
        ListStatus=0
#ifdef MPI_PARALL_GRID
      END IF
#endif
      nbMatch=0
      DO IP=1,np_write
        IE=ListNestInfo(iGrid) % HOT_IE(IP)
        IF (IE .gt. 0) THEN
          nbMatch = nbMatch + 1
        END IF
      END DO
      allocate(Listmatch(nbMatch), ACsend(MSC,MDC,nbMatch), VAR_ONEDsend(nbOned, nbMatch), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
      ACsend=0
      VAR_ONEDsend=0
      idx=0
      DO IP=1,np_write
        IE=ListNestInfo(iGrid) % HOT_IE(IP)
        IF (IE .gt. 0) THEN
          idx=idx+1
          ListMatch(idx)=IP
          DO I=1,3
            eW=ListNestInfo(iGrid) % HOT_W(I,IP)
            IP2=INE(I,IE)
            eVect(1)=WINDXY(IP2,1)
            eVect(2)=WINDXY(IP2,2)
            eVect(3)=PRESSURE(IP2)
            eVect(4)=DVWIND(IP2,1)
            eVect(5)=DVWIND(IP2,2)
            eVect(6)=CURTXY(IP2,1)
            eVect(7)=CURTXY(IP2,2)
            eVect(8)=DVCURT(IP2,1)
            eVect(9)=DVCURT(IP2,2)
            eVect(10)=DDEP(IP2,1)
            eVect(11)=DDEP(IP2,2)
            eVect(12)=DCUX(IP2,1)
            eVect(13)=DCUX(IP2,2)
            eVect(14)=DCUY(IP2,1)
            eVect(15)=DCUY(IP2,2)
            eVect(16)=WATLEV(IP2)
            eVect(17)=WATLEVOLD(IP2)
            eVect(18)=DVWALV(IP2)
            eVect(19)=WLDEP(IP2)
            eVect(20)=DEPDT(IP2)
            eVect(21)=QBLOCAL(IP2)
            eVect(22)=DISSIPATION(IP2)
            eVect(23)=AIRMOMENTUM(IP2)
            eVect(24)=UFRIC(IP2)
            eVect(25)=ALPHA_CH(IP2)
            eVect(26)=TAUW(IP2)
            eVect(27)=TAUTOT(IP2)
            eVect(28)=TAUWX(IP2)
            eVect(29)=TAUWY(IP2)
            eVect(30)=TAUHF(IP2)
            eVect(31)=Z0(IP2)
            eVect(32)=CD(IP2)
            eVect(33)=USTDIR(IP2)
            eVect(34)=RSXX(IP2)
            eVect(35)=RSXY(IP2)
            eVect(36)=RSYY(IP2)
            eVect(37)=FORCEXY(IP2,1)
            eVect(38)=FORCEXY(IP2,2)
            ACsend(:,:,idx) = ACsend(:,:,idx) + eW * AC2(:,:,IP2)
            VAR_ONEDsend(:,idx) = VAR_ONEDsend(:,idx) + eW * eVect
          END DO
        END IF
      END DO
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        DO idx=1,nbMatch
          IP=ListMatch(idx)
          ACwrite(:,:,IP)     = ACsend(:,:,idx)
          VAR_ONEDwrite(:,IP) = VAR_ONEDsend(:,idx)
          ListStatus(IP)=1
        END DO
        deallocate(ACsend, VAR_ONEDsend, ListMatch)
        DO iProc=2,nproc
          CALL MPI_RECV(eInt, 1, itype, iProc-1, 2401, comm, istatus, ierr)
          nbMatchLoc=eInt(1)
          allocate(ListMatch(nbMatchLoc), ACsend(MSC,MDC,nbMatchLoc), VAR_ONEDsend(nbOned, nbMatchLoc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_nesting, allocate error 1')
          CALL MPI_RECV(ListMatch, nbMatchLoc, itype, iProc-1, 2402, comm, istatus, ierr)
          CALL MPI_RECV(ACsend, MSC*MDC*nbMatchLoc, itype, iProc-1, 2403, comm, istatus, ierr)
          CALL MPI_RECV(VAR_ONEDsend, nbOned*nbMatchLoc, itype, iProc-1, 2404, comm, istatus, ierr)
          DO idx=1,nbMatchLoc
            IP=ListMatch(idx)
            ACwrite(:,:,IP) = ACsend(:,:,idx)
            VAR_ONEDwrite(:,IP) = VAR_ONEDsend(:,idx)
            ListStatus(IP)=1
          END DO
          deallocate(ListMatch, ACsend, VAR_ONEDsend)
        END DO
      ELSE
        eInt(1)=nbMatch
        CALL MPI_SEND(eInt, 1, itype, 0, 2401, comm, ierr)
        CALL MPI_SEND(ListMatch, nbMatch, itype, 0, 2402, comm, ierr)
        CALL MPI_SEND(ACsend, MSC*MDC*nbMatch, rtype, 0, 2403, comm, ierr)
        CALL MPI_SEND(VAR_ONEDsend, nbOned*nbMatch, rtype, 0, 2404, comm, ierr)
      END IF
#else
      ACwrite = ACsend
      VAR_ONEDwrite = VAR_ONEDsend
#endif      
#ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
#endif
        eTimeDay=ListNestInfo(iGrid) % eTime % BMJD
        POS=1
        CALL WRITE_HOTFILE_PART_2(FILERET, eTimeDay, POS, np_write, ACwrite, VAR_ONEDwrite)
#ifdef MPI_PARALL_GRID
      END IF
#endif
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
