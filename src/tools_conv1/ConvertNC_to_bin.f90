      PROGRAM MERGE_WWM_NC_FILE
      USE NETCDF
      IMPLICIT NONE
      CHARACTER(LEN = 100) :: HisPrefix
      CHARACTER(LEN = 100) :: VarName
      CHARACTER(LEN = 100) :: BinFile
      character(len = 256) :: FILE_NAME
      character(len = 728) :: CHRERR
      integer, dimension(nf90_max_var_dims) :: dimids

      integer nbArg, nbProc, iProc, FHNDL
      integer sizBlock
      integer IP, len, nb1, nb2, eMNP
      integer ISTAT, idim_id, ncid, var_id, itime_id
      integer iTime, nbTime, iEnt
      integer np_global, len1, len2
      REAL*8, PARAMETER   :: DAY2SEC  = 86400
      character (len = *), parameter :: CallFct="MERGE_WWM_NC_FILE"

      real*4 :: TIME_4
      real*4, allocatable :: XP(:), YP(:)
      real*4, allocatable :: TheMat(:,:)
      real*4, allocatable :: TheR(:,:)
      real*4, allocatable :: TheMax(:)
      real*4, allocatable :: TheWrite(:)
      real*8, allocatable :: LTimeDay(:)
      integer, allocatable :: ListMNP(:)
      integer, allocatable :: iplg(:)


      nbArg=command_argument_count()
      IF (nbArg == 0) THEN
        PRINT *, 'Reading from multiple netcdf file'
        Print *, 'and create a simple bin file'
        Print *, 'Use as'
        Print *, 'ConvertNC_to_bin [HisPrefix] [Var] [BinFile]'
        Print *, ''
        Print *, 'If file is WWM_output_0001.nc then'
        Print *, 'HisPrefix=WWM_output_'
        Print *, 'Var is HS, TP, etc.'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1, HisPrefix)
      CALL GET_COMMAND_ARGUMENT(2, VarName)
      CALL GET_COMMAND_ARGUMENT(3, BinFile)
      FHNDL=10

      OPEN(FHNDL, FILE  = TRIM(BinFile), FORM = 'UNFORMATTED')
      CALL DETERMINE_NUMBER_PROC(HisPrefix, nbProc)
      Print *, 'nbProc=', nbProc
      WRITE (FILE_NAME,10) TRIM(HisPrefix),0
!
! Determine nbTime and np_global
! Also read ocean_time_day
!
      Print *, 'FILE_NAME=', TRIM(FILE_NAME)
      ISTAT = NF90_OPEN(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ISTAT = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

      ISTAT = NF90_INQUIRE_VARIABLE(ncid, ITIME_ID, dimids = dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

      ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = nbTime)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

      ALLOCATE(LTimeDay(nbTime))
      ISTAT = NF90_GET_VAR(ncid, itime_id, LTimeDay)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

      ISTAT = nf90_inq_dimid(ncid, 'np_global', idim_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

      ISTAT = nf90_inquire_dimension(ncid, idim_id, len = np_global)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

      ISTAT = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
!
!  Build XP, YP global
!
      ALLOCATE(ListMNP(nbProc))
      ALLOCATE(XP(np_global))
      ALLOCATE(YP(np_global))
      ALLOCATE(TheMax(np_global))
      WRITE (FILE_NAME,10) TRIM(HisPrefix),0

      ISTAT = NF90_OPEN(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

      ISTAT = nf90_inq_varid(ncid, 'x', var_id)
      IF (ISTAT /= 0) THEN
        ISTAT = nf90_inq_varid(ncid, 'lon', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)
      END IF

      ISTAT = NF90_GET_VAR(ncid, var_id, XP)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

      ISTAT = nf90_inq_varid(ncid, 'y', var_id)
      IF (ISTAT /= 0) THEN
        ISTAT = nf90_inq_varid(ncid, 'lat', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)
      END IF

      ISTAT = NF90_GET_VAR(ncid, var_id, YP)
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

      ISTAT = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)
      Print *, 'XP, min=', minval(XP), ' max=', maxval(XP)
      Print *, 'YP, min=', minval(YP), ' max=', maxval(YP)

      TheMax=0
      DO iProc=1,nbProc
        WRITE (FILE_NAME,10) TRIM(HisPrefix),iProc-1
        ISTAT = NF90_OPEN(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

        ISTAT = nf90_inq_dimid(ncid, 'mnp', idim_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

        ISTAT = nf90_inquire_dimension(ncid, idim_id, len = eMNP)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)
        ListMNP(iProc)=eMNP
        Print *, 'iProc=', iProc, ' eMNP=', eMNP

        ISTAT = NF90_CLOSE(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)

      END DO
!
!  Now loop
!
      allocate(TheWrite(np_global))
      sizBlock=60
      iEnt=1
      DO
        nb1=1+(iEnt-1)*sizBlock
        IF (nb1 .gt. nbTime) THEN
          EXIT
        END IF
        nb2=iEnt*sizBlock
        IF (nb2 .gt. nbTime) THEN
          nb2=nbTime
        END IF
        len=nb2+1-nb1
        ALLOCATE(TheMat(len,np_global))
        DO iProc=1,nbProc
          Print *, 'iEnt=', iEnt, ' iProc=', iProc
          WRITE (FILE_NAME,10) TRIM(HisPrefix),iProc-1
          ISTAT = NF90_OPEN(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 19, ISTAT)

          ISTAT = nf90_inq_varid(ncid, TRIM(VarName), var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 20, ISTAT)

          eMNP=ListMNP(iProc)

          ISTAT = nf90_inq_varid(ncid, 'iplg', var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 21, ISTAT)

          ALLOCATE(iplg(eMNP))
          ISTAT = NF90_GET_VAR(ncid, var_id, iplg)
          CALL GENERIC_NETCDF_ERROR(CallFct, 22, ISTAT)

          ISTAT = nf90_inq_varid(ncid, TRIM(VarName), var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 23, ISTAT)

          ISTAT = NF90_INQUIRE_VARIABLE(ncid, var_id, dimids = dimids)
          CALL GENERIC_NETCDF_ERROR(CallFct, 24, ISTAT)

          ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = len1)
          CALL GENERIC_NETCDF_ERROR(CallFct, 25, ISTAT)

          ISTAT = nf90_inquire_dimension(ncid, dimids(2), len = len2)
          CALL GENERIC_NETCDF_ERROR(CallFct, 26, ISTAT)

          ALLOCATE(TheR(eMNP, len))
!          Print *, 'nb1=', nb1, 'len=', len
          ISTAT = NF90_GET_VAR(ncid, var_id, TheR, start=(/1,nb1/), count=(/eMNP,len/))
          CALL GENERIC_NETCDF_ERROR(CallFct, 27, ISTAT)

          ISTAT = NF90_CLOSE(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 28, ISTAT)

          DO IP=1,eMNP
            TheMat(:,iplg(IP))=TheR(IP,:)
          END DO
          DEALLOCATE(iplg)
          DEALLOCATE(TheR)
        END DO
        DO iTime=nb1,nb2
          Print *, 'iTime=', iTime, '/', nbTime
          TIME_4=SNGL(DAY2SEC*(LTimeDay(iTime) - LTimeDay(1)))
          TheWrite=TheMat(iTime+1-nb1,:)
          Print *, 'min=', minval(TheWrite), ' max=', maxval(TheWrite), 'avg=', sum(TheWrite)/REAL(np_global)
          WRITE(FHNDL) TIME_4
          WRITE(FHNDL) (XP(IP), YP(IP), TheMat(iTime+1-nb1,IP), IP=1,np_global)
          DO IP=1,np_global
            IF (TheMat(iTime+1-nb1,IP) .gt. TheMax(IP)) THEN
              TheMax(IP)=TheMat(iTime+1-nb1,IP)
            END IF
          END DO
        END DO
        DEALLOCATE(TheMat)
        iEnt=iEnt+1
      END DO
      WRITE(FHNDL) TIME_4
      WRITE(FHNDL) (XP(IP), YP(IP), TheMax(IP), IP=1,np_global)

      close(FHNDL)
      deallocate(ListMNP)
      deallocate(TheWrite)
      deallocate(XP, YP)
  10  FORMAT (a,i4.4,'.nc')
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE DETERMINE_NUMBER_PROC(HisPrefix, nbProc)
      IMPLICIT NONE
      character(len=*), intent(in) :: HisPrefix
      integer, intent(out) :: nbProc
      character(len = 256) :: FILE_NAME
      integer :: iRankTest
      logical :: test
      iRankTest=1
      DO
        WRITE (FILE_NAME,10) TRIM(HisPrefix),iRankTest-1
        INQUIRE(FILE=TRIM(FILE_NAME), EXIST=test)
        IF (test.eqv..false.) THEN
          nbProc=iRankTest-1
          EXIT
        END IF
        iRankTest=iRankTest+1
      END DO
  10  FORMAT (a,i4.4,'.nc')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE NETCDF
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        Print *, TRIM(CallFct), ' - ', idx, ' - ', TRIM(CHRERR)
        STOP
      ENDIF
      END SUBROUTINE
