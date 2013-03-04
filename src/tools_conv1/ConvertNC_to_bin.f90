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

      real*4 :: TIME_4
      real*8, allocatable :: XPloc(:), YPloc(:)
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
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -1-', TRIM(CHRERR)
        STOP
      ENDIF

      ISTAT = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -2-', TRIM(CHRERR)
        STOP
      ENDIF

      ISTAT = NF90_INQUIRE_VARIABLE(ncid, ITIME_ID, dimids = dimids)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -3-', TRIM(CHRERR)
        STOP
      ENDIF

      ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = nbTime)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -4-', TRIM(CHRERR)
        STOP
      ENDIF

      ALLOCATE(LTimeDay(nbTime))
      ISTAT = NF90_GET_VAR(ncid, itime_id, LTimeDay)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -5-', TRIM(CHRERR)
        STOP
      ENDIF

      ISTAT = nf90_inq_dimid(ncid, 'np_global', idim_id)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -6-', TRIM(CHRERR)
        STOP
      ENDIF


      ISTAT = nf90_inquire_dimension(ncid, idim_id, len = np_global)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -7-', TRIM(CHRERR)
        STOP
      ENDIF

      ISTAT = NF90_CLOSE(ncid)
      IF (ISTAT .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(ISTAT)
        Print *, 'ConvertNC_to_bin -8-', TRIM(CHRERR)
        STOP
      ENDIF
!
!  Build XP, YP global
!
      ALLOCATE(ListMNP(nbProc))
      ALLOCATE(XP(np_global))
      ALLOCATE(YP(np_global))
      ALLOCATE(TheMax(np_global))
      TheMax=0
      DO iProc=1,nbProc
        WRITE (FILE_NAME,10) TRIM(HisPrefix),iProc-1
        ISTAT = NF90_OPEN(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -9-', TRIM(CHRERR)
          STOP
        ENDIF

        ISTAT = nf90_inq_dimid(ncid, 'mnp', idim_id)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -10-', TRIM(CHRERR)
          STOP
        ENDIF

        ISTAT = nf90_inquire_dimension(ncid, idim_id, len = eMNP)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -11-', TRIM(CHRERR)
          STOP
        ENDIF
        ListMNP(iProc)=eMNP
        Print *, 'iProc=', iProc, ' eMNP=', eMNP

        ISTAT = nf90_inq_varid(ncid, 'x', var_id)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -12-', TRIM(CHRERR)
          STOP
        ENDIF

        ALLOCATE(XPloc(eMNP))
        ISTAT = NF90_GET_VAR(ncid, var_id, XPloc)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -13-', TRIM(CHRERR)
          STOP
        ENDIF


        ISTAT = nf90_inq_varid(ncid, 'y', var_id)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -14-', TRIM(CHRERR)
          STOP
        ENDIF

        ALLOCATE(YPloc(eMNP))
        ISTAT = NF90_GET_VAR(ncid, var_id, YPloc)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -15-', TRIM(CHRERR)
          STOP
        ENDIF


        ISTAT = nf90_inq_varid(ncid, 'iplg', var_id)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -16-', TRIM(CHRERR)
          STOP
        ENDIF

        ALLOCATE(iplg(eMNP))
        ISTAT = NF90_GET_VAR(ncid, var_id, iplg)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -17-', TRIM(CHRERR)
          STOP
        ENDIF
        DO IP=1,eMNP
          XP(iplg(IP))=SNGL(XPloc(IP))
          YP(iplg(IP))=SNGL(YPloc(IP))
        END DO
        deallocate(XPloc)
        deallocate(YPloc)
        deallocate(iplg)

        ISTAT = NF90_CLOSE(ncid)
        IF (ISTAT .NE. nf90_noerr) THEN
          CHRERR = nf90_strerror(ISTAT)
          Print *, 'ConvertNC_to_bin -18-', TRIM(CHRERR)
          STOP
        ENDIF

      END DO
      Print *, 'XP, min=', minval(XP), ' max=', maxval(XP)
      Print *, 'YP, min=', minval(YP), ' max=', maxval(YP)
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
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -19-', TRIM(CHRERR)
            STOP
          ENDIF

          ISTAT = nf90_inq_varid(ncid, TRIM(VarName), var_id)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -20-', TRIM(CHRERR)
            STOP
          ENDIF

          eMNP=ListMNP(iProc)

          ISTAT = nf90_inq_varid(ncid, 'iplg', var_id)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -21-', TRIM(CHRERR)
            STOP
          ENDIF

          ALLOCATE(iplg(eMNP))
          ISTAT = NF90_GET_VAR(ncid, var_id, iplg)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -22-', TRIM(CHRERR)
            STOP
          ENDIF

          ISTAT = nf90_inq_varid(ncid, TRIM(VarName), var_id)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -23-', TRIM(CHRERR)
            STOP
          ENDIF

          ISTAT = NF90_INQUIRE_VARIABLE(ncid, var_id, dimids = dimids)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -24-', TRIM(CHRERR)
            STOP
          ENDIF

          ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = len1)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -25-', TRIM(CHRERR)
            STOP
          ENDIF
!          Print *, 'len1=', len1

          ISTAT = nf90_inquire_dimension(ncid, dimids(2), len = len2)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -26-', TRIM(CHRERR)
            STOP
          ENDIF
!          Print *, 'len2=', len2



          ALLOCATE(TheR(eMNP, len))
!          Print *, 'nb1=', nb1, 'len=', len
          ISTAT = NF90_GET_VAR(ncid, var_id, TheR, start=(/1,nb1/), count=(/eMNP,len/))
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -27-', TRIM(CHRERR)
            STOP
          ENDIF
!          Print *, 'TheR min=', minval(TheR), ' max=', maxval(TheR)

          ISTAT = NF90_CLOSE(ncid)
          IF (ISTAT .NE. nf90_noerr) THEN
            CHRERR = nf90_strerror(ISTAT)
            Print *, 'ConvertNC_to_bin -28-', TRIM(CHRERR)
            STOP
          ENDIF

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

