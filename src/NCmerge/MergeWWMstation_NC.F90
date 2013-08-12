      PROGRAM MERGE_WWM_HISTORY_NC
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN = 100) :: HisPrefix
      CHARACTER(LEN = 100) :: VarName
      CHARACTER(LEN = 100) :: BinFile
      character(len = 1000) :: FILE_NAME_SPLIT
      character(len = 1000) :: FILE_NAME_MERGE
      character(len = 728) :: CHRERR
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      integer, dimension(nf90_max_var_dims) :: dimids

      integer nbArg, nbProc, iProc, FHNDL
      integer sizBlock
      integer IP, len, nb1, nb2, eMNP
      integer iret, idim_id, ncid, ncidB, var_id, itime_id
      integer nnode_dims, ntime_dims, istat
      integer iTime, nbTime, iEnt, IVAR
      integer len1, len2
      character (len = *), parameter :: CallFct="MergeWWMstation_NC"

      integer, allocatable :: ListORIGIN(:), ISum(:), IFound(:), ISMAX(:)
      integer ifile, I, nproc
      integer :: eInt(1)
     
      logical test
      real(rkind), allocatable :: LTimeDay(:)
      real(rkind) :: eTimeDay
      real(rkind), allocatable :: OUTPAR_STATIONS_R(:,:), OUTPAR_STATIONS_W(:,:)
      real(rkind), allocatable :: WK_STATIONS_R(:,:), WK_STATIONS_W(:,:)
      real(rkind), allocatable :: AC_STATIONS_R(:,:,:), AC_STATIONS_W(:,:,:)
      real(rkind), allocatable :: ACOUT_1D_STATIONS_R(:,:,:), ACOUT_1D_STATIONS_W(:,:,:)
      real(rkind), allocatable :: ACOUT_2D_STATIONS_R(:,:,:), ACOUT_2D_STATIONS_W(:,:,:)
      nbArg=command_argument_count()
      IF (nbArg .ge. 2) THEN
        PRINT *, 'MergeWWMstation_NC merges netcdf station files'
        Print *, 'It is supposed to be used just like WWM'
        Print *, ''
        Print *, 'i.e. MergeWWMhistory_NC'
        Print *, 'or   MergeWWMhistory_NC wwminput.nml'
        Print *, 'with the wwminput the same as the one of the run'
        STOP
      END IF
      IF (nbArg.eq.0) THEN
        INP%FNAME  = 'wwminput.nml'
      ELSE
        CALL GET_COMMAND_ARGUMENT(1, INP%FNAME)
      ENDIF
!
! Determine nbTime and np_global
! Also read ocean_time_day
!

!      INP%FHNDL=440
!      open(INP%FHNDL, FILE=TRIM(INP%FNAME))
!      GRD%FHNDL=230
!      STAT%FHNDL=235
!      open(STAT%FHNDL, FILE='wwmstat.out', status='unknown')
      CALL INIT_FILE_HANDLES()
      CALL READ_WWMINPUT
!      CALL INIT_ARRAYS
      CLOSE(INP%FHNDL)
      open(GRD%FHNDL, FILE=TRIM(GRD%FNAME))
!
      Print *, 'GRD%FNAME=', TRIM(GRD%FNAME)

      ALLOCATE(XP(MNP))
      ALLOCATE(YP(MNP))
      ALLOCATE(DEP(MNP))
      ALLOCATE(INVSPHTRANS(MNP,2))
      ALLOCATE(INE(3,MNE))
      ALLOCATE(IEN(6,MNE))
      ALLOCATE(TRIA(MNE))
      np_total=MNP
      ne_total=MNE
!
      CALL READ_SPATIAL_GRID
      CALL SPATIAL_GRID
      close(GRD%FHNDL)
      MULTIPLEOUT_STAT=0
      CALL INIT_SPECTRAL_GRID()
!
! Now the allocations
!
      allocate(OUTPAR_STATIONS_R(IOUTS,OUTVARS_COMPLETE), WK_STATIONS_R(IOUTS,MSC), AC_STATIONS_R(IOUTS,MSC,MDC),        OUTPAR_STATIONS_W(IOUTS,OUTVARS_COMPLETE), WK_STATIONS_W(IOUTS,MSC), AC_STATIONS_W(IOUTS,MSC,MDC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 1')
      OUTPAR_STATIONS_W=0
      WK_STATIONS_W=0
      AC_STATIONS_W=0
      IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
        allocate(ACOUT_1D_STATIONS_R(IOUTS, MSC, 3), ACOUT_2D_STATIONS_R(IOUTS, MSC, MDC), ACOUT_1D_STATIONS_W(IOUTS, MSC, 3), ACOUT_2D_STATIONS_W(IOUTS, MSC, MDC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 3')
        ACOUT_1D_STATIONS_W=0
        ACOUT_2D_STATIONS_W=0
      ENDIF
!
! Creating the first FILE_NAME
!
      ifile=1
      iProc=1
      CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
!
! Reading nproc
!
      iret = NF90_OPEN(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)

      iret = nf90_inq_varid(ncid, 'nproc', var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)

      iret = NF90_GET_VAR(ncid, var_id, eInt)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      nproc=eInt(1)

      iret = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
!
! Read the "ifound" array
!
      allocate(ListORIGIN(IOUTS), IFound(IOUTS), ISmax(IOUTS), ISum(IOUTS), stat=istat)
      DO I=1,IOUTS
        STATION(I)%IFOUND = 0
        STATION(I)%ISMAX  = 0
      END DO
      IF (istat/=0) CALL WWM_ABORT('wwm_output, allocate error 3')
      ListORIGIN=0
      ISum=0
      DO iProc=1,nproc
        ifile=1
        CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
        iret = NF90_OPEN(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)

        iret = nf90_inq_varid(ncid, 'ifound', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)

        iret = NF90_GET_VAR(ncid, var_id, IFound)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)

        iret = nf90_inq_varid(ncid, 'ismax', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)

        iret = NF90_GET_VAR(ncid, var_id, ISMAX)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)

        DO I=1,IOUTS
          IF (IFound(I) .eq. 1) THEN
            ListORIGIN(I)=iProc
            ISum(I)=1
            STATION(I)%IFOUND=1
            STATION(I)%ISMAX = ISMAX(I)
          END IF
        END DO
        iret = NF90_CLOSE(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      END DO
      DEALLOCATE(IFound)
!
! Now looping ...
!
      ifile=1
      DO
        iProc=1
        CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
        INQUIRE(FILE=TRIM(FILE_NAME_SPLIT),EXIST=test)
        IF (test .eqv. .FALSE.) THEN
          EXIT
        END IF
        Print *, 'ifile=', ifile
        !
        ! Reading the time
        !
        iret = nf90_open(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        iret = NF90_INQUIRE_VARIABLE(ncid, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
        iret = nf90_inquire_dimension(ncid, dimids(1), len = nbTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
        ALLOCATE(LTimeDay(nbTime))
        iret = NF90_GET_VAR(ncid, itime_id, LTimeDay)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
        !
        ! Now creating the merged file.
        !
        CALL GET_FILE_NAME_MERGE(FILE_NAME_MERGE, ifile)
        Print *, 'FILE_NAME_MERGE=', TRIM(FILE_NAME_MERGE)
        CALL DEFINE_STATION_NC(FILE_NAME_MERGE, MULTIPLEOUT_STAT)


        DO iTime=1,nbTime
          Print *, 'iTime=', iTime, '/', nbTime
          DO iProc=1,nproc
!            Print *, 'iProc=', iProc
            CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
            iret=nf90_open(TRIM(FILE_NAME_SPLIT), nf90_write, ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
            IF (VAROUT_STATION%AC) THEN
              iret=nf90_get_var(ncid,var_id,AC_STATIONS_R, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              DO I=1,IOUTS
                IF (ListORIGIN(I) .eq. iProc) THEN
                  AC_STATIONS_W(I,:,:)=AC_STATIONS_R(I,:,:)
                END IF
              END DO
            END IF
            IF (VAROUT_STATION%WK) THEN
              iret=nf90_inq_varid(ncid, 'WK', var_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
              iret=nf90_get_var(ncid,var_id,WK_STATIONS_R, start = (/1,1,iTime/), count=(/IOUTS,MSC,1/))
              DO I=1,IOUTS
                IF (ListORIGIN(I) .eq. iProc) THEN
                  WK_STATIONS_W(I,:)=WK_STATIONS_R(I,:)
                END IF
              END DO
            END IF
            IF (VAROUT_STATION%ACOUT_1D) THEN
              iret=nf90_inq_varid(ncid, 'ACOUT_1D', var_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
              iret=nf90_get_var(ncid,var_id,ACOUT_1D_STATIONS_R, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,3,1/))
              DO I=1,IOUTS
                IF (ListORIGIN(I) .eq. iProc) THEN
                  ACOUT_1D_STATIONS_W(I,:,:)=ACOUT_1D_STATIONS_R(I,:,:)
                END IF
              END DO
            END IF
            IF (VAROUT_STATION%ACOUT_2D) THEN
              iret=nf90_inq_varid(ncid, 'ACOUT_2D', var_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
              iret=nf90_get_var(ncid,var_id,ACOUT_2D_STATIONS_R, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
            END IF
            DO IVAR=1,OUTVARS_COMPLETE
              IF (VAROUT_STATION%LVAR(IVAR)) THEN
                CALL NAMEVARIABLE(IVAR, eStr, eStrFullName, eStrUnit)
                iret=nf90_inq_varid(ncid, TRIM(eStr), var_id)
                CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
                iret=nf90_get_var(ncid,var_id,OUTPAR_STATIONS_R(:,IVAR), start = (/1, iTime/), count = (/ IOUTS, 1 /))
                CALL GENERIC_NETCDF_ERROR(CallFct, 21, iret)
              END IF
            END DO
            DO I=1,IOUTS
              IF (ListORIGIN(I) .eq. iProc) THEN
                OUTPAR_STATIONS_W(I,:)=OUTPAR_STATIONS_R(I,:)
              END IF
            END DO
            iret=nf90_close(ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 22, iret)
          END DO
          !
          ! Now doing the writing
          !
          Print *, 'FILE_NAME_MERGE=', TRIM(FILE_NAME_MERGE)
          iret=nf90_open(TRIM(FILE_NAME_MERGE), nf90_write, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 23, iret)
          IF (VAROUT_STATION%AC) THEN
            iret=nf90_inq_varid(ncid, 'AC', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 24, iret)
            IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
              iret=nf90_put_var(ncid,var_id,AC_STATIONS_W, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 25, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(AC_STATIONS_W), start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 26, iret)
            ENDIF
          END IF
          IF (VAROUT_STATION%WK) THEN
            iret=nf90_inq_varid(ncid, 'WK', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)
            IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
              iret=nf90_put_var(ncid,var_id,WK_STATIONS_W, start = (/1,1,iTime/), count=(/IOUTS,MSC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(WK_STATIONS_W), start = (/1,1,iTime/), count=(/IOUTS,MSC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 29, iret)
            ENDIF
          END IF
          IF (VAROUT_STATION%ACOUT_1D) THEN
            iret=nf90_inq_varid(ncid, 'ACOUT_1D', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 30, iret)
            IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
              iret=nf90_put_var(ncid,var_id,ACOUT_1D_STATIONS_W, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,3,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(ACOUT_1D_STATIONS_W), start = (/1,1,1,iTime/), count=(/IOUTS,MSC,3,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 32, iret)
            ENDIF
          END IF
          IF (VAROUT_STATION%ACOUT_2D) THEN
            iret=nf90_inq_varid(ncid, 'ACOUT_2D', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 33, iret)
            IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
              iret=nf90_put_var(ncid,var_id,ACOUT_2D_STATIONS_W, start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 34, iret)
            ELSE
              iret=nf90_put_var(ncid,var_id,SNGL(ACOUT_2D_STATIONS_W), start = (/1,1,1,iTime/), count=(/IOUTS,MSC,MDC,1/))
              CALL GENERIC_NETCDF_ERROR(CallFct, 35, iret)
            ENDIF
          END IF
          DO IVAR=1,OUTVARS_COMPLETE
            IF (VAROUT_STATION%LVAR(IVAR)) THEN
              CALL NAMEVARIABLE(IVAR, eStr, eStrFullName, eStrUnit)
              iret=nf90_inq_varid(ncid, TRIM(eStr), var_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 36, iret)
              IF (NF90_RUNTYPE == NF90_OUTTYPE_STAT) THEN
                iret=nf90_put_var(ncid,var_id,OUTPAR_STATIONS_W(:,IVAR), start = (/1, iTime/), count = (/ IOUTS, 1 /))
                CALL GENERIC_NETCDF_ERROR(CallFct, 37, iret)
              ELSE
                iret=nf90_put_var(ncid,var_id,SNGL(OUTPAR_STATIONS_W(:,IVAR)), start = (/1, iTime/), count = (/ IOUTS, 1 /))
                CALL GENERIC_NETCDF_ERROR(CallFct, 38, iret)
              ENDIF
            END IF
          END DO
          eTimeDay=LTimeDay(iTime)
          CALL WRITE_NETCDF_TIME(ncid, iTime, eTimeDay)
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 39, iret)
        END DO
        DEALLOCATE(LTimeDay)
        IF (OUT_STATION%IDEF .le. 0) THEN
          EXIT
        END IF
        ifile=ifile+1
      END DO
      DEALLOCATE(OUTPAR_STATIONS_R, WK_STATIONS_R, AC_STATIONS_R, OUTPAR_STATIONS_W, WK_STATIONS_W, AC_STATIONS_W)
      IF (VAROUT_STATION%ACOUT_1D.or.VAROUT_STATION%ACOUT_2D) THEN
        deallocate(ACOUT_1D_STATIONS_R, ACOUT_2D_STATIONS_R, ACOUT_1D_STATIONS_W, ACOUT_2D_STATIONS_W)
      ENDIF
      DEALLOCATE(ListORIGIN)
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_SPLIT(FILE_NAME, ifile, iProc)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ifile, iProc
      character(len=1000), intent(out) :: FILE_NAME
      character(len =256) :: PRE_FILE_NAME
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(OUT_STATION%FNAME)
      IF (OUT_STATION%IDEF.gt.0) THEN
         WRITE (PRE_FILE_NAME,10) OUT_STATION%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4)
      ELSE
         WRITE (PRE_FILE_NAME,20) OUT_STATION%FNAME(1:LPOS)
  20     FORMAT (a)
      ENDIF
      WRITE (FILE_NAME,40) TRIM(PRE_FILE_NAME),iProc-1
  40  FORMAT (a,'_',i4.4,'.nc')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_MERGE(FILE_NAME, ifile)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ifile
      character(len=1000), intent(out) :: FILE_NAME
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(OUT_STATION%FNAME)
      IF (OUT_STATION%IDEF.gt.0) THEN
         WRITE (FILE_NAME,10) OUT_STATION%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4,'.nc')
      ELSE
         WRITE (FILE_NAME,20) OUT_STATION%FNAME(1:LPOS)
  20     FORMAT (a,'.nc')
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_HANDLES()
#ifdef MPI_PARALL_GRID
         use elfe_msgp
#endif
         USE DATAPOOL
         IMPLICIT NONE
         CHARACTER (LEN = 30) :: FDB
         INTEGER              :: LFDB

#ifdef SELFE
            INP%FNAME  = 'wwminput.nml'
#endif
            CHK%FNAME  = 'wwmcheck.nml'
           QSTEA%FNAME = 'qstea.out'
         WINDBG%FNAME  = 'winddbg.out'
         IOBPOUT%FNAME = 'iobp.out'
        IOBPDOUT%FNAME = 'iobpd.out'
!
!2do ... dinstinguish between binary and ascii stuff ...
!
             BND%FHNDL  = STARTHNDL + 1
             WIN%FHNDL  = STARTHNDL + 2
             CUR%FHNDL  = STARTHNDL + 3
             WAT%FHNDL  = STARTHNDL + 4
             WAV%FHNDL  = STARTHNDL + 5
             CHK%FHNDL  = STARTHNDL + 6
           HOTIN%FHNDL  = STARTHNDL + 7
          HOTOUT%FHNDL  = STARTHNDL + 8
             INP%FHNDL  = STARTHNDL + 9
             GRD%FHNDL  = STARTHNDL + 10
          GRDCOR%FHNDL  = STARTHNDL + 11

           QSTEA%FHNDL  = STARTHNDL + 12

         IOBPOUT%FHNDL  = STARTHNDL + 13
         IOBPDOUT%FHNDL = STARTHNDL + 14

         DBG%FHNDL      = STARTHNDL + 15 
         STAT%FHNDL     = STARTHNDL + 16 
         WINDBG%FHNDL   = STARTHNDL + 17 

#ifndef MPI_PARALL_GRID
         open(DBG%FHNDL,file='wwmdbg.out',status='unknown') !non-fatal errors
         open(STAT%FHNDL,file='wwmstat.out',status='unknown') !non-fatal errors
         open(WINDBG%FHNDL,file='windbg.out',status='unknown') !non-fatal errors
#else
# ifdef SELFE
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file='outputs/'//fdb,status='replace') 
# else
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file=fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file=fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file=fdb,status='replace')
# endif
#endif
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
#ifdef MPIP_PARALL_GRID
         IF (myrank == 0) THEN
#endif
         WRITE(STAT%FHNDL,*) 'Input Filename   =', TRIM(INP%FNAME)
         WRITE(STAT%FHNDL,*) 'Check Filename   =', TRIM(CHK%FNAME)
         WRITE(STAT%FHNDL,*) 'Qstea Filename   =', TRIM(QSTEA%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobp Filename    =', TRIM(IOBPOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobpd Filename   =', TRIM(IOBPDOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'WindDbg Filename =', TRIM(WINDBG%FNAME)
#ifdef MPIP_PARALL_GRID
         ENDIF
#endif

         OPEN( INP%FHNDL,      FILE = TRIM(INP%FNAME))
         OPEN( CHK%FHNDL,      FILE = TRIM(CHK%FNAME))
         OPEN( QSTEA%FHNDL,    FILE = TRIM(QSTEA%FNAME))
         OPEN( IOBPOUT%FHNDL,  FILE = TRIM(IOBPOUT%FNAME))
         OPEN( IOBPDOUT%FHNDL, FILE = TRIM(IOBPDOUT%FNAME))

           OUT1D%FHNDL = STARTHNDL + 18 
            MISC%FHNDL = STARTHNDL + 19 
         OUTSP1D%FHNDL = STARTHNDL + 20 
         OUTPARM%FHNDL = STARTHNDL + 21 
         OUTSP2D%FHNDL = STARTHNDL + 22 

         OUT%FHNDL     = STARTHNDL + 23 

      END SUBROUTINE
